/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/*=======================================================================
   //Class Problem: Handle the data and their functions for time stepping

                 and coupled processes within each time step, and finally
                 solve a problem.
   Design and implementation:  WW
   Start time:  09.07.2008
   Modification:
            12.2008 WW Incorporate the changes from previous versions.
   ========================================================================*/
#include "problem.h"

#if defined(USE_MPI)
#include <mpi.h>
#include "SplitMPI_Communicator.h"
#endif

#include <cfloat>
#include <iostream>
#include <sstream>
// kg44: max size for size_t (system_dependent) is set normally here
#include <limits>
// WW
//
/*------------------------------------------------------------------------*/
/* Pre-processor definitions */
#include "makros.h"
#include "display.h"
#include "MemWatch.h"
/*------------------------------------------------------------------------*/
// MSHLib
#include "msh_lib.h"
#include "msh_node.h"

/*------------------------------------------------------------------------*/
// Data file
// OK411
extern int ReadData(char*, GEOLIB::GEOObjects& geo_obj, std::string& unique_name);
/* PCS */
#include "pcs_dm.h"
#include "rf_pcs.h"
// 16.12.2008.WW #include "rf_apl.h"

#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
#include "par_ddc.h"
#endif

#include "rf_react.h"
#include "rf_react_int.h"
#include "rf_react_cap.h" //CB merge CAP 0311

#include "rf_st_new.h"
#include "rf_tim_new.h"
#include "rfmat_cp.h"
//#include "rf_vel_new.h"
#include "rf_fluid_momentum.h"
#include "rf_random_walk.h"
// Finite element
#include "DUMUX.h" // BG 02/2011
#include "Eclipse.h" //BG 09/2009
#include "Output.h"
#include "fem_ele_std.h"
#include "fem_ele_vec.h"
#include "files0.h" // GetLineFromFile1
#include "rf_bc_new.h"
#include "rf_node.h"
#include "rf_out_new.h"
#include "tools.h"
#include "timer.h"
#include "rf_msp_new.h" //WX:01.2013
//
#ifdef CHEMAPP
#include "eqlink.h"
#endif
#ifdef UDE_REACT
#include "rf_REACT_ude.h"
#endif
#ifdef GEM_REACT
#include "rf_REACT_GEM.h"
#endif
#ifdef BRNS
// BRNS dll link; HB 02.11.2007
#include "rf_REACT_BRNS.h"
#endif
#include "rf_kinreact.h"

#include "ShapeFunctionPool.h"

#if defined(USE_PETSC) // || defined(other parallel libs)//03.3012. WW
#include "PETSC/PETScLinearSolver.h"
#endif
namespace process
{
class CRFProcessDeformation;
}
using process::CRFProcessDeformation;

/**************************************************************************
   GeoSys - Function: Constructor
   Task:
   Programing:
   07/2008 WW Set it as an constructor of class problem based on the
            PreTimeloop
   Modification:
 ***************************************************************************/
Problem::Problem (char* filename) :
	dt0(0.), print_result(true),
	_line_shapefunction_pool(NULL), _quadr_shapefunction_pool(NULL),
	_geo_obj (new GEOLIB::GEOObjects), _geo_name (filename),
	mrank(0), msize(0)
{
	if (filename != NULL)
	{
		// read data
		ReadData(filename, *_geo_obj, _geo_name);
#if !defined(USE_PETSC) // &&  !defined(other parallel libs)//03~04.3012. WW
		DOMRead(filename);
#endif
	}
#if !defined(USE_PETSC) && !defined(NEW_EQS) // && defined(other parallel libs)//03~04.3012. WW
	//#ifndef NEW_EQS
	ConfigSolverProperties(); //_new. 19.10.2008. WW
#endif

	// set the link to Problem instance in CRFProcess objects
	for (size_t i = 0; i < pcs_vector.size(); i++)
		pcs_vector[i]->setProblemObjectPointer(this);

	for (size_t i = 0; i < pcs_vector.size(); i++)
	{
		hasAnyProcessDeactivatedSubdomains = (pcs_vector[i]->NumDeactivated_SubDomains > 0);
		if (hasAnyProcessDeactivatedSubdomains)
			break;
	}
	//----------------------------------------------------------------------
	// Create ST
	// OK STCreateFromPNT();
	//----------------------------------------------------------------------
	GetHeterogeneousFields(); // OK/MB
	//----------------------------------------------------------------------
	// Test MSH-MMP //OK
	bool validMatID = true;
	if (fem_msh_vector.size() == 1)
	{
		size_t max_matId = MSHGetMaxPatchIndex(fem_msh_vector[0]);
		validMatID = (max_matId + 1 <= mmp_vector.size());
	}
	else
	{
		int g_max_mmp_groups = MSHSetMaxMMPGroups();
		validMatID = !(g_max_mmp_groups > (int)mmp_vector.size());
	}

	if (!validMatID)
	{
		std::cout << "Error: not enough MMP data. please check MMP and material IDs in a mesh." << std::endl;
		print_result = false; // OK
		return;
	}
	//----------------------------------------------------------------------
	// Create PCS processes
	PCSCreate();
	if (!PCSCheck()) // OK4910 reactivated
	{
		print_result = false; // OK
		return;
	}

	initializeConstrainedProcesses(pcs_vector);

	CRFProcess* pcs_fluid_momentum = PCSGet("FLUID_MOMENTUM");
	if (pcs_fluid_momentum)
	{
		pcs_fluid_momentum->_idxVx = pcs_fluid_momentum->GetNodeValueIndex("VELOCITY1_X", true);
		pcs_fluid_momentum->_idxVy = pcs_fluid_momentum->GetNodeValueIndex("VELOCITY1_Y", true);
		pcs_fluid_momentum->_idxVz = pcs_fluid_momentum->GetNodeValueIndex("VELOCITY1_Z", true);
	}
	// delete pcs_fluid_momentum;

	//
	// JT: Set to true to force node copy at end of loop
	force_post_node_copy = true;
	//
	// JT: Certain restrictions might be made if an external simulator is being used
	external_coupling_exists = false;
	for (size_t i = 0; i < pcs_vector.size(); i++)
	{
		if (pcs_vector[i]->simulator.compare("GEOSYS") != 0)
			external_coupling_exists = true;
	}
#ifdef GEM_REACT
	external_coupling_exists = true;
#endif
#ifdef LIBPHREEQC
	external_coupling_exists = true;
#endif
#ifdef BRNS
	external_coupling_exists = true;
#endif
#ifdef CHEMAPP
	external_coupling_exists = true;
#endif
	//
	//......................................................................
	//#ifdef RESET_4410
	//  //if(pcs_vector[0]->pcs_type_name.compare("TWO_PHASE_FLOW")==0) //OK
	//  if(total_processes[3])  // 3: TWO_PHASE_FLOW. 12.12.2008. WW
	//    PCSCalcSecondaryVariables(); //OK
	//#endif
	//......................................................................
	// 09.07.2008 WW
	SetActiveProcesses();
	// OK if (!Check()) return; //OK
	//----------------------------------------------------------------------
	// REACTIONS
	// CB before the first time step
	if (REACTINT_vec.size() == 0)
	{
		for (size_t i = 0; i < mmp_vector.size(); i++)
		{
			if (mmp_vector[i]->porosity_model == 13)
			{
				std::cout << " Error in Model setup: Porosity model 13 is used, " << "\n";
				std::cout << " but no reaction interface is defined! Exiting now..." << "\n";
				exit(0);
			}
		}
	}
	// if(MASS_TRANSPORT_Process) // if(MASS_TRANSPORT_Process&&NAPL_Dissolution) //CB Todo
	CreateClockTime(); // CB time
	if (transport_processes.size() > 0) // 12.12.2008. WW
	{
		// set the id variable flow_pcs_type for Saturation and velocity calculation
		// in mass transport element matrices
		SetFlowProcessType();
		//----------------------------------------------------------------------
		KRConfig(*_geo_obj, _geo_name);

		// initialyse the reaction interface if  not done yet
		if (REACTINT_vec.size() > 0)
		{
			if (REACTINT_vec[0]->unitconversion)
			{
				CRFProcess* flow_pcs = NULL;
				flow_pcs = PCSGetFlow();
				if (flow_pcs->type == 1212) // in case of mutlltiphase flow, sat water must be calculated here, required
				                            // by pgc interface
					flow_pcs->CalcSecondaryVariables(true);
			}
			REACTINT_vec[0]->InitREACTINT();
		}
		//----------------------------------------------------------------------
		if (KinReactData_vector.size() > 0)
		{
			// Configure Data for Blobs (=>NAPL dissolution)
			KBlobConfig(*_geo_obj, _geo_name);
			KBlobCheck();
			// in case of Twophaseflow before the first time step
			if (total_processes[3] || total_processes[4])
				if (KNaplDissCheck()) // 3: TWO_PHASE_FLOW. 12.12.2008. WW
					KNaplCalcDensity(); // PCSCalcSecondaryVariables();
			// CB _drmc_ data for microbes
			if (MicrobeData_vector.size() > 0)
				MicrobeConfig();
		}
	}
//----------------------------------------------------------------------
// REACTIONS
// Initialization of REACT structure for rate exchange between MTM2 and Reactions

//--------------------------------------------------
// HB, for the GEM chemical reaction engine 05.2007
//--------------------------------------------------
#ifdef GEM_REACT
	m_vec_GEM = new REACT_GEM();
	GEMRead(FileName, m_vec_GEM);

	string path = ""; // to get the path of the file;
	path = FileName; // first get full path and project name;
	int pos, npos;
	pos = 0;
	npos = (int)path.size();

// Get path
#ifdef _WIN32
	pos = (int)path.rfind("\\"); // HS keep this on windows
#else
	pos = (int)path.rfind("/"); // HS keep this on linux
#endif // _WIN32
	if (pos < npos)
		path = path.substr(0, pos + 1);

	// now start initialization of GEMS
	if (m_vec_GEM->Init_Nodes(path) == 0)
	{
		if (m_vec_GEM->Init_RUN(path) == 0)
		{
			m_vec_GEM->initialized_flag = 1;
		}
		else // something is wrong and we stop execution
		{
			cout << " GEMS: Error in Init_Nodes..check input "
			     << "\n";
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
			MPI_Finalize(); // make sure MPI exits
#endif

			exit(1);
		}
	}
	else // something is wrong and we stop execution
	{
		cout << " GEMS: Error in Init_RUN..check input "
		     << "\n";
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
		MPI_Finalize(); // make sure MPI exits
#endif

		exit(1);
	}
#else // GEM_REACT
	//---------------------------------------------------
	REACT* rc = NULL; // SB
	//  rc->TestPHREEQC(); // Test if *.pqc file is present
	rc = rc->GetREACT();
	if (rc) // OK
	{
		if (rc->flag_pqc)
		{
			if (cp_vec.size() > 0)
			{ // OK
#ifdef REACTION_ELEMENT
				rc->CreateREACT(); // SB
				rc->InitREACT0();
				rc->ExecuteReactionsPHREEQC0();
				REACT_vec.clear();
				REACT_vec.push_back(rc);
#else

				rc->CreateREACT(); // SB
				rc->InitREACT();
// SB4501        rc->ExecuteReactions();
#ifdef OGS_FEM_IPQC
				rc->ExecutePQCString();
#elif LIBPHREEQC // MDL: new functions with built-in phreeqc
				rc->ExecuteReactionsPHREEQCNewLib();
#else
				rc->ExecuteReactionsPHREEQCNew();
#endif // LIBPHREEQC
				REACT_vec.clear();
				REACT_vec.push_back(rc);
#endif // REACTION_ELEMENT
			}
		}
		//  delete rc;
	}
	// CB merge CAP 0311
	// Initialize using ChemApp
	if (REACT_CAP_vec.size() > 0)
	{
		// SB 10/2009 do a first equilibrium calculation
		REACT_CAP_vec[0]->ExecuteReactionsChemApp(0, -1); // DL/SB 11/2008 //DL 2011.11.24 comment for AGU
		// Copy new timelevel to old time level
		REACT_CAP_vec[0]->ConvertIC2BC(*_geo_obj, _geo_name);
	}
#endif // GEM_REACT

#ifdef BRNS
	// Here to test BRNS; HB 02.11.2007
	// REACT_BRNS* pBRNS;
	// pBRNS = new REACT_BRNS();
	m_vec_BRNS = new REACT_BRNS();
	m_vec_BRNS->InitBRNS(this);
#endif

#ifdef CHEMAPP
	CEqlink* eq = NULL;
	eq = eq->GetREACTION();
	if (cp_vec.size() > 0 && eq) // MX
	{
		eq->TestCHEMAPPParameterFile(pcs_vector[0]->file_name_base);
		if (eq->flag_chemapp)
			eq->callCHEMAPP(pcs_vector[0]->file_name_base);
	}
#endif
	//  delete rc;

	if (REACTINT_vec.size() > 0)
		REACTINT_vec[0]->ReactionPostProcessing(true);
	//----------------------------------------------------------------------
	// DDC
	size_t no_processes = pcs_vector.size();
	CRFProcess* m_pcs = NULL;
#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
	//----------------------------------------------------------------------
	// DDC
	if (dom_vector.size() > 0)
	{
		DOMCreate();
		//
		for (size_t i = 0; i < no_processes; i++)
		{
			m_pcs = pcs_vector[i];
			m_pcs->CheckMarkedElement();
			CountDoms2Nodes(m_pcs);
			// Config boundary conditions for domain decomposition
			m_pcs->SetBoundaryConditionSubDomain(); // WW
		}
		//
		node_connected_doms.clear();
// Release some memory. WW
#if defined(USE_MPI) // TEST_MPI WW
		// Release memory of other domains. WW
		for (size_t i = 0; i < dom_vector.size(); i++)
		{
			if (i != (size_t)myrank)
			{
// If shared memory, skip the following line
#if defined(NEW_BREDUCE2)
				dom_vector[i]->ReleaseMemory();
#else
				// If MPI__Allreduce is used for all data conlection, activate following
				delete dom_vector[i];
				dom_vector[i] = NULL;
#endif
			}
		}
#endif
	}
#endif //#if !defined(USE_PETSC) // && !defined(other parallel libs)//03.3012. WW
	//----------------------------------------------------------------------
	PCSRestart(); // SB
	OUTCheck(); // new SB
	//========================================================================
	// Controls for coupling. WW
	cpl_overall_max_iterations = 1;
	cpl_overall_min_iterations = 1; // JT2012
	//========================================================================
	// WW
	char line[MAX_ZEILE];
	std::string line_string;
	std::ios::pos_type position;
	std::stringstream in_num;
	// File handling
	std::string num_file_name = FileName + NUM_FILE_EXTENSION;
	std::ifstream num_file(num_file_name.data(), std::ios::in);
	if (num_file.good())
	{
		num_file.seekg(0L, std::ios::beg);
		while (!num_file.eof())
		{
			num_file.getline(line, MAX_ZEILE);
			line_string = line;
			if (line_string.find("#STOP") != std::string::npos)
				break;
			if (line_string.find("$OVERALL_COUPLING") != std::string::npos)
			{
				in_num.str(GetLineFromFile1(&num_file));
				// JT// in_num >> max_coupling_iterations >> coupling_tolerance;
				// JT: coupling_tolerance is process dependent, cannot be here. See m_num->cpl_tolerance (with
				// $COUPLING_CONTROL)
				in_num >> cpl_overall_min_iterations >> cpl_overall_max_iterations;
				break;
			}
		}
		num_file.close();
	}
	//========================================================================
	// For time stepping. WW
	CTimeDiscretization* m_tim = NULL;
	start_time = 1.0e+25; // 1.e+8;  kg44 I need a really big time, as I have starting times bigger than 3.e+13 (1
	                      // Million years!!!)
	end_time = 0.;
	max_time_steps = 0;
	bool time_ctr = false;
	// Determine the start and end times from all available process related data.
	for (size_t i = 0; i < time_vector.size(); i++)
	{
		m_tim = time_vector[i];
		m_tim->FillCriticalTime();
		if (m_tim->time_start < start_time)
			start_time = m_tim->time_start;
		if (m_tim->time_end > end_time)
			end_time = m_tim->time_end;
		if (max_time_steps < m_tim->time_step_vector.size())
			max_time_steps = m_tim->time_step_vector.size();
		if (m_tim->GetPITimeStepCrtlType() > 0)
			time_ctr = true;
		m_tim->last_active_time = start_time; // NW

		// check maximum number of coupling iterations against maximum time step increase
		if (m_tim->time_control_type == TimeControlType::SELF_ADAPTIVE
		    && (m_tim->adapt_itr_type == IterationType::COUPLED
		        || m_tim->adapt_itr_type == IterationType::COUPLED_STABLE_ERROR)
		    && cpl_overall_max_iterations < m_tim->time_adapt_tim_vector.back())
			std::cout << "Warning: Maximum number of coupling iterations is smaller than maximum time step increase!!!"
			          << std::endl;
	}
	if (max_time_steps == 0)
		max_time_steps = std::numeric_limits<std::size_t>::max()
		                 - 1; // ULONG_MAX-1;  //kg44 increased the number to maximum number (size_t)
	current_time = start_time;
	if (time_ctr)
	{
		// Modified on 14.02.2011. WW
		long maxi_eqs_dim = 0;
		for (size_t i = 0; i < no_processes; i++)
		{
			m_pcs = pcs_vector[i];
			if (m_pcs->size_unknowns > maxi_eqs_dim)
				maxi_eqs_dim = m_pcs->size_unknowns;
		}
		buffer_array = new double[maxi_eqs_dim];
		buffer_array1 = new double[maxi_eqs_dim];
	}
	else
	{
		buffer_array = NULL;
		buffer_array1 = NULL;
	}
	//========================================================================
	CRFProcessDeformation* dm_pcs = NULL;

	//  //WW
	for (size_t i = 0; i < no_processes; i++)
	{
		m_pcs = pcs_vector[i];
		m_pcs->CalcSecondaryVariables(true); // WW
		m_pcs->Extropolation_MatValue(); // WW
	}
	// Calculation of the initial stress and released load for excavation simulation
	// 07.09.2007  WW
	// Excavation for defromation
	dm_pcs = (CRFProcessDeformation*)total_processes[12];
	if (dm_pcs)
		dm_pcs->CreateInitialState4Excavation();

#ifdef OGS_DELETE_EDGES_AFTER_INIT
	// Free memory occupied by edges. 09.2012. WW
	bool fluid_mom_pcs = false;
	for (size_t i = 0; i < no_processes; i++)
	{
		m_pcs = pcs_vector[i];
		if (m_pcs->getProcessType() == FiniteElement::FLUID_MOMENTUM) // 09.2012 WW
		{
			fluid_mom_pcs = true;
			break;
		}
	}
	if (!fluid_mom_pcs)
	{
		for (size_t k = 0; k < fem_msh_vector.size(); k++)
		{
			fem_msh_vector[k]->FreeEdgeMemory();
		}
	}
#endif
}

/**************************************************************************
   GeoSys - Function: Desstructor
   Task:
   Programing:
   08/2008 WW Set it as an constructor of class problem based on the
            PreTimeloop

   Modification:
   12.2008  WW
 ***************************************************************************/
Problem::~Problem()
{
	if (_geo_obj)
		delete _geo_obj;
	delete[] active_processes;
	delete[] exe_flag;
	if (buffer_array)
		delete[] buffer_array;
	if (buffer_array1)
		delete[] buffer_array1;
	buffer_array = NULL;
	buffer_array1 = NULL;
	active_processes = NULL;
	exe_flag = NULL;
	//
	PCSDestroyAllProcesses();
	//
	if (GetRFProcessProcessingAndActivation("MT") && GetRFProcessNumComponents() > 0)
	{
		DestroyREACT(); // SB
		cp_vec.clear(); // Destroy component properties vector
	}
//
#ifdef CHEMAPP
	if (Eqlink_vec.size() > 0)
	{
		Eqlink_vec[0]->DestroyMemory();
		Eqlink_vec.clear();
	}
#endif
// WW ClockTimeVec[0]->PrintTimes();
#ifdef GEM_REACT
	// HS:
	delete m_vec_GEM;
#endif

#ifdef BRNS
	// Here to delete BRNS instance; HB 12.11.2007
	// delete m_vec_BRNS.at(0);
	delete m_vec_BRNS;
#endif

	if (_line_shapefunction_pool)
	{
		if (_line_shapefunction_pool == _quadr_shapefunction_pool)
		{
			delete _line_shapefunction_pool;
			_line_shapefunction_pool = NULL;
			_quadr_shapefunction_pool = NULL;
		}
		else
		{
			delete _line_shapefunction_pool;
			_line_shapefunction_pool = NULL;
		}
	}
	if (_quadr_shapefunction_pool)
		delete _quadr_shapefunction_pool;

#if defined(USE_PETSC) ||defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL) || defined(USE_MPI_GEMS)
	if(mrank == 0)
#endif
		std::cout << "\n^O^: Your simulation is terminated normally ^O^ "
		          << "\n";
}

/*-------------------------------------------------------------------------
   GeoSys - Function: SetActiveProcesses
   Task:
   total_processes:
    0: LIQUID_FLOW     | 1: GROUNDWATER_FLOW  | 2: RICHARDS_FLOW
    3: PS_GLOBAL   | 4: MULTI_PHASE_FLOW  | 5: COMPONENTAL_FLOW
    6: OVERLAND_FLOW   | 7: AIR_FLOW          | 8: HEAT_TRANSPORT
    9: FLUID_MOMENTUM  |10: RANDOM_WALK       |11: MASS_TRANSPORT
   12: DEFORMATION     | 14: TNEQ             |15: TES
   Return:
   Programming:
   07/2008 WW
   03/2009 PCH added PS_GLOBAL
   Modification:
   -------------------------------------------------------------------------*/
inline int Problem::AssignProcessIndex(CRFProcess* m_pcs, bool activefunc)
{
	//	if (m_pcs->pcs_type_name.compare("OVERLAND_FLOW") == 0) {
	if (m_pcs->getProcessType() == FiniteElement::OVERLAND_FLOW)
	{
		if (!activefunc)
			return 0;
		total_processes[0] = m_pcs;
		active_processes[0] = &Problem::OverlandFlow;
		return 0;
		//	} else if (m_pcs->pcs_type_name.compare("AIR_FLOW") == 0) {
	}
	else if (m_pcs->getProcessType() == FiniteElement::GROUNDWATER_FLOW)
	{
		if (!activefunc)
			return 1;
		total_processes[1] = m_pcs;
		active_processes[1] = &Problem::GroundWaterFlow;
		return 1;
		//	} else if (m_pcs->pcs_type_name.compare("RICHARDS_FLOW") == 0) {
	}
	else if (m_pcs->getProcessType() == FiniteElement::RICHARDS_FLOW)
	{
		if (!activefunc)
			return 2;
		total_processes[2] = m_pcs;
		active_processes[2] = &Problem::RichardsFlow;
		return 2;
		//	} else if (m_pcs->pcs_type_name.compare("TWO_PHASE_FLOW") == 0) {
	}
	else if (m_pcs->getProcessType() == FiniteElement::TWO_PHASE_FLOW)
	{
		if (!activefunc)
			return 3;
		total_processes[3] = m_pcs;
		active_processes[3] = &Problem::TwoPhaseFlow;
		return 3;
		//	} else if (m_pcs->pcs_type_name.compare("MULTI_PHASE_FLOW") == 0) {
	}
	else if (m_pcs->getProcessType() == FiniteElement::MULTI_PHASE_FLOW)
	{
		if (!activefunc)
			return 4;
		total_processes[4] = m_pcs;
		active_processes[4] = &Problem::MultiPhaseFlow;
		return 4;
		//	} else if (m_pcs->pcs_type_name.compare("COMPONENTAL_FLOW") == 0) {
		//	} else if (m_pcs->getProcessType() == COMPONENTAL_FLOW) {
		//		if (!activefunc)
		//			return 5;
		//		total_processes[5] = m_pcs;
		//		active_processes[5] = &Problem::ComponentalFlow;
		//		return 5;
		//	} else if (m_pcs->pcs_type_name.compare("OVERLAND_FLOW") == 0) {
	}
	else if (m_pcs->getProcessType() == FiniteElement::LIQUID_FLOW)
	{
		if (!activefunc)
			return 6;
		total_processes[6] = m_pcs;
		active_processes[6] = &Problem::LiquidFlow;
		return 6;
		//	} else if (m_pcs->pcs_type_name.compare("GROUNDWATER_FLOW") == 0) {
	}
	else if (m_pcs->getProcessType() == FiniteElement::AIR_FLOW)
	{
		if (!activefunc)
			return 7;
		total_processes[7] = m_pcs;
		active_processes[7] = &Problem::AirFlow;
		return 7;
		//	} else if (m_pcs->pcs_type_name.compare("HEAT_TRANSPORT") == 0) {
	}
	else if (m_pcs->getProcessType() == FiniteElement::HEAT_TRANSPORT)
	{
		if (!activefunc)
			return 8;
		total_processes[8] = m_pcs;
		active_processes[8] = &Problem::HeatTransport;
		return 8;
		//	} else if (m_pcs->pcs_type_name.compare("FLUID_MOMENTUM") == 0) {
	}
	else if (m_pcs->getProcessType() == FiniteElement::FLUID_MOMENTUM)
	{
		if (!activefunc)
			return 9;
		total_processes[9] = m_pcs;
		active_processes[9] = &Problem::FluidMomentum;
		return 9;
		//	} else if (m_pcs->pcs_type_name.compare("RANDOM_WALK") == 0) {
	}
	else if (m_pcs->getProcessType() == FiniteElement::RANDOM_WALK)
	{
		if (!activefunc)
			return 10;
		total_processes[10] = m_pcs;
		active_processes[10] = &Problem::RandomWalker;
		DATWriteParticleFile(0); // YS
		return 10;
		//	} else if (m_pcs->pcs_type_name.compare("MASS_TRANSPORT") == 0) {
	}
	else if (m_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT)
	{
		if (!activefunc)
			return 11;
		total_processes[11] = m_pcs;
		active_processes[11] = &Problem::MassTrasport;
		return 11;
		//	} else if (m_pcs->pcs_type_name.find("DEFORMATION") != string::npos) {
	}
	else if (isDeformationProcess(m_pcs->getProcessType()))
	{
		if (!activefunc)
			return 12;
		total_processes[12] = m_pcs;
		active_processes[12] = &Problem::Deformation;
		return 12;
		//	} else if (m_pcs->pcs_type_name.find("PS_GLOBAL") != string::npos) {
	}
	else if (m_pcs->getProcessType() == FiniteElement::PS_GLOBAL)
	{
		//    if(!activefunc) return 13;
		if (!activefunc)
			return 3;
		total_processes[3] = m_pcs;
		active_processes[3] = &Problem::PS_Global;
		return 3;
	}
	else if (m_pcs->getProcessType() == FiniteElement::MULTI_COMPONENTIAL_FLOW)
	{
		if (!activefunc)
			return 5;
		total_processes[5] = m_pcs;
		active_processes[5] = &Problem::MULTI_COMPONENTIAL_FLOW;
		return 5;
	}
	else if (m_pcs->getProcessType() == FiniteElement::TNEQ)
	{
		if (!activefunc)
			return 14;
		total_processes[14] = m_pcs;
		active_processes[14] = &Problem::TNEQ;
		return 14;
	}
	else if (m_pcs->getProcessType() == FiniteElement::TES)
	{
		if (!activefunc)
			return 15;
		total_processes[15] = m_pcs;
		active_processes[15] = &Problem::TES;
		return 15;
	}
	std::cout << "Error: no process is specified. " << '\n';
	return -1;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: SetActiveProcesses
   Task:
   total_processes:
    0: LIQUID_FLOW     | 1: GROUNDWATER_FLOW  | 2: RICHARDS_FLOW
    3: TWO_PHASE_FLOW  | 4: MULTI_PHASE_FLOW  | 5: COMPONENTAL_FLOW
    6: OVERLAND_FLOW   | 7: AIR_FLOW          | 8: HEAT_TRANSPORT
    9: FLUID_MOMENTUM  |10: RANDOM_WALK       |11: MASS_TRANSPORT
   12: DEFORMATION     |13: PS_GLOBAL         |14: TNEQ              | 15: TES
   Return:
   Programming:
   07/2008 WW
   03/2009 PCH add PS_GLOBAL
   Modification:
   --------------------------------------------------------------------*/
void Problem::SetActiveProcesses()
{
	CRFProcess* m_pcs = NULL;
	total_processes.resize(max_processes);
	active_processes = new ProblemMemFn[max_processes];
	coupled_process_index.resize(max_processes);
	exe_flag = new bool[max_processes];
	//
	for (size_t i = 0; i < max_processes; i++)
	{
		total_processes[i] = NULL;
		active_processes[i] = NULL;
		coupled_process_index[i] = -1;
	}
	//
	for (size_t i = 0; i < pcs_vector.size(); i++)
	{
		m_pcs = pcs_vector[i];
		AssignProcessIndex(m_pcs);
	}
	//
	for (size_t i = 0; i < max_processes; i++)
		if (total_processes[i])
		{
			// JT: Check for a coupled variable or process (variable is probably necessary for multiple component
			// transport situations.
			// First try for a variable, because this is more general for mass transport
			m_pcs = PCSGet(total_processes[i]->m_num->cpl_variable, true);
			if (m_pcs)
			{
				m_pcs->pcs_is_cpl_underling = true;
				total_processes[i]->pcs_is_cpl_overlord = true;
				//
				coupled_process_index[i] = AssignProcessIndex(m_pcs, false);
				m_pcs->cpl_overlord = total_processes[i];
				total_processes[i]->cpl_underling = m_pcs;
			}
			else
			{ // Try for a process, because it may have been assigned this way
				m_pcs = PCSGet(total_processes[i]->m_num->cpl_process);
				if (m_pcs)
				{
					m_pcs->pcs_is_cpl_underling = true;
					total_processes[i]->pcs_is_cpl_overlord = true;
					//
					coupled_process_index[i] = AssignProcessIndex(m_pcs, false);
					m_pcs->cpl_overlord = total_processes[i];
					total_processes[i]->cpl_underling = m_pcs;
				}
			}
			active_process_index.push_back(i);
		}
	// Transport  porcesses
	for (size_t k = 0; k < pcs_vector.size(); k++)
	{
		m_pcs = pcs_vector[k];
		//		if (m_pcs->pcs_type_name.compare("MASS_TRANSPORT") == 0)
		// TF
		if (m_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT)
			transport_processes.push_back(m_pcs);
		//		if (m_pcs->pcs_type_name.compare("TWO_PHASE_FLOW") == 0) //09.01.2008. WW
		// TF
		if ((m_pcs->getProcessType() == FiniteElement::TWO_PHASE_FLOW)
		    || (m_pcs->getProcessType() == FiniteElement::MULTI_PHASE_FLOW)) // BG 04/2011
			multiphase_processes.push_back(m_pcs);

		if ((m_pcs->getProcessType() == FiniteElement::GROUNDWATER_FLOW)
		    || (m_pcs->getProcessType() == FiniteElement::LIQUID_FLOW)) // BG 04/2011
			singlephaseflow_process.push_back(m_pcs);
	}
}

/**************************************************************************     <
   ROCKFLOW - Function: PCSCreate
   Task:
   Programing:
   02/2003 OK Implementation
   03/2003 OK H processes
   04/2003 OK C processes
   05/2003 OK T processes
   05/2003 OK TH processes
   08/2004 OK PCS2
   08/2004 WW The new creation of the deformation process
   10/2004 OK H gas processes
   01/2005 WW New element calculation
   01/2005 OK H unsatutated process
   02/2005 MB switch case in config()
   06/2005 OK MMP2PCSRelation
   07/2008 WW Capsulated into class Problem
   Modification:
   05/2010 TF formated source code
 ***************************************************************************/
void Problem::PCSCreate()
{
#if defined(USE_PETSC) || defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL) \
    || defined(USE_MPI_GEMS)
	if (mrank == 0)
	{
#endif
		std::cout << "---------------------------------------------"
		          << "\n";
		std::cout << "Create PCS processes"
		          << "\n";
#if defined(USE_PETSC) || defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL) \
    || defined(USE_MPI_GEMS)
	}
#endif

	size_t no_processes = pcs_vector.size();
	// OK_MOD if(pcs_deformation>0) Init_Linear_Elements();
	for (size_t i = 0; i < no_processes; i++)
	{
		pcs_vector[i]->pcs_type_number = i;
		pcs_vector[i]->Config(); // OK
	}

#if defined(NEW_EQS)
	CreateEQS_LinearSolver();
#endif

	for (size_t i = 0; i < no_processes; i++)
	{
#if defined(USE_PETSC) || defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL) \
    || defined(USE_MPI_GEMS)
		if (mrank == 0)
		{
#endif
			std::cout << "............................................."
			          << "\n";
			FiniteElement::ProcessType pcs_type(pcs_vector[i]->getProcessType());
			std::cout << "Create: " << FiniteElement::convertProcessTypeToString(pcs_type) << "\n";
			//		if (!pcs_vector[i]->pcs_type_name.compare("MASS_TRANSPORT")) {
			// YS   // TF
			// if (pcs_type != FiniteElement::MASS_TRANSPORT && pcs_type != FiniteElement::FLUID_MOMENTUM
			//				&& pcs_type != FiniteElement::RANDOM_WALK)
			if (pcs_type == FiniteElement::MASS_TRANSPORT) // SB
			{
				std::cout << " for " << pcs_vector[i]->pcs_primary_function_name[0] << " ";
				std::cout << " pcs_component_number " << pcs_vector[i]->pcs_component_number;
			}
			std::cout << "\n";
#if defined(USE_PETSC) || defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL) \
    || defined(USE_MPI_GEMS)
		}
#endif

		pcs_vector[i]->Create();
	}

	createShapeFunctionPool(); // WW

	for (size_t i = 0; i < no_processes; i++)
	{ // WW
		CRFProcess* pcs = pcs_vector[i];
		pcs->SetBoundaryConditionAndSourceTerm();

		if (   pcs->getProcessType() == FiniteElement::DEFORMATION_DYNAMIC
			|| pcs->getProcessType() == FiniteElement::DEFORMATION
			|| pcs->getProcessType() == FiniteElement::DEFORMATION_FLOW
			|| pcs->getProcessType() == FiniteElement::DEFORMATION_H2)
		{
			CRFProcessDeformation* dm_pcs = dynamic_cast<CRFProcessDeformation*>(pcs);
			dm_pcs->InitGauss();
		}
	}

#if defined(USE_PETSC) // || defined(other solver libs)//03.3012. WW
	CreateEQS_LinearSolver();
#endif

	for (size_t i = 0; i < no_processes; i++)
		MMP2PCSRelation(pcs_vector[i]);

	for (size_t i = 0; i < no_processes; i++) // WW

		pcs_vector[i]->ConfigureCouplingForLocalAssemblier();

	for (size_t i = 0; i < out_vector.size(); i++)
		// initialize process and mesh attributes of COutput objects
		out_vector[i]->init();
}

/*-------------------------------------------------------------------------
   ROCKFLOW - Function: PCSRestart
   Task: Insert process to list
   Programming:
   06/2003 OK Implementation
   07/2008 WW Capsulated into class Problem
   Modification:
   -------------------------------------------------------------------------*/
void Problem::PCSRestart()
{
	const size_t no_processes(pcs_vector.size());
	if (no_processes == 0)
		return; // OK41

	int ok = 0; // = ReadRFRRestartData(file_name_base);

	if (ok == 0)
	{
		std::cout << "RFR: no restart data"
		          << "\n";
		return;
	}

	// WW int nidx0; //, nidx1;
	for (size_t i = 0; i < no_processes; i++)
	{
		CRFProcess* m_pcs = pcs_vector[i];
		for (size_t j = 0; j < m_pcs->GetPrimaryVNumber(); j++)
		{
			// timelevel=0;
			// WW nidx0 = m_pcs->GetNodeValueIndex(m_pcs->GetPrimaryVName(j));
			// timelevel= 1;
			// WW nidx1 = nidx0 + 1;
			// OK411      CopyNodeVals(nidx1,nidx0);
		}
	}
}

/*-------------------------------------------------------------------------
   GeoSys - Function: SetTimeActiveProcesses
   Task: For independent time stepping. Set processes to active in time.
   Note: Minimum time step allowance in handled in CalcTimeStep()
   Return:
   Programming:
   03/2012 JT
   Modification:
   --------------------------------------------------------------------*/
void Problem::SetTimeActiveProcesses()
{
	size_t ii;
	double tval, next_time, lowest_next_active;
	CTimeDiscretization* m_tim = NULL;
	CTimeDiscretization* inactive_tim = NULL;
	next_time = current_time + dt;
	//
	lowest_next_active = DBL_MAX;
	for (ii = 0; ii < active_process_index.size(); ii++)
	{
		m_tim = total_processes[active_process_index[ii]]->Tim;
		m_tim->time_active = true; // activate
		if (m_tim->time_independence && m_tim->next_active_time > next_time)
		{ // Process is then not active this time step
			m_tim->time_active = false; // deactivate
			//
			// store the lowest next active time of inactive processes
			if (m_tim->next_active_time < lowest_next_active)
			{
				lowest_next_active = m_tim->next_active_time;
				inactive_tim = m_tim;
			}
		}
	}
	//
	// Check if we should shift time step slightly to hit a non-active process
	if (inactive_tim)
	{
		tval = next_time + dt / 1.0e3; // a small dt increase is better than a miniscule dt on the next step
		if (tval > lowest_next_active)
		{ // allow this slight increase
			inactive_tim->time_active = true;
			dt = lowest_next_active - current_time;
			next_time = current_time + dt;
		}
		else if (tval + dt > lowest_next_active)
		{ // Try to smooth to the target time from 2 time steps away
			dt = (lowest_next_active - current_time) / 2.0;
			next_time = current_time + dt;
		}
	}
	//
	// Set times for all active processes
	for (ii = 0; ii < active_process_index.size(); ii++)
	{
		m_tim = total_processes[active_process_index[ii]]->Tim;
		if (m_tim->time_active)
		{
			m_tim->time_step_length = next_time - m_tim->last_active_time;
			m_tim->last_active_time = next_time;
		}
	}
}

/**************************************************************************
   FEMLib-Method:
   07/2008 WW Implementation
   01/2009 WW Update
   03/2012 JT Many changes. Allow independent time stepping.
**************************************************************************/
void Problem::Euler_TimeDiscretize()
{
	long accepted_times = 0;
	long rejected_times = 0;
	double dt_rec;
	int i;
	bool force_output;
	last_dt_accepted = false; // JT: false first. Thus copy node values after first dt.
	//
	CTimeDiscretization* m_tim = NULL;
	aktueller_zeitschritt = 0;
	ScreenMessage("\n\n***Start time steps\n");
//
// Output zero time initial values
#if defined(USE_MPI) || defined(USE_MPI_KRC)
	if (mrank == 0)
	{
#endif
		OUTData(current_time, aktueller_zeitschritt, true);
#if defined(USE_MPI) || defined(USE_MPI_KRC)
	}
#endif

	// check if this is a steady state simulation
	bool isSteadySimulation = true;
	for (i = 0; i < (int)active_process_index.size(); i++)
	{
		if (total_processes[active_process_index[i]]->tim_type != TimType::STEADY)
		{
			isSteadySimulation = false;
			break;
		}
	}

	//
	// ------------------------------------------
	// PERFORM TRANSIENT SIMULATION
	// ------------------------------------------
	double previous_rejected_dt = .0;
	while (end_time > current_time)
	{
		// Get time step
		dt = dt_rec = DBL_MAX;
		for (i = 0; i < (int)active_process_index.size(); i++)
		{
			m_tim = total_processes[active_process_index[i]]->Tim;
			if (!m_tim->time_active)
				continue; // JT
			dt = MMin(dt, m_tim->CalcTimeStep(current_time));
			dt_rec = MMin(dt_rec, m_tim->recommended_time_step); // to know if we have a critical time alteration
		}

		if (!last_dt_accepted && dt == previous_rejected_dt)
		{
			ScreenMessage("Stop this simulation. New time step size is same as the rejected one.\n");
			break;
		}

		SetTimeActiveProcesses(); // JT2012: Activate or deactivate processes with independent time stepping
//
#if defined(USE_MPI)
		MPI_Bcast(
		    &dt, 1, MPI_DOUBLE, 0, comm_DDC); // all processes use the same time stepping (JT->WW. Must they always?)
#endif
		//
		// Update time settings
		aktueller_zeitschritt++;
		current_time += dt;
		aktuelle_zeit = current_time;
//
// Print messsage
#if defined(USE_PETSC) || defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL) \
    || defined(USE_MPI_GEMS)
		if (mrank == 0)
		{
#endif
			std::cout << "\n\n#############################################################\n";
			std::cout << "Time step: " << aktueller_zeitschritt << "|  Time: " << current_time
			          << "|  Time step size: " << dt << "\n";
			if (dt_rec > dt)
			{
				std::cout << "This time step size was modified to match a critical time!"
				          << "\n";
			}
#if defined(USE_PETSC) || defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL) \
    || defined(USE_MPI_GEMS)
		}
#endif
		if (CouplingLoop())
		{
			// ---------------------------------
			// TIME STEP ACCEPTED
			// ---------------------------------
			last_dt_accepted = true;
#if defined(USE_PETSC) || defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL) \
    || defined(USE_MPI_GEMS)
			if (mrank == 0)
#endif
				ScreenMessage("This step is accepted.\n");
			PostCouplingLoop();
			if (print_result)
			{
				if (current_time < end_time)
					force_output = false;
				else // JT: Make sure we printout on last time step
					force_output = true;
#if defined(USE_MPI) || defined(USE_MPI_KRC)
				if (myrank == 0)
				{
#endif
					//
					OUTData(current_time, aktueller_zeitschritt, force_output);
#if defined(USE_MPI)
				}
#endif
			}
			accepted_times++;
			for (i = 0; i < (int)active_process_index.size(); i++)
			{
				m_tim = total_processes[active_process_index[i]]->Tim;
				if (m_tim->time_active)
					m_tim->accepted_step_count++;
			}
		}
		else if (isSteadySimulation)
		{
			ScreenMessage(
			    "This time step is rejected. We stop the simulation because this is steady state simulation.\n");
			break;
		}
		else
		{
			// ---------------------------------
			// TIME STEP FAILED
			// ---------------------------------
			last_dt_accepted = false;
			ScreenMessage("This step is rejected: Redo, with a new time step.\n");
			rejected_times++;
			current_time -= dt;
			aktuelle_zeit = current_time;
			aktueller_zeitschritt--;
			previous_rejected_dt = dt;
			//
			// decrement active dt, and increment count
			for (i = 0; i < (int)active_process_index.size(); i++)
			{
				m_tim = total_processes[active_process_index[i]]->Tim;
				if (!m_tim->time_active)
					continue;
				m_tim->rejected_step_count++;
				m_tim->last_active_time -= dt;
				m_tim->step_current--;
				m_tim->repeat = true;
				m_tim->last_rejected_timestep = aktueller_zeitschritt + 1;
				//
				// Copy nodal values in reverse
				if (isDeformationProcess(total_processes[active_process_index[i]]->getProcessType()))
					continue;
				total_processes[active_process_index[i]]->CopyTimestepNODValues(false);
				// JT: This wasn't done before. Is it needed? //
				// total_processes[active_process_index[i]]->CopyTimestepELEValues(false);
			}
			for (i = 0; i < (int)total_processes.size(); i++)
			{
				if (!active_processes[i] && total_processes[i] && total_processes[i]->tim_type == TimType::STEADY)
				{
					m_tim = total_processes[i]->Tim;
					m_tim->step_current--;
					m_tim->repeat = true;
				}
			}
		}
		ScreenMessage("\n#############################################################\n");
		if (aktueller_zeitschritt >= max_time_steps)
			break;

		//		// executing only one time step for profiling
		//		current_time = end_time;
	}
#if defined(USE_PETSC) || defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL) \
    || defined(USE_MPI_GEMS) || defined(USE_MPI_KRC)
	if (mrank == 0)
	{
#endif
		std::cout << "\n----------------------------------------------------\n";
		for (i = 0; i < (int)active_process_index.size(); i++) // JT2012
		{
			m_tim = total_processes[active_process_index[i]]->Tim;
			std::cout << "\nFor process: "
			          << convertProcessTypeToString(total_processes[active_process_index[i]]->getProcessType()) << "\n";
			if (m_tim->time_control_type == TimeControlType::FIXED_STEPS)
			{
				std::cout << "No time control for this process."
				          << "\n";
			}
			else
			{
				std::cout << "Accepted time steps:                " << m_tim->accepted_step_count << "\n";
				std::cout << "Rejected time steps:                " << m_tim->rejected_step_count << "\n";
			}
			if (total_processes[active_process_index[i]]->m_num->nls_max_iterations > 1)
			{
				std::cout << "Number of non-converged iterations: "
				          << total_processes[active_process_index[i]]->num_notsatisfied << "\n";
				std::cout << "Number of stagnated iterations:     "
				          << total_processes[active_process_index[i]]->num_diverged << "\n";
			}
		}
		std::cout << "\n----------------------------------------------------\n";
#if defined(USE_PETSC) || defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL) \
    || defined(USE_MPI_GEMS) || defined(USE_MPI_KRC)

#if defined(USE_PETSC) // 05.2014. WW
		for (size_t i = 0; i < out_vector.size(); i++)
		{
			out_vector[i]->NODDomainWriteBinary_Header();
		}
#endif
	}
#endif
}

/*-----------------------------------------------------------------------
   GeoSys - Function: Coupling loop
   Task:
   Return: error
   Programming:
   07/2008 WW
   Modification:
   12.2008 WW Update
   03.2012 JT All new. Different strategy, generalized tolerance criteria
   -------------------------------------------------------------------------*/
bool Problem::CouplingLoop()
{
	int i, index, cpl_index;
	double max_outer_error, max_inner_error; //, error;
	bool transient_bc = false;
	bool run_flag[max_processes];
	int outer_index, inner_index, inner_max; //, inner_min;
	//
	CRFProcess* a_pcs = NULL;
	CRFProcess* b_pcs = NULL;
	CTimeDiscretization* m_tim = NULL;
	//
	print_result = false;
	int acounter = 0;
	//
	for (i = 0; i < (int)pcs_vector.size(); i++)
	{
		pcs_vector[i]->UpdateTransientBC();
		if (pcs_vector[i]->bc_transient_index.size() != 0)
			transient_bc = true;
	}
	if (transient_bc)
		pcs_vector[0]->WriteBC();

	for (i = 0; i < (int)total_processes.size(); i++)
	{
		if (active_processes[i] && total_processes[i]->selected)
		{
			exe_flag[i] = true;
			m_tim = total_processes[i]->Tim;
			total_processes[i]->SetDefaultTimeStepAccepted();
			acounter++;
			m_tim->step_current++;
			// reset
			total_processes[i]->iter_nlin_max = 0;
			total_processes[i]->iter_lin_max = 0;
		}
		else
		{ // 21.05.2010.  WW
			if (total_processes[i] && total_processes[i]->tim_type == TimType::STEADY)
			{
				acounter++;
				m_tim = total_processes[i]->Tim;
				m_tim->step_current++; // NW increment needed to get correct time step length in
				// CTimeDiscretization::CalcTimeStep()
			}
			exe_flag[i] = false;
		}
	}
	int num_processes = (int)active_process_index.size();
	//
	// JT: All active processes must run on the overall loop. Strange this wasn't the case before.
	for (i = 0; i < (int)total_processes.size(); i++)
	{
		run_flag[i] = exe_flag[i];
	}
	for (i = 0; i < num_processes; i++)
	{
		index = active_process_index[i];
		total_processes[index]->first_coupling_iteration = true;
	}
	//
	// To do
	// SB->WW I do not understand this condition, why switch off output?
	// WW Reason:
	/// Make output when all defined processes are activated.
	// JT->WW->SB:  I agree with SB. Just b/c one process is deactivated doesn't mean we don't want output for the
	// others.
	// if(acounter == num_processes)
	print_result = true;
	//
	bool accept = true;
	max_outer_error = 0.0;
	for (outer_index = 0; outer_index < cpl_overall_max_iterations; outer_index++)
	{
		// JT: All active processes must run on the overall loop. Strange this wasn't the case before.
		for (i = 0; i < num_processes; i++)
		{
			index = active_process_index[i];
			run_flag[index] = exe_flag[index];
		}
		for (i = 0; i < num_processes; i++)
		{
			index = active_process_index[i];
			m_tim = total_processes[index]->Tim;
			if (!m_tim->time_active)
				run_flag[index] = false;
		}

		max_outer_error = 0.0; // NW reset error for each iteration
		for (i = 0; i < num_processes; i++)
		{
			index = active_process_index[i];
			if (!run_flag[index])
				continue; // JT: may have been turned off after an inner loop!
			cpl_index = coupled_process_index[index];
			//
			// PERFORM AN INNER COUPLING
			// ---------------------------------------
			if (cpl_index >= 0 && run_flag[cpl_index])
			{
				a_pcs = total_processes[index];
				b_pcs = total_processes[cpl_index];
				//
				inner_max = a_pcs->m_num->cpl_max_iterations;
				//				inner_min = a_pcs->m_num->cpl_min_iterations; // variable set but never used
				//
				a_pcs->iter_outer_cpl = outer_index;
				b_pcs->iter_outer_cpl = outer_index;
				//
				max_inner_error = 0.0;
				for (inner_index = 0; inner_index < a_pcs->m_num->cpl_max_iterations; inner_index++)
				{
					a_pcs->iter_inner_cpl = inner_index;
					b_pcs->iter_inner_cpl = inner_index;
					//
					// FIRST PROCESS
					loop_process_number = i;
					if (a_pcs->first_coupling_iteration)
						PreCouplingLoop(a_pcs);
					//					error = Call_Member_FN(this, active_processes[index])();
					Call_Member_FN(this, active_processes[index])();
					if (!a_pcs->TimeStepAccept())
					{
						accept = false;
						break;
					}
					//
					// COUPLED PROCESS
					loop_process_number = i + 1;
					if (b_pcs->first_coupling_iteration)
						PreCouplingLoop(b_pcs);
					//					error = Call_Member_FN(this, active_processes[cpl_index])();
					Call_Member_FN(this, active_processes[cpl_index])();
					if (!b_pcs->TimeStepAccept())
					{
						accept = false;
						break;
					}
					//
					// Check for break criteria
					max_inner_error = MMax(a_pcs->cpl_max_relative_error, b_pcs->cpl_max_relative_error);
					a_pcs->first_coupling_iteration
					    = false; // No longer true (JT: these are important, and are also used elswhere).
					b_pcs->first_coupling_iteration = false; // No longer true.
					//
					// Store the outer loop error
					if (inner_index == 0)
						max_outer_error = MMax(max_outer_error, max_inner_error);
					//
					std::cout << "\n======================================================\n";
					std::cout << "Inner coupling loop " << inner_index + 1 << "/" << inner_max << " complete."
					          << "\n";
					std::cout << "Max coupling error (relative to tolerance): " << max_inner_error << "\n";
					std::cout << "======================================================\n";
					//
					// Coupling convergence criteria (use loop minimum from a_pcs because this is where the coupled
					// process was called)
					if (max_inner_error <= 1.0
					    && inner_index + 2
					           > a_pcs->m_num->cpl_min_iterations) // JT: error is relative to the tolerance.
						break;
				}
				run_flag[cpl_index] = false; // JT: CRUCIAL!!
			}
			else
			{
				// PERFORM AN OUTER COUPLING
				// ---------------------------------------
				a_pcs = total_processes[index];
				a_pcs->iter_outer_cpl = outer_index;
				a_pcs->iter_inner_cpl = 0;
				//
				loop_process_number = i;
				if (a_pcs->first_coupling_iteration)
					PreCouplingLoop(a_pcs);
				//				error = Call_Member_FN(this, active_processes[index])(); // TF: error set, but never
				//used
				Call_Member_FN(this, active_processes[index])();
				if (!a_pcs->TimeStepAccept())
				{
					accept = false;
					break;
				}
				a_pcs->first_coupling_iteration = false; // No longer true.
				// Check for break criteria
				max_outer_error = MMax(max_outer_error, a_pcs->cpl_max_relative_error);

// Reapply BCs if constrained BC
#if defined(USE_MPI) || defined(USE_PETSC)
				int has_constrained_bc_i = a_pcs->hasConstrainedBC() ? 1 : 0;
				int n_has_constrained_bc = 0;
#ifdef USE_PETSC
				MPI_Comm comm = MPI_COMM_WORLD;
#else
				MPI_Comm comm = comm_DDC;
#endif
				MPI_Allreduce(&has_constrained_bc_i, &n_has_constrained_bc, 1, MPI_INT, MPI_SUM, comm);
				if (n_has_constrained_bc > 0)
/* Commented out the following codes because MPI_C_BOOL is supported from the MPI 2.2
bool has_constrained_bc_i = a_pcs->hasConstrainedBC();
bool has_constrained_bc = false;
MPI_Allreduce(&has_constrained_bc_i, &has_constrained_bc, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);
if(has_constrained_bc > 0)
*/
#else
				if (a_pcs->hasConstrainedBC())
#endif
				{
#if defined(USE_MPI) || defined(USE_PETSC)
					const int rank = mrank;
#else
					const int rank = -1;
#endif
					if (rank < 1)
						std::cout << "this process has constrained BCs. reapply BCs." << std::endl;
					a_pcs->IncorporateBoundaryConditions(rank);
				}
			}
			if (!accept)
				break;
		}
		if (!accept)
		{
			std::cout << "\n";
			break;
		}
		//
		if (cpl_overall_max_iterations > 1)
		{
			std::cout << "\n======================================================\n";
			std::cout << "Outer coupling loop " << outer_index + 1 << "/" << cpl_overall_max_iterations << " complete."
			          << "\n";
			std::cout << "Max coupling error (relative to tolerance): " << max_outer_error << "\n";
			std::cout << "======================================================\n";

			// update max cpl error variable in all processes
			for (std::size_t i(0); i < active_process_index.size(); i++)
			{
				std::size_t thisPCSIndex(active_process_index[i]);
				CRFProcess* thisPCS(total_processes[thisPCSIndex]);
				thisPCS->cpl_max_relative_error_overall = max_outer_error;
			}
		}
		else
		{
			std::cout << "\n";
		}
		// Coupling convergence criteria
		// if(max_outer_error <= 1.0 && outer_index+2 > cpl_overall_min_iterations) // JT: error is relative to the
		// tolerance.	//MW outer_index + 2 is hard to follow - is it faster to write it like this?
		if (max_outer_error <= 1.0
		    && outer_index + 1 >= cpl_overall_min_iterations) // JT: error is relative to the tolerance.
			break;

		// MW
		if (max_outer_error > 1 && outer_index + 1 == cpl_overall_max_iterations
		    && cpl_overall_max_iterations > 1) // m_tim->step_current>1 &&
		{
			accept = false;
			break;
		}
	}
	//
	return accept;
}

/*-----------------------------------------------------------------------
   GeoSys - Function: pre Coupling loop
   Task: Process solution is beginning. Perform any pre-loop configurations
   Programming:
   03/2012 JT
   Modification:
-------------------------------------------------------------------------*/
void Problem::PreCouplingLoop(CRFProcess* m_pcs)
{
	if (!last_dt_accepted || force_post_node_copy) // if last time step not accepted or values were already copied.
		return;
	//
	/*For mass transport this routine is only called once (for the overall transport process)
	  and so we need to copy for all transport components*/
	if (m_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT)
	{
		CRFProcess* c_pcs = NULL;
		for (size_t i = 0; i < pcs_vector.size(); i++)
		{
			c_pcs = pcs_vector[i];
			if (c_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT)
			{
				c_pcs->CopyTimestepNODValues();
				c_pcs->CopyTimestepELEValues();
			}
		}
	}
	else
	{ // Otherwise, just copy this process
		m_pcs->CopyTimestepNODValues();
		m_pcs->CopyTimestepELEValues();
	}
}

/*-----------------------------------------------------------------------
   GeoSys - Function: post Coupling loop
   Task:
   Return: error
   Programming:
   08/2008 WW
   Modification:
   -------------------------------------------------------------------------*/
void Problem::PostCouplingLoop()
{
	CRFProcess* m_pcs = NULL;
	bool prqvec = false;
	bool capvec = false;
	if (total_processes[12])
	{
		CRFProcessDeformation* dm_pcs = (CRFProcessDeformation*)(total_processes[12]);
		bool doPostExcav = false; // WX
		for (size_t l = 0; l < msp_vector.size(); l++)
		{
			if (msp_vector[l]->GetBoolExcavated())
				doPostExcav = true;
		}
		if (dm_pcs->ExcavMaterialGroup >= 0 || doPostExcav)
			dm_pcs->PostExcavation(); // WX:07.2011
		if (dm_pcs->UpdateIniState == 1) // WX:10.2011
			dm_pcs->UpdateIniStateValue();

		if (H_Process && dm_pcs->type / 10 != 4) // HM partitioned scheme
			dm_pcs->ResetTimeStep();
		dm_pcs->Extropolation_GaussValue();
		dm_pcs->WriteGaussPointStress();
	}

	// Reaction postprocessing
	if (REACT_CAP_vec.size() > 0)
		capvec = true;

	if ((KinReactData_vector.size() > 0) || (REACT_vec.size() > 0) || capvec || prqvec)
	{
		// map concentrations in radial model
		if (KinReactData_vector.size() > 0)
			if (KinReactData_vector[0]->copy_concentrations)
				KinReactData_vector[0]->CopyConcentrations();
		if (transport_processes.size() > 0)
		{
			if (total_processes[3] || total_processes[4])
				if (KNaplDissCheck())
				{ // Check if NAPLdissolution is modeled
					// return P_new, and Phase_volumina
					CalcNewPhasePressure();
				}
		}
		// update porosities and permeabilities from reactions
		if (REACTINT_vec.size() > 0)
			REACTINT_vec[0]->ReactionPostProcessing(false);
		if (KinReactData_vector.size() > 0)
		{
			if (KinReactData_vector[0]->NumberMineralkinetics > 0)
				KinReactData_vector[0]
				    ->PostprocessMinKin(); // Set Mineral surface areas for next time step based on deltaC
		}
		m_pcs = PCSGetFlow();
		m_pcs->Extropolation_GaussValue();
	}

	//  Update the results
	for (int i = 0; i < (int)pcs_vector.size(); i++)
	{
		m_pcs = pcs_vector[i];
		if (hasAnyProcessDeactivatedSubdomains && m_pcs->ExcavMaterialGroup < 0) // NW
			// WX:11.2012 when excavated, not do CheckMarkedElement
			m_pcs->CheckMarkedElement();
#if defined(USE_MPI) // 18.10.2007 WW
		if (myrank == 0)
		{
#endif
			m_pcs->WriteSolution(); // WW
#ifdef GEM_REACT
			if (i == 0) // for GEM_REACT we also need information on porosity (node porosity internally stored in Gems
			            // process)!....do it only once and it does not matter for which process ! ....we assume that
			            // the first pcs process is the flow process...if reload not defined for every process,
			            // restarting with gems will not work in any case

				if ((m_pcs->reload == 1 || m_pcs->reload == 3)
				    && !((aktueller_zeitschritt % m_pcs->nwrite_restart) > 0))
					m_vec_GEM->WriteReloadGem();

#endif
#if defined(USE_MPI) // 18.10.2007 WW
		}
#endif

		m_pcs->Extropolation_MatValue(); // WW
		if (m_pcs->cal_integration_point_value) // WW
			m_pcs->Extropolation_GaussValue();
		// BG
		if ((m_pcs->simulator == "ECLIPSE") || (m_pcs->simulator == "DUMUX"))
		{
			m_pcs->Extropolation_GaussValue();
		}
		// JT: Now done in PreCouplingLoop() // m_pcs->CopyTimestepNODValues(); //MB
		if (force_post_node_copy)
		{ // JT: safety valve. Set this value to true (in Problem()) and values will be copied here.
			m_pcs->CopyTimestepNODValues();
			m_pcs->CopyTimestepELEValues();
		}
	}
// WW
#if !defined(USE_PETSC) && !defined(NEW_EQS) // && defined(other parallel libs)//03~04.3012. WW
	//#ifndef NEW_EQS                                //WW. 07.11.2008
	if (total_processes[1])
		total_processes[1]->AssembleParabolicEquationRHSVector();
#endif
	LOPCalcELEResultants();
}

const GEOLIB::GEOObjects* Problem::getGeoObj() const
{
	return _geo_obj;
}

const std::string& Problem::getGeoObjName() const
{
	return _geo_name;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: LiquidFlow
   Task: Similate liquid flow
   Return: error
   Programming:
   07/2008 WW Extract from LOPTimeLoop_PCS();
   Modification:
   -------------------------------------------------------------------------*/
inline double Problem::LiquidFlow()
{
	int success;
	double error = 0.;
	CRFProcess* m_pcs = total_processes[6];
	if (!m_pcs->selected)
		return error;
	//  error = m_pcs->Execute();
	// Cases: decide, weather to use GEOSYS, ECLIPSE or DuMux; BG 10/2010
	if ((m_pcs->simulator.compare("GEOSYS")
	     == 0)) //         ||(m_pcs->simulator.compare("ECLIPSE")==0)){ // standard: use GeoSys
	{
		error = m_pcs->ExecuteNonLinear(loop_process_number);
#ifdef RESET_4410
		PCSCalcSecondaryVariables(); // PCS member function
#endif
		m_pcs->CalIntegrationPointValue(); // WW
		if (m_pcs->tim_type == TimType::STEADY)
			m_pcs->selected = false;
	}

	if (m_pcs->simulator.compare("ECLIPSE") == 0) // use ECLIPSE to calculate one phase liquid flow, BG
	{
		if (m_pcs->EclipseData == NULL) // SBG if this is the first call, make a new instance
			m_pcs->EclipseData = new CECLIPSEData();
		// call ECLIPSE interface
		success = m_pcs->EclipseData->RunEclipse(m_pcs->Tim->step_current, m_pcs);
		if (success == 0)
			std::cout << "Error running Eclipse!"
			          << "\n";
		if (m_pcs->tim_type == TimType::STEADY)
			m_pcs->selected = false;
	}

	else if (m_pcs->simulator.compare("DUMUX") == 0)
	{
		if (m_pcs->DuMuxData == NULL) // SBG if this is the first call, make a new instance
			m_pcs->DuMuxData = new CDUMUXData();
		// call DUMUX interface
		success = m_pcs->DuMuxData->RunDuMux(m_pcs->Tim->step_current, m_pcs);
		if (success == 0)
			std::cout << "Error running DuMux!"
			          << "\n";
	}

	return error;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: RichardsFlow
   Task: Similate Richards flow
   Return: error
   Programming:
   07/2008 WW Extract from LOPTimeLoop_PCS();
   Modification:
   -------------------------------------------------------------------------*/
inline double Problem::RichardsFlow()
{
	//-------  WW
	double error = 0.;
	CRFProcess* m_pcs = total_processes[2];
	if (!m_pcs->selected)
		return error;
	bool twoflowcpl = false;
	// if(GROUNDWATER_FLOW|| OVERLAND_FLOW) WW
	if (total_processes[1] || total_processes[6])
		twoflowcpl = true;
	if (twoflowcpl)
	{ //-------  WW
		// JOD coupling
		lop_coupling_iterations = m_pcs->m_num->cpl_max_iterations;
		if (pcs_vector.size() > 1 && lop_coupling_iterations > 1)
			m_pcs->CopyCouplingNODValues();
		// WW if(m_pcs->adaption) PCSStorage();
		CFEMesh* m_msh = FEMGet("RICHARDS_FLOW");
		if (m_msh->geo_name.compare("REGIONAL") == 0)
			LOPExecuteRegionalRichardsFlow(m_pcs, loop_process_number);
		else
			error = m_pcs->ExecuteNonLinear(loop_process_number);
		if (m_pcs->saturation_switch == true)
			m_pcs->CalcSaturationRichards(1, false); // JOD
		else
			// WW
			m_pcs->CalcSecondaryVariablesUnsaturatedFlow();
		// WW#ifndef NEW_EQS //WW. 07.11.2008
		// WW      if(lop_coupling_iterations > 1) // JOD  coupling
		// WW         pcs_coupling_error = m_pcs->CalcCouplingNODError();
		// WW#endif
		conducted = true; // WW
	}
	else // WW
	{
		CFEMesh* m_msh = FEMGet("RICHARDS_FLOW"); // WW
		if (m_msh->geo_name.compare("REGIONAL") == 0)
			LOPExecuteRegionalRichardsFlow(m_pcs, loop_process_number);
		else
			error = m_pcs->ExecuteNonLinear(loop_process_number);
		if (m_pcs->TimeStepAccept())
		{
			// WW
			m_pcs->CalcSecondaryVariablesUnsaturatedFlow();
			CalcVelocities = true;
			conducted = true; // WW
		}
	}
	if (m_pcs->TimeStepAccept())
		m_pcs->CalIntegrationPointValue(); // WW
	return error;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: TwoPhaseFlow
   Task: Similate twp-phase flow
   Return: error
   Programming:
   07/2008 WW Extract from LOPTimeLoop_PCS();
   Modification:
   12.2008 WW Update
   01.2008 WW Add phases
   -------------------------------------------------------------------------*/
inline double Problem::TwoPhaseFlow()
{
	double error = 1.0e+8;
	CRFProcess* m_pcs = total_processes[3];
	if (!m_pcs->selected)
		return error; // 12.12.2008 WW
	//
	// 08.01.2009. WW
	for (int i = 0; i < (int)multiphase_processes.size(); i++)
	{
		m_pcs = multiphase_processes[i];
		error = m_pcs->ExecuteNonLinear(loop_process_number);
		if (m_pcs->TimeStepAccept())
		{
			PCSCalcSecondaryVariables();
			m_pcs->CalIntegrationPointValue();
			// CB 12/09 (first time added on 010808) Velocity at CenterOfGravity, required for NAPL dissolution
			if (i == 0) // is 0 in all cases the correct index?
				m_pcs->CalcELEVelocities();
		}
	}
	return error;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: MultiPhaseFlow()
   Task: Similate multi-phase flow by p-p scheme
   Return: error
   Programming:
   07/2008 WW Extract from LOPTimeLoop_PCS();
   Modification:
   12.2008 WW Update
   -------------------------------------------------------------------------*/
inline double Problem::MultiPhaseFlow()
{
	double error = 1.0e+8;
	int success = 0; // BG
	CRFProcess* m_pcs = total_processes[4];
	if (!m_pcs->selected)
		return error; // 12.12.2008 WW
	// initialize density and viscosity if the CO2 phase transition is used
	if (m_pcs->Phase_Transition_Model == 1)
	{
		if (m_pcs->Tim->step_current == 1)
		{
			std::cout << " The Viscosity is not calculated yet!!!"
			          << "\n";
			m_pcs->CalculateFluidDensitiesAndViscositiesAtNodes(m_pcs);
		}
		else
			m_pcs->Phase_Transition_CO2(m_pcs, 1);
	}

	// m_pcs->CalculateFluidDensitiesAndViscositiesAtNodes(m_pcs);

	// Cases: decide, weather to use GEOSYS, ECLIPSE or DuMux; BG 10/2010
	if ((m_pcs->simulator.compare("GEOSYS")
	     == 0)) // ||(m_pcs->simulator.compare("ECLIPSE")==0)){ // standard: use GeoSys
	{
		// if((m_pcs->simulator.compare("GEOSYS")==0) ||(m_pcs->simulator.compare("ECLIPSE")==0)){ // standard: use
		// GeoSys
		error = m_pcs->ExecuteNonLinear(loop_process_number);
		if (m_pcs->TimeStepAccept())
			m_pcs->CalIntegrationPointValue(); // WW
	}

	if (m_pcs->simulator.compare("ECLIPSE") == 0) // use ECLIPSE to calculate multi-phase flow, BG
	{
		if (m_pcs->EclipseData == NULL) // SBG if this is the first call, make a new instance
			m_pcs->EclipseData = new CECLIPSEData();
		// redo WTP
		m_pcs->EclipseData->dissolved_co2_pcs_name_ECL
		    = m_pcs->dissolved_co2_pcs_name; // hand over process name for storing dissolved CO2
		m_pcs->EclipseData->dissolved_co2_ingas_pcs_name_ECL = m_pcs->dissolved_co2_ingas_pcs_name;

		// call ECLIPSE interface
		success = m_pcs->EclipseData->RunEclipse(m_pcs->Tim->step_current, m_pcs);
		if (success == 0)
		{
			std::cout << "Error running Eclipse!"
			          << "\n";
			system("Pause");
			exit(0);
		}
	}
	else if (m_pcs->simulator.compare("DUMUX") == 0)
	{
		if (m_pcs->DuMuxData == NULL) // SBG if this is the first call, make a new instance
		{
			m_pcs->DuMuxData = new CDUMUXData();
			m_pcs->DuMuxData->dissolved_co2_pcs_name_DUMUX
			    = m_pcs->dissolved_co2_pcs_name; // hand over process name for storing dissolved CO2
		}
		// call DUMUX interface
		success = m_pcs->DuMuxData->RunDuMux(m_pcs->Tim->step_current, m_pcs);
		if (success == 0)
			std::cout << "Error running DuMux!"
			          << "\n";
	}
	// CO2-Phase_Transition BG, NB
	if ((m_pcs->Phase_Transition_Model == 1) && ((m_pcs->simulator.compare("GEOSYS") == 0)))
	{
		// check if mfp-model for density and viscosity is 18
		if (m_pcs->Tim->step_current == 1)
		{
			CFluidProperties* FluidProp;

			FluidProp = MFPGet("LIQUID");
			if ((FluidProp->density_model != 18) || (FluidProp->viscosity_model != 18))
			{
				std::cout
				    << "If the Phase_Transition_Model is used the density model and the viscosity model should be 18!"
				    << "\n";
				std::cout << "The run is terminated now ..."
				          << "\n";
				system("Pause");
				exit(0);
			}
			FluidProp = MFPGet("GAS");
			if ((FluidProp->density_model != 18) || (FluidProp->viscosity_model != 18))
			{
				std::cout
				    << "If the Phase_Transition_Model is used the density model and the viscosity model should be 18!"
				    << "\n";
				std::cout << "The run is terminated now ..."
				          << "\n";
				system("Pause");
				exit(0);
			}
		}
		if (m_pcs->Phase_Transition_Model == 1)
		{
			m_pcs->CalculateFluidDensitiesAndViscositiesAtNodes(m_pcs);
			m_pcs->Phase_Transition_CO2(m_pcs, 0);
		}
	}

	if (m_pcs->tim_type == TimType::STEADY)
		m_pcs->selected = false;

	// TestOutputEclipse(m_pcs);
	// TestOutputDuMux(m_pcs);

	if (m_pcs->OutputMassOfGasInModel == true) // 05/2012 BG
		OutputMassOfGasInModel(m_pcs);

	return error;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: TestOutputDuMux()
   Task: provides the sum of CO2 in the model domain in output files
   Return: nothing
   Programming:
   02/2011 BG
   -------------------------------------------------------------------------*/
void Problem::TestOutputDuMux(CRFProcess* m_pcs)
{
	// Testoutput amount of co2 in model domain
	CFEMesh* m_msh = fem_msh_vector[0]; // SB: ToDo hart gesetzt
	MeshLib::CElem* m_ele = NULL;
	MeshLib::CNode* m_node = NULL;
	CMediumProperties* m_mat_mp = NULL;
	std::ostringstream temp;
	double mass_CO2_gas, mass_CO2_water, mass_CO2;
	int index;
	double saturation_CO2;
	double saturation_water;
	double node_volume;
	double time;
	std::string tempstring;
	std::vector<std::string> vec_string;
	// int position;
	std::string path;
	double density_CO2;
	double porosity = 0.0;
	int variable_index;
	double concentration_CO2_water;
	int indexConcentration_CO2;
	// CRFProcess *n_pcs = NULL;
	int group;

	path = m_pcs->file_name_base;
	int position = int(path.find_last_of("\\"));
	std::string path_new;
	path_new = path.substr(0, position);
	// position = int(path_new.find_last_of("\\"));
	// path_new = path_new.substr(0,position);
	if (m_pcs->DuMuxData->Windows_System == true)
		tempstring = path_new + "\\Sum_CO2_nodes.csv";
	else
		tempstring = path_new + "Sum_CO2_nodes.csv";

	if (m_pcs->Tim->step_current == 1)
		// Header of the file
		vec_string.push_back("Time, massCO2_gas, massCO2_water, massCO2, porosity");
	else
	{
		// read file and store data
		CReadTextfiles_DuMux* TextFile;
		TextFile = new CReadTextfiles_DuMux;
		TextFile->Read_Text(tempstring);

		for (int i = 0; i < TextFile->NumberOfRows; i++)
			vec_string.push_back(TextFile->Data[i]);
	}

	mass_CO2 = mass_CO2_gas = mass_CO2_water = 0;
	// +1: new timelevel
	indexConcentration_CO2 = pcs_vector[m_pcs->DuMuxData->ProcessIndex_CO2inLiquid]->GetNodeValueIndex(
	                             pcs_vector[m_pcs->DuMuxData->ProcessIndex_CO2inLiquid]->pcs_primary_function_name[0])
	                         + 1;

	for (long i = 0; i < (long)m_msh->nod_vector.size(); i++)
	{
		m_node = m_msh->nod_vector[i]; // get element
		node_volume = 0;
		saturation_CO2 = 0;
		if (mfp_vector[1]->density_model == 18)
		{
			variable_index = m_pcs->GetNodeValueIndex("DENSITY2");
			density_CO2 = m_pcs->GetNodeValue(i, variable_index);
		}
		else
			density_CO2 = mfp_vector[1]->Density();

		for (int j = 0; j < int(m_node->getConnectedElementIDs().size()); j++)
		{
			m_ele = m_msh->ele_vector[m_node->getConnectedElementIDs()[j]];

			// get the phase volume of current element elem
			group = m_ele->GetPatchIndex();
			m_mat_mp = mmp_vector[group];
			// CB Now provides also heterogeneous porosity, model 11
			porosity = m_mat_mp->Porosity(m_ele->GetIndex(), 1);
			node_volume = node_volume + m_ele->GetVolume() / 8 * porosity;
		}

		//+1... new time level
		index = m_pcs->GetNodeValueIndex("SATURATION1") + 1;
		saturation_water = m_pcs->GetNodeValue(i, index);
		// if (saturation_water < 1)
		saturation_CO2 = 1 - saturation_water;
		concentration_CO2_water
		    = pcs_vector[m_pcs->DuMuxData->ProcessIndex_CO2inLiquid]->GetNodeValue(i, indexConcentration_CO2);

		mass_CO2_gas = mass_CO2_gas + node_volume * saturation_CO2 * density_CO2;
		mass_CO2_water = mass_CO2_water
		                 + node_volume * saturation_water * concentration_CO2_water * m_pcs->DuMuxData->Molweight_CO2
		                       * 0.001;
		// cout << " Node: " << i << " saturation: " << saturation_water << " Density_CO2: " << density_CO2 << "
		// node_volume: " << node_volume << "\n";
	}
	mass_CO2 = mass_CO2_gas + mass_CO2_water;
	// calculating time
	time = 0;
	for (int k = 0; k < m_pcs->Tim->step_current; k++)
		time += m_pcs->Tim->time_step_vector[k];
	temp.str("");
	temp.clear();
	temp << time;
	tempstring = temp.str();
	tempstring += ", ";
	temp.str("");
	temp.clear();
	temp.precision(12);
	temp << mass_CO2_gas;
	tempstring += temp.str();
	tempstring += ", ";
	temp.str("");
	temp.clear();
	temp.precision(12);
	temp << mass_CO2_water;
	tempstring += temp.str();
	tempstring += ", ";
	temp.str("");
	temp.clear();
	temp.precision(12);
	temp << mass_CO2;
	tempstring += temp.str();
	tempstring += ", ";
	temp.str("");
	temp.clear();
	temp.precision(12);
	temp << porosity;
	tempstring += temp.str();

	vec_string.push_back(tempstring);

	// within the first timestep create file and write header
	CWriteTextfiles_DuMux* TextFile;
	TextFile = new CWriteTextfiles_DuMux;
	if (m_pcs->DuMuxData->Windows_System == true)
		tempstring = path_new + "\\Sum_CO2_nodes.csv";
	else
		tempstring = path_new + "Sum_CO2_nodes.csv";
	TextFile->Write_Text(tempstring, vec_string);

	// Testoutput amount of co2 in model domain calculated at nodes-DuMux
	// double element_volume;
	Math_Group::vec<MeshLib::CNode*> ele_nodes(8);

	path = m_pcs->file_name_base;
	position = int(path.find_last_of("\\"));
	path_new = path.substr(0, position);
	// position = int(path_new.find_last_of("\\"));
	// path_new = path_new.substr(0,position);
	if (m_pcs->DuMuxData->Windows_System == true)
		tempstring = path_new + "\\Sum_CO2_nodes_DuMux.csv";
	else
		tempstring = path_new + "Sum_CO2_nodes_DuMux.csv";

	vec_string.clear();
	if (m_pcs->Tim->step_current == 1)
		// Header of the file
		vec_string.push_back("Time, massCO2_gas, massCO2_water, massCO2, porosity");
	else
	{
		// read file and store data
		CReadTextfiles_DuMux* TextFile;
		TextFile = new CReadTextfiles_DuMux;
		TextFile->Read_Text(tempstring);

		for (int i = 0; i < TextFile->NumberOfRows; i++)
			vec_string.push_back(TextFile->Data[i]);
	}

	mass_CO2 = mass_CO2_gas = mass_CO2_water = 0;
	// +1: new timelevel
	indexConcentration_CO2 = pcs_vector[m_pcs->DuMuxData->ProcessIndex_CO2inLiquid]->GetNodeValueIndex(
	                             pcs_vector[m_pcs->DuMuxData->ProcessIndex_CO2inLiquid]->pcs_primary_function_name[0])
	                         + 1;

	for (long i = 0; i < (long)m_pcs->DuMuxData->NodeData.size(); i++)
	{
		m_node = m_msh->nod_vector[i]; // get element
		saturation_CO2 = 0;
		node_volume = 0;

		density_CO2 = m_pcs->DuMuxData->NodeData[i]->getPhaseDensity()[1];

		for (int j = 0; j < int(m_node->getConnectedElementIDs().size()); j++)
		{
			m_ele = m_msh->ele_vector[m_node->getConnectedElementIDs()[j]];

			// get the phase volume of current element elem
			group = m_ele->GetPatchIndex();
			m_mat_mp = mmp_vector[group];
			// CB Now provides also heterogeneous porosity, model 11
			porosity = m_mat_mp->Porosity(m_ele->GetIndex(), 1);
			node_volume = node_volume + m_ele->GetVolume() / 8 * porosity;
		}

		saturation_water = m_pcs->DuMuxData->NodeData[i]->getPhaseSaturation()[0];
		saturation_CO2 = 1 - saturation_water;

		// concentration_CO2_water = pcs_vector[m_pcs->DuMuxData->ProcessIndex_CO2inLiquid]->GetNodeValue(i,
		// indexConcentration_CO2);
		concentration_CO2_water = m_pcs->DuMuxData->NodeData[i]->getCO2InLiquid()
		                          * m_pcs->DuMuxData->NodeData[i]->getPhaseDensity()[0]
		                          / (m_pcs->DuMuxData->Molweight_CO2 * 1e-3);

		mass_CO2_gas = mass_CO2_gas + node_volume * saturation_CO2 * density_CO2;
		mass_CO2_water = mass_CO2_water
		                 + node_volume * saturation_water * concentration_CO2_water * m_pcs->DuMuxData->Molweight_CO2
		                       * 0.001;
	}
	mass_CO2 = mass_CO2_gas + mass_CO2_water;
	// calculating time
	time = 0;
	for (int k = 0; k < m_pcs->Tim->step_current; k++)
		time += m_pcs->Tim->time_step_vector[k];
	temp.str("");
	temp.clear();
	temp << time;
	tempstring = temp.str();
	tempstring += ", ";
	temp.str("");
	temp.clear();
	temp.precision(12);
	temp << mass_CO2_gas;
	tempstring += temp.str();
	tempstring += ", ";
	temp.str("");
	temp.clear();
	temp.precision(12);
	temp << mass_CO2_water;
	tempstring += temp.str();
	tempstring += ", ";
	temp.str("");
	temp.clear();
	temp.precision(12);
	temp << mass_CO2;
	tempstring += temp.str();
	tempstring += ", ";
	temp.str("");
	temp.clear();
	temp.precision(12);
	temp << porosity;
	tempstring += temp.str();

	vec_string.push_back(tempstring);

	// within the first timestep create file and write header
	TextFile = new CWriteTextfiles_DuMux;
	if (m_pcs->DuMuxData->Windows_System == true)
		tempstring = path_new + "\\Sum_CO2_nodes_DuMux.csv";
	else
		tempstring = path_new + "Sum_CO2_nodes_DuMux.csv";
	TextFile->Write_Text(tempstring, vec_string);
}

/*-------------------------------------------------------------------------
   GeoSys - Function: TestOutputEclipse()
   Task: provides the sum of CO2 in the model domain in output files
   Return: nothing
   Programming:
   02/2011 BG
   -------------------------------------------------------------------------*/
void Problem::TestOutputEclipse(CRFProcess* m_pcs)
{
	// Testoutput amount of co2 in model domain calculated at nodes
	CFEMesh* m_msh = fem_msh_vector[0]; // SB: ToDo hart gesetzt
	MeshLib::CElem* m_ele = NULL;
	MeshLib::CNode* m_node = NULL;
	CMediumProperties* m_mat_mp = NULL;
	std::ostringstream temp;
	double mass_CO2_gas, mass_CO2_water, mass_CO2;
	int index;
	double saturation_CO2;
	double saturation_water;
	double node_volume;
	double time;
	std::string tempstring;
	std::vector<std::string> vec_string;
	// int position;
	std::string path;
	double density_CO2;
	double porosity = 0.0;
	int variable_index;
	double concentration_CO2_water;
	int indexConcentration_CO2;
	// CRFProcess *n_pcs = NULL;
	int group;

	path = m_pcs->file_name_base;
	int position = int(path.find_last_of("\\"));
	std::string path_new;
	path_new = path.substr(0, position);
	// position = int(path_new.find_last_of("\\"));
	// path_new = path_new.substr(0,position);
	tempstring = path_new + "\\Sum_CO2_nodes.csv";

	if (m_pcs->Tim->step_current == 1)
		// Header of the file
		vec_string.push_back("Time, massCO2_gas, massCO2_water, massCO2, porosity");
	else
	{
		// read file and store data
		CReadTextfiles_ECL* TextFile;
		TextFile = new CReadTextfiles_ECL;
		TextFile->Read_Text(tempstring);

		for (int i = 0; i < TextFile->NumberOfRows; i++)
			vec_string.push_back(TextFile->Data[i]);
	}

	mass_CO2 = mass_CO2_gas = mass_CO2_water = 0;
	// +1: new timelevel
	indexConcentration_CO2 = pcs_vector[m_pcs->EclipseData->ProcessIndex_CO2inLiquid]->GetNodeValueIndex(
	                             pcs_vector[m_pcs->EclipseData->ProcessIndex_CO2inLiquid]->pcs_primary_function_name[0])
	                         + 1;

	for (long i = 0; i < (long)m_msh->nod_vector.size(); i++)
	{
		m_node = m_msh->nod_vector[i]; // get element
		node_volume = 0;
		saturation_CO2 = 0;
		if (mfp_vector[1]->density_model == 18)
		{
			variable_index = m_pcs->GetNodeValueIndex("DENSITY2");
			density_CO2 = m_pcs->GetNodeValue(i, variable_index);
		}
		else
			density_CO2 = mfp_vector[1]->Density();

		for (int j = 0; j < int(m_node->getConnectedElementIDs().size()); j++)
		{
			m_ele = m_msh->ele_vector[m_node->getConnectedElementIDs()[j]];

			// get the phase volume of current element elem
			group = m_ele->GetPatchIndex();
			m_mat_mp = mmp_vector[group];
			// CB Now provides also heterogeneous porosity, model 11
			porosity = m_mat_mp->Porosity(m_ele->GetIndex(), 1);
			node_volume = node_volume + m_ele->GetVolume() / 8 * porosity;
		}

		//+1... new time level
		index = m_pcs->GetNodeValueIndex("SATURATION1") + 1;
		saturation_water = m_pcs->GetNodeValue(i, index);
		// if (saturation_water < 1)
		saturation_CO2 = 1 - saturation_water;
		concentration_CO2_water
		    = pcs_vector[m_pcs->EclipseData->ProcessIndex_CO2inLiquid]->GetNodeValue(i, indexConcentration_CO2);

		mass_CO2_gas = mass_CO2_gas + node_volume * saturation_CO2 * density_CO2;
		mass_CO2_water = mass_CO2_water
		                 + node_volume * saturation_water * concentration_CO2_water * m_pcs->EclipseData->Molweight_CO2
		                       * 0.001;
		// cout << " Node: " << i << " saturation: " << saturation_water << " Density_CO2: " << density_CO2 << "
		// node_volume: " << node_volume << "\n";
	}
	mass_CO2 = mass_CO2_gas + mass_CO2_water;
	// calculating time
	time = 0;
	for (int k = 0; k < m_pcs->Tim->step_current; k++)
		time += m_pcs->Tim->time_step_vector[k];
	temp.str("");
	temp.clear();
	temp << time;
	tempstring = temp.str();
	tempstring += ", ";
	temp.str("");
	temp.clear();
	temp.precision(12);
	temp << mass_CO2_gas;
	tempstring += temp.str();
	tempstring += ", ";
	temp.str("");
	temp.clear();
	temp.precision(12);
	temp << mass_CO2_water;
	tempstring += temp.str();
	tempstring += ", ";
	temp.str("");
	temp.clear();
	temp.precision(12);
	temp << mass_CO2;
	tempstring += temp.str();
	tempstring += ", ";
	temp.str("");
	temp.clear();
	temp.precision(12);
	temp << porosity;
	tempstring += temp.str();

	vec_string.push_back(tempstring);

	// within the first timestep create file and write header
	CWriteTextfiles_ECL* TextFile;
	TextFile = new CWriteTextfiles_ECL;
	tempstring = path_new + "\\Sum_CO2_nodes.csv";
	TextFile->Write_Text(tempstring, vec_string);

	// Testoutput amount of co2 in model domain calculated at elements
	double element_volume;
	Math_Group::vec<MeshLib::CNode*> ele_nodes(8);

	path = m_pcs->file_name_base;
	position = int(path.find_last_of("\\"));
	path_new = path.substr(0, position);
	// position = int(path_new.find_last_of("\\"));
	// path_new = path_new.substr(0,position);
	tempstring = path_new + "\\Sum_CO2_elements.csv";

	vec_string.clear();
	if (m_pcs->Tim->step_current == 1)
		// Header of the file
		vec_string.push_back("Time, massCO2_gas, massCO2_water, massCO2, porosity");
	else
	{
		// read file and store data
		CReadTextfiles_ECL* TextFile;
		TextFile = new CReadTextfiles_ECL;
		TextFile->Read_Text(tempstring);

		for (int i = 0; i < TextFile->NumberOfRows; i++)
			vec_string.push_back(TextFile->Data[i]);
	}

	mass_CO2 = mass_CO2_gas = mass_CO2_water = 0;
	// +1: new timelevel
	indexConcentration_CO2 = pcs_vector[m_pcs->EclipseData->ProcessIndex_CO2inLiquid]->GetNodeValueIndex(
	                             pcs_vector[m_pcs->EclipseData->ProcessIndex_CO2inLiquid]->pcs_primary_function_name[0])
	                         + 1;

	for (long i = 0; i < long(m_pcs->EclipseData->eclgrid.size()); i++)
	{
		m_ele = m_msh->ele_vector[i];
		m_node = m_msh->nod_vector[i]; // get element
		saturation_CO2 = 0;
		element_volume = 0;

		if (m_pcs->EclipseData->E100 == true)
			density_CO2 = m_pcs->EclipseData->Data[i][m_pcs->EclipseData->GetVariableIndex("GAS_DEN")];
		else
			density_CO2 = m_pcs->EclipseData->Data[i][m_pcs->EclipseData->GetVariableIndex("DENG")];

		group = m_ele->GetPatchIndex();
		m_mat_mp = mmp_vector[group];
		// CB Now provides also heterogeneous porosity, model 11
		porosity = m_mat_mp->Porosity(m_ele->GetIndex(), 1);
		element_volume = m_ele->GetVolume() * porosity;
		if (m_pcs->EclipseData->E100 == true)
			saturation_water = 1 - m_pcs->EclipseData->Data[i][m_pcs->EclipseData->GetVariableIndex("SGAS")];
		else
			saturation_water = m_pcs->EclipseData->Data[i][m_pcs->EclipseData->GetVariableIndex("SWAT")];
		saturation_CO2 = m_pcs->EclipseData->Data[i][m_pcs->EclipseData->GetVariableIndex("SGAS")];

		m_ele->GetNodes(ele_nodes);
		concentration_CO2_water = 0;
		for (int j = 0; j < int(ele_nodes.Size()); j++)
			concentration_CO2_water = concentration_CO2_water
			                          + pcs_vector[m_pcs->EclipseData->ProcessIndex_CO2inLiquid]->GetNodeValue(
			                                ele_nodes[j]->GetIndex(), indexConcentration_CO2) / ele_nodes.Size();

		mass_CO2_gas = mass_CO2_gas + element_volume * saturation_CO2 * density_CO2;
		mass_CO2_water = mass_CO2_water
		                 + element_volume * saturation_water * concentration_CO2_water
		                       * m_pcs->EclipseData->Molweight_CO2 * 0.001;
	}
	mass_CO2 = mass_CO2_gas + mass_CO2_water;
	// calculating time
	time = 0;
	for (int k = 0; k < m_pcs->Tim->step_current; k++)
		time += m_pcs->Tim->time_step_vector[k];
	temp.str("");
	temp.clear();
	temp << time;
	tempstring = temp.str();
	tempstring += ", ";
	temp.str("");
	temp.clear();
	temp.precision(12);
	temp << mass_CO2_gas;
	tempstring += temp.str();
	tempstring += ", ";
	temp.str("");
	temp.clear();
	temp.precision(12);
	temp << mass_CO2_water;
	tempstring += temp.str();
	tempstring += ", ";
	temp.str("");
	temp.clear();
	temp.precision(12);
	temp << mass_CO2;
	tempstring += temp.str();
	tempstring += ", ";
	temp.str("");
	temp.clear();
	temp.precision(12);
	temp << porosity;
	tempstring += temp.str();

	vec_string.push_back(tempstring);

	// within the first timestep create file and write header
	TextFile = new CWriteTextfiles_ECL;
	tempstring = path_new + "\\Sum_CO2_elements.csv";
	TextFile->Write_Text(tempstring, vec_string);
}

/*-------------------------------------------------------------------------
   GeoSys - Function: OutputMassOfGasInModel()
   Task: provides the sum of CO2 in the model domain in output files
   Return: nothing
   Programming:
   02/2011 BG
   -------------------------------------------------------------------------*/
void Problem::OutputMassOfGasInModel(CRFProcess* m_pcs)
{
	// Testoutput amount of co2 in model domain calculated at nodes
	CFEMesh* m_msh = fem_msh_vector[0]; // SB: ToDo hart gesetzt
	MeshLib::CElem* m_ele = NULL;
	MeshLib::CNode* m_node = NULL;
	CMediumProperties* m_mat_mp = NULL;
	std::ostringstream temp;
	double mass_Gas_gas, mass_Gas_water, mass_Gas;
	int index;
	double saturation_Gas;
	double saturation_water;
	double node_volume;
	double time;
	std::string tempstring;
	std::vector<std::string> vec_string;
	// int position;
	std::string path;
	double density_Gas;
	double porosity = 0.0;
	int variable_index;
	double concentration_Gas_water;
	int indexConcentration_Gas = 0;
	int ProcessIndex_GasInLiquid;
	CRFProcess* n_pcs = NULL;
	int group;
	std::string transport_process_name, Processname;
	double Molweight_Gas;
	double V_model;
	// bool Windows_System;
	int position;
	std::string filename;

	path = m_pcs->file_name_base;

#ifdef _WIN32
	{
		position = int(path.find_last_of("\\"));
		if (position < 0)
			filename = "Sum_Gas_nodes.csv";
		else
		{
			std::string path_new;
			path_new = path.substr(0, position);
			filename = path_new + "\\Sum_Gas_nodes.csv";
		}
	}
#else
	{
		position = int(path.find_last_of("/"));
		if (position < 0)
			filename = "Sum_Gas_nodes.csv";
		else
		{
			std::string path_new;
			path_new = path.substr(0, position);
			filename = path_new + "/Sum_Gas_nodes.csv";
		}
	}
#endif

	transport_process_name = "H2";
	Molweight_Gas = 2;
	Processname = "PS_GLOBAL";

	ProcessIndex_GasInLiquid = -1;
	if (ProcessIndex_GasInLiquid == -1)
		for (int i = 0; i < int(pcs_vector.size()); i++)
		{
			n_pcs = pcs_vector[i];
			// identify your process and store idx of pcs-vector
			if (n_pcs->nod_val_name_vector[0] == transport_process_name)
				ProcessIndex_GasInLiquid = i;
		}

	if (m_pcs->Tim->step_current == 1)
		// Header of the file
		vec_string.push_back("Time, massGas_gas, massGas_water, massGas, porosity, Modelvolume, Pressure2, Density2");
	else
	{
		CReadTextfiles_ECL TextFile;
		TextFile.Read_Text(filename);

		for (int i = 0; i < TextFile.NumberOfRows; i++)
			vec_string.push_back(TextFile.Data[i]);
	}

	V_model = mass_Gas = mass_Gas_gas = mass_Gas_water = 0;
	// +1: new timelevel
	if (ProcessIndex_GasInLiquid > -1)
		indexConcentration_Gas = pcs_vector[ProcessIndex_GasInLiquid]->GetNodeValueIndex(
		                             pcs_vector[ProcessIndex_GasInLiquid]->pcs_primary_function_name[0]) + 1;

	for (long i = 0; i < (long)m_msh->nod_vector.size(); i++)
	{
		m_node = m_msh->nod_vector[i]; // get element
		node_volume = 0;
		saturation_Gas = 0;
		if (mfp_vector[1]->density_model == 18)
		{
			variable_index = m_pcs->GetNodeValueIndex("DENSITY2");
			density_Gas = m_pcs->GetNodeValue(i, variable_index);
		}
		else
		{
			variable_index = m_pcs->GetNodeValueIndex("PRESSURE2");
			double p2 = m_pcs->GetNodeValue(i, variable_index);
			double dens_arg[3];
			dens_arg[0] = p2;
			density_Gas = mfp_vector[1]->Density(dens_arg);
			// cout << "Knoten: " << i << " Dichte: " << density_Gas << " P2: " << p2 << "\n";
		}

		for (int j = 0; j < int(m_node->getConnectedElementIDs().size()); j++)
		{
			m_ele = m_msh->ele_vector[m_node->getConnectedElementIDs()[j]];

			// get the phase volume of current element elem
			group = m_ele->GetPatchIndex();
			m_mat_mp = mmp_vector[group];
			// CB Now provides also heterogeneous porosity, model 11
			porosity = m_mat_mp->Porosity(m_ele->GetIndex(), 1);
			// node_volume = node_volume + m_ele->GetVolume() / 8 * porosity;
			if (m_ele->GetElementType() == 2)
				node_volume = node_volume + m_ele->GetVolume() / m_ele->GetNodesNumber(true) * porosity;
			else
				node_volume = node_volume + m_ele->GetVolume() / m_ele->GetNodesNumber(false) * porosity;
		}

		V_model = V_model + node_volume;

		//+1... new time level
		if (Processname == "PS_GLOBAL")
			index = m_pcs->GetNodeValueIndex("SATURATION1");
		else
			index = m_pcs->GetNodeValueIndex("SATURATION1") + 1;
		saturation_water = m_pcs->GetNodeValue(i, index);
		// if (saturation_water < 1)
		saturation_Gas = 1 - saturation_water;
		if (ProcessIndex_GasInLiquid > -1)
		{
			concentration_Gas_water = pcs_vector[ProcessIndex_GasInLiquid]->GetNodeValue(i, indexConcentration_Gas);
			mass_Gas_water = mass_Gas_water
			                 + node_volume * saturation_water * concentration_Gas_water * Molweight_Gas * 0.001;
		}
		else
			mass_Gas_water = 0;

		mass_Gas_gas = mass_Gas_gas + node_volume * saturation_Gas * density_Gas;
		// cout << " Node: " << i << " saturation: " << saturation_water << " Density_CO2: " << density_CO2 << "
		// node_volume: " << node_volume << "\n";
	}
	mass_Gas = mass_Gas_gas + mass_Gas_water;
	// calculating time
	time = 0;
	for (int k = 0; k < m_pcs->Tim->step_current; k++)
		time += m_pcs->Tim->time_step_vector[k];
	temp.str("");
	temp.clear();
	temp << time;
	tempstring = temp.str();
	tempstring += ", ";
	temp.str("");
	temp.clear();
	temp.precision(12);
	temp << mass_Gas_gas;
	tempstring += temp.str();
	tempstring += ", ";
	temp.str("");
	temp.clear();
	temp.precision(12);
	temp << mass_Gas_water;
	tempstring += temp.str();
	tempstring += ", ";
	temp.str("");
	temp.clear();
	temp.precision(12);
	temp << mass_Gas;
	tempstring += temp.str();
	tempstring += ", ";
	temp.str("");
	temp.clear();
	temp.precision(12);
	temp << porosity;
	tempstring += temp.str();
	tempstring += ", ";
	temp.str("");
	temp.clear();
	temp.precision(12);
	temp << V_model;
	tempstring += temp.str();

	// Test output pressure and density
	variable_index = m_pcs->GetNodeValueIndex("PRESSURE2");
	double p2 = m_pcs->GetNodeValue(1, variable_index);
	double dens_arg[3];
	dens_arg[0] = p2;
	density_Gas = mfp_vector[1]->Density(dens_arg);
	tempstring += ", ";
	temp.str("");
	temp.clear();
	temp.precision(12);
	temp << p2;
	tempstring += temp.str();
	tempstring += ", ";
	temp.str("");
	temp.clear();
	temp.precision(12);
	temp << density_Gas;
	tempstring += temp.str();
	// Test output pressure and density

	vec_string.push_back(tempstring);

	// within the first timestep create file and write header
	CWriteTextfiles_ECL TextFile;
	TextFile.Write_Text(filename, vec_string);
}

/*-------------------------------------------------------------------------
GeoSys - Function: OutputMassOfComponentInModel()
Task: provides the sum of each transport component in the model domain in output files
Return: nothing
Programming:
05/2011 BG
-------------------------------------------------------------------------*/
void Problem::OutputMassOfComponentInModel(std::vector<CRFProcess*> flow_pcs, CRFProcess* transport_pcs)
{
	// Testoutput mass of component in model domain calculated at nodes
	MeshLib::CFEMesh* m_msh = fem_msh_vector[0]; // SB: ToDo hart gesetzt
	MeshLib::CElem* m_ele = NULL;
	MeshLib::CNode* m_node = NULL;
	CMediumProperties* m_mat_mp = NULL;
	std::ostringstream temp;
	double ComponentMass;
	// //int index;
	double saturation_water;
	double node_volume;
	double time;
	std::string tempstring;
	std::vector<std::string> vec_string;
	// int position;
	std::string path;
	//	double density_water; // 2012-08 TF, variable ‘density_water’ set but not used [-Wunused-but-set-variable]
	double porosity;
	//	int variable_index; // 2012-08 TF, variable set but not used
	double ComponentConcentration, TotalVolume, SourceTerm;
	int indexComponentConcentration;
	// CRFProcess *n_pcs = NULL;
	int group;
	std::string filename;
	// bool Windows_System;
	int position;

	path = flow_pcs[0]->file_name_base;
// position = int(path.find_last_of("\\"));
// if (position >= 0)
//	Windows_System = true;
// else
//	Windows_System = false;
// if (Windows_System == true)
#ifdef _WIN32
	{
		position = int(path.find_last_of("\\"));
		if (position < 0)
			filename = "Mass_" + transport_pcs->nod_val_name_vector[0] + "_nodes.csv";
		else
		{
			std::string path_new;
			path_new = path.substr(0, position);
			filename = path_new + "\\Mass_" + transport_pcs->nod_val_name_vector[0] + "_nodes.csv";
		}
	}
#else
	{
		position = int(path.find_last_of("/"));
		if (position < 0)
			filename = "Mass_" + transport_pcs->nod_val_name_vector[0] + "_nodes.csv";
		else
		{
			std::string path_new;
			path_new = path.substr(0, position);
			filename = path_new + "/Mass_" + transport_pcs->nod_val_name_vector[0] + "_nodes.csv";
		}
	}
#endif
	if (flow_pcs[0]->Tim->step_current == 1)
	{
		// Header of the file
		vec_string.push_back("Time, ComponentMass, TotalVolume");
	}
	else
	{
		// read file and store data
		CReadTextfiles_ECL TextFile;
		TextFile.Read_Text(filename);

		for (int i = 0; i < TextFile.NumberOfRows; i++)
		{
			vec_string.push_back(TextFile.Data[i]);
		}
	}

	ComponentMass = TotalVolume = SourceTerm = 0;

	for (long i = 0; i < (long)m_msh->nod_vector.size(); i++)
	{
		m_node = m_msh->nod_vector[i]; // get element
		node_volume = 0;
		if (mfp_vector[0]->density_model == 18)
		{
			//			variable_index = flow_pcs[0]->GetNodeValueIndex("DENSITY1"); // // 2012-08 TF, variable set but not
			//used
			//			density_water = flow_pcs[0]->GetNodeValue(i, variable_index); // 2012-08 TF, variable
			//‘density_water’ set but not used [-Wunused-but-set-variable]
		}
		//		else
		//			density_water = mfp_vector[0]->Density(); // 2012-08 TF, variable ‘density_water’ set but not used
		//[-Wunused-but-set-variable]

		for (int j = 0; j < int(m_node->getConnectedElementIDs().size()); j++)
		{
			m_ele = m_msh->ele_vector[m_node->getConnectedElementIDs()[j]];

			// get the phase volume of current element elem
			group = m_ele->GetPatchIndex();
			m_mat_mp = mmp_vector[group];
			porosity
			    = m_mat_mp->Porosity(m_ele->GetIndex(), 1); // CB Now provides also heterogeneous porosity, model 11
			node_volume = node_volume
			              + porosity * m_ele->GetVolume() * m_ele->GetFluxArea() / m_ele->GetNodesNumber(false);
			// cout << m_ele->GetNodesNumber(false) << " " << "\n";
		}

		TotalVolume += node_volume;
		if (flow_pcs[0]->getProcessType() == FiniteElement::GROUNDWATER_FLOW)
			saturation_water = 1; // saturation does not exist for example in groundwater flow
		else
		{
			int index = flow_pcs[0]->GetNodeValueIndex("SATURATION1") + 1; //+1... new time level
			saturation_water = flow_pcs[0]->GetNodeValue(i, index);
		}
		indexComponentConcentration = transport_pcs->GetNodeValueIndex(transport_pcs->pcs_primary_function_name[0])
		                              + 1; // +1: new timelevel
		ComponentConcentration = transport_pcs->GetNodeValue(i, indexComponentConcentration);

		ComponentMass = ComponentMass + node_volume * saturation_water * ComponentConcentration;
		// cout << " Node: " << i << " saturation: " << saturation_water << " Density_water: " << density_water << "
		// node_volume: " << node_volume  << " Concentration: " << ComponentConcentration << "\n";
	}

	// calculating source term
	for (long i = 0; i < (long)transport_pcs->st_node_value.size(); i++)
	{
		SourceTerm = transport_pcs->st_node_value[i]->node_value;
	}
	SourceTerm *= transport_pcs->Tim->this_stepsize;

	// calculating time
	time = 0;
	for (int k = 0; k < transport_pcs->Tim->step_current; k++)
	{
		time += transport_pcs->Tim->time_step_vector[k];
	}
	temp.str("");
	temp.clear();
	temp.precision(14);
	temp << time;
	tempstring = temp.str();
	tempstring += ", ";
	temp.str("");
	temp.clear();
	temp.precision(14);
	temp << ComponentMass;
	tempstring += temp.str();
	tempstring += ", ";
	temp.str("");
	temp.clear();
	temp.precision(14);
	temp << TotalVolume;
	tempstring += temp.str();

	vec_string.push_back(tempstring);

	// within the first timestep create file and write header
	CWriteTextfiles_ECL TextFile;
	TextFile.Write_Text(filename, vec_string);
}

/*-------------------------------------------------------------------------
   GeoSys - Function: PS_Global()
   Task: Similate multi-phase flow by p-p scheme
   Return: error
   Programming:
   03/2009 PCH Implementation
   Modification:
   -------------------------------------------------------------------------*/
inline double Problem::PS_Global()
{
	double error = 1.0e+8;
	CRFProcess* m_pcs = total_processes[3];
	if (!m_pcs->selected)
		return error;
	error = m_pcs->ExecuteNonLinear(loop_process_number);
	if (m_pcs->TimeStepAccept())
		m_pcs->CalIntegrationPointValue();

	if (m_pcs->OutputMassOfGasInModel == true) // 05/2012 BG
		OutputMassOfGasInModel(m_pcs);

	return error;
}
/*-------------------------------------------------------------------------
   GeoSys - Function: MULTI_COMPONENTIAL_FLOW()
   Task: Multi-componential flow with global approach
   Return: error
   Programming:
   02/2011 AKS/NB Implementation
   Modification:
   -------------------------------------------------------------------------*/
inline double Problem::MULTI_COMPONENTIAL_FLOW()
{
	double error = 1.0e+8;
	CRFProcess* m_pcs = total_processes[5];
	if (!m_pcs->selected)
		return error;
	error = m_pcs->ExecuteNonLinear(loop_process_number);
	if (m_pcs->TimeStepAccept())
		m_pcs->CalIntegrationPointValue();
	return error;
}
/*-------------------------------------------------------------------------
GeoSys - Function: TNEQ()
Task: Simulate Reactive thermal nonequilibrium flow
Return: error
Programming:
07/2013 HS,TN Implementation
Modification:
-------------------------------------------------------------------------*/
inline double Problem::TNEQ()
{
	double error = 1.0e+8;
	CRFProcess* m_pcs = total_processes[14];
	if (!m_pcs->selected)
		return error;
	error = m_pcs->ExecuteNonLinear(loop_process_number);
	if (m_pcs->TimeStepAccept())
		m_pcs->CalIntegrationPointValue();
	return error;
}

/*-------------------------------------------------------------------------
GeoSys - Function: TES()
Task: Simulate Reactive thermal nonequilibrium flow
Return: error
Programming:
07/2013 HS,TN Implementation
Modification:
-------------------------------------------------------------------------*/
inline double Problem::TES()
{
	double error = 1.0e+8;
	CRFProcess* m_pcs = total_processes[15];
	if (!m_pcs->selected)
		return error;
	error = m_pcs->ExecuteNonLinear(loop_process_number);
	if (m_pcs->TimeStepAccept())
		m_pcs->CalIntegrationPointValue();
	return error;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: GroundWaterFlow()
   Task:
   Return: error
   Programming:
   07/2008 WW Extract from LOPTimeLoop_PCS();
   Modification:
   12.2008 WW Update
   05.2009 WW For surface-soil-ground coupled model
   -------------------------------------------------------------------------*/
inline double Problem::GroundWaterFlow()
{
	double error = 1.0e+8;
	CRFProcess* m_pcs = total_processes[1];
	if (!m_pcs->selected)
		return error; // 12.12.2008 WW
	ClockTimeVec[0]->StartTime(); // CB time

	//----- For the coupling with the soil column approach. 05.2009. WW
	MeshLib::GridsTopo* neighb_grid = NULL;
	MeshLib::GridsTopo* neighb_grid_this = NULL;
	CRFProcess* neighb_pcs = total_processes[2];
	std::vector<double> border_flux;
	int idx_flux = 0, idx_flux_this;
	// WW int no_local_nodes;

	long i;
	if (neighb_pcs)
	{
		for (i = 0; i < (long)neighb_pcs->m_msh->grid_neighbors.size(); i++)
		{
			neighb_grid = neighb_pcs->m_msh->grid_neighbors[i];
			if (neighb_grid->getNeighbor_Name().find("SECTOR_GROUND") != std::string::npos)
				break;
			neighb_grid = NULL;
		}
		for (i = 0; i < (long)m_pcs->m_msh->grid_neighbors.size(); i++)
		{
			neighb_grid_this = m_pcs->m_msh->grid_neighbors[i];
			if (neighb_grid_this->getNeighbor_Name().find("SECTOR_SOIL") != std::string::npos)
				break;
			neighb_grid_this = NULL;
		}
	}
	//------------
	if (neighb_grid)
	{
		CSourceTerm* m_st = NULL;
		CNodeValue* m_nod_val = NULL;
		long l = 0, m = 0;
		// WW no_local_nodes = neighb_pcs->m_msh->getNumberOfMeshLayers()+1;

		if (aktueller_zeitschritt == 1)
		{
			m_st = new CSourceTerm();
			for (i = 0; i < neighb_grid->getBorderNodeNumber(); i++)
			{
				m_nod_val = new CNodeValue();
				m_pcs->st_node_value.push_back(m_nod_val);
				m_pcs->st_node.push_back(m_st);
			}
		}
		l = (long)m_pcs->st_node_value.size();
		long* local_indxs = neighb_grid->getBorderNodeIndicies();
		long* local_indxs_this = neighb_grid_this->getBorderNodeIndicies();
		border_flux.resize(neighb_grid->getBorderNodeNumber());
		idx_flux = neighb_pcs->GetNodeValueIndex("VELOCITY_Z1");

		for (i = 0; i < neighb_grid->getBorderNodeNumber(); i++)
			border_flux[local_indxs_this[i]] = neighb_pcs->GetNodeValue(local_indxs[i], idx_flux)
			                                   / neighb_pcs->time_unit_factor;

		m_pcs->Integration(border_flux);

		m = l - neighb_grid->getBorderNodeNumber();
		for (i = m; i < l; i++)
		{
			m_nod_val = m_pcs->st_node_value[i];
			m_nod_val->msh_node_number = local_indxs_this[i - m];
			m_nod_val->geo_node_number = local_indxs_this[i - m];
			m_nod_val->node_value = -border_flux[local_indxs_this[i - m]];
		}
	}
	//-------------------------------

	error = m_pcs->ExecuteNonLinear(loop_process_number);
	//................................................................
	// Calculate secondary variables
	// NOD values
	conducted = true; // WW
	// std::cout << "      Calculation of secondary NOD values" << "\n";
	if (m_pcs->TimeStepAccept())
	{
#ifdef RESET_4410
		PCSCalcSecondaryVariables(); // PCS member function
#endif
		// std::cout << "      Calculation of secondary GP values" << "\n";
		m_pcs->CalIntegrationPointValue(); // WW
		m_pcs->cal_integration_point_value = true; // WW Do not extropolate Gauss velocity

		if (neighb_grid)
		{
			m_pcs->Extropolation_GaussValue();
			m_pcs->cal_integration_point_value = false;
			idx_flux_this = m_pcs->GetNodeValueIndex("VELOCITY_Z1");
			for (i = 0; i < neighb_grid->getBorderNodeNumber(); i++)
				border_flux[i] = m_pcs->GetNodeValue(i, idx_flux_this);
			m_pcs->Integration(border_flux);
			for (i = 0; i < neighb_grid->getBorderNodeNumber(); i++)
				neighb_pcs->SetNodeValue(i, idx_flux, border_flux[i]);
			//         neighb_pcs->SetNodeValue(i+no_local_nodes, idx_flux, border_flux[i]);
		}
	}
// ELE values
#if !defined(USE_PETSC) && !defined(NEW_EQS) // && defined(other parallel libs)//03~04.3012. WW
	//#ifndef NEW_EQS                                //WW. 07.11.2008
	if (m_pcs->tim_type == TimType::STEADY) // CMCD 05/2006
	{
		// std::cout << "      Calculation of secondary ELE values" << "\n";
		m_pcs->AssembleParabolicEquationRHSVector(); // WW LOPCalcNODResultants();
		m_pcs->CalcELEVelocities();
		m_pcs->selected = false;
	}
#endif

	ClockTimeVec[0]->StopTime("Flow"); // CB time
	ClockTimeVec[0]->StartTime(); // CB time

	return error;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: ComponentalFlow();
   Task:
   Return: error
   Programming:
   07/2008 WW Extract from LOPTimeLoop_PCS();
   Modification:
   -------------------------------------------------------------------------*/
inline double Problem::ComponentalFlow()
{
	double error = 1.e8;
	CRFProcess* m_pcs = total_processes[5];
	if (!m_pcs->selected)
		return error; // 12.12.2008 WW
	//
	error = m_pcs->ExecuteNonLinear(loop_process_number);
	m_pcs->CalIntegrationPointValue(); // WW
	return error;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: OverlandFlow()
   Task:
   Return: error
   Programming:
   07/2008 WW Extract from LOPTimeLoop_PCS();
   Modification:
   -------------------------------------------------------------------------*/
inline double Problem::OverlandFlow()
{
	double error = 1.e8;
	CRFProcess* m_pcs = total_processes[0];
	if (!m_pcs->selected)
		return error; // 12.12.2008 WW

	error = m_pcs->ExecuteNonLinear(loop_process_number);
	if (m_pcs->TimeStepAccept())
		PCSCalcSecondaryVariables();
	return error;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: AirFlow()
   Task:
   Return: error
   Programming:
   07/2008 WW Extract from LOPTimeLoop_PCS();
   Modification:
   -------------------------------------------------------------------------*/
inline double Problem::AirFlow()
{
	double error = 1.e8;
	CRFProcess* m_pcs = total_processes[7];
	if (!m_pcs->selected)
		return error; // 12.12.2008 WW

	error = m_pcs->ExecuteNonLinear(loop_process_number);
	m_pcs->CalIntegrationPointValue(); // WW
	m_pcs->cal_integration_point_value = false; // AKS
	m_pcs->CalcELEVelocities(); // OK
	//
	return error;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: HeatTransport
   Task: Similate heat transport
   Return: error
   Programming:
   07/2008 WW Extract from LOPTimeLoop_PCS();
   Modification:
   -------------------------------------------------------------------------*/
inline double Problem::HeatTransport()
{
	double error = 1.0e+8;
	CRFProcess* m_pcs = total_processes[8];
	if (!m_pcs->selected)
		return error; // 12.12.2008 WW
	// CB This is a cheat to map a 2D horizontal heat pump distribution on a vertical model
	if (REACTINT_vec.size() > 0)
		if (REACTINT_vec[0]->heatpump_2DhTO2Dv)
			REACTINT_vec[0]->Heatpump_2DhTO2Dv_Mapping(false); // Restore the old solution

	error = m_pcs->ExecuteNonLinear(loop_process_number);
	// if(m_pcs->non_linear)
	//  error = m_pcs->ExecuteNonLinear();
	// else
	//  error = m_pcs->Execute();

	// CB This is a cheat to map a 2D horizontal heat pump distribution on a vertical model
	if (REACTINT_vec.size() > 0)
		if (REACTINT_vec[0]->heatpump_2DhTO2Dv)
			REACTINT_vec[0]->Heatpump_2DhTO2Dv_Mapping(true); // Map Heat plume CL to 2D vertical

	return error;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: MassTrasport
   Task: Similate heat transport
   Return: error
   Programming:
   07/2008 WW Extract from LOPTimeLoop_PCS();
   Modification:
   12.2008 WW Immigrtate the new functionalities  from loop_pcs.cpp
   -------------------------------------------------------------------------*/
inline double Problem::MassTrasport()
{
	double error = 1.0e+8;
	bool capvec = false;
	// bool prqvec = false;
	CRFProcess* m_pcs = total_processes[11];
	//
	if (!m_pcs->selected)
		return error; // 12.12.2008 WW
	// ClockTimeVec[0]->StartTime(); // CB time, remove if flow time is stopped

	for (int i = 0; i < (int)transport_processes.size(); i++)
	{
		m_pcs = transport_processes[i]; // 18.08.2008 WW
		// Component Mobile ?

		// MW reduce CONCENTRATION1 to non-negative values above water table for stability for Sugio approach with
		// RICHARDS
		if (mmp_vector[0]->permeability_saturation_model[0] == 10)
		{
			int nidx0 = m_pcs->GetNodeValueIndex("CONCENTRATION1", 0);
			CRFProcess* local_richards_flow = PCSGet("PRESSURE1", true);
			if (local_richards_flow != NULL)
			{
				int nidx1 = local_richards_flow->GetNodeValueIndex("PRESSURE1", 0);
				for (size_t j = 0; j < m_pcs->m_msh->GetNodesNumber(false); j++)
				{
					double local_conc = m_pcs->GetNodeValue(j, nidx0 + 1);
					double local_pressure = local_richards_flow->GetNodeValue(j, nidx1 + 1);
					if (local_pressure < 0 && local_conc < 0)
						m_pcs->SetNodeValue(j, nidx0 + 1, 0);
				}
			}
		}

		if (CPGetMobil(m_pcs->GetProcessComponentNumber()) > 0)
			error = m_pcs->ExecuteNonLinear(
			    loop_process_number); // NW. ExecuteNonLinear() is called to use the adaptive time step scheme

		int component = m_pcs->pcs_component_number;
		CompProperties* m_cp = cp_vec[component];

		if (m_cp->OutputMassOfComponentInModel == 1)
		{ // 05/2012 BG
			// TODO: Check
			if (singlephaseflow_process.size() > 0)
				OutputMassOfComponentInModel(singlephaseflow_process, m_pcs);
			else
				OutputMassOfComponentInModel(multiphase_processes, m_pcs);
		}
	}

	// Calculate Chemical reactions, after convergence of flow and transport
	// Move inside iteration loop if couplingwith transport is implemented SB:todo
	// SB:todo move into Execute Reactions	  if((aktueller_zeitschritt % 1) == 0)
	// REACT *rc = NULL; //OK
	// rc = REACT_vec[0]; //OK
	//				if(rc->flag_pqc) rc->ExecuteReactions();
	//				delete rc;
	// CB2406 #ifdef OGS_FEM_CAP // CAP_REACT  //CB merge CAP 0311
	if (REACT_CAP_vec.size() > 0)
		capvec = true;

	// if( (KinReactData_vector.size() > 0) || (REACT_vec.size()>0) ||  capvec || (REACT_PRQ_vec.size()>0) )  //CB merge
	// CAP 0311
	if (REACTINT_vec.size() > 0)
		REACTINT_vec[0]->ReactionPreProcessing();

	if (KinReactData_vector.size() > 0)
	{ // WW moved the following lines into this curly braces. 12.12.2008
		ClockTimeVec[0]->StopTime("Transport"); // CB time
		ClockTimeVec[0]->StartTime(); // CB time
		// Calculate Chemical reactions, after convergence of flow and transport
		// Move inside iteration loop if couplingwith transport is implemented SB:todo
		// First calculate kinetic reactions

		// DL --> for liquid phase only,
		// this will provide activities and speciation for kinetics
		if (capvec && KinReactData_vector[0]->NumberMineralkinetics > 0)
			REACT_CAP_vec[0]->ExecuteReactionsChemApp(-1, -1);

		KinReactData_vector[0]->ExecuteKinReact();
		ClockTimeVec[0]->StopTime("KinReactions"); // CB time
		ClockTimeVec[0]->StartTime(); // CB time
	}
	if (REACT_vec.size() > 0) // OK
	{
		if (REACT_vec[0]->flag_pqc)
		{
#ifdef REACTION_ELEMENT
			REACT_vec[0]->ExecuteReactionsPHREEQC0();
#else
// REACT_vec[0]->ExecuteReactions();

#ifdef OGS_FEM_IPQC
			REACT_vec[0]->ExecutePQCString();
#elif LIBPHREEQC
			// MDL: built-in phreeqc
			REACT_vec[0]->ExecuteReactionsPHREEQCNewLib();
#else
			REACT_vec[0]->ExecuteReactionsPHREEQCNew();
#endif // LIBPHREEQC
#endif // REACTION_ELEMENT
		}
	}

	if (REACT_CAP_vec.size() > 0)
	{ // use ChemApp
		REACT_CAP_vec[0]->ExecuteReactionsChemApp(1, -1); // DL
		// REACT_CAP_vec[0]->ConvertIC2BC();
	}

#ifdef GEM_REACT
	else // WW moved these pare of curly braces inside  ifdef GEM_REACT
	    if (m_vec_GEM->initialized_flag == 1) // when it was initialized.
	{
		int m_time = 1; // 0-previous time step results; 1-current time step results

		// Check if the Sequential Iterative Scheme needs to be intergrated
		if (m_pcs->m_num->cpl_max_iterations > 1)
			m_vec_GEM->flag_iterative_scheme = 1; // set to standard iterative scheme;
		// Move current xDC to previous xDC
		m_vec_GEM->CopyCurXDCPre();
		// Get info from MT
		// m_vec_GEM->ConvPorosityNodeValue2Elem(); //
		// second arguments should be one if we work with concentrations
		m_vec_GEM->GetReactInfoFromMassTransport(m_time);
		// m_vec_GEM->ConcentrationToMass();
		m_vec_GEM->Run_MainLoop(); // Run GEM
		m_vec_GEM->CopyCurBPre();

		// m_vec_GEM->MassToConcentration();
		// Calculate the different of xDC
		m_vec_GEM->UpdateXDCChemDelta();
		// Set info in MT
		m_vec_GEM->SetReactInfoBackMassTransport(m_time);
		// m_vec_GEM->ConvPorosityNodeValue2Elem(); // update element porosity and push back values
	}
#endif // GEM_REACT

#ifdef CHEMAPP
	if (Eqlink_vec.size() > 0)
		Eqlink_vec[0]->ExecuteEQLINK();
#endif
#ifdef BRNS
	if (m_vec_BRNS->init_flag == true)
		m_vec_BRNS->RUN(dt /*time value in seconds*/);
#endif

	// if(KinReactData_vector.size() > 0)  //12.12.2008 WW
	ClockTimeVec[0]->StopTime("EquiReact"); // CB time

	return error;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: FluidMomentum()
   Task:
   Return: error
   Programming:
   07/2008 WW Extract from LOPTimeLoop_PCS();
   Modification:
   -------------------------------------------------------------------------*/
inline double Problem::FluidMomentum()
{
	double error = 1.0e+8;
	CRFProcess* m_pcs = total_processes[9];
	//
	if (!m_pcs->selected)
		return error; // 12.12.2008 WW

	CFluidMomentum* fm_pcs = NULL; // by PCH
	//
	CFEMesh* m_msh = fem_msh_vector[0]; // Something must be done later on here.

	fm_pcs = m_msh->fm_pcs;
	fm_pcs->Execute(loop_process_number);

	// Switch off rechard flow if
	if (m_pcs->num_type_name.compare("STEADY") == 0 && aktueller_zeitschritt > 1)
	{
		// Turn off FLUID_MOMENTUM
		m_pcs->selected = false;
		// Turn off RICHARDS_FLOW
		m_pcs = PCSGet("RICHARDS_FLOW");
		if (m_pcs)
			m_pcs->selected = false;
		// Turn off LIQUID_FLOW
		m_pcs = PCSGet("LIQUID_FLOW");
		if (m_pcs)
			m_pcs->selected = false;
		// Turn off GROUNDWATER_FLOW
		m_pcs = PCSGet("GROUNDWATER_FLOW");
		if (m_pcs)
			m_pcs->selected = false;
	}
	//
	// error = 0.0 // JT... in unsteady flow, setting error=0.0 corresponds to error_cpl=0.0, and the coupling loop
	// ceases before RWPT is performed
	//            // What is the correct way to handle this, rather than setting error=1.e8???
	return error;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: RandomWalker()
   Task:
   Return: error
   Programming:
   07/2008 WW Extract from LOPTimeLoop_PCS();
   Modification:
   12.2008 WW
   -------------------------------------------------------------------------*/
inline double Problem::RandomWalker()
{
	double error = 1.0e+8;
	//
	CRFProcess* m_pcs = total_processes[10];
	//
	if (!m_pcs->selected)
		return error; // 12.12.2008 WW
	//
	// CFEMesh* m_msh = NULL;

	if (m_pcs && m_pcs->selected)
	{
		lop_coupling_iterations = 1;

		// Mount the proper mesh
		CFEMesh* m_msh = NULL;
		for (int i = 0; i < (int)pcs_vector.size(); ++i)
		{
			m_pcs = pcs_vector[i];

			// Select the mesh whose process name has the mesh for Fluid_Momentum
			//			if( m_pcs->pcs_type_name.find("RICHARDS_FLOW")!=string::npos)
			// TF
			if (m_pcs->getProcessType() == FiniteElement::RICHARDS_FLOW)
				m_msh = FEMGet("RICHARDS_FLOW");
			//			else if( m_pcs->pcs_type_name.find("LIQUID_FLOW")!=string::npos)
			// TF
			else if (m_pcs->getProcessType() == FiniteElement::LIQUID_FLOW)
				m_msh = FEMGet("LIQUID_FLOW");
			//			else if( m_pcs->pcs_type_name.find("GROUNDWATER_FLOW")!=string::npos)
			// TF
			else if (m_pcs->getProcessType() == FiniteElement::GROUNDWATER_FLOW)
				m_msh = FEMGet("GROUNDWATER_FLOW");
		}

		if (!m_msh)
			m_msh = FEMGet("FLUID_MOMENTUM");
		RandomWalk* rw_pcs = NULL; // By PCH
		rw_pcs = m_msh->PT;
		if (rw_pcs->getFlowPCS())
		{
			if (rw_pcs->getFlowPCS()->cal_integration_point_value) // WW
				rw_pcs->getFlowPCS()->Extropolation_GaussValue();
		}

		// Do I need velocity fileds solved by the FEM?
		if (m_pcs->tim_type == TimType::PURERWPT)
		{
			rw_pcs->PURERWPT = 1;
			char* dateiname = NULL;
			int sizeOfWord = 100;
			dateiname = (char*)malloc(sizeOfWord * sizeof(char));

			std::string filename = FileName;
			for (int i = 0; i <= (int)filename.size(); ++i)
				dateiname[i] = filename[i];

			rw_pcs->ReadInVelocityFieldOnNodes(dateiname);

			delete[] dateiname;
		}

		// Set the mode of the RWPT method
		if (m_pcs->num_type_name.compare("HETERO") == 0)
		{
			rw_pcs->RWPTMode = 1; // Set it for heterogeneous media
			std::cout << "RWPT is on " << m_pcs->num_type_name << " mode."
			          << "\n";
		}
		else if (m_pcs->num_type_name.compare("HOMO_ADVECTION") == 0)
		{
			rw_pcs->RWPTMode = 2;
			std::cout << "RWPT is on " << m_pcs->num_type_name << " mode."
			          << "\n";
		}
		else if (m_pcs->num_type_name.compare("HETERO_ADVECTION") == 0)
		{
			rw_pcs->RWPTMode = 3;
			std::cout << "RWPT is on " << m_pcs->num_type_name << " mode."
			          << "\n";
		}
		else if (m_pcs->num_type_name.compare("HOMO_DISPERSION") == 0)
		{
			rw_pcs->RWPTMode = 4;
			std::cout << "RWPT is on " << m_pcs->num_type_name << " mode."
			          << "\n";
		}
		else if (m_pcs->num_type_name.compare("HETERO_DISPERSION") == 0)
		{
			rw_pcs->RWPTMode = 5;
			std::cout << "RWPT is on " << m_pcs->num_type_name << " mode."
			          << "\n";
		}
		else if (m_pcs->num_type_name.compare("HETERO_FDM") == 0)
		{
			rw_pcs->RWPTMode = 1; // Set it for heterogeneous media
			std::cout << "RWPT is on " << m_pcs->num_type_name << " mode."
			          << "\n";
		}
		else if (m_pcs->num_type_name.compare("HOMO_ADVECTION_FDM") == 0)
		{
			rw_pcs->RWPTMode = 2;
			std::cout << "RWPT is on " << m_pcs->num_type_name << " mode."
			          << "\n";
		}
		else if (m_pcs->num_type_name.compare("HETERO_ADVECTION_FDM") == 0)
		{
			rw_pcs->RWPTMode = 3;
			std::cout << "RWPT is on " << m_pcs->num_type_name << " mode."
			          << "\n";
		}
		else if (m_pcs->num_type_name.compare("HOMO_DISPERSION_FDM") == 0)
		{
			rw_pcs->RWPTMode = 4;
			std::cout << "RWPT is on " << m_pcs->num_type_name << " mode."
			          << "\n";
		}
		else if (m_pcs->num_type_name.compare("HETERO_DISPERSION_FDM") == 0)
		{
			rw_pcs->RWPTMode = 5;
			std::cout << "RWPT is on " << m_pcs->num_type_name << " mode."
			          << "\n";
		}
		else // HOMO Advection + Dispersion
		{
			rw_pcs->RWPTMode = 0;
			std::cout << "RWPT is on HOMO_ADVECTION_DISPERSION mode."
			          << "\n";
		}

		if (m_pcs->num_type_name.find("FDM") != std::string::npos)
		{
			rw_pcs->PURERWPT = 2;
			if (rw_pcs->FDMIndexSwitch == 0)
			{
				rw_pcs->buildFDMIndex();
				// Switch off
				rw_pcs->FDMIndexSwitch = 1;
			}
		}

		if (rwpt_numsplits < 0)
			rwpt_numsplits = 10; // JT 2010 set default value, unless specified in .tim input file

		rw_pcs->AdvanceBySplitTime(dt, rwpt_numsplits);
		//	rw_pcs->TraceStreamline(); // JT, no longer needed
		print_result = true;
		rw_pcs->RandomWalkOutput(aktuelle_zeit, aktueller_zeitschritt);
	}

	return 0.0;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: Deformation
   Task: Similate deformation
   Return: error
   Programming:
   07/2008 WW Extract from LOPTimeLoop_PCS();
   Modification:
   -------------------------------------------------------------------------*/
inline double Problem::Deformation()
{
	CRFProcessDeformation* dm_pcs = NULL;
	double error = 1.0e+8;
	CRFProcess* m_pcs = total_processes[12];
	//
	dm_pcs = (CRFProcessDeformation*)(m_pcs);
	error = dm_pcs->Execute(loop_process_number);

	if (dm_pcs->pcs_type_name_vector.size() > 0 && dm_pcs->pcs_type_name_vector[0].find("DYNAMIC") != std::string::npos)
	{
		return error;
	}
	// Error
	if (dm_pcs->type / 10 == 4)
	{
		if (!dm_pcs->isDynamic())
		{
			m_pcs->cal_integration_point_value = true;
			dm_pcs->CalIntegrationPointValue();
		}
		if (dm_pcs->type == 42) // H2M. 07.2011. WW
			dm_pcs->CalcSecondaryVariablesUnsaturatedFlow();
	}
	return error;
}

/**************************************************************************
   FEMLib-Method:
   02/2005 OK Implementation
   08/2005 WW Changes due to geometry objects applied
   08/2005 MB Changes ... (OK to what ?)
   04/2006 OK and once again ...
   07/2008 WW Extract from LOPTimeLoop_PCS();
   05.2009 WW For surface-soil-ground coupled model
**************************************************************************/
inline void Problem::LOPExecuteRegionalRichardsFlow(CRFProcess* m_pcs_global, int loop_process_number)
{
	int j, k;
	long i;
	MeshLib::CElem* m_ele = NULL;
	// TF not used	MeshLib::CNode* m_nod = NULL;
	int no_local_elements = m_pcs_global->m_msh->getNumberOfMeshLayers();
	int no_local_nodes = no_local_elements + 1;
	long g_element_number, g_node_number;
	CFEMesh* m_msh_local = NULL;
	CRFProcess* m_pcs_local = NULL;
	Math_Group::vec<MeshLib::CNode*> ele_nodes(20);
	double value;
#if defined(USE_MPI_REGSOIL)
	double* values;
	int rp;
	int num_parallel_blocks;
	int l;
#endif
	int timelevel = 1;
	int idxp = m_pcs_global->GetNodeValueIndex("PRESSURE1") + timelevel;
	// WW int idxcp = m_pcs_global->GetNodeValueIndex("PRESSURE_CAP") + timelevel;
	int idxS = m_pcs_global->GetNodeValueIndex("SATURATION1") + timelevel;
	MeshLib::CElem* m_ele_local = NULL;
	MeshLib::CNode* m_nod_local = NULL;

#if defined(USE_MPI_REGSOIL)
	values = new double[no_local_nodes]; // Should be more sophisticated
#endif

	//======================================================================
	if (aktueller_zeitschritt == 1)
	{
		//--------------------------------------------------------------------
		// Create local RICHARDS process
		std::cout << "    Create local RICHARDS process"
		          << "\n";
		m_pcs_local = new CRFProcess();
		m_pcs_local->isRSM = true; // WW

		//    m_pcs_local->pcs_type_name = m_pcs_global->pcs_type_name;
		// TF
		m_pcs_local->setProcessType(m_pcs_global->getProcessType());
		m_pcs_local->num_type_name = m_pcs_global->num_type_name;
		m_pcs_local->cpl_type_name = m_pcs_global->cpl_type_name;
		m_pcs_local->Write_Matrix = m_pcs_global->Write_Matrix;
		m_pcs_local->pcs_type_number = (int)pcs_vector.size();
		m_pcs_local->Config();
		m_pcs_local->setProblemObjectPointer(m_pcs_global->getProblemObjectPointer()); // WW 18.08.2011. WW
		//--------------------------------------------------------------------
		// Create local MSH
		m_msh_local = new CFEMesh(m_pcs_global->m_msh->getGEOObjects(), m_pcs_global->m_msh->getProjectName());
		m_msh_local->geo_name = "RICHARDS_FLOW_LOCAL";
		m_msh_local->setElementType(MshElemType::LINE);
		m_msh_local->setNumberOfMeshLayers(m_pcs_global->m_msh->getNumberOfMeshLayers());
		//....................................................................
		m_msh_local->ele_vector.resize(no_local_elements);
		for (j = 0; j < no_local_elements; j++)
		{
			m_ele = m_pcs_global->m_msh->ele_vector[j];
			m_ele_local = new MeshLib::CElem(j, m_ele);
			for (k = 0; k < 2; k++) // ele_type
				m_ele_local->getNodeIndices()[k] = j + k;
			m_msh_local->ele_vector[j] = m_ele_local;
		}
		m_msh_local->nod_vector.resize(no_local_nodes);
		m_msh_local->Eqs2Global_NodeIndex.resize(no_local_nodes);
		for (j = 0; j < no_local_nodes; j++)
		{
			// TF        m_nod = m_pcs_global->m_msh->nod_vector[j];
			m_nod_local = new MeshLib::CNode(j, m_pcs_global->m_msh->nod_vector[j]->getData());
			// m_nod_local = m_nod;
			m_msh_local->nod_vector[j] = m_nod_local;
			m_msh_local->nod_vector[j]->SetEquationIndex(j);
			m_msh_local->Eqs2Global_NodeIndex[j] = m_msh_local->nod_vector[j]->GetIndex();
		}
		m_msh_local->ConstructGrid();
		m_msh_local->FillTransformMatrix();
		m_msh_local->FaceNormal();
		//....................................................................
		fem_msh_vector.push_back(m_msh_local);
		//....................................................................
		m_pcs_local->m_msh = m_msh_local;
		m_pcs_local->Create();

		m_pcs_local->fem->setShapeFunctionPool(m_pcs_global->fem->getShapeFunctionPool(0),
		                                       m_pcs_global->fem->getShapeFunctionPool(1));

		m_pcs_local->SetBoundaryConditionAndSourceTerm();
		//....................................................................
		// BC
		//....................................................................
		// ST
		/*
		    for(s=0;s<(int)m_pcs_global->st_node_value.size();s++)
		    {
		      m_nod_value = new CNodeValue(); //OK
		      //m_nod_value = m_pcs_global->st_node_value[s];
		      m_nod_value->node_value = m_pcs_global->st_node_value[s]->node_value;
		      m_nod_value->msh_node_number = m_pcs_global->st_node_value[s]->msh_node_number;
		      m_nod_value->geo_node_number =  m_nod_value->msh_node_number; //WW
		      m_pcs_local->st_node_value.push_back(m_nod_value);
		    }
		 */
		m_pcs_local->st_node_value.clear();
		m_pcs_local->st_node.clear();
		for (j = 0; j < (int)m_pcs_global->st_node_value.size(); j++)
			m_pcs_local->st_node_value.push_back(m_pcs_global->st_node_value[j]);
		for (j = 0; j < (int)m_pcs_global->st_node.size(); j++)
			m_pcs_local->st_node.push_back(m_pcs_global->st_node[j]);
		//....................................................................
		pcs_vector.push_back(m_pcs_local);
	}
	//======================================================================
	else
		for (i = 0; i < (long)m_pcs_global->m_msh->nod_vector.size(); i++)
		{
			value = m_pcs_global->GetNodeValue(i, idxp);
			m_pcs_global->SetNodeValue(i, idxp - 1, value);
			// WW value = m_pcs_global->GetNodeValue(i,idxcp);
			// WW m_pcs_global->SetNodeValue(i,idxcp-1,value);
			value = m_pcs_global->GetNodeValue(i, idxS);
			m_pcs_global->SetNodeValue(i, idxS - 1, value);
		}
	//======================================================================
	std::cout << "\n      ================================================"
	          << "\n";
	std::cout << "    ->Process " << loop_process_number << ": "
	          << "REGIONAL_" << convertProcessTypeToString(m_pcs_global->getProcessType()) << "\n";
	std::cout << "      ================================================"
	          << "\n";
	int no_richards_problems = (int)(m_pcs_global->m_msh->ele_vector.size() / no_local_elements);

	//--- For couping with ground flow process. WW
	int idx_v;
	MeshLib::GridsTopo* neighb_grid = NULL;
	CRFProcess* neighb_pcs = total_processes[1];

	if (neighb_pcs)
		for (i = 0; i < (int)neighb_pcs->m_msh->grid_neighbors.size(); i++)
		{
			neighb_grid = neighb_pcs->m_msh->grid_neighbors[i];
			if (neighb_grid->getNeighbor_Name().find("SECTOR_SOIL") != std::string::npos)
				break;
			neighb_grid = NULL;
		}
//------------
#ifndef USE_MPI_REGSOIL
	for (i = 0; i < no_richards_problems; i++)
	// for(i=0;i<2;i++)
	{
		if (i > 0)
			std::cout << "      ================================================"
			          << "\n";
		std::cout << "      ->Column number: " << i << "\n";
		std::cout << "      ================================================"
		          << "\n";
		m_pcs_local = pcs_vector[(int)pcs_vector.size() - 1];
		m_pcs_local->pcs_number = i;
		m_msh_local = fem_msh_vector[(int)fem_msh_vector.size() - 1];
		//....................................................................
		// Set local NODs
		for (j = 0; j < no_local_nodes; j++)
		{
			g_node_number = j + (i * no_local_nodes);
			// TF not used			m_nod = m_pcs_global->m_msh->nod_vector[g_node_number];
			m_nod_local = m_msh_local->nod_vector[j];
			// m_nod_local = m_nod;
			m_nod_local->getConnectedElementIDs().push_back(i);
			// m_nod_local->ok_dummy = i;
		}
		//....................................................................
		// Set local ELEs
		for (j = 0; j < no_local_elements; j++)
		{
			g_element_number = j + (i * no_local_elements);
			m_ele = m_pcs_global->m_msh->ele_vector[g_element_number];
			m_ele_local = m_msh_local->ele_vector[j];
			m_ele_local->SetPatchIndex(m_ele->GetPatchIndex());
		}
		//....................................................................
		// Set ICs
		if (aktueller_zeitschritt == 1) // YD

			for (j = 0; j < no_local_nodes; j++)
			{
				g_node_number = j + (i * no_local_nodes);
				value = m_pcs_global->GetNodeValue(g_node_number, idxp);
				m_pcs_local->SetNodeValue(j, idxp - 1, value);
				m_pcs_local->SetNodeValue(j, idxp, value);
				// WW value = m_pcs_global->GetNodeValue(g_node_number,idxcp);
				// WW m_pcs_local->SetNodeValue(j,idxcp-1,value);
				// WW  m_pcs_local->SetNodeValue(j,idxcp,value);
				value = m_pcs_global->GetNodeValue(g_node_number, idxS);
				m_pcs_local->SetNodeValue(j, idxS - 1, value);
				m_pcs_local->SetNodeValue(j, idxS, value);
			}
		else
			for (j = 0; j < no_local_nodes; j++)
			{
				g_node_number = j + (i * no_local_nodes);
				value = m_pcs_global->GetNodeValue(g_node_number, idxp);
				m_pcs_local->SetNodeValue(j, idxp - 1, value);
				// value = m_pcs_global->GetNodeValue(g_node_number,idxcp);
				// m_pcs_local->SetNodeValue(j,idxcp-1,value);
				value = m_pcs_global->GetNodeValue(g_node_number, idxS);
				m_pcs_local->SetNodeValue(j, idxS - 1, value);
			}
		//....................................................................
		// Set local BCs
		// WW m_pcs_local->CreateBCGroup();
		//....................................................................
		// Set local STs
		// WW m_pcs_local->CreateSTGroup();
		// look for corresponding OF-triangle
		// m_ele_of = m_msh_of->GetElement(m_pnt_sf);
		//....................................................................
		//---- Source term from neighbor process. 25.05.2009. WW
		idx_v = m_pcs_local->GetNodeValueIndex("VELOCITY_Z1");
		if (neighb_grid)
		{
			CSourceTerm* m_st = NULL;
			CNodeValue* m_nod_val = NULL;
			if (aktueller_zeitschritt == 1)
			{
				m_st = new CSourceTerm();
				m_nod_val = new CNodeValue();
				m_pcs_local->st_node_value.push_back(m_nod_val);
				m_pcs_local->st_node.push_back(m_st);
			}
			else
			{
				m_st = m_pcs_local->st_node[(int)m_pcs_local->st_node.size() - 1];
				m_nod_val = m_pcs_local->st_node_value[(int)m_pcs_local->st_node_value.size() - 1];
			}
			// WW long *local_indxs = NULL;
			// WW local_indxs = neighb_grid->getBorderNodeIndicies();
			m_nod_val->msh_node_number = no_local_nodes - 1;
			m_nod_val->geo_node_number = no_local_nodes - 1;
			//       m_nod_val->node_value = -neighb_pcs->GetNodeValue(local_indxs[i],
			//       neighb_pcs->GetNodeValueIndex("VELOCITY_Z1"))
			//                               /neighb_pcs->time_unit_factor;
			m_nod_val->node_value = -m_pcs_global->GetNodeValue(i, m_pcs_global->GetNodeValueIndex("VELOCITY_Z1"))
			                        / neighb_pcs->time_unit_factor;

			m_pcs_global->SetNodeValue(i, m_pcs_global->GetNodeValueIndex("VELOCITY_Z1"), 0.);
		}
		//-----------------------------------------------------

		m_pcs_local->ExecuteNonLinear(loop_process_number, false);
		// Velocity. 22.05.2009. WW
		m_pcs_local->CalIntegrationPointValue();
		m_pcs_local->Extropolation_GaussValue();
		m_pcs_local->cal_integration_point_value = false;
		//....................................................................
		// Store results in global PCS tables
		for (j = 0; j < no_local_nodes; j++)
		{
			g_node_number = j + (i * no_local_nodes);
			value = m_pcs_local->GetNodeValue(j, idxp);
			m_pcs_global->SetNodeValue(g_node_number, idxp, value);
			// value = m_pcs_local->GetNodeValue(j,idxcp);
			// m_pcs_global->SetNodeValue(g_node_number,idxcp,value);
			value = m_pcs_local->GetNodeValue(j, idxS);
			m_pcs_global->SetNodeValue(g_node_number, idxS, value);
			// WW. 22.05.2009
			value = m_pcs_local->GetNodeValue(j, idx_v);
			m_pcs_global->SetNodeValue(g_node_number, idx_v, value);
		}
		// m_pcs->m_msh = FEMGet("RICHARDS_FLOW");
		// pcs->m_msh->RenumberNodesForGlobalAssembly(); related to Comment 1  // MB/WW
		if (!m_pcs_local->TimeStepAccept())
		{
			m_pcs_global->accepted = false;
			break;
		}
	}
#else // ifndef USE_MPI_REGSOIL
	num_parallel_blocks = no_richards_problems / size;

#ifdef TRACE
	std::cout << "Num parallel blocks: " << num_parallel_blocks << "\n";
#endif

	for (i = 0; i < num_parallel_blocks + 1; i++)
	// for(i=0;i<2;i++)
	{
		if (i * size + myrank < no_richards_problems) // Do a parallel block
		{
			rp = i * size + myrank;
			m_pcs_local = pcs_vector[(int)pcs_vector.size() - 1];
			m_pcs_local->pcs_number = rp;
			m_msh_local = fem_msh_vector[(int)fem_msh_vector.size() - 1];
			//....................................................................
			// Set local NODs
			for (j = 0; j < no_local_nodes; j++)
			{
				g_node_number = j + (rp * no_local_nodes);
				m_nod = m_pcs_global->m_msh->nod_vector[g_node_number];
				m_nod_local = m_msh_local->nod_vector[j];
				// m_nod_local = m_nod;
				// ????
				m_nod_local->connected_elements.push_back(rp);
			}
			//....................................................................
			// Set local ELEs
			for (j = 0; j < no_local_elements; j++)
			{
				g_element_number = j + (rp * no_local_elements);
				m_ele = m_pcs_global->m_msh->ele_vector[g_element_number];
				m_ele_local = m_msh_local->ele_vector[j];
				m_ele_local->SetPatchIndex(m_ele->GetPatchIndex());
			}
			//....................................................................
			// Set ICs
			if (aktueller_zeitschritt == 1) // YD

				for (j = 0; j < no_local_nodes; j++)
				{
					g_node_number = j + (rp * no_local_nodes);
					value = m_pcs_global->GetNodeValue(g_node_number, idxp);
					m_pcs_local->SetNodeValue(j, idxp - 1, value);
					m_pcs_local->SetNodeValue(j, idxp, value);
					value = m_pcs_global->GetNodeValue(g_node_number, idxcp);
					m_pcs_local->SetNodeValue(j, idxcp - 1, value);
					m_pcs_local->SetNodeValue(j, idxcp, value);
					value = m_pcs_global->GetNodeValue(g_node_number, idxS);
					m_pcs_local->SetNodeValue(j, idxS - 1, value);
					m_pcs_local->SetNodeValue(j, idxS, value);
				}
			else
				for (j = 0; j < no_local_nodes; j++)
				{
					g_node_number = j + (rp * no_local_nodes);
					value = m_pcs_global->GetNodeValue(g_node_number, idxp);
					m_pcs_local->SetNodeValue(j, idxp - 1, value);
					value = m_pcs_global->GetNodeValue(g_node_number, idxcp);
					m_pcs_local->SetNodeValue(j, idxcp - 1, value);
					value = m_pcs_global->GetNodeValue(g_node_number, idxS);
					m_pcs_local->SetNodeValue(j, idxS - 1, value);
				}
			//....................................................................
			// Set local BCs
			m_pcs_local->CreateBCGroup();
//....................................................................
// Set local STs
// m_pcs_local->CreateSTGroup();
// look for corresponding OF-triangle
// m_ele_of = m_msh_of->GetElement(m_pnt_sf);
//....................................................................
#ifdef TRACE
			std::cout << "Executing local process."
			          << "\n";
#endif
			m_pcs_local->ExecuteNonLinear(loop_process_number);
		} // End of parallel block
		//....................................................................
		// Store results in global PCS tables
		for (int k = 0; k < size; k++)
		{
			rp = i * size + k;
			if (rp < no_richards_problems)
			{
				if (myrank == k)
				{
					// idxp
					for (l = 0; l < no_local_nodes; l++)
						values[l] = m_pcs_local->GetNodeValue(l, idxp);
					MPI_Bcast((void*)values, no_local_nodes, MPI_DOUBLE, k, comm_DDC);
					for (l = 0; l < no_local_nodes; l++)
						m_pcs_global->SetNodeValue(l + rp * no_local_nodes, idxp, values[l]);
					// idxcp
					for (l = 0; l < no_local_nodes; l++)
						values[l] = m_pcs_local->GetNodeValue(l, idxcp);
					MPI_Bcast((void*)values, no_local_nodes, MPI_DOUBLE, k, comm_DDC);
					for (l = 0; l < no_local_nodes; l++)
						m_pcs_global->SetNodeValue(l + rp * no_local_nodes, idxcp, values[l]);
					// idxS
					for (l = 0; l < no_local_nodes; l++)
						values[l] = m_pcs_local->GetNodeValue(l, idxS);
					MPI_Bcast((void*)values, no_local_nodes, MPI_DOUBLE, k, comm_DDC);
					for (l = 0; l < no_local_nodes; l++)
						m_pcs_global->SetNodeValue(l + rp * no_local_nodes, idxS, values[l]);
				}
				else
				{
					// idxp
					MPI_Bcast((void*)values, no_local_nodes, MPI_DOUBLE, k, comm_DDC);
					for (l = 0; l < no_local_nodes; l++)
						m_pcs_global->SetNodeValue(l + rp * no_local_nodes, idxp, values[l]);
					// idxcp
					MPI_Bcast((void*)values, no_local_nodes, MPI_DOUBLE, k, comm_DDC);
					for (l = 0; l < no_local_nodes; l++)
						m_pcs_global->SetNodeValue(l + rp * no_local_nodes, idxcp, values[l]);
					// idxS
					MPI_Bcast((void*)values, no_local_nodes, MPI_DOUBLE, k, comm_DDC);
					for (l = 0; l < no_local_nodes; l++)
						m_pcs_global->SetNodeValue(l + rp * no_local_nodes, idxS, values[l]);
				}
			}
		}
	}
#endif
//----------------------------------------------------------------------

#if defined(USE_MPI_REGSOIL)
	delete[] values;
#endif
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   08/2003 SB Implementation
   01/2004 MX k_eff
   11/2004 OK fluid mass fluxes
   08/2008 WW Extract from LOPTimeLoop_PCS();
   last modification:
**************************************************************************/
void Problem::LOPCalcELEResultants()
{
	size_t no_processes = pcs_vector.size();
	CRFProcess* m_pcs = NULL;

	for (size_t p = 0; p < no_processes; p++)
	{
		m_pcs = pcs_vector[p];
		if (!m_pcs->selected) // OK4108
			continue;
		// cout << "LOPCalcELEResultants: " << m_pcs->pcs_type_name << "\n";
		// TF
		std::string pcs_type_name(convertProcessTypeToString(m_pcs->getProcessType()));
		//		switch (m_pcs->pcs_type_name[0]) {
		switch (pcs_type_name[0])
		{
			default:
				break;
			case 'L': // Liquid flow
				m_pcs->CalcELEVelocities();
				break;
			case 'G': // Groundwater flow
				m_pcs->CalcELEVelocities();
				break;
			case 'A': // Gas flow
				m_pcs->CalcELEVelocities();
				// m_pcs->CalcELEMassFluxes();			// BG
				break;
			case 'T': // Two-phase flow
				break;
			case 'C': // Componental flow
				break;
			case 'H': // Heat transport
				break;
			case 'M': // Mass transport
				break;
			case 'D': // Deformation
				break;
			case 'O': // Overland flow
				m_pcs->CalcELEVelocities();
				break;
			case 'R': // Richards flow
				m_pcs->CalcELEVelocities();
				break;
			case 'F': // Fluid Momentum
				break;
		}
	}
}

/**************************************************************************
   ROCKFLOW - Funktion: ASMCalcNodeWDepth

   Task:
   Berechnung und Speichern der Knotenfl?se
   Parameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E: long i: node index
   Result:
   - void -

   Programmaenderungen:
   11/2002   MB/OK  Implementation
   10/2004   MB     PCS
   12/2008   WW     Encapsulate to this class
**************************************************************************/
inline void Problem::ASMCalcNodeWDepth(CRFProcess* m_pcs)
{
	int nidx, nidy, nidz;
	// OK411 int timelevel = 1;
	double WDepth;

	nidx = m_pcs->GetNodeValueIndex("HEAD") + 1;
	nidy = m_pcs->GetNodeValueIndex("WDEPTH");
	nidz = m_pcs->GetNodeValueIndex("COUPLING");

	const size_t n_nodes(m_pcs->m_msh->nod_vector.size());
	for (size_t nn = 0; nn < n_nodes; nn++)
	{
		WDepth = m_pcs->GetNodeValue(nn, nidx) - m_pcs->m_msh->nod_vector[nn]->getData()[2];
		// JOD only needed for GREEN_AMPT source term
		m_pcs->SetNodeValue(nn, nidz, m_pcs->GetNodeValue(nn, nidz + 1));
		if (WDepth < 0.0)
			WDepth = 0.0;
		m_pcs->SetNodeValue(nn, nidy, WDepth);
	}
}

/**************************************************************************/
/* ROCKFLOW - Funktion: PCSCalcSecondaryVariables
 */
/* Aufgabe:
   Berechung von secondary variables w?rend der Zeitschleife
   Abfrage je nach Prozess, der aktiv ist
 */
/* Programmaenderungen:
   08/2003   SB   Implementation
   01/2006   YD   add dual porosity
   01/2007 OK Two-phase flow
 */
/**************************************************************************/
void Problem::PCSCalcSecondaryVariables()
{
	// WW  long j;

	int i, ptype;
	CRFProcess* m_pcs = NULL;
	// OK411 CRFProcess* m_pcs_phase_1 = NULL;
	// OK411 CRFProcess* m_pcs_phase_2 = NULL;
	// WW int ndx_p_gas_old,ndx_p_gas_new,ndx_p_liquid_old,ndx_p_liquid_new,ndx_p_cap_old;
	//----------------------------------------------------------------------
	// OK411 bool pcs_cpl = true;
	//----------------------------------------------------------------------
	// Check if NAPLdissolution is modeled, required by MMPCalcSecondaryVariablesNew
	bool NAPLdiss = false;
	NAPLdiss = KNaplDissCheck();
	/* go through all processes */
	int no_processes = (int)pcs_vector.size();
	for (i = 0; i < no_processes; i++)
	{
		/* get process */
		// pcs = pcs->GetProcessByNumber(i+1);
		m_pcs = pcs_vector[i]; // JOD
		if (m_pcs != NULL)
		{
			ptype = m_pcs->GetObjType();
			switch (ptype)
			{
				case 1: /* Flow process */
					// do nothing
					break;
				case 66:
					// Temp mit pcs, only for test MB
					ASMCalcNodeWDepth(m_pcs);
					break;
				case 11: /* Non-isothermal flow process */
					break;
				case 13: /* Non-isothermal flow process */
					break;
				case 2: /* Mass transport process */
					break;
				case 3: /* Heat transport */
					// do nothing
					break;
				case 4: /* Deformation */
					// do nothing
					break;
				case 41: /* Deformation-flow coupled system in monolithic scheme */
					// do nothing
					break;
				case 12: /* Multi-phase flow process */
					MMPCalcSecondaryVariablesNew(m_pcs, NAPLdiss);
					// MMPCalcSecondaryVariablesNew(m_pcs);
					break;
				default:
					// std::cout << "PCSCalcSecondaryVariables - nothing to do" <<
					//"\n";
					break;
			}
		} // If
	} // while
}

/**************************************************************************
   FEMLib-Method:
   05/2009 OK Implementation
**************************************************************************/
bool Problem::Check()
{
	CRFProcess* m_pcs = NULL;
	for (int i = 0; i < (int)total_processes.size(); i++)
	{
		m_pcs = total_processes[i];
		if (!m_pcs->Check())
			return false;
	}
	return true;
}

void Problem::createShapeFunctionPool()
{
	CRFProcess* pcs_0c_fem = NULL;
	CRFProcess* pcs_1c_fem = NULL;
	for (std::size_t i = 0; i < pcs_vector.size(); i++)
	{
		CRFProcess* pcs = pcs_vector[i];
		if (pcs->getProcessType() == FiniteElement::FLUID_MOMENTUM
		    || pcs->getProcessType() == FiniteElement::RANDOM_WALK)
			continue;

		if (pcs->getProcessType() == FiniteElement::DEFORMATION
		    || pcs->getProcessType() == FiniteElement::DEFORMATION_DYNAMIC
		    || pcs->getProcessType() == FiniteElement::DEFORMATION_FLOW
		    || pcs->getProcessType() == FiniteElement::DEFORMATION_H2)
		{
			pcs_1c_fem = pcs;
		}
		else
		{
			pcs_0c_fem = pcs;
		}
	}

	CFiniteElementStd* lin_fem_assembler = NULL;
	CFiniteElementVec* fem_assembler = NULL;
	if (pcs_0c_fem)
	{
		lin_fem_assembler = pcs_0c_fem->getLinearFEMAssembler();
		lin_fem_assembler->setOrder(1);
	}
	if (pcs_1c_fem)
	{
		if (!lin_fem_assembler)
		{
			lin_fem_assembler = pcs_1c_fem->getLinearFEMAssembler();
			if (lin_fem_assembler)
				lin_fem_assembler->setOrder(1);
		}

		CRFProcessDeformation* dm_pcs = dynamic_cast<CRFProcessDeformation*>(pcs_1c_fem);
		fem_assembler = dm_pcs->GetFEMAssembler();
		fem_assembler->setOrder(2);
	}

	// Check element types of meshes
	std::vector<MshElemType::type> elem_types;
	elem_types.reserve(MshElemType::LAST);

	for (std::size_t i = 0; i < static_cast<std::size_t>(MshElemType::LAST); i++)
	{
		elem_types.push_back(MshElemType::INVALID);
	}

	for (std::size_t i = 0; i < fem_msh_vector.size(); i++)
	{
		MeshLib::CFEMesh* mesh = fem_msh_vector[i];
		if (mesh->getNumberOfLines() > 0)
			elem_types[static_cast<int>(MshElemType::LINE) - 1] = MshElemType::LINE;
		if (mesh->getNumberOfTris() > 0)
		{
			elem_types[static_cast<int>(MshElemType::TRIANGLE) - 1] = MshElemType::TRIANGLE;
			elem_types[static_cast<int>(MshElemType::LINE) - 1] = MshElemType::LINE;
		}
		if (mesh->getNumberOfQuads() > 0)
		{
			elem_types[static_cast<int>(MshElemType::QUAD) - 1] = MshElemType::QUAD;
			elem_types[static_cast<int>(MshElemType::LINE) - 1] = MshElemType::LINE;
		}
		if (mesh->getNumberOfHexs() > 0)
		{
			elem_types[static_cast<int>(MshElemType::HEXAHEDRON) - 1] = MshElemType::HEXAHEDRON;
			elem_types[static_cast<int>(MshElemType::QUAD8) - 1] = MshElemType::QUAD8;
			elem_types[static_cast<int>(MshElemType::LINE) - 1] = MshElemType::LINE;
		}
		if (mesh->getNumberOfTets() > 0)
		{
			elem_types[static_cast<int>(MshElemType::TETRAHEDRON) - 1] = MshElemType::TETRAHEDRON;
			elem_types[static_cast<int>(MshElemType::TRIANGLE) - 1] = MshElemType::TRIANGLE;
			elem_types[static_cast<int>(MshElemType::LINE) - 1] = MshElemType::LINE;
		}
		if (mesh->getNumberOfPrisms() > 0)
		{
			elem_types[static_cast<int>(MshElemType::PRISM) - 1] = MshElemType::PRISM;
			elem_types[static_cast<int>(MshElemType::QUAD8) - 1] = MshElemType::QUAD8;
			elem_types[static_cast<int>(MshElemType::TRIANGLE) - 1] = MshElemType::TRIANGLE;
			elem_types[static_cast<int>(MshElemType::LINE) - 1] = MshElemType::LINE;
		}
		if (mesh->getNumberOfPyramids() > 0)
		{
			elem_types[static_cast<int>(MshElemType::PYRAMID) - 1] = MshElemType::PYRAMID;
			elem_types[static_cast<int>(MshElemType::QUAD8) - 1] = MshElemType::QUAD8;
			elem_types[static_cast<int>(MshElemType::TRIANGLE) - 1] = MshElemType::TRIANGLE;
			elem_types[static_cast<int>(MshElemType::LINE) - 1] = MshElemType::LINE;
		}
	}

	const int num_gauss_sample_pnts = num_vector[0]->getNumIntegrationSamplePoints();
	if (lin_fem_assembler)
		_line_shapefunction_pool
		    = new FiniteElement::ShapeFunctionPool(elem_types, *lin_fem_assembler, num_gauss_sample_pnts);
	if (fem_assembler)
	{
		_quadr_shapefunction_pool
		    = new FiniteElement::ShapeFunctionPool(elem_types, *fem_assembler, num_gauss_sample_pnts);
		if (!_line_shapefunction_pool)
		{
			fem_assembler->setOrder(1);
			_line_shapefunction_pool
			    = new FiniteElement::ShapeFunctionPool(elem_types, *fem_assembler, num_gauss_sample_pnts);
			fem_assembler->setOrder(2);
		}
	}

	// Set ShapeFunctionPool
	for (std::size_t i = 0; i < pcs_vector.size(); i++)
	{
		CRFProcess* pcs = pcs_vector[i];
		if (pcs->getProcessType() == FiniteElement::FLUID_MOMENTUM
		    || pcs->getProcessType() == FiniteElement::RANDOM_WALK)
			continue;

		if (pcs->getProcessType() == FiniteElement::DEFORMATION)
		{
			CRFProcessDeformation* dm_pcs = dynamic_cast<CRFProcessDeformation*>(pcs);
			CFiniteElementVec* fem_assem_h = dm_pcs->GetFEMAssembler();
			fem_assem_h->setShapeFunctionPool(_line_shapefunction_pool, _quadr_shapefunction_pool);
		}
		else if (pcs->getProcessType() == FiniteElement::DEFORMATION_DYNAMIC
		         || pcs->getProcessType() == FiniteElement::DEFORMATION_FLOW
		         || pcs->getProcessType() == FiniteElement::DEFORMATION_H2)
		{
			CRFProcessDeformation* dm_pcs = dynamic_cast<CRFProcessDeformation*>(pcs);
			CFiniteElementVec* fem_assem_h = dm_pcs->GetFEMAssembler();
			fem_assem_h->setShapeFunctionPool(_line_shapefunction_pool, _quadr_shapefunction_pool);
			CFiniteElementStd* fem_assem = dm_pcs->getLinearFEMAssembler();
			fem_assem->setShapeFunctionPool(_line_shapefunction_pool, _quadr_shapefunction_pool);
		}
		else
		{
			CFiniteElementStd* fem_assem = pcs->getLinearFEMAssembler();
			if (!_quadr_shapefunction_pool)
				_quadr_shapefunction_pool = _line_shapefunction_pool;
			fem_assem->setShapeFunctionPool(_line_shapefunction_pool, _quadr_shapefunction_pool);
		}
	}
}

/**************************************************************************
   FEMLib-Method:
   06/2009 OK Implementation
**************************************************************************/
bool MODCreate()
{
	PCSConfig(); // OK
	if (!PCSCheck()) // OK
	{
		std::cout << "Not enough data for MOD creation.\n";
		return false;
	}
	else
		return true;
}

#ifdef BRNS

// BRNS-Coupling: For writing spatially resolved reaction rates at the final iteration,
// we need to get the timing information.

double Problem::getCurrentTime()
{
	return current_time;
}

double Problem::getEndTime()
{
	return end_time;
}
#endif // BRNS
