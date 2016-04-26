/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   FEMLib - Object: NUM
   Task:
   Programing:
   11/2004 OK Implementation
   last modified:
**************************************************************************/
#include "rf_num_new.h"

// There is a name conflict between stdio.h and the MPI C++ binding
// with respect to the names SEEK_SET, SEEK_CUR, and SEEK_END.  MPI
// wants these in the MPI namespace, but stdio.h will #define these
// to integer values.  #undef'ing these can cause obscure problems
// with other include files (such as iostream), so we instead use
// #error to indicate a fatal error.  Users can either #undef
// the names before including mpi.h or include mpi.h *before* stdio.h
// or iostream.
#ifdef USE_MPI // WW
#include "mpi.h"
#include "par_ddc.h"
//#undef SEEK_SET
//#undef SEEK_END
//#undef SEEK_CUR
#endif

// C++ STL
#include <cfloat>
#include <cmath>
#include <fstream>
#include <iostream>
#include <list>
#include <string>

#include "makros.h"
#include "memory.h"
#include "display.h"

// FEM-Makros
#include "files0.h"
#include "makros.h"
// GeoSys-GeoLib
// GeoSys-FEMLib
#ifndef NEW_EQS // WW. 06.11.2008
#include "matrix_routines.h"
#endif
#include "StringTools.h"
#include "mathlib.h"
#include "rf_pcs.h"
#include "tools.h"
// GeoSys-MSHLib

using namespace std;

extern std::ios::pos_type GetNextSubKeyword(ifstream* file, string* line, bool* keyword);
extern size_t max_dim; // OK411 todo

//==========================================================================
vector<CNumerics*> num_vector;
/**************************************************************************
   FEMLib-Method:
   Task: constructor
   Programing:
   11/2004 OK Implementation
   10/2005 OK pcs_type_name
   07/2007 OK DEFORMATION
**************************************************************************/
CNumerics::CNumerics(string name)
{
	pcs_type_name = name; // OK
	// GLOBAL
	renumber_method = 0;
	//
	// LS - Linear Solver
	ls_method = 2; // OK41
	ls_max_iterations = 1000;
	ls_error_method = 1;
	ls_error_tolerance = 1e-12;
	ls_theta = 1.0;
	ls_precond = 1;
	ls_storage_method = 2; // OK41
	m_cols = 5; // 06.2010. WW
	ls_extra_arg = ""; // NW
	//
	// NLS - Nonlinear Solver
	nls_method_name = "PICARD";
	nls_method = -1; // Default linear, 0: Picard. 1: Newton. 2:JFNK
	nls_error_method = 1; // JT2012
	nls_max_iterations = 1; // OK
	nls_relaxation = 0.0;
	for (size_t i = 0; i < DOF_NUMBER_MAX; i++) // JT2012
		nls_error_tolerance[i] = -1.0; // JT2012: should not default this. Should always be entered by user!
	//
	// CPL - Coupled processes
	cpl_error_specified = false;
	cpl_master_process = false;
	cpl_process = "INVALID_PROCESS"; // JT2012: do not couple with any process, unless indicated
	cpl_variable = "NONE";
	cpl_variable_JOD = "FLUX";
	cpl_max_iterations = 1; // OK
	cpl_min_iterations = 1; // JT2012
	// Local picard1                                //NW
	local_picard1_tolerance = 1.0e-3;
	local_picard1_max_iterations = 1;
	update_velocity_within_nonlinear = 0;

	for (size_t i = 0; i < DOF_NUMBER_MAX; i++) // JT2012
		cpl_error_tolerance[i] = -1.0; // JT2012: should not default this. Should always be entered by user!
	//
	// ELE
	ele_gauss_points = 2;
	ele_mass_lumping = 0;
	ele_upwind_method = 0; // CB
	ele_upwinding = 0;
	ele_supg_method = 0; // NW
	ele_supg_method_length = 0; // NW
	ele_supg_method_diffusivity = 0; // NW
	fct_method = -1; // NW
	fct_prelimiter_type = 0; // NW
	fct_const_alpha = -1.0; // NW
	//----------------------------------------------------------------------
	// Deformation
	GravityProfile = 0;
	DynamicDamping = NULL; // WW
	if (pcs_type_name.compare("DEFORMATION") == 0)
	{
		ls_method = 2;
		ls_error_method = 2;
		ls_error_tolerance = 1e-12;
		ls_max_iterations = 2000;
		ls_precond = 100;
		ls_storage_method = 4;
	}
	//----------------------------------------------------------------------
	if (pcs_type_name.compare("RICHARDS_FLOW") == 0)
	{
		ele_mass_lumping = 1;
		ele_upwinding = 0.5;
		ls_max_iterations = 2000;
		ls_error_method = 2;
		ls_error_tolerance = 1e-10;
		ls_precond = 4;
		ls_storage_method = 4;
		nls_max_iterations = 25;
	}

#ifdef USE_PETSC
	lsover_name = "bcgs";
	pres_name = "bjacobi";
#endif
}

/**************************************************************************
   FEMLib-Method:
   Task: deconstructor
   Programing:
   11/2004 OK Implementation
**************************************************************************/
CNumerics::~CNumerics(void)
{
	if (DynamicDamping)
		delete[] DynamicDamping;
	DynamicDamping = NULL;
}

/**************************************************************************
   FEMLib-Method:
   Task: OBJ read function
   Programing:
   11/2004 OK Implementation
**************************************************************************/
bool NUMRead(string file_base_name)
{
	//----------------------------------------------------------------------
	// OK  NUMDelete();
	//----------------------------------------------------------------------
	CNumerics* m_num = NULL;
	char line[MAX_ZEILE];
	bool overall_coupling_exists = false; // JT
	string sub_line;
	string line_string;
	ios::pos_type position;
	//========================================================================
	// File handling
	string num_file_name = file_base_name + NUM_FILE_EXTENSION;
	ifstream num_file(num_file_name.data(), ios::in);
	if (!num_file.good())
		return false;
	num_file.seekg(0L, ios::beg);
	//========================================================================
	// Keyword loop
	cout << "NUMRead" << "\n";
	int max_num_integration_pnts = 0;
	while (!num_file.eof())
	{
		num_file.getline(line, MAX_ZEILE);
		line_string = line;
		if (line_string.find("#STOP") != string::npos)
		{
			// Unify the number of integration points.
			if (max_num_integration_pnts > 3)
				max_num_integration_pnts = 3;
			for (std::size_t i = 0; i < num_vector.size(); i++)
			{
				num_vector[i]->ele_gauss_points = max_num_integration_pnts;
			}

			return true;
		}
		//
		if (line_string.find("$OVERALL_COUPLING") != string::npos)
		{
			overall_coupling_exists = true; // JT: for error checking
		}
		//----------------------------------------------------------------------
		// keyword found
		if (line_string.find("#NUMERICS") != string::npos)
		{
			m_num = new CNumerics("default");
			position = m_num->Read(&num_file);

			max_num_integration_pnts = std::max(max_num_integration_pnts, m_num->ele_gauss_points);

			num_vector.push_back(m_num);
			num_file.seekg(position, ios::beg);
			m_num->NumConfigure(overall_coupling_exists); // JT2012
		} // keyword found
	} // eof

	return true;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   05/2005 WW Implementation
**************************************************************************/
bool CNumerics::CheckDynamic()
{
	if (DynamicDamping)
		return true;
	else
		return false;
}

/**************************************************************************
FEMLib-Method:
Task: After a read of each m_num, configure any defaults here
Programing:
05/2011 JT Implementation
**************************************************************************/
void CNumerics::NumConfigure(bool overall_coupling_exists)
{
	// Overall coupling check
	if (overall_coupling_exists && !cpl_error_specified)
	{
		if (this->nls_method < 0)
		{
			std::cout << "ERROR in NUMRead. Overall coupling requested, but ";
			std::cout << this->pcs_type_name << " was not\n";
			std::cout << "supplied with coupling tolerance. See $COUPLING_CONTROL keyword to enter this.\n";
			exit(1);
		}
		else
		{ // Can take the non-linear tolerance as an emergency backup
			std::cout << "WARNING in NUMRead. Overall coupling requested, but ";
			std::cout << this->pcs_type_name << " was not\n";
			std::cout << "supplied with coupling tolerance. Adopting 10*non_linear_tolerance.\n";
			setCouplingErrorMethod(getNonLinearErrorMethod());
			for (size_t i = 0; i < DOF_NUMBER_MAX; i++)
			{
				cpl_error_tolerance[i] = 10.0 * nls_error_tolerance[i];
			}
			cpl_error_specified = true;
		}
	}
	//
	// Check master processes
	if (cpl_master_process && !cpl_error_specified)
	{
		std::cout << "ERROR in NUMRead. Process coupling requested, but ";
		std::cout << this->pcs_type_name << " was not\n";
		std::cout << "supplied with coupling tolerance. See $COUPLING_CONTROL keyword to enter this.\n";
		exit(1);
	}
	//
	// We are ok. Now check the tolerances.
	if (this->nls_method < 0)
	{ // linear solution
		if (cpl_error_specified)
		{ // A coupling error was entered. Adopt this for error calculations.
			for (size_t i = 0; i < DOF_NUMBER_MAX; i++)
			{
				nls_error_tolerance[i] = cpl_error_tolerance[i];
			}
			setNonLinearErrorMethod(getCouplingErrorMethod());
		}
		else
		{ // We have no error tolerances for non-linear or coupled simulations. Force some defaults.
			setNonLinearErrorMethod(FiniteElement::LMAX);
			setCouplingErrorMethod(FiniteElement::LMAX);
			nls_error_tolerance[0] = cpl_error_tolerance[0] = 1.0;
		}
	}
	// Default CPL error method to NLS method. Just so error is not checked twice
	if (!cpl_error_specified)
	{
		setCouplingErrorMethod(getNonLinearErrorMethod());
	}
	//
	// Default all NLS tolerances to the previous DOF, if they were not entered.
	for (size_t i = 1; i < DOF_NUMBER_MAX; i++)
	{
		if (nls_error_tolerance[i] < 0.0)
			nls_error_tolerance[i] = nls_error_tolerance[i - 1];
	}
	//
	// Default all CPL tolerances to the previous DOF, if they were not entered.
	for (size_t i = 1; i < DOF_NUMBER_MAX; i++)
	{
		if (cpl_error_tolerance[i] < 0.0)
			cpl_error_tolerance[i] = cpl_error_tolerance[i - 1];
	}
}

/**************************************************************************
   FEMLib-Method:
   Task: OBJ read function
   Programing:
   11/2004 OK Implementation
**************************************************************************/
ios::pos_type CNumerics::Read(ifstream* num_file)
{
	string line_string;
	std::string error_method_name;
	std::string coupling_target;
	bool new_keyword = false;
	bool new_subkeyword = false;
	ios::pos_type position;
	ios::pos_type position_subkeyword;
	std::stringstream line;
	//========================================================================
	// Schleife ueber alle Phasen bzw. Komponenten
	while (!new_keyword)
	{
		if (new_subkeyword)
			num_file->seekg(position, ios::beg);
		new_subkeyword = false;
		position = GetNextSubKeyword(num_file, &line_string, &new_keyword);
		if (new_keyword)
			return position;
		//....................................................................
		// subkeyword found
		if (line_string.find("$PCS_TYPE") != string::npos)
		{
			line.str(GetLineFromFile1(num_file));
			line >> pcs_type_name;
			line.clear();
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$RENUMBER") != string::npos)
		{
			line.str(GetLineFromFile1(num_file));
			line >> renumber_method;
			if (renumber_method == 2)
				line >> renumber_parameter;
			line.clear();
			continue;
		}
		//....................................................................
		// JT->WW: Local tolerance previously found in $NON_LINEAR_SOLVER for NEWTON. Moved here for now.
		if (line_string.find("$PLASTICITY_TOLERANCE") != string::npos)
		{
			line.str(GetLineFromFile1(num_file));
			line >> nls_plasticity_local_tolerance;
		}
		//....................................................................
		// subkeyword found ($NON_LINEAR_ITERATION  -or-  $NON_LINEAR_ITERATIONS)
		if (line_string.find("$NON_LINEAR_ITERATION") != string::npos)
		{
			// JT:	in >> nls_method_name
			//		in >> error_method_name
			//		in >> max iter
			//		in >> relaxation
			//		in >> tolerance[1:dof]
			//
			line.str(GetLineFromFile1(num_file));
			line >> nls_method_name;
			line >> error_method_name;
			line >> nls_max_iterations;
			line >> nls_relaxation;
			//
			setNonLinearErrorMethod(FiniteElement::convertErrorMethod(error_method_name));
			switch (getNonLinearErrorMethod())
			{
				case FiniteElement::ENORM: // only 1 tolerance required
					line >> nls_error_tolerance[0];
					break;
				//
				case FiniteElement::ERNORM: // only 1 tolerance required
					line >> nls_error_tolerance[0];
					break;
				//
				case FiniteElement::EVNORM: // 1 tolerance for each primary variable (for Deformation, only 1 tolerance
				                            // required. Applies to x,y,z)
					for (int i = 0; i < DOF_NUMBER_MAX; i++)
						line >> nls_error_tolerance[i];
					break;
				//
				case FiniteElement::LMAX: // 1 tolerance for each primary variable (for Deformation, only 1 tolerance
				                          // required. Applies to x,y,z)
					for (int i = 0; i < DOF_NUMBER_MAX; i++)
						line >> nls_error_tolerance[i];
					break;
				//
				case FiniteElement::BNORM: // only 1 tolerance required
					line >> nls_error_tolerance[0];
					break;
				//
				default:
					ScreenMessage("ERROR in NUMRead. Invalid non-linear iteration error method selected.\n");
					exit(1);
					break;
			}

			nls_method = 0;
			if (nls_method_name.find("NEWTON") != string::npos)
				nls_method = 1;
			else if (nls_method_name.find("JFNK") != string::npos) //  Jacobian free Newton-Krylov method
				nls_method = 2;
			//
			line.clear();
			continue;
		}
		else if (line_string.find("$NON_LINEAR_SOLVER") != string::npos)
		{
			ScreenMessage(" --\n Using old $NON_LINEAR_SOLVER keyword.\n");
			ScreenMessage(" Eventually this will be obsolete. Consider switching to\n");
			ScreenMessage(" $NON_LINEAR_ITERATIONS for better results and greater flexibility.\n");
			ScreenMessage(" --\n");
			//
			// JT:	in >> method_name
			//		in >> tolerance
			//		if(NEWTON) in >> tolerance_local
			//		in >> max iter
			//		in >> relaxation
			//
			//
			line.str(GetLineFromFile1(num_file));
			line >> nls_method_name;
			//
			nls_method = 0;
			if (nls_method_name.find("NEWTON") != string::npos)
				nls_method = 1;
			else if (nls_method_name.find("JFNK") != string::npos) //  Jacobian free Newton-Krylov method
				nls_method = 2;
			//
			if (nls_method > 0)
			{
				line >> nls_error_tolerance[0];
				line >> nls_plasticity_local_tolerance;
				error_method_name = "BNORM"; // JT: this is hardwired in old version
			}
			else
			{
				line >> nls_error_tolerance[0];
				error_method_name = "LMAX"; // JT: this is hardwired in old version
			}
			setNonLinearErrorMethod(FiniteElement::convertErrorMethod(error_method_name));
			//
			line >> nls_max_iterations;
			line >> nls_relaxation;
			line.clear();
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$LINEAR_SOLVER") != string::npos)
		{
			std::string str_buf = GetLineFromFile1(num_file); // WW
			line.str(str_buf);
#ifdef USE_PETSC
			if (str_buf.find("petsc") != string::npos) // 03.2012. WW
			{
				line >> str_buf >> lsover_name >> pres_name >> ls_error_tolerance >> ls_max_iterations >> ls_theta;
			}
			else
#endif
			{
				line >> ls_method;
				line >> ls_error_method;
				line >> ls_error_tolerance;
				line >> ls_max_iterations;
				line >> ls_theta;
				line >> ls_precond;
				line >> ls_storage_method;
				/// For GMRES. 06.2010. WW
				if (ls_method == 13)
					line >> m_cols;
			}
			line.clear();
			continue;
		}
		//....................................................................
		// JT subkeyword found
		if (line_string.find("$COUPLING_ITERATIONS") != string::npos)
		{
			ScreenMessage("$COUPLING_ITERATIONS keyword obsolete.\n");
			ScreenMessage("Use $COUPLING_CONTROL and $COUPLED_PROCESS for process couplings.\n");
			exit(1);
		}
		//....................................................................
		// JT subkeyword found
		if (line_string.find("$COUPLING_CONTROL")
		    != string::npos) // JT: For this process, how is coupling error handled?
		{
			// JT:	in >> error_method_name
			//		in >> tolerance[1:dof]
			//
			line.str(GetLineFromFile1(num_file));
			line >> error_method_name;
			//
			cpl_error_specified = true;
			setCouplingErrorMethod(FiniteElement::convertErrorMethod(error_method_name));
			switch (getCouplingErrorMethod())
			{
				case FiniteElement::ENORM: // only 1 tolerance required
					line >> cpl_error_tolerance[0];
					break;
				//
				case FiniteElement::ERNORM: // only 1 tolerance required
					line >> cpl_error_tolerance[0];
					break;
				//
				case FiniteElement::EVNORM: // 1 tolerance for each primary variable (for Deformation, only 1 tolerance
				                            // required. Applies to x,y,z)
					for (int i = 0; i < DOF_NUMBER_MAX; i++)
						line >> cpl_error_tolerance[i];
					break;
				//
				case FiniteElement::LMAX: // 1 tolerance for each primary variable (for Deformation, only 1 tolerance
				                          // required. Applies to x,y,z)
					for (int i = 0; i < DOF_NUMBER_MAX; i++)
						line >> cpl_error_tolerance[i];
					break;
				//
				case FiniteElement::BNORM:
					ScreenMessage("ERROR in NUMRead. BNORM not configured for process couplings.\n");
					ScreenMessage("We suggest ENORM as a valid companion for NEWTON couplings.\n");
					exit(1);
					break;
				//
				default:
					ScreenMessage("ERROR in NUMRead. Invalid coupling error method selected.\n");
					exit(1);
					break;
			}
			//
			line.clear();
			continue;
		}
		//....................................................................
		// JT subkeyword found
		if (line_string.find("$COUPLED_PROCESS")
		    != string::npos) // JT: Is this process coupled to another process in an inner loop?
		{
			// in >> process name >> min iter >> max iter
			//
			line.str(GetLineFromFile1(num_file));
			line >> coupling_target; // name of coupled process -OR- process variable
			line >> cpl_min_iterations;
			line >> cpl_max_iterations;
			//
			cpl_master_process = true;
			//
			// Is coupling through a process name or a primary variable?
			if (FiniteElement::convertPrimaryVariable(coupling_target) != FiniteElement::INVALID_PV)
			{ // Then a valid process VARIABLE is entered. Use this.
				cpl_variable = coupling_target;
			}
			else if (PCSGet(coupling_target))
			{ // Then a valid process is entered
				cpl_process = coupling_target;
			}
			else
			{
				ScreenMessage(
				    "WARNING. $COUPLED_PROCESS keyword encountered, but a valid process OR primary variable was not "
				    "found.\n");
				cpl_master_process = false;
			}
			//
			line.clear();
			continue;
		}
		//....................................................................
		if (line_string.find("$EXTERNAL_SOLVER_OPTION") != string::npos) // subkeyword found
		{
			ls_extra_arg = GetLineFromFile1(num_file);
			trim(ls_extra_arg);
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$ELE_GAUSS_POINTS") != string::npos)
		{
			line.str(GetLineFromFile1(num_file));
			line >> ele_gauss_points; // probably element-type-wise
			line.clear();
			continue;
		}
		// subkeyword found
		if (line_string.find("$ELE_MASS_LUMPING") != string::npos)
		{
			line.str(GetLineFromFile1(num_file));
			line >> ele_mass_lumping;
			line.clear();
			cout << "-> Mass lumping selected for " << pcs_type_name << "\n"; // JOD 2014-11-10
			continue;
		}
		// subkeyword found
		if (line_string.find("$ELE_UPWINDING") != string::npos)
		{
			line.str(GetLineFromFile1(num_file));
			// CB now read also upwinding method
			line >> ele_upwinding >> ele_upwind_method;
			line.clear();
			continue;
		}
		// subkeyword found
		if (line_string.find("$ELE_SUPG") != string::npos)
		{
			line.str(GetLineFromFile1(num_file));
			// NW
			line >> ele_supg_method >> ele_supg_method_length >> ele_supg_method_diffusivity;
			line.clear();
			cout << "->SUPG method is selected."
			     << "\n";
			continue;
		}
		// subkeyword found
		if (line_string.find("$GRAVITY_PROFILE") != string::npos)
		{
			line.str(GetLineFromFile1(num_file)); // WW
			line >> GravityProfile;
			line.clear();
			continue;
		}
		// subkeyword found
		if (line_string.find("$DYNAMIC_DAMPING") != string::npos)
		{
			line.str(GetLineFromFile1(num_file)); // WW
			DynamicDamping = new double[3];
			// Default
			DynamicDamping[0] = 0.515;
			DynamicDamping[1] = 0.51;
			DynamicDamping[2] = 0.51;
			line >> DynamicDamping[0] >> DynamicDamping[1] >> DynamicDamping[2];
			line.clear();
			continue;
		}
		// Flux corrected transport by Kuzmin (2009)
		if (line_string.find("$LOCAL_PICARD1") != string::npos) // NW
		{
			line.str(GetLineFromFile1(num_file));
			line >> local_picard1_max_iterations;
			line >> local_picard1_tolerance;
			line.clear();
			continue;
		}
		if (line_string.find("$NON_LINEAR_UPDATE_VELOCITY") != string::npos) // NW
		{
			line.str(GetLineFromFile1(num_file));
			line >> update_velocity_within_nonlinear;
			line.clear();
			continue;
		}
		// NW
		if (line_string.find("$FEM_FCT") != string::npos)
		{
			line.str(GetLineFromFile1(num_file));
			line >> fct_method; // 1: linearized FCT
			line >> fct_prelimiter_type; // 0: just cancel, 1: minmod, 2: superbee
			line >> fct_const_alpha; //-1: off, [0.0,1.0] 0: Upwind, 1: Galerkin
			line.clear();
			cout << "->FEM_FCT method is selected."
			     << "\n";
			continue;
		}

		//....................................................................
		/*
		    if(line_string.find("$TIME_STEPS")!=string::npos) { // subkeyword found
		      while((!new_keyword)||(!new_subkeyword)||(!num_file->eof())){
		        position = num_file->tellg();
		        line_string = GetLineFromFile1(num_file);
		        if(line_string.find("#")!=string::npos){
		          return position;
		        }
		        if(line_string.find("$")!=string::npos){
		          new_subkeyword = true;
		          break;
		   }
		   line.str(line_string);
		   line >> no_time_steps;
		   line >> time_step_length;
		   for(i=0;i<no_time_steps;i++)
		   time_step_vector.push_back(time_step_length);
		   line.clear();
		   }
		   }
		 */
		//....................................................................
	}
	return position;
}

/**************************************************************************
   FEMLib-Method:
   Task: master write function
   Programing:
   11/2004 OK Implementation
   last modification:
**************************************************************************/
void NUMWrite(string base_file_name)
{
	CNumerics* m_num = NULL;
	string sub_line;
	string line_string;
	//========================================================================
	// File handling
	string num_file_name = base_file_name + NUM_FILE_EXTENSION;
	fstream num_file(num_file_name.data(), ios::trunc | ios::out);
	num_file.setf(ios::scientific, ios::floatfield);
	num_file.precision(12);
	if (!num_file.good())
		return;
	num_file.seekg(0L, ios::beg);
	//========================================================================
	num_file << "GeoSys-NUM: Numerics ------------------------------------------------\n";
	//========================================================================
	// OUT vector
	int num_vector_size = (int)num_vector.size();
	int i;
	for (i = 0; i < num_vector_size; i++)
	{
		m_num = num_vector[i];
		m_num->Write(&num_file);
	}
	num_file << "#STOP";
	num_file.close();
}

/**************************************************************************
   FEMLib-Method:
   Task: write function
   Programing:
   11/2004 OK Implementation
   last modification:
**************************************************************************/
void CNumerics::Write(fstream* num_file)
{
	// KEYWORD
	*num_file << "#NUMERICS"
	          << "\n";
	//--------------------------------------------------------------------
	/*OK
	   *num_file << " $METHOD" << "\n";
	   *num_file << method_name << "\n";
	   if(method_name.find("LAGRANGE")!=string::npos){
	   *num_file << lag_quality << " " << lag_max_steps << " " << lag_local_eps << " ";
	   *num_file << lag_time_weighting << " " << lag_min_weight << " ";
	   *num_file << lag_use_matrix << " " << lag_vel_method;
	   *num_file << "\n";
	   }
	 */
	//--------------------------------------------------------------------
	*num_file << " $PCS_TYPE"
	          << "\n";
	*num_file << "  " << pcs_type_name << "\n";
	//--------------------------------------------------------------------
	*num_file << " $NON_LINEAR_SOLVER"
	          << "\n";
	*num_file << "  " << nls_method_name;
	*num_file << " " << nls_error_tolerance;
	*num_file << " " << nls_max_iterations;
	*num_file << " " << nls_relaxation;
	*num_file << "\n";
	//--------------------------------------------------------------------
	*num_file << " $LINEAR_SOLVER"
	          << "\n";
	*num_file << "  " << ls_method;
	*num_file << " " << ls_error_method;
	*num_file << " " << ls_error_tolerance;
	*num_file << " " << ls_max_iterations;
	*num_file << " " << ls_theta;
	*num_file << " " << ls_precond;
	*num_file << " " << ls_storage_method;
	*num_file << "\n";
	//--------------------------------------------------------------------
	*num_file << " $ELE_GAUSS_POINTS"
	          << "\n";
	*num_file << "  " << ele_gauss_points;
	*num_file << "\n";
	//--------------------------------------------------------------------
	*num_file << " $ELE_MASS_LUMPING"
	          << "\n";
	*num_file << "  " << ele_mass_lumping;
	*num_file << "\n";
	//--------------------------------------------------------------------
	*num_file << " $ELE_UPWINDING"
	          << "\n";
	*num_file << "  " << ele_upwinding;
	*num_file << "\n";
	//--------------------------------------------------------------------
}

//////////////////////////////////////////////////////////////////////////
// LINEAR_SOLVER
//////////////////////////////////////////////////////////////////////////
//#ifndef NEW_EQS                                   //WW. 06.11.2008
#if !defined(NEW_EQS) && !defined(USE_PETSC) // && !defined(other parallel solver) //WW. 04.10.2012

/*************************************************************************
   ROCKFLOW - Funktion: SetZeroLinearSolver

   Aufgabe:
   Setzt Matrix, rechte Seite und oesungsvektor auf NULL.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E LINEAR_SOLVER *ls: Zeiger auf eine Instanz vom Typ LINEAR_SOLVER.

   Ergebnis:
   - Adresse des LS's im Erfolgsfall, ansonsten NULL-Zeiger -

   Programmaenderungen:
   02/1999     AH         Erste Version

*************************************************************************/
LINEAR_SOLVER* SetZeroLinearSolver(LINEAR_SOLVER* ls)
{
	/* Initialisieren der Felder */
	if (!ls)
		return NULL;
	SetLinearSolver(ls);
	if (ls->matrix)
		MXInitMatrix();
	if (ls->b)
		MNulleVec(ls->b, ls->dim);
	if (ls->x)
		MNulleVec(ls->x, ls->dim);
	/*if (ls->xx) MNulleVec(ls->xx,ls->dim); */
	/*if (ls->memory) MNulleVec(ls->memory,ls->dim); */

	return ls;
}

/*************************************************************************
   ROCKFLOW - Funktion: SetLinearSolver

   Aufgabe:
   Konstruktor fuer LINEAR_SOLVER

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E LINEAR_SOLVER *ls: Zeiger auf eine Instanz vom Typ LINEAR_SOLVER.

   Ergebnis:
   - Adresse des LS's im Erfolgsfall, ansonsten NULL-Zeiger -

   Programmaenderungen:
   02/1999     AH         Erste Version
   7/2000     CT         Implizite Kenntnis ueber Matrixspeichermodell

*************************************************************************/
LINEAR_SOLVER* SetLinearSolver(LINEAR_SOLVER* ls)
{
	if (!ls)
		return NULL;
	if (ls->matrix)
		MXSetMatrixPointer(ls->matrix);

	return ls;
}

/*************************************************************************
   ROCKFLOW - Funktion: GetUnknownVectorDimensionLinearSolver

   Aufgabe:
   Liefert die Dimension der Unbekannten im Vektors des lin. Loesers.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E LINEAR_SOLVER *ls: Zeiger auf eine Instanz vom Typ LINEAR_SOLVER.

   Ergebnis:
   - s. o. -

   Programmaenderungen:
   08/2000     MK         Erste Version

*************************************************************************/
int GetUnknownVectorDimensionLinearSolver(LINEAR_SOLVER* ls)
{
	if (!ls)
		return -1;
	else
		return ls->unknown_vector_dimension;
}

/*************************************************************************
   ROCKFLOW - Funktion: DestroyLinearSolver

   Aufgabe:
   Destruktor fuer LINEAR_SOLVER

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E LINEAR_SOLVER *ls: Zeiger auf eine Instanz vom Typ LINEAR_SOLVER.

   Ergebnis:
   - void -

   Programmaenderungen:
   02/1999     AH         Erste Version
   08/1999     AH         Erweiterung memory  (neq)
   10/1999     AH         Systemzeit
   07/2000     C.Thorenz  Bugfix fuer MXDestroy

*************************************************************************/
LINEAR_SOLVER* DestroyLinearSolver(LINEAR_SOLVER* ls)
{
	if (!ls)
		return NULL;
	SetLinearSolver(ls);

	if (ls->system_time_name)
		ls->system_time_name = (char*)Free(ls->system_time_name);
	if (ls->system_time_assemble_function_name)
		ls->system_time_assemble_function_name = (char*)Free(ls->system_time_assemble_function_name);

	// OK    if (ls->name)
	// OK        ls->name = (char *) Free(ls->name);
	if (ls->b)
		ls->b = (double*)Free(ls->b);
	if (ls->x)
		ls->x = (double*)Free(ls->x);
	/* //WW
	   if (ls->bc_name)
	    ls->bc_name = (char *) Free(ls->bc_name);
	   if (ls->bc_name2)
	    ls->bc_name2 = (char *) Free(ls->bc_name2);
	   if (ls->bc_name3)
	    ls->bc_name3 = (char *) Free(ls->bc_name3);
	   if (ls->sousin_name1)
	    ls->sousin_name1 = (char *) Free(ls->sousin_name1);
	   if (ls->sousin_name2)
	    ls->sousin_name2 = (char *) Free(ls->sousin_name2);
	 */
	if (ls->lsp_name)
		ls->lsp_name = (char*)Free(ls->lsp_name);
	if (ls->r)
		ls->r = (double*)Free(ls->r);
	if (ls->xx)
		ls->xx = (double*)Free(ls->xx);
	if (ls->memory)
		ls->memory = (double*)Free(ls->memory);
	if (ls->new_memory)
		ls->new_memory = (double**)DestroyMemoryLinearSolver(ls);
	if (ls->matrix)
	{
		MXSetMatrixPointer(ls->matrix);
		ls->matrix = MXDestroyMatrix();
	}
	/* Multiple unknowns (Dim) / multiple solvers (MultSolvers) */
	if (ls->name_ls)
		ls->name_ls = (char*)Free(ls->name_ls);

	if (ls->name_group_ls)
		ls->name_group_ls = (char*)Free(ls->name_group_ls);
	if (ls)
		ls = (LINEAR_SOLVER*)Free(ls);
	return NULL;
}

/*************************************************************************
   ROCKFLOW - Funktion: DestroyMemoryLinearSolver

   Aufgabe:
   Extra-Speicher-Destruktor einer Instanz vom Typ LINEAR_SOLVER.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E LINEAR_SOLVER *ls: Zeiger auf eine Instanz vom Typ LINEAR_SOLVER.

   Ergebnis:
   - Ungleich NULL im Erfolgsfall -

   Programmaenderungen:
   03/2000    AH      Erste Version

*************************************************************************/
LINEAR_SOLVER* DestroyMemoryLinearSolver(LINEAR_SOLVER* ls)
{
	int i;

	if (!ls)
		return NULL;
	SetLinearSolver(ls);

	for (i = 0; i < ls->memory_number; i++)
		if (ls->new_memory[i])
			ls->new_memory[i] = (double*)Free(ls->new_memory[i]);
	if (ls->new_memory)
		ls->new_memory = (double**)Free(ls->new_memory);

	return ls;
}

/*************************************************************************
   ROCKFLOW - Funktion: InitLinearSolver

   Aufgabe:
   Konstruktor fuer LINEAR_SOLVER

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E LINEAR_SOLVER *ls: Zeiger auf eine Instanz vom Typ LINEAR_SOLVER.

   Ergebnis:
   - Adresse des LS's im Erfolgsfall, ansonsten NULL-Zeiger -

   Programmaenderungen:
   02/1999     AH         Erste Version
   08/1999     AH         Erweiterung memory  (neq)
   7/2000     CT         Implizite Kenntnis ueber Matrixspeichermodell
   Fehlerbereinigung

*************************************************************************/
LINEAR_SOLVER* InitLinearSolver(LINEAR_SOLVER* ls)
{
	if (!ls)
		return NULL;

	//#ifdef USE_MPI //WW
	//    ls->matrix = NULL;
	//#else
	if (ls->matrix)
	{
		MXSetMatrixPointer(ls->matrix);
		ls->matrix = MXDestroyMatrix();
	}
	ls->matrix = MXSetupMatrix(ls->dim, ls->store, 0l);
	//#endif
	if (ls->dim == 0)
		return ls;

	if (ls->b)
		ls->b = (double*)Free(ls->b);
	ls->b = (double*)Malloc(ls->dim * sizeof(double));
	if (ls->x)
		ls->x = (double*)Free(ls->x);
	ls->x = (double*)Malloc(ls->dim * sizeof(double));
	if (ls->r)
		ls->r = (double*)Free(ls->r);
	ls->r = (double*)Malloc(ls->dim * sizeof(double));

	if (ls->xx)
		ls->xx = (double*)Free(ls->xx);
	ls->xx = (double*)Malloc(ls->dim * sizeof(double));
	if (ls->memory)
		ls->memory = (double*)Free(ls->memory);
	ls->memory = (double*)Malloc(ls->dim * sizeof(double));

	return ls;
}

#ifdef USE_MPI // WW
/////////////////////////////////////////////////////////////////////
/**************************************************************************
   FEMLib-Method:
   Task:
      Based on LINEAR_SOLVER *InitLinearSolver(LINEAR_SOLVER * ls)
   07/2006 WW
   last modified:
**************************************************************************/

LINEAR_SOLVER* InitVectorLinearSolver(LINEAR_SOLVER* ls)
{
	if (!ls)
		return NULL;
	if (ls->dim == 0)
		return ls;

	if (ls->b)
		ls->b = (double*)Free(ls->b);
	ls->b = (double*)Malloc(ls->dim * sizeof(double));
	if (ls->x)
		ls->x = (double*)Free(ls->x);
	ls->x = (double*)Malloc(ls->dim * sizeof(double));
	if (ls->r)
		ls->r = (double*)Free(ls->r);
	ls->r = (double*)Malloc(ls->dim * sizeof(double));

	if (ls->xx)
		ls->xx = (double*)Free(ls->xx);
	ls->xx = (double*)Malloc(ls->dim * sizeof(double));
	if (ls->memory)
		ls->memory = (double*)Free(ls->memory);
	ls->memory = (double*)Malloc(ls->dim * sizeof(double));

	return ls;
}
#endif
/////////////////////////////////////////////////////////////////////

/*************************************************************************
   ROCKFLOW - Funktion: ConfigSolverProperties

   Aufgabe:
   Konfiguration von Loeser-Parametern

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)

   Ergebnis:
   - void -

   Programmaenderungen:
   12/1999   OK   Implementierung

*************************************************************************/
void ConfigSolverProperties(void)
{
	/* Speichertechnik automatisch optimieren */
	switch (max_dim)
	{
		case 0:
			sp2_start = 3;
			sp2_inc = 1;
			break;
		case 1:
			sp2_start = 11;
			sp2_inc = 2;
			break;
		case 2:
			sp2_start = 33;
			sp2_inc = 6;
			break;
	}
	/* Umnummerierer-Methode auswaehlen */
	// WW ConfigRenumberProperties();
	/* umnummerierer muss definiert sein ! */
}

/*************************************************************************
   ROCKFLOW - Funktion: InitializeLinearSolver

   Aufgabe:
   Loesungsvektor vorbelegen.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E LINEAR_SOLVER *ls: Zeiger auf eine Instanz vom Typ LINEAR_SOLVER.

   Ergebnis:
   - Adresse des LS's im Erfolgsfall, ansonsten NULL-Zeiger -

   Programmaenderungen:
   02/1999     AH         Erste Version
   08/1999     AH         Erweiterung memory  (neq)
   09/2004      Remove BC/Sink incoorperation
*************************************************************************/
// WW
void SetLinearSolverType(LINEAR_SOLVER* ls, CNumerics* m_num)
{
	if (!ls)
		return;
	/* Gruppeneingenschaften setzen */
	/*
	   SetPropertiesLinearSolver(ls, prop_name);
	   // Benutzer-Eigenschaften
	   lsp = GetTLinearSolverProperties(ls->lsp_name, lsp, 0);
	   if (lsp) {
	    ls->lsp = lsp;
	    ls->store = get_lsp_store(ls->lsp);
	   }
	   // Entwickler-Eigenschaften
	   if (ls->lsp_name && !ls->lsp) {
	    sprintf(string, "DEFAULT_");
	   strcat(string, ls->lsp_name);
	   lsp = GetLinearSolverProperties(string);
	   if (lsp) {
	   ls->lsp = lsp;
	   ls->store = get_lsp_store(ls->lsp);
	   }
	   }
	   // Master Default-Eigenschaften
	   if (!ls->lsp) {
	   ls->lsp = GetLinearSolverProperties(MASTER_DEFAULT_LINEAR_SOLVER_PROPERTIES);
	   if (ls->lsp) {
	   ls->store = get_lsp_store(ls->lsp);
	   }
	   }
	 */
	ls->store = m_num->ls_storage_method;
	switch (m_num->ls_method) // OK (ls->lsp->type){
	{
		case 1:
			ls->LinearSolver = SpGauss;
			break;
		case 2:
#ifndef USE_MPI
			ls->LinearSolver = SpBICGSTAB;
#endif
			break;
		case 3:
			ls->LinearSolver = SpBICG;
			break;
		case 4:
			ls->LinearSolver = SpQMRCGSTAB;
			break;
		case 5:
			ls->LinearSolver = SpCG;
			break;
		case 6:
			ls->LinearSolver = SpCGNR;
			break;
		case 7:
			ls->LinearSolver = SpCGS;
			break;
		case 8:
			ls->LinearSolver = SpRichardson;
			break;
		case 9:
			ls->LinearSolver = SpJOR;
			break;
		case 10:
			ls->LinearSolver = SpSOR;
			break;
		case 11:
			ls->LinearSolver = SpAMG1R5;
			break;
		case 12:
			ls->LinearSolver = SpUMF;
			break;
		default:
			cout << "***ERROR in SetLinearSolverType(): Specified linear solver type (" << m_num->ls_method
			     << ") is not supported. "
			     << "\n";
			exit(1);
	}
}

// WW
LINEAR_SOLVER* InitializeLinearSolver(LINEAR_SOLVER* ls, CNumerics* m_num)
{
	SetLinearSolverType(ls, m_num);
	//#ifdef USE_MPI                                 //WW
	//	InitVectorLinearSolver(ls);
	//#else
	InitLinearSolver(ls);
	//#endif
	/* Internen Speicher allokieren */
	InitMemoryLinearSolver(ls, ls->memory_number);
	/* Speicher initialisieren */
	SetMemoryZeroLinearSolver(ls);
	return ls;
}

/*************************************************************************
   ROCKFLOW - Funktion: InitMemoryLinearSolver

   Aufgabe:
   Extra-Speicher-Konstruktor einer Instanz vom Typ LINEAR_SOLVER.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E LINEAR_SOLVER *ls: Zeiger auf eine Instanz vom Typ LINEAR_SOLVER.

   Ergebnis:
   - Ungleich NULL im Erfolgsfall -

   Programmaenderungen:
   03/2000    AH      Erste Version

*************************************************************************/
LINEAR_SOLVER* InitMemoryLinearSolver(LINEAR_SOLVER* ls, int memory_number)
{
	int i;

	if (!ls || memory_number <= 0 || ls->dim <= 0)
		return NULL;
	SetLinearSolver(ls);

	if (ls->new_memory)
		DestroyMemoryLinearSolver(ls);

	ls->new_memory = (double**)Malloc(memory_number * sizeof(double*));
	if (ls->new_memory == NULL)
		return NULL;
	for (i = 0; i < memory_number; i++)
		ls->new_memory[i] = NULL;

	for (i = 0; i < memory_number; i++)
	{
		ls->new_memory[i] = (double*)Malloc(ls->dim * sizeof(double));
		if (ls->new_memory[i] == NULL)
			return DestroyMemoryLinearSolver(ls);
	}

	ls->memory_number = memory_number;
	return ls;
}

/*************************************************************************
   ROCKFLOW - Funktion: SetMemoryZeroLinearSolver

   Aufgabe:
   Setzt Speichervektor auf NULL.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E LINEAR_SOLVER *ls: Zeiger auf eine Instanz vom Typ LINEAR_SOLVER.

   Ergebnis:
   - Adresse des LS's im Erfolgsfall, ansonsten NULL-Zeiger -

   Programmaenderungen:
   08/1999     AH         Erste Version (ah neq)

*************************************************************************/
LINEAR_SOLVER* SetMemoryZeroLinearSolver(LINEAR_SOLVER* ls)
{
	if (!ls)
		return NULL;
	SetLinearSolver(ls);
	if (ls->xx)
		MNulleVec(ls->xx, ls->dim);

	return ls;
}

/*************************************************************************
   ROCKFLOW - Funktion: CreateLinearSolver

   Aufgabe:
   Konstruktor fuer LINEAR_SOLVER

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E char *name: Zeiger auf den Namen des LS's.

   Ergebnis:
   - Adresse des LS's im Erfolgsfall, ansonsten NULL-Zeiger -

   Programmaenderungen:
   02/1999     AH         Erste Version
   08/1999     AH         Erweiterung memory  (neq)
   10/1999     AH         Systemzeit
   10/1999     AH         Systemzeit
   05/2000     OK         Erweiterung fuer vektorielle Groessen
   08/2000     MK         Erweiterung fuer variable Dimension der Unbekannten
   04/2003     MK         Erweiterung fuer variable Anzahl der Loeserobjekte

*************************************************************************/
LINEAR_SOLVER* CreateLinearSolver(long store, long dim)
{
	LINEAR_SOLVER* ls;

	ls = (LINEAR_SOLVER*)Malloc(sizeof(LINEAR_SOLVER));
	if (ls == NULL)
		return NULL;

	// OK    ls->name = NULL;
	ls->store = store;
	ls->dim = dim;
	ls->matrix = NULL;
	ls->b = NULL;
	ls->x = NULL;
	ls->lsp = NULL;
	ls->unknown_vector_dimension = 1;
	/* BC */
	//   ls->bc_names_dirichlet = NULL;
	//   ls->bc_names_neumann = NULL;
	ls->unknown_vector_indeces = NULL;
	ls->unknown_node_numbers = NULL;
	ls->unknown_update_methods = NULL;
	// wW   ls->bc_name = NULL;
	//    ls->bc_name2 = NULL;
	//    ls->bc_name3 = NULL;
	//    ls->boundary_conditions_function = NULL;
	/* SS */
	//    ls->sousin_name1 = NULL;
	//    ls->sousin_name2 = NULL;
	//    ls->sousin_name3 = NULL;
	//    ls->source_multiplicator_function = NULL;
	//    ls->source_multiplicator_function1 = NULL;
	//    ls->source_multiplicator_function2 = NULL;

	ls->lsp_name = NULL;

	ls->r = NULL;
	ls->xx = NULL;
	ls->memory = NULL;

	ls->memory_number = 0;
	ls->new_memory = NULL;

	ls->assemble_function = NULL;

	ls->system_time_name = NULL;
	ls->system_time_assemble_function_name = NULL;

	/* Multiple unknowns (Dim) / multiple solvers (MultSolvers) */
	ls->name_group_ls = NULL;
	ls->name_ls = NULL;
	ls->number_ls = -1;
	ls->num_of_unknowns_ls = 0;
	return ls;
}

/*************************************************************************
   ROCKFLOW - Funktion: CreateLinearSolverDim

   Aufgabe:
   Konstruktor fuer LINEAR_SOLVER

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E char *name: Zeiger auf den Namen des LS's.

   Ergebnis:
   - Adresse des LS's im Erfolgsfall, ansonsten NULL-Zeiger -

   Programmaenderungen:
   02/1999     AH         Erste Version
   08/1999     AH         Erweiterung memory  (neq)
   10/1999     AH         Systemzeit
   10/1999     AH         Systemzeit
   05/2000     OK         Erweiterung fuer vektorielle Groessen
   08/2000     MK         CreateLinearSolver->CreateLinearSolverDim
   04/2003     MK         wird durch CreateLinearSolver wieder abgeloest
*************************************************************************/
LINEAR_SOLVER* CreateLinearSolverDim(long store, int unknown_vector_dimension, long dim)
{
	LINEAR_SOLVER* ls;

	ls = (LINEAR_SOLVER*)Malloc(sizeof(LINEAR_SOLVER));
	if (ls == NULL)
		return NULL;

	// OK    ls->name = NULL;
	ls->store = store;
	ls->dim = dim;
	ls->matrix = NULL;
	ls->b = NULL;
	ls->x = NULL;
	ls->lsp = NULL;
	ls->unknown_vector_dimension = unknown_vector_dimension;
	/* BC */
	/* //WW
	   ls->bc_names_dirichlet = NULL;
	   ls->bc_names_neumann = NULL;
	 */
	ls->unknown_vector_indeces = NULL;
	ls->unknown_node_numbers = NULL;
	ls->unknown_update_methods = NULL;
	/* //WW
	   ls->bc_name = NULL;
	   ls->bc_name2 = NULL;
	   ls->bc_name3 = NULL;
	   ls->boundary_conditions_function = NULL;
	   // SS
	   ls->sousin_name1 = NULL;
	   ls->sousin_name2 = NULL;
	   ls->sousin_name3 = NULL;
	   ls->source_multiplicator_function = NULL;
	   ls->source_multiplicator_function1 = NULL;
	   ls->source_multiplicator_function2 = NULL;
	 */

	ls->lsp_name = NULL;

	ls->r = NULL;
	ls->xx = NULL;
	ls->memory = NULL;

	ls->memory_number = 0;
	ls->new_memory = NULL;

	ls->assemble_function = NULL;

	ls->system_time_name = NULL;
	ls->system_time_assemble_function_name = NULL;

	/* Multiple unknowns (Dim) / multiple solvers (MultSolvers) */
	ls->name_group_ls = NULL;
	ls->name_ls = NULL;
	ls->number_ls = -1;
	ls->num_of_unknowns_ls = 0;
	return ls;
}

//////////////////////////////////////////////////////////////////////////
// NUM
//////////////////////////////////////////////////////////////////////////

double GetNumericalTimeCollocation(char* name)
{
	CRFProcess* m_pcs = NULL;
	m_pcs = PCSGet(name);
	if (m_pcs)
		return m_pcs->m_num->ls_theta;
	else
		cout << "Fatal error in GetNumericalTimeCollocation: No valid PCS"
		     << "\n";
	// OK return np_array[GetNumericsIndex(name)].time_collocation;
	return 1.0;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: CalcIterationError
 */
/* Aufgabe:
   Ermittelt den Fehler bei Iterationen

 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E double *new_iteration : Vektor des neuen Iterationsschritts
   E double *old_iteration : Vektor des alten Iterationsschritts
   E double *reference     : Vektor des alten Zeitschritts (als Referenz)
   E longlength            : Laenge der Vektoren
   E int method            : Methode der Fehlerermittlung
 */
/* Ergebnis:
   - void -
 */
/* Programmaenderungen:
   1/1999     C.Thorenz  Zweite Version                                                                          */
/**************************************************************************/
double NUMCalcIterationError(double* new_iteration, double* old_iteration, double* reference, long length, int method)
{
	static long i;
	static double error, change, max_c, min_c;

	error = 0.;
	change = 0.;

	max_c = 0.;
	min_c = 1.e99;

	switch (method)
	{
		default:
		case 0:
			return 0.;
		/* Max. Unterschied zwischen altem und neuem Iterationsschritt */
		case 1:
			for (i = 0l; i < length; i++)
				error = max(error, fabs(new_iteration[i] - old_iteration[i]));
			return error;

		/* Max. Unterschied zwischen altem und neuem Iterationsschritt,
		   jeweils normiert mit dem Mittelwert der Groesse des Wertes  */
		case 2:
			for (i = 0l; i < length; i++)
				error = max(error, 2. * fabs(new_iteration[i] - old_iteration[i])
				                       / (fabs(new_iteration[i]) + fabs(old_iteration[i]) + MKleinsteZahl));
			return error;

		/* Max. Unterschied zwischen altem und neuem Iterationsschritt,
		   normiert mit dem groessten Wert */
		case 3:
			for (i = 0l; i < length; i++)
			{
				error = max(error, fabs(new_iteration[i] - old_iteration[i]));
				max_c = max(max(max_c, fabs(new_iteration[i])), fabs(old_iteration[i]));
			}
			return error / (max_c + MKleinsteZahl);

		/* Max. Unterschied zwischen altem und neuem Iterationsschritt,
		   normiert mit der Spanne der Werte */
		case 4:
			for (i = 0l; i < length; i++)
			{
				error = max(error, fabs(new_iteration[i] - old_iteration[i]));
				min_c = min(min_c, fabs(new_iteration[i]));
				max_c = max(max_c, fabs(new_iteration[i]));
			}
			return error / (max_c - min_c + MKleinsteZahl);

		/* Max. Unterschied zwischen altem und neuem Iterationsschritt,
		   normiert mit dem Unterschied zum alten Zeitschritt. Die
		   genaueste Methode, da die Fehlerberechnung dann Zeitschritt-
		   unabhaengig wird! */
		case 5:
			for (i = 0l; i < length; i++)
				error = max(error,
				            fabs(new_iteration[i] - old_iteration[i])
				                / (fabs(new_iteration[i] - reference[i]) + MKleinsteZahl));
			return error;

		/* Max. Unterschied zwischen altem und neuem Iterationsschritt,
		   normiert mit dem maximalen Unterschied zum alten Zeitschritt */
		case 6:
			for (i = 0l; i < length; i++)
			{
				error = max(error, fabs(new_iteration[i] - old_iteration[i]));
				change = max(change, fabs(new_iteration[i] - reference[i]));
			}
			return error / (change + MKleinsteZahl);

		/* Der Vektorabstand */
		case 7:
			return MVekDist(old_iteration, new_iteration, length);

		/* Der Vektorabstand der Iteration, normiert mit dem Vektorabstand
		   zur alten Zeitebene */
		case 8:
			return MVekDist(old_iteration, new_iteration, length)
			       / (MVekDist(reference, new_iteration, length) + MKleinsteZahl);
	}
}
#endif // ifndef NEW_EQS //WW. 06.11.2008

int GetNumericsGaussPoints(int element_dimension)
{
	int m_gaussian_points = 3;
	int g_gaussian_points = 3;
	switch (element_dimension)
	{
		case 1:
			m_gaussian_points = 1;
			break;
		case 2:
			m_gaussian_points = g_gaussian_points;
			break;
		case 3:
			m_gaussian_points = g_gaussian_points;
			break;
		case 4:
			m_gaussian_points = g_gaussian_points;
			break;
	}
	return m_gaussian_points;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   05/2005 OK Implementation
   last modification:
**************************************************************************/
CNumerics* NUMGet(string num_name)
{
	CNumerics* m_num = NULL;
	for (int i = 0; i < (int)num_vector.size(); i++)
	{
		m_num = num_vector[i];
		if (m_num->pcs_type_name.compare(num_name) == 0)
			return m_num;
	}
	return NULL;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   01/2005 OK Implementation
   last modified:
**************************************************************************/
void NUMDelete()
{
	long i;
	int no_num = (int)num_vector.size();
	for (i = 0; i < no_num; i++)
		delete num_vector[i];
	num_vector.clear();
}
