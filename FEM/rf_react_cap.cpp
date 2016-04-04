/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/*

    rf_react_cap.cpp

Reaction package to go with CHEMAPP

Dedong Li and Sebastian Bauer

Kiel, 11/2008


*/

//#include "stdafx.h" /* MFC */
#include <signal.h>
#include "display.h"
#include "makros.h"
#include "memory.h"
#include "rf_pcs.h"
//#include "nodes.h"
#include "rf_pcs.h"
#include "mathlib.h"
#include "rfmat_cp.h"
#include "rf_react_cap.h"
#include "rf_react.h"
#include "rf_react_int.h" // new reaction interface

#include "rf_ic_new.h"
#include "stdio.h"
//#include "geo_strings.h"

#include "rf_pcs.h" //OK_MOD"
#include "files0.h"
//#include "files.h"
//#include "nodes.h"
//#include "elements.h"           /* fuer ElGetElementGroupNumber */
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include "rf_tim_new.h"
#include "rf_mmp_new.h"
#include "rf_kinreact.h"
// Elem object
#include "fem_ele_std.h"

#ifdef OGS_FEM_CAP // CAP_REACT
#include "cacint.h" //DL
#else
/* Length of a TQ String */
#define TQSTRLEN 25
#endif
#include "CAP_IO.h"
#include "VLE.h"
#include "HKF.h"
// CB2406 #endif

using namespace std;
/*
#define PHREEQC
int NR_INORG_COMP_EQU;
int NR_MIN_EQU;
extern char *file_name;
extern char *crdat; // MX
//REACTION_MODEL *rcml=NULL;
extern double gravity_constant;
*/

vector<REACT_CAP*> REACT_CAP_vec;

/**************************************************************************/
/* Constructor */
REACT_CAP::REACT_CAP(void)
{
	flag_cap = false; /* DL 28,10,08*/
	data_file = "false";
	data_format = "false";
	mass_type = "false";
	check_no_reaction_nodes = false; // ToDo
	nodenumber = 0;

	CAP_MODE = 0; // default case, when Chemapp IS available
}
/* Destructor */
REACT_CAP::~REACT_CAP(void)
{
}
/**************************************************************************/

/**************************************************************************
   Function: REACT_CAP::Read

   Task:
   Read ChemApp interface file *.cap

   Programming:
   11/2008     DL/SB         First Version
**************************************************************************/
ios::pos_type REACT_CAP::Read(ifstream* rfd_file, const GEOLIB::GEOObjects& geo_obj, const std::string& unique_name)
{
	bool new_keyword = false;
	std::string hash("#");
	std::string line_string;
	std::stringstream in;
	ios::pos_type position;

	int ipx, ik = 0, ikh = 0, ieh = 0, j;
	std::string ncomp_x, idphase;
	char* ncomp;
	long value;
	std::string geometry_type, geometry_name;
	std::string line_str1, line_str2;
	double diff_value;
	vector<std::string> pcs_name1;
	vector<double> pcs_stoi1;

	this->species_relative_activity.clear();
	this->species_relative_activity_name.clear();

	this->show_data_file = false;
	this->Node_All_SI.clear();

	// go through input file
	while (!new_keyword)
	{
		line_string = GetLineFromFile1(rfd_file);
		// if(line_string.size() < 1) break;

		if (line_string.find(hash) != string::npos)
		{
			new_keyword = true;
			break;
		}

		/* read keywords */
		//....................................................................
		if (line_string.find("$CAP_MODE") != string::npos)
		{ // subkeyword found
			in.str(GetLineFromFile1(rfd_file));
			in >> CAP_MODE; // sub_line
			in.clear();
		}
		//....................................................................
		if (line_string.find("$DATA_FILE") != string::npos)
		{ // subkeyword found
			in.str(GetLineFromFile1(rfd_file));
			in >> data_file; // sub_line
			in.clear();
		}
		//....................................................................
		if (line_string.find("$DATA_FORMAT") != string::npos)
		{ // subkeyword found
			in.str(GetLineFromFile1(rfd_file));
			in >> data_format; // sub_line
			in.clear();
		}
		if (line_string.find("$SHOW_DATA_FILE") != string::npos) // subkeyword found
			this->show_data_file = true;
		//....................................................................
		if (line_string.find("$MASS_NUM") != string::npos)
		{ // subkeyword found
			in.str(GetLineFromFile1(rfd_file));
			in >> mass_num; // sub_line
			in.clear();
		}
		//....................................................................
		if (line_string.find("$MASS_TYPE") != string::npos)
		{ // subkeyword found
			in.str(GetLineFromFile1(rfd_file));
			in >> mass_type; // sub_line
			in.clear();
		}
		//....................................................................
		if (line_string.find("$PHASE_LIST") != string::npos)
		{ // subkeyword found
			this->phase_list.clear();
			this->phase_post.clear();
			while (1)
			{
				in.str(GetLineFromFile1(rfd_file));
				in >> idphase;
				if (!idphase.find("END"))
				{
					in.clear();
					break;
				}
				this->phase_list.push_back(atoi(idphase.c_str()));
				ncomp_x = "";

				in >> ncomp_x;
				ncomp = new char[(int)ncomp_x.size() + 1];
				strcpy(ncomp, ncomp_x.c_str());
				this->phase_post.push_back(ncomp);
				in.clear();
			}
		}
		//....................................................................
		if (line_string.find("$SPECIES_DEFINE") != string::npos)
		{ // subkeyword found
			in.str(GetLineFromFile1(rfd_file));
			in >> species_define; // sub_line
			in.clear();
		}
		//....................................................................
		if (line_string.find("$SPECIES_LIST") != string::npos)
		{ // subkeyword found
			this->species_list_name.clear();
			this->species_list_phase.clear();
			while (1)
			{
				in.str(GetLineFromFile1(rfd_file));
				in >> ncomp_x;
				in.clear();
				if (!ncomp_x.find("END"))
					break;
				this->CompNamePhase(ncomp_x, ncomp, ipx);
				species_list_name.push_back(ncomp);
				species_list_phase.push_back(ipx);
			}
		}
		//....................................................................
		if (line_string.find("$KINETIC_SPECIES") != string::npos)
		{ // subkeyword found
			this->species_kin_name.clear();
			this->species_kin_phase.clear();
			while (1)
			{
				in.clear();
				in.str(GetLineFromFile1(rfd_file));
				in >> ncomp_x;
				in.clear();
				if (!ncomp_x.find("END"))
					break;
				this->CompNamePhase(ncomp_x, ncomp, ipx);
				species_kin_name.push_back(ncomp);
				species_kin_phase.push_back(ipx);
			}
		}

		//....................................................................
		if (line_string.find("$KINETIC_REACTIONS") != string::npos)
		{ // subkeyword found
			this->Kin_Reactions.clear();
			while (1)
			{
				in.str(GetLineFromFile1(rfd_file));
				in >> ncomp_x;
				if (!ncomp_x.find("END"))
				{
					if (ik > 0)
						this->Kin_Reactions.push_back(this->React);
					in.clear();
					break;
				}
				else if (!ncomp_x.find("REACTION"))
				{
					if (ik > 0)
						this->Kin_Reactions.push_back(this->React);
					in >> this->React.type;
					this->React.species_stoi.clear();
					this->React.species_idx.clear();
					this->React.species_name.clear();
					this->React.species_phase.clear();
					in.clear();
					ik++;
				}
				else
				{
					React.species_stoi.push_back((double)atof(ncomp_x.c_str()));
					in >> ncomp_x;
					in.clear();
					this->CompNamePhase(ncomp_x, ncomp, ipx);
					React.species_name.push_back(ncomp);
					React.species_phase.push_back(ipx);
					React.species_idx.push_back(-1);
				}
			}
		}
		//....................................................................
		if (line_string.find("$KINETIC_HKF_REACTIONS") != string::npos)
		{ // subkeyword found
			this->Kin_HKF_Reactions.clear();
			while (1)
			{
				in.str(GetLineFromFile1(rfd_file));
				in >> ncomp_x;
				if (!ncomp_x.find("END"))
				{
					if (ikh > 0)
						this->Kin_HKF_Reactions.push_back(this->React);
					in.clear();
					break;
				}
				else if (!ncomp_x.find("REACTION"))
				{
					if (ikh > 0)
						this->Kin_HKF_Reactions.push_back(this->React);
					in >> this->React.type;
					this->React.species_stoi.clear();
					this->React.species_idx.clear();
					this->React.species_name.clear();
					this->React.species_phase.clear();
					in.clear();
					ikh++;
				}
				else
				{
					React.species_stoi.push_back((double)atof(ncomp_x.c_str()));
					in >> ncomp_x;
					in.clear();
					this->CompNamePhase(ncomp_x, ncomp, ipx);
					React.species_name.push_back(ncomp);
					React.species_phase.push_back(ipx);
					React.species_idx.push_back(-1);
				}
			}
		}
		if (line_string.find("$NLOG_AC") != string::npos)
		{ // subkeyword found
			this->species_nlog_name.clear();
			this->species_nlog_phase.clear();
			while (1)
			{
				in.str(GetLineFromFile1(rfd_file));
				in >> ncomp_x;
				if (!ncomp_x.find("END"))
				{
					in.clear();
					break;
				}
				ncomp = new char[(int)ncomp_x.size()];
				strcpy(ncomp, ncomp_x.c_str());
				this->nlog_name.push_back(ncomp);
				in >> ncomp_x;
				in.clear();
				this->CompNamePhase(ncomp_x, ncomp, ipx);
				species_nlog_name.push_back(ncomp);
				species_nlog_phase.push_back(ipx);
			}
		}

		if (line_string.find("$PCS_RENAME") != string::npos)
		{ // subkeyword found
			this->pcs_rename0.clear();
			this->pcs_rename1.clear();
			while (1)
			{
				pcs_name1.clear();
				pcs_stoi1.clear();
				in.str(GetLineFromFile1(rfd_file));
				line_str2 = in.str();
				in >> ncomp_x;
				if (!ncomp_x.find("END"))
				{
					in.clear();
					break;
				}
				this->pcs_rename0.push_back(ncomp_x);
				// for(j=1;j<(int)REACT_PRQ::string2vector(line_str2).size();j+=2){
				for (j = 1; j < (int)REACTINT::string2vector(line_str2).size(); j += 2)
				{
					in >> ncomp_x;
					if (ncomp_x == "*")
						pcs_stoi1.push_back(999999);
					else if (ncomp_x == "/")
						pcs_stoi1.push_back(-999999);
					else
						pcs_stoi1.push_back((double)atof(ncomp_x.c_str()));
					in >> ncomp_x;
					pcs_name1.push_back(ncomp_x);
				}
				this->pcs_rename1.push_back(pcs_name1);
				this->pcs_rename_stoi.push_back(pcs_stoi1);
				in.clear();
			}
		}

		if (line_string.find("$PCS_PRE_RENAME") != string::npos)
		{ // subkeyword found
			this->pcs_rename0_pre.clear();
			this->pcs_rename1_pre.clear();
			while (1)
			{
				pcs_name1.clear();
				pcs_stoi1.clear();
				in.str(GetLineFromFile1(rfd_file));
				line_str2 = in.str();
				in >> ncomp_x;
				if (!ncomp_x.find("END"))
				{
					in.clear();
					break;
				}
				this->pcs_rename0_pre.push_back(ncomp_x);
				// for(j=1;j<(int)REACT_PRQ::string2vector(line_str2).size();j+=2){
				for (j = 1; j < (int)REACTINT::string2vector(line_str2).size(); j += 2)
				{
					in >> ncomp_x;
					if (ncomp_x == "*")
						pcs_stoi1.push_back(999999);
					else if (ncomp_x == "/")
						pcs_stoi1.push_back(-999999);
					else
						pcs_stoi1.push_back((double)atof(ncomp_x.c_str()));
					in >> ncomp_x;
					pcs_name1.push_back(ncomp_x);
				}
				this->pcs_rename1_pre.push_back(pcs_name1);
				this->pcs_rename_stoi_pre.push_back(pcs_stoi1);
				in.clear();
			}
		}

		if (line_string.find("$RELATIVE_ACTIVITY") != string::npos)
		{ // subkeyword found
			this->species_relative_activity_name.clear();
			this->species_relative_activity.clear();
			while (1)
			{
				in.str(GetLineFromFile1(rfd_file));
				line_str2 = in.str();
				in >> ncomp_x;
				if (!ncomp_x.find("END"))
				{
					in.clear();
					break;
				}

				this->CompNamePhase(ncomp_x, ncomp, ipx);
				this->species_relative_activity_name.push_back(ncomp);
				in >> ncomp_x;
				this->species_relative_activity.push_back((double)atof(ncomp_x.c_str()));
				in.clear();
			}
		}

		if (line_string.find("$NODE_ALL_SI") != string::npos)
		{ // subkeyword found
			this->Node_All_SI.clear();
			while (1)
			{
				in.str(GetLineFromFile1(rfd_file));
				line_str2 = in.str();
				in >> ncomp_x;
				if (!ncomp_x.find("END"))
				{
					in.clear();
					break;
				}
				this->Node_All_SI.push_back((int)atoi(ncomp_x.c_str()));
				in.clear();
			}
		}
		//....................................................................
		if (line_string.find("$REDOX_EH") != string::npos)
		{ // subkeyword found
			this->species_redox_name.clear();
			this->species_redox_phase.clear();
			this->redox_stoi.clear();
			while (1)
			{
				in.str(GetLineFromFile1(rfd_file));
				in >> ncomp_x;
				if (!ncomp_x.find("END"))
				{
					in.clear();
					break;
				}
				if (ieh < 2)
				{
					ncomp = new char[(int)ncomp_x.size()];
					strcpy(ncomp, ncomp_x.c_str());
					ipx = 0;
				}
				else
					this->CompNamePhase(ncomp_x, ncomp, ipx);
				species_redox_name.push_back(ncomp);
				species_redox_phase.push_back(ipx);
				in >> value;
				in.clear();
				this->redox_stoi.push_back(value);
				ieh++;
			}
		}
		//....................................................................
		if (line_string.find("$SET_BC_EQ_GEO") != string::npos)
		{ // subkeyword found
			this->ic_2_bc_geometry_type.clear();
			this->ic_2_bc_geometry_name.clear();
			this->ic_2_bc_GeoID.clear();
			while (1)
			{
				line_str1 = GetLineFromFile1(rfd_file);
				if ((line_str1.find("$") != string::npos) || (line_str1.find("#STOP") != string::npos)
				    || (line_str1.find("END") != string::npos))
				{
					break;
				}
				in.str(line_str1);
				in >> geometry_type >> geometry_name;
				ic_2_bc_geometry_type.push_back(geometry_type);
				ic_2_bc_geometry_name.push_back(geometry_name);

				// TF 06/2010 - for the change string-ID to size_t-ID
				size_t geo_obj_idx(std::numeric_limits<size_t>::max());

				if (geometry_type.find("POINT") != std::string::npos)
				{
					// get the point vector and set the geo_obj_idx
					if (!((geo_obj.getPointVecObj(unique_name))->getElementIDByName(geometry_name, geo_obj_idx)))
					{
						std::cerr << "error in CKinReactData::Read: (type=" << geometry_type << "): " << geometry_name
						          << " point name not found!"
						          << "\n";
						exit(1);
					}
				}
				if (geometry_type.find("POLYLINE") != std::string::npos)
				{
					// get the point vector and set the geo_obj_idx
					if (!((geo_obj.getPolylineVecObj(unique_name))->getElementIDByName(geometry_name, geo_obj_idx)))
					{
						std::cerr << "error in CKinReactData::Read: polyline name " << geometry_name << " not found!"
						          << "\n";
						exit(1);
					}
				}

				ic_2_bc_GeoID.push_back(geo_obj_idx);
				in.clear();
				// break;
			}
		}
		//....................................................................
		if (line_string.find("$SET_BC_EQ_SPECIES") != std::string::npos)
		{ // subkeyword found
			this->ic_2_bc_species.clear();
			this->ic_2_bc_species_value.clear();
			while (1)
			{
				line_str1 = GetLineFromFile1(rfd_file);
				if ((line_str1.find("$") != string::npos) || (line_str1.find("#STOP") != string::npos)
				    || (line_str1.find("END") != string::npos))
				{
					break;
				}
				in.str(line_str1);
				in >> geometry_name; // actually, this reads in a species name ...
				ic_2_bc_species.push_back(geometry_name);
				diff_value = 0.0;
				in >> diff_value;
				this->ic_2_bc_species_value.push_back(diff_value);
				in.clear();
				// break;
			}
		}
	} // end while
	return position;
}

void REACT_CAP::RecoverChemSystem(void)
{
	// int i,ii;
	// LI nphase, npcon, noerr;
	if (!this->species_define.compare("TRANS"))
		this->DefSpeciesInChemApp();
	else if (!this->species_define.compare("LIST"))
		this->DefSpeciesInChemAppList();
	else
		this->DefSpeciesInChemAppAll();
}

/* ChemApp subroutine for chemical reaction calculation   DL 28,10,08 */
void REACT_CAP::ExecuteReactionsChemApp(int f, int nodeflag)
{
	// CAP_MODE=2; // now from input file
	CAP_icount = 1;
	CAP_Time = 0;
	CAP_Node = 0;

	int ii, ok = 0, ik = 0;

	if (f == 0)
	{
		this->CreateREACT(); // set nodenumber and rateflag
		this->InitChemApp(); // read DATAFILE

		if (this->show_data_file)
		{
			cout << " ==> Show Data File"
			     << "\n";
			this->ShowDataFile();
		}

		this->LoadMassPCS(); // read process and get MASS TRANSPORT list
		this->CmpSpeciesListPCS(); // consistent checking between TRANS and DATAFILE
		this->CmpSpeciesListKIN(); // DL 03.2010

		if (!this->species_define.compare("TRANS"))
		{
			cout << " Define TRANS as Chemical Reaction Species "
			     << "\n";
			this->DefSpeciesInChemApp();
		}
		else if (!this->species_define.compare("LIST"))
		{
			cout << " Define LIST as Chemical Reaction Species "
			     << "\n";
			this->CmpSpeciesListDAT(); // consistent checking between .cap species_list and DATAFILE
			if (this->CmpSpeciesList()) // consistent checking between .cap species_list and TRANS
				cout << " The two species lists are compatible. OK..."
				     << "\n";
			else
				cout << " Incompatibility checking ERROR..."
				     << "\n";
			this->DefSpeciesInChemAppList();
		}
		else
		{
			// this->UndefOutSpecies();
			cout << " UNDEF SPECIES, USING THE TOTAL SPECIES LIST "
			     << "\n";
			this->DefSpeciesInChemAppAll();
		}

		if (!this->mass_type.compare("ELEMENT"))
		{
			cout << " SET MASS TRANSPORT AS ELEMENT "
			     << "\n";
			if (CAP_MODE == 2)
				cout << " File Mode "
				     << "\n";
			else
				this->ResetElement();
		}
		else
			cout << " SET MASS TRANSPORT AS SPECIES (defualt) "
			     << "\n";

		this->CmpSpeciesListKinReact(); // DL 03.2010

		this->SetKinHKFspecies(); // DL 11.2010
		if (this->Kin_HKF_species.size() > 0)
			this->LoadHKFparam();

		this->CmpSpeciesListNLOG();
		this->CmpSpeciesListREDOX();

		this->CmpSpeciesListRelativeActivity();

		this->warning_out.clear();
		this->warning_out.open("ChemApp_Warning.txt", ios::out);

		cout << "-------------------------------"
		     << "\n";
		cout << "      END CHEMAPP INITIAL      "
		     << "\n";
		cout << "-------------------------------"
		     << "\n";

		// this->KinInit();

		cout.flush() << " Calculating geochemical equilibrium at all nodes using initial values."
		             << "\n";
		this->LoopNodeReact(0, nodeflag);

	} // f == 0

	if (f == 1 || f == -1)
	{
		// Check for nodes without reactions
		if ((int)this->check_no_reaction_nodes == false)
		{
			ok = this->CheckNoReactionNodes(); // default rateflag=1 for reaction calc
			if (!ok)
				cout << "Error when checking for nodes without reactions"
				     << "\n"
				     << "\n";
		}

		for (ii = 0; ii < this->nodenumber; ii++)
			if (this->rateflag[ii] != 0)
				ik++;
		cout << " Calculating geochemical equilibrium at " << ik << " nodes."
		     << "\n";

		if (f == 1)
			this->LoopNodeReact(1, nodeflag);
		else if (f == -1)
			this->LoopNodeReact_Liquid_Vapor(1, nodeflag);
	}

	// if(f==2){
	//	for(ii=0;ii<this->nodenumber;ii++) if(this->rateflag[ii] != 0) ik++;
	//	cout << " Calculating geochemical kinetic at " << ik << " nodes." << "\n";
	//}
}

/***********************************************************************************************

************************************************************************************************/
int REACT_CAP::CheckNoReactionNodes(void)
{
	int ok = 0;
	long l;
	cout << " CheckNoReactionNodes "
	     << "\n";

	CFEMesh* m_msh = fem_msh_vector[0]; // SB: ToDo hart gesetzt
	if (m_msh == NULL)
	{
		cout << "No mesh in CheckNoReactionNodes"
		     << "\n";
		exit(0);
	}

	// CB 19.1.2011
	// Get the reaction interface data for checking for dried out nodes from eclipse coupling --> rateflag = 0
	REACTINT* m_rei = NULL;
	if (REACTINT_vec.size() > 0)
		m_rei = REACTINT_vec[0];

	if (aktueller_zeitschritt < 2)
	{ // do not in the very first calculation before first time step and in the first time step
		this->check_no_reaction_nodes = false;
	}
	else
	{
		CKinReactData* m_krd = NULL;
		if (KinReactData_vector.size() > 0)
			m_krd = KinReactData_vector[0];
		if (m_krd == NULL)
		{
			// no KinReactData specified in *.krc file
			cout << "No CKinReactData available in CheckNoReactionNodes"
			     << "\n";
		}
		else
		{
			if (m_krd->is_a_CCBC.size() > 0)
			{ // reaction nodes specified in krc input file
				// Initialize vector is_a_CCBC
				for (l = 0; l < (long)m_msh->nod_vector.size(); l++)
				{ // node 1 needed for phreeqc-input
					this->rateflag[l] = 1;
					if (m_krd->is_a_CCBC[l] == true)
						this->rateflag[l] = 0; // rateflag == 0 means no reactions calculated
					// cout << " Node " << l << " is " << this->rateflag[l] << "\n";
				}
			}
		}

		// also, switch off reactions for nodes with Sat Water < WaterSatLimit to avoid dryout problem in Eclipse
		// coupling
		if (m_rei)
		{
			if (m_rei->s_water_limit)
				for (l = 0; l < (long)m_msh->nod_vector.size(); l++)
					if (m_rei->dried_out_nodes[l])
						rateflag[l] = 0;
		}

		this->check_no_reaction_nodes = true;
	}

	ok = 1;
	return ok;
}

void REACT_CAP::CreateREACT(void)
{
	int i, vector_size;
	long l;
	CRFProcess* m_pcs = NULL;

	vector_size = (int)pcs_vector.size();
	for (i = 0; i < vector_size; i++)
	{
		m_pcs = pcs_vector[i];
		// if(m_pcs->pcs_type_name.compare("MASS_TRANSPORT")==0){
		if (m_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT)
		{
			nodenumber = (long)m_pcs->m_msh->GetNodesNumber(false);
			// elenumber = (long) m_pcs->m_msh->ele_vector.size();
		}
	}
	rateflag = (int*)Malloc(sizeof(int) * nodenumber);
	for (l = 0; l < nodenumber; l++)
		rateflag[l] = 1;
	CAP_count = 0;
}

/**************************************************************
   Function: REACT_CAP::InitChemApp

   Task:
   Loading DATABASE file.
                                                    DL 12.2008
**************************************************************/
void REACT_CAP::InitChemApp(void)
{
	LI noerr, version;
	// CRFProcess* m_pcs = NULL;
	char* cap_datafile;

#ifndef OGS_FEM_CAP // CAP_REACT
	cout << " Warning in REACT_CAP::InitChemApp():"
	     << "\n";
	cout << "    ChemApp source files are not included in this configuration."
	     << "\n";
	cout << "    Installation of ChemApp and usage of CMake option OGS_FEM_CAP"
	     << "\n";
	cout << "    is required for coupled OGS_CAP simulation."
	     << "\n";

	if (CAP_MODE != 2)
	{
		cout << "    CAP_MODE = " << CAP_MODE << " was defined in *.cap input file:"
		     << "\n";
		if (CAP_MODE == 0)
		{
			cout << "      (normal ChemApp calculation, this is the default case, "
			     << "\n";
			cout << "      also when $CAP_MODE Keyword is not used in *.cap input file."
			     << "\n";
		}
		else if (CAP_MODE == 1)
		{
			cout << "      (normal ChemApp calculation with storage of results in text files."
			     << "\n";
		}
		cout << "      Options for $CAP_MODE Keyword: "
		     << "\n";
		cout << "      0: default, 1: write result files, 2: read result files)"
		     << "\n";
		CAP_MODE = 2; // default case, when Chemapp IS NOT available
		cout << " Setting CAP_MODE = " << CAP_MODE << ": Read ChemApp results from precalculated files."
		     << "\n";
	}
	else
		cout << " CAP_MODE = " << CAP_MODE << ": Read ChemApp results from precalculated files."
		     << "\n";
#endif

	if (CAP_MODE == 2)
	{
		// check for results file folder
		cout << "\n"
		     << "  CAP_MODE = " << CAP_MODE << "\n";
		cout << "  Checking OGS input files..."
		     << "\n";
		bool files = CAP_check_file();

		if (files == false)
		{
			cout << "\n"
			     << "  Warning: No ChemApp results files available. Exiting now...";
			exit(0);
		}
		cout << "  Read ChemApp results from precalculated files"
		     << "\n";
	}

	CAP_tqini(&noerr);
	CAP_tqvers(&version, &noerr);
	cout << "\n"
	     << "---------------------------------------------------------------"
	     << "\n";
	cout << "   This program has been compiled with ChemApp version " << version << "\n";
	cout << "   ChemApp Module for Chemical Reaction Calculation...Start... "
	     << "\n";
	cout << "---------------------------------------------------------------"
	     << "\n";
	// CGSProject* m_gsp = NULL;
	// m_gsp = GSPGetMember("mmp");
	// if(m_gsp)
	//  cap_datafile = new char [m_gsp->path.length() + this->data_file.length()];
	// strcpy(cap_datafile, m_gsp->path.c_str());
	// strcat(cap_datafile, this->data_file.c_str());
	cap_datafile = new char[this->data_file.length()];
	strcpy(cap_datafile, this->data_file.c_str());
	cout << "THERMOCHEM DATA FILE --> " << cap_datafile << "\n";
	noerr = 1;
	if (data_format == "DAT")
	{
		CAP_tqopna(cap_datafile, 10, &noerr);
		CAP_tqrfil(&noerr);
	}
	if (data_format == "BIN")
	{
		CAP_tqopnb(cap_datafile, 10, &noerr);
		CAP_tqrbin(&noerr);
	}
	if (data_format == "CST")
	{
		CAP_tqopnt(cap_datafile, 10, &noerr);
		CAP_tqrcst(&noerr);
	}
	if (!noerr)
		cout << " READ DATA FILE OK ..."
		     << "\n"
		     << "\n";
	else
	{
		cout << " READ DATA FILE ERROR !!! STOP."
		     << "\n"
		     << "\n";
		exit(1);
	} // STOP and EXIT
	CAP_tqclos(10, &noerr);

	// pH Eh save as species output
	// for(ii=this->mass_num; ii<(int)pcs_mass_idx.size(); ii++){
	//	m_pcs = pcs_vector[pcs_mass_idx[ii]];
	//	if(strcmp(m_pcs->pcs_primary_function_name[0], "pH")==0 || strcmp(m_pcs->pcs_primary_function_name[0],
	//"Eh")==0){
	//		cout << "     --> " << m_pcs->pcs_primary_function_name[0] << "\n";
	//		pcs_ospecies_idx.push_back(pcs_mass_idx[ii]);
	//	}
	//}
}

void REACT_CAP::ShowDataFile(void)
{
	LI noerr, nscom, npcon, i;
	char scname[TQSTRLEN];

	CAP_tqnosc(&nscom, &noerr);
	printf("--> %li elements\n", nscom);
	for (i = 1; i <= nscom; i++)
	{ // elements list
		CAP_tqgnsc(i, scname, &noerr);
		printf("%li: %s\n", i, scname);
	}

	CAP_tqnop(&nscom, &noerr);
	printf("--> %li phases\n", nscom);
	for (i = 1; i <= nscom; i++)
	{ // phases list
		CAP_tqgnp(i, scname, &noerr);
		printf("%li: %s\n", i, scname);
	}

	CAP_tqnopc(1, &npcon, &noerr);
	printf("--> %li constituents in gas phase\n", npcon);
	for (i = 1; i <= npcon; i++)
	{ // gas species
		CAP_tqgnpc(1, i, scname, &noerr);
		printf("%li: %s\n", i, scname);
	}

	CAP_tqnopc(2, &npcon, &noerr);
	printf("--> %li constituents in aqueous solution\n", npcon);
	for (i = 1; i <= npcon; i++)
	{ // aqueous solution species
		CAP_tqgnpc(2, i, scname, &noerr);
		printf("%li: %s\n", i, scname);
	}
}

void REACT_CAP::CompNamePhase(std::string ncomp_x, char*& ncomp, int& ipx)
{
	size_t i, np;
	int flag = 0;
	std::string postfix;
	for (i = 0; i < this->phase_list.size(); i++)
	{
		postfix = phase_post[i];
		np = ncomp_x.size() - postfix.size();
		if (ncomp_x.rfind(postfix) == np && postfix != "")
		{
			flag = 1;
			ipx = phase_list[i];
			ncomp = (char*)malloc(25); // new char [np];  //CB-DL??
			// ncomp = new char [np]; // (char *)malloc(25);  // CB200412: I'm not sure, mallloc was in my version
			strcpy(ncomp, ncomp_x.substr(0, np).c_str());
		}
		if (flag == 0 && postfix == "")
		{
			flag = 1;
			ipx = phase_list[i];
			ncomp = new char[np];
			strcpy(ncomp, ncomp_x.substr(0, np).c_str());
		}
	}
	if (flag == 0)
	{
		cout << " ERROR!...  " << ncomp_x << "  out of PHASE LIST. STOP!"
		     << "\n";
		exit(0);
	}
}

/**************************************************************
    Function: REACT_CAP::LoadMassPCS

    Task:
    Loading Mass_Transport from Processes list.
    Get index of phases by postfix of name
    check name of species validity
                                                    DL 01.2009
**************************************************************/
void REACT_CAP::LoadMassPCS(void)
{
	int i, ims, no_pcs, ipx; // number of process and transport species  ipx 0-solid 1-gas 2-aqueous
	CRFProcess* m_pcs = NULL;
	char* ncomp = new char[TQSTRLEN];
	std::string ncomp_x;

	int j, jx, nrs, is_idx, nrs_pre, is_idx_pre;
	vector<int> idx0, idx1, idx0_pre, idx1_pre;
	vector<vector<int> > pcs_idx1, pcs_idx1_pre;

	cout << "\n"
	     << " => LoadMassPCS() "
	     << "\n";
	cout << " REACTIVE_MASS_NUM " << this->mass_num << "\n";
	pcs_mass_idx.clear();
	ims = 0;
	// loading in mass_transport from processes list
	species_idx.clear();
	species_phase.clear();
	species_name.clear();
	species_mobil.clear();
	no_pcs = (int)pcs_vector.size();

	idx0.clear();
	pcs_idx1.clear();
	nrs = (int)pcs_rename0.size();
	for (j = 0; j < nrs; j++)
	{
		idx0.push_back(-1);
		idx1.clear();
		for (jx = 0; jx < (int)pcs_rename1[j].size(); jx++)
			idx1.push_back(-1);
		pcs_idx1.push_back(idx1);
	}

	idx0_pre.clear();
	pcs_idx1_pre.clear();
	nrs_pre = (int)pcs_rename0_pre.size();
	for (j = 0; j < nrs_pre; j++)
	{
		idx0_pre.push_back(-1);
		idx1_pre.clear();
		for (jx = 0; jx < (int)pcs_rename1_pre[j].size(); jx++)
			idx1_pre.push_back(-1);
		pcs_idx1_pre.push_back(idx1_pre);
	}

	for (i = 0; i < no_pcs; i++)
	{
		m_pcs = pcs_vector[i];
		cout << left << " Process " << setw(6) << i << setw(18) << m_pcs->getProcessType() << " "
		     << m_pcs->pcs_primary_function_name[0] << "\n";
		// if(!m_pcs->pcs_type_name.compare("MASS_TRANSPORT")){
		if (m_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT)
		{
			pcs_mass_idx.push_back(i); // mass index storage
			species_mobil.push_back(cp_vec[m_pcs->pcs_component_number]->mobil);
			if (ims < mass_num)
			{
				ims++;
				ncomp_x = m_pcs->pcs_primary_function_name[0];
				this->CompNamePhase(ncomp_x, ncomp, ipx);
				species_name.push_back(ncomp);
				species_phase.push_back(ipx);
			}
		}

		for (j = 0; j < nrs; j++)
		{
			if (m_pcs->pcs_primary_function_name[0] == pcs_rename0[j])
				idx0[j] = i;
			for (jx = 0; jx < (int)pcs_rename1[j].size(); jx++)
				if (m_pcs->pcs_primary_function_name[0] == pcs_rename1[j][jx])
					pcs_idx1[j][jx] = i;
		}

		for (j = 0; j < nrs_pre; j++)
		{
			if (m_pcs->pcs_primary_function_name[0] == pcs_rename0_pre[j])
				idx0_pre[j] = i;
			for (jx = 0; jx < (int)pcs_rename1_pre[j].size(); jx++)
				if (m_pcs->pcs_primary_function_name[0] == pcs_rename1_pre[j][jx])
					pcs_idx1_pre[j][jx] = i;
		}
	}

	// set pcs rename
	cout << "\n"
	     << " PCS_RENAME "
	     << "\n";
	this->pcs_rename_idx0.clear();
	this->pcs_rename_idx1.clear();
	for (i = 0; i < nrs; i++)
	{
		is_idx = 1;
		if (idx0[i] == -1)
			is_idx = 0;
		for (jx = 0; jx < (int)pcs_rename1[i].size(); jx++)
			if (pcs_idx1[i][jx] == -1)
				is_idx = 0;
		if (is_idx == 0)
		{
			cout << " Warning!!!, line " << i + 1 << " can not be found in PCS name list ! "
			     << "\n";
			exit(0);
		}
		pcs_rename_idx0.push_back(idx0[i]);
		pcs_rename_idx1.push_back(pcs_idx1[i]);
	}
	for (i = 0; i < int(pcs_rename_idx0.size()); i++)
	{
		cout << " " /*<< setw(4) << pcs_rename_idx0[i]*/ << setw(12) << pcs_rename0[i] << " <--- ";
		for (jx = 0; jx < (int)pcs_rename1[i].size(); jx++)
			cout << "    " << pcs_rename_stoi[i][jx] << "  " << pcs_rename1[i][jx];
		cout << "\n";
	}

	// set pcs rename pre
	cout << "\n"
	     << " PCS_PRE_RENAME "
	     << "\n";
	this->pcs_rename_idx0_pre.clear();
	this->pcs_rename_idx1_pre.clear();
	for (i = 0; i < nrs_pre; i++)
	{
		is_idx_pre = 1;
		if (idx0_pre[i] == -1)
			is_idx_pre = 0;
		for (jx = 0; jx < (int)pcs_rename1_pre[i].size(); jx++)
			if (pcs_idx1_pre[i][jx] == -1)
				is_idx_pre = 0;
		if (is_idx_pre == 0)
		{
			cout << " Warning!!!, line " << i + 1 << " can not be found in PCS name list ! "
			     << "\n";
			exit(0);
		}
		pcs_rename_idx0_pre.push_back(idx0_pre[i]);
		pcs_rename_idx1_pre.push_back(pcs_idx1_pre[i]);
	}
	for (i = 0; i < int(pcs_rename_idx0_pre.size()); i++)
	{
		cout << " " /*<< setw(4) << pcs_rename_idx0_pre[i]*/ << setw(12) << pcs_rename0_pre[i] << " <--- ";
		for (jx = 0; jx < (int)pcs_rename1_pre[i].size(); jx++)
			cout << "    " << pcs_rename_stoi_pre[i][jx] << "  " << pcs_rename1_pre[i][jx];
		cout << "\n";
	}
}

/**************************************************************
    Function: REACT_CAP::CmpSpeciesListREDOX

    Task:
    Compare the species REDOX_EH list with DATAFILE.
    Get index of phases and species.
    Get pcs_No. for Eh output species.
                                                    DL 01.2009
**************************************************************/
bool REACT_CAP::CmpSpeciesListREDOX(void)
{
	bool sr = true;
	LI ipc, noerr = 0;
	int i, ii, ns; // number of process and transport species
	CRFProcess* m_pcs = NULL;
	std::string ncomp_x;

	cout << "\n"
	     << "\n"
	     << " => CmpSpeciesListREDOX() "
	     << "\n";
	// compare REDOX_EH LIST with DATAFILE and get redox_idx
	this->pcs_redox = -1;
	ns = (int)species_redox_name.size();
	for (i = 0; i < ns; i++)
	{
		if (i == 0)
		{
			for (ii = this->mass_num; ii < (int)pcs_mass_idx.size(); ii++)
			{
				m_pcs = pcs_vector[pcs_mass_idx[ii]];
				if (this->species_mobil[ii] == 0
				    && strcmp(m_pcs->pcs_primary_function_name[0], species_redox_name[i]) == 0)
				{
					cout << " --> " << species_redox_name[i] << "\n";
					pcs_redox = pcs_mass_idx[ii];
				}
			}
			species_redox_idx.push_back(0);
		}
		else if (i == 1) // electron element
			species_redox_idx.push_back(0);
		else
		{
			if (species_redox_phase[i] == 0)
				CAP_tqinp(species_redox_name[i], &ipc, &noerr);
			else
				CAP_tqinpc(species_redox_name[i], species_redox_phase[i], &ipc, &noerr);
			if (noerr)
			{
				cout << " ERROR !!! SPECIES REDOX INPUT out of DATAFILE. STOP."
				     << "\n";
				exit(1);
			} // STOP and EXIT
			species_redox_idx.push_back(ipc);
			cout << " " << setw(4) << redox_stoi[i] << setw(4) << species_redox_phase[i] << setw(4)
			     << species_redox_idx[i] << species_redox_name[i] << "\n";
		}
	}
	cout << " REDOX_Eh OF SPECIES DEFINE OK ..."
	     << "\n";
	return sr;
}

/**************************************************************
    Function: REACT_CAP::CmpSpeciesListNLOG

    Task:
    Compare the species NLOG_AC list with DATAFILE.
    Get index of phases and species.
    Get pcs_No. for nlog output species.
                                                    DL 01.2009
**************************************************************/
bool REACT_CAP::CmpSpeciesListNLOG(void)
{
	bool sr = true;
	LI ipc, noerr = 0;
	int i, ii, ns; // number of process and transport species
	CRFProcess* m_pcs = NULL;
	std::string ncomp_x;

	cout << "\n"
	     << " => CmpSpeciesListNLOG() "
	     << "\n";
	// compare NLOG_AC LIST with DATAFILE and get list_idx
	pcs_nspecies_idx.clear(); // clear up nlog-pcs-idx vector
	nspecies_idx.clear(); // clear up nlog-idx vector
	ns = (int)species_nlog_name.size();
	for (i = 0; i < ns; i++)
	{
		noerr = 0;
		if (species_nlog_phase[i] == 0)
			CAP_tqinp(species_nlog_name[i], &ipc, &noerr);
		else
			CAP_tqinpc(species_nlog_name[i], species_nlog_phase[i], &ipc, &noerr);
		if (noerr)
		{
			cout << " ERROR !!! SPECIES NLOG INPUT out of DATAFILE. STOP."
			     << "\n";
			exit(1);
		} // STOP and EXIT
		species_nlog_idx.push_back(ipc);
		cout << " " << setw(8) << nlog_name[i] << setw(4) << species_nlog_phase[i] << setw(6) << species_nlog_idx[i]
		     << species_nlog_name[i] << "\n";

		// DL 02.02.2009 nlog value output
		for (ii = this->mass_num; ii < (int)pcs_mass_idx.size(); ii++)
		{
			m_pcs = pcs_vector[pcs_mass_idx[ii]];
			if (this->species_mobil[ii] == 0 && strcmp(m_pcs->pcs_primary_function_name[0], nlog_name[i]) == 0)
			{
				cout << " --> " << species_nlog_name[i] << "\n";
				pcs_nspecies_idx.push_back(pcs_mass_idx[ii]);
				nspecies_idx.push_back(i);
			}
		}
	}
	cout << " NLOG(ACTIVITY) OF SPECIES DEFINE OK ..."
	     << "\n";
	return sr;
}

/**************************************************************
    Function: REACT_CAP::CmpSpeciesListPCS

    Task:
    Compare the species MASS_TRANSPORT list with DATAFILE.
    Get index of phases and species.
                                                    DL 01.2009
**************************************************************/
bool REACT_CAP::CmpSpeciesListPCS(void)
{
	bool sr = true;
	LI ipc, noerr = 0;
	int i, ns; // number of process and transport species
	// CRFProcess* m_pcs = NULL;
	std::string ncomp_x;

	cout << "\n"
	     << "\n"
	     << " => CmpSpeciesListPCS() "
	     << "\n";
	// check up consistency between pcs_mass_transport and DATAFILE
	cout << " Pcs Phase Index Species List for Chemical Reactive Mass"
	     << "\n";
	ns = (int)species_name.size();
	for (i = 0; i < ns; i++)
	{
		if (strcmp(species_name[i], "EA") == 0)
			ipc = 0;
		else
		{
			noerr = 0;
			// std::string name_st(species_name[i]);
			// if(name_st.size()>24){
			//	species_name[i] = new char [24];
			//	strcpy(species_name[i], name_st.substr(0,24).c_str());
			//}
			// cout << species_name[i] << " " << name_st << " " << name_st.size() << "\n";
			if (species_phase[i] == 0)
				CAP_tqinp(species_name[i], &ipc, &noerr);
			else
				CAP_tqinpc(species_name[i], species_phase[i], &ipc, &noerr);
			if (noerr)
			{
				cout << " SPECIES INPUT ERROR !!! STOP."
				     << "\n";
				exit(1);
			} // STOP and EXIT
		}
		species_idx.push_back(ipc);
		species_dormant.push_back(0);
		cout << " " << setw(4) << pcs_mass_idx[i] << setw(6) << species_phase[i] << setw(6) << species_idx[i]
		     << species_name[i] << "\n";
	}
	cout << " TRANSPORT SPECIES INPUT OK ..."
	     << "\n";
	return sr;
}

/**************************************************************
    Function: REACT_CAP::CmpSpeciesListDAT

    Task:
    Compare the two species .cap list with DATAFILE.
    Get index of phases and species.
    Get pcs_No. for output species.
    respectively, find out incompatible species.
                                                    DL 01.2009
**************************************************************/
bool REACT_CAP::CmpSpeciesListDAT(void)
{
	bool sr = true;
	LI ipc, noerr = 0;
	int i, ii, ns; // number of process and transport species
	CRFProcess* m_pcs = NULL;
	std::string ncomp_x;

	cout << "\n"
	     << "\n"
	     << " => CmpSpeciesListDAT() "
	     << "\n";
	// compare LIST with DATAFILE and get list_idx
	pcs_ospecies_idx.clear(); // clear up species output vector
	cout << " SPECIES LIST FOR CHEMICAL REACTION CALCULATION "
	     << "\n";
	cout << " Phase Index Species List for Chemical Reaction"
	     << "\n";
	ns = (int)species_list_name.size();
	for (i = 0; i < ns; i++)
	{
		if (species_list_phase[i] == 0)
			CAP_tqinp(species_list_name[i], &ipc, &noerr);
		else
			CAP_tqinpc(species_list_name[i], species_list_phase[i], &ipc, &noerr);
		if (noerr)
		{
			cout << " ERROR !!! SPECIES LIST INPUT out of DATAFILE. STOP."
			     << "\n";
			exit(1);
		} // STOP and EXIT
		species_list_idx.push_back(ipc);
		cout << " " << setw(6) << species_list_phase[i] << setw(6) << species_list_idx[i] << species_list_name[i]
		     << "\n";

		// DL 21.01.2009 species value output
		if (species_list_phase[i] == 2)
			for (ii = this->mass_num; ii < (int)pcs_mass_idx.size(); ii++)
			{
				m_pcs = pcs_vector[pcs_mass_idx[ii]];
				if (/*this->species_mobil[ii]==0 &&*/ strcmp(m_pcs->pcs_primary_function_name[0], species_list_name[i])
				    == 0)
				{
					cout << "     --> " << species_list_name[i] << "\n";
					pcs_ospecies_idx.push_back(pcs_mass_idx[ii]);
				}
			}
	}

	cout << " REACTION SPECIES LIST INPUT OK ..."
	     << "\n";
	return sr;
}

/**************************************************************
    Function: REACT_CAP::CmpSpeciesListKIN

    Task:
    Compare the species KINETIC list with DATAFILE.
    Get index of phases and species.
                                                    DL 03.2010
**************************************************************/
bool REACT_CAP::CmpSpeciesListKIN(void)
{
	bool sr = true;
	LI ipc, noerr = 0;
	int i, j;
	// CRFProcess* m_pcs = NULL;

	cout << "\n"
	     << " => CmpSpeciesListKIN() "
	     << "\n"; // compare KIN LIST with DATAFILE and get list_idx
	for (i = 0; i < (int)this->species_kin_name.size(); i++)
	{
		if (this->species_kin_phase[i] == 0)
			CAP_tqinp(this->species_kin_name[i], &ipc, &noerr);
		else
			CAP_tqinpc(this->species_kin_name[i], this->species_kin_phase[i], &ipc, &noerr);
		if (noerr)
		{
			cout << " ERROR !!! SPECIES KIN INPUT out of DATAFILE. STOP."
			     << "\n";
			exit(1);
		} // STOP and EXIT
		this->species_kin_idx.push_back(ipc);
		cout << " " << setw(4) << right << this->species_kin_phase[i];
		cout << setw(4) << this->species_kin_idx[i] << "  ";
		cout << setw(10) << left << this->species_kin_name[i] << "\n";
	}

	for (i = 0; i < (int)this->species_name.size(); i++)
		for (j = 0; j < (int)this->species_kin_name.size(); j++)
			if (strcmp(species_name[i], species_kin_name[j]) == 0 && species_phase[i] == species_kin_phase[j])
				species_dormant[i] = 1;

	cout << " KINETIC SPECIES DEFINITION OK ..."
	     << "\n";
	return sr;
}

// DL 04.2013

bool REACT_CAP::CmpSpeciesListRelativeActivity(void)
{
	bool sr = true, is_pcs;
	LI ipc, noerr = 0;
	int i, j;
	// CRFProcess* m_pcs = NULL;
	cout << "\n"
	     << " => CmpSpeciesListRelativeActivity() "
	     << "\n";

	for (i = 0; i < (int)this->species_relative_activity_name.size(); i++)
	{
		CAP_tqinp(this->species_relative_activity_name[i], &ipc, &noerr);
		if (noerr)
		{
			cout << " ERROR !!! SPECIES RELATIVE ACTIVITY INPUT out of DATAFILE. STOP."
			     << "\n";
			exit(1);
		} // STOP and EXIT
		this->species_relative_activity_idx.push_back(ipc);

		cout << " " << setw(8) << left << this->species_relative_activity_idx[i];
		cout << setw(10) << left << this->species_relative_activity_name[i];
		cout << "  " << right << this->species_relative_activity[i] << "\n";

		is_pcs = false;
		for (j = 0; j < (int)species_name.size(); j++)
		{
			if (strcmp(species_name[j], species_relative_activity_name[i]) == 0)
			{
				this->species_relative_activity_idx_pcs.push_back(j);
				is_pcs = true;
			}
		}
		if (!is_pcs)
		{
			cout << " ERROR !!! SPECIES RELATIVE ACTIVITY INPUT out of PCS-MT list. STOP."
			     << "\n";
			exit(1);
		} // STOP and EXIT
	}
	return sr;
}

/**************************************************************
    Function: REACT_CAP::CmpSpeciesListKinReact

    Task:
    Compare the species Kin React list with DATAFILE.
    Get index of phases and species.
                                                    DL 03.2010
**************************************************************/
bool REACT_CAP::CmpSpeciesListKinReact(void)
{
	bool sr = true;
	LI ipc = 0, noerr = 0;
	int i, ii;
	// CRFProcess* m_pcs = NULL;

	cout << "\n"
	     << " => CmpSpeciesListKinReact() "
	     << "\n"; // compare Kin React LIST with DATAFILE and get list_idx
	for (ii = 0; ii < (int)this->Kin_Reactions.size(); ii++)
	{
		cout << " Reaction " << Kin_Reactions[ii].type << "\n";
		for (i = 0; i < (int)this->Kin_Reactions[ii].species_name.size(); i++)
		{
			if (this->Kin_Reactions[ii].species_phase[i] == 0)
				CAP_tqinp(this->Kin_Reactions[ii].species_name[i], &ipc, &noerr);
			else
				CAP_tqinpc(
				    this->Kin_Reactions[ii].species_name[i], this->Kin_Reactions[ii].species_phase[i], &ipc, &noerr);
			if (noerr)
			{
				cout << " ERROR !!! SPECIES LOGK INPUT out of DATAFILE. STOP."
				     << "\n";
				exit(1);
			} // STOP and EXIT
			this->Kin_Reactions[ii].species_idx[i] = ipc;
			cout << " " << setw(4) << right << this->Kin_Reactions[ii].species_phase[i];
			cout << setw(4) << this->Kin_Reactions[ii].species_idx[i];
			cout << setw(4) << this->Kin_Reactions[ii].species_stoi[i] << "  ";
			cout << setw(10) << left << this->Kin_Reactions[ii].species_name[i] << "\n";
		}
	}
	cout << " KINETIC REACTIONS DEFINITION OK ..."
	     << "\n";
	return sr;
}

/**************************************************************
    Function: REACT_CAP::UndefOutSpecies
    Task:
    Search and find the species in PCS list for output of value.
                                                 DL 10.09.2009
**************************************************************/
void REACT_CAP::UndefOutSpecies(void)
{
	LI ipc = 0, noerr = 0;
	int ii; // number of process and transport species
	CRFProcess* m_pcs = NULL;

	for (ii = this->mass_num; ii < (int)pcs_mass_idx.size(); ii++)
	{
		m_pcs = pcs_vector[pcs_mass_idx[ii]];
		CAP_tqinpc((char*)m_pcs->pcs_primary_function_name[0], 2, &ipc, &noerr);
		if (this->species_mobil[ii] == 0 && (!noerr))
		{
			cout << "     --> " << m_pcs->pcs_primary_function_name[0] << "\n";
			pcs_ospecies_idx.push_back(pcs_mass_idx[ii]);
		}
	}
}

/**************************************************************
    Function: REACT_CAP::CmpSpeciesList

    Task:
    Compare the two species lists from MASS_TRANSPORT and .cap
    respectively, find out incompatible species.
                                                    DL 12.2008
**************************************************************/
bool REACT_CAP::CmpSpeciesList(void)
{
	int i, j, ns, nls, sf;
	bool sr = true;
	ns = (int)this->species_name.size();
	nls = (int)this->species_list_name.size();

	cout << "\n"
	     << "\n"
	     << " => CmpSpeciesList() "
	     << "\n";
	for (i = 0; i < ns; i++)
	{
		sf = 0; // flag sf->0 can not find
		for (j = 0; j < nls; j++)
			if (strcmp(species_name[i], species_list_name[j]) == 0)
				sf = 1;
		if (sf != 1)
		{
			cout << " WARNING... " << species_name[i] << " is incompatible with SPECIES LIST in .cap"
			     << "\n";
			sr = false;
		}
	}
	return sr;
}

/**************************************************************
    Function: REACT_CAP::DefSpeciesInChemApp

    Task:
    Eliminate redundant species from DataBase for reaction
                                                    DL 12.2008
**************************************************************/
void REACT_CAP::DefSpeciesInChemApp(void)
{
	LI i, ii, ns, nphase = 0, npcon = 0, noerr = 0; //, numcon=0;
	// cout << "\n" << "\n" << "\n" << " => DefSpeciesInChemApp() " << "\n";
	ns = (LI)species_name.size();
	CAP_tqnop(&nphase, &noerr); // get the number of phase --> nphase
	if (!noerr)
	{
		for (i = 1; i <= nphase; i++)
		{
			CAP_tqnopc(i, &npcon, &noerr);
			if (npcon == 1)
				CAP_tqcsp(i, (char*)"eliminated", &noerr);
			else
				for (ii = 1; ii <= npcon; ii++)
					CAP_tqcspc(i, ii, (char*)"eliminated", &noerr);
		}
		for (i = 0; i < ns; i++)
		{
			if (species_phase[i] == 0)
				CAP_tqcsp(species_idx[i], (char*)"entered", &noerr);
			if (species_phase[i] > 0)
				CAP_tqcspc(species_phase[i], species_idx[i], (char*)"entered", &noerr);
		}
		for (i = 0; i < (int)species_kin_name.size(); i++)
		{
			if (species_kin_phase[i] == 0)
				CAP_tqcsp(species_kin_idx[i], (char*)"dormant", &noerr);
			if (species_kin_phase[i] > 0)
				CAP_tqcspc(species_kin_phase[i], species_kin_idx[i], (char*)"dormant", &noerr);
		}
	}
}

void REACT_CAP::DefSpeciesInChemAppList(void)
{
	LI i, ii, ns, nphase = 0, npcon = 0, noerr = 0;
	// cout << "\n" << "\n" << " => DefSpeciesInChemAppList() " << "\n";
	ns = (LI)species_list_name.size();
	CAP_tqnop(&nphase, &noerr); // get the number of phase --> nphase
	for (i = 1; i <= nphase; i++)
	{
		CAP_tqnopc(i, &npcon, &noerr);
		if (npcon == 1)
			CAP_tqcsp(i, (char*)"eliminated", &noerr);
		else
			for (ii = 1; ii <= npcon; ii++)
				CAP_tqcspc(i, ii, (char*)"eliminated", &noerr);
	}
	for (i = 0; i < ns; i++)
	{
		if (species_list_phase[i] == 0)
			CAP_tqcsp(species_list_idx[i], (char*)"entered", &noerr);
		if (species_list_phase[i] > 0)
			CAP_tqcspc(species_list_phase[i], species_list_idx[i], (char*)"entered", &noerr);
	}
	for (i = 0; i < (int)species_kin_name.size(); i++)
	{
		if (species_kin_phase[i] == 0)
			CAP_tqcsp(species_kin_idx[i], (char*)"dormant", &noerr);
		if (species_kin_phase[i] > 0)
			CAP_tqcspc(species_kin_phase[i], species_kin_idx[i], (char*)"dormant", &noerr);
	}
}

void REACT_CAP::DefSpeciesInChemAppAll(void)
{
	int i, ii;
	LI nphase, npcon, noerr;
	CAP_tqnop(&nphase, &noerr);
	for (i = 1; i <= nphase; i++)
	{
		CAP_tqnopc(i, &npcon, &noerr);
		if (npcon == 1)
			CAP_tqcsp(i, (char*)"entered", &noerr);
		else
			for (ii = 1; ii <= npcon; ii++)
				CAP_tqcspc(i, ii, (char*)"entered", &noerr);
	}
}

void REACT_CAP::SetAllSolidAsDormant(void)
{
	LI i, nphase = 0, npcon = 0, noerr = 0;
	// LI ii, ns ;
	// cout  << "\n"  << "\n" << " => SetAllSolidAsDormant() " << "\n";
	CAP_tqnop(&nphase, &noerr); // get the number of phase --> nphase
	for (i = 1; i <= nphase; i++)
	{
		CAP_tqnopc(i, &npcon, &noerr);
		if (npcon == 1)
			CAP_tqcsp(i, (char*)"dormant", &noerr);
	}
}

/***********************************************************
      Function: REACT_CAP::ResetElement

      Task:
      Change Elements as Component for ChemApp
                                             DL   12/08
***********************************************************/
void REACT_CAP::ResetElement(void)
{
	int i, j, ipw, ns, nm = 0;
	LI noerr, nscom;
	DB *stoi, *zstoi, wmass;
	typedef char charsc[TQSTRLEN];
	charsc *newsc, scname;
	vector<std::string> newscv;

	cout << "\n"
	     << "\n"
	     << " => ResetElement() "
	     << "\n";
	cout << " RESET ELEMENT..."
	     << "\n";
	ipw = 0;
	newscv.clear();
	ns = (int)species_name.size();
	for (i = 0; i < ns; i++)
		if (species_phase[i] == 2 && species_mobil[i] == 1)
			ipw++;
	CAP_tqnosc(&nscom, &noerr);
	if (nscom < ipw)
	{
		cout << " ERROR... DEFINED TOO MANY Mass Transport Elements "
		     << "\n";
		exit(1);
	}
	// if(nscom==ipw) // mass transport element composition is same with that of DATABASE
	//	for(i=0;i<ns;i++)
	//	    if(species_phase[i]==2)	newscv.push_back(species_name[i]);

	if (nscom >= ipw)
	{
		stoi = new DB[nscom];
		zstoi = new DB[nscom];
		for (i = 0; i < nscom; i++)
			zstoi[i] = 0;
		for (i = 0; i < ns; i++)
			if (strcmp(species_name[i], "EA") != 0 && species_phase[i] == 2 && species_mobil[i] == 1)
			{
				newscv.push_back(species_name[i]);
				CAP_tqstpc(species_phase[i], species_idx[i], stoi, &wmass, &noerr);
				for (j = 0; j < nscom; j++)
					zstoi[j] = zstoi[j] + stoi[j]; // find which element out of new element-list
			}

		for (i = 0; i < nscom; i++) // sum the account of other elements
			if (zstoi[i] == 0)
				nm++;

		if (newscv.size() + nm + 1 == (size_t)nscom)
			newscv.push_back("EA"); // add electron as first element

		for (i = 0; i < nscom; i++) // add previous other elements into new element-system
			if (zstoi[i] == 0)
			{
				CAP_tqgnsc(i + 1, scname, &noerr);
				newscv.push_back(scname);
			}

		if (newscv.size() != (size_t)nscom)
		{
			cout << " ERROR... DEFINED Elements with potential CONFLICT."
			     << "\n";
			exit(1);
		}
	}
	cout << " OK "
	     << "\n";
	newsc = new charsc[nscom];
	for (i = 0; i < nscom; i++)
		strcpy(newsc[i], newscv[i].c_str());
	strcpy(newsc[nscom], ""); // statement for end of element list
	cout << " New Element"
	     << "\n";
	for (i = 0; i < nscom; i++)
		cout << " " << setw(4) << i << newsc[i] << "\n";

	CAP_tqcsc((CHP)newsc, &noerr);
	if (!noerr)
		cout << " ELEMENT INPUT OK ..."
		     << "\n"
		     << "\n";
	else
	{
		cout << " ELEMENT INPUT ERROR !!! STOP."
		     << "\n"
		     << "\n";
		exit(1);
	} // STOP and EXIT
}

/***********************************************************
      Function: REACT_CAP::LoopNodeReact

      Task:
      Chemical Reaction Calculation at all nodes
                                             DL   12/08
***********************************************************/
void REACT_CAP::LoopNodeReact(int f, int nodeflag)
{
	int ff;
	int ii, ns, isc, ix, widx = 0, fg; // idx,
	size_t i, iv;

	LI ipc, iEA, noerr, numcon = 0, nscom = 0;
	DB TT = 298.15, PP = 1.0, vals[2], value = 0.0, value1 = 0.0, *stoi,
	   wmass; //, zEA; // zEA -> sum of electronic charge in solution
	CRFProcess* m_pcs = NULL;

	vector<double> species_value_d, kin_value, residual_value;
	vector<double> species_value, species_value_b; // back up the value  DL 2012.2.12

	vector<int> residual_idx;
	bool first_time = true;
	// bool is_2nd_try;
	int number_UnderSat;

	bool MinKinReact = KMinKinCheck(); //  CB 3-12-2010 Check, if KinReact module is used
	// double mv[8], cv_CO2, cv[8], mCO2, dens ;
	double unitfactor_l = 1, unitfactor_s = 1;
	// CB / DL 21.1.2011: Attention: Gas phase volume fraction : (1-S)*n

	stringstream ss;
	// int i_r, i_s;
	int i_np; // for residual value return back pcs
	int np, precision_type;
	double sp_value;

	double delta, err_vle = 1.0e-12; // VLE
	bool is_equilibrium;
	bool is_VLE;
	int iter_eq;
	vector<bool> iv_eq;
	vector<double> a, b, x;
	vector<double> species_value_s;

	ff = f; // if ff=0-->initial calculation,  ff=1->full geochemical system, ff=-1->for liquid system
	fg = f;
	f = 1; // old time level or new time level

	cout << " ChemApp Reaction node loop."
	     << "\n";
	cout.flush();

	// CB 19.1.2011
	// Get the reaction interface data
	REACTINT* m_rei = NULL;
	if (REACTINT_vec.size() > 0)
	{
		m_rei = REACTINT_vec[0];
		if (ff == 0 && m_rei->icSolidUpdate)
			ff = 1; // allow update of solid matrix concentrations
	}

	ns = (int)species_name.size();

	CAP_tqinsc((char*)"EA", &iEA, &noerr);
	CAP_tqnosc(&nscom, &noerr);
	stoi = new DB[nscom];
	// cout << " iEA= " << iEA << " " << nscom << "\n";

	// if(ff != -1){
	if (nodeflag < 0)
	{
		node_logK.clear();
		node_ac.clear();
		node_HKF_logK.clear();
	}
	//}
	CAP_tqcio((char*)"ERROR", 0, &noerr);

	// node loop
	for (ii = 0; ii < this->nodenumber; ii++)
	{ // ii==0 as boundary point without reaction calc

		CAP_icount = 0;
		CAP_Time = aktueller_zeitschritt;
		CAP_Node = ii;

		this->RecoverChemSystem();

		if (nodeflag >= 0)
			if (nodeflag != ii)
				continue;

		CAP_tqremc(0, &noerr); // removes all input conditions relating to incoming amounts

		// return pcs rename pre
		// in new version, set in problem.cpp ?
		for (i = 0; i < pcs_rename_idx0_pre.size(); i++)
		{
			value = 0.0;
			for (ix = 0; ix < (int)pcs_rename_idx1_pre[i].size(); ix++)
			{
				m_pcs = pcs_vector[pcs_rename_idx1_pre[i][ix]];
				if (pcs_rename_stoi_pre[i][ix] == 999999)
					value *= m_pcs->GetNodeValue(ii, f);
				else if (pcs_rename_stoi_pre[i][ix] == -999999)
					value /= m_pcs->GetNodeValue(ii, f);
				else
					value += m_pcs->GetNodeValue(ii, f) * pcs_rename_stoi_pre[i][ix];
			}
			m_pcs = pcs_vector[pcs_rename_idx0_pre[i]];
			m_pcs->SetNodeValue(ii, f, value);
		}

		// only for nodes, where Cnew != Cold
		if (rateflag[ii] != 0)
		{
			// CAP_tqremc(-2, &noerr); // remove all condition and targets set previously
			// CAP_tqstrm("inputs", &noerr); // remove stream input

			// get species conc. after transport mol/m2 Liquid
			species_value.clear();
			for (i = 0; i < (size_t)ns; i++)
			{
				m_pcs = pcs_vector[this->pcs_mass_idx[i]];
				species_value.push_back(m_pcs->GetNodeValue(ii, f));
				// if(ii==1) cout << i << " " << species_name[i] << " " << species_value[i] << "\n";
			}

			// CB 19.1.2011
			// based on porosity, calculate TOTALS Ti,w Tj,s before coputing equilirium chemistry
			// mol (/m3aquifer)
			if (m_rei)
			{
				if (m_rei->unitconversion)
				{
					m_rei->CalcUnitConversionFactors(ii, &unitfactor_l, &unitfactor_s, false);
					// unitfactor_l = m_rei->node_porosity[ii] * m_rei->GetWaterSaturation(ii);
					// if(unitfactor_l==0.0) unitfactor_l = m_rei->node_porosity[ii]* 1;
					//  unitfactor_s = 1 - m_rei->node_porosity[ii];
					for (i = 0; i < (size_t)ns; i++)
					{
						// idx = pcs_vector[pcs_mass_idx[i]]->GetProcessComponentNumber();
						if (species_phase[i] == 2)
						{ // liquid phase
							if (strcmp(species_name[i], "H2O") == 0 || strcmp(species_name[i], "H2O_liquid") == 0
							    || strcmp(species_name[i], "water_liquid") == 0)
							{
								// the water pcs node value is now updated after preprocessing,
								// and after kinreact anyway, so use the node value
								// species_value[i] =  m_rei->water_conc[ii] * unitfactor_l ; // set the total amount of
								// water
								widx = i; // save the species index of watre species
							}
							// else
							species_value[i] *= unitfactor_l; // Ti,w = Ci,w * n * S
						}
						else if (species_phase[i] == 0) // solid phase
							species_value[i] *= unitfactor_s; // Tj,s = Cj,s * (1-n)
						// CB / DL 21.1.2011: Attention: Gas phase volume fraction : (1-S)*n
						// else if (species_phase[i]==1)// gas phase
					}
				}
				// set P, T input values
				TT = m_rei->GetTemperature(ii);
				PP = m_rei->GetPressure(ii);
			}

			// for(i=0;i<ns;i++)	if(ii==1) cout << i << " " << species_name[i] << " " << species_value[i] << "\n";

			//----VLE--init--
			species_value_s.clear();
			for (i = 0; i < species_value.size(); i++)
				species_value_s.push_back(species_value[i]);
			iv_eq.clear();
			a.clear();
			b.clear();
			x.clear();

			if (m_rei)
				for (iv = 0; iv < m_rei->VLE_conditions.size(); iv++)
				{
					iv_eq.push_back(false);
					a.push_back(0);
					b.push_back(0);
					x.push_back(0);
					m_pcs = pcs_vector[m_rei->VLE_conditions[iv].vp_idx];
					m_rei->VLE_conditions[iv].vp_value = m_pcs->GetNodeValue(ii, f);
					for (i = 0; i < (size_t)ns; i++)
						if (this->pcs_mass_idx[i] == m_rei->VLE_conditions[iv].aq_idx)
						{
							m_rei->VLE_conditions[iv].idx_aq_species = i;
							m_rei->VLE_conditions[iv].aq_value = species_value[i];
						}
				}

			//----END--init--

			//----VLE_P--init----
			if (m_rei)
				for (iv = 0; iv < m_rei->VLE_pressure.size(); iv++)
				{
					iv_eq.push_back(false);
					a.push_back(0);
					b.push_back(0);
					x.push_back(0);
					m_pcs = pcs_vector[m_rei->VLE_pressure[iv].vp_idx];
					m_rei->VLE_pressure[iv].vp_value = m_pcs->GetNodeValue(
					    ii, f); // vp_value is the partial pressure of the gas, from current node value (IC or BC ??)
					for (i = 0; i < (size_t)ns; i++)
						if (this->pcs_mass_idx[i] == m_rei->VLE_pressure[iv].aq_idx)
						{
							m_rei->VLE_pressure[iv].idx_aq_species = i;
							// m_rei->VLE_pressure[iv].aq_value=species_value[i];
							// m_rei->VLE_pressure[iv].aq_value=   m_rei->VLE_pressure[iv].vp_value; //
							// species_value[i]; //should get value from geochemcalc
						}
				}

			iter_eq = 0;
			is_equilibrium = false;
			// record input species_values

			while (!is_equilibrium)
			{
				//====RECOVERING LOOP START====

				for (i_np = -1; i_np < 4; i_np++)
				{ // ChemApp eq calc recovering

					// set P, T input values
					// TT = 273.15;
					// PP = 100.0;
					// CALL ChemAppCalc(species_name, species_value, T, P)
					CAP_tqsetc((char*)"T", 0, 0, TT, &numcon, &noerr);
					CAP_tqsetc((char*)"P", 0, 0, PP, &numcon, &noerr);
					// cout << ii << " T " << TT << " P " << PP << "\n";

					// for(i=0;i<ns;i++)
					//	  cout << i << " " << species_value[i] << "\n";
					// setting input concentration values for Chemapp

					// zEA=0;
					residual_value.clear();
					residual_idx.clear();

					for (i = 0; i < (size_t)ns; i++)
					{
						precision_type = 0;
						if (species_dormant[i] == 0)
						{ // dormant = 1 --> kinetic species, does not take part in eq reactions

							if (strcmp(species_name[i], "EA") == 0)
							{
								CAP_tqinsc((char*)"EA", &iEA, &noerr); // for electrons, e-
								CAP_tqsetc((char*)"ia", 0, iEA, species_value[i], &numcon, &noerr);
							}
							else
							{
								if (species_phase[i] == 0 && species_value[i] < 0.0e0 && species_value[i] > -50)
								{ // DL to set different solid composition using IC flag "-100"
									residual_value.push_back(species_value[i]);
									residual_idx.push_back(i);
									species_value[i] = 0; // solid, negative c
								}

								if (m_rei)
								{ // set precision type for old version, in new version it is no used
									if (m_rei->unitconversion)
									{
										if (species_phase[i] != 0 && species_value[i] < 1.0e-7 && species_value[i] != 0)
											precision_type = 1;
										// H and O element
										if (species_value[i] < 1.0e-5 && species_value[i] != 0)
											if (m_rei->formula2index(species_name[i])[8] != 0
											    || m_rei->formula2index(species_name[i])[9] != 0)
												precision_type = 1;
									}
									else
									{
										if (species_phase[i] != 0 && species_value[i] < 1.0e-10
										    && species_value[i] != 0)
											precision_type = 1;
										// H and O element
										if (species_value[i] < 1.0e-8 && species_value[i] != 0)
											if (m_rei->formula2index(species_name[i])[8] != 0
											    || m_rei->formula2index(species_name[i])[9] != 0)
												precision_type = 1;
									}
								}
								else
								{
									if (species_phase[i] != 0 && species_value[i] < 1.0e-10 && species_value[i] != 0)
										precision_type = 1;
									// H and O element
									if (species_value[i] < 1.0e-8 && species_value[i] != 0)
										if (m_rei->formula2index(species_name[i])[8] != 0
										    || m_rei->formula2index(species_name[i])[9] != 0)
											precision_type = 1;
								}

								// new version for precision setting
								if (species_phase[i] == 0)
								{
									if (ff != -1) // for full system, ff==-1 for the liquid system
										if (species_value[i] > -50)
										{
											CAP_tqsetc((char*)"ia",
											           species_idx[i],
											           0,
											           species_value[i],
											           &numcon,
											           &noerr); // single comp phase (solid)
										}
								}
								else
								{
									if (i_np == -1)
										sp_value = species_value[i];
									else if ((precision_type == 0 || precision_type == 1) && i_np >= 0 && i_np < 3)
									{
										if (species_value[i] > 0)
										{
											np = (int)log10(species_value[i]);
											if (np > 0)
												np = 7;
											else if (np > -2)
												np = 6;
											else if (np > -4)
												np = 5;
											else if (np > -5)
												np = 4;
											else if (np <= -5 && np >= -7)
												np = 3;
											else if (np < -7 && np > -16)
												np = 2;
											else if (np <= -16)
												np = 1;

											np -= i_np;
											if (np < 0)
												np = 0;
											ss.clear();
											ss << setprecision(np) << species_value[i];
											ss >> sp_value;
										}
										else
											sp_value = 0;
									}
									else if (precision_type == 1 && i_np == 3)
										sp_value = 0;
									else
										sp_value = species_value[i];

									residual_value.push_back(species_value[i] - sp_value);
									residual_idx.push_back(i);
									CAP_tqsetc((char*)"ia",
									           species_phase[i],
									           species_idx[i],
									           sp_value,
									           &numcon,
									           &noerr); // any mixture phase
									CAP_tqstpc(species_phase[i], species_idx[i], stoi, &wmass, &noerr);
									// zEA=zEA+stoi[iEA-1]*species_value[i];		//set charge balance automatically
								}
							} // if strcmp
						} // if dormant
					} // for ns

					// cout << " zEA= " << zEA << "\n";
					// CAP_tqsetc("ia", 0, iEA, zEA, &numcon, &noerr); // set electronic charge balance
					// if(zEA>0) CAP_tqsetc("ia", 2, 34, zEA, &numcon, &noerr);
					// if(zEA<0) CAP_tqsetc("ia", 2, 71, -zEA, &numcon, &noerr);
					// cout  << "\n" << ii << "-----------------------------"  << "\n" << "Input      P Id   Name" <<
					// "\n";
					// if(ii==2)
					// for(i=0;i<ns;i++){
					// cout << setprecision(20) << species_value[i] << " " << 0 << " " << species_phase[i] << " " <<
					// setw(4) << species_idx[i] << " " << species_name[i] << "\n";
					//}
					// CAP_tqshow(&noerr);

					if (this->species_relative_activity_name.size() == 0)
					// if(!this->species_relative_activity_name.size()>0)
					{ // no relative activity

						for (i = 0; i < (size_t)ns; i++)
						{ // DL set the solid species as eliminated when the value is -100
							if (species_phase[i] == 0 && species_value[i] < -50 && species_value[i] > -200)
								CAP_tqcsp(this->species_idx[i], (char*)"eliminated", &noerr);
						}
						// here, calculate equilibrium geochemistry for a node
						noerr = 0;
						if (first_time && fg == 0)
						{ // do when the first node is called
							CAP_tqce(const_cast<char*>(" "), 0, 0, vals, &noerr);
							first_time = false;
						}
						else
						{
							if (ii == 1 /*|| ii==1 || ii==8 || ii==12 || ii==16 || ii==20 */)
								CAP_tqce(const_cast<char*>(" "),
								         0,
								         0,
								         vals,
								         &noerr); // this is faster, using previuos result as start for iteration
							else
								CAP_tqce(const_cast<char*>(" "),
								         0,
								         0,
								         vals,
								         &noerr); // this is faster, using previuos result as start for iteration
						}
					}
					// DL ---------------relative activity---------------start
					else if (this->species_relative_activity_name.size() > 0)
					{
						this->species_relative_activity_state.clear();
						for (i = 0; i < species_relative_activity_name.size(); i++)
							species_relative_activity_state.push_back(0);
						for (ix = 0; ix < 5; ix++)
						{
							// pre setting
							for (i = 0; i < species_relative_activity_name.size(); i++)
							{
								if (species_value[species_relative_activity_idx_pcs[i]] > -50)
								{
									if (species_value[species_relative_activity_idx_pcs[i]] > 1.0e-16)
									{ // precision
										if (species_relative_activity_state[i] == 1
										    || species_relative_activity_state[i] == 0)
										{
											CAP_tqcsp(species_relative_activity_idx[i], (char*)"entered", &noerr); // DL
											CAP_tqsetc((char*)"ac",
											           species_relative_activity_idx[i],
											           0,
											           species_relative_activity[i],
											           &numcon,
											           &noerr); // DL
											species_relative_activity_state[i] = 1;
										}
										if (species_relative_activity_state[i] == 2)
										{
											CAP_tqcsp(species_relative_activity_idx[i], (char*)"dormant", &noerr); // DL
											CAP_tqsetc((char*)"ia",
											           species_relative_activity_idx[i],
											           0,
											           species_value[this->species_relative_activity_idx_pcs[i]],
											           &numcon,
											           &noerr);
										}
									}
									else if (species_value[this->species_relative_activity_idx_pcs[i]] <= 1.0e-16)
									{
										if (this->species_relative_activity_state[i] == 3
										    || this->species_relative_activity_state[i] == 0)
										{
											CAP_tqcsp(
											    this->species_relative_activity_idx[i], (char*)"dormant", &noerr); // DL
											CAP_tqsetc((char*)"ia",
											           this->species_relative_activity_idx[i],
											           0,
											           0.0,
											           &numcon,
											           &noerr);
											this->species_relative_activity_state[i] = 3;
										}
										if (this->species_relative_activity_state[i] == 4)
										{
											CAP_tqcsp(
											    this->species_relative_activity_idx[i], (char*)"entered", &noerr); // DL
											CAP_tqsetc((char*)"ac",
											           this->species_relative_activity_idx[i],
											           0,
											           this->species_relative_activity[i],
											           &numcon,
											           &noerr); // DL
										}
									}
								}
							}

							for (i = 0; i < (size_t)ns; i++)
							{ // DL set the solid species as eliminated when the value is -100
								if (species_phase[i] == 0 && species_value[i] < -50 && species_value[i] > -200)
									CAP_tqcsp(this->species_idx[i], (char*)"eliminated", &noerr);
							}
							// CAP_tqshow(&noerr);

							// calculate
							noerr = 0;
							if (first_time && fg == 0)
							{ // do when the first node is called
								CAP_tqce(const_cast<char*>(" "), 0, 0, vals, &noerr);
								first_time = false;
							}
							else
							{
								if (ii == 1 || ii == 6 /*|| ii==8 || ii==12 || ii==16 || ii==20 */)
									CAP_tqcen((char*)" ",
									          0,
									          0,
									          vals,
									          &noerr); // this is faster, using previuos result as start for iteration
								else
									CAP_tqcen((char*)" ",
									          0,
									          0,
									          vals,
									          &noerr); // this is faster, using previuos result as start for iteration
							}

							// post checking
							if (noerr)
								ix = 100;
							else
							{
								this->species_relative_activity_ia.clear();
								this->species_relative_activity_calc.clear();
								number_UnderSat = 0;
								for (i = 0; i < species_relative_activity_name.size(); i++)
								{
									if (species_value[this->species_relative_activity_idx_pcs[i]] > -50)
									{ // DL====
										CAP_tqgetr((char*)"ia",
										           this->species_relative_activity_idx[i],
										           0,
										           &value,
										           &noerr); // get income amount, "ia"
										this->species_relative_activity_ia.push_back(value);
										CAP_tqgetr((char*)"ac",
										           this->species_relative_activity_idx[i],
										           0,
										           &value,
										           &noerr); // get relative activity after calc, "ac"
										this->species_relative_activity_calc.push_back(value);
									}
									else
									{
										this->species_relative_activity_ia.push_back(-100);
										this->species_relative_activity_calc.push_back(-100);
									}
									if (species_value[this->species_relative_activity_idx_pcs[i]] > -50)
									{ // DL====
										if (this->species_relative_activity_state[i] == 1
										    || this->species_relative_activity_state[i] == 4)
										{ // ra case
											if (this->species_relative_activity_ia[i] > 0.0
											    && this->species_relative_activity_ia[i]
											           > species_value[this->species_relative_activity_idx_pcs[i]])
											{
												this->species_relative_activity_state[i] = 2;
												number_UnderSat++;
											}
											else
											{
												this->species_relative_activity_state[i] = 1;
												this->species_relative_activity_ia[i]
												    = species_value[this->species_relative_activity_idx_pcs[i]]
												      - this->species_relative_activity_ia[i];
											}
										}
										else if (this->species_relative_activity_state[i] == 2
										         || this->species_relative_activity_state[i] == 3)
										{ // dormant case
											if (this->species_relative_activity_calc[i]
											    > this->species_relative_activity[i])
											{
												this->species_relative_activity_state[i] = 4;
												number_UnderSat++;
											}
											else
											{
												this->species_relative_activity_state[i] = 3;
												this->species_relative_activity_ia[i] = 0.0;
											}
										}
										// cout << ix << " state " << this->species_relative_activity_state[i] << "\n";
									} // > -50
								}
								if (number_UnderSat == 0)
									ix = 100;
							}
						} // for ix< 5
					}
					// DL --------------------------------------------end ra

					if (noerr)
						warning_out << " Time Step " << aktueller_zeitschritt << "   Node " << ii
						            << ". Recovering ... i_np = " << i_np << " --> " << i_np + 1 << "\n";
					else
						break;
				} // end ChemApp Eq calc Recovering

				//====RECOVERING LOOP END====

				if (noerr)
				{
					cout << "\n"
					     << "\n"
					     << " Time Step " << aktueller_zeitschritt << "\n";
					cout << " WARNING... At node " << ii << " ChemApp ERROR... NO CHEMICAL REACTION CALCULATION!"
					     << "\n"; // exit(1)
					cout << ii << "-----------------------------"
					     << "\n"
					     << "         Input P   Id Name"
					     << "\n";
					for (i = 0; i < (size_t)ns; i++)
					{
						if (species_dormant[i] == 0)
							cout << setprecision(8) << setw(14) << species_value[i] << " " << species_phase[i] << " "
							     << setw(4) << species_idx[i] << " " << species_name[i] << "\n";
					}

					// CAP_tqshow(&noerr);
					// Kinetic Parameters Updata, in case of error, an empty vector will be stored
					if (MinKinReact && ff != 1 && nodeflag < 0)
					{
						// if(MinKinReact && ff != -1){
						species_value_d.clear();
						for (i = 0; i < species_value.size(); i++)
							species_value_d.push_back(species_value[i]);
						this->KinParamUpdata(ii, noerr, species_value_d);
						this->KinParamUpdataHKF(ii, TT, PP, noerr); // CB include noerr
					}
				}
				else
				{
					// show sateration index list
					if (fg == 0 && this->Node_All_SI.size() > 0)
					{
						for (i = 0; i < Node_All_SI.size(); i++)
						{
							if (ii == this->Node_All_SI[i])
							{
								cout << " => show all SI at Node : " << ii << "\n";
								this->SetAllSolidAsDormant();
								CAP_tqcel((char*)" ", 0, 0, vals, &noerr);
								this->RecoverChemSystem();
							}
						}
					}

					species_value_b.clear();
					for (i = 0; i < species_value.size(); i++)
						species_value_b.push_back(species_value[i]); // back up the input value
					// cout << "\n" << ii << "-----------------------------" << "\n" << "Output     P Id   Name" <<
					// "\n";
					isc = 0;
					// get results from Chemapp mol
					for (i = 0; i < (size_t)ns; i++)
					{
						if (species_value[i] > -50)
						{
							if (species_phase[i] == 0)
								CAP_tqgetr((char*)"a",
								           species_idx[i],
								           0,
								           &value,
								           &noerr); // a = amount, phase 0, value = C, error flag
							else
							{
								if (species_phase[i] == 2)
								{ // aqueous
									if (this->mass_type.compare("ELEMENT")) // element  //DL 2013
										CAP_tqgetr(
										    (char*)"a", 2, species_idx[i], &value, &noerr); // A Equilibrium amount
									else
									{ // species
										isc++;
										CAP_tqgetr((char*)"ap", 2, isc, &value, &noerr); // AP Equilibrium amount /
										// Equilibrium amount of system
										// component in a phase , isc
										// counter
									}
								}
								else // other phases
									CAP_tqgetr((char*)"a", species_phase[i], species_idx[i], &value, &noerr);
							}
							if (species_dormant[i] == 0)
								species_value[i] = value; // get eq value from chemapp, set to input vector

							// DL 04.2013 for relative activity ---------------------
							for (ix = 0; ix < (int)this->species_relative_activity_name.size(); ix++)
								if (this->species_relative_activity_idx_pcs[ix]
								    == (int)i) // found the pcs mt in relative activity list
									species_value[i] = this->species_relative_activity_ia[ix];
						}
						//---------------------------------------
						// if(ii==1)
						//	cout << setprecision(20) << species_value[i] << " " << isc << " " << species_phase[i] << " "
						//<< setw(4) << species_idx[i] << " " << species_name[i] << "\n";

					} // for ns

				} // end ChemApp calc

				//----VLE--calc-- direct iteration method
				// for(i=0;i<m_rei->VLE_conditions.size();i++){
				//	//cout << " VLE calc ... " << iter_eq << "\n";
				//	if(m_rei->VLE_conditions[i].vp_value > 1.0e-8 ){
				//		delta= m_rei->VLE_conditions[i].aq_value -
				// species_value[m_rei->VLE_conditions[i].idx_aq_species];
				//		delta *= 1000.0;
				//		if(abs(delta)>1.0e-8){
				//			if(delta>0)
				//				if(delta<=m_rei->VLE_conditions[i].vp_value){
				//					m_rei->VLE_conditions[i].vp_value -= delta;
				//					species_value[m_rei->VLE_conditions[i].idx_aq_species] += delta;
				//				}
				//				else{
				//					m_rei->VLE_conditions[i].vp_value =0;
				//					species_value[m_rei->VLE_conditions[i].idx_aq_species] +=
				// m_rei->VLE_conditions[i].vp_value;
				//				}
				//			else{
				//				m_rei->VLE_conditions[i].vp_value -= delta;
				//				species_value[m_rei->VLE_conditions[i].idx_aq_species] += delta;
				//			}
				//		}
				//		else
				//			iv_eq[i]=true;
				//	}
				//	else
				//		iv_eq[i]=true;
				//}

				//----bisection method--  VLE
				if (m_rei)
					if (m_rei->VLE_conditions.size() > 0)
					{
						for (i = 0; i < m_rei->VLE_conditions.size(); i++)
						{
							m_rei->VLE_conditions[i].delta = 0.0;
							if (m_rei->VLE_conditions[i].vp_value > 0)
							{
								delta = m_rei->VLE_conditions[i].aq_value
								        - species_value[m_rei->VLE_conditions[i].idx_aq_species];
								if (iter_eq == 0)
								{
									a[i] = delta;
									b[i] = 1.0; // 10.0*m_rei->VLE_conditions[i].aq_value;
								}
								if (abs(m_rei->VLE_conditions[i].aq_value) > 0)
									err_vle = 1.0e-6 * abs(m_rei->VLE_conditions[i].aq_value);
								else
									err_vle = 1.0e-12;
								if (abs(delta) < err_vle)
									iv_eq[i] = true;
								else
								{
									if (delta >= m_rei->VLE_conditions[i].vp_value)
									{
										m_rei->VLE_conditions[i].delta = m_rei->VLE_conditions[i].vp_value;
										// m_rei->VLE_conditions[i].vp_value = 0;
									}
									else
									{
										if (iter_eq > 0)
										{
											if (delta > 0)
												a[i] = x[i];
											else
												b[i] = x[i];
										}
										x[i] = 0.5 * (a[i] + b[i]);
										m_rei->VLE_conditions[i].delta = x[i];
									}
								}
							}
							else
								iv_eq[i] = true;
						}

						is_equilibrium = true;
						for (i = 0; i < m_rei->VLE_conditions.size(); i++)
						{
							if (!iv_eq[i])
								is_equilibrium = false;
						}
						iter_eq++;
						if (iter_eq >= 20)
						{
							is_equilibrium = true;
							cout << " warning..., at node " << ii << ", iter_eq max is 20. ";
						}
						if (!is_equilibrium)
						{
							for (i = 0; i < species_value.size(); i++)
								species_value[i] = species_value_s[i];
							for (i = 0; i < m_rei->VLE_conditions.size(); i++)
								species_value[m_rei->VLE_conditions[i].idx_aq_species]
								    += m_rei->VLE_conditions[i].delta;
						}
					}

				//----bisection method--  VLE pressure
				if (m_rei)
					if (m_rei->VLE_pressure.size() > 0)
					{
						for (i = 0; i < m_rei->VLE_pressure.size(); i++)
						{
							m_rei->VLE_pressure[i].aq_value = exp(VLE::LnPHI_CO2(TT, PP))
							                                  * m_rei->VLE_pressure[i].vp_value
							                                  / exp(VLE::Henry_const_CO2(TT));

							std::cout << m_rei->VLE_pressure[i].aq_value << "\n";
							std::cout << (VLE::LnPHI_CO2(TT, PP)) << "\n";
							std::cout << m_rei->VLE_pressure[i].vp_value << "\n";
							std::cout << VLE::Henry_const_CO2(TT) << "\n";

							// calc the real solubility of gas from geochemcalc, and store the value
							m_rei->VLE_pressure[i].delta = 0.0;
							if (m_rei->VLE_pressure[i].vp_value > 0)
							{
								delta = m_rei->VLE_pressure[i].aq_value
								        - species_value[m_rei->VLE_pressure[i].idx_aq_species]
								              / (species_value[widx] / 55.51);
								// compare with the new value after reaction
								if (iter_eq == 0)
								{
									a[i] = 0.0; // delta;
									// b[i]=2.0; //up limiting
									b[i] = 10.0 * species_value[widx] / 55.51;
								}
								if (abs(m_rei->VLE_pressure[i].aq_value) > 0)
									err_vle = 1.0e-3 * abs(m_rei->VLE_pressure[i].aq_value);
								else
									err_vle = 1.0e-6;
								if (abs(delta) < err_vle)
									iv_eq[i] = true;
								else
								{
									if (iter_eq > 0)
									{
										if (delta > 0)
											a[i] = x[i];
										else
											b[i] = x[i];
									}
									x[i] = 0.5 * (a[i] + b[i]);
									// cout << " aq " << m_rei->VLE_pressure[i].aq_value << " sp " <<
									// species_value[m_rei->VLE_pressure[i].idx_aq_species]/(species_value[widx]/55.51)
									// << " delta " << delta<< " x[i] " << x[i] << "\n";
									m_rei->VLE_pressure[i].delta = x[i];
								}
							}
							else
								iv_eq[i] = true;
						}
						is_equilibrium = true;
						for (i = 0; i < m_rei->VLE_pressure.size(); i++)
						{
							if (!iv_eq[i])
								is_equilibrium = false;
						}

						// cout << " iter_eq " << iter_eq << " water " << species_value[widx] << "\n";
						iter_eq++;
						if (iter_eq >= 20)
						{
							is_equilibrium = true;
							cout << " warning..., at node " << ii << ", iter_eq max is 20. for VLE_P ";
						}

						if (!is_equilibrium)
						{
							for (i = 0; i < species_value.size(); i++)
								species_value[i] = species_value_s[i];
							for (i = 0; i < m_rei->VLE_pressure.size(); i++)
								species_value[m_rei->VLE_pressure[i].idx_aq_species] += m_rei->VLE_pressure[i].delta;
						}
					}
				//-----------------

				if (m_rei)
				{
					if (m_rei->VLE_pressure.size() == 0 && m_rei->VLE_conditions.size() == 0)
						is_equilibrium = true;
				}
				else
					is_equilibrium = true;

			} // END  while(!is_equilibrium)
			//====Final Chemical Equilibria at the Node====

			is_VLE = false;

			if (m_rei)
			{
				if (m_rei->VLE_pressure.size() > 0)
				{
					for (i = 0; i < m_rei->VLE_pressure.size(); i++)
					{
						if (m_rei->VLE_pressure[i].vp_value > 0)
							is_VLE = true;
					}
				}
			}

			if (is_VLE)
				cout << "at node " << ii << " VLE_P iterations " << iter_eq << "\n";

			if (!noerr)
			{
				if (m_rei)
					// revise VLE vp_value
					for (i = 0; i < m_rei->VLE_conditions.size(); i++)
					{
						m_rei->VLE_conditions[i].vp_value -= m_rei->VLE_conditions[i].delta;
						m_pcs = pcs_vector[m_rei->VLE_conditions[i].vp_idx];
						m_pcs->SetNodeValue(ii, f, m_rei->VLE_conditions[i].vp_value);
					}

				// for(i=0;i<VLS_idx_aq.size();i++){
				//	delta_value_aq.push_back(0.0);
				//	delta_VLS_aq[i] =species_value_b[VLS_idx_aq[i]]-species_value[VLS_idx_aq[i]];
				//	if(abs( delta_VLS_aq[i]/(species_value[VLS_idx_aq[i]]+VLS_value_vs[i]) ) > 1.0e-8 )
				//}
				// pass parameters from pcs-eclipse

				////delete amount for charge balance
				// for(i=0;i<ns;i++){
				//	if(species_phase[i]==2){
				//		if(zEA >=0 && species_idx[i]==34)
				//			species_value[i] -= zEA;
				//		if(zEA <=0 && species_idx[i]==71)
				//			species_value[i] += zEA;
				//	}
				//}

				//-----------------------------------
				// Kinetic Parameters Updata
				if (MinKinReact && ff != 1 && nodeflag < 0)
				{
					// if(MinKinReact && ff != -1){

					species_value_d.clear();
					for (i = 0; i < species_value.size(); i++)
						species_value_d.push_back(species_value[i]);
					this->KinParamUpdata(
					    ii, noerr, species_value_d); // in case of  no error, a filled vector will be stored
					this->KinParamUpdataHKF(ii, TT, PP, noerr); // CB include noerr
					//-----------------------------------
					// species_value_d=this->KineticReact(species_value_d);
					//-----------------------------------
				}

				//-----------------------------------
				// return back the residual value
				if (m_rei)
					if (m_rei->residual_back)
						for (i = 0; i < residual_idx.size(); i++)
							species_value[residual_idx[i]] += residual_value[i];
				//-----------------------------------
				// species_value=this->KineticReact(species_value);
				//-----------------------------------

				// if(ii==1) cout << " value " << species_value.back() << "\n";
				// unit factors are only different from 1 if(m_rei && m_rei->unitconversion)
				for (i = 0; i < (size_t)ns; i++)
				{ // mol/m3
					if (species_dormant[ns - i - 1] == 0)
					{ // return only values of non-dormant species
						m_pcs = pcs_vector[pcs_mass_idx[ns - i - 1]];
						if (this->species_phase[ns - i - 1] == 0)
						{
							if (ff == 1) // full geochemistry calculation
								m_pcs->SetNodeValue(
								    ii, f, species_value.back() / unitfactor_s); // old f==0 or new f==1 time step
						}
						else if (this->species_phase[ns - i - 1] == 2)
						{
							m_pcs->SetNodeValue(ii, f, species_value.back() / unitfactor_l); // old or new time step
							if (ns - (int)i - 1 == widx && m_rei) // update water concentration vector as well
								if (m_rei->unitconversion)
									m_rei->water_conc[ii] = species_value.back() / unitfactor_l;
						}
						else // CB / DL 21.1.2011: Attention: Gas phase volume fraction : (1-S)*n
							m_pcs->SetNodeValue(ii, f, species_value.back()); // new time step
					}
					species_value.pop_back();
				} // push back chemical equilibrium values

				// TEST CB
				// species_value.clear();
				//	for(i=0;i<ns;i++){
				//		m_pcs = pcs_vector[this->pcs_mass_idx[i]];
				//		species_value.push_back(m_pcs->GetNodeValue(ii,1));
				//	}

				// return species value at this node and time step, if pcs_ospecies_idx.size>0  22.01.2009
				for (i = 0; i < pcs_ospecies_idx.size(); i++)
				{
					m_pcs = pcs_vector[pcs_ospecies_idx[i]];
					CAP_tqinpc((char*)m_pcs->pcs_primary_function_name[0], 2, &ipc, &noerr);
					CAP_tqgetr((char*)"a", 2, ipc, &value, &noerr);
					value /= unitfactor_l;
					m_pcs->SetNodeValue(ii, f, value);
				} // end return species value

				// return nlog value of species at this node, if pcs_nspecies_idx.size>0
				for (i = 0; i < pcs_nspecies_idx.size(); i++)
				{
					m_pcs = pcs_vector[pcs_nspecies_idx[i]];
					if (species_nlog_phase[nspecies_idx[i]] == 0)
						CAP_tqgetr((char*)"ac", species_nlog_idx[nspecies_idx[i]], 0, &value, &noerr);
					else
						CAP_tqgetr((char*)"ac",
						           species_nlog_phase[nspecies_idx[i]],
						           species_nlog_idx[nspecies_idx[i]],
						           &value,
						           &noerr);
					// if(value==0) exit(1);
					if (abs(value - 1.0) < 1.0e-8)
						value = 1.0;
					value = 0.0 - log10(value);
					m_pcs->SetNodeValue(ii, f, value);
				}
				// return redox Eh value of given reaction at this node
				if (pcs_redox > -1)
				{
					value1 = 0.0;
					for (i = 2; i < species_redox_name.size(); i++)
					{
						if (species_redox_phase[i] == 0)
							CAP_tqgetr((char*)"ac", species_redox_idx[i], 0, &value, &noerr);
						else
							CAP_tqgetr((char*)"ac", species_redox_phase[i], species_redox_idx[i], &value, &noerr);
						if (value > 0)
							value1 = value1 + redox_stoi[i] * log(value);
					}
					if (redox_stoi[1] != 0)
						value1 = redox_stoi[0] + 8.31451 * TT / 96485.3 * value1 / redox_stoi[1];
					else
						value1 = redox_stoi[0];
					m_pcs = pcs_vector[pcs_redox];
					m_pcs->SetNodeValue(ii, f, value1);
				}

			} // if(noerr)

		} // if(rateflag[ii]!=0)
		else if (MinKinReact && ff != 1 && nodeflag < 0)
		{ // Kinetic Parameters Updata				  //
			// when rateflag[ii] == 0, noerr is not defined. Do so for Kinetic Parameters Updating
			this->KinParamUpdata(ii, 1, species_value_d);
			this->KinParamUpdataHKF(ii, TT, PP, 1); // CB include noerr
		}

		// return pcs rename
		for (i = 0; i < pcs_rename_idx0.size(); i++)
		{
			value = 0.0;
			for (ix = 0; ix < (int)pcs_rename_idx1[i].size(); ix++)
			{
				//--return back residual
				// for(i_s=0;i_s<(int)pcs_mass_idx.size();i_s++){
				//	if(pcs_rename_idx1[i][ix]==pcs_mass_idx[i_s])
				//		for(i_r=0;i_r<(int)residual_value.size();i_r++){
				//			if(residual_idx[i_r]==i_s)
				//				value += residual_value[i_r];
				//		}
				//}

				m_pcs = pcs_vector[pcs_rename_idx1[i][ix]];
				if (pcs_rename_stoi[i][ix] == 999999)
					value *= m_pcs->GetNodeValue(ii, f);
				else if (pcs_rename_stoi[i][ix] == -999999)
					value /= m_pcs->GetNodeValue(ii, f);
				else
					value += m_pcs->GetNodeValue(ii, f) * pcs_rename_stoi[i][ix];
				;
			}
			m_pcs = pcs_vector[pcs_rename_idx0[i]];
			m_pcs->SetNodeValue(ii, f, value);
		}

	} // for (ii node loop)

} // end of function loopnodereact

void REACT_CAP::LoopNodeReact_Liquid_Vapor(int /*f*/, int nodeflag)
{
	LI i, noerr = 0;
	cout << "\n"
	     << " --> Liquid Vapor Reactions "
	     << "\n";

	for (i = 0; i < (int)species_name.size(); i++)
		if (species_phase[i] == 0)
			CAP_tqcsp(species_idx[i], (char*)"dormant", &noerr);

	this->LoopNodeReact(-1, nodeflag);

	for (i = 0; i < (int)species_name.size(); i++)
		if (species_phase[i] == 0)
			CAP_tqcsp(species_idx[i], (char*)"entered", &noerr);
	for (i = 0; i < (int)species_kin_name.size(); i++)
		if (species_kin_phase[i] == 0)
			CAP_tqcsp(species_kin_idx[i], (char*)"dormant", &noerr);
}

void REACT_CAP::LoopNodeReact_Liquid_Solid(int /*f*/, int /*nodeflag*/)
{
	// LI i, noerr=0;
	cout << "\n"
	     << " --> Liquid Solid Reactions "
	     << "\n";

	// for(i=0;i<(int)species_name.size();i++)
	//	if(species_phase[i]==0) CAP_tqcsp(species_idx[i], "dormant", &noerr);

	// this->LoopNodeReact(-1, nodeflag);

	// for(i=0;i<(int)species_name.size();i++)
	//	if(species_phase[i]==0) CAP_tqcsp(species_idx[i], "entered", &noerr);
	// for(i=0;i<(int)species_kin_name.size();i++)
	//	if(species_kin_phase[i]==0) CAP_tqcsp(species_kin_idx[i], "dormant", &noerr);
}

/***********************************************************
      Function: REACT_CAP::KinParamUpdat(void)

      Task:
      kinetic parameters updata at the ChemApp state
                                             DL   10/10
***********************************************************/
void REACT_CAP::KinParamUpdata(int /*ii*/, int err, vector<double>& spvc)
{
	int i, i_re, i_resp;
	long noerr = 0;
	double value_ac, value_dG, logK;
	vector<double> species_ac, reaction_logK;

	species_ac.clear();
	reaction_logK.clear();

	if (err == 0)
	{ // value was calculated
		for (i = 0; i < (int)species_name.size(); i++)
		{ // loop over species
			if (species_phase[i] == 0) // mineral phase
				species_ac.push_back(1.0);
			// CAP_tqgetr("ac", species_idx[i], 0, &value_ac, &noerr);	// this is activity!! not a.coeff
			else
			{
				CAP_tqgetr((char*)"ac", species_phase[i], species_idx[i], &value_ac, &noerr); // liquid phase
				if (strcmp(species_name[i], "H2O") == 0 || strcmp(species_name[i], "H2O_liquid") == 0
				    || strcmp(species_name[i], "water_liquid") == 0) // for H2O use value_ac directly
					species_ac.push_back(value_ac);
				else if (value_ac /*spvc[i]*/ > 0) // otherwise divide by molality to obtain gamma, if species conc > 0
					species_ac.push_back(value_ac / spvc[i] * spvc[0] / MOLH2OPERKG);
				else // store activity coefficient = 1, of species concentration = 0
					species_ac.push_back(1.0);
			}
		}

		for (i_re = 0; i_re < (int)this->Kin_Reactions.size(); i_re++)
		{ // loop over reactions
			logK = 0;
			for (i_resp = 0; i_resp < (int)this->Kin_Reactions[i_re].species_name.size(); i_resp++)
			{ // loop over species of reaction
				if (Kin_Reactions[i_re].species_phase[i_resp] == 0)
				{
					CAP_tqgdpc((char*)"G",
					           Kin_Reactions[i_re].species_idx[i_resp],
					           0,
					           &value_dG,
					           &noerr); // Get Gibbs energy for species
				}
				else
				{
					CAP_tqgdpc((char*)"G",
					           Kin_Reactions[i_re].species_phase[i_resp],
					           Kin_Reactions[i_re].species_idx[i_resp],
					           &value_dG,
					           &noerr);
				}
				logK += -value_dG * Kin_Reactions[i_re].species_stoi[i_resp] / 2.302585; // calculate logK for reaction
			}
			reaction_logK.push_back(logK); // store logK of a reaction
		}
		// node_logK.push_back(reaction_logK);
	}
	// pushback species_ac and reaction_logK vectors,
	// or empty vectors instead when noerr = 1 or rateflag == 0
	node_ac.push_back(species_ac);
	node_logK.push_back(reaction_logK);
}

/***********************************************************
      Function: REACT_CAP::KinParamUpdataHKF(void)

      Task:
      kinetic parameters updata HKF
                                             DL   11/10
***********************************************************/
void REACT_CAP::KinParamUpdataHKF(int /*ii*/, double T, double P, int err)
{
	int i = -1, i_re, i_resp;
	// long noerr = 0;
	// double value_ac;
	double value_dG, logK;
	vector<double> species_ac, reaction_logK;

	// cout << " KinParamUpdataHKF " << "\n";

	reaction_logK.clear();

	if (err == 0)
	{ // value was calculated
		HKFcalc(T, P);
		for (i_re = 0; i_re < (int)this->Kin_HKF_Reactions.size(); i_re++)
		{
			logK = 0;
			for (i_resp = 0; i_resp < (int)this->Kin_HKF_Reactions[i_re].species_name.size(); i_resp++)
			{
				i++;
				value_dG = Kin_HKF_species[Kin_HKF_index[i]].G;
				logK += Kin_HKF_Reactions[i_re].species_stoi[i_resp] * value_dG * 4.18 / 2.302585 / 8.314 / T;
			}
			reaction_logK.push_back(logK);
		}
		// node_HKF_logK.push_back(reaction_logK);
	}
	// pushback reaction logK vector,
	// or empty vector instead when noerr = 1 or rateflag == 0
	node_HKF_logK.push_back(reaction_logK);
}

/***********************************************************
      Function: REACT_CAP::KineticReact(int)

      Task:
      calc Kinetic reaction and updata the pcs_vector
                                             DL   10/10
***********************************************************/
vector<double> REACT_CAP::KineticReact(vector<double> species_value0)
{
	int i, iter_max = 5;
	double t_step, h_calc;
	vector<double> h, species_value1;
	vector<double> dcdt0;

	t_step = 0.01;

	h.clear();
	for (i = 0; i < (int)species_value0.size(); i++)
	{
		h_calc = 1 / (log10(species_value0[i]) + 5);

		h.push_back(h_calc + 10);
	}

	species_value1 = this->ODE(species_value0, h, iter_max, t_step);
	for (i = 0; i < (int)species_value0.size(); i++)
	{
		cout << setw(32) << species_name[i] << setw(16) << species_value0[i] << setw(10) << species_value1[i] << "\n";
	}

	// dcdt0=derivs(species_value0);
	// for(i=0;i<(int)species_value0.size();i++){
	//	cout << setw(32) << species_name[i] << setw(16) << species_value0[i] /*<< setw(10) << species_value1[i] */<< "
	// dcdt " << dcdt0[i] << "\n";
	//}

	cout << "\n"
	     << "\n";

	// this->KinRate();
	return species_value0; // species_value1
}

/***********************************************************
      Function: REACT_CAP::ODE
                                             DL   10/10
***********************************************************/
vector<double> REACT_CAP::ODE(vector<double> c0, vector<double> h, int iter_max, double t_step)
{
	vector<double> dcdt0, dcdt1;
	double k = 0, t0, t1, t_sum = 0.0;
	int i, j = 0, n;

	while (j < iter_max)
	{
		j++;
		t0 = t_step;
		t1 = t_step;
		n = 0;

		if (j > 1)
			dcdt1 = dcdt0;

		dcdt0 = derivs(c0);

		if (j > 1)
			for (i = 0; i < (int)c0.size(); i++)
				if (dcdt0[i] * dcdt1[i] <= 0)
					n++;

		for (i = 0; i < (int)c0.size(); i++)
		{
			if (h[i] > 0)
			{
				if (dcdt0[i] > 0)
					k = 1 - pow(10, h[i]);
				if (dcdt0[i] < 0)
					k = 1 - pow(10, -h[i]);
				if (dcdt0[i] == 0)
					k = 0.0;
				if (dcdt0[i] != 0 && c0[i] != 0)
					t1 = abs(c0[i] * k / dcdt0[i]);
				cout << " dcdt0 " << setw(16) << dcdt0[i] << " c0 " << setw(16) << c0[i] << " h " << h[i] << " t1 "
				     << setw(16) << t1 /*<< "\n"*/;
				cout << " i " << i << " k " << k << " t0 " << t0 << "\n";
				if (t1 < t0)
					t0 = t1;
			}
		}

		t_sum += t0;
		cout << " t sum " << t_sum << "\n";

		if (t_sum >= t_step)
		{
			for (i = 0; i < (int)c0.size(); i++)
				c0[i] += dcdt0[i] * (t_step - t_sum + t0);
			cout << " Kinetic... "
			     << "\n";
			break;
		}

		if (n == (int)c0.size())
		{
			cout << " Equilibrium... "
			     << "\n";
			break;
		}

		for (i = 0; i < (int)c0.size(); i++)
			c0[i] += dcdt0[i] * t0;
	}

	if (j == iter_max)
		cout << " WARNING... MAX ITERATION"
		     << "\n";
	cout << " iter " << j << "\n";
	return c0;
}

/***********************************************************
      Function: REACT_CAP::derivs

      Task:
      dc/dt=f(c) for Kinetic rate
                                             DL   10/10
***********************************************************/
vector<double> REACT_CAP::derivs(vector<double> c)
{
	vector<double> dcdt;
	double *stoi, wmass;
	int i, ii, i_sp, i_re, i_resp;
	long nscom, noerr = 0;

	// initial output vector
	dcdt.clear();
	for (i = 0; i < (int)c.size(); i++)
		dcdt.push_back(0);

	// for(i=0;i<(int)species_idx.size();i++){
	//	cout << " " << setw(6) << species_phase[i] << setw(6) << species_idx[i] << species_name[i] << "\n";
	//}
	// for(i=0;i<(int)species_kin_idx.size();i++){
	//	cout << " " << setw(6) << species_kin_phase[i] << setw(6) << species_kin_idx[i] << species_kin_name[i] << "\n";
	//}

	this->KinRate();

	// for(i_re=0;i_re<(int)Kin_Reactions.size();i_re++){

	// cout << Kin_Reactions[i_re].type << " " << Kin_Reactions[i_re].rate;
	// for(i_resp=0;i_resp<(int)Kin_Reactions[i_re].species_name.size();i_resp++){
	//	cout << " " << Kin_Reactions[i_re].species_stoi[i_resp] << " " << Kin_Reactions[i_re].species_name[i_resp];
	//}
	// cout << "\n";
	//}

	if (this->mass_type == "SPECIES")
	{
		for (i_sp = 0; i_sp < (int)this->species_name.size(); i_sp++)
		{
			for (i_re = 0; i_re < (int)this->Kin_Reactions.size(); i_re++)
			{
				if (this->Kin_Reactions[i_re].is_Kin)
				{
					for (i_resp = 0; i_resp < (int)this->Kin_Reactions[i_re].species_name.size(); i_resp++)
					{
						if (strcmp(species_name[i_sp], Kin_Reactions[i_re].species_name[i_resp]) == 0
						    && species_phase[i_sp] == Kin_Reactions[i_re].species_phase[i_resp])
						{
							dcdt[i_sp] += Kin_Reactions[i_re].species_stoi[i_resp] * Kin_Reactions[i_re].rate;
						}
					}
				}
			}
		}
		// cout << " SPECIES " << "\n";
	}

	if (this->mass_type == "ELEMENT")
	{
		CAP_tqnosc(&nscom, &noerr);
		stoi = new double[nscom];

		for (i_re = 0; i_re < (int)this->Kin_Reactions.size(); i_re++)
		{
			for (i_resp = 0; i_resp < (int)this->Kin_Reactions[i_re].species_name.size(); i_resp++)
			{
				if (this->Kin_Reactions[i_re].is_Kin)
				{
					// for non-aq phase species, renew value directly
					if (Kin_Reactions[i_re].species_phase[i_resp] != 2)
						for (i_sp = 0; i_sp < (int)this->species_name.size(); i_sp++)
						{
							if (strcmp(species_name[i_sp], Kin_Reactions[i_re].species_name[i_resp]) == 0
							    && species_phase[i_sp] == Kin_Reactions[i_re].species_phase[i_resp])
							{
								dcdt[i_sp] += Kin_Reactions[i_re].species_stoi[i_resp] * Kin_Reactions[i_re].rate;
								// cout << " " << species_name[i_sp] << "\n";
							}
						}

					// for aq species, need to calc new stoichiometry of the species in system
					if (Kin_Reactions[i_re].species_phase[i_resp] == 2)
					{
						CAP_tqstpc(Kin_Reactions[i_re].species_phase[i_resp],
						           Kin_Reactions[i_re].species_idx[i_resp],
						           stoi,
						           &wmass,
						           &noerr);
						ii = 0;
						for (i_sp = 0; i_sp < (int)this->species_name.size(); i_sp++)
						{
							if (species_phase[i_sp] == 2)
							{
								dcdt[i_sp]
								    += stoi[ii] * Kin_Reactions[i_re].species_stoi[i_resp] * Kin_Reactions[i_re].rate;
								// cout << species_name[i_sp] << " " << Kin_Reactions[i_re].type << " " << stoi[ii] << "
								// " << Kin_Reactions[i_re].species_stoi[i_resp] << " " <<
								// Kin_Reactions[i_re].species_name[i_resp] << "\n";
								ii++;
							}
						}
					}
				}
			}
		}
		// cout << " ELEMENT " << "\n";
	}

	// for(i=0;i<(int)c.size();i++)
	//	cout << dcdt[i] << " " << species_name[i] << "\n";

	return dcdt;
}

void REACT_CAP::KinRate(void)
{
	int i, ix, n_r, ii; // j,
	double Kin_rate, T; //, P;

	ii = current_node;
	T = current_TT;
	// P=current_PP;
	// cout << " kinetic minerals " << "\n";
	n_r = (int)this->Kin_Reactions.size();

	for (i = 0; i < n_r; i++)
	{
		if (Kin_Reactions[i].is_Kin)
		{
			Kin_rate = 0;
			// cout << " " << Kin_Reactions[i].type;
			for (ix = 0; ix < (int)Kin_Reactions[i].EA.size(); ix++)
			{
				if (Kin_Reactions[i].ref_species[ix] == "")
					Kin_rate
					    += Kin_Reactions[i].K25[ix] * exp(-Kin_Reactions[i].EA[ix] / 8.314 * (1.0 / T - 1.0 / 298.15));
				else
					Kin_rate += Kin_Reactions[i].K25[ix]
					            * exp(-Kin_Reactions[i].EA[ix] / 8.314 * (1.0 / T - 1.0 / 298.15))
					            * pow(node_ac[ii][Kin_Reactions[i].ref_idx[ix]], Kin_Reactions[i].ref_expo[ix]);
			}
			Kin_rate *= 1.0;
			// cout << " " << Kin_rate << "\n";
			Kin_Reactions[i].rate = Kin_rate;
		}
		else
			Kin_Reactions[i].rate = 0.0;
	}
}

void REACT_CAP::KinInit(void)
{
	int i, ix, jx, n_r;
	size_t j;
	CKinReact* m_kr = NULL;
	string ref_name_kr;
	bool is_species;

	n_r = (int)this->Kin_Reactions.size();
	for (i = 0; i < n_r; i++)
		Kin_Reactions[i].is_Kin = false;

	for (i = 0; i < (int)KinReact_vector.size(); i++)
	{
		m_kr = KinReact_vector[i];
		for (ix = 0; ix < n_r; ix++)
		{
			if (m_kr->chemapp_name.compare(Kin_Reactions[ix].type) == 0)
			{ // found
				Kin_Reactions[ix].is_Kin = true;
				if (m_kr->Km_uniform)
				{
					Kin_Reactions[ix].logK_type = 0;
					Kin_Reactions[ix].logK_kin = m_kr->Km_default;
				}
				else if (m_kr->Km_CHEMAPP)
				{
					Kin_Reactions[ix].logK_type = 1;
					Kin_Reactions[ix].logK_kin = -99;
				}
				else if (m_kr->Km_HKF)
				{
					Kin_Reactions[ix].logK_type = 2;
					Kin_Reactions[ix].logK_kin = -99;
				}
				else
				{
					Kin_Reactions[ix].logK_type = -1;
					Kin_Reactions[ix].logK_kin = -99;
				}
				Kin_Reactions[ix].eta = m_kr->Eta;
				Kin_Reactions[ix].theta = m_kr->Theta;
				Kin_Reactions[ix].SA = m_kr->Am_ini;
				Kin_Reactions[ix].PF = m_kr->precipfactor;

				Kin_Reactions[ix].EA.clear(); // energy activity
				Kin_Reactions[ix].K25.clear(); // logK_25C
				Kin_Reactions[ix].ref_species.clear();
				Kin_Reactions[ix].ref_idx.clear();
				Kin_Reactions[ix].ref_expo.clear();

				for (j = 0; j < m_kr->mechvec.size(); j++)
				{
					Kin_Reactions[ix].EA.push_back(m_kr->mechvec[j]->Eact);
					Kin_Reactions[ix].K25.push_back(m_kr->mechvec[j]->k25);

					if (m_kr->mechvec[j]->mechSpeciesNames.size() > 0)
					{
						ref_name_kr = m_kr->mechvec[j]->mechSpeciesNames[0];
						Kin_Reactions[ix].ref_species.push_back(ref_name_kr);
						Kin_Reactions[ix].ref_expo.push_back(m_kr->mechvec[j]->mechSpeciesExpo[0]);

						is_species = false;
						for (jx = 0; jx < (int)species_name.size(); jx++)
						{
							if (ref_name_kr == species_name[jx] && species_phase[jx] == 2)
							{
								is_species = true;
								Kin_Reactions[ix].ref_idx.push_back(jx);
							}
						}
						if (is_species == false)
						{
							cout << " WARNING ... " << ref_name_kr << " is NOT ChemApp species "
							     << "\n";
							exit(0);
						}
					}
					else
					{
						Kin_Reactions[ix].ref_species.push_back("");
						Kin_Reactions[ix].ref_expo.push_back(-99);
						Kin_Reactions[ix].ref_idx.push_back(-99);
					}
				}
			}
		}
	}

	// for(i=0;i<n_r;i++){
	//	cout << "\n" << Kin_Reactions[i].type << "\n";
	//	if(Kin_Reactions[i].is_Kin){
	//		cout << " Kin true " << "\n";
	//		cout << Kin_Reactions[i].is_Kin << "\n";
	//		cout << Kin_Reactions[i].logK_kin << "\n";
	//		cout << Kin_Reactions[i].logK_type << "\n";
	//		cout << Kin_Reactions[i].theta << "\n";
	//		cout << Kin_Reactions[i].eta << "\n";
	//		cout << Kin_Reactions[i].SA << "\n";
	//		cout << Kin_Reactions[i].PF << "\n";
	//		for(ix=0;ix<(int)Kin_Reactions[i].EA.size();ix++){
	//			cout << Kin_Reactions[i].EA[ix] << "\n";
	//			cout << Kin_Reactions[i].K25[ix] << "\n";
	//			cout << Kin_Reactions[i].ref_species[ix] << "\n";
	//			cout << Kin_Reactions[i].ref_expo[ix] << "\n";
	//			cout << Kin_Reactions[i].ref_idx[ix] << "\n";
	//		}
	//	}
	//	else
	//		cout << " Kin false " << "\n";
	//}
}

bool REACT_CAP::SetKinHKFspecies()
{
	bool sr = true;
	int i = -1, i_re, i_resp, i_sp, is_new = 1;
	SpeciesData HKFspecies;

	this->Kin_HKF_index.clear();
	this->Kin_HKF_species.clear();

	cout << "\n"
	     << " => SetKinHKFspecies() "
	     << "\n";

	for (i_re = 0; i_re < (int)this->Kin_HKF_Reactions.size(); i_re++)
	{
		cout << " Reaction " << Kin_HKF_Reactions[i_re].type << "\n";
		for (i_resp = 0; i_resp < (int)this->Kin_HKF_Reactions[i_re].species_name.size(); i_resp++)
		{
			for (i_sp = 0; i_sp < (int)Kin_HKF_species.size(); i_sp++)
			{
				if (strcmp(this->Kin_HKF_species[i_sp].name0.c_str(),
				           this->Kin_HKF_Reactions[i_re].species_name[i_resp])
				    == 0)
				{
					is_new = 0;
					this->Kin_HKF_index.push_back(i_sp);
				}
			}
			if (is_new == 1)
			{
				i++;
				HKFspecies.name0 = Kin_HKF_Reactions[i_re].species_name[i_resp];
				// cout << " HKFspecies.name " << HKFspecies.name <<  " "
				//<< Kin_HKF_Reactions[i_re].species_name[i_resp] <<"\n";
				// HKFspecies.type=
				this->Kin_HKF_species.push_back(HKFspecies);
				this->Kin_HKF_index.push_back(i);
			}
			is_new = 1;
			cout << " " << setw(4) << right << this->Kin_HKF_Reactions[i_re].species_phase[i_resp];
			cout << setw(4) << this->Kin_HKF_Reactions[i_re].species_stoi[i_resp] << "  ";
			cout << setw(10) << left << this->Kin_HKF_Reactions[i_re].species_name[i_resp] << "\n";
		}
	}
	// for(i=0;i<(int)this->Kin_HKF_index.size();i++)
	//	cout << i << " KIN index " << this->Kin_HKF_index[i] << "\n";
	return sr;
}

bool REACT_CAP::LoadHKFparam()
{
	cout << "\n"
	     << " => LoadHKFparam "
	     << "\n";
	bool sr = true;
	int i, ix, iy, is_load, type = -1;
	double charge = -99, HKFparam[9][4];
	std::string species_name0;

	// for(ix=0;ix<9;ix++)
	//	for(iy=0;iy<4;iy++)
	//		HKFparam[ix][iy]=0.0;
	// cout << " size " << this->Kin_HKF_species.size() << "\n";
	for (i = 0; i < (int)this->Kin_HKF_species.size(); i++)
	{
		// cout << i << " " << Kin_HKF_species[i].name0 << " " << "\n";
		if (Kin_HKF_species[i].name0.size() > 14)
		{
			cout << " Warning ... " << Kin_HKF_species[i].name0 << " is not Kin_HKF_species in slop98.dat. EXIT"
			     << "\n";
			exit(0);
		}
		if (Kin_HKF_species[i].name0 == "H2O" || Kin_HKF_species[i].name0 == "H2O_liquid"
		    || Kin_HKF_species[i].name0 == "water_liquid")
			Kin_HKF_species[i].type = -1;
		else
		{
			species_name0 = Kin_HKF_species[i].name0;
			is_load = HKF::HKF_OGS_loadparam(species_name0, type, charge, HKFparam);

			if (is_load == 1)
			{
				Kin_HKF_species[i].type = type;
				Kin_HKF_species[i].charge = charge;
				for (ix = 0; ix < 9; ix++)
					for (iy = 0; iy < 4; iy++)
						Kin_HKF_species[i].param[ix][iy] = HKFparam[ix][iy];
			}
			else
			{
				Kin_HKF_species[i].type = -1;
				cout << " Warning ... load " << Kin_HKF_species[i].name0 << " parameters failed ... "
				     << "\n";
			}
		}
	}

	// is_load=HKF_load_param(spec,name_type,phase_type);

	// if(HKF_load_param(spec,name_type,phase_type)==1){
	//	cout << " HKF load param...... OK! " << "\n";
	//	i_counter=0;
	//	for(i=0;i<Kin_HKF_species.size();i++)
	//		if(Kin_HKF_species[i].name!="H2O" && Kin_HKF_species[i].name!="H2O_liquid"){
	//			Kin_HKF_species[i]=spec[i_counter];
	//			i_counter++;
	//		}
	//}

	return sr;
}

bool REACT_CAP::HKFcalc(double T, double P)
{
	bool sr = true;
	int i; // ix,iy
	// int type;
	double G, H, S; // ,param[9][4],charge

	for (i = 0; i < (int)this->Kin_HKF_species.size(); i++)
	{
		if (Kin_HKF_species[i].type != -1)
		{
			HKF::HKF_OGS_calc(
			    T, P, G, H, S, Kin_HKF_species[i].type, Kin_HKF_species[i].charge, Kin_HKF_species[i].param);
			// cout << Kin_HKF_species[i].name0 << " " << G << " " << H << " " << S << "\n";
			Kin_HKF_species[i].G = G;
			Kin_HKF_species[i].H = H;
			Kin_HKF_species[i].S = S;
		}
		else if (Kin_HKF_species[i].name0 == "H2O" || Kin_HKF_species[i].name0 == "H2O_liquid"
		         || Kin_HKF_species[i].name0 == "water_liquid")
		{
			Kin_HKF_species[i].G = HKF::HKFcalcw(T, P, 0);
			Kin_HKF_species[i].H = HKF::HKFcalcw(T, P, 1);
			Kin_HKF_species[i].S = HKF::HKFcalcw(T, P, 2);
		}
	}
	return sr;
}

/**************************************************************************/
/*  Global function  */
/**************************************************************************/
bool REACT_CAP_Read(std::string file_base_name, const GEOLIB::GEOObjects& geo_obj, const std::string& unique_name)
{
	char line[MAX_ZEILE];
	std::string file_name_cap, line_string;
	ios::pos_type position;

	REACT_CAP* rc_cap = new REACT_CAP();
	// look if file is there
	file_name_cap = file_base_name + REACTION_EXTENSION_CHEMAPP;
	ifstream cap_file(file_name_cap.data(), ios::in);
	if (!cap_file.good())
	{
		delete rc_cap;
		rc_cap = NULL;
	}
	else
	{ // file is there - use ChemApp
		if (REACT_vec.capacity() > 0)
		{ // Test, if PHREEQC is used also
			cout << "\n"
			     << " Warning!  ChemApp and PHREEQC both actived. ChemApp will NOT be used ! "
			     << "\n";
			delete rc_cap;
			rc_cap = NULL;
		}
		else
		{
			rc_cap->flag_cap = true;
			// Read input file *.cap
			cap_file.clear();
			cap_file.seekg(0, ios::beg);
			cout << "ChemApp_Read"
			     << "\n";
			cap_file.getline(line, MAX_ZEILE); // first line
			cap_file.getline(line, MAX_ZEILE); // second line ToDo
			line_string = line;
			if (line_string.find("#STOP") != string::npos)
				return true;

			// Call the object read function
			position = rc_cap->Read(&cap_file, geo_obj, unique_name);
			// store data in global vector
			REACT_CAP_vec.push_back(rc_cap);
			// close file
			cap_file.close();
		}
	}
	return true;
}
/***********************************************************
      Function: REACT_CAP::ConvertIC2BC

      Task:
      Convert calculated equilibrium values
      at all mass transport processes to be used as
      boundary condition for the remainsing simulation
      at a specified geometry

      SB   10/2009   First Implementation
***********************************************************/
void REACT_CAP::ConvertIC2BC(const GEOLIB::GEOObjects& geo_obj, const std::string& unique_name)
{
	double species_value;
	int ns = 0, i, ii, j, l, ll;
	CRFProcess* m_pcs = NULL;
	vector<long> node_nr;
	ns = (int)species_name.size();
	std::string s_geo_name, s_geo_type, testname;
	vector<long> nodes_vector;
	bool replace = false;

	// Do only if flag is specified
	if (this->ic_2_bc_geometry_name.size() > 0)
	{
		cout << "\n"
		     << "\n"
		     << " => ConvertIC2BC() "
		     << "\n";

		// Get nodes at specified geometry

		/********************************************************/
		CFEMesh* m_msh = fem_msh_vector[0]; // SB: ToDo hart gesetzt
		if (m_msh == NULL)
		{
			cout << "ConvertIC2BC"
			     << "\n";
			exit(0);
		}
		// Go through specified geometry elements
		for (j = 0; j < (int)this->ic_2_bc_geometry_name.size(); j++)
		{
			s_geo_name = this->ic_2_bc_geometry_name[j];
			s_geo_type = this->ic_2_bc_geometry_type[j];
			//------------------------------------------------------------------
			if (s_geo_type.compare("POINT") == 0)
			{
				// 06/2010 TF switch to new GEOLIB
				// CGLPoint* m_geo_point = NULL; // make new GEO point
				// m_geo_point = GEOGetPointByName(s_geo_name);//Get GEO point by name
				// if(m_geo_point)
				//  l = m_msh->GetNODOnPNT(m_geo_point); // + ShiftInNodeVector; // find MSH point number stored in l
				const std::vector<GEOLIB::Point*>* pnt_vec(geo_obj.getPointVec(unique_name));
				l = m_msh->GetNODOnPNT((*pnt_vec)[this->ic_2_bc_GeoID[j]]);
				node_nr.push_back(l);
			}
			//------------------------------------------------------------------
			if (s_geo_type.compare("POLYLINE") == 0)
			{
				// CGLPolyline *m_polyline = NULL;
				// m_polyline = GEOGetPLYByName(s_geo_name);// get Polyline by name
				// CGLPolyline *m_polyline (polyline_vector[this->ic_2_bc_GeoID[j]]);
				std::vector<GEOLIB::Polyline*> const* const ply_vec(geo_obj.getPolylineVec(unique_name));
				GEOLIB::Polyline const* const m_polyline((*ply_vec)[this->ic_2_bc_GeoID[j]]);

				if (m_polyline)
				{
					// if(m_polyline->getType()==100) //WW
					// m_msh->GetNodesOnArc(m_polyline,nodes_vector);
					// else
					m_msh->GetNODOnPLY(m_polyline, nodes_vector);
					for (i = 0; i < (long)nodes_vector.size(); i++)
					{
						ll = nodes_vector[i];
						l = ll; //+ShiftInNodeVector;
						node_nr.push_back(l);
					}
				}
			} // if(POLYLINE)
			//------------------------------------------------------------------
			if (s_geo_type.compare("SURFACE") == 0)
			{
				Surface* m_surface = NULL;
				m_surface = GEOGetSFCByName(s_geo_name);
				if (m_surface)
				{
					m_msh->GetNODOnSFC(m_surface, nodes_vector);
					for (i = 0; i < (long)nodes_vector.size(); i++)
					{
						ll = nodes_vector[i];
						l = ll; //+ShiftInNodeVector;
						node_nr.push_back(l);
					}
				}
			}
		} // end of for(j=0;j<m_krd->NoReactGeoName.size()...
		// test output
		/* */ cout << " Vector node_nr: "
		           << "\n";
		for (l = 0; l < (long)node_nr.size(); l++)
			cout << " Node number: " << node_nr[l] << ", ";
		cout << "\n";

		// Transfer simulated value to bc data structure
		// First get them from results vectors of mass transport processes and store new values in species_value
		for (ii = 0; ii < (int)node_nr.size(); ii++)
		{ // count through nodes on geometry
			for (i = 0; i < ns; i++)
			{
				m_pcs = pcs_vector[this->pcs_mass_idx[i]];
				// species_value = m_pcs->GetNodeValue(ii,0);  // CB DL todo check, which one is the correct
				species_value = m_pcs->GetNodeValue(ii, 1);
				// check if bc replacement is wanted for this species = process
				replace = false;
				for (int jj = 0; jj < (int)this->ic_2_bc_species.size(); jj++)
				{
					testname = this->ic_2_bc_species[jj];
					if (testname.compare(m_pcs->pcs_primary_function_name[0]) == 0)
					{
						replace = true;
						species_value += this->ic_2_bc_species_value[jj];
						break;
					}
				}

				if (ic_2_bc_species.size() == 0)
					replace = true;

				if (replace)
				{
					cout << " Species : " << setw(4) << this->pcs_mass_idx[i] << " " << setw(16)
					     << m_pcs->pcs_primary_function_name[0] << " ";
					// cout << old/new c in m_pcs: " << m_pcs->GetNodeValue(ii,0) <<" / " << m_pcs->GetNodeValue(ii,1)
					// << "\n";

					// Store new value in boundary condition
					for (j = 0; j < (int)m_pcs->bc_node_value.size(); j++)
					{
						if (m_pcs->bc_node_value[j]->msh_node_number == node_nr[ii])
						{
							cout << "      Replaced bc value  " << m_pcs->bc_node_value[j]->node_value << "  by  "
							     << species_value << "\n";
							pcs_vector[this->pcs_mass_idx[i]]->bc_node_value[j]->node_value = species_value;
						}
					}
				} // end if(replace)
			}
		} // end for(ii
	}
}
