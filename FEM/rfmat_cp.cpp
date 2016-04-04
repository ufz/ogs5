/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************/
/* ROCKFLOW - Modul: rfmat_cp.c
 */
/* Task:
   Methods for ComponentProperties
 */
/* Programming:
   10/2004   SB  First Implemented
 */
/**************************************************************************/
#ifdef WIN32
#pragma warning(disable : 4786) /*Visual C++ 6.0*/
#endif
// C
#include <math.h>
// C++
#include <cfloat>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

#include "files0.h"
#include "makros.h"
#include "display.h"
#include "rf_mfp_new.h"
#include "rf_mmp_new.h"
#include "rf_msp_new.h"
#include "rf_pcs.h"
#include "rf_pcs.h" //OK_MOD" //GetRFProcess..()
#include "rf_tim_new.h"
#include "rfmat_cp.h"
#include "tools.h"
#include "rf_react.h"
#include "rf_react_int.h"

#ifdef GEM_REACT
#include "rf_REACT_GEM.h"
#endif
using SolidProp::CSolidProperties;
/* Vector auf CompProperties , globale Zugriffe */
// vector <CompProperties*> cp_vec;
// do not need this anymore, use global map structure instead.
std::map<int, CompProperties*> cp_vec;
std::map<std::string, int> cp_name_2_idx;

/*========================================================================*/
/* Component Properties                                                  */
/*========================================================================*/

/**************************************************************************
   FEMLib-Method:
   Task: CompProperties Constructor
   Programing:
   02/2004 SB Implementation
**************************************************************************/
CompProperties::CompProperties(/* int n // HS we do not need this. */) : idx(std::numeric_limits<size_t>::max())
{
	// if ( idx != std::numeric_limits<size_t>::max() ) // this means idx is set.

	// compname = name1; // HS it will be loaded later.
	iupac_formula = ""; // CB 0311 check
	mobil = 1; // by default, set to mobile species.
	transport_phase = 0; // by default, set to the 1st phase.
	fluid_phase = 0; // by default, set to water

	diffusion_model = -1;
	count_of_diffusion_model_values = 0;

	decay_model = -1;
	count_of_decay_model_values = 0;

	isotherm_model = -1;
	count_of_isotherm_model_values = 0;

	bubble_velocity_model = -1;
	bubble_velocity[0] = bubble_velocity[1] = bubble_velocity[2] = 0.0;
	file_base_name = "nix";
	valence = 0;
	a_zero = 0;
	mineral_density = 0;
	molar_density = 0;
	molar_weight = 0;
	max_solubility = 0;

	//// CB _ctx_
	//   _ctx_ = false;
	//   ct_substratename = "";
	//   _ctx_Coefficient = 0;

	this->setProcessType(FiniteElement::MASS_TRANSPORT);
	this->setProcessPrimaryVariable(FiniteElement::CONCENTRATION);

	OutputMassOfComponentInModel = 0; // 05/2012 BG
}

/**************************************************************************
   FEMLib-Method:
   Task: CompProperties Destructor
   Programing:
   02/2004 SB Implementation
**************************************************************************/
CompProperties::~CompProperties(void)
{
}

/**************************************************************************
   FEMLib-Method:
   Task: CompProperties read function
   Programing:
   02/2004 SB Implementation - adapted from OK rf_bc_new
   10/2004 SB Adapted to new file structure
   01/2005 OK boolean type
   01/2005 OK Destruct before read
**************************************************************************/
bool CPRead(std::string file_base_name)
{
	//----------------------------------------------------------------------
	// OK  MCPDelete();
	//----------------------------------------------------------------------
	CompProperties* m_cp = NULL;
	char line[MAX_ZEILE];
	std::string sub_line;
	std::string line_string;
	std::ios::pos_type position;
	//========================================================================
	// File handling
	std::string cp_file_name = file_base_name + CP_FILE_EXTENSION;
	ifstream cp_file(cp_file_name.data(), ios::in);
	if (!cp_file.good())
	{
		if (pcs_vector[0]->getProcessType() == FiniteElement::TNEQ)
		{
			DisplayMsgLn("ERROR. TNEQ requires specification of inert and reactive components in mcp file.");
			exit(1);
		}
		if (pcs_vector[0]->getProcessType() == FiniteElement::TES)
		{
			DisplayMsgLn("ERROR. TES requires specification of inert and reactive components in mcp file.");
			exit(1);
		}
		return false;
	}
	cp_file.seekg(0L, ios::beg);

	//========================================================================
	cp_vec.clear();
	cout << "CPRead"
	     << "\n";
	// Schleife ueber alle Phasen bzw. Komponenten
	while (!cp_file.eof())
	{
		cp_file.getline(line, MAX_ZEILE);
		line_string = line;
		if (line_string.find("#STOP") != std::string::npos)
			break;
		//----------------------------------------------------------------------
		// keyword found
		if (line_string.find("#COMPONENT_PROPERTIES") != std::string::npos)
		{
			m_cp = new CompProperties();
			m_cp->file_base_name = file_base_name;
			position = m_cp->Read(&cp_file);
			// HS the index of this component is filled
			// one after another by its sequence in mcp file.
			m_cp->idx = cp_vec.size();
			cp_name_2_idx[m_cp->compname] = m_cp->idx;
			cp_vec[m_cp->idx] = m_cp;
			m_cp = NULL;
			cp_file.seekg(position, ios::beg);
		} // keyword found
	} // eof
	// immediately check if enough PCS objects are available
	size_t pcs_mt_count = 0;
	size_t pcs_rwpt_count = 0;
	size_t i;
	for (i = 0; i < pcs_vector.size(); i++)
	{
		if (pcs_vector[i]->getProcessType() == FiniteElement::MASS_TRANSPORT)
			pcs_mt_count++;
		if (pcs_vector[i]->getProcessType() == FiniteElement::RANDOM_WALK)
			pcs_rwpt_count++;
	}
	if (pcs_rwpt_count == 0) // HS, no random walk detected.
	{
		if ((pcs_mt_count != cp_vec.size() || pcs_mt_count != cp_name_2_idx.size())
		    && pcs_vector[0]->getProcessType() != FiniteElement::TNEQ
		    && pcs_vector[0]->getProcessType() != FiniteElement::TES)
		{
			DisplayMsgLn("Mass transport components and Mass transport processes do not fit!");
			exit(1);
		}
		else if (cp_vec.size() < 2 && (pcs_vector[0]->getProcessType() == FiniteElement::TNEQ
		                               || pcs_vector[0]->getProcessType() == FiniteElement::TES))
		{
			DisplayMsgLn("ERROR. TNEQ/TES requires specification of inert and reactive components in mcp file.");
			exit(1);
		}
		else
		{
			// and then link MCP with the PCS.
			std::map<int, CompProperties*>::iterator cp_iter = cp_vec.begin();
			for (i = 0; i < pcs_vector.size(); i++)
			{
				if (pcs_vector[i]->getProcessType() == FiniteElement::MASS_TRANSPORT)
				{
					cp_iter->second->setProcess(pcs_vector[i]);
					++cp_iter;
				}
				if (pcs_vector[i]->getProcessType() == FiniteElement::TNEQ
				    || pcs_vector[i]->getProcessType() == FiniteElement::TES)
					std::cout
					    << "Warning! For TNEQ/TES, the component order in the mcp file needs to be INERT, REACTIVE!\n";
			}
		} // end of else

		// Assign fluid id's for use in fluid property calculation
		if (pcs_vector[0]->getProcessType() == FiniteElement::TNEQ
		    || pcs_vector[0]->getProcessType() == FiniteElement::TES)
			for (i = 0; i < cp_vec.size(); i++)
			{
				if (cp_vec[i]->compname == "N2")
					cp_vec[i]->fluid_id = 3;
				else if (cp_vec[i]->compname == "H2O")
					cp_vec[i]->fluid_id = 1;
				else if (cp_vec[i]->compname == "O2")
				{
					cp_vec[i]->fluid_id = 5;
					std::cout << "Warning: The thermal conductivity critical enhancement model for oxygen has not been "
					             "implemented. Dilute and residual parts only.\n";
				}
				else
					std::cout << "Warning: Component name unknown to TNEQ/TES.\n";
			}
	}
	return true;
}

/**************************************************************************
   FEMLib-Method:
   Task: CompProperties read function
   Programing:
   02/2004 SB Implementation - adapted from OK rf_bc_new
   05/2007 PCH: Anisotropic diffusion coefficient added
**************************************************************************/
ios::pos_type CompProperties::Read(ifstream* rfd_file)
{
	//  char line[MAX_ZEILE];
	std::string sub_line;
	std::string line_string;
	std::string delimiter(" ");
	bool new_keyword = false;
	std::string hash("#");
	std::stringstream in;
	int j;
	//  double *read_help;
	//  long index;
	ios::pos_type position;

	//========================================================================
	// Schleife ueber alle Phasen bzw. Komponenten
	while (!new_keyword)
	{
		//    new_subkeyword = false;
		position = rfd_file->tellg();

		//    if(!GetLineFromFile(line,rfd_file)) break;
		//    line_string = line;
		line_string = GetLineFromFile1(rfd_file);
		if (line_string.size() < 1)
			break;

		if (line_string.find(hash) != std::string::npos)
		{
			new_keyword = true;
			break;
		}

		/* Keywords nacheinander durchsuchen */
		//....................................................................
		// subkeyword found
		if (line_string.find("$NAME") != std::string::npos)
		{
			in.str(GetLineFromFile1(rfd_file));
			in >> compname; // sub_line
			in.clear();
			//	  compname = (char *) sub_line.c_str();
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$FORMULA") != std::string::npos)
		{ // subkeyword found
			in.str(GetLineFromFile1(rfd_file));
			in >> iupac_formula; // sub_line
			in.clear();
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$MOBILE") != std::string::npos)
		{
			//      rfd_file->getline(line,MAX_ZEILE);
			//      line_string = line;
			//      mobil = atoi(line_string.substr(0).data());
			in.str(GetLineFromFile1(rfd_file));
			in >> mobil;
			in.clear();
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$TRANSPORT_PHASE") != std::string::npos)
		{
			//      rfd_file->getline(line,MAX_ZEILE);
			//      line_string = line;
			//	  transport_phase = atoi(line_string.substr(0).data());
			in.str(GetLineFromFile1(rfd_file));
			in >> transport_phase;
			in.clear();
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$FLUID_PHASE") != std::string::npos)
		{
			in.str(GetLineFromFile1(rfd_file));
			in >> fluid_phase;
			in.clear();
			// critical constant and corresponding components fluid Property // AKS
			if (mfp_vector[fluid_phase])
				// found corresponding fluid phase.
				mfp_vector[fluid_phase]->component_vector.push_back(this);
			else
			{
				cout << "Error! Corresponding fluid phase not found!!! Terminate program..."
				     << "\n";
				exit(0);
			}
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$MOL_MASS") != std::string::npos)
		{
			in.str(GetLineFromFile1(rfd_file));
			in >> molar_mass;
			in.clear();
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$CRITICAL_PRESSURE") != std::string::npos)
		{
			in.str(GetLineFromFile1(rfd_file));
			in >> pc;
			in.clear();
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$CRITICAL_TEMPERATURE") != std::string::npos)
		{
			in.str(GetLineFromFile1(rfd_file));
			in >> Tc;
			in.clear();
		}
		//....................................................................
		if (line_string.find("$ACENTRIC_FACTOR") != std::string::npos)
		{
			in.str(GetLineFromFile1(rfd_file));
			in >> omega;
			in.clear();
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$FLUID_ID") != std::string::npos)
		{
			in.str(GetLineFromFile1(rfd_file));
			in >> fluid_id;
			in.clear();
		}
		if (line_string.find("$MOLAR_VOLUME") != std::string::npos)
		{
			in.str(GetLineFromFile1(rfd_file));
			in >> Vm;
			in.clear();
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$VOLUME_DIFFUSION") != std::string::npos)
		{
			in.str(GetLineFromFile1(rfd_file));
			in >> Vd;
			in.clear();
		}
		//....................................................................
		if (line_string.find("$MINERAL_DENSITY") != std::string::npos)
		{ // subkeyword found [g/cm3]:AB
			in.str(GetLineFromFile1(rfd_file));
			in >> mineral_density;
			in.clear();
		}

		//....................................................................
		// subkeyword found
		if (line_string.find("$DIFFUSION") != std::string::npos)
		{
			in.str(GetLineFromFile1(rfd_file));
			in >> diffusion_model;
			if (diffusion_model == 0) // read curve

				in >> diffusion_function_name;
			else
			{
				count_of_diffusion_model_values = GetNumberDiffusionValuesCompProperties(diffusion_model);
				// unknown parameter model
				if (count_of_diffusion_model_values < 0)
				{
					DisplayMsgLn(" Unknown Diffusion model - program stops !");
					exit(1);
				}

				if (diffusion_model > 0)
					//		read_help = (double *) Malloc(count_of_diffusion_model_values * sizeof(double));

					for (j = 0; j < count_of_diffusion_model_values; j++)
					{
						if (in.peek() > 0)

							/*
							            if(StrTestInv(&sub[p_sub += pos],&pos)) {
							                    StrReadString (&name,&sub[p_sub += pos],f,TFString,&pos);
							                    sprintf(group_name,"%s%d",name_component_properties,j);
							                    ptr_d=get_tp_diffusion_model_value_ptr(cp1,j);
							                    SetInverseVariable(name,group_name,1,(void *)ptr_d);
							                    name=(char *)Free(name);
							                }
							                else  */
							//					in >> read_help[j];
							in >> diffusion_model_values[j];
						else
							cout << "Warning: Missing diffusion model values for component " << this->compname << "\n";
						if ((diffusion_model == 1) && (j == 0))
							if ((diffusion_model_values[j] < 0.0) || (diffusion_model_values[j] > 100000.0))
								cout << "Warning: Funny diffusion model values specified"
								     << "\n";
					} // end for(j...)

				//		diffusion_model_values = read_help;

				if (diffusion_model < 0)
					DisplayMsgLn("Error: Diffusion model must be larger than or 0");
			}
			in.clear();
		} // subkeyword found
		// subkeyword found
		if (line_string.find("$DECAY") != std::string::npos)
		{
			in.str(GetLineFromFile1(rfd_file));
			in >> decay_model;
			if (decay_model == 0) // curve

				in >> decay_function_name;
			else
			{
				count_of_decay_model_values = GetNumberDecayValuesCompProperties(decay_model);

				if (count_of_decay_model_values < 0) // unknown parameter model
				{
					DisplayMsgLn(" Unknown Aqueous Decay model - program stops !");
					exit(1);
				}

				if (decay_model > 0)
					//		read_help = (double *) Malloc(count_of_decay_model_values * sizeof(double));

					for (j = 0; j < count_of_decay_model_values; j++)
						/*
						            if(StrTestInv(&sub[p_sub += pos],&pos)) {
						                    StrReadString (&name,&sub[p_sub += pos],f,TFString,&pos);
						                    sprintf(group_name,"%s%d",name_component_properties,j);
						                    ptr_d=get_tp_diffusion_model_value_ptr(cp1,j);
						                    SetInverseVariable(name,group_name,1,(void *)ptr_d);
						                    name=(char *)Free(name);
						                }
						                else  */
						//					in >> read_help[j];
						in >> decay_model_values[j];

				// end for(j...)

				//		decay_model_values = read_help;

				if (decay_model < -1)
					DisplayMsgLn("Error: Aqueous decay model must be larger than or 0");
			}
			in.clear();
		} // subkeyword found
		// subkeyword found
		if (line_string.find("$ISOTHERM") != std::string::npos)
		{
			in.str(GetLineFromFile1(rfd_file));
			in >> isotherm_model;
			if (isotherm_model == 0)
				in >> isotherm_function_name;
			else
			{
				count_of_isotherm_model_values = GetNumberIsothermValuesCompProperties(isotherm_model);

				if (count_of_isotherm_model_values < 0) // unknown parameter model
				{
					DisplayMsgLn(" Unknown Isotherm model - program stops !");
					exit(1);
				}

				if (isotherm_model > 0)
					//		read_help = (double *) Malloc(count_of_isotherm_model_values * sizeof(double));

					for (j = 0; j < count_of_isotherm_model_values; j++)
						/*
						            if(StrTestInv(&sub[p_sub += pos],&pos)) {
						                    StrReadString (&name,&sub[p_sub += pos],f,TFString,&pos);
						                    sprintf(group_name,"%s%d",name_component_properties,j);
						                    ptr_d=get_tp_diffusion_model_value_ptr(cp1,j);
						                    SetInverseVariable(name,group_name,1,(void *)ptr_d);
						                    name=(char *)Free(name);
						                }
						                else  */
						//					in >> read_help[j];
						in >> isotherm_model_values[j];

				// end for(j...)

				//		isotherm_model_values = read_help;
				if (isotherm_model < 0)
					DisplayMsgLn("Error: Isotherm model must be larger than or 0");
			}
			in.clear();
		} // subkeyword found

		// subkeyword found
		if (line_string.find("$BUBBLE_VELOCITY") != std::string::npos)
		{
			in.str(GetLineFromFile1(rfd_file));
			in >> bubble_velocity_model;
			if (bubble_velocity_model == 1)
				in >> bubble_velocity[0] >> bubble_velocity[1] >> bubble_velocity[2];
			in.clear();
		} // subkeyword found

		//....................................................................

		// parameters for NAPL dissolution CB140708
		// subkeyword found
		if (line_string.find("$MOLAR_DENSITY") != std::string::npos)
		{
			in.str(GetLineFromFile1(rfd_file));
			in >> molar_density;
			in.clear();
			if (molar_density <= 0) // unphysical entry
			{
				DisplayMsgLn("Error in MOLAR_DENSITY - setting molar_density to 1.0!");
				molar_density = 1.0;
			}
		}
		// subkeyword found
		if (line_string.find("$MOLAR_WEIGHT") != std::string::npos)
		{
			in.str(GetLineFromFile1(rfd_file));
			in >> molar_weight;
			in.clear();
			if (molar_weight <= 0) // unphysical entry
			{
				DisplayMsgLn("Error in MOLAR_WEIGHT - setting molar_weight to 1.0!");
				molar_weight = 1.0;
			}
		}
		// subkeyword found
		if (line_string.find("$MAXIMUM_AQUEOUS_SOLUBILITY") != std::string::npos)
		{
			in.str(GetLineFromFile1(rfd_file));
			in >> max_solubility;
			in.clear();
			if (max_solubility <= 0) // unphysical entry
			{
				DisplayMsgLn(": Error in MAXIMUM_AQUEOUS_SOLUBILITY - setting max_solubility to 1.0!");
				max_solubility = 1.0;
			}
		}
		//....................................................................
		// subkeyword found
		// Output of the Mass of the component within the model if subkeyword is found, BG 05/2012
		if (line_string.find("$OutputMassOfComponentInModel") != std::string::npos)
		{
			OutputMassOfComponentInModel = 1;
		}
		//....................................................................
		// parameters for MINERAL KINETICS CB0410
		if (line_string.find("$VALENCE") != std::string::npos)
		{ // subkeyword found
			in.str(GetLineFromFile1(rfd_file));
			in >> valence;
			in.clear();
			if (fabs(double(valence)) >= 7)
			{ // unphysical entry
				DisplayMsgLn("Error in VALENCE - setting valence to 0!");
				valence = 0;
			}
		}
		if (line_string.find("$A_ZERO") != std::string::npos)
		{ // subkeyword found
			in.str(GetLineFromFile1(rfd_file));
			in >> a_zero;
			in.clear();
			if (a_zero >= 10)
			{ // unphysical entry
				// DisplayMsgLn("Error in A_ZERO - setting valence to 0!");
				// a_zero = 0.0;
			}
		}
	} // end while
	return position;
}

/**************************************************************************
   CompProperties - Method: Write
   Task: ComponentProperties write function - echo of input values to rfe - file
   Programing:
   05/2004 SB Implementation
**************************************************************************/
void CPWrite(std::string base_file_name, int flag)
{
	CompProperties* m_kr = NULL;
	std::string rfe_file_name;
	int i, length;
	std::ofstream rfe_file;

	//========================================================================
	// File handling
	if (flag == 1)
	{
		//	rfe_file_name = base_file_name + ".rfe";
		rfe_file_name = base_file_name;
		rfe_file.open(rfe_file_name.data(), ios::app);
	}
	else
	{
		rfe_file_name = base_file_name + CP_FILE_EXTENSION;
		rfe_file.open(rfe_file_name.data(), ios::trunc | ios::out);
	}

	if (!rfe_file.good())
		return;
	if (flag == 1)
		rfe_file.seekp(0L, ios::end); // go to end
	//========================================================================
	rfe_file << "\n"
	         << "; CompProperties ----------------------------------------------------------------- "
	         << "\n";
	// Output all components
	length = (int)cp_vec.size();
	for (i = 0; i < length; i++)
	{
		m_kr = cp_vec[i];
		m_kr->Write(&rfe_file);
	}
	rfe_file << "\n"
	         << "#STOP "
	         << "\n";
	rfe_file.close();
	//  delete rfe_file;
}

/**************************************************************************
   Component Properties:
   Task: Component class write function
   Programing:
   05/2004 SB Implementation - adapted from OK rf_bc_new
**************************************************************************/
void CompProperties::Write(ofstream* rfe_file)
{
	int i;

	// Write Keyword
	*rfe_file << "#COMPONENT_PROPERTIES"
	          << "\n";
	// Name of component
	*rfe_file << "$NAME"
	          << "\n"
	          << compname << "\n";
	// chemical formula
	*rfe_file << "$FORMULA"
	          << "\n"
	          << iupac_formula << "\n";
	// mobile or not?
	*rfe_file << "$MOBILE"
	          << "\n"
	          << mobil << "\n";
	// TRANSPORT_PHASE
	*rfe_file << "$TRANSPORT_PHASE"
	          << "\n"
	          << transport_phase << "\n";
	// FLUID_PHASE
	*rfe_file << "$FLUID_PHASE"
	          << "\n"
	          << transport_phase << "\n";
	// Diffusion
	if (diffusion_model > -1)
	{
		*rfe_file << "$DIFFUSION"
		          << "\n";
		*rfe_file << diffusion_model << "  ";
		for (i = 0; i < count_of_diffusion_model_values; i++)
			*rfe_file << diffusion_model_values[i] << " ";
		*rfe_file << "\n";
	}
	// DECAY
	if (decay_model > -1)
	{
		*rfe_file << "$DECAY"
		          << "\n";
		*rfe_file << decay_model << "  ";
		for (i = 0; i < count_of_decay_model_values; i++)
			*rfe_file << decay_model_values[i] << " ";
		*rfe_file << "\n";
	}
	// Isotherm
	if (isotherm_model > -1)
	{
		*rfe_file << "$ISOTHERM"
		          << "\n";
		*rfe_file << isotherm_model << "  ";
		for (i = 0; i < count_of_isotherm_model_values; i++)
			*rfe_file << isotherm_model_values[i] << " ";
		*rfe_file << "\n";
	}
	// NAPL Dissolution CB140708
	if (molar_density > 0)
	{
		*rfe_file << "$MOLAR_DENSITY"
		          << "\n";
		*rfe_file << molar_density << "  ";
		*rfe_file << "\n";
	}
	if (molar_weight > 0)
	{
		*rfe_file << "$MOLAR_WEIGHT"
		          << "\n";
		*rfe_file << molar_weight << "  ";
		*rfe_file << "\n";
	}
	if (max_solubility > 0)
	{
		*rfe_file << "$MAXIMUM_AQUEOUS_SOLUBILITY"
		          << "\n";
		*rfe_file << max_solubility << "  ";
		*rfe_file << "\n";
	}
	if (a_zero > 0)
	{
		*rfe_file << "$A_ZERO"
		          << "\n";
		*rfe_file << a_zero << "  ";
		*rfe_file << "\n";
	}
	if (valence)
	{
		*rfe_file << "$VALENCE"
		          << "\n";
		*rfe_file << valence << "  ";
		*rfe_file << "\n";
	}
	if (mineral_density > 0)
	{
		*rfe_file << "$MINERAL_DENSITY"
		          << "\n";
		*rfe_file << mineral_density << "  ";
		*rfe_file << "\n";
	}

	*rfe_file << "\n";
}

/*************************************************************************
   ROCKFLOW - Funktion: CalcDiffusionCoefficientCP

   Aufgabe:
   Liefert die Diffusion in Abhaengigkeit von den uebergebenen
   Argumenten

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)

   Ergebnis:
   -1.0: Fehler
   Sonst Wert der berechneten Diffusion

   Programmaenderungen:
   06/2000    AH      Erste Version
   09/2003    MX      Cases 1, 2 von Methode 1 in die Funktion gekapselt
   02/2004    SB      Angepasst an CompProperties
   09/2004	SB		Fitted to CP Class Structure, now CompProperties Member function

*************************************************************************/
double CompProperties::CalcDiffusionCoefficientCP(long index, double theta, CRFProcess* m_pcs)
{
	(void)theta;
	int p_idx = -1, t_idx = -1;
	/*int dependence = 0; */
	double diffusion_coefficient = -1.0;
	double pressure_average = 1e5;
	double temperature_average = 293.;
	double diffusion_average = 0.0;
	double* k = NULL;
	double Dm;
	// static long *element_nodes;
	CRFProcess* t_pcs = NULL;
#ifdef GEM_REACT
	static int count_nodes;
	// TF unused variable - fix a compiler warning
	//   static double eta = 0.0;                       //, theta = 1.0;
	double porosity = 1.0, saturation = 1.0, dummy = 0.0;
	int i, flowflag;
	CRFProcess* f_pcs = NULL;
#endif
	// static int p_ind, t_ind;
	// OK411
	diffusion_average = diffusion_average;
	pressure_average = pressure_average;
	temperature_average = temperature_average;
	p_idx = p_idx;
	t_idx = t_idx;

#ifdef GEM_REACT
	const long group(m_pcs->m_msh->ele_vector[index]->GetPatchIndex());
	CMediumProperties* m_mat_mp(mmp_vector[group]);
#endif

	k = diffusion_model_values;

	switch (diffusion_model)
	{
		case -1:
			return 0.0; // no diffusion specified
		case 0: /* curve value */
		{
			// DisplayMsgLn("Not implemented");
			return 0.0;
		}
		case 1: /* Konstanter Diffusionswert */
		{
			if (count_of_diffusion_model_values < 1)
				return 0.0;
			Dm = k[0];
			return Dm;
		}
		/* all following cases are handled in CalcTracerDiffusionCoefficient1, and the needed parameters are gathered
		 * below */
		case 2:
		case 3:
		case 4:
		case 5:
		case 6:
		case 7:
/*OK411
        {
        count_nodes = ElNumberOfNodes[ElGetElementType(index) - 1];
        element_nodes = ElGetElementNodes(index);
      p_ind = 1;
      sprintf(name_pres,"PRESSURE%ld",transport_phase+1);
   // SB:GS4 todo               p_idx = m_pcs->PCSGetNODValueIndexNew(name_pres,1);
        if (p_ind) {
            pressure_average = 0.0;
            for (i = 0; i < count_nodes; i++)
                pressure_average += GetNodeVal(element_nodes[i], p_idx);
   pressure_average /= count_nodes;
   }
   t_ind = 1;
   if (! GetRFProcessHeatReactModel()){
   DisplayMsgLn("No heat transport active, no temperature available for determination of diffusion coefficient ");
   break;
   };

   sprintf(name_heat,"TEMPERATURE%ld",transport_phase+1);
   // SB:GS4 todo		m_pcs->PCSGetNODValueIndexNew(name_heat,1);
   if (t_ind) {
   temperature_average = 0.0;
   for (i = 0; i < count_nodes; i++)
   temperature_average += GetNodeVal(element_nodes[i], t_idx);
   temperature_average /= count_nodes;
   }
   eta = mfp_vector[transport_phase]->Viscosity(); //GetFluidViscosity(transport_phase, index, 0., 0., 0., theta);
   diffusion_coefficient = CalcDiffusionCoefficientCP_Method1(index, temperature_average, pressure_average, eta);
   element_nodes = NULL;
   break;
   }
 */
#ifdef GEM_REACT
		case 8: /* Archies law De = Dp * tau * tau_zero* poros   with tau = poros^m/poros0^m   ....tortuosity = tau_zero
		           * tau depends on porosity and tau_zero is multiplicated  (like porosity) outside this function  */
		{
			MeshLib::CElem* m_Elem;

			if (count_of_diffusion_model_values < 2)
				return 0.0;

			porosity = m_mat_mp->Porosity(index, theta);

			m_Elem = m_pcs->m_msh->ele_vector[index];
			// we average the diffusion coefficient directly
			count_nodes = m_Elem->GetNodesNumber(false);
			diffusion_average = 0.0;
			// we need saturation for multi-phase flow: ***  attention: ***  only coupling to Richards flow implemented
			// so far
			flowflag = m_vec_GEM->REACT_GEM::GetFlowType_MT();
			if (flowflag == 3)
				f_pcs = PCSGet("RICHARDS_FLOW");

			for (i = 0; i < count_nodes; i++) // calculate harmonic mean of node based diffusion coefficients
			{
				// then get the values from nodes
				dummy = m_vec_GEM->REACT_GEM::GetNodePorosityValue(m_Elem->GetNodeIndex(i));
				// for Multi-Phase flow we also have to account for the saturation......the following line would account
				// for element saturatuions all flow types
				if (flowflag == 3) // RICHARDS_FLOW ...get node based saturations
				{
					saturation = f_pcs->GetNodeValue(m_Elem->GetNodeIndex(i),
					                                 f_pcs->GetNodeValueIndex("SATURATION1")
					                                     + 1); // current volume of water phase after hydraulic step
				}
				else
				{
					saturation = 1.0;
				}
				// Attention....this does not work for gas phase/non wetting phase ....tr_phase=10
				dummy *= saturation; // account for "mean element" saturation in Archies law ....better would be to get
				// node saturations!
				Dm = k[0] * pow(dummy, k[1])
				     * dummy; // node based diffusion coefficient---last porosity is for effective diffusion coefficient
				// then get the values for initial porosity from nodes
				dummy = m_vec_GEM->REACT_GEM::GetNodePorosityValueInitial(
				    m_Elem->GetNodeIndex(i)); // ATTENTION: This does not include initial saturation different from 1
				Dm /= pow(dummy, k[1]); // node based diffusion coefficient with modified Archies relation

				diffusion_average += 1.0 / Dm;
				//	cout << "debug: " << Dm << " porosity: " << GetNodePorosityValue_MT(m_Elem->GetNodeIndex ( i ), 0)
				//<< endl;
			}
			Dm = count_nodes / diffusion_average; // This is now harmonic mean of node diffusion coefficients
			// end calculation of diffusion coefficient
			if (flowflag == 3) // RICHARDS_FLOW ...get node based saturations
			{
				saturation = PCSGetEleMeanNodeSecondary_2(
				    index, m_pcs->flow_pcs_type, "SATURATION1", 1); // this line accounts for element saturatuions for
				// all flow types....of course we have a problem if
				// we have saturation ne 1 and not Richards flow
				if (saturation <= 1.0e-20)
					saturation = 1.0e-20; // set to an arbitrary small number to avoid divsion by zero...
			}
			else
			{
				saturation = 1.0;
			}
			Dm = Dm / (porosity * saturation); // correct for multiplication with element porosities & saturations ->
			// Pore diffusion coefficient
			//			cout << " CalcDiffusionCoefficientCP: De: " << Dm*porosity*saturation << " saturation " <<
			// saturation << " porosity " << endl;

			return Dm;
		}
		case 9: /*  De is calculated independently from element porosity. We use node porosity values with with Archies
		           law De = Dp * poros^m   and do a harmonic average of the node diffusion coefficients!   as Dp is part
		           of the dispersion tensor, we  modify Dp -> Dp=Dp0*poro^(m-1)*/
		{ /* Remark: tortuosity is considered as constant factor and does not change with porosity */
			MeshLib::CElem* m_Elem;

			if (count_of_diffusion_model_values < 2)
				return 0.0;

			porosity = m_mat_mp->Porosity(index, theta);

			m_Elem = m_pcs->m_msh->ele_vector[index];
			// we average the diffusion coefficient directly
			count_nodes = m_Elem->GetNodesNumber(false);
			diffusion_average = 0.0;
			// we need saturation for multi-phase flow: ***  attention: ***  only coupling to Richards flow implemented
			// so far
			flowflag = m_vec_GEM->REACT_GEM::GetFlowType_MT();
			if (flowflag == 3)
				f_pcs = PCSGet("RICHARDS_FLOW");

			for (i = 0; i < count_nodes; i++) // calculate harmonic mean of node based diffusion coefficients
			{
				// then get the values from nodes
				dummy = m_vec_GEM->REACT_GEM::GetNodePorosityValue(m_Elem->GetNodeIndex(i));
				// for Multi-Phase flow we also have to account for the saturation......the following line would account
				// for element saturatuions all flow types
				if (flowflag == 3) // RICHARDS_FLOW ...get node based saturations
				{
					saturation = f_pcs->GetNodeValue(m_Elem->GetNodeIndex(i),
					                                 f_pcs->GetNodeValueIndex("SATURATION1")
					                                     + 1); // current volume of water phase after hydraulic step
				}
				else
				{
					saturation = 1.0;
				}
				// Attention....this does not work for gas phase/non wetting phase ....tr_phase=10
				dummy *= saturation; // account for "mean element" saturation in Archies law ....better would be to get
				// node saturations!
				Dm = k[0] * pow(dummy, k[1]); // node based diffusion coefficient

				diffusion_average += 1.0 / Dm;
				//	cout << "debug: " << Dm << " porosity: " << GetNodePorosityValue_MT(m_Elem->GetNodeIndex ( i ), 0)
				//<< "\n";
			}
			Dm = count_nodes / diffusion_average; // This is now harmonic mean of node diffusion coefficients
			// end calculation of diffusion coefficient
			if (flowflag == 3) // RICHARDS_FLOW ...get node based saturations
			{
				saturation = PCSGetEleMeanNodeSecondary_2(
				    index, m_pcs->flow_pcs_type, "SATURATION1", 1); // this line accounts for element saturatuions for
				// all flow types....of course we have a problem if
				// we have saturation ne 1 and not Richards flow
				if (saturation <= 1.0e-20)
					saturation = 1.0e-20; // set to an arbitrary small number to avoid divsion by zero...
			}
			else
			{
				saturation = 1.0;
			}
			Dm = Dm / (porosity * saturation); // correct for multiplication with element porosities & saturations ->
			// Pore diffusion coefficient
			//			cout << " CalcDiffusionCoefficientCP: De: " << Dm*porosity*saturation << " saturation " <<
			// saturation << " porosity " << "\n";

			return Dm;
		}
#endif
		case 10: // Temperature dependence Yaws
		{
			if (count_of_diffusion_model_values < 2)
				return 0.0;
			temperature_average = 0;
			size_t nnodes = m_pcs->m_msh->ele_vector[index]->GetNodesNumber(false);
			if (REACTINT_vec.size() > 0)
			{
				for (size_t i = 0; i < nnodes; i++)
				{
					long node = m_pcs->m_msh->ele_vector[index]->GetNodeIndex(i);
					temperature_average += REACTINT_vec[0]->GetTemperature(node) / double(nnodes);
				}
			}
			else if ((t_pcs = PCSGet("HEAT_TRANSPORT")))
			{
				// t_pcs = PCSGet("HEAT_TRANSPORT");
				int idx_t = t_pcs->GetNodeValueIndex("TEMPERATURE1") + 1;
				for (size_t i = 0; i < nnodes; i++)
				{
					long node = m_pcs->m_msh->ele_vector[index]->GetNodeIndex(i);
					temperature_average += t_pcs->GetNodeValue(node, idx_t) / double(nnodes);
				}
			}
			// else
			//  return 0;
			if (temperature_average > 0)
				Dm = pow(10, (k[0] + k[1] / temperature_average)) / 10000; // cm^2/s -> m^2/s;
			else
			{
				Dm = 0;
				DisplayMsgLn("Something wrong in diffusion model 10! T = 0.");
			}
			// A = -1.4943; //Daq(T) PCE
			// B = -1059.0; //Daq(T) PCE
			return Dm;
		}

		default:
			DisplayMsgLn("Unknown diffusion model!");
			break;
	}
	/* switch */

	return diffusion_coefficient;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: CalcDiffusionCoefficientCP_Method1
 */
/* Aufgabe:
   Berechnet in Abhaengigkeit der Konzentration die zur Isotherme
   gehoerige sorbierten Menge.
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E long index     :
   E long phase     :
   E long component :

   E double T       : Temperatur in Kelvin [K]
   E double p       : Druck in Pascal [Pa]
   E double Eta     : Dynamische Viskositaet des Wassers [K]
   Richtwerte:
   = 0.894 cP (25 Grad Celsius) [*1.e-3 N s m^(-2) bzw. J s m^(-3)]
   = 1.002 cP (20 Grad Celsius) [*1.e-3 N s m^(-2) bzw. J s m^(-3)]
   = 1.306 cP (10 Grad Celsius) [*1.e-3 N s m^(-2) bzw. J s m^(-3)]

   E double k1  : Faktor k1 der Isotherme
   E double k2  : Faktor k2 der Isotherme
   E int    iso : Isotherme-Identifikator
 */
/* Ergebnis:
   - double - wert der Sorption
 */
/* Programmaenderungen:
   06/2000    AH      Erste Version
   09/2004	   SB		  MOved to Class CompProperties as member function
 */
/**************************************************************************/
double CompProperties::CalcDiffusionCoefficientCP_Method1(long index, double T, double P, double eta)
{
	int diffusion_model = 0; // SB:cp GetTracerDiffusionModel(index, phase, component);
	double* k = NULL;
	int count = 0;
	double t = aktuelle_zeit;
	double msol, Vs;
	double Dm, Daq; // WW, Dg;
	double Kb, Rm;
	// WW int curve;
	double mg = 28.97; /* Mol.-Gew. Luft [g/mol] */
	double Vg = 20.1; /* molares Volumen von Luft [cm^3/mol] */

	double m; /* Molekulargewicht der diffundierenden Substanz [g/mol] */
	double V; /* molares Volumen der diffundierenden Substanz [cm^3/mol] */
	/* In erster Naeherung: Quotient Mol.-Gew./Dichte */

	double p = P / 101325.; /* Umrechnen P in Pascal in p in atm */

	index = index;
	if (t < 0.0 || eta < 0.0)
		return -1.;

	k = diffusion_model_values;

	switch (diffusion_model)
	{
		case 0: /* Keine Diffusion */
		{
			return 0.0;
		}
		case 1: /* Konstanter Diffusionswert */
		{
			if (count < 1)
				return 0.0;
			Dm = k[0];
			return Dm;
		}
		case 2: /* Variabler Diffusionswert (Zeitabhaengig). Abhaegigkeit ueber Kurve */
		{
			if (count < 1)
				return 0.0;
			Dm = k[0];
			// WW curve = (int) k[1];
			return Dm;
		}
		case 3: /* Worch, 1993 */
		{
			if (count < 1)
				return 0.0;
			m = k[0];
			Daq = 3.595e-7 * T / eta / pow(m, 0.53) * 1.e-4;
			return Daq;
		}
		case 4: /* Hayduk und Laudie, 1974 */
		{
			if (count < 1)
				return 0.0;
			V = k[0];
			Daq = 13.26e-5 / pow(eta, 1.14) / pow(V, 0.589) * 1.e-4;
			return Daq;
		}
		case 5: /* Wilke und Chang, 1955 */
		{
			if (count < 2)
				return 0.0;
			msol = k[0];
			Vs = k[0];
			Daq = 7.4e-8 * T * sqrt(msol) / eta / pow(Vs, 0.6) * 1.e-4;
			return Daq;
		}
		case 6: /* Stokes-Einstein (fuer Partikel/Makromolekuele) */
		{
			if (count < 1)
				return 0.0;
			Rm = k[0];
			Kb = 1.38066e-23; /* Boltzmann Konstante [J/K] */
			Daq = Kb * T / 6. / PI / Rm / eta * 1.e-4;
			return Daq;
		}
		case 7: /*  FSG-Methode, Lyman et al., 1990 */
		{
			if (count < 2)
				return 0.0;
			m = k[0];
			V = k[1];
			Vg = 20.1; /* Molares Volumen von Luft [cm?/mol] */
			// WW Dg = (0.001 * pow(T, 1.75) * sqrt(1. / mg + 1. / m)) / (p * pow(pow(Vg, 1. / 3.) + pow(V, 1. / 3.),
			// 2.)) * 1.e-4;
			return (0.001 * pow(T, 1.75) * sqrt(1. / mg + 1. / m))
			       / (p * MathLib::fastpow(pow(Vg, 1. / 3.) + pow(V, 1. / 3.), 2)) * 1.e-4; // WW
		}
	}
	DisplayMsgLn("Unknown diffusion model specified!");
	return 0.;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: GetNumberDiffusionValuesCompProperties
 */
/* Aufgabe:
   Liefert in Abhaengigkeit des Diffusionsmodells die Anzahl der
   benoetigten Koeffizienten.
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E int diffusion_model : Diffusionsmodell
 */
/* Ergebnis:
   - int - Anzahl der Koeffizienten
 */
/* Programmaenderungen:
   06/2000    AH      Erste Version  (Fall 0 bis 7)
   02/2004    SB      Adapted to CompProperties
 */
/**************************************************************************/
int CompProperties::GetNumberDiffusionValuesCompProperties(int diffusion_model)
{
	int n = -1;

	if (diffusion_model < -1)
		return n;

	switch (diffusion_model)
	{
		case -1:
			n = 0;
			break; /* Keine Diffusion */
		case 0: /* curve */
			n = 1;
		case 1:
			n = 1;
			break; /* Konstanter Diffusionswert */
		case 2:
			n = 1;
			break; /* Variabler Diffusionswert (Zeitabhaengig) */

		/* Diffusionskoeffizienten in Wasser Daq [m^2/s] */
		case 3:
			n = 1;
			break; /* Worch, 1993 */
		case 4:
			n = 1;
			break; /* Hayduk und Laudie, 1974 */
		case 5:
			n = 2;
			break; /* Wilke und Chang, 1955 */
		case 6:
			n = 1;
			break; /* Stokes-Einstein (fuer Partikel/Makromolekuele) */

		/* Diffusionskoeffizienten in Luft Dg [m^2/s] */
		case 7:
			n = 2;
			break; /* FSG-Methode, Lyman et al., 1990 */
		case 8:
			n = 2;
			break; /* Archies Law */
		case 9:
			n = 2;
			break; /* Archies Law */
		case 10:
			n = 2;
			break; /* Yaws empirical model*/

	} /* switch */

	/* switch */

	return n;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: GetNumberDecayValuesCompProperties
 */
/* Aufgabe:
   Liefert in Abhaengigkeit des in der sorbierten Phase Zerfallsmodells
   die Anzahl der benoetigten Koeffizienten.
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E int isotherm : Isotherme-Identifikator
 */
/* Ergebnis:
   - int - Anzahl der Isothermen-Koeffizienten
 */
/* Programmaenderungen:
   06/2000    AH      Erste Version (Fall 0 bis 7)
   09/2004	  SB	  Adapted to CompProperties
 */
/**************************************************************************/
int CompProperties::GetNumberDecayValuesCompProperties(int decay_model)
{
	int n = -1;

	if (decay_model < -1)
		return n;

	switch (decay_model)
	{
		case -1:
			n = 0;
			break;
		case 0:
			n = 1;
			break; /* No Decay */
		case 1:
			n = 2;
			break; /* Any Order Decay with constant Decay Rate */
		case 2:
			n = 2;
			break; /* Monod or Michaelis-Menten kinetics with constant rate coefficients */
		default:
			DisplayMsgLn(" Error: Unknown model for decay  ");
			break;
	} /* switch */

	return n;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: GetNumberIsothermValuesCompProperties
 */
/* Aufgabe:
   Liefert in Abhaengigkeit der Isotherme die Anzahl der Koeffizienten.
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E int isotherm : Isotherme-Identifikator
 */
/* Ergebnis:
   - int - Anzahl der Isothermen-Koeffizienten
 */
/* Programmaenderungen:
   10/1999     AH         Erste Version  (Fall 0 bis 11)
   02/2004     SB         Adapted to CompProperties - nur method 0-3+15
   04/2006     CMCD       Included for analytical fracture surface method, Tang et al. 1981, WRR.
 */
/**************************************************************************/
int CompProperties::GetNumberIsothermValuesCompProperties(int isotherm)
{
	int n = 0;

	if (isotherm < -1)
		return -isotherm;

	switch (isotherm)
	{
		case -1:
			n = 0; /* no isotherm */
		case 0:
			n = 1;
			break; /* get KD from curve derivative */
		case 1:
			n = 1;
			break; /* Henry Isotherm */
		case 2:
			n = 2;
			break; /* Freundlich Isotherm */
		case 3:
			n = 2;
			break; /* Langmuir Isotherm */
		case 4:
			n = 1;
			break; /* Linear Isotherm, for fracture surface CMCD */
		case 5:
			n = 3;
			break; /* two-rate model */
		/*
		    case 4:
		        n = 3;  break;                  // Freundlich Langmuir Isotherm
		    case 5:
		        n = 4;  break;                  // Double Langmuir Isotherm
		    case 6:
		        n = 3;  break;                  // Extented Freundlich Isotherm
		    case 7:
		        n = 3;  break;                  // Gunary Isotherm
		    case 8:
		        n = 3;  break;                  // Fitter-Sutton Isotherm
		   case 9:
		   n = 4;  break;                  // Barry Isotherm
		   case 10:
		   n = 2;  break;                  // Power Isotherm
		   case 11:
		   n = 2;  break;                  // Modified Kielland Isotherm
		 */
		/*	case 15:
		      n=1;    break;
		 */
		default:
			DisplayMsgLn(" Error - this ISOTHERM model found ");
			break;
	} /* switch */

	return n;
}

/**************************************************************************
   ROCKFLOW - Funktion: CalcElementRetardationFactorNew

   Aufgabe:
   Berechnet den Retardationsfaktor des angegebenen Elements (fuer
   1D, 2D und 3D identisch)

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E long number: Index des Elements, dessen Retardationsfaktor berechnet
                  werden soll

   Ergebnis:
   Retardationsfaktor

   Programmaenderungen:
   04/2003     SB         Erste Version
   09/2004	   SB		  Moved to Class CompProperties as member function
**************************************************************************/
double CompProperties::CalcElementRetardationFactorNew(long index, double* gp, CRFProcess* m_pcs)
{
	static double porosity, density_rock, isotherm;
	static double conc, retard = 0.0;
	double theta = m_pcs->m_num->ls_theta;
	int gueltig;
	long group = 0; // SB4200 ToDO
	double fracture_width = 0.0; // CMCD needed for wall retardation in a fracture
	gp = gp; // OK411

	group = m_pcs->m_msh->ele_vector[index]->GetPatchIndex();

	CMediumProperties* m_mat_mp = NULL;
	m_mat_mp = mmp_vector[group];

	CSolidProperties* m_msp = NULL;
	m_msp = msp_vector[group];

	// porosity = GetSoilPorosity(index);
	porosity = m_mat_mp->Porosity(index, theta);
	// porosity = 0.5; //SB- set porosity
	// density_rock = 2000.0; // GetSolidDensity(index);
	density_rock = fabs(m_msp->Density());
	// SB-todo: CHeck saturatiuon dependence

	/* Get mean element concentration from last time step */
	conc = CalcElementMeanConcNew(index, m_pcs);
	/* DisplayMsg(" Mean conc: "); DisplayDouble(conc,0,0); DisplayMsgLn(" "); */

	switch (isotherm_model)
	{
		case -1: /* no sorption */
			isotherm = 0.0;
			// CMCD moved here
			retard = 1. + (1. - porosity) * density_rock * isotherm / porosity;
			break;
		case 0: /* from curve */
			isotherm = 0.0;
			isotherm = GetCurveDerivative((int)isotherm_function_name, 0, fabs(conc), &gueltig);
			// CB added here 31.10.07
			retard = 1. + (1. - porosity) * density_rock * isotherm / porosity;
			break;
		case 1: /* linear isotherm   */
			isotherm = isotherm_model_values[0];
			//  retard = 1. + density_rock*isotherm/porosity;//CMCD moved here
			// CMCD moved here
			retard = 1. + (1. - porosity) * density_rock * isotherm / porosity;
			break;
		/* rausgenommen, solange der PRozess linear ist */
		case 2: /* Freundlich Isotherm */
			/* Freundlich is linear */
			if (fabs(isotherm_model_values[1] - 1.0) < MKleinsteZahl)
				isotherm = isotherm_model_values[0];
			else if (fabs(conc) < MKleinsteZahl)
				isotherm = isotherm_model_values[0];
			else
				isotherm = isotherm_model_values[0] * isotherm_model_values[1]
				           * pow(fabs(conc), isotherm_model_values[1] - 1.0);
			// CMCD moved here
			retard = 1. + (1. - porosity) * density_rock * isotherm / porosity;
			break;
		case 3: /* Langmuir Isotherm */
			isotherm = isotherm_model_values[0] * isotherm_model_values[1]
			           / MathLib::fastpow((1 + isotherm_model_values[0] * fabs(conc)), 2);
			// CMCD moved here
			retard = 1. + (1. - porosity) * density_rock * isotherm / porosity;
			break;
		case 4: /* Face retardation model*/
			isotherm = isotherm_model_values[0];
			fracture_width = m_pcs->m_msh->ele_vector[index]->GetFluxArea();
			// Tang et al. 1981, Contaminant transport in fractured porous media, anaylitacl solution for a single
			// fracture, WRR 17, 3, 555-564.
			retard = 1 + (isotherm / (fracture_width / 2.0));
			break;
		case 5: /* two-rate model */
			break; // YS: the two-rate model is described in the RWPT
		// case 15: /* Input by curve */
		//	 isotherm = 0.0;
		//	 isotherm = GetCurveDerivative((int) isotherm_model_values[0], 0, fabs(conc), &gueltig);
		//	 break;
		default:
			DisplayMsgLn("Unknown sorption isotherm type. Assuming no sorption");
			isotherm = 0.0;
			// CMCD moved here
			retard = 1. + (1. - porosity) * density_rock * isotherm / porosity;
			break;
	}
	// DisplayMsg(" conc: "); DisplayDouble(conc,0,0);DisplayMsg(" isotherm: "); DisplayDouble(isotherm,0,0);
	// DisplayMsgLn(" ");
	// if(conc < 0.0) isotherm = 0.0;

	//  retard = 1. + (1.-porosity)*density_rock*isotherm/porosity; Case 4 doesn't use this function.
	//  if(index < 1) {DisplayMsg(" Retardation factor: "); DisplayDouble(retard,0,0); DisplayMsgLn(" ");}

	return retard;
}

/**************************************************************************
   ROCKFLOW - Funktion: CalcElementMeanConcentration

   Aufgabe:
   Berechnet den mittleren Wert eines Elements aus den Werten der
   angrenzenden Knoten

   Programmaenderungen:
   04/2003     SB         Erste Version
   09/2004	   SB		  Moved to Class CompProperties as member function
**************************************************************************/
double CompProperties::CalcElementMeanConcNew(long index, CRFProcess* m_pcs)
{
	static long i, nn, idx;
	//  static long *element_nodes;
	static double val, val1, val2;
	double theta; // GetNumericalTimeCollocation("TRANSPORT");
	//  CFEMesh* m_msh = m_pcs->m_msh; // Get mesh from Process
	MeshLib::CElem* elem = NULL;

	elem = m_pcs->m_msh->ele_vector[index];
	nn = elem->GetVertexNumber();
	idx = m_pcs->GetNodeValueIndex((char*)compname.data());
	theta = m_pcs->m_num->ls_theta;
	/* get Value after last iteration */
	/*
	   val1=0.0;
	   if (ElGetElement(index)!=NULL)   // wenn Element existiert
	     if (ElGetElementActiveState(index)){  // nur aktive Elemente
	        nn = ElGetElementNodesNumber(index);
	        element_nodes = ElGetElementNodes(index);
	        for (i=0;i<nn;i++) {
	   //			val1 += PCSGetNODConcentration(element_nodes[i],component,timelevel);
	         val1 += PCSGetNODValue(element_nodes[i], (char *) name.data(),1); //name is class member, new timelevel
	      }
	        val1 /= (double)nn;
	   element_nodes = NULL;
	   }
	 */
	val1 = 0.0;
	for (i = 0; i < nn; i++)
		val1 += m_pcs->GetNodeValue(elem->GetNodeIndex(i), idx);
	val1 /= (double)nn;

	/* get value from beginning of timestep */
	/*
	   val2=0.0;
	   if (ElGetElement(index)!=NULL)   // wenn Element existiert
	     if (ElGetElementActiveState(index)){  // nur aktive Elemente
	        nn = ElGetElementNodesNumber(index);
	        element_nodes = ElGetElementNodes(index);
	        for (i=0;i<nn;i++) {
	   //			val2 += ergebnis[element_nodes[i]];
	         val2 += PCSGetNODValue(element_nodes[i],(char *) name.data(),0);
	      }
	        val2 /= (double)nn;
	   element_nodes = NULL;
	   }
	 */
	val2 = 0.0;
	for (i = 0; i < nn; i++)
		val2 += m_pcs->GetNodeValue(elem->GetNodeIndex(i), idx + 1);
	val2 /= (double)nn;
	/* calculate mean value */

	//	 DisplayMsgLn(" "); DisplayMsg(" val1: "); DisplayDouble(val1,0,0);DisplayMsg(", val2: ");
	// DisplayDouble(val2,0,0); DisplayMsgLn("");
	val = theta * val1 + (1.0 - theta) * val2;

	return val;
}

/**************************************************************************
   ROCKFLOW - Funktion: MTM2CalcElementDecayRate

   Aufgabe:
   Berechnet die Zerfallsrate des angegebenen Elements (fuer
   1D, 2D und 3D identisch)

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E long number: Index des Elements, dessen Retardationsfaktor berechnet
                  werden soll

   Ergebnis:
   Zerfallsrate

   Programmaenderungen:
   04/2003     SB         Erste Version
   09/2004	   SB		  Moved to Class CompProperties as member function
**************************************************************************/
double CompProperties::CalcElementDecayRateNew(long index, CRFProcess* m_pcs)
{
	static double conc, lambda = 0.0;
	int gueltig;

	// if(index >= anz_active_elements){DisplayMsgLn(" Too many elements "); return 0.0; }
	/* Get mean element concentration from last time step */
	conc = CalcElementMeanConcNew(index, m_pcs);

	switch (decay_model)
	{
		case -1: /* no decay */
			lambda = 0.0;
			break;
		case 0: /* curve - linear interpolation*/
			lambda = 0.0;
			lambda = GetCurveValue((int)decay_function_name, 0, fabs(conc), &gueltig);
			break;
		case 1: /* Any order decay  with n > 1 */
			if (fabs(decay_model_values[1] - 1.0) < MKleinsteZahl)
			{
				/* First order decay */
				lambda = decay_model_values[0];
				break;
			}
			/* zero - order decay */
			else if (decay_model_values[1] < MKleinsteZahl)
			{
				if (fabs(conc) < MKleinsteZahl)
					lambda = decay_model_values[0];
				else
					lambda = decay_model_values[0] / conc;
			}
			else
			{
				/* Any order decay with n not equal to one or zero*/
				if (fabs(conc) < MKleinsteZahl)
					lambda = decay_model_values[0];
				else
					lambda = decay_model_values[0] * pow(fabs(conc), decay_model_values[1] - 1.0);
			}
			break;
		case 2: /* Monod-Kinetics */
			lambda = decay_model_values[0] / (decay_model_values[1] + conc);
			break;
		default:
			lambda = 0.0;
			break;
	}
	// if(index < 0){DisplayMsg(" Decay Rate lambda: "); DisplayDouble(lambda,0,0); DisplayMsgLn(" ");}
	return lambda;
}

// SB:todo Wie kann ich die gut ersetzen (wird nur in loop_pcs gebraucht, um zu schauen ob der process mobil ist ??
int CPGetMobil(long comp)
{
	int mobil = -1;
	CompProperties* cp_single = NULL;

	cp_single = cp_vec[comp];
	if (cp_single == NULL)
	{
		DisplayMsgLn("The requested component properties are not defined!");
		return 0;
	}
	mobil = cp_single->mobil;
	if (mobil < 0)
		DisplayMsgLn("Error getting component property");

	return mobil;
}

/* Kommentar:

   Eingabe for Schl?selw?ter:

   DECAY

   Method 0: curve
   Give curve number - not yet implemented
     0 ; no decay

   Method 1: decay with any order: includes orders out of [0, inf] positive
   Give Method  decay_rate  order_of_decay
   1  1.0e-6  1.0  ; decay with first order kinetics
   1  1.0e-06 2.0  ; decay with second order kinetics

   Method 2:  Monod-Kinetics
   Give Method degradation_rate half_saturation_concentration
   2  1.0e-6  0.5 ; Monod Kinetics with dC/dt = C * (1.0e-6/ 0.5 + c)

   $ISOTHERM
   Method 0: Input by curve
   Give: curve_number
   0  7 ; curve derivative is taken at the element_concentration

   Method 1: linear isotherm
   Give: Method linear_distribution_coefficient
   1  1e-3;

   Method 2: Freundlich Isotherm
   Give: Method Freundlich_coefficient freundlich_exponent
   2  1.0e-4  0.8

   Method 3: Langmuir Isotherm
   Give: param1  param2
   3   1e-3 0.5

 */

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   01/2005 OK Implementation
   last modified:
**************************************************************************/
void MCPDelete()
{
	long i;
	int no_mcp = (int)cp_vec.size();
	for (i = 0; i < no_mcp; i++)
		delete cp_vec[i];
	cp_vec.clear();
}
