/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   FEMLib - Object: MFP Fluid Properties
   Task:
   Programing:
   08/2004 OK Implementation
   last modified:
**************************************************************************/
#include "makros.h"

#include <cfloat>

// FEM-Makros
//#include "mathlib.h"
#include "eos.h" //NB
// GeoSys-GeoLib
#include "files0.h"
// GeoSys-FEMLib
#include "fem_ele_std.h"
//
#include "rf_mfp_new.h"
//#include "rf_mmp_new.h"
extern double InterpolValue(long number, int ndx, double r, double s, double t);
//#include "rf_pcs.h"
#include "rfmat_cp.h"
extern double GetCurveValue(int, int, double, int*);
#include "tools.h" //GetLineFromFile
// KG44 used for GEMS part
#ifdef GEM_REACT
#include "rf_REACT_GEM.h"
#endif

#include "PhysicalConstant.h"

/* Umrechnungen SI - Amerikanisches System */
// WW #include "steam67.h"
#define PSI2PA 6895.
#define PA2PSI 1.4503263234227701232777374909355e-4

double gravity_constant = 9.81; // TEST for FEBEX OK 9.81;

using namespace PhysicalConstant;
using namespace std;

//==========================================================================
std::vector<CFluidProperties*> mfp_vector;

double TemperatureUnitOffset()
{
	return process::isTemperatureUnitCesius() ? PhysicalConstant::CelsiusZeroInKelvin : 0.0;
}

/**************************************************************************
   FEMLib-Method:
   Task: OBJ constructor
   Programing:
   08/2004 OK Implementation
**************************************************************************/
CFluidProperties::CFluidProperties() : name("WATER"),
	_reference_temperature(PhysicalConstant::CelsiusZeroInKelvin + 20.0)
{
	phase = 0;
	// Density
	density_model = 1;
	rho_0 = 1000.;
	drho_dp = 0.;
	drho_dT = 0.;
	drho_dC = 0.;
	// Viscosity
	viscosity_model = 1;
	my_0 = 1e-3;
	dmy_dp = 0.;
	dmy_dT = 0.;
	dmy_dC = 0.;
	// Thermal properties
	heat_capacity_model = 1;
	specific_heat_capacity = 4680.; // CMCD we should give this as SHC not HC GeoSys 4 9/2004
	heat_conductivity_model = 1;
	heat_conductivity = 0.6;
	// Electrical properties
	// Chemical properties
	diffusion_model = 1;
	diffusion = 2.13e-6;
	// State variables
	p_0 = 101325.;
	T_0 = PhysicalConstant::CelsiusZeroInKelvin + 20.0;
	C_0 = 0.;
	Z = 1.;
	cal_gravity = true;
	drho_dT_unsaturated = false; // considering fluid expansion due to temperature in the unsaturated case? (Richards)
	// Data
	mode = 0; // Gauss point values
	Fem_Ele_Std = NULL;
	// WW
	molar_mass = MolarMass::Air;

	compressibility_model_pressure = -1;
	specific_heat_source = 0.0;
	beta_T = 0.0;
	cmpN = 0;

	fluid_id = 1; // Water

#ifdef MFP_TEST // WW
	scatter_data = NULL;
#endif
}

/**************************************************************************
   FEMLib-Method:
   Task: OBJ deconstructor
   Programing:
   05/2010 OK/AKS Implementation
**************************************************************************/
CFluidProperties::~CFluidProperties(void)
{
	for (int i = 0; i < (int)component_vector.size(); i++)
		component_vector[i] = NULL;
	component_vector.clear();

#ifdef MFP_TEST
	if (scatter_data) // WW
		delete scatter_data;
#endif
}

/**************************************************************************
   FEMLib-Method:
   Task: OBJ read function
   Programing:
   08/2004 OK Implementation
   11/2004 SB string streaming
**************************************************************************/
std::ios::pos_type CFluidProperties::Read(std::ifstream* mfp_file)
{
	std::string sub_line;
	std::string line_string;
	std::string delimiter(" ");
	bool new_keyword = false;
	std::string hash("#");
	std::ios::pos_type position;
	std::string sub_string;
	// WW bool new_subkeyword = false;
	std::string dollar("$");
	std::string delimiter_type(":");
	std::stringstream in;
	//========================================================================
	// Schleife ueber alle Phasen bzw. Komponenten
	while (!new_keyword)
	{
		// WW new_subkeyword = false;
		position = mfp_file->tellg();
		// SB    mfp_file->getline(buffer,MAX_ZEILE);
		// SB    line_string = buffer;
		line_string = GetLineFromFile1(mfp_file);
		if (line_string.size() < 1)
			break;
		if (line_string.find(hash) != string::npos)
		{
			new_keyword = true;
			break;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$FLUID_TYPE") != string::npos)
		{
			in.str(GetLineFromFile1(mfp_file));
			in >> name; // sub_line
			switch (name[0])
			{
				case 'C': fluid_id = 0; break; // CARBON DIOXIDE
				case 'W': fluid_id = 1; break; // WATER
				case 'M': fluid_id = 2; break; // METHANE
				case 'N': fluid_id = 3; break; // NITROGEN
				case 'H': fluid_id = 4; break; // HYDROGEN
				case 'O': fluid_id = 5; break; // OXYGEN
				default:
				{
					std::cout << "Fluid type of " << name << " is not valid. Take the default type of WATER" << std::endl;
					fluid_id = 1;
				}
			}

			in.clear();
			continue;
		}
		//....................................................................
		// AKS
		if (line_string.find("$COMPONENTS") != string::npos)
		{
			in.str(GetLineFromFile1(mfp_file));
			in >> cmpN >> cmpNm1 >> cmpNm2 >> cmpNm3 >> cmpNm4; // get up to four component
			in.clear();
			continue;
		}
		//....................................................................
		// NB 4.8.01
		if (line_string.find("$FLUID_NAME") != string::npos)
		{
			in.str(GetLineFromFile1(mfp_file));
			in >> fluid_name; // sub_line
			therm_prop(fluid_name); // NB 4.9.05 (getting thermophysical constants of specified substance)
			// TODO: add choosing of property functions in input file (NB)
			in.clear();
			continue;
		}
		//....................................................................
		if (line_string.find("$EOS_TYPE") != string::npos)
		{
			in.str(GetLineFromFile1(mfp_file));
			in >> eos_name; // sub_line
			in.clear();
			continue;
		}
		//....................................................................
		// NB Oct-2009
		if (line_string.find("$COMPRESSIBILITY") != string::npos)
		{
			in.str(GetLineFromFile1(mfp_file));
			in >> compressibility_model_pressure; // sub_line 1 for first phase
			if (compressibility_model_pressure == 15) // components constant density
			{
				if (eos_name == "CONSTANT")
					in >> beta_p[0] >> beta_p[1] >> beta_p[2] >> beta_p[3];
			}
			else
				in >> compressibility_pressure; // sub_line 1
			in.clear();

			in.str(GetLineFromFile1(mfp_file));
			in >> compressibility_model_temperature; // sub_line 2 for second phase
			if (compressibility_model_temperature == 15) // components constant density
			{
				if (eos_name == "CONSTANT")
					in >> alpha_T[0] >> alpha_T[1] >> alpha_T[2] >> alpha_T[3];
			}
			else
				in >> compressibility_temperature; // sub_line 2
			in.clear();

			in.str(GetLineFromFile1(mfp_file));
			in >> solutal_expansivity_model; // sub_line 2 for second phase
			in >> solutal_expansivity; // sub_line 2
			in.clear();

			// available models see CFluidProperties::drhodP and CFluidProperties::drhodT
			// 0 incompressible
			// 1 constant slope
			// 2 slope from fct_table
			// 3 difference quotient
			// 4 analytical derivation

			in.clear();
			continue;
		}
		//....................................................................
		if (line_string.find("$JTC") != string::npos)
		{
			in.str(GetLineFromFile1(mfp_file));
			in >> mu_JT; // sub_line
			in.clear();
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$DAT_TYPE") != string::npos)
		{
			in.str(GetLineFromFile1(mfp_file));
			in >> name; // sub_line
			in.clear();
			continue;
		}
		// YD/WW subkeyword found
		if (line_string.find("$NON_GRAVITY") != string::npos)
		{
			cal_gravity = false;
			continue;
		}
		if (line_string.find("$DRHO_DT_UNSATURATED")
		    != string::npos) // JM considering drho/dT for the unsaturated case (richards)?
		{
			drho_dT_unsaturated = true; // if keyword found, drho/dT will be considered for unsaturated case (richards)
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$DENSITY") != string::npos)
		{
			// WW new_subkeyword = false;
			in.str(GetLineFromFile1(mfp_file));
			in >> density_model;

			// TF - _rho_fct_name is only used for writing it back to the file
			if (density_model == 0) // rho = f(x)
			{ // JOD 2014-11-10
				in >> density_curve_number;
				density_pcs_name_vector.push_back("TEMPERATURE1");
			}
			if (density_model == 1) // rho = const
				in >> rho_0;

			if (density_model == 2) // rho(p) = rho_0*(1+beta_p*(p-p_0))
			{
				in >> rho_0;
				in >> p_0;
				in >> drho_dp;
				density_pcs_name_vector.push_back("PRESSURE1");
			}
			if (density_model == 3) // rho(C) = rho_0*(1+beta_C*(C-C_0))
			{
				in >> rho_0;
				in >> C_0;
				in >> drho_dC;
				density_pcs_name_vector.push_back(
				    "CONCENTRATION1"); // JOD 2014-11-10, consistent with output DELTA_CONCENTRATION1
				// density_pcs_name_vector.push_back("Isochlor");  // PCH
			}
			if (density_model == 4) // rho(T) = rho_0*(1+beta_T*(T-T_0))
			{
				in >> rho_0;
				in >> T_0;
				T_0 += TemperatureUnitOffset();
				in >> drho_dT;
				density_pcs_name_vector.push_back("TEMPERATURE1");
			}
			if (density_model == 5) // rho(C,T) = rho_0*(1+beta_C*(C-C_0)+beta_T*(T-T_0))
			{
				in >> rho_0;
				in >> C_0;
				in >> drho_dC;
				in >> T_0;
				T_0 += TemperatureUnitOffset();
				in >> drho_dT;
				density_pcs_name_vector.push_back("CONCENTRATION1");
				density_pcs_name_vector.push_back("TEMPERATURE1");
			}
			if (density_model == 6 || density_model == 14) // rho(p,T) = rho_0*(1+beta_p*(p-p_0)+beta_T*(T-T_0))
			{
				in >> rho_0;
				in >> p_0;
				in >> drho_dp;
				in >> T_0;
				T_0 += TemperatureUnitOffset();
				in >> drho_dT;
				density_pcs_name_vector.push_back("PRESSURE1");
				density_pcs_name_vector.push_back("TEMPERATURE1");
			}
			if (density_model == 7) // rho(p,p_v,T)
			{
				// no input data required
			}
			if (density_model == 8) // rho(p,T,C)
			{
				in >> C_0;
				density_pcs_name_vector.push_back("PRESSURE1");
				density_pcs_name_vector.push_back("TEMPERATURE1");
			}
			if (density_model == 9) // WW
				// Molar mass
				in >> molar_mass;

			if ((density_model == 10) // NB 4.8.01  read density from a rho-P-T table
			    || (density_model == 11) // NB 4.9.05  Peng-Robinson Equation of State
			    || (density_model == 12) // NB 4.9.05  Redlich-Kwong Equation of State
			    || (density_model == 13)) // NB JUN 09  Fundamental equation
			{
				std::string arg1, arg2, arg3;
				in >> arg1 >> arg2 >> arg3; // get up to three arguments for density model
				if (arg1.length() > 0)
				{
					if (isdigit(arg1[0]) != 0) // first argument is reference temperature
					{
						T_0 = atof(arg1.c_str());
						T_0 += TemperatureUnitOffset();
						arg1 = arg2;
						arg2 = arg3;
					}
				}
				else
					T_0 = TemperatureUnitOffset();

				if (arg1.length() == 0) // if no arguments are given use standard
				{
					arg1 = "PRESSURE1";
					arg2 = "TEMPERATURE1";
				}
				else if (arg2.length() == 0) // if only PRESSURE argument is given
					arg2 = "TEMPERATURE1";

				density_pcs_name_vector.push_back(arg1);
				if (T_Process)
					density_pcs_name_vector.push_back(arg2);
			}
			// AKS
			if (density_model == 15) // components constant density
			{
				if (eos_name == "CONSTANT")
					in >> rho[0] >> rho[1] >> rho[2] >> rho[3];
			}

			if (density_model == 18) // BG, NB calculated node densities from the phase transition model
			{
			}
			/// Density model no. 19: Extracts fluid densities from GEM3K calculations
			/// The densites depend on the thermodynamic data base (including pressure and temperature corrections)
			/// in addition various gas/fluid models (&EOS) are defined in GEMS3K
			if (density_model == 19) // KG44 get node densities from GEMS calculations
			{
			}
			if (density_model == 20) // rho(C) = rho (p, T)*(1+beta_C*(C-C_0))  for water, range p < 100 MPa, 0 <= T <=
			// 350 °C   Magri GFZ thesis
			{ // JOD 2014-11-10
				in >> C_0;
				in >> drho_dC;
				density_pcs_name_vector.push_back("PRESSURE1");
				density_pcs_name_vector.push_back("TEMPERATURE1");
				density_pcs_name_vector.push_back("CONCENTRATION1");
			}
			if (density_model == 26) // && pcs_vector[0]->getProcessType() != FiniteElement::TNEQ)
			{
				std::cout << "Warning: This density model requires two components and their molar masses defined in "
				             "the mcp file!\n";
				density_pcs_name_vector.push_back("PRESSURE1");
				density_pcs_name_vector.push_back("TEMPERATURE1");
				density_pcs_name_vector.push_back("CONCENTRATION1");
			}
			//      mfp_file->ignore(MAX_ZEILE,'\n');
			in.clear();
			continue;
		}
		if (line_string.find("$TEMPERATURE") != string::npos) // subkeyword found 11/2010, BG, NB, DL, SB
		{
			in.str(GetLineFromFile1(mfp_file));
			in >> T_0 >> T_0;
			T_0 += TemperatureUnitOffset();
			in.clear();
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$VISCOSITY") != string::npos)
		{
			in.str(GetLineFromFile1(mfp_file));
			in >> viscosity_model;
			// TF 11/2011 - used only in read- and write-method
			if (viscosity_model == 0) // my = fct(x)
				in >> viscosity_curve_number; // JOD 2014-11-10
			//				in >> _my_fct_name;
			if (viscosity_model == 1) // my = const
				in >> my_0;

			if (viscosity_model == 2) // my(p) = my_0*(1+gamma_p*(p-p_0))
			{
				in >> my_0;
				in >> p_0;
				in >> dmy_dp;
				viscosity_pcs_name_vector.push_back("PRESSURE1");
			}
			if (viscosity_model == 3) // my(T), Yaws et al. (1976)
			{ // optional: read reference temperature for viscosity model
				viscosity_pcs_name_vector.push_back("PRESSURE1"); // JM dummy wird benoetigt!
				// OK4704
				viscosity_pcs_name_vector.push_back("TEMPERATURE1");
			}
			if (viscosity_model == 4) // my(T), ???
			{
				viscosity_pcs_name_vector.push_back("PRESSURE1"); // JM dummy wird benoetigt!
				viscosity_pcs_name_vector.push_back("TEMPERATURE1"); // added by CB
			}
			if (viscosity_model == 5) // my(p,T), Reichenberg (1971)
			{
			}
			if (viscosity_model == 6) // my(C,T),
			{
			}
			if (viscosity_model == 7) // my(p,T,C)
			{
				in >> C_0;
				viscosity_pcs_name_vector.push_back("PRESSURE1");
				viscosity_pcs_name_vector.push_back("TEMPERATURE1");
			}
			if (viscosity_model == 9) // my(rho,T)
			{
				std::string arg1, arg2;
				in >> arg1 >> arg2; // get up to three arguments for density model

				if (arg1.length() == 0) // if no arguments are given use standard
				{
					arg1 = "PRESSURE1";
					arg2 = "TEMPERATURE1";
				}
				else if (arg2.length() == 0) // if only PRESSURE argument is given

					arg2 = "TEMPERATURE1";

				viscosity_pcs_name_vector.push_back(arg1);
				if (T_Process)
					viscosity_pcs_name_vector.push_back(arg2);
			}
			// AKS
			if (density_model == 15) // components constant viscosity
			{
				if (eos_name == "CONSTANT")
					in >> mu[0] >> mu[1] >> mu[2] >> mu[3];
			}

			if (viscosity_model == 18) // BG, NB calculated node viscosities from the phase transition model
			{
			}
			/// Viscosity model no. 19: Extracts viscosities from GEM3K calculations (not yet implemented)
			if (viscosity_model == 19) // KG44 extract viscosity from GEMS
			{
			}
			if (viscosity_model == 26) // && pcs_vector[0]->getProcessType() != FiniteElement::TNEQ)
			{
				std::cout << "Warning: This viscosity model requires two components and their molar masses defined in "
				             "the mcp file!\n";
				viscosity_pcs_name_vector.push_back("PRESSURE1");
				viscosity_pcs_name_vector.push_back("TEMPERATURE1");
				viscosity_pcs_name_vector.push_back("CONCENTRATION1");
			}
			if (viscosity_model == 30) // Reynolds model
			{
				in >> viscosity0 >> T_0 >> viscosity_T_star;
				viscosity_pcs_name_vector.push_back("PRESSURE1");
				viscosity_pcs_name_vector.push_back("TEMPERATURE1");
			}

			//    mfp_file->ignore(MAX_ZEILE,'\n');
			in.clear();
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$SPECIFIC_HEAT_CAPACITY") != string::npos)
		{
			in.str(GetLineFromFile1(mfp_file));
			in >> heat_capacity_model;
			// TF 11/2011 - used only in read- and write-method
			//			if(heat_capacity_model == 0) // c = fct(x)
			//				in >> heat_capacity_fct_name;
			if (heat_capacity_model == 1) // c = const
			{
				in >> specific_heat_capacity;
				specific_heat_capacity_pcs_name_vector.push_back("PRESSURE1");
				specific_heat_capacity_pcs_name_vector.push_back("TEMPERATURE1");
			}
			if (heat_capacity_model == 2) // my(p,T,C)
			{
				in >> C_0;
				specific_heat_capacity_pcs_name_vector.push_back("PRESSURE1");
				specific_heat_capacity_pcs_name_vector.push_back("TEMPERATURE1");
			}
			if (heat_capacity_model == 3) // YD: improved phase change
			{
				in >> T_Latent1; // Tmin for phase change
				T_Latent1 += TemperatureUnitOffset();
				in >> T_Latent2; // Tmax for phase change
				T_Latent2 += TemperatureUnitOffset();
				in >> heat_phase_change_curve;
				specific_heat_capacity_pcs_name_vector.push_back("PRESSURE1");
				specific_heat_capacity_pcs_name_vector.push_back("TEMPERATURE1");
				specific_heat_capacity_pcs_name_vector.push_back("SATURATION1");
				enthalpy_pcs_name_vector.push_back("TEMPERATURE1");
			}
			if (heat_capacity_model == 4) // YD: improved phase change, function
			{
				in >> T_Latent1; // Tmin for phase change
				T_Latent1 += TemperatureUnitOffset();
				in >> T_Latent2; // Tmax for phase change
				T_Latent1 += TemperatureUnitOffset();
				in >> specific_heat_capacity; // ^c
				in >> latent_heat; // L
				specific_heat_capacity_pcs_name_vector.push_back("PRESSURE1");
				specific_heat_capacity_pcs_name_vector.push_back("TEMPERATURE1");
				specific_heat_capacity_pcs_name_vector.push_back("SATURATION1");
				enthalpy_pcs_name_vector.push_back("TEMPERATURE1");
			}

			// AKS
			if (density_model == 15) // components constant density
			{
				if (eos_name == "CONSTANT")
					in >> cp[0] >> cp[1] >> cp[2] >> cp[3];
			}

			if (heat_capacity_model == 11
			    || heat_capacity_model == 12
			    || heat_capacity_model == 13) // && pcs_vector[0]->getProcessType() != FiniteElement::TNEQ)
			{
				std::cout << "Warning: This heat capacity model requires two components and their molar masses defined "
				             "in the mcp file!\n";
				specific_heat_capacity_pcs_name_vector.push_back("PRESSURE1");
				specific_heat_capacity_pcs_name_vector.push_back("TEMPERATURE1");
				specific_heat_capacity_pcs_name_vector.push_back("CONCENTRATION1");
			}

			in.clear();
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$HEAT_CONDUCTIVITY") != string::npos)
		{
			in.str(GetLineFromFile1(mfp_file));
			in >> heat_conductivity_model;
			// TF 11/2011 - used only in read- and write-method
			//			if(heat_conductivity_model == 0) // my = fct(x)
			//				in >> heat_conductivity_fct_name;
			if (heat_conductivity_model == 1) // my = const

				in >> heat_conductivity;
			if (heat_conductivity_model == 2) // my = f(p,T,C)
			{
				in >> C_0;
				heat_conductivity_pcs_name_vector.push_back("PRESSURE1");
				heat_conductivity_pcs_name_vector.push_back("TEMPERATURE1");
			}
			if (heat_conductivity_model == 3) // my = f(p,T) NB
			{
				std::string arg1, arg2;
				in >> arg1 >> arg2; // get up to three arguments for density model

				if (arg1.length() == 0) // if no arguments are given use standard
				{
					arg1 = "PRESSURE1";
					arg2 = "TEMPERATURE1";
				}
				else if (arg2.length() == 0) // if only PRESSURE argument is given

					arg2 = "TEMPERATURE1";

				heat_conductivity_pcs_name_vector.push_back(arg1);
				heat_conductivity_pcs_name_vector.push_back(arg2);
			}
			// AKS
			if (density_model == 15) // components constant density
			{
				if (eos_name == "CONSTANT")
					in >> kappa[0] >> kappa[1] >> kappa[2] >> kappa[3];
			}

			if (heat_conductivity_model == 11) // && pcs_vector[0]->getProcessType() != FiniteElement::TNEQ)
			{
				std::cout << "Warning: This heat conductivity model requires two components and their molar masses "
				             "defined in the mcp file!\n";
				heat_conductivity_pcs_name_vector.push_back("PRESSURE1");
				heat_conductivity_pcs_name_vector.push_back("TEMPERATURE1");
				heat_conductivity_pcs_name_vector.push_back("CONCENTRATION1");
			}

			in.clear(); // OK
			continue;
		}
		// subkeyword found
		if (line_string.find("$PHASE_DIFFUSION") != string::npos)
		{
			in.str(GetLineFromFile1(mfp_file));
			in >> diffusion_model;
			// TF 11/2011 - only read - not used
			//			if(diffusion_model == 0) // D = fct(x)
			//				in >> dif_fct_name;
			if (diffusion_model == 1) // D = const //MX
				in >> diffusion;
			if (diffusion_model == 2) // D = const //MX
			{
				in >> A_Daq >> B_Daq; // SP: Daq calculated after Yaws; fct(T)
				phase_diffusion_pcs_name_vector.push_back("TEMPERATURE1");
			}
			in.clear();
			continue;
		}
		//....................................................................
		if (line_string.find("$DIFFUSION") != string::npos)
		{
			in.str(GetLineFromFile1(mfp_file));
			in >> diffusion_model; // sub_line 1 for first phase
			if (diffusion_model == 15 && eos_name == "CONSTANT")
				in >> D0[0] >> D0[1] >> D0[2] >> D0[3];
			else
				in >> diffusion_coef; // sub_line 1
			in.clear();
		}
		//....................................................................
		if (line_string.find("$DECAY") != string::npos)
		{
			in.str(GetLineFromFile1(mfp_file));
			in >> decay_model >> lambda[0] >> lambda[1] >> lambda[2] >> lambda[3];
			in.clear();
		}
		//....................................................................
		if (line_string.find("$ISOTHERM") != string::npos)
		{
			in.str(GetLineFromFile1(mfp_file));
			in >> isotherm_model >> Kd[0] >> Kd[1] >> Kd[2] >> Kd[3];
			in.clear();
		}
		//....................................................................
		// CMCD outer space version
		if (line_string.find("$GRAVITY") != string::npos)
		{
			in.str(GetLineFromFile1(mfp_file));
			in >> gravity_constant;
			in.clear();
			continue;
		}
		if (line_string.find("$SPECIFIC_HEAT_SOURCE") != string::npos)
		{
			in.str(GetLineFromFile1(mfp_file));
			in >> specific_heat_source;
			in.clear();
			continue;
		}
	}
	return position;
}

/**************************************************************************
   FEMLib-Method:
   Task: Master read function
   Programing:
   08/2004 OK Implementation
   01/2005 OK Boolean type
   01/2005 OK Destruct before read
**************************************************************************/
bool MFPRead(std::string file_base_name)
{
	//----------------------------------------------------------------------
	// OK  MFPDelete();
	//----------------------------------------------------------------------
	CFluidProperties* m_mfp = NULL;
	char line[MAX_ZEILE];
	std::string sub_line;
	std::string line_string;
	std::ios::pos_type position;
	//========================================================================
	// File handling
	std::string mfp_file_name = file_base_name + MFP_FILE_EXTENSION;
	std::ifstream mfp_file(mfp_file_name.data(), std::ios::in);
	if (!mfp_file.good())
		return false;
	mfp_file.seekg(0L, std::ios::beg);
	//========================================================================
	// Keyword loop
	std::cout << "MFPRead"
	          << "\n";
	while (!mfp_file.eof())
	{
		mfp_file.getline(line, MAX_ZEILE);
		line_string = line;
		if (line_string.find("#STOP") != std::string::npos)
			return true;
		//----------------------------------------------------------------------
		// keyword found
		if (line_string.find("#FLUID_PROPERTIES") != std::string::npos)
		{
			m_mfp = new CFluidProperties();
			position = m_mfp->Read(&mfp_file);
			m_mfp->phase = (int)mfp_vector.size(); // OK4108
			mfp_vector.push_back(m_mfp);
			mfp_file.seekg(position, std::ios::beg);
		} // keyword found
	} // eof
	//========================================================================
	// Configuration
	int i;
	int no_fluids = (int)mfp_vector.size();
	if (no_fluids == 1)
	{
		m_mfp = mfp_vector[0];
		m_mfp->phase = 0;
	}
	else if (no_fluids == 2)
		for (i = 0; i < no_fluids; i++)
		{
			m_mfp = mfp_vector[i];
			if (m_mfp->name.find("GAS") != string::npos)
				m_mfp->phase = 0;
			else
				m_mfp->phase = 1;
		}
	//----------------------------------------------------------------------
	// Test
	if (mfp_vector.size() == 0)
	{
		std::cout << "Error in MFPRead: no MFP data"
		          << "\n";
		abort();
	}
	//----------------------------------------------------------------------
	return true;
}

/**************************************************************************
   FEMLib-Method:
   Task: write function
   Programing:
   11/2004 SB Implementation
   last modification:
**************************************************************************/
void CFluidProperties::Write(std::ofstream* mfp_file) const
{
	// KEYWORD
	*mfp_file << "#FLUID_PROPERTIES"
	          << "\n";
	*mfp_file << " $FLUID_TYPE"
	          << "\n";
	*mfp_file << "  " << name << "\n";
	*mfp_file << " $DAT_TYPE"
	          << "\n";
	*mfp_file << "  " << name << "\n";
	*mfp_file << " $DENSITY"
	          << "\n";
	// TF 11/2011 - _rho_fct_name is used only here and in the read-method
	//	if(density_model == 0)
	//		*mfp_file << "  " << density_model << " " << _rho_fct_name << "\n";
	if (density_model == 1)
		*mfp_file << "  " << density_model << " " << rho_0 << "\n";
	// todo
	*mfp_file << " $VISCOSITY"
	          << "\n";
	// TF 11/2011 - used only in read- and write-method
	//	if(viscosity_model == 0)
	//		*mfp_file << "  " << viscosity_model << " " << _my_fct_name << "\n";
	if (viscosity_model == 1)
		*mfp_file << "  " << viscosity_model << " " << my_0 << "\n";
	*mfp_file << " $PHASE_DIFFUSION"
	          << "\n";
	// if(diffusion_model == 0) *mfp_file << "  " << diffusion_model << " " << dif_fct_name << "\n";
	if (diffusion_model == 1)
		*mfp_file << "  " << diffusion_model << " " << diffusion << "\n";
	if (diffusion_model == 2)
		*mfp_file << "  " << diffusion_model << " "
		          << "\n";
	// todo
	*mfp_file << " $SPECIFIC_HEAT_CAPACITY"
	          << "\n";
	// TF 11/2011 - used only in read- and write-method
	//	if(heat_capacity_model == 0)
	//		*mfp_file << "  " << heat_capacity_model << " " << heat_capacity_fct_name <<
	//		"\n";
	if (heat_capacity_model == 1)
		*mfp_file << "  " << heat_capacity_model << " " << specific_heat_capacity << "\n";
	*mfp_file << " $SPECIFIC_HEAT_CONDUCTIVITY"
	          << "\n";
	// TF 11/2011 - used only in read- and write-method
	//	if(heat_conductivity_model == 0)
	//		*mfp_file << "  " << heat_conductivity_model << " " <<
	//		heat_conductivity_fct_name << "\n";
	if (heat_conductivity_model == 1)
		*mfp_file << "  " << heat_conductivity_model << " " << heat_conductivity << "\n";
	//--------------------------------------------------------------------
}

/**************************************************************************
   FEMLib-Method:
   Task: Master write function
   Programing:
   08/2004 OK Implementation
   last modification:
**************************************************************************/
void MFPWrite(std::string base_file_name)
{
	CFluidProperties* m_mfp = NULL;
	string sub_line;
	string line_string;
	ofstream mfp_file;
	//========================================================================
	// File handling
	std::string mfp_file_name = base_file_name + MFP_FILE_EXTENSION;
	mfp_file.open(mfp_file_name.data(), std::ios::trunc | std::ios::out);
	mfp_file.setf(std::ios::scientific, std::ios::floatfield);
	mfp_file.precision(12);
	if (!mfp_file.good())
		return;
	//  mfp_file.seekg(0L,ios::beg);
	//========================================================================
	mfp_file << "GeoSys-MFP: Material Fluid Properties -------------"
	         << "\n";
	//========================================================================
	// OUT vector
	int mfp_vector_size = (int)mfp_vector.size();
	int i;
	for (i = 0; i < mfp_vector_size; i++)
	{
		m_mfp = mfp_vector[i];
		m_mfp->Write(&mfp_file);
	}
	mfp_file << "#STOP";
	mfp_file.close();
	//  delete mfp_file;
}

////////////////////////////////////////////////////////////////////////////
// Properties functions
////////////////////////////////////////////////////////////////////////////
/**************************************************************************
   FEMLib-Method:
   Task: Master calc function
   Programing:
   09/2005 WW implementation
   11/2005 YD modification
   11/2005 CMCD Inclusion current and previous time step quantities
   05/2007 PCH improvement for density-dependent flow
   last modification:
**************************************************************************/
void CFluidProperties::CalPrimaryVariable(std::vector<std::string>& pcs_name_vector)
{
	CRFProcess* m_pcs = NULL;

	int nidx0, nidx1;
	if (!Fem_Ele_Std) // OK
		return;

	primary_variable[0] = 0;
	primary_variable[1] = 0.;
	primary_variable[2] = 0;

	for (int i = 0; i < (int)pcs_name_vector.size(); i++)
	{
		// MX  m_pcs = PCSGet("HEAT_TRANSPORT");
		m_pcs = PCSGet(pcs_name_vector[i], true);
		if (!m_pcs)
			return; // MX
		nidx0 = m_pcs->GetNodeValueIndex(pcs_name_vector[i]);
		nidx1 = nidx0 + 1;

		if (mode == 0) // Gauss point values
		{
			primary_variable_t0[i] = Fem_Ele_Std->interpolate(nidx0, m_pcs);
			primary_variable_t1[i] = Fem_Ele_Std->interpolate(nidx1, m_pcs);
			primary_variable[i] = (1. - Fem_Ele_Std->pcs->m_num->ls_theta) * Fem_Ele_Std->interpolate(nidx0, m_pcs)
			                      + Fem_Ele_Std->pcs->m_num->ls_theta * Fem_Ele_Std->interpolate(nidx1, m_pcs);
		}
		else if (mode == 2) // Element average value
		{
			primary_variable[i] = (1. - Fem_Ele_Std->pcs->m_num->ls_theta) * Fem_Ele_Std->elemnt_average(nidx0, m_pcs)
			                      + Fem_Ele_Std->pcs->m_num->ls_theta * Fem_Ele_Std->elemnt_average(nidx1, m_pcs);
			primary_variable_t0[i] = Fem_Ele_Std->elemnt_average(nidx0, m_pcs);
			primary_variable_t1[i] = Fem_Ele_Std->elemnt_average(nidx1, m_pcs);
		}
		else if (mode == 3) // NB, just testing
		{
			primary_variable[i] = Fem_Ele_Std->interpolate(nidx0, m_pcs);
		}
		else
		{
			if (pcs_name_vector[i].compare("TEMPERATURE1") == 0)
			{
				primary_variable[i] = _reference_temperature;
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////
// Fluid density
/**************************************************************************
   FEMLib-Method:
   Task: Master calc function
   Programing:
   08/2004 OK MFP implementation
           based on MATCalcFluidDensity by OK/JdJ,AH,MB
   11/2005 YD Modification
   05/2008 WW Add an argument: double* variables: P, T, C
   last modification:
   NB 4.9.05
**************************************************************************/
double CFluidProperties::Density(double* variables)
{
	double pressure;
	double Rho = 0.0;
	static double density;
	// static double air_gas_density,vapour_density,vapour_pressure;

	int gueltig;

	//----------------------------------------------------------------------
	if (variables) // This condition is added by WW
	{
		//----------------------------------------------------------------------
		// Duplicate the following lines just to enhance computation. WW
		switch (density_model)
		{
			case 0: // rho = f(x)
				density = GetCurveValue(density_curve_number, 0, variables[1], &gueltig);
				break;
			case 1: // rho = const
				density = rho_0;
				break;
			case 2: // rho(p) = rho_0*(1+beta_p*(p-p_0))
				density = rho_0 * (1. + drho_dp * (max(variables[0], 0.0) - p_0));
				break;
			case 3: // rho(C) = rho_0*(1+beta_C*(C-C_0))
				density = rho_0 * (1. + drho_dC * (max(variables[2], 0.0) - C_0));
				break;
			case 4: // rho(T) = rho_0*(1+beta_T*(T-T_0))
				density = rho_0 * (1. + drho_dT * (max(variables[1], 0.0) - T_0));
				break;
			case 5: // rho(C,T) = rho_0*(1+beta_C*(C-C_0)+beta_T*(T-T_0))
				density = rho_0
				          * (1. + drho_dC * (max(variables[2], 0.0) - C_0) + drho_dT * (max(variables[1], 0.0) - T_0));
				break;
			case 6: // rho(p,T) = rho_0*(1+beta_p*(p-p_0)+beta_T*(T-T_0))
				density = rho_0
				          * (1. + drho_dp * (max(variables[0], 0.0) - p_0) + drho_dT * (max(variables[1], 0.0) - T_0));
				break;
			case 7: // Pefect gas. WW
				density = variables[0] * molar_mass / (PhysicalConstant::IdealGasConstant * variables[1]);
				break;
			case 8: // M14 von JdJ // 25.1.12 Added by CB for density output AB-model
				density = MATCalcFluidDensityMethod8(variables[0], variables[1], variables[2]);
				break;
			case 10: // Get density from temperature-pressure values from fct-file	NB 4.8.01
				density = GetMatrixValue(variables[1], variables[0], fluid_name, &gueltig);
				break;
			case 11: // Redlich-Kwong EOS for different fluids NB 4.9.05
				density = rkeos(variables[1], variables[0], fluid_id);
				break;
			case 12: // Peng-Robinson EOS for different fluids NB 4.9.05
				// NB
				density = preos(this, variables[1], variables[0]);
				break;
			case 13: // Helmholtz free Energy NB JUN 09
				// NB
				density = zero(variables[1], variables[0], fluid_id, 1e-8);
				break;
			case 14: // #Exponential law#
				density = rho_0 * exp(drho_dp * (max(variables[0], 0.0) - p_0) + drho_dT * (max(variables[2], 0.0))
				                      + drho_dC * max(variables[2], 0.0));
				break;

			case 15: // mixture 1/rho= sum_i x_i/rho_i #p, T, x:->Amagat's law#
				for (int CIndex = 2; CIndex < cmpN + 2; CIndex++)
					Rho += variables[CIndex] / ComponentDensity(CIndex, variables);
				density = 1 / Rho;
				break;

			case 18: // using calculated densities at nodes from the phase transition model, BG, NB 11/2010
				variables[2] = phase;
				density = GetElementValueFromNodes(int(variables[0]),
				                                   int(variables[1]),
				                                   int(variables[2]),
				                                   0); // hand over element index, Gauss point index and phase index
				break;
			case 19: // KG44 get the density from GEMS calculations
				// seems complicated, as we probably have to call GEMS.....or take values from last GEMS calculation
				// ---> update during iterations is not possible
				//  long elem = Fem_Ele_Std->GetMeshElement()->GetIndex();
				density = 1000.0;
#ifdef GEM_REACT
				if (!Fem_Ele_Std) // for Richards flow (when saturation is needed initially in GEMS setup) we have to
				// make
				{
					density = 1000.0;
				}
				else
				{
					// elem = Fem_Ele_Std->GetMeshElement()->GetIndex(); //kg44 need element index or node index for
					// GEMS ...currently we stick to arithmetric average of node data (even for gauss points!)
					density = m_vec_GEM->REACT_GEM::FluidDensity(
					    long(variables[0]), int(variables[1])); // hand over element index and gauss point index
					// here we can interpolate values from nodes to elemnt or to gauss points
				}
#endif
				// insert call for GEMS densities..
				break;
			case 20: // rho(p,T, C) for water, range p < 100 MPa, 0 <= T <= 350 °C   Magri GFZ thesis

				pressure = variables[0] / 1e5;
				density = 9.99792877961606e+02 + 5.07605113140940e-04 * max(pressure, 0.0)
				          - 5.28425478164183e-10 * pow(max(pressure, 0.0), 2.)
				          + (5.13864847162196e-02 - 3.61991396354483e-06 * max(pressure, 0.0)
				             + 7.97204102509724e-12 * pow(max(pressure, 0.0), 2.))
				                * max(variables[1], 0.0)
				          + (-7.53557031774437e-03 + 6.32712093275576e-08 * max(pressure, 0.0)
				             - 1.66203631393248e-13 * pow(max(pressure, 0.0), 2.))
				                * pow(max(variables[1], 0.0), 2.)
				          + (4.60380647957350e-05 - 5.61299059722121e-10 * max(pressure, 0.0)
				             + 1.80924436489400e-15 * pow(max(pressure, 0.0), 2.))
				                * pow(max(variables[1], 0.0), 3.)
				          + (-2.26651454175013e-07 + 3.36874416675978e-12 * max(pressure, 0.0)
				             - 1.30352149261326e-17 * pow(max(pressure, 0.0), 2.))
				                * pow(max(variables[1], 0.0), 4.)
				          + (6.14889851856743e-10 - 1.06165223196756e-14 * max(pressure, 0.0)
				             + 4.75014903737416e-20 * pow(max(pressure, 0.0), 2.))
				                * pow(max(variables[1], 0.0), 5.)
				          + (-7.39221950969522e-13 + 1.42790422913922e-17 * max(pressure, 0.0)
				             - 7.13130230531541e-23 * pow(max(pressure, 0.0), 2.))
				                * pow(max(variables[1], 0.0), 6.);

				if (fabs(drho_dC) > 1.e-20)
					density *= 1. + drho_dC * (max(variables[2], 0.0) - C_0);

				break;
			case 26: // Dalton's law + ideal gas for use with TNEQ/TES
			{
				const double M0 = cp_vec[0]->molar_mass; // molar mass of component 0
				const double M1 = cp_vec[1]->molar_mass;
				const double p = variables[0];
				const double T = variables[1];
				const double x = variables[2]; // gas mass fraction of component 1
				// assert(0.0 <= x && x <= 1.0);

				// gas molar fraction of component 1
				const double xn = M0 * x / (M0 * x + M1 * (1.0 - x));

				density = p / (PhysicalConstant::IdealGasConstant * T) * (M1 * xn + M0 * (1.0 - xn)); // R_uni in mNs
			}
			break;

			default:
				std::cout << "Error in CFluidProperties::Density: no valid model"
				          << "\n";
				break;
		}
	}
	else
	{
		CalPrimaryVariable(density_pcs_name_vector);

		//----------------------------------------------------------------------
		switch (density_model)
		{
			case 0: // rho = f(x)
				density = GetCurveValue(density_curve_number, 0, primary_variable[0], &gueltig);
				break;
			case 1: // rho = const
				density = rho_0;
				break;
			case 2: // rho(p) = rho_0*(1+beta_p*(p-p_0))
				density = rho_0 * (1. + drho_dp * (max(primary_variable[0], 0.0) - p_0));
				break;
			case 3: // rho(C) = rho_0*(1+beta_C*(C-C_0))
				density = rho_0 * (1. + drho_dC * (max(primary_variable[0], 0.0) - C_0));
				break;
			case 4: // rho(T) = rho_0*(1+beta_T*(T-T_0))
				density = rho_0 * (1. + drho_dT * (max(primary_variable[0], 0.0) - T_0));
				break;
			case 5: // rho(C,T) = rho_0*(1+beta_C*(C-C_0)+beta_T*(T-T_0))
				density = rho_0 * (1. + drho_dC * (max(primary_variable[0], 0.0) - C_0)
				                   + drho_dT * (max(primary_variable[1], 0.0) - T_0));
				break;
			case 6: // rho(p,T) = rho_0*(1+beta_p*(p-p_0)+beta_T*(T-T_0))
				density = rho_0 * (1. + drho_dp * (max(primary_variable[0], 0.0) - p_0)
				                   + drho_dT * (max(primary_variable[1], 0.0) - T_0));
				break;
			case 7: // rho_w^l(p,T) for gas phase
				/* //WW
				{
				    const double vapour_pressure = MFPCalcVapourPressure(primary_variable[0]);
				    air_gas_density = (MolarMass::Air * (primary_variable[1]-vapour_pressure)) /
				(PhysicalConstant::IdealGasConstant*(primary_variable[0]+0.0));
				    vapour_density = (MolarMass::Water * vapour_pressure) /
				(PhysicalConstant::IdealGasConstant*(primary_variable[0]+0.0));
				    density = vapour_density + air_gas_density;
				}
				*/
				break;
			case 8: // M14 von JdJ
				density = MATCalcFluidDensityMethod8(primary_variable[0], primary_variable[1], primary_variable[2]);
				break;
			case 10: // Get density from temperature-pressure values from fct-file NB
				density = GetMatrixValue(primary_variable[1], primary_variable[0], fluid_name, &gueltig);
				break;
			case 11: // Peng-Robinson equation of state NB
				density = rkeos(primary_variable[1], primary_variable[0], fluid_id);
				break;
			case 12: // Redlich-Kwong equation of state NB
				density = preos(this, primary_variable[1], primary_variable[0]);
				break;
			case 13: // Helmholtz free Energy NB JUN 09
				// NB
				density = zero(primary_variable[1], primary_variable[0], fluid_id, 1e-8);
				break;

			case 14: // #Exponential law#
				density = rho_0 * exp(drho_dp * (max(variables[0], 0.0) - p_0) + drho_dT * (max(variables[2], 0.0))
				                      + drho_dC * max(variables[2], 0.0));
				break;
			case 15: // mixture 1/rho= sum_i x_i/rho_i #p, T, x:-> Amagat's law#
				for (int CIndex = 2; CIndex < cmpN + 2; CIndex++)
					Rho += variables[CIndex] / ComponentDensity(CIndex, variables);
				density = 1 / Rho;
				break;
			case 18: // using calculated densities at nodes from the phase transition model, BG, NB 11/2010
				// Just dummy function, so density is not 0 ;
				std::cout << " Error - Density Model 18 not implemented, usind dummy density of 1000."
				          << "\n";
				density = 1000; // Achtung - dummy
                                break;

			case 19: // KG44 get the density from GEMS calculations
				// seems complicated, as we probably have to call GEMS.....or take values from last GEMS calculation
				// ---> update during iterations is not possible
				//  long elem = Fem_Ele_Std->GetMeshElement()->GetIndex();
				density = 1000.0;
#ifdef GEM_REACT
				if (!Fem_Ele_Std) // for Richards flow (when saturation is needed initially in GEMS setup) we have to
				// make
				{
					density = 1000.0;
				}
				else
				{
					long elem = Fem_Ele_Std->GetMeshElement()->GetIndex(); // kg44 need element index or node index for
					// GEMS ...currently we stick to arithmetric
					// average of node data (even for gauss
					// points!)
					density = m_vec_GEM->REACT_GEM::FluidDensity(
					    elem, -1); // hand over element index and set gauss point to -1...
					// Remark: Interpolation from Gauss points to element is not possible here: Fluid densities are only
					// calculated at nodes, therefore interpolation from gauss points onto element will not work
				}
#endif
				// insert call for GEMS densities..
				break;
			case 20: // rho(p,T, C) for water, range p < 100 MPa, 0 <= T <= 350 °C   Magri GFZ thesis
				pressure = primary_variable[0] / 1e5;

				density = 9.99792877961606e+02 + 5.07605113140940e-04 * max(pressure, 0.0)
					- 5.28425478164183e-10 * pow(max(pressure, 0.0), 2.)
					+ (5.13864847162196e-02 - 3.61991396354483e-06 * max(pressure, 0.0)
							+ 7.97204102509724e-12 * pow(max(pressure, 0.0), 2.))
					* max(primary_variable[1], 0.0)
					+ (-7.53557031774437e-03 + 6.32712093275576e-08 * max(pressure, 0.0)
							- 1.66203631393248e-13 * pow(max(pressure, 0.0), 2.))
					* pow(max(primary_variable[1], 0.0), 2.)
					+ (4.60380647957350e-05 - 5.61299059722121e-10 * max(pressure, 0.0)
							+ 1.80924436489400e-15 * pow(max(pressure, 0.0), 2.))
					* pow(max(primary_variable[1], 0.0), 3.)
					+ (-2.26651454175013e-07 + 3.36874416675978e-12 * max(pressure, 0.0)
							- 1.30352149261326e-17 * pow(max(pressure, 0.0), 2.))
					* pow(max(primary_variable[1], 0.0), 4.)
					+ (6.14889851856743e-10 - 1.06165223196756e-14 * max(pressure, 0.0)
							+ 4.75014903737416e-20 * pow(max(pressure, 0.0), 2.))
					* pow(max(primary_variable[1], 0.0), 5.)
					+ (-7.39221950969522e-13 + 1.42790422913922e-17 * max(pressure, 0.0)
							- 7.13130230531541e-23 * pow(max(pressure, 0.0), 2.))
					* pow(max(primary_variable[1], 0.0), 6.);

				if (fabs(drho_dC) > 1.e-20)
					density *= 1. + drho_dC * (max(primary_variable[2], 0.0) - C_0);
				break;
			case 26: // Dalton's law + ideal gas for use with TNEQ/TES
			{
				const double M0 = cp_vec[0]->molar_mass; // molar mass of component 0
				const double M1 = cp_vec[1]->molar_mass;
				const double p = primary_variable[0];
				const double T = primary_variable[1];
				const double x = primary_variable[2]; // gas mass fraction of component 1
				// assert(0.0 <= x && x <= 1.0);

				// gas molar fraction of component 1
				const double xn = M0 * x / (M0 * x + M1 * (1.0 - x));

				density = p / (PhysicalConstant::IdealGasConstant * T) * (M1 * xn + M0 * (1.0 - xn)); // R_uni in mNs
				break;
			}
			default:
				std::cout << "Error in CFluidProperties::Density : no valid model (variables definded?)"
					<< "\n";
				break;
		}
	}
	return density;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: GetElementValueFromNodes
   Task: Interpolates node values like density or viscosity to the elements (if GPIndex < 0) or to element gauss points
   Return: interpolated variable
   Programming: 11/2010 BG
   Modification:
   -------------------------------------------------------------------------*/
double CFluidProperties::GetElementValueFromNodes(long ElementIndex, int GPIndex, int PhaseIndex, int VariableIndex)
{
	CFEMesh* m_msh = fem_msh_vector[0]; // SB: ToDo hart gesetzt
	MeshLib::CElem* m_ele = NULL;
	MeshLib::CNode* m_node = NULL;
	CRFProcess* m_pcs;
	double var, variable;
	int variable_index = 0, nNodes;
	double distance, weight, sum_weights;
	Math_Group::vec<long> vec_nod_index(8);

	variable = 0;

	m_pcs = PCSGet("MULTI_PHASE_FLOW");
	// check if PHASE_TRANSITION is used for the process
	if ((density_model == 18) || (viscosity_model == 18))
		if (m_pcs->Phase_Transition_Model != 1)
		{
			cout << "The Phase_Transition_Model should be used together with the density and viscosity model 18 !"
			     << "\n";
			cout << "The run is terminated now ..."
			     << "\n";
			//system("Pause");
			exit(0);
		}

	m_ele = m_msh->ele_vector[ElementIndex]; // get element
	// if GPIndex > 0 -> interpolation to GP if not then interpolation to the Element centre
	if (PhaseIndex == 0)
	{
		switch (VariableIndex)
		{
			case 0: // Density
				variable_index = m_pcs->GetNodeValueIndex("DENSITY1");
				break;
			case 1: // Viscosity
				variable_index = m_pcs->GetNodeValueIndex("VISCOSITY1");
				break;
		}
	}
	else
	{
		switch (VariableIndex)
		{
			case 0: // Density
				variable_index = m_pcs->GetNodeValueIndex("DENSITY2");
				break;
			case 1: // Viscosity
				variable_index = m_pcs->GetNodeValueIndex("VISCOSITY2");
				break;
		}
	}

	if (GPIndex > -1)
	{
		// interpolate density to gauss point
		if (m_ele->GetMark()) // Marked for use
			// Configure Element for interpolation of node velocities to GP velocities
			// Fem_Ele_Std->ConfigElement(m_ele);
			variable = Fem_Ele_Std->InterpolatePropertyToGausspoint(GPIndex, m_pcs, variable_index);
	}
	else
	{
		distance = weight = sum_weights = 0.0;
		if (m_ele->GetElementType() == 0)
			nNodes = m_ele->GetNodesNumber(true);
		else
			nNodes = m_ele->GetNodesNumber(false);
		m_ele->GetNodeIndeces(vec_nod_index);

		for (long i = 0; i < int(nNodes); i++) // go through list of connected nodes
		{ // Get the connected node
			m_node = m_msh->nod_vector[vec_nod_index[i]];
			// calculate distance between the node and the barycentre
			double const* gravity_centre(m_ele->GetGravityCenter());
			double const* const pnt(m_node->getData());
			distance = (gravity_centre[0] - pnt[0]) * (gravity_centre[0] - pnt[0]);
			distance += (gravity_centre[1] - pnt[1]) * (gravity_centre[1] - pnt[1]);
			distance += (gravity_centre[2] - pnt[2]) * (gravity_centre[2] - pnt[2]);
			distance = sqrt(distance);

			// Weight of each face depending on distance
			weight = (1.0 / distance);
			// Sum of weights
			sum_weights += weight;
			// Density
			var = m_pcs->GetNodeValue(int(m_node->GetIndex()), variable_index);
			variable += var * weight;
		}
		variable = variable / sum_weights;
	}
	// cout << "Variable: " << variable << "\n";
	return variable;
}

/*************************************************************************
   ROCKFLOW - Funktion: MATCalcFluidDensityMethod8

   Task:
   Density of a geothermal fluid as a function of temperature, pressure
   and concentration.

   Fluid-Density function according to IAPWS-IF97

   Programmaenderungen:
   09/2003   CMCD  ARL  First implementation
   09/2004   CMCD  Inclusion in GeoSys vs. 4

*************************************************************************/
double CFluidProperties::MATCalcFluidDensityMethod8(double Press, double TempK, double Conc)
{
	Conc = Conc;
	/*int c_idx;*/
	double rho_0;
	double GammaPi, Pressurevar, Tau, pressure_average, temperature_average;
	double Tstar, Pstar, GazConst;
	double L[35], J[35], n[35];
	int i;
	double salinity;

	pressure_average = Press;
	temperature_average = TempK;
	salinity = C_0;
	Tstar = 1386;
	Pstar = 16.53e6; // MPa
	GazConst = 0.461526e3; //

	n[0] = 0.0;
	n[1] = 0.14632971213167;
	n[2] = -0.84548187169114;
	n[3] = -0.37563603672040e1;
	n[4] = 0.33855169168385e1;
	n[5] = -0.95791963387872;
	n[6] = 0.15772038513228;
	n[7] = -0.16616417199501e-1;
	n[8] = 0.81214629983568e-3;
	n[9] = 0.28319080123804e-3;
	n[10] = -0.60706301565874e-3;
	n[11] = -0.18990068218419e-1;
	n[12] = -0.32529748770505e-1;
	n[13] = -0.21841717175414e-1;
	n[14] = -0.52838357969930e-4;
	n[15] = -0.47184321073267e-3;
	n[16] = -0.30001780793026e-3;
	n[17] = 0.47661393906987e-4;
	n[18] = -0.44141845330846e-5;
	n[19] = -0.72694996297594e-15;
	n[20] = -0.31679644845054e-4;
	n[21] = -0.28270797985312e-5;
	n[22] = -0.85205128120103e-9;
	n[23] = -0.22425281908000e-5;
	n[24] = -0.65171222895601e-6;
	n[25] = -0.14341729937924e-12;
	n[26] = -0.40516996860117e-6;
	n[27] = -0.12734301741641e-8;
	n[28] = -0.17424871230634e-9;
	n[29] = -0.68762131295531e-18;
	n[30] = 0.14478307828521e-19;
	n[31] = 0.26335781662795e-22;
	n[32] = -0.11947622640071e-22;
	n[33] = 0.18228094581404e-23;
	n[34] = -0.93537087292458e-25;

	L[0] = 0.;
	L[1] = 0.;
	L[2] = 0.;
	L[3] = 0.;
	L[4] = 0.;
	L[5] = 0.;
	L[6] = 0.;
	L[7] = 0.;
	L[8] = 0.;
	L[9] = 1.;
	L[10] = 1.;
	L[11] = 1.;
	L[12] = 1.;
	L[13] = 1.;
	L[14] = 1.;
	L[15] = 2.;
	L[16] = 2.;
	L[17] = 2.;
	L[18] = 2.;
	L[19] = 2.;
	L[20] = 3.;
	L[21] = 3.;
	L[22] = 3.;
	L[23] = 4.;
	L[24] = 4.;
	L[25] = 4.;
	L[26] = 5.;
	L[27] = 8.;
	L[28] = 8.;
	L[29] = 21.;
	L[30] = 23.;
	L[31] = 29.;
	L[32] = 30.;
	L[33] = 31.;
	L[34] = 32.;

	J[0] = -2.;
	J[1] = -2.;
	J[2] = -1.;
	J[3] = 0.;
	J[4] = 1.;
	J[5] = 2.;
	J[6] = 3.;
	J[7] = 4.;
	J[8] = 5.;
	J[9] = -9.;
	J[10] = -7.;
	J[11] = -1.;
	J[12] = 0.;
	J[13] = 1.;
	J[14] = 3.;
	J[15] = -3.;
	J[16] = 0.;
	J[17] = 1.;
	J[18] = 3.;
	J[19] = 17.;
	J[20] = -4.;
	J[21] = 0.;
	J[22] = 6.;
	J[23] = -5.;
	J[24] = -2.;
	J[25] = 10.;
	J[26] = -8.;
	J[27] = -11.;
	J[28] = -6.;
	J[29] = -29.;
	J[30] = -31.;
	J[31] = -38.;
	J[32] = -39.;
	J[33] = -40.;
	J[34] = -41.;

	Pressurevar = pressure_average / Pstar;
	Tau = Tstar / temperature_average;

	/*BEGIN:Calculation of GammaPi*/
	GammaPi = 0.;

	for (i = 1; i < 35; i++)

		GammaPi = GammaPi - (n[i]) * (L[i]) * (pow((7.1 - Pressurevar), (L[i] - 1))) * (pow((Tau - 1.222), J[i]));

	/*END: Calculation of GammaPi*/

	/*BEGIN: Calculation of density*/
	rho_0 = pressure_average / (GazConst * temperature_average * Pressurevar * GammaPi);
	/*END: Calculation of density*/

	/*  return rho_0 + drho_dC * (concentration_average - c0);   */
	/*printf("%f", rho_0 + salinity);*/
	return rho_0 + salinity;
}

////////////////////////////////////////////////////////////////////////////
// Fluid viscosity
/**************************************************************************
   FEMLib-Method:
   Task: Master calc function
   Programing:
   08/2004 OK Implementation
   11/2005 YD Modification
   10/2008 OK Faster data access
   03/2009 NB Viscosity depending on Density()
**************************************************************************/
// OK4709
double CFluidProperties::Viscosity(double* variables)
{
	double viscosity = 0.0;

	//----------------------------------------------------------------------
	if (variables) // OK4709: faster data access
	{
		primary_variable[0] = variables[0]; // p (single phase)
		primary_variable[1] = variables[1]; // T (temperature)
		primary_variable[2] = variables[2]; // C (salinity)//index in case of AIR_FLOW
	}
	else
		CalPrimaryVariable(viscosity_pcs_name_vector);
	//----------------------------------------------------------------------
	switch (viscosity_model)
	{
		case 0: // rho = f(x)
		{
			int gueltig;
			viscosity = GetCurveValue(viscosity_curve_number, 0, primary_variable[1], &gueltig);
			break;
		}
		case 1: // my = const
			viscosity = my_0;
			break;
		case 2: // my(p) = my_0*(1+gamma_p*(p-p_0))
			viscosity = my_0 * (1. + dmy_dp * (max(primary_variable[0], 0.0) - p_0));
			break;
		case 3: // my^l(T), Yaws et al. (1976)
		{
			CRFProcess* m_pcs = PCSGet("MULTI_COMPONENTIAL_FLOW");

			if (mode == 1) // OK4704 for nodal output
			{
				m_pcs = PCSGet("HEAT_TRANSPORT");
				// if(!m_pcs) return 0.0;
				primary_variable[1] = m_pcs->GetNodeValue(node, m_pcs->GetNodeValueIndex("TEMPERATURE1") + 1);
			}
			// ToDo pcs_name
			// used within this model
			viscosity = LiquidViscosity_Yaws_1976(primary_variable[1]);
			break;
		}
		case 4: // my^g(T), Marsily (1986)
			viscosity = LiquidViscosity_Marsily_1986(primary_variable[1]);
			// viscosity = LiquidViscosity_Marsily_1986(primary_variable[0]);
			break;
		case 5: // my^g(p,T), Reichenberg (1971)
			viscosity = GasViscosity_Reichenberg_1971(primary_variable[0], primary_variable[1]);
			break;
		case 6: // my(C,T),
			viscosity = LiquidViscosity_NN(primary_variable[0], primary_variable[1]);
			break;
		case 8: // my(p,C,T),
			viscosity = LiquidViscosity_CMCD(primary_variable[0], primary_variable[1], primary_variable[2]);
			break;
		case 9: // viscosity as function of density and temperature, NB
		{
			double mfp_arguments[2];
			// TODO: switch case different models...
			// Problem.cpp 3 PS_GLOBAL, 1212,1313 pcs_type
			// TODO: default fluid_ID, if not specified
			mfp_arguments[0] = primary_variable[0]; // rescue primary_variable before its destroyed by Density();
			mfp_arguments[1] = primary_variable[1];

			const double density
			    = Density(mfp_arguments); // TODO: (NB) store density (and viscosity) as secondary variable
			// NB
			viscosity = Fluid_Viscosity(density, mfp_arguments[1], mfp_arguments[0], fluid_id);
			break;
		}
		case 15: // mixture 1/�= sum_i y_i/�_i:: VTPR-EoS
		{
			CRFProcess* m_pcs = PCSGet("MULTI_COMPONENTIAL_FLOW");
			double my = 0.0;

			for (int CIndex = 2; CIndex < cmpN + 2; CIndex++)
			{
				if (eos_name == "CONSTANT")
				{
					my += variables[CIndex] / mu[CIndex - 2];
				}
				else
				{
					therm_prop(m_pcs->pcs_primary_function_name[CIndex]);
					my += variables[CIndex]
					      / Fluid_Viscosity(ComponentDensity(CIndex, variables), variables[1], variables[0], fluid_id);
				}
			}
			viscosity = 1.0 / my;
			break;
		}

		case 18: // BG, NB using calculated viscosities at nodes from the phase transition model
			variables[2] = phase;
			viscosity = GetElementValueFromNodes(
			    int(variables[0]), int(variables[1]), int(variables[2]),
			    1); // hand over element index, Gauss point index, phase index and variable index
			break;
		case 19: // reserved for GEMS coupling
			break;
		case 26: // Wilke (see Poling, B. E.; Prausnitz, J. M.; John Paul, O. & Reid, R. C. The properties of gases and
			// liquids McGraw-Hill New York, 2001, 5: page 9.21)
			{
				double x[2], Vs[2];

				const double M1 = cp_vec[0]->molar_mass;
				const double M2 = cp_vec[1]->molar_mass;

				const double p = primary_variable[0];
				const double T = primary_variable[1];
				const double X = primary_variable[2];

				// reactive component
				x[0] = M1 * X / (M1 * X + M2 * (1.0 - X)); // mass in mole fraction
				Vs[0] = Fluid_Viscosity(M2 * p / (PhysicalConstant::IdealGasConstant * T), T, p, cp_vec[1]->fluid_id);
				// inert component
				x[1] = 1.0 - x[0];
				Vs[1] = Fluid_Viscosity(M1 * p / (PhysicalConstant::IdealGasConstant * T), T, p,
				                        cp_vec[0]->fluid_id); // R_uni in mNs

				const double M1_over_M2(M2 / M1); // reactive over inert
				const double V1_over_V2(Vs[0] / Vs[1]);

				const double phi_12 = (1.0 + pow(V1_over_V2, 0.5) * pow(1.0 / M1_over_M2, 0.25))
				                      * (1.0 + pow(V1_over_V2, 0.5) * pow(1.0 / M1_over_M2, 0.25))
				                      / pow(8.0 * (1.0 + M1_over_M2), 0.5);
				const double phi_21 = phi_12 * M1_over_M2 / V1_over_V2;

				viscosity = Vs[0] * x[0] / (x[0] + x[1] * phi_12) + Vs[1] * x[1] / (x[1] + x[0] * phi_21);

				break;
			}
		case 30: // exp(T) e.g Reynolds
			viscosity = LiquidViscosity_expo(primary_variable[1]);
			break;
		default:
			cout << "Error in CFluidProperties::Viscosity: no valid model"
			     << "\n";
			break;
	}
	//----------------------------------------------------------------------

	return viscosity;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Dynamische Gasviskositaet nach Reichenberg (1971)
   als Funktion von Temperatur und Druck
   in Reid et al. (1988), S. 420
   Programing:
   08/2004 OK MFP implementation based on CalcFluidViscosityMethod7 by OK
   last modification:
**************************************************************************/
double CFluidProperties::GasViscosity_Reichenberg_1971(double p, double T)
{
	double my, my0;
	double A, B, C, D;
	double Q, Pr, Tr;
	double alpha1, alpha2, beta1, beta2, gamma1, gamma2, delta1, delta2;
	double a, c, d;
	double pc, Tc;

	alpha1 = 1.9824e-3;
	alpha2 = 5.2683;
	beta1 = 1.6552;
	beta2 = 1.2760;
	gamma1 = 0.1319;
	gamma2 = 3.7035;
	delta1 = 2.9496;
	delta2 = 2.9190;
	a = -0.5767;
	c = -79.8678;
	d = -16.6169;
	Q = 1.0;
	Tc = 126.2;
	Tr = T / Tc;
	pc = 33.9 * 10000.; /* bar->Pascal*/
	Pr = p / pc;

	A = alpha1 / Tr * exp(alpha2 * pow(Tr, a));
	B = A * (beta1 * Tr - beta2);
	C = gamma1 / Tr * exp(gamma2 * pow(Tr, c));
	D = delta1 / Tr * exp(delta2 * pow(Tr, d));

	/* Poise->Pa*s */
	my0 = 26.69 * sqrt(28.96) * sqrt(T) / (3.7 * 3.7) * 1.e-6 * 0.1;
	my = my0 * (1.0 + Q * (A * pow(Pr, 1.5)) / (B * Pr + (1 / (1 + C * pow(Pr, D)))));
	return my;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Dynamische Fluessigkeits-Viskositaet nach Yaws et al. (1976)
   als Funktion von Temperatur
   in Reid et al. (1988), S. 441/455
   Eqn.(3): ln(my) = A + B/T + CT + DT^2
   Programing:
   08/2004 OK MFP implementation
           based on CalcFluidViscosityMethod8 by OK (06/2001)
   last modification:
**************************************************************************/
double CFluidProperties::LiquidViscosity_Yaws_1976(double T)
{
	double ln_my, my;
	double A, B, C, D;

	A = -2.471e+01;
	B = 4.209e+03;
	C = 4.527e-02;
	D = -3.376e-5;

	// A = -11.6225;     // CB values for water viscosity after YAWS et al. (1976)
	// B =  1.9490e+03;
	// C =  2.1641e-02;
	// D = -1.5990e-05;

	if (T > 0)
	{
		ln_my = A + B / T + C * T + D * T * T;
		my = exp(ln_my); /* in cP */
		// my = pow(10, ln_my);
		my = my * 1.e-3; /* in Pa s */
	}
	else
		my = 1e-3;
	return my;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Fluessigkeits-Viskositaet in Abhaengigkeit von der Temperatur
   (nach Marsily 1986)
   Programing:
   08/2004 OK MFP implementation
           based on CalcFluidViscosityMethod9 by OK (05/2001)
   last modification:
**************************************************************************/
double CFluidProperties::LiquidViscosity_Marsily_1986(double T)
{
	const double A = 2.29E-03, B = -1.01E-03;

	// MW:	coefficients are wrong!
	//		log(T) gives natural logarithm, here decadal log is needed.
	//		There is no benchmark for this, so nobody noticed...
	// my = 2.285e-5 + 1.01e-3 * log(T);

	if (T > 0)
		return A + B * log10(T); // T in Celsius needed
	else
		return 0.001758784;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Liefert die Viskositaet in Abhaengigkeit von der Konzentration
   und der Temperatur.
   Programing:
   02/2000 AH Erste Version
   08/2004 OK MFP implementation
   last modification:
**************************************************************************/
double CFluidProperties::LiquidViscosity_NN(double c, double T)
{
	double f1, f2, mu0 = 0.001, mu;
	double omega0, omega, sigma0, sigma;

	if (rho_0 < MKleinsteZahl || T_0 < MKleinsteZahl)
		return 0.;

	omega = c / rho_0;
	omega0 = C_0 / rho_0;
	sigma = (T - T_0) / T_0;
	sigma0 = 0.;

	f1 = (1. + 1.85 * omega0 - 4.1 * omega0 * omega0 + 44.5 * omega0 * omega0 * omega0)
	     / (1. + 1.85 * omega - 4.1 * omega * omega + 44.5 * omega * omega * omega);
	f2 = (1 + 0.7063 * sigma - 0.04832 * sigma * sigma * sigma)
	     / (1 + 0.7063 * sigma0 - 0.04832 * sigma0 * sigma * sigma0);
	mu = mu0 / (f1 + f2);
	return mu;
}

////////////////////////////////////////////////////////////////////////////
// Fluid diffusion: required for Phase transition rate coefficients (NAPL)
/**************************************************************************
FEMLib-Method:
Task: Master calc function
Programming:
06/2013 SP Implementation
**************************************************************************/
double CFluidProperties::PhaseDiffusion(double* variables)
{
	static double diff;
	// int fct_number = 0;
	// int gueltig;
	// double mfp_arguments[2];
	// double density;

	// double TTT=0, PPP=0;
	// long Element_Index;
	// int nod_index;

	// string val_name;
	// CRFProcess* m_pcs = NULL;

	//----------------------------------------------------------------------
	// bool New = false;                              // To be
	// if(fem_msh_vector.size()>0) New = true;

	//----------------------------------------------------------------------
	if (variables) // OK4709: faster data access
	{
		primary_variable[0] = variables[0]; // p (single phase)// SP: has to be T
		primary_variable[1] = variables[1]; // T (temperature) // SP: has to be p
		primary_variable[2] = variables[2]; // C (salinity)//index in case of AIR_FLOW
	}
	else
	{
		CalPrimaryVariable(phase_diffusion_pcs_name_vector);
	}
	//----------------------------------------------------------------------
	switch (diffusion_model)
	{
		// case 0:                                     // rho = f(x)
		//   diff = GetCurveValue(fct_number,0,primary_variable[0],&gueltig);
		//   break;
		case 1: // my = const
			diff = diffusion;
			break;
		case 2: // my(p) = my_0*(1+gamma_p*(p-p_0))
			diff = PhaseDiffusion_Yaws_1976(primary_variable[1]);
			break;
		default:
			cout << "Error in CFluidProperties::PhaseDiffusion: no valid model" << endl;
			break;
	}
	//----------------------------------------------------------------------

	return diff;
}

/**************************************************************************
FEMLib-Method:
Task:
   Aquatic diffusion coefficient after Yaws et al. (1976)
   (as function of temperature)
   in YAWS et al. (1976), p. 505
   Eqn.: log10 Daq = A + B/T [cm2/sec]
Programming:
06/2013 SP Implementation
**************************************************************************/
double CFluidProperties::PhaseDiffusion_Yaws_1976(double T)
{
	double A, B;

	A = A_Daq;
	B = B_Daq;

	diff = pow(10, (A + (B / T))) / 10000; // cm�/s -> m�/s
	return diff;
}

////////////////////////////////////////////////////////////////////////////
// Fluid thermal properties
/**************************************************************************
   FEMLib-Method:
   Task: Master calc function
   Programing:
   08/2004 OK MFP implementation based on MATCalcFluidHeatCapacity (OK)
   10/2005 WW/YD Case 3: simplified phase change
   10/2005 YD Case 4: improved phase change
   11/2005 CMCD edited cases and expanded case 3 & 4
   last modification: NB Jan 09 4.9.05
**************************************************************************/
// NB
double CFluidProperties::SpecificHeatCapacity(double* variables)
{
	CRFProcess* m_pcs;
	m_pcs = PCSGet("MULTI_COMPONENTIAL_FLOW");
	double Cp = 0.0;
	double x[2], Cp_c[2];

	int gueltig = -1;
	double pressure, saturation, temperature;

	if (variables) // NB Jan 09
	{
		primary_variable[0] = variables[0]; // p (single phase)
		primary_variable[1] = variables[1]; // T (temperature)
		primary_variable[2] = variables[2]; // ELE index
	}
	else
		CalPrimaryVariable(specific_heat_capacity_pcs_name_vector);

	pressure = primary_variable[0];
	temperature = primary_variable[1];
	saturation = primary_variable[2];
	//......................................................................
	//
	switch (heat_capacity_model)
	{
		case 0: // c = f(x)
			specific_heat_capacity = GetCurveValue(0, 0, temperature, &gueltig);
			break;
		case 1: // c = const, value already read in to specific_heat_capacity
			break;
		case 2: // c = f(p,T,Conc)
			specific_heat_capacity = MATCalcFluidHeatCapacityMethod2(pressure, temperature, saturation);
			break;
		case 5:
			specific_heat_capacity = GetCurveValue(heat_phase_change_curve, 0, temperature_buffer, &gueltig);
			break;
		case 9:
			specific_heat_capacity = isobaric_heat_capacity(Density(primary_variable), primary_variable[1], fluid_id);
			break;
		case 11: // mole fraction weighted average of molar isochoric specific heat capacities converted into isobaric
			// mixture specific heat capacity (ideal conversion)
			{
				// reactive component
				const double M1 = cp_vec[0]->molar_mass;
				const double M2 = cp_vec[1]->molar_mass;
				const double p = primary_variable[0];
				const double T = primary_variable[1];
				const double X = primary_variable[2];

				x[0] = M1 * X / (M1 * X + M2 * (1.0 - X)); // mass in mole fraction
				therm_prop("W");
				const double R = PhysicalConstant::IdealGasConstant;
				Cp_c[0] = isochoric_heat_capacity(M2 * p / (R * T), T, cp_vec[1]->fluid_id);
				// inert component
				x[1] = 1.0 - x[0];
				therm_prop("N");
				Cp_c[1] = isochoric_heat_capacity(M1 * p / (R * T), T, cp_vec[0]->fluid_id);
				specific_heat_capacity
				    = Cp_c[0] * M2 * x[0] + Cp_c[1] * M1 * x[1]; // mixture isochoric molar heat capacities
				specific_heat_capacity += R; // isochoric in isobaric
				specific_heat_capacity /= (M1 * x[1] + M2 * x[0]); // molar in specific of mixture value
				break;
			}
		case 12: // mass fraction weighted average of isobaric specific heat capacities using a linearised model
			// reactive component
			x[0] = primary_variable[2]; // mass fraction
			Cp_c[0] = linear_heat_capacity(variables[1], cp_vec[1]->fluid_id);
			// inert component
			x[1] = 1.0 - x[0];
			Cp_c[1] = linear_heat_capacity(variables[1], cp_vec[0]->fluid_id);
			specific_heat_capacity = Cp_c[0] * x[0] + Cp_c[1] * x[1]; // mixture isobaric specific heat capacities
			break;
		case 13: //mass fraction weighted average of isobaric specific heat capacities using a polynomial model
			//reactive component
			x[0] = primary_variable[2]; //mass fraction
			Cp_c[0] = polynomial_heat_capacity(variables[1],cp_vec[1]->fluid_id);
			//inert component
			x[1] = 1.0 - x[0];
			Cp_c[1] = polynomial_heat_capacity(variables[1],cp_vec[0]->fluid_id);
			specific_heat_capacity = Cp_c[0]*x[0] + Cp_c[1]*x[1]; //mixture isobaric specific heat capacities
			break;
		case 15: // mixture cp= sum_i y_i*cp:: P, T, x dependent
			for (int CIndex = 2; CIndex < cmpN + 2; CIndex++)
			{
				if (eos_name == "CONSTANT")
				{
					Cp += variables[CIndex] * cp[CIndex - 2];
				}
				else
				{
					therm_prop(m_pcs->pcs_primary_function_name[CIndex]);
					Cp += variables[CIndex]
					      * isobaric_heat_capacity(ComponentDensity(CIndex, variables), variables[1], fluid_id);
				}
			}
			specific_heat_capacity = Cp;

			break;
	}
	return specific_heat_capacity;
}

/**************************************************************************
   FEMLib-Method:
   Task: calculate heat capacity for phase change
   Programing:
   02/2008 JOD moved from CFluidProperties::SpecificHeatCapacity()
   last modification:
**************************************************************************/
double CFluidProperties::PhaseChange()
{
	int gueltig = -1;
	double heat_capacity_phase_change = 0;
	double pressure; // WW , saturation, temperature;
	double density_vapor, humi, drdT;
	double H1, H0, T0, T_1; // T1 defined at 662 in steam67???

	//......................................................................
	CalPrimaryVariable(specific_heat_capacity_pcs_name_vector);
	pressure = primary_variable[0];
	// WW temperature =  primary_variable[1];
	// WW saturation = primary_variable[2];

	if (heat_capacity_model == 3)
	{
		T_1 = primary_variable_t1[1];
		if (T_1 <= T_Latent1 || T_1 >= T_Latent2)
			heat_capacity_phase_change = GetCurveValue(heat_phase_change_curve, 0, temperature_buffer, &gueltig);
		else
		{
			heat_capacity_model = 5; // ??? JOD
			H1 = CalcEnthalpy(T_1);
			T0 = primary_variable_t0[1];
			if (fabs(T_1 - T0) < 1.0e-8)
				T_1 += 1.0e-8;
			H0 = CalcEnthalpy(T0);
			heat_capacity_phase_change = (H1 - H0) / (T_1 - T_0);
			heat_capacity_model = 3;
		}
	}
	else if (heat_capacity_model == 4)
	{
		T_1 = primary_variable_t1[1];
		if (T_1 <= T_Latent1 || T_1 >= T_Latent2)
		{
			humi
			    = exp(pressure / (PhysicalConstant::SpecificGasConstant::WaterVapour * temperature_buffer * Density()));
			density_vapor = humi * Density();
			drdT = (vaporDensity_derivative(temperature_buffer) * humi
			        - density_vapor * pressure / (PhysicalConstant::SpecificGasConstant::WaterVapour * Density()
			                                      * (temperature_buffer * temperature_buffer)))
			       / Density();
			H1 = latent_heat + specific_heat_capacity * (temperature_buffer - T_Latent1);
			heat_capacity_phase_change = H1 * drdT;
		}
		else
		{
			heat_capacity_model = 5;
			H1 = CalcEnthalpy(T_1);
			T0 = primary_variable_t0[1];
			if (fabs(T_1 - T0) < 1.0e-8)
				T_1 += 1.0e-8;
			H0 = CalcEnthalpy(T0);
			heat_capacity_phase_change = (H1 - H0) / (T_1 - T0);
			heat_capacity_model = 4;
		}
	}

	return heat_capacity_phase_change;
}

/**************************************************************************
   FEMLib-Method:
   Task: Master calc function
   Programing:
   02/2007 WW MFP implementation based on MATCalcFluidHeatCapacity (OK)
**************************************************************************/
double MFPCalcFluidsHeatCapacity(CFiniteElementStd* assem)
{
	double heat_capacity_fluids = 0.0;
	double PG = 0.0, Sw = 0.0, TG, rhow, rho_gw, p_gw, dens_aug[3], rho_g;
	CFluidProperties* m_mfp = NULL;
	CRFProcess* m_pcs = assem->cpl_pcs;
	// if (m_pcs->pcs_type_name.find("MULTI_PHASE_FLOW")!=string::npos)
	if (m_pcs && m_pcs->type == 1212) // non-isothermal multi-phase flow
	{
		// Capillary pressure
		PG = assem->interpolate(assem->NodalValC1);
		Sw = assem->MediaProp->SaturationCapillaryPressureFunction(PG);
		double PG2 = assem->interpolate(assem->NodalVal_p2);
		TG = assem->interpolate(assem->NodalVal1);
		rhow = assem->FluidProp->Density();
		rho_gw = assem->FluidProp->vaporDensity(TG) * exp(-PG / (rhow * SpecificGasConstant::WaterVapour * TG));
		p_gw = rho_gw * SpecificGasConstant::WaterVapour * TG;
		dens_aug[0] = PG2 - p_gw;
		dens_aug[1] = TG;
		m_mfp = mfp_vector[1];
		// 2 Dec 2010 AKS
		rho_g = rho_gw + m_mfp->Density(dens_aug);
		//
		m_mfp = mfp_vector[0];
		heat_capacity_fluids = Sw * m_mfp->Density() * m_mfp->SpecificHeatCapacity();
		m_mfp = mfp_vector[1];
		heat_capacity_fluids += (1.0 - Sw) * rho_g * m_mfp->SpecificHeatCapacity();
	}

	else
	{
		heat_capacity_fluids = assem->FluidProp->Density() * assem->FluidProp->SpecificHeatCapacity();

		if (m_pcs && m_pcs->type != 1) // neither liquid nor ground water flow
		{
			//  pressure
			PG = assem->interpolate(assem->NodalValC1);

			if (PG < 0.0)
			{
				Sw = assem->MediaProp->SaturationCapillaryPressureFunction(-PG);
				heat_capacity_fluids *= Sw;
				if (assem->GasProp != 0)
					heat_capacity_fluids
					    += (1. - Sw) * assem->GasProp->Density() * assem->GasProp->SpecificHeatCapacity();
				heat_capacity_fluids += (1. - Sw) * assem->FluidProp->PhaseChange();
			}
		}
	}
	return heat_capacity_fluids;
}

/**************************************************************************
   FEMLib-Method:
   Task: Master calc function
   Programing:
   08/2004 OK MFP implementation based on MATCalcFluidHeatCapacity (OK)
   11/2005 YD Modification
   last modification:
   NB 4.9.05
**************************************************************************/
// NB Dec 08 4.9.05
double CFluidProperties::HeatConductivity(double* variables)
{
	if (variables) // NB Dec 08
	{
		primary_variable[0] = variables[0]; // p (single phase)
		primary_variable[1] = variables[1]; // T (temperature)
		primary_variable[2] = variables[2]; // ELE index
	}
	else
		CalPrimaryVariable(heat_conductivity_pcs_name_vector);

	switch (heat_conductivity_model)
	{
		case 0: // rho = f(x)
		{
			int fct_number = 0;
			int gueltig;
			heat_conductivity = GetCurveValue(fct_number, 0, primary_variable[0], &gueltig);
		}
		break;
		case 1: // c = const
			break;
		case 2:
			heat_conductivity
			    = MATCalcHeatConductivityMethod2(primary_variable[0], primary_variable[1], primary_variable[2]);
			break;
		case 3: // NB
			heat_conductivity = Fluid_Heat_Conductivity(Density(), primary_variable[1], fluid_id);
			// if (heat_conductivity<0.03) // not sure about this
			break;
		case 9:
			heat_conductivity = Fluid_Heat_Conductivity(Density(primary_variable), primary_variable[1], fluid_id);
			break;

		case 11: // Wassilijewa, Maso&Saxena (see Poling, B. E.; Prausnitz, J. M.; John Paul, O. & Reid, R. C. The
			// properties of gases and liquids McGraw-Hill New York, 2001, 5: page 10.30f.)
			{
				double x[2], k[2];
				const double M0 = cp_vec[0]->molar_mass;
				const double M1 = cp_vec[1]->molar_mass;
				const double p = primary_variable[0];
				const double T = primary_variable[1];
				const double X = primary_variable[2];

				// TODO [CL] max() is redundant if the fraction is guaranteed to be between 0 and 1.
				// reactive component
				const double R = PhysicalConstant::IdealGasConstant;
				x[0] = max(M0 * X / (M0 * X + M1 * (1.0 - X)), 0.); // convert mass to mole fraction
				k[0] = Fluid_Heat_Conductivity(M1 * p / (R * T), T, cp_vec[1]->fluid_id);
				// inert component
				x[1] = 1.0 - x[0];
				k[1] = Fluid_Heat_Conductivity(M0 * p / (R * T), T, cp_vec[0]->fluid_id);

				const double M1_over_M2 = M1 / M0; // reactive over inert
				const double V1_over_V2 = Fluid_Viscosity(M1 * p / (R * T), T, p, cp_vec[1]->fluid_id)
				                          / Fluid_Viscosity(M0 * p / (R * T), T, p, cp_vec[0]->fluid_id);
				const double L1_over_L2 = V1_over_V2 / M1_over_M2;

				const double phi_12 = (1.0 + pow(L1_over_L2, 0.5) * pow(M1_over_M2, -0.25))
				                      * (1.0 + pow(V1_over_V2, 0.5) * pow(M1_over_M2, -0.25))
				                      / pow(8.0 * (1.0 + M1_over_M2), 0.5);
				const double phi_21 = phi_12 * M1_over_M2 / V1_over_V2;

				heat_conductivity = k[0] * x[0] / (x[0] + x[1] * phi_12);
				heat_conductivity += k[1] * x[1] / (x[1] + x[0] * phi_21);
				break;
			}
		case 15: // mixture k_m= sum_i y_i*k_i:: p, T, x
		{
			CRFProcess* m_pcs = PCSGet("MULTI_COMPONENTIAL_FLOW");
			double Kappa = 0.0;

			for (int CIndex = 2; CIndex < cmpN + 2; CIndex++)
			{
				if (eos_name == "CONSTANT")
				{
					Kappa += variables[CIndex] * kappa[CIndex - 2];
				}
				else
				{
					therm_prop(m_pcs->pcs_primary_function_name[CIndex]);
					Kappa += variables[CIndex]
					         * Fluid_Heat_Conductivity(ComponentDensity(CIndex, variables), variables[1], fluid_id);
				}
			}
			heat_conductivity = Kappa;
		}
		break;
	}

	return heat_conductivity;
}

/**************************************************************************
   FEMLib-Method:
   Task: Master calc function
   Programing:
   08/2004 OK MFP implementation based on MATCalcFluidHeatCapacity (OK)
   last modification:
   10/2010 TF changed access to process type
**************************************************************************/
double MFPCalcFluidsHeatConductivity(long index, double* gp, double theta, CFiniteElementStd* assem)
{
	gp = gp; // OK411
	index = index;
	double saturation_phase = 0.0; // OK411
	double heat_conductivity_fluids = 0.0;
	int nidx0, nidx1;
	bool New = false; // To be removed. WW
	if (fem_msh_vector.size() > 0)
		New = true;

	//--------------------------------------------------------------------
	//----------------------------------------------------------------------
	CFluidProperties* m_mfp = NULL;
	int no_fluids = (int)mfp_vector.size();
	// YD-----------
	CRFProcess* m_pcs = NULL;
	for (int i = 0; i < (int)pcs_vector.size(); i++)
	{
		m_pcs = pcs_vector[i];
		//    if(m_pcs->pcs_type_name.find("RICHARDS_FLOW"))
		if (m_pcs->getProcessType() == FiniteElement::RICHARDS_FLOW)
			no_fluids = 1;
	}
	// YD-----------
	switch (no_fluids)
	{
		case 1:
			m_mfp = mfp_vector[0];
			heat_conductivity_fluids = m_mfp->HeatConductivity();
			break;
		case 2:
			m_mfp = mfp_vector[0];
			if (New) // WW
			{
				nidx0 = m_pcs->GetNodeValueIndex("SATURATION1");
				nidx1 = nidx0 + 1;
				saturation_phase
				    = (1. - theta) * assem->interpolate(nidx0, m_pcs) + theta * assem->interpolate(nidx1, m_pcs);
			}
			heat_conductivity_fluids = saturation_phase * m_mfp->HeatConductivity();
			m_mfp = mfp_vector[1];
			heat_conductivity_fluids += (1.0 - saturation_phase) * m_mfp->HeatConductivity();
			break;
		default:
			cout << "Error in MFPCalcFluidsHeatConductivity: no fluid phase data"
			     << "\n";
	}
	return heat_conductivity_fluids;
}

double CFluidProperties::LiquidViscosity_expo(double T)
{
	return viscosity0 * exp(-(T - T_0) / viscosity_T_star);
}

/**************************************************************************
   FEMLib-Method:
   Task:

   Programing:

   last modification:
**************************************************************************/
double CFluidProperties::LiquidViscosity_CMCD(double Press, double TempK, double C)
{
	C = C;
	/*CMcD variables for 20 ALR*/
	double A1, A2, A3, A4, A5, A6, A7, A8; /*constants*/
	double TempC, TempF, Pbar, Salinity; /* Temperature [K], Temperature [F], Pressure [bar]*/
	double my_Zero, PsatBar,
	    PsatKPa; /*my pure water, my saline water, Saturation pressure [bar], Saturation pressure [KPa]*/
	/*intermediate values*/
	double sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, exponent;
	/*CMcD end variables for 20 ALR*/
	/* //Prepared for introduction of solute transport in PCS version
	   //Average Concentration
	   comp=0; // nur fuer Einkomponenten-Systeme
	   timelevel=1;
	   concentration_average = 0.0;
	   for (i = 0; i < count_nodes; i++)
	    concentration_average += GetNodeVal(element_nodes[i], c_idx);
	   concentration_average /= count_nodes;
	   fp->salinity=fp->rho_0+concentration_average*fp->drho_dC;
	 */
	// Link to function from ALR
	Salinity = C_0 / 1000.;
	/***constant values*******/
	A1 = -7.419242;
	A2 = -0.29721;
	A3 = -0.1155286;
	A4 = -0.008685635;
	A5 = 0.001094098;
	A6 = 0.00439993;
	A7 = 0.002520658;
	A8 = 0.0005218684;

	/*Unit conversions*/
	TempC = TempK - PhysicalConstant::CelsiusZeroInKelvin;
	TempF = TempC * 1.8 + 32.0;
	Pbar = Press / 100000.0;
	/*end of units conversion*/

	/*Calculation of the saturation pressure*/
	//   sum1=pow((0.65-0.01*TempK),0)*A1;
	sum1 = 1.0 * A1;
	//   sum2=pow((0.65-0.01*TempK),1)*A2;
	sum2 = (0.65 - 0.01 * TempK) * A2;
	sum3 = (0.65 - 0.01 * TempK) * (0.65 - 0.01 * TempK) * A3;
	sum4 = MathLib::fastpow((0.65 - 0.01 * TempK), 3) * A4;
	sum5 = MathLib::fastpow((0.65 - 0.01 * TempK), 4) * A5;
	sum6 = MathLib::fastpow((0.65 - 0.01 * TempK), 5) * A6;
	sum7 = MathLib::fastpow((0.65 - 0.01 * TempK), 6) * A7;
	sum8 = MathLib::fastpow((0.65 - 0.01 * TempK), 7) * A8;

	/*intermediate value*/
	exponent = sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7 + sum8;
	exponent = exponent * (374.136 - TempC) / TempK; /*intermediate value*/

	PsatKPa = exp(exponent) * 22088; /*saturation pressure in kPa*/
	PsatBar = PsatKPa / (1000 * 100000); /*Saturation pressure in bar*/

	/*Viscosity of pure water in Pa-S*/
	my_Zero = 243.18e-7 * (pow(10., (247.8 / (TempK - 140)))) * (1 + (Pbar - PsatBar) * 1.0467e-6 * (TempK - 305));

	/*Viscosity of saline water in Pa-S*/
	viscosity = my_Zero * (1 - 0.00187 * (sqrt(Salinity)) + 0.000218 * (MathLib::fastpow(sqrt(Salinity), 5))
	                       + (sqrt(TempF) - 0.0135 * TempF)
	                             * (0.00276 * Salinity - 0.000344 * (MathLib::fastpow(sqrt(Salinity), 3))));
	return viscosity;
}

/**************************************************************************/
/* ROCKFLOW - Function: MATCalcFluidHeatConductivityMethod2
 */
/* Task:
   Calculate heat conductivity of all fluids
 */
/* Parameter: (I: Input; R: Return; X: Both)
   I double temperature
 */
/* Return:
   Value of heat conductivity of fluids as a function of (p,C)
 */
/* Programming:
   09/2003   CMCD ARL   Implementation
   08/2004   CMCD inclusion in GeoSys v. 4.
 */
/**************************************************************************/
double CFluidProperties::MATCalcHeatConductivityMethod2(double Press, double TempK, double Conc)
{
	Conc = Conc;
	int i, j;
	double TauTC, Nabla, Delta, Nabla0, Nabla1, Nabla2;
	double heat_conductivity, Rho, temperature_average, pressure_average, viscosity;
	double Rhostar, TstarTC, Lambdastar; // WW, Pstar1;
	double nZero[4];
	double n[5][6];
	double A1, A2, A3, A4, A5, A6, A7, A8; /*constants*/
	double TempC, TempF, Pbar, Salinity; /* Temperature [K], Temperature [F], Pressure [bar]*/
	double my_Zero, PsatBar,
	    PsatKPa; /*my pure water, my saline water, Saturation pressure [bar], Saturation pressure [KPa]*/
	/*intermediate values*/
	double sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, exponent;

	/*************************************************************************************************/
	/*************************************************************************************************/
	/*************************Partial derivatives calculation*****************************************/
	double GammaPi, GammaPiTau, GammaPiPi, Pii, Tau, GazConst;
	double LGamma[35];
	double JGamma[35];
	double nGamma[35];
	double Tstar, Pstar;
	double TstarTilda, PstarTilda, RhostarTilda;
	double First_derivative, Second_derivative;

	pressure_average = Press;
	temperature_average = TempK;
	Salinity = C_0;
	Tstar = 1386;
	Pstar = 16.53e6; // MPa
	GazConst = 0.461526e3; //!!!!Given by equation (2.1)
	TstarTilda = 647.226;
	PstarTilda = 22.115e6;
	RhostarTilda = 317.763;
	Lambdastar = 0.4945;
	Rhostar = 317.763;
	TstarTC = 647.226;
	// WW Pstar1 = 22.115e-6;

	/*BEGIN: reduced dimensions*/
	TauTC = TstarTC / temperature_average;
	// Delta = Rho / Rhostar;
	// WW PiiTC = pressure_average / Pstar1;
	/*END: reduced dimensions*/

	nGamma[1] = 0.14632971213167;
	nGamma[2] = -0.84548187169114;
	nGamma[3] = -0.37563603672040e1;
	nGamma[4] = 0.33855169168385e1;
	nGamma[5] = -0.95791963387872;
	nGamma[6] = 0.15772038513228;
	nGamma[7] = -0.16616417199501e-1;
	nGamma[8] = 0.81214629983568e-3;
	nGamma[9] = 0.28319080123804e-3;
	nGamma[10] = -0.60706301565874e-3;
	nGamma[11] = -0.18990068218419e-1;
	nGamma[12] = -0.32529748770505e-1;
	nGamma[13] = -0.21841717175414e-1;
	nGamma[14] = -0.52838357969930e-4;
	nGamma[15] = -0.47184321073267e-3;
	nGamma[16] = -0.30001780793026e-3;
	nGamma[17] = 0.47661393906987e-4;
	nGamma[18] = -0.44141845330846e-5;
	nGamma[19] = -0.72694996297594e-15;
	nGamma[20] = -0.31679644845054e-4;
	nGamma[21] = -0.28270797985312e-5;
	nGamma[22] = -0.85205128120103e-9;
	nGamma[23] = -0.22425281908000e-5;
	nGamma[24] = -0.65171222895601e-6;
	nGamma[25] = -0.14341729937924e-12;
	nGamma[26] = -0.40516996860117e-6;
	nGamma[27] = -0.12734301741641e-8;
	nGamma[28] = -0.17424871230634e-9;
	nGamma[29] = -0.68762131295531e-18;
	nGamma[30] = 0.14478307828521e-19;
	nGamma[31] = 0.26335781662795e-22;
	nGamma[32] = -0.11947622640071e-22;
	nGamma[33] = 0.18228094581404e-23;
	nGamma[34] = -0.93537087292458e-25;

	LGamma[1] = 0.;
	LGamma[2] = 0.;
	LGamma[3] = 0.;
	LGamma[4] = 0.;
	LGamma[5] = 0.;
	LGamma[6] = 0.;
	LGamma[7] = 0.;
	LGamma[8] = 0.;
	LGamma[9] = 1.;
	LGamma[10] = 1.;
	LGamma[11] = 1.;
	LGamma[12] = 1.;
	LGamma[13] = 1.;
	LGamma[14] = 1.;
	LGamma[15] = 2.;
	LGamma[16] = 2.;
	LGamma[17] = 2.;
	LGamma[18] = 2.;
	LGamma[19] = 2.;
	LGamma[20] = 3.;
	LGamma[21] = 3.;
	LGamma[22] = 3.;
	LGamma[23] = 4.;
	LGamma[24] = 4.;
	LGamma[25] = 4.;
	LGamma[26] = 5.;
	LGamma[27] = 8.;
	LGamma[28] = 8.;
	LGamma[29] = 21.;
	LGamma[30] = 23.;
	LGamma[31] = 29.;
	LGamma[32] = 30.;
	LGamma[33] = 31.;
	LGamma[34] = 32.;

	JGamma[1] = -2.;
	JGamma[2] = -1.;
	JGamma[3] = 0.;
	JGamma[4] = 1.;
	JGamma[5] = 2.;
	JGamma[6] = 3.;
	JGamma[7] = 4.;
	JGamma[8] = 5.;
	JGamma[9] = -9.;
	JGamma[10] = -7.;
	JGamma[11] = -1.;
	JGamma[12] = 0.;
	JGamma[13] = 1.;
	JGamma[14] = 3.;
	JGamma[15] = -3.;
	JGamma[16] = 0.;
	JGamma[17] = 1.;
	JGamma[18] = 3.;
	JGamma[19] = 17.;
	JGamma[20] = -4.;
	JGamma[21] = 0.;
	JGamma[22] = 6.;
	JGamma[23] = -5.;
	JGamma[24] = -2.;
	JGamma[25] = 10.;
	JGamma[26] = -8.;
	JGamma[27] = -11.;
	JGamma[28] = -6.;
	JGamma[29] = -29.;
	JGamma[30] = -31.;
	JGamma[31] = -38.;
	JGamma[32] = -39.;
	JGamma[33] = -40.;
	JGamma[34] = -41.;

	Pii = pressure_average / Pstar;
	Tau = Tstar / temperature_average;

	/*BEGIN:Calculation of GammaPi*/
	GammaPi = 0;

	for (i = 1; i < 35; i++)
		GammaPi = GammaPi
		          - (nGamma[i]) * (LGamma[i]) * (pow((7.1 - Pii), (LGamma[i] - 1))) * (pow((Tau - 1.222), JGamma[i]));
	/*END: Calculation of GammaPi*/

	/*BEGIN:Calculation of GammaPiTau*/
	GammaPiTau = 0;
	for (i = 1; i < 35; i++)
		GammaPiTau = GammaPiTau
		             - (nGamma[i]) * (LGamma[i]) * (pow((7.1 - Pii), (LGamma[i] - 1))) * (JGamma[i])
		                   * (pow((Tau - 1.222), (JGamma[i] - 1)));
	/*END: Calculation of GammaPiTau*/

	/*BEGIN:Calculation of GammaPiPi*/
	GammaPiPi = 0;
	for (i = 1; i <= 34; i++)
		GammaPiPi = GammaPiPi
		            + (nGamma[i]) * (LGamma[i]) * (LGamma[i] - 1) * (pow((7.1 - Pii), (LGamma[i] - 2)))
		                  * (pow((Tau - 1.222), (JGamma[i])));
	/*END: Calculation of GammaPiPi*/

	/*BEGIN:Calculation of derivative*/
	First_derivative = ((TstarTilda) * (Pstar) * ((GammaPiTau) * (Tstar) - (GammaPi) * (temperature_average)))
	                   / (PstarTilda * temperature_average * temperature_average * GammaPiPi),
	Second_derivative = ((-1) * (PstarTilda) * (GammaPiPi))
	                    / ((RhostarTilda) * (temperature_average) * (GazConst) * ((GammaPi * GammaPi)));
	/*End:Calculation of derivative*/

	/*BEGIN: Calculation of density*/
	Rho = pressure_average / (GazConst * (temperature_average) * (Pii) * (GammaPi));
	/*END: Calculation of density*/

	/*************************Partial derivatives calculation*****************************************/
	/*************************************************************************************************/
	/*************************************************************************************************/

	/*BEGIN: Constant values*/

	Lambdastar = 0.4945;
	Rhostar = 317.763;
	TstarTC = 647.226;
	// WW Pstar1 = 22.115e6;

	nZero[0] = 0.1e1;
	nZero[1] = 0.6978267e1;
	nZero[2] = 0.2599096e1;
	nZero[3] = -0.998254;

	n[0][0] = 0.13293046e1;
	n[0][1] = -0.40452437;
	n[0][2] = 0.24409490;
	n[0][3] = 0.18660751e-1;
	n[0][4] = -0.12961068;
	n[0][5] = 0.44809953e-1;

	n[1][0] = 0.17018363e1;
	n[1][1] = -0.22156845e1;
	n[1][2] = 0.16511057e1;
	n[1][3] = -0.76736002;
	n[1][4] = 0.37283344;
	n[1][5] = -0.11203160;

	n[2][0] = 0.52246158e1;
	n[2][1] = -0.10124111e2;
	n[2][2] = 0.49874687e1;
	n[2][3] = -0.27297694;
	n[2][4] = -0.43083393;
	n[2][5] = 0.13333849;

	n[3][0] = 0.87127675e1;
	n[3][1] = -0.95000611e1;
	n[3][2] = 0.43786606e1;
	n[3][3] = -0.91783782;
	n[3][4] = 0;
	n[3][5] = 0;

	n[4][0] = -0.18525999e1;
	n[4][1] = 0.93404690;
	n[4][2] = 0;
	n[4][3] = 0;
	n[4][4] = 0;
	n[4][5] = 0;
	/*END: Constant values*/

	/*BEGIN: reduced dimensions*/
	TauTC = TstarTC / temperature_average;
	Delta = Rho / Rhostar;
	// WW PiiTC = pressure_average / Pstar1;
	/*END: reduced dimensions*/

	/*BEGIN: Nabla0*/
	Nabla0 = 0;
	for (i = 0; i <= 3; i++)
		Nabla0 = Nabla0 + (nZero[i]) * (MathLib::fastpow(TauTC, i));

	Nabla0 = Nabla0 * (sqrt(TauTC));
	Nabla0 = 1 / Nabla0;
	/*END: Nabla0*/

	/*BEGIN: Nabla1*/
	Nabla1 = 0;
	for (i = 0; i <= 4; i++)
	{
		const double t(MathLib::fastpow(TauTC - 1, i)); // TF
		for (j = 0; j <= 5; j++)
			Nabla1 = Nabla1 + (n[i][j]) * (t * MathLib::fastpow(Delta - 1, j));
	}

	Nabla1 = Delta * (Nabla1);
	Nabla1 = exp(Nabla1);
	/*END: Nabla1*/

	/*Calculate Viscosity*/
	/*Link to function from ALR*/

	TempK = temperature_average;
	Press = pressure_average;
	Salinity = C_0 / 1000.;

	/***constant values*******/
	A1 = -7.419242;
	A2 = -0.29721;
	A3 = -0.1155286;
	A4 = -0.008685635;
	A5 = 0.001094098;
	A6 = 0.00439993;
	A7 = 0.002520658;
	A8 = 0.0005218684;

	/*Unit conversions*/
	TempC = TempK - PhysicalConstant::CelsiusZeroInKelvin;
	TempF = TempC * 1.8 + 32.0;
	if (TempF < 0.0)
		TempF = 0.0;
	Pbar = Press / 100000.0;
	/*end of units conversion*/

	/*Calculation of the saturation pressure*/
	// TF   sum1=pow((0.65-0.01*TempK),0)*A1;
	sum1 = A1;
	// TF  sum2=pow((0.65-0.01*TempK),1)*A2;
	sum2 = (0.65 - 0.01 * TempK) * A2;
	sum3 = (0.65 - 0.01 * TempK) * (0.65 - 0.01 * TempK) * A3;
	sum4 = MathLib::fastpow((0.65 - 0.01 * TempK), 3) * A4;
	sum5 = MathLib::fastpow((0.65 - 0.01 * TempK), 4) * A5;
	sum6 = MathLib::fastpow((0.65 - 0.01 * TempK), 5) * A6;
	sum7 = MathLib::fastpow((0.65 - 0.01 * TempK), 6) * A7;
	sum8 = MathLib::fastpow((0.65 - 0.01 * TempK), 7) * A8;

	/*intermediate value*/
	exponent = sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7 + sum8;
	exponent = exponent * (374.136 - TempC) / TempK; /*intermediate value*/

	PsatKPa = exp(exponent) * 22088; /*saturation pressure in kPa*/
	PsatBar = PsatKPa / (1000 * 100000); /*Saturation pressure in bar*/

	/*Viscosity of pure water in Pa-S*/
	my_Zero = 243.18e-7 * (pow(10., (247.8 / (TempK - 140)))) * (1 + (Pbar - PsatBar) * 1.0467e-6 * (TempK - 305));

	/*Viscosity of saline water in Pa-S*/
	viscosity = my_Zero * (1 - 0.00187 * (sqrt(Salinity)) + 0.000218 * (MathLib::fastpow(sqrt(Salinity), 5))
	                       + (sqrt(TempF) - 0.0135 * TempF)
	                             * (0.00276 * Salinity - 0.000344 * (MathLib::fastpow(sqrt(Salinity), 3))));

	/* End of viscosity function*/

	/*BEGIN: Nabla2*/
	Nabla2 = 0.0013848 / ((viscosity) / 55.071e-6) * (1.0 / (TauTC * Delta) * (TauTC * Delta))
	         * (First_derivative * First_derivative) * (pow((Delta * (Second_derivative)), 0.4678)) * (sqrt(Delta))
	         * exp(-18.66 * ((1 / TauTC - 1) * (1 / TauTC - 1)) - (MathLib::fastpow(Delta - 1, 4)));
	/*END: Nabla2*/

	/*BEGIN: Nabla => heat_conductivity*/
	Nabla = Nabla0 * (Nabla1) + Nabla2;
	heat_conductivity = Nabla * (Lambdastar);
	/*END: Nabla => Lambda*/

	return heat_conductivity;
}

/**************************************************************************
   ROCKFLOW - Funktion:
   Task:
   Programming:
   08/2005 WW Implementation
**************************************************************************/
double CFluidProperties::vaporDensity(const double T_abs)
{
	return 1.0e-3 * exp(19.81 - 4975.9 / T_abs);
}

double CFluidProperties::vaporDensity_derivative(const double T_abs)
{
	return 4.9759 * exp(19.81 - 4975.9 / T_abs) / (T_abs * T_abs);
}

/**************************************************************************/
/* ROCKFLOW - Function: MATCalcFluidHeatCapacity
 */
/* Task:
   Calculate heat capacity of all fluids
 */
/* Parameter: (I: Input; R: Return; X: Both)
   I double temperature
 */
/* Return:
   Value of heat capacity of all fluids
 */
/* Programming:
   09/2003   CMCD ARL   Implementation
   09/2004   CMCD included in Geosys ver. 4

 */
/**************************************************************************/
double CFluidProperties::MATCalcFluidHeatCapacityMethod2(double Press, double TempK, double Conc)
{
	Conc = Conc;
	double Pressurevar, Tau, pressure_average, temperature_average, Tstar, Pstar, GazConst;
	double GammaPi, GammaPiTau, GammaPiPi, GammaTauTau;
	double L[35], J[35], n[35];
	int i;
	// WW double salinity;
	double Cp; // WW , Cv;

	pressure_average = Press;
	temperature_average = TempK;

	// WW salinity=C_0;

	Tstar = 1386;

	Pstar = 16.53e6; // MPa
	GazConst = 0.461526e3; //

	n[0] = 0.0;
	n[1] = 0.14632971213167;
	n[2] = -0.84548187169114;
	n[3] = -0.37563603672040e1;
	n[4] = 0.33855169168385e1;
	n[5] = -0.95791963387872;
	n[6] = 0.15772038513228;
	n[7] = -0.16616417199501e-1;
	n[8] = 0.81214629983568e-3;
	n[9] = 0.28319080123804e-3;
	n[10] = -0.60706301565874e-3;
	n[11] = -0.18990068218419e-1;
	n[12] = -0.32529748770505e-1;
	n[13] = -0.21841717175414e-1;
	n[14] = -0.52838357969930e-4;
	n[15] = -0.47184321073267e-3;
	n[16] = -0.30001780793026e-3;
	n[17] = 0.47661393906987e-4;
	n[18] = -0.44141845330846e-5;
	n[19] = -0.72694996297594e-15;
	n[20] = -0.31679644845054e-4;
	n[21] = -0.28270797985312e-5;
	n[22] = -0.85205128120103e-9;
	n[23] = -0.22425281908000e-5;
	n[24] = -0.65171222895601e-6;
	n[25] = -0.14341729937924e-12;
	n[26] = -0.40516996860117e-6;
	n[27] = -0.12734301741641e-8;
	n[28] = -0.17424871230634e-9;
	n[29] = -0.68762131295531e-18;
	n[30] = 0.14478307828521e-19;
	n[31] = 0.26335781662795e-22;
	n[32] = -0.11947622640071e-22;
	n[33] = 0.18228094581404e-23;
	n[34] = -0.93537087292458e-25;

	L[0] = 0.;
	L[1] = 0.;
	L[2] = 0.;
	L[3] = 0.;
	L[4] = 0.;
	L[5] = 0.;
	L[6] = 0.;
	L[7] = 0.;
	L[8] = 0.;
	L[9] = 1.;
	L[10] = 1.;
	L[11] = 1.;
	L[12] = 1.;
	L[13] = 1.;
	L[14] = 1.;
	L[15] = 2.;
	L[16] = 2.;
	L[17] = 2.;
	L[18] = 2.;
	L[19] = 2.;
	L[20] = 3.;
	L[21] = 3.;
	L[22] = 3.;
	L[23] = 4.;
	L[24] = 4.;
	L[25] = 4.;
	L[26] = 5.;
	L[27] = 8.;
	L[28] = 8.;
	L[29] = 21.;
	L[30] = 23.;
	L[31] = 29.;
	L[32] = 30.;
	L[33] = 31.;
	L[34] = 32.;

	J[0] = -2.;
	J[1] = -2.;
	J[2] = -1.;
	J[3] = 0.;
	J[4] = 1.;
	J[5] = 2.;
	J[6] = 3.;
	J[7] = 4.;
	J[8] = 5.;
	J[9] = -9.;
	J[10] = -7.;
	J[11] = -1.;
	J[12] = 0.;
	J[13] = 1.;
	J[14] = 3.;
	J[15] = -3.;
	J[16] = 0.;
	J[17] = 1.;
	J[18] = 3.;
	J[19] = 17.;
	J[20] = -4.;
	J[21] = 0.;
	J[22] = 6.;
	J[23] = -5.;
	J[24] = -2.;
	J[25] = 10.;
	J[26] = -8.;
	J[27] = -11.;
	J[28] = -6.;
	J[29] = -29.;
	J[30] = -31.;
	J[31] = -38.;
	J[32] = -39.;
	J[33] = -40.;
	J[34] = -41.;

	Pressurevar = pressure_average / Pstar;
	Tau = Tstar / temperature_average;

	/*BEGIN:Calculation of GammaPi*/
	GammaPi = 0;

	for (i = 1; i < 35; i++)
		GammaPi = GammaPi - (n[i]) * (L[i]) * (pow((7.1 - Pressurevar), (L[i] - 1.))) * (pow((Tau - 1.222), J[i]));
	/*END: Calculation of GammaPi*/

	/*BEGIN:Calculation of GammaPiTau*/
	GammaPiTau = 0;
	for (i = 1; i < 35; i++)
		GammaPiTau
		    = GammaPiTau
		      - (n[i]) * (L[i]) * (pow((7.1 - Pressurevar), (L[i] - 1.))) * (J[i]) * (pow((Tau - 1.222), (J[i] - 1.)));
	/*END: Calculation of GammaPiTau*/

	/*BEGIN:Calculation of GammaTauTau*/
	GammaTauTau = 0;
	for (i = 1; i < 35; i++)
		GammaTauTau
		    = GammaTauTau
		      + (n[i]) * (pow((7.1 - Pressurevar), (L[i]))) * (J[i]) * (J[i] - 1.) * (pow((Tau - 1.222), (J[i] - 2)));
	/*END: Calculation of GammaTauTau*/

	/*BEGIN:Calculation of GammaPiPi*/
	GammaPiPi = 0;
	for (i = 1; i < 35; i++)
		GammaPiPi
		    = GammaPiPi
		      + (n[i]) * (L[i]) * (L[i] - 1) * (pow((7.1 - Pressurevar), (L[i] - 2.))) * (pow((Tau - 1.222), (J[i])));
	/*END: Calculation of GammaPiPi*/

	/*************************Partial derivatives calculation*****************************************/
	/*************************************************************************************************/
	/*************************************************************************************************/

	/*BEGIN: Fluid isobaric heat capacity*/
	Cp = -(Tau * Tau) * (GammaTauTau) * (GazConst);
	/*END: Fluid isobaric heat capacity*/

	/*BEGIN: Fluid isochoric heat capacity*/
	/* Cv is not used currently 9.2003*/
	// WW Cv = (- (pow(Tau,2))* (GammaTauTau) + pow(GammaPi - Tau * (GammaPiTau),2) / GammaPiPi) * GazConst;
	/*BEGIN: Fluid isochoric heat capacity*/

	return Cp;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   01/2005 OK Implementation
   last modified:
**************************************************************************/
void MFPDelete()
{
	long i;
	int no_mfp = (int)mfp_vector.size();
	for (i = 0; i < no_mfp; i++)
		delete mfp_vector[i];
	mfp_vector.clear();
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   10/2005 OK/YD Implementation
**************************************************************************/
CFluidProperties* MFPGet(const string& name)
{
	CFluidProperties* m_mfp = NULL;
	for (int i = 0; i < (int)mfp_vector.size(); i++)
	{
		m_mfp = mfp_vector[i];
		if (m_mfp->name.compare(name) == 0)
			return m_mfp;
	}
	return NULL;
}

CFluidProperties* MFPGet(int fluid) // NB
{
	CFluidProperties* m_mfp = NULL;
	for (int i = 0; i < (int)mfp_vector.size(); i++)
	{
		m_mfp = mfp_vector[i];
		if (m_mfp->fluid_id == fluid)
			return m_mfp;
	}
	return NULL;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   10/2005 YD Calculate Enthalpy
**************************************************************************/
double CFluidProperties::CalcEnthalpy(double temperature)
{
	//----------------------------------------------------------------------
	CalPrimaryVariable(enthalpy_pcs_name_vector);
	//----------------------------------------------------------------------
	// double temperature = primary_variable[0];
	double val = 0.0;
	double T0_integrate = 0.0;
	double T1_integrate = 0.0;
	double heat_capacity_all;

	switch (GetHeatCapacityModel())
	{
		case 5:
			MFPGet("LIQUID");
			//------------PART 1--------------------------------
			T0_integrate = PhysicalConstant::CelsiusZeroInKelvin;
			T1_integrate = T_Latent1;
			int npoint = 100; // Gauss point
			double DT = (T1_integrate - T0_integrate) / npoint;
			for (int i = 0; i < npoint; i++)
			{
				temperature_buffer = T0_integrate + i * DT;
				heat_capacity_all = Fem_Ele_Std->FluidProp->SpecificHeatCapacity();
				temperature_buffer = T0_integrate + (i + 1) * DT;
				heat_capacity_all += Fem_Ele_Std->FluidProp->SpecificHeatCapacity();
				val += 0.5 * DT * heat_capacity_all;
			}

			//------------PART 2--------------------------------
			npoint = 500;
			T0_integrate = T_Latent1;
			T1_integrate = temperature;
			DT = (T1_integrate - T0_integrate) / npoint;
			for (int i = 0; i < npoint; i++)
			{
				temperature_buffer = T0_integrate + i * DT;
				heat_capacity_all = Fem_Ele_Std->FluidProp->SpecificHeatCapacity();
				temperature_buffer = T0_integrate + (i + 1) * DT;
				heat_capacity_all += Fem_Ele_Std->FluidProp->SpecificHeatCapacity();
				val += 0.5 * DT * heat_capacity_all;
			}
			break;
			//  case 5:
	}
	return val;
}

/**************************************************************************
   PCSLib-Method:
   08/2008 OK
   last change: 11/2008 NB
**************************************************************************/
double MFPGetNodeValue(long node, const string& mfp_name, int phase_number)
{
	CFluidProperties* m_mfp = mfp_vector[max(phase_number, 0)];
	const int restore_mode = m_mfp->mode;
	m_mfp->mode = 0;
	m_mfp->node = node;

	int mfp_id = -1;
	std::vector<std::string>* vec_var_names;
	switch (mfp_name[0])
	{
		case 'V':
			mfp_id = 0; // VISCOSITY
			vec_var_names = &m_mfp->viscosity_pcs_name_vector;
			break;
		case 'D':
			mfp_id = 1; // DENSITY
			vec_var_names = &m_mfp->density_pcs_name_vector;
			break;
		case 'H':
			mfp_id = 2; // HEAT_CONDUCTIVITY
			vec_var_names = &m_mfp->heat_conductivity_pcs_name_vector;
			break;
		case 'S':
			mfp_id = 3; // SPECIFIC HEAT CAPACITY
			vec_var_names = &m_mfp->specific_heat_capacity_pcs_name_vector;
			break;
		default:
			mfp_id = -1;
			static std::vector<std::string> default_var_names;
			if (default_var_names.empty())
			{
				default_var_names.push_back("PRESSURE1");
				default_var_names.push_back("TEMPERATURE1");
				default_var_names.push_back("CONCENTRATION1");
			}
			vec_var_names = &default_var_names;
			break;
	}

	// std::vector<double> arguments(vec_var_names->size());
	std::vector<double> arguments(std::max(size_t(3ul), vec_var_names->size()));
	for (unsigned i = 0; i < vec_var_names->size(); i++)
	{
		CRFProcess* pcs = PCSGet((*vec_var_names)[i], true);
		if (pcs)
		{
			int var_idx = pcs->GetNodeValueIndex((*vec_var_names)[i], true);
			if ((*vec_var_names)[i] == "PRESSURE1")
				arguments[0] = pcs->GetNodeValue(node, var_idx);
			else if ((*vec_var_names)[i] == "TEMPERATURE1")
				arguments[1] = pcs->GetNodeValue(node, var_idx);
			else if ((*vec_var_names)[i] == "CONCENTRATION1")
				arguments[2] = pcs->GetNodeValue(node, var_idx);
			else
				std::cout << "The variable " << (*vec_var_names)[i] << " is not supported in MFPGetNodeValue."
				          << std::endl;
		}
		else
		{
			arguments[i] = 0.0;
		}
	}

	if (m_mfp->cmpN > 0)
	{
		CRFProcess* m_pcs = PCSGet("MULTI_COMPONENTIAL_FLOW");
		arguments.resize(m_mfp->cmpN + 2);
		for (int PVIndex = 0; PVIndex < m_mfp->cmpN + 2; PVIndex++)
			arguments[PVIndex]
			    = m_pcs->GetNodeValue(node, m_pcs->GetNodeValueIndex(m_pcs->pcs_primary_function_name[PVIndex]));
	}

	//......................................................................
	double mfp_value = .0;
	switch (mfp_id)
	{
		case 0:
			mfp_value = m_mfp->Viscosity(arguments.data());
			break;
		// NB 4.8.01
		case 1:
			mfp_value = m_mfp->Density(arguments.data());
			break;
		case 2:
			mfp_value = m_mfp->HeatConductivity(arguments.data());
			break;
		// NB AUG 2009
		case 3:
			mfp_value = m_mfp->SpecificHeatCapacity(arguments.data());
			break;
		default:
			cout << "MFPGetNodeValue: no MFP data"
			     << "\n";
			break;
	}
	//......................................................................
	m_mfp->mode = restore_mode; // NB changeback
	return mfp_value;
}

/**************************************************************************
   Task: Derivative of density with respect to pressure at constant T and mass fraction
   Programing: 09/2009 NB
**************************************************************************/
double CFluidProperties::drhodP(double* variables)
{
	const double p = variables[0];
	const double T = variables[1];

	double arguments[2];
	double rho1, rho2, drhodP = 0.0;

	if (p < 0)
		return 0;

	switch (compressibility_model_pressure)
	{
		case 0: // incompressible
			drhodP = 0;
			break;

		case 1: // constant slope compressibility
			drhodP = compressibility_pressure;
			break;

		case 2: // use of fct file
			drhodP = 0; // to be done
			break;
		case 3: // use of difference quotient
			arguments[1] = T;
			// in case 3, compressibility_pressure acts as delta P
			arguments[0] = p + (compressibility_pressure / 2.);
			rho1 = Density(arguments);

			arguments[0] = p - (compressibility_pressure / 2.);
			rho2 = Density(arguments);

			// drhodP = (rho(p+dP/2)-rho(P-dP/2))/dP
			// in case 3, compressibility_pressure acts as delta P
			drhodP = (rho1 - rho2) / compressibility_pressure;

			break; // use of difference quotient
		case 7: // use of fct file
			drhodP = 1.0 / p; // to be done
			break;

		case 15: // volume translated Peng-Robinson
			if (eos_name == "VTPR" || eos_name == "PR" || eos_name == "IDEAL")
			{
				CRFProcess* m_pcs = PCSGet("MULTI_COMPONENTIAL_FLOW");
				std::vector<double> roots;

				double beta_m = 0.0;
				double v_m = 0.0;
				const double R = 1000.0 * PhysicalConstant::IdealGasConstant;

				for (int CIndex = 2; CIndex < cmpN + 2; CIndex++)
				{
					double A, B, c, C, dvdp, fctA, fctB, fctC, beta, Tr, z, z1, z2, z3, v, Trr, fct;
					therm_prop(m_pcs->pcs_primary_function_name[CIndex]);
					Tr = T / Tc;
					// double a0 = 1.0 + m0*(1.0 - Tr) + n0*(1.0 - Tr)*(0.7 - Tr);
					Trr = 1.0 - pow(Tr, 0.6667);
					fct = -0.25 + pc / p;
					if (fluid_id == 1)
						fct = (1.0780 + p / 6.8e8);
					c = fct * (k1 + k2 * Trr + k3 * Trr * Trr) * R * Tc / pc;
					if (eos_name == "PR")
						c = 0.0;
					A = a * p * pow(R * T, -2.0);
					B = b * p * pow(R * T, -1.0);
					C = c * p * pow(R * T, -1.0);
					z1 = B - 1.0 + 3.0 * C;
					z2 = (A - 3.0 * B * B + 3.0 * C * C + 2.0 * B * C - 2.0 * B - 2.0 * C);
					z3 = (B * B * B + C * C * C + B * B - C * C + B * C * C - 3.0 * B * B * C - 2.0 * B * C + C * A
					      - A * B);
					NsPol3(z1, z2, z3, &roots);
					z = FindMax(roots);
					if (fluid_id == 1)
						z = FindMin(roots);
					if (eos_name == "IDEAL")
					{
						z = 1.0;
						c = 0.0;
						b = 0.0;
						a = 0.0;
					}
					v = (z * R * T / p);
					fctB = (v + c) * (v + c + b) + b * (v + c - b);
					fctA = v + c - b;
					fctC = v + c + b;
					dvdp = fctA * fctB / (p * fctB + 2.0 * p * fctA * fctC - 2.0 * R * T * fctC + a);
					beta = dvdp;
					beta_m += beta;
					v_m += v;
				}

				drhodP = beta_m / v_m; // beta_t=beta_g*(vg/vt)+beta_l*(vl/vt)
			}
			else if (eos_name == "CONSTANT")
			{
				for (int CIndex = 2; CIndex < cmpN + 2; CIndex++)
				{
					drhodP = variables[CIndex - 2] * beta_p[CIndex - 2];
				}
			}

			break;

		default:
			drhodP = drho_dp;
	}

	return drhodP;
}

/**************************************************************************
   Task: derivative of density with respect to temperature at constant p and mass fraction
   Programing: 09/2009 NB
**************************************************************************/
double CFluidProperties::drhodT(double* variables)
{
	double a0, A, B, c, C, da0, da, alpha, alpha_m, dvdT, fctA, fctB, fctC, p, R, T, Tr, z, z1, z2, z3, v, v_m, Trr,
	    fct;
	std::vector<double> roots;
	CRFProcess* m_pcs;
	m_pcs = PCSGet("MULTI_COMPONENTIAL_FLOW");
	p = variables[0];
	T = variables[1];
	R = 1000.0 * PhysicalConstant::IdealGasConstant;
	alpha_m = 0.0;
	v_m = 0.0;
	double arguments[2];
	double rho1, rho2, drhodT = 0.0;

	if (!drho_dT_unsaturated) // fluid expansion (drho/dT) for unsaturated case activated?
	{
		if (p < 0)
			return 0;
	}

	switch (compressibility_model_temperature)
	{
		case 0: // fluid is incompressible
			drhodT = 0;
			break;
		case 1: // constant slope compressibility, (for test cases)
			drhodT = compressibility_temperature;
			break;

		case 2: // use of fct file
			drhodT = 0; // to be done
			break;

		case 3: // use of difference quotient
			// in case 3, compressibility_temperature acts as delta T
			arguments[0] = p;

			arguments[1] = T + (compressibility_temperature / 2.);
			rho1 = Density(arguments);

			arguments[1] = T - (compressibility_temperature / 2.);
			rho2 = Density(arguments);

			drhodT = (rho1 - rho2) / compressibility_temperature;
			break;
		case 6: // rho(p,T) = rho_0*(1+beta_p*(p-p_0)+beta_T*(T-T_0))
			return drho_dT;
		case 7: // use of fct file
			drhodT = 1.0 / T;
			break;

		case 15: // volume translated Peng-Robinson
			if (eos_name == "VTPR" || eos_name == "PR" || eos_name == "IDEAL")
			{
				for (int CIndex = 2; CIndex < cmpN + 2; CIndex++)
				{
					therm_prop(m_pcs->pcs_primary_function_name[CIndex]);
					Tr = T / Tc;
					a0 = 1.0 + m0 * (1.0 - Tr) + n0 * (1.0 - Tr) * (0.7 - Tr);
					da0 = -(m0 + n0 * (1.0 - Tr) + n0 * (0.7 - Tr)) / Tc;
					da = 0.457235 * 2.0 * a0 * da0 * pow(R * Tc, 2.0) / pc;
					Trr = 1.0 - pow(Tr, 0.6667);
					fct = -0.25 + pc / p;
					if (fluid_id == 1)
						fct = (1.0780 + p / 6.8e8);
					c = fct * (k1 + k2 * Trr + k3 * Trr * Trr) * R * Tc / pc;
					if (eos_name == "PR")
						c = 0.0;
					A = a * p * pow(R * T, -2.0);
					B = b * p * pow(R * T, -1.0);
					C = c * p * pow(R * T, -1.0);
					z1 = B - 1.0 + 3.0 * C;
					z2 = (A - 3.0 * B * B + 3.0 * C * C + 2.0 * B * C - 2.0 * B - 2.0 * C);
					z3 = (B * B * B + C * C * C + B * B - C * C + B * C * C - 3.0 * B * B * C - 2.0 * B * C + C * A
					      - A * B);
					NsPol3(z1, z2, z3, &roots);
					z = FindMax(roots);
					if (fluid_id == 1)
						z = FindMin(roots);
					if (eos_name == "IDEAL")
					{
						z = 1.0;
						c = 0.0;
						b = 0.0;
						a = 0.0, da = 0.0;
					}
					v = (z * R * T / p);
					fctB = (v + c) * (v + c + b) + b * (v + c - b);
					fctA = v + c - b;
					fctC = v + c + b;
					dvdT = (R * fctB - da * fctA) / (p * fctB + 2.0 * p * fctA * fctC - 2.0 * R * T * fctC + a);
					alpha = dvdT;
					alpha_m += alpha;
					v_m += v;
				}
				drhodT = -alpha_m / v_m; // alpha_t=alpha_g*(vg/vt)+alpha_l*(vl/vt)
			}

			if (eos_name == "CONSTANT")
				for (int CIndex = 2; CIndex < cmpN + 2; CIndex++)
				{
					drhodT = variables[CIndex - 2] * alpha_T[CIndex - 2];
				}

			break;

		default:
			drhodT = 0;
			break;
	}
	return drhodT;
}

/**************************************************************************
   Task: derivative of density with respect to mass fraction at constant p and T
   Programing: 09/2012 AKS
**************************************************************************/
double CFluidProperties::drhodX(int CIndex, double* variables)
{
	static double drhodX;
	switch (solutal_expansivity_model)
	{
		case 1:
			drhodX = solutal_expansivity;
			break;

		case 15:
			if (eos_name == "CONSTANT")
				drhodX = -Density(variables) / rho[CIndex - 2];
			else
				drhodX = -Density(variables) / ComponentDensity(CIndex, variables);
			break;

		default:
			drhodX = 0.0;
			break;
	}
	return drhodX;
}

/**************************************************************************
   Task: Calculation of component density based on VTPR, PR and IDEAL EoS
   Programing: 09/2012 AKS
**************************************************************************/
double CFluidProperties::ComponentDensity(int CIndex, double* variables)
{
	double A, B, z, z1, z2, z3, Tr, c, C, Trr, T, p, density, R, fct;
	CRFProcess* m_pcs;
	m_pcs = PCSGet("MULTI_COMPONENTIAL_FLOW");
	std::vector<double> roots;
	p = variables[0];
	T = variables[1];
	R = 1000.0 * PhysicalConstant::IdealGasConstant;

	if (eos_name == "VTPR" || eos_name == "PR" || eos_name == "IDEAL")
	{
		therm_prop(m_pcs->pcs_primary_function_name[CIndex]);
		Tr = T / Tc;
		// double a0 = 1.0 + m0*(1.0 - Tr) + n0*(1.0 - Tr)*(0.7 - Tr);
		Trr = 1.0 - pow(Tr, 0.6667);
		fct = -0.25 + pc / p;
		if (fluid_id == 1)
			fct = (1.0780 + p / 6.8e8);
		c = fct * (k1 + k2 * Trr + k3 * Trr * Trr) * R * Tc / pc;
		if (eos_name == "PR")
			c = 0.0;
		A = a * p * pow(R * T, -2.0);
		B = b * p * pow(R * T, -1.0);
		C = c * p * pow(R * T, -1.0);
		z1 = B - 1.0 + 3.0 * C;
		z2 = (A - 3.0 * B * B + 3.0 * C * C + 2.0 * B * C - 2.0 * B - 2.0 * C);
		z3 = (B * B * B + C * C * C + B * B - C * C + B * C * C - 3.0 * B * B * C - 2.0 * B * C + C * A - A * B);
		NsPol3(z1, z2, z3, &roots);
		z = FindMax(roots);
		if (fluid_id == 1)
			z = FindMin(roots);
		if (eos_name == "IDEAL")
			z = 1.0;
		density = p * molar_mass / (z * R * T);
	}
	else if (eos_name == "CONSTANT")
	{
		density = rho[CIndex - 2];
	}
	else
	{
		density = 0.0;
	}

	return density;
}
/**************************************************************************
   Task: Calculation of component's effective diffusion coefficient in a mixture
   Programing: 05/2010 AKS
**************************************************************************/
double CFluidProperties::EffectiveDiffusionCoef(int CIndex, double* variables)
{
	CRFProcess* m_pcs;
	m_pcs = PCSGet("MULTI_COMPONENTIAL_FLOW");
	double MI = 0.0, VdI = 0.0, effective_diffusion_coef = 0.0, T, p;
	T = variables[1];
	p = variables[0];

	switch (diffusion_model)
	{
		case 1:
			effective_diffusion_coef = diffusion_coef;
			break;

		case 15:
			if (eos_name == "VTPR" || eos_name == "PR" || eos_name == "IDEAL")
			{
				for (int i = 0; i < cmpN; i++)
				{
					therm_prop(m_pcs->pcs_primary_function_name[i + 2]);
					MI += 1.0 / molar_mass;
					VdI += pow(Vd, 1.0 / 3.0);
				}
				therm_prop(m_pcs->pcs_primary_function_name[CIndex]);
				MI = MI - 1.0 / molar_mass;
				VdI = VdI - pow(Vd, 1.0 / 3.0);
				effective_diffusion_coef = 0.0143 * pow(T, 1.75) / (p * pow(2.0 * pow(MI, -1), 0.5) * pow(VdI, 2));
			}
			if (eos_name == "CONSTANT")
				effective_diffusion_coef = D0[CIndex - 2];
			break;
	}
	return effective_diffusion_coef;
}
#ifdef MFP_TEST
//-----------------------------------------------------
//
/*!
   \brief constructor of class Hash_Table

    Read data and fill the member data

    WW 07.2011
 */
Hash_Table::Hash_Table(string f_name)
{
	string aline;
	std::stringstream ss;
	int i_buff = 0.;
	int i, length;

	f_name = FilePath + f_name;
	ifstream ins(f_name.c_str());
	if (!ins.good())
	{
		cout << "File " << f_name << " cannot openned. Program exit ";
		exit(1);
	}

	getline(ins, aline);
	if (aline.find("Varaiable_Number") != string::npos)
	{
		ss.str(aline);
		ss >> aline >> num_var;
		ss.clear();
	}
	else
	{
		cout << "Number of varaiables are not defined. Program exit ";
		exit(1);
	}

	num_par = 2;
	num_var -= 2;
	length = num_par + num_var;

	names.resize(length);
	getline(ins, aline);
	ss.str(aline);
	for (i = 0; i < length; i++)
		ss >> names[i];
	ss.clear();

	int counter = 0;
	double par = 0.;
	table_section_ends.push_back(0);
	while (!ins.eof())
	{
		getline(ins, aline);
		if (aline.find("...") != string::npos)
		{
			table_section_ends.push_back(counter + 1);
			hash_row_par.push_back(par);
			continue;
		}
		if (aline.find("---") != string::npos)
			break;

		double* data_vp = new double[length - 1];
		hash_table_data.push_back(data_vp);

		ss.str(aline);
		ss >> par;
		for (i = 0; i < length - 1; i++)
			ss >> data_vp[i];
		ss.clear();

		counter++;
	}
}
/*!
   \brief desconstructor of class Hash_Table


    WW 07.2011
 */
Hash_Table::~Hash_Table()
{
	while (hash_table_data.size() > 0)
	{
		delete[] hash_table_data[hash_table_data.size() - 1];
		hash_table_data.pop_back();
	}
}

/*!
   \brief desconstructor of class Hash_Table


    WW 07.2011
 */
double Hash_Table::CalcValue(double* var, const int var_id) const
{
	size_t i, j, k;
	double *data_0, *data_1;
	double val_e[4];

	if (var[0] < hash_row_par[0])
		var[0] = hash_row_par[0];
	if (var[0] > hash_row_par[hash_row_par.size() - 1])
		var[0] = hash_row_par[hash_row_par.size() - 1];

	for (i = 0; i < hash_row_par.size() - 1; i++)
		if ((var[0] >= hash_row_par[i]) && (var[0] < hash_row_par[i + 1]))
		{
			for (j = 0; j < 2; j++)
			{
				data_0 = hash_table_data[i + j];
				for (k = table_section_ends[i + j]; k < table_section_ends[i + j + 1]; k++)
				{
				}
			}
			break;
		}
	return 0.;
}
#endif //#ifdef MFP_TEST
