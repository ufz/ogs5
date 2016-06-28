/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   FEMLib - Object: MSP Solid Properties
   Task:
   Programing:
   08/2004 WW Implementation
   last modified:
**************************************************************************/

// C++ STL
//#include <math.h>
//#include <string>
//#include <fstream>
//#include <iostream>
//#include <sstream>
#include <cfloat>

// FEM-Makros
#include "makros.h"
#include "rf_pcs.h"

// Time
#include "fem_ele_std.h"
#include "fem_ele_vec.h"
#include "rf_msp_new.h"
#include "rf_tim_new.h"
//#include "rf_mmp_new.h"
#include "pcs_dm.h"

#include "StringTools.h"
#include "files0.h" // GetLineFromFile1
#include "tools.h" // GetLineFromFile
#include "PhysicalConstant.h"

#include "Eigen/Dense"
#include "minkley.h"
#include "burgers.h"

std::vector<SolidProp::CSolidProperties*> msp_vector;
std::vector<std::string> msp_key_word_vector; // OK

namespace SolidProp
{
using namespace std;
using FiniteElement::ElementValue_DM;
using Math_Group::Matrix;

/**************************************************************************
   FEMLib-Method:
   Task: OBJ read function
   Programing:
   08/2004 OK Implementation for fuild properties
   08/2004 WW Modification for solid properties
   12/2005 WW Creep properties
**************************************************************************/
std::ios::pos_type CSolidProperties::Read(std::ifstream* msp_file)
{
	char buffer[MAX_ZEILE];
	std::string sub_line;
	std::string line_string;
	std::string delimiter(" ");
	bool new_keyword = false;
	std::string hash("#");
	std::ios::pos_type position;
	//  ios::pos_type position0;
	std::string sub_string;
	// WW bool new_subkeyword = false;
	std::string dollar("$");
	std::string delimiter_type(":");
	std::stringstream in_sd;

	int i = 0, Size = 0;

	//========================================================================
	// Schleife ueber alle Phasen bzw. Komponenten
	while (!new_keyword)
	{
		// WW new_subkeyword = false;

		position = msp_file->tellg();
		if (!GetLineFromFile(buffer, msp_file))
			break;
		line_string = buffer;
		trim(line_string); // NW
		trim(line_string, ':'); // NW
		if (line_string.find(hash) != string::npos)
		{
			new_keyword = true;
			break;
		}
		//....................................................................
		// NAME //OK
		// subkeyword found
		if (line_string.find("$NAME") != string::npos)
		{
			in_sd.str(GetLineFromFile1(msp_file));
			in_sd >> name; // sub_line
			in_sd.clear();
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$SWELLING_PRESSURE_TYPE") != string::npos)
		{
			in_sd.str(GetLineFromFile1(msp_file));
			in_sd >> SwellingPressureType;
			if (SwellingPressureType == 1 || SwellingPressureType == 2)
			{
				in_sd >> Max_SwellingPressure;
				in_sd.clear();
			}
			// 10.03.2008 WW
			else if (SwellingPressureType == 3 || SwellingPressureType == 4)
			{
				if (!PCSGet("MULTI_PHASE_FLOW") || !PCSGet("RICHARDS_FLOW"))
				{
					data_Youngs = new Matrix(9);
					//! 0: \f$ \kappa_i0     \f$
					//! 1: \f$ \alpha_i     \f$
					//! 2: \f$ \kappa_{s0}  \f$
					//! 3: \f$ \alpha_{sp}  \f$
					//! 4: \f$ \alpha_{ss}  \f$
					//! 5: \f$ p_ref        \f$
					//! 6: \f$ buffer       \f$
					if (SwellingPressureType == 3)
						in_sd >> (*data_Youngs)(0) >> (*data_Youngs)(1) >> (*data_Youngs)(2) >> (*data_Youngs)(3)
						    >> (*data_Youngs)(4) >> (*data_Youngs)(5);
					else if (SwellingPressureType == 4)
						in_sd >> (*data_Youngs)(0) >> (*data_Youngs)(1) >> (*data_Youngs)(2);
					in_sd.clear();
				}
				else
				{
					std::cout << "No multi-phase flow coupled. The thermal elatic model can only be used in H2 coupled "
					             "proccess."
					          << "\n";
					std::cout << "Quit the simulation now!"
					          << "\n";
					abort();
				}
			}
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$DENSITY") != string::npos)
		{
			in_sd.str(GetLineFromFile1(msp_file));
			in_sd >> Density_mode;
			if (Density_mode == 0) // rho = f(x)
			{
				in_sd >> Size;
				in_sd.clear();
				data_Density = new Matrix(Size, 2);
				for (i = 0; i < Size; i++)
				{
					in_sd.str(GetLineFromFile1(msp_file));
					in_sd >> (*data_Density)(i, 0) >> (*data_Density)(i, 1);
					in_sd.clear();
				}
			}
			else if (Density_mode == 1) // rho = const
			{
				data_Density = new Matrix(1);
				in_sd >> (*data_Density)(0);
				in_sd.clear();
			}
		}
		//....................................................................
		if (line_string.find("$THERMAL") != string::npos)
		{
			in_sd.str(GetLineFromFile1(msp_file));
			in_sd >> sub_line >> ThermalExpansion;
			in_sd.clear();
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("CAPACITY") != string::npos)
		{
			in_sd.str(GetLineFromFile1(msp_file));
			in_sd >> Capacity_mode;
			switch (Capacity_mode)
			{
				case 0: //  = f(x)
					in_sd >> Size;
					in_sd.clear();
					data_Capacity = new Matrix(Size, 2);
					for (i = 0; i < Size; i++)
					{
						in_sd.str(GetLineFromFile1(msp_file));
						in_sd >> (*data_Capacity)(i, 0) >> (*data_Capacity)(i, 1);
						in_sd.clear();
					}
					break;
				case 1: //  = const
					data_Capacity = new Matrix(1);
					in_sd >> (*data_Capacity)(0);
					in_sd.clear();
					break;
				case 2: // boiling model for rock. WW
					// 0. Wet capacity
					// 1. Dry capacity
					// 2. Boiling temperature
					// 3. Boiling temperature range
					// 4. Latent of vaporization
					data_Capacity = new Matrix(5);
					for (i = 0; i < 5; i++)
						in_sd >> (*data_Capacity)(i);
					in_sd.clear();
					break;
				case 3: // DECOVALEX THM1, Bentonite
					in_sd.clear();
					break;
				case 4:
					// 0. Capacity at density 1
					// 1. Capacity at density 2
					// 2. density 1
					// 3. density 2
					data_Capacity = new Matrix(5);
					for (i = 0; i < 4; i++)
						in_sd >> (*data_Capacity)(i);
					in_sd.clear();
					break;
				// TES
				case 5: // Capacity depending on solid conversion
					// 0. Capacity at lower_density_limit (reactive system property)
					// 1. Capacity at upper_density_limit (reactive system property)
					data_Capacity = new Matrix(3);
					for (i = 0; i < 2; i++)
						in_sd >> (*data_Capacity)(i);
					in_sd.clear();
					break;
				case 6: // Capacity depending on loading with adsorbate
					// 0. Capacity at desorbed state (lower density limit)
					// 1. Capacity of adsorbate
					data_Capacity = new Matrix(3);
					for (i = 0; i < 2; i++)
						in_sd >> (*data_Capacity)(i);
					in_sd.clear();
					break;
				case 7: //Capacity depending on solid conversion (density-based average)
					//0. Capacity at lower_density_limit (reactive system property)
					//1. Capacity at upper_density_limit (reactive system property)
					data_Capacity = new Matrix(3);
					for(i=0; i<2; i++)
						in_sd>> (*data_Capacity)(i);
					in_sd.clear();
					break;
			}
		}

		//....................................................................
		// subkeyword found
		if (line_string.compare("CONDUCTIVITY") == 0)
		{
			in_sd.str(GetLineFromFile1(msp_file));
			in_sd >> Conductivity_mode;
			switch (Conductivity_mode)
			{
				case 0: //  = f(T) //21.12.2009 WW
					in_sd >> Size;
					in_sd.clear();
					data_Conductivity = new Matrix(Size, 2);
					for (i = 0; i < Size; i++)
					{
						in_sd.str(GetLineFromFile1(msp_file));
						in_sd >> (*data_Conductivity)(i, 0) >> (*data_Conductivity)(i, 1);
						in_sd.clear();
					}
					// WW
					conductivity_pcs_name_vector.push_back("TEMPERATURE1");
					break;
				case 1: //  = const
					data_Conductivity = new Matrix(1);
					in_sd >> (*data_Conductivity)(0);
					in_sd.clear();
					break;
				case 2: // boiling model for rock. WW
					// 0. Wet conductivity
					// 1. Dry conductivity
					// 2. Boiling temperature
					// 3. Boiling temperature range
					data_Conductivity = new Matrix(4);
					for (i = 0; i < 4; i++)
						in_sd >> (*data_Conductivity)(i);
					in_sd.clear();
					capacity_pcs_name_vector.push_back("TEMPERATURE1");
					capacity_pcs_name_vector.push_back("SATURATION1");
					break;
				case 3: // DECOVALEX THM1, Bentonite
					in_sd.clear();
					capacity_pcs_name_vector.push_back("TEMPERATURE1");
					capacity_pcs_name_vector.push_back("SATURATION1");
					break;
				case 30: // another model for bentonite. WW
					// 0. maximum conductivity
					// 1. minimum conductivity
					// 2. saturation
					data_Conductivity = new Matrix(3);
					for (i = 0; i < 3; i++)
						in_sd >> (*data_Conductivity)(i);
					in_sd.clear();
					capacity_pcs_name_vector.push_back("SATURATION1");
					break;
				case 4: //  = f(S) //21.12.2009 WW
					in_sd >> Size;
					in_sd.clear();
					data_Conductivity = new Matrix(Size, 2);
					for (i = 0; i < Size; i++)
					{
						in_sd.str(GetLineFromFile1(msp_file));
						in_sd >> (*data_Conductivity)(i, 0) >> (*data_Conductivity)(i, 1);
						in_sd.clear();
					}
					break;
				case 5: // DECOVALEX2015, Task B2, Buffer: f(S,T) by matrix function
					in_sd >> T_0;
					in_sd.clear();
					conductivity_pcs_name_vector.push_back("TEMPERATURE1");
					conductivity_pcs_name_vector.push_back("SATURATION1");
					break;
			}
			in_sd.clear();
		}

		//....................................................................
		if (line_string.find("CONDUCTIVITY_TENSOR") != string::npos) // subkeyword found
		{
			in_sd.str(GetLineFromFile1(msp_file));
			in_sd >> thermal_conductivity_tensor_type_name;
			thermal_conductivity_tensor_dim = 0; // NW
			switch (thermal_conductivity_tensor_type_name[0])
			{
				case 'I': // isotropic
					thermal_conductivity_tensor_type = 0;
					in_sd >> thermal_conductivity_tensor[0];
					thermal_conductivity_tensor[1] = thermal_conductivity_tensor[2] = thermal_conductivity_tensor[0];
					break;
				case 'O': // orthotropic
					thermal_conductivity_tensor_type = 1;
					in_sd >> thermal_conductivity_tensor_dim;
					if (thermal_conductivity_tensor_dim == 0)
						std::cout << "Error in CSolidProperties::Read: no tensor dimension"
						          << "\n";
					if (thermal_conductivity_tensor_dim == 2)
					{
						in_sd >> thermal_conductivity_tensor[0];
						in_sd >> thermal_conductivity_tensor[1];
					}
					if (thermal_conductivity_tensor_dim == 3)
					{
						in_sd >> thermal_conductivity_tensor[0];
						in_sd >> thermal_conductivity_tensor[1];
						in_sd >> thermal_conductivity_tensor[2];
					}
					break;
				case 'A': // anisotropic
					thermal_conductivity_tensor_type = 2;
					in_sd >> thermal_conductivity_tensor_dim;
					if (thermal_conductivity_tensor_dim == 0)
						std::cout << "Error in CSolidProperties::Read: no tensor dimension"
						          << "\n";
					if (thermal_conductivity_tensor_dim == 2)
					{
						in_sd >> thermal_conductivity_tensor[0];
						in_sd >> thermal_conductivity_tensor[1];
						in_sd >> thermal_conductivity_tensor[2];
						in_sd >> thermal_conductivity_tensor[3];
					}
					if (thermal_conductivity_tensor_dim == 3)
					{
						in_sd >> thermal_conductivity_tensor[0];
						in_sd >> thermal_conductivity_tensor[1];
						in_sd >> thermal_conductivity_tensor[2];
						in_sd >> thermal_conductivity_tensor[3];
						in_sd >> thermal_conductivity_tensor[4];
						in_sd >> thermal_conductivity_tensor[5];
						in_sd >> thermal_conductivity_tensor[6];
						in_sd >> thermal_conductivity_tensor[7];
						in_sd >> thermal_conductivity_tensor[8];
					}
					break;
				default:
					cout << "Error in CSolidProperties::Read: no valid thermal conductivity tensor type"
					     << "\n";
					break;
			}
			in_sd.clear();
		}

		//....................................................................
		// subkeyword found
		if (line_string.find("$ELASTICITY") != string::npos)
		{
			in_sd.str(GetLineFromFile1(msp_file));
			in_sd >> sub_line >> PoissonRatio;
			in_sd.clear();
		}
		//....................................................................
		// 12.2009. WW
		if (line_string.find("$EXCAVATION") != string::npos)
		{
			in_sd.str(GetLineFromFile1(msp_file));
			in_sd >> excavation;
			in_sd.clear();
		}
		//....................................................................
		if (line_string.find("YOUNGS_MODULUS") != string::npos)
		{ // subkeyword found
			in_sd.str(GetLineFromFile1(msp_file));
			in_sd >> Youngs_mode;
			int type = Youngs_mode;
			if (type > 9 && type < 14)
				type = 1000;
			switch (type) // 15.03.2008 WW
			{
				case 0: //  = f(x)
					in_sd >> Size;
					in_sd.clear();
					data_Youngs = new Matrix(Size, 2);
					for (i = 0; i < Size; i++)
					{
						in_sd.str(GetLineFromFile1(msp_file));
						in_sd >> (*data_Youngs)(i, 0) >> (*data_Youngs)(i, 1);
						in_sd.clear();
					}
					break;
				case 1: //  = const
					data_Youngs = new Matrix(1);
					in_sd >> (*data_Youngs)(0);
					in_sd.clear();
					break; // UJG 24.11.2009
				case 2: //  = const
					// data_Youngs Lubby1 model
					//  0: E_0
					//  1: a (factor)
					//  2: n (exponent)
					data_Youngs = new Matrix(3);
					in_sd >> (*data_Youngs)(0) >> (*data_Youngs)(1) >> (*data_Youngs)(2);
					in_sd.clear();
					break;
				case 1000: // case 10-13: transverse isotropic linear elasticity (UJG 24.11.2009)
					// data_Youngs transverse isotropic linear elasticity
					//  0: E_i     (Young's modulus of the plane of isotropy)
					//  1: E_a     (Young's modulus w.r.t. the anisotropy direction)
					//  2: nu_{ia} (Poisson's ratio w.r.t. the anisotropy direction)
					//  3: G_a     (shear modulus w.r.t. the anisotropy direction)
					//  4: n_x     (x-coefficient of the local axis of anisotropy (2D case: -\sin\phi))
					//  5: n_y     (y-coefficient of the local axis of anisotropy (2D case: \cos\phi))
					//  6: n_z     (z-coefficient of the local axis of anisotropy (2D case: 0))
					data_Youngs = new Matrix(7);
					in_sd >> (*data_Youngs)(0) >> (*data_Youngs)(1) >> (*data_Youngs)(2) >> (*data_Youngs)(3)
					    >> (*data_Youngs)(4) >> (*data_Youngs)(5) >> (*data_Youngs)(6);
					in_sd.clear();
					break;
			}
		}
		//....................................................................
		// WX:06.2012 E_Function
		if (line_string.find("$E_Function") != string::npos)
		{
			in_sd.str(GetLineFromFile1(msp_file));
			in_sd >> E_Function_Model >> E_Function_Model_Value[0];
			in_sd.clear();
		}
		// Time dependent YOUNGS POISSON
		if (line_string.find("$TIME_DEPENDENT_YOUNGS_POISSON") != string::npos)
		{
			in_sd.str(GetLineFromFile1(msp_file));
			in_sd >> Time_Dependent_E_nv_mode;
			switch (Time_Dependent_E_nv_mode)
			{
				case 1: // isotropic
					in_sd >> Time_Dependent_E_nv_value[0] >> Time_Dependent_E_nv_value[1]; //[0] for E, [1] for nv
					break;
				case 2: // transversely isotropic
					in_sd >> Time_Dependent_E_nv_value[0] >> Time_Dependent_E_nv_value[1]
					    >> Time_Dependent_E_nv_value[2] >> Time_Dependent_E_nv_value[3] >> Time_Dependent_E_nv_value[4];
					//[0]: Ei, [1]: Ea, [2]: ni, [3]: Eia, [4]: Ga
					break;
				default:
					cout << "!ERROR in msp file, no valid TIME_DEPENDENT_YOUNG_POISSON mode!" << endl;
					break;
			}
			in_sd.clear();
		}
		//....................................................................
		if (line_string.find("$CREEP") != string::npos)
		{
			if (line_string.find("NORTON") != string::npos)
			{
				Creep_mode = 1;
				/*! \subsection Norton creep model */
				/*! \f$\dot\epsilon_s=A \left(\frac{\sigma_v}{\sigma^\ast}\right)^n\f$ */
				// data_Creep:
				//  0: A,   coefficient
				//  1: n,   exponential
				data_Creep = new Matrix(2);
				in_sd.str(GetLineFromFile1(msp_file));
				in_sd >> (*data_Creep)(0);
				in_sd >> (*data_Creep)(1);
				in_sd.clear();
			}
			if (line_string.find("BGRA") != string::npos)
			{
				Creep_mode = 2;
				/*! \subsection Temperature dependent creep model by BGR */
				/*! \f$\dot\epsilon_s=A\exp^{-Q/RT}\left(\frac{\sigma_v}{\sigma^\ast}\right)^n\f$ */
				// data_Creep:
				//  0: A,   coefficient
				//  1: n,   exponential
				//  2: Q,   activation energy
				data_Creep = new Matrix(3);
				in_sd.str(GetLineFromFile1(msp_file));
				in_sd >> (*data_Creep)(0);
				in_sd >> (*data_Creep)(1);
				in_sd >> (*data_Creep)(2);
				in_sd.clear();
			}
			// TN..................................................................
			if (line_string.find("BGRB") != string::npos)
			{
				Creep_mode = 3;
				// data_Creep:
				//  0: A1,  coefficient A1
				//  1: n1,   exponential n1
				//  2: Q1,   activation energy Q1
				//  3: A2,  coefficient A2
				//  4: n2,  coefficient n2
				//  5: Q2,  activation energy Q2
				data_Creep = new Matrix(6);
				in_sd.str(GetLineFromFile1(msp_file));
				in_sd >> (*data_Creep)(0);
				in_sd >> (*data_Creep)(1);
				in_sd >> (*data_Creep)(2);
				in_sd >> (*data_Creep)(3);
				in_sd >> (*data_Creep)(4);
				in_sd >> (*data_Creep)(5);
				in_sd.clear();
			}
			if (line_string.find("BGRSF") != string::npos)
			{
				Creep_mode = 4;
				// data_Creep:
				//  0: A,   coefficient A
				//  1: n,   exponential
				//  2: Q,   activation energy
				//  3: C,   coefficient C
				data_Creep = new Matrix(4);
				in_sd.str(GetLineFromFile1(msp_file));
				in_sd >> (*data_Creep)(0);
				in_sd >> (*data_Creep)(1);
				in_sd >> (*data_Creep)(2);
				in_sd >> (*data_Creep)(3);
				in_sd.clear();
			}
			//....................................................................
			if (line_string.find("LUBBY2") != string::npos)
			{
				Creep_mode = 1000;
				// data_Creep:
				//  0: eta_m
				//  1: m
				//  2: l
				//  3: eta_k
				//  4: k1
				//  5: k2
				//  6: G_k
				data_Creep = new Matrix(7, 2);
				in_sd.str(GetLineFromFile1(msp_file));
				for (i = 0; i < 7; i++)
					in_sd >> (*data_Creep)(i, 0);
				in_sd.clear();
			}
			if(line_string.find("BURGERS") != string::npos)
			{
				Creep_mode = 1001;
				// data_Creep:
				//  0: G_K0
				//  1: m_K
				//  2: mu_K
				//  3: m_vK
				//  4: G_M0
				//  5: K_M0
				//  6: mu_M
				//  7: m_vM
				//  8: m_GM // slope of elesticity temperature dependence of G_M
				//  9: m_KM // slope of elesticity temperature dependence of K_M
				// 10: T_ref // reference temperature;
				// 11: B // constant factor for Arrhenius term
				// 12: Q // activation energy in Arrhenius term
				data_Creep = new Matrix(14);
				in_sd.str(GetLineFromFile1(msp_file));
				for(i = 0; i < 14; i++)
					in_sd >> (*data_Creep)(i);
				in_sd.clear();

				//Local Newton scheme for Burgers model. TN 06.06.2014
				//initialised fully 3D
				material_burgers = new Burgers::SolidBurgers(data_Creep);

				SolidMath::InitialiseProjectionTensors();
			}
			if(line_string.find("MINKLEY") != string::npos)
			{
				Creep_mode = 1002;
				// data_Creep:
				//  0: Kelvin shear modulus
				//  1: Kelvin shear viscosity
				//  2: Maxwell shear modulus
				//  3: Maxwell bulk modulus
				//  4: Maxwell shear viscosity
				//  5: m -- effective stress factor in viscosity relationship
				//  6: n -- effective stress exponent in viscosity relationship
				//  7: initial cohesion
				//  8: hardening/softening modulus
				//  9: friction angle
				// 10: dilatancy angle
				// 11: transition angle for corner smoothing
				// 12: viscosity for viscoplastic regularisation
				// 13: l_0 // for temperature dependency of maxwell viscosity

				data_Creep = new Matrix(15);
				in_sd.str(GetLineFromFile1(msp_file));
				for(i = 0; i < 16; i++)
					in_sd >> (*data_Creep)(i);
				in_sd.clear();

				//Local Newton scheme for Burgers model. TN 06.06.2014
				//initialised fully 3D
				material_minkley = new Minkley::SolidMinkley(data_Creep);

				SolidMath::InitialiseProjectionTensors();
			}
		}
		// WX:10.2012, threshold dev. stress for Lubby2
		if (line_string.find("$THRESHOLD_DEV_STR") != string::npos)
		{
			in_sd.str(GetLineFromFile1(msp_file));
			in_sd >> threshold_dev_str;
			in_sd.clear();
		}
		//....................................................................
		if (line_string.find("$BIOT_CONSTANT") != string::npos)
		{
			in_sd.str(GetLineFromFile1(msp_file));
			in_sd >> biot_const;
			in_sd.clear();
		}
		//....................................................................
		if (line_string.find("$SOLID_BULK_MODULUS") != string::npos) // WX: 04.2013
		{
			in_sd.str(GetLineFromFile1(msp_file));
			in_sd >> Ks;
			in_sd.clear();
		}
		//....................................................................
		if (line_string.find("BISHOP_COEFFICIENT") != string::npos) // WX
		{
			in_sd.str(GetLineFromFile1(msp_file));
			in_sd >> bishop_model;
			switch (bishop_model)
			{
				case 1: // constant
					in_sd >> bishop_model_value;
					break;
				case 2: // pow(Se, parameter)
					in_sd >> bishop_model_value;
					break;
				case 3: // JM model 3:    if p<bishop_model_value -> bishop_parameter=0.0;  else -> bishop_parameter=1.0
					in_sd >> bishop_model_value;
					break;
				default:
					break;
			}
			in_sd.clear();
		}
		//....................................................................
		if (line_string.find("$STRESS_INTEGRATION_TOLERANCE") != string::npos)
		{
			in_sd.str(GetLineFromFile1(msp_file));
			in_sd >> f_tol >> s_tol;
			in_sd.clear();
		}
		if (line_string.find("$STRESS_UNIT") != string::npos)
		{
			in_sd.str(GetLineFromFile1(msp_file));
			in_sd >> sub_line;
			if (sub_line.compare("MegaPascal") == 0)
				grav_const = 9.81e-6;
			else if (sub_line.find("KiloPascal") == 0)
				grav_const = 9.81e-3;
			in_sd.clear();
		}
		// WX:08.2011 define gravity constant
		if (line_string.find("$GRAVITY_CONSTANT") != string::npos)
		{
			in_sd.str(GetLineFromFile1(msp_file));
			in_sd >> grav_const;
			in_sd.clear();
		}
		if(line_string.find("$GRAVITY_RAMP")!=string::npos)
		{
			in_sd.str(GetLineFromFile1(msp_file));
			in_sd>>grav_curve_id;
			in_sd>>grav_const;
			in_sd.clear();
			gravity_ramp = 1;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$PLASTICITY") != string::npos)
		{
			in_sd.str(GetLineFromFile1(msp_file));
			in_sd >> sub_line;
			in_sd.clear();
			if (sub_line.find("DRUCKER-PRAGER") != string::npos)
			{
				devS = new double[6];
				Plasticity_type = 1;
				// No return mapping
				if (sub_line.find("NORETURNMAPPING") != string::npos)
				{
					Plasticity_type = 10;
					if (sub_line.find("TENSIONCUTOFF") != string::npos) // WX: 08.2010
					{
						Plasticity_type = 11;
						dFtds = new double[6]; // WX: 08.2010 Tensile yield function
						dGtds = new double[6];
						ConstitutiveMatrix = new Matrix(6, 6); // WX: 08.2010
					}
					dFds = new double[6];
					dGds = new double[6];
					D_dFds = new double[6];
					D_dGds = new double[6];
				}
				Size = 5;
				if (Plasticity_type == 11)
					Size = 6;
				/*
				   Material parameters for Cam-Clay model
				   i : parameter
				   0 : The initial yield stress
				   1 : Plastic hardening parameter
				   2 : Internal frection angle
				   3 : Dilatancy angle
				   4 : Localized softening modulus
				   5 : Tensile strength //WX
				 */
			}
			else if (sub_line.find("SINGLE_YIELD_SURFACE") != string::npos)
			{
				Plasticity_type = 2;
				Size = 23;
				AllocateMemoryforSYS();
				/*
				   Material parameters for Single yield surface model
				   i: parameter
				   0: alpha0
				   1: beta0
				   2: delta0
				   3: epsilon0
				   4: kappa0
				   5: gamma0
				   6: m0

				   7: alpha1
				   8: beta1
				   9: delta1
				   10: epsilon1
				   11: kappa1
				   12: gamma1
				   13: m1

				   14: Psi1
				   15: Psi2

				   16: Ch
				   17: Cd
				   18: br
				   19: mr

				   20: Initial stress_xx
				   21: Initial stress_yy
				   22: Initial stress_zz
				 */
			}
			else if (sub_line.find("CAM-CLAY") != string::npos)
			{
				Plasticity_type = 3;
				Size = 10;
				/*
				   Material parameters for Cam-Clay model
				   i: parameter
				   0 : M: slope of the critical line
				   1 : Lamda, the virgin compression index
				   2 : Kappa, swelling index
				   3 : p0, preconsolidation pressure
				   4 : e0, initial void ratio
				   5 : OCR
				   6 : Initial stress_xx
				   7 : Initial stress_yy
				   8 : Initial stress_zz
				   9 : Mimimum p: ( stress_xx + stress_yy + stress_zz )/3
				 */
			}
			else if (sub_line.find("MOHR-COULOMB") != string::npos) // WX
			{
				devS = new double[6];
				TransMatrixA = new Matrix(6, 6); // WX:08.2011
				TransMatrixA_T = new Matrix(6, 6);
				Inv_TransA = new Matrix(6, 6);
				Inv_TransA_T = new Matrix(6, 6);
				TmpDe = new Matrix(6, 6);
				Inv_De = new Matrix(6, 6);
				TMatrix = new Matrix(6, 6);
				TmpMatrix = new Matrix(6, 6);
				TmpMatrix2 = new Matrix(6, 6);
				Dep_l = new Matrix(3, 3);
				dDep_l = new Matrix(3, 3);
				dGds_dFds = new Matrix(6, 6); // WX:08.2011
				ConstitutiveMatrix = new Matrix(6, 6);
				Plasticity_type = 4;
				Size = 6;
				/*
				   i	parameter
				   0	cohesion
				   1	friction angle
				   2	dilatance angle
				   3	tension strength
				   4   hardening curve for friction angle
				   5   hardening curve for cohesion
				 */
			}
			else if (sub_line.find("HOEK-BROWN") != string::npos) // WX
			{
				devS = new double[6];
				ConstitutiveMatrix = new Matrix(6, 6);
				Plasticity_type = 5;
				Size = 4;
				/*
				   i   parameter
				   0   a
				   1   s
				   2   mb
				   3   sigci
				 */
			}
			if (Plasticity_type == 1)
				data_Plasticity = new Matrix(Size + 1);
			else
				data_Plasticity = new Matrix(Size);
			for (i = 0; i < Size; i++)
			{
				in_sd.str(GetLineFromFile1(msp_file));
				in_sd >> (*data_Plasticity)(i);
				in_sd.clear();
			}
		}
		// Solid reaction system
		if (line_string.find("$REACTIVE_SYSTEM") != string::npos)
		{
			in_sd.str(GetLineFromFile1(msp_file));
			in_sd >> reaction_system;
			this->setSolidReactiveSystem(FiniteElement::convertSolidReactiveSystem(reaction_system));
			in_sd.clear();
			if (reaction_system.compare("SINUSOIDAL") == 0)
			{ // For Benchmarks
				in_sd.str(GetLineFromFile1(msp_file));
				in_sd >> reaction_enthalpy; // in J/kg, negative for exothermic composition reaction
				in_sd.clear();
			}
			SetSolidReactiveSystemProperties();
		}

		// Solid non reactive fraction
		if (line_string.find("$NON_REACTIVE_FRACTION") != string::npos)
		{
			in_sd.str(GetLineFromFile1(msp_file));
			in_sd >> non_reactive_solid_volume_fraction >> non_reactive_solid_density;
			in_sd.clear();
		}

		if (line_string.find("$SPECIFIC_HEAT_SOURCE") != string::npos)
		{
			in_sd.str(GetLineFromFile1(msp_file));
			in_sd >> specific_heat_source;
			in_sd.clear();
		}

		if (line_string.find("$MICRO_STRUCTURE_PLAS") != string::npos) // WX:09.2011 for anisotropic plasticity
		{
			if (line_string.find("NORETURNMAPPING") != string::npos)
				Plasticity_type = 44; // direct stress intergration anisotropic plas.
			Plasticity_Bedding = true;
			MicroStruTensor = new double[3];
			Bedding_Norm = new double[3];
			TransMicroStru = new Matrix(6, 6);
			TransMicroStru_T = new Matrix(6, 6);
			TransMicroStru_TInv = new Matrix(6, 6);
			for (i = 0; i < 4; i++)
			{
				in_sd.str(GetLineFromFile1(msp_file));
				in_sd >> line_string;
				if (line_string.find("MICRO_STRUCTURE_TENSOR") != string::npos)
				{
					in_sd >> MicroStruTensor[0] >> MicroStruTensor[1] >> MicroStruTensor[2];
				}
				else if (line_string.find("BEDDING_NORM") != string::npos)
				{
					in_sd >> Bedding_Norm[0] >> Bedding_Norm[1] >> Bedding_Norm[2];
				}
				/*else if(line_string.find("FRICTION_CURVE")!=string::npos)
				{
				in_sd>>bedding_fric_curve;
				}*/
				else if (line_string.find("UNIAXI_COMP_CURVE") != string::npos)
				{
					in_sd >> bedding_uc_curve_order;
					comp_para = new double[bedding_uc_curve_order + 1];
					for (int ii = 0; ii < bedding_uc_curve_order + 1; ii++)
						in_sd >> comp_para[ii];
				}
				else if (line_string.find("TENSION_CURVE") != string::npos)
				{
					in_sd >> bedding_tens_curve_order;
					tens_para = new double[bedding_tens_curve_order + 1];
					for (int ii = 0; ii < bedding_tens_curve_order + 1; ii++)
						in_sd >> tens_para[ii];
				}
				in_sd.clear();
			}
			CalTransMatrixMicroStru(TransMicroStru, Bedding_Norm);
			TransMicroStru->GetTranspose(*TransMicroStru_T);
			Cal_Inv_Matrix(6, TransMicroStru_T, TransMicroStru_TInv);
			// TransMicroStru_T->Write();
			// TransMicroStru_TInv->Write();
			// TransMicroStru->Write();
		}
		in_sd.clear();
	}
	return position;
}

//==========================================================================

/**************************************************************************
   FEMLib-Method:
   Task: Set values for solid reactive system
**************************************************************************/
void CSolidProperties::SetSolidReactiveSystemProperties() // Definition auch in conversion_rate::conversion_rate
{
	if (reaction_system.compare("CaOH2") == 0)
	{
		lower_solid_density_limit = 1665.1; //density Calciumoxid
		upper_solid_density_limit = 2200.0; //density Calciumhydroxid
		reaction_enthalpy = -1.083e+05/PhysicalConstant::MolarMass::Water; //in J/kg (Molar heat of reaction divided by molar mass of water)
		//negative for exothermic composition reaction
		reaction_entropy = -143.5/PhysicalConstant::MolarMass::Water; //in J/kgK
	}
	else if (reaction_system.compare("Mn3O4") == 0)
	{
		lower_solid_density_limit = 4500.0; // density Mn2O3 (Mangan(III)-oxid)
		upper_solid_density_limit = 4860.0; // density Mn3O4 (Mangan(II,III)-oxid)
		reaction_enthalpy = -1.376e+05 / 0.032; // in J/kg (Molar heat of reaction divided by molar mass of oxygen)
		// negative for exothermic composition reaction
		reaction_entropy = -114.1 / 0.032; // in J/kgK
	}
	else if (reaction_system.compare("Z13XBF") == 0)
	{
		lower_solid_density_limit = 1150.0; // Dehydrated Zeolite pellet
		upper_solid_density_limit = 1.e6; // state dependent; not needed for calculation
		reaction_enthalpy = 0.0; // state dependent; calculated where needed
		reaction_entropy = 0.0; // state dependent; calculated where needed
	}
	return;
}

//==========================================================================

/**************************************************************************
FEMLib-Method:
Task: Constructor and destructor
Programing:
08/2004 WW Implementation
04/2005 WW Curve dependency
**************************************************************************/
CSolidProperties::CSolidProperties()
    : data_Youngs(NULL), data_Density(NULL), data_Capacity(NULL), data_Conductivity(NULL), data_Plasticity(NULL),
      data_Creep(NULL)
{
	PoissonRatio = 0.2;
	ThermalExpansion = 0.0;
	biot_const = 1.0;
	// Default, all data are constant
	Density_mode = -1;
	Youngs_mode = -1;
	Capacity_mode = -1;
	Conductivity_mode = -1;
	T_0 = 0.0;
	Creep_mode = -1;
	grav_const = 9.81; // WW
	gravity_ramp = 0;
	excavation = -1; // 12.2009. WW
	excavated = false; // To be .....  12.2009. WW
	E_Function_Model = -1; // WX:06.2012 E dependence
	threshold_dev_str = -1.; // WX:12.2012
	Time_Dependent_E_nv_mode = -1; // WX:01.2013
	Time_Dependent_E_nv_value[0] = -1; // WX:01.2013

	SwellingPressureType = -1;
	Max_SwellingPressure = 0.0;
	// Default, elasticity
	Plasticity_type = -1;

	E = Lambda = G = K = 0.0;
	Ks = 0.; // WX: 04.2013
	devS = NULL;
	axisymmetry = false;
	dl2 = 0.0;
	// SYS
	d2G_dSdS = NULL;
	d2G_dSdM = NULL;
	LocalJacobi = NULL;
	inv_Jac = NULL;
	sumA_Matrix = NULL;
	rhs_l = NULL;
	x_l = NULL;
	Li = NULL;
	// Drucker-Prager
	dFds = NULL;
	dGds = NULL;
	D_dFds = NULL;
	D_dGds = NULL;
	dFtds = NULL; // WX
	dGtds = NULL; // WX
	ConstitutiveMatrix = NULL; // WX
	// WX:08.2011 for mohr coulomb
	TransMatrixA = NULL;
	TransMatrixA_T = NULL;
	Inv_TransA = NULL;
	Inv_TransA_T = NULL;
	TmpDe = NULL;
	Inv_De = NULL;
	TMatrix = NULL;
	TmpMatrix = NULL;
	TmpMatrix2 = NULL;
	Dep_l = NULL;
	dDep_l = NULL;
	dGds_dFds = NULL;
	Plasticity_Bedding = false; // WX:09.2011
	TransMicroStru = NULL; // WX:11.2011
	MicroStruTensor = NULL; // WX:11.2011
	TransMicroStru_T = NULL;
	TransMicroStru_TInv = NULL;
	comp_para = NULL;
	tens_para = NULL;
	// Curve variable type
	// 0: Time
	// 1: ...
	CurveVariableType_Conductivity = -1;
	mode = 0; // Gauss point values //OK
	//
	s_tol = 1e-9;
	f_tol = 1e-9;

	Crotm = NULL; // rotation matrix for matrices: UJG 25.11.2009
	D_tran = NULL; // UJG/WW

	// Thermal conductivity tensor (default: iso)
	thermal_conductivity_tensor_type = 0;
	thermal_conductivity_tensor_dim = 1;
	thermal_conductivity_tensor[0] = 1.0;

	// Reactive system
	reaction_system = "INERT";
	_reactive_system = FiniteElement::INERT;
	lower_solid_density_limit = 0.0;
	upper_solid_density_limit = 0.0;
	reaction_enthalpy = 0.0;
	reaction_entropy = 0.0;
	non_reactive_solid_volume_fraction = 0.0;
	non_reactive_solid_density = 0.0;
	material_minkley = NULL;
	material_burgers = NULL;

	specific_heat_source = 0.0;

	bishop_model = -1; // 15.08.2011. WW
}
CSolidProperties::~CSolidProperties()
{
	if (data_Density)
		delete data_Density;
	if (data_Density)
		delete data_Youngs;
	if (data_Plasticity)
		delete data_Plasticity;
	if (data_Capacity)
		delete data_Capacity;
	if (data_Conductivity)
		delete data_Conductivity;
	if (data_Creep)
		delete data_Creep;
	if (devS)
		delete[] devS;

	if (Crotm)
		delete Crotm; // rotation matrix for matrices: UJG 25.11.2009
	if (D_tran)
		delete D_tran; // rotation matrix for matrices: UJG 25.11.2009
	data_Density = NULL;
	data_Youngs = NULL;
	data_Plasticity = NULL;
	data_Capacity = NULL;
	data_Conductivity = NULL;
	data_Creep = NULL;
	devS = NULL;
	Crotm = NULL;
	D_tran = NULL;

	if (d2G_dSdS)
		delete d2G_dSdS;
	if (d2G_dSdM)
		delete d2G_dSdM;
	if (LocalJacobi)
		delete LocalJacobi; // To store local Jacobi matrix
	if (inv_Jac)
		delete inv_Jac; // To store the inverse of the  Jacobi matrix
	if (sumA_Matrix)
		delete sumA_Matrix;
	if (rhs_l)
		delete[] rhs_l; // To store local unknowns of 15
	if (x_l)
		delete[] x_l; // To store local unknowns of 15
	if (Li)
		delete[] Li;

	if (dFds)
		delete[] dFds;
	if (dGds)
		delete[] dGds;
	if (D_dFds)
		delete[] D_dFds;
	if (D_dGds)
		delete[] D_dGds;
	if (dFtds)
		delete[] dFtds; // WX:
	if (dGtds)
		delete[] dGtds; // WX:
	if (ConstitutiveMatrix)
		delete ConstitutiveMatrix; // WX:
	if (TransMatrixA)
		delete TransMatrixA; // WX:08.2011
	if (TransMatrixA_T)
		delete TransMatrixA_T;
	if (Inv_TransA)
		delete Inv_TransA;
	if (Inv_TransA_T)
		delete Inv_TransA_T;
	if (TmpDe)
		delete TmpDe;
	if (Inv_De)
		delete Inv_De;
	if (TMatrix)
		delete TMatrix;
	if (TmpMatrix)
		delete TmpMatrix;
	if (TmpMatrix2)
		delete TmpMatrix2;
	if (Dep_l)
		delete Dep_l;
	if (dDep_l)
		delete dDep_l;
	if (dGds_dFds)
		delete dGds_dFds; // WX:08.2011
	if (TransMicroStru)
		delete TransMicroStru; // WX:11.2011
	if (MicroStruTensor)
		delete MicroStruTensor; // WX:11.2011
	if (TransMicroStru_T)
		delete TransMicroStru_T;
	if (TransMicroStru_TInv)
		delete TransMicroStru_TInv;
	if (comp_para)
		delete comp_para;
	if (tens_para)
		delete tens_para;

	dFds = NULL;
	dGds = NULL;
	D_dFds = NULL;
	D_dGds = NULL;
	dFtds = NULL; // WX:
	dGtds = NULL; // WX:
	TransMatrixA = NULL; // WX:08.2011
	TransMatrixA_T = NULL;
	Inv_TransA = NULL;
	Inv_TransA_T = NULL;
	TmpDe = NULL;
	Inv_De = NULL;
	TMatrix = NULL;
	TmpMatrix = NULL;
	TmpMatrix2 = NULL;
	Dep_l = NULL;
	dDep_l = NULL;
	dGds_dFds = NULL; // WX:08.2011
	ConstitutiveMatrix = NULL; // WX:
	TransMicroStru = NULL; // WX:11.2011
	MicroStruTensor = NULL; // WX:11.2011
	TransMicroStru_T = NULL;
	TransMicroStru_TInv = NULL;
	comp_para = NULL;
	tens_para = NULL; // WX
	d2G_dSdS = NULL;
	d2G_dSdM = NULL;
	LocalJacobi = NULL;
	inv_Jac = NULL;
	sumA_Matrix = NULL;
	rhs_l = NULL;
	x_l = NULL;
	Li = NULL;
	material_minkley = NULL;
	material_burgers = NULL;
	smath = NULL;
}
//----------------------------------------------------------------------------

/**************************************************************************
   FEMLib-Method: CSolidProperties::CalulateValue
   Task: Linear interpolation
   Programing:
   08/2004 WW Implementation
**************************************************************************/
double CSolidProperties::CalulateValue(const Matrix* data, const double x) const
{
	int i;
	double Val = 0.0;

	if (!data)
		return 0.0;

	const int Size = data->Rows();
	// If given varial is not in the range of data
	if (x < (*data)(0, 0))
		Val = (*data)(0, 1);
	else if (x >= (*data)(Size - 1, 0))
		Val = (*data)(Size - 1, 1);
	else
		for (i = 0; i < Size - 1; i++)
			if ((x >= (*data)(i, 0)) && (x < (*data)(i + 1, 0)))
				Val = (x - (*data)(i, 0)) * ((*data)(i + 1, 1) - (*data)(i, 1)) / ((*data)(i + 1, 0) - (*data)(i, 0))
				      + (*data)(i, 1);

	return Val;
}

/**************************************************************************
   FEMLib-Method: CSolidProperties::Density(double refence = 0.0)
   Task: Get density
   Programing:
   08/2004 WW Implementation
**************************************************************************/
double CSolidProperties::Density(double refence)
{
	double val = 0.0;
	switch (Density_mode)
	{
		case 0:
			val = CalulateValue(data_Density, refence);
			break;
		case 1:
			val = (*data_Density)(0);
			break;
	}
	return val;
}

/**************************************************************************
   FEMLib-Method: CSolidProperties::Heat_Capacity(const double refence = 0.0) const
   Task: Get heat capacity
   Programing:
   08/2004 WW Implementation
   **************************************************************************/
double CSolidProperties::Heat_Capacity(double refence)
{
	double val = 0.0;
	switch (Capacity_mode)
	{
		case 0:
			val = CalulateValue(data_Capacity, refence);
			break;
		case 1:
			val = (*data_Capacity)(0);
			break;
		case 3:
			// WW        val=1.38*(273.15+refence)+732.5;
			val = 1.38 * refence + 732.5;
			break;
		case 4: // solid capacity depending on solid density (for thermochemical heat storage) - TN
			// refence contains value of solid density (current)
			val = (*data_Capacity)(0)
			      + ((*data_Capacity)(1) - (*data_Capacity)(0)) / ((*data_Capacity)(3) - (*data_Capacity)(2))
			            * (refence - (*data_Capacity)(2));
			break;
		case 5:
			// TODO [CL]: TES heat capacity model number changed 4-->5
			val = (*data_Capacity)(0)
			      + ((*data_Capacity)(1) - (*data_Capacity)(0))
			            / (upper_solid_density_limit - lower_solid_density_limit)
			            * (refence - lower_solid_density_limit);
			break;
		case 6: // in a sorption system
		{
			// TODO [CL]: TES heat capacity model number changed 5-->6
			double C = refence / lower_solid_density_limit - 1.0;
			val = lower_solid_density_limit / refence * ((*data_Capacity)(0) + C * (*data_Capacity)(1));
		}
		break;
		case 7:
			//linear density average (compare model 5)
			val = lower_solid_density_limit*(*data_Capacity)(0);
			val += (upper_solid_density_limit*(*data_Capacity)(1)-lower_solid_density_limit*(*data_Capacity)(0))/
				(upper_solid_density_limit  - lower_solid_density_limit)*(refence - lower_solid_density_limit);
			val /= refence;
		break;
		default:
			val = (*data_Capacity)(0);
			break;
	}
	return val;
}

/*************************************************************************
FEMLib-Method:
Task: Get heat phase change temperature
Programing:
09/2005 WW Implementation
**************************************************************************/
bool CSolidProperties::CheckTemperature_in_PhaseChange(const double T0, const double T1)
{
	bool stat = false;
	double T_a = (*data_Capacity)(2);
	double T_b = (*data_Capacity)(2) + (*data_Capacity)(3);
	switch (Capacity_mode)
	{
		case 0:
			break;
		case 1:
			break;
		case 2:
			if ((T1 > T0) && (*data_Capacity)(4) > 0.0)
			{
				if ((T0 < T_a) && (T1 > T_a))
					stat = true;
				else if ((T0 < T_a) && (T1 > T_b))
					stat = true;
				else if ((T0 < T_b) && (T1 > T_b))
					stat = true;
				else if ((T1 >= T_a) && (T1 <= T_b))
					stat = true;
			}
			break;
	}
	return stat;
}

/*************************************************************************
FEMLib-Method:
Task: Get heat capacity with boiling model
      latent_factor=density_w*porosity*saturation
Programing:
09/2005 WW Implementation
**************************************************************************/
double CSolidProperties::Heat_Capacity(double temperature, double porosity, double Sat)
{
	double val = 0.0;
	double sign = 1;
	double dens = fabs(Density());
	// If sign =1, temperature increases in evolution. Otherwise, decrease.
	if (fabs(temperature) > 1.e-9)
		sign = fabs(temperature) / temperature;
	temperature *= sign;
	CFluidProperties* m_mfp = NULL;
	m_mfp = mfp_vector[0];

	// 0. Wet capacity
	// 1. Dry capacity
	// 2. Boiling temperature
	// 3. Boiling temperature range
	// 4. Latent of vaporization
	if (temperature < (*data_Capacity)(2)) // Wet
		val = dens * (*data_Capacity)(0);
	else if ((temperature >= (*data_Capacity)(2)) && (temperature < ((*data_Capacity)(2) + (*data_Capacity)(3))))
	{
		if (sign > 0.0) // Temperature increase
			val = dens * (*data_Capacity)(0)
			      + dens * ((*data_Capacity)(1) - (*data_Capacity)(0)) * (temperature - (*data_Capacity)(2))
			            / (*data_Capacity)(3)
			      + porosity * Sat * m_mfp->Density() * (*data_Capacity)(4) / (*data_Capacity)(3);
		else
			val = dens * (*data_Capacity)(0)
			      + dens * ((*data_Capacity)(1) - (*data_Capacity)(0)) * (temperature - (*data_Capacity)(2))
			            / (*data_Capacity)(3);
	}
	else
		val = (*data_Capacity)(1);
	return val;
}

/**************************************************************************
   FEMLib-Method:
   Task: Get heat capacity with boiling model
      latent_factor=density_w*porosity*saturation
     temperature increases in evolution
   Programing:
   09/2005 WW Implementation
**************************************************************************/
double CSolidProperties::Enthalpy(double temperature, const double latent_factor)
{
	double val = 0.0;
	double T0 = 0.0; //�C, a reference temperature
	double dT = 0.0;

	// 0. Wet capacity
	// 1. Dry capacity
	// 2. Boiling temperature
	// 3. Boiling temperature range
	// 4. Latent of vaporization
	if (temperature < (*data_Capacity)(2)) // Wet
		val = (*data_Capacity)(0) * (temperature - T0);
	else if ((temperature >= (*data_Capacity)(2)) && (temperature < ((*data_Capacity)(2) + (*data_Capacity)(3))))
	{
		dT = temperature - (*data_Capacity)(2) - T0;
		val = (*data_Capacity)(0) * ((*data_Capacity)(2) - T0);
		val += dT * latent_factor * (*data_Capacity)(4) / (*data_Capacity)(3) / Density();
		val += (*data_Capacity)(0) * dT
		       + 0.5 * dT * dT * ((*data_Capacity)(1) - (*data_Capacity)(0)) / (*data_Capacity)(3);
	}
	else
	{
		dT = (*data_Capacity)(3);
		val = (*data_Capacity)(0) * ((*data_Capacity)(2) - T0);
		val += dT * latent_factor * (*data_Capacity)(4) / (*data_Capacity)(3) / Density();
		val += (*data_Capacity)(0) * dT
		       + 0.5 * dT * dT * ((*data_Capacity)(1) - (*data_Capacity)(0)) / (*data_Capacity)(3);
		val += (*data_Capacity)(1) * (temperature - (*data_Capacity)(3) - (*data_Capacity)(2) - T0);
	}
	return val * fabs(Density());
}

/**************************************************************************
   FEMLib-Method: CSolidProperties::Heat_Capacity(const double reference = 0.0) const
   Task: Get density
   Programing:
   08/2004 WW Implementation
**************************************************************************/
double CSolidProperties::Heat_Conductivity(double reference)
{
	double val = 0.0;
	int gueltig;
	switch (Conductivity_mode)
	{
		case 0:
			val = CalulateValue(data_Conductivity, reference);
			break;
		case 1:
			val = (*data_Conductivity)(0);
			break;
		case 2:
		{
			const double* k_T = data_Conductivity->getEntryArray();

			// 0. Wet conductivity
			// 1. Dry conductivity
			// 2. Boiling temperature
			// 3. Boiling temperature range
			if (reference < k_T[2]) // Wet
				val = k_T[0];
			else if ((reference >= k_T[2]) && (reference < (k_T[2] + k_T[3])))
				val = k_T[0] + (k_T[0] - k_T[1]) * (reference - k_T[2]) / k_T[3];
			else
				val = k_T[1];
		}
		break;
		case 3: // reference: saturation
			// val = 1.28-0.71/(1+10.0*exp(reference-0.65));  //MX
			val = 1.28 - 0.71 / (1 + exp(10.0 * (reference - 0.65)));
			break;
		case 30: // Another model for bentonite. 10.2013. WW
		{
			// val = k_max-k_min/(1+10.0*exp(reference-S0));
			const double* k_T = data_Conductivity->getEntryArray();
			//		val = k_T[0] - (k_T[0]-k_T[1]) / (1 + exp(10.0 * (reference - k_T[2])));
			val = k_T[0] + k_T[1] * (reference - k_T[2]);
		}
		break;
		case 4: // 21.12.2009. WW
			val = CalulateValue(data_Conductivity, reference);
			break;
		case 5: // DECOVALEX2015, TaskB2 JM
			CalPrimaryVariable(capacity_pcs_name_vector);
			val = GetMatrixValue(primary_variable[0] + T_0, primary_variable[1], name, &gueltig);
			break;
	}
	return val;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   09/2004 OK MSP implementation
   03/2005 WW Conductivity from input file
   04/2011 NW Conductivity tensor
   last modification:
   ToDo: geo_dimension
**************************************************************************/
void CSolidProperties::HeatConductivityTensor(const int dim, double* tensor, int group)
{
	group = group; // OK411
	// static double tensor[9];
	double temperature = 0.0;
	double saturation = 0.0;
	int i = 0;
	// WW
	CalPrimaryVariable(conductivity_pcs_name_vector);
	//--------------------------------------------------------------------
	// MMP medium properties
	// WW CMediumProperties *m_mmp = NULL;
	// WW m_mmp = mmp_vector[group];  //MX
	// Test for DECOVALEX
	// thermal_conductivity_tensor[0] =0.5+0.8*PCSGetELEValue(number,NULL,theta,"SATURATION1");

	//--------------------------------------------------------------------
	// There are a number of cases where the heat conductivity tensor is defined by the capacity model;
	double base_thermal_conductivity = .0;
	int gueltig = 1;
	switch (Conductivity_mode)
	{
		case 0:
			// WW
			base_thermal_conductivity = Heat_Conductivity(primary_variable[0]);
			break;
		case 1:
			// WW
			base_thermal_conductivity = Heat_Conductivity(0);
			break;
		case 2: // Boiling model. DECOVALEX THM2
			temperature = primary_variable[0];
			base_thermal_conductivity = Heat_Conductivity(temperature);
			break;
		case 3: // DECOVALEX THM1
			saturation = primary_variable[1];
			base_thermal_conductivity = Heat_Conductivity(saturation);
			break;
		case 5:
			base_thermal_conductivity = GetMatrixValue(primary_variable[1], primary_variable[0] + T_0, name, &gueltig);
			break;
		default: // Normal case
			cout << "***Error in CSolidProperties::HeatConductivityTensor(): conductivity mode is not supported "
			     << "\n";
			// base_thermal_conductivity = Heat_Conductivity();
	}

	//--------------------------------------------------------------------
	// Set unit tensor
	// check
	if (thermal_conductivity_tensor_type > 0 && dim != thermal_conductivity_tensor_dim)
		cout << "***Error in CSolidProperties::HeatConductivityTensor(): problem dimension and the given tensor "
		        "dimension are not same."
		     << "\n";
	// reset
	for (i = 0; i < 9; i++)
		tensor[i] = 0.0;
	// set
	switch (thermal_conductivity_tensor_type)
	{
		case 0: // Iso
			for (i = 0; i < dim; i++)
				tensor[i * dim + i] = 1.0;
			break;
		case 1: // Ortho
			switch (thermal_conductivity_tensor_dim)
			{
				case 1: // 1-D
					tensor[0] = thermal_conductivity_tensor[0];
					break;
				case 2: // 2-D
					tensor[0] = thermal_conductivity_tensor[0];
					tensor[3] = thermal_conductivity_tensor[1];
					break;
				case 3: // 3-D
					tensor[0] = thermal_conductivity_tensor[0];
					tensor[4] = thermal_conductivity_tensor[1];
					tensor[8] = thermal_conductivity_tensor[2];
					break;
			}
			break;
		case 2: // Aniso
			for (i = 0; i < thermal_conductivity_tensor_dim * thermal_conductivity_tensor_dim; i++)
				tensor[i] = thermal_conductivity_tensor[i];
			break;
	}

	//--------------------------------------------------------------------
	// tensor form
	for (i = 0; i < dim * dim; i++)
		tensor[i] *= base_thermal_conductivity;
}

/**************************************************************************
   FEMLib-Method: CSolidProperties::Youngs_Modulus(const double reference = 0.0) const
   Task: Get density
   Programing:
   08/2004 WW Implementation
**************************************************************************/
double CSolidProperties::Youngs_Modulus(double reference)
{
	double val = 0.0;
	switch (Youngs_mode)
	{
		case 0:
			val = CalulateValue(data_Youngs, reference);
			break;
		case 1:
			val = (*data_Youngs)(0);
			break;
	}
	return val;
}

//-------------------------------------------------------------------------
// Kronecker function
double CSolidProperties::Kronecker(const int ii, const int jj)
{
	double delta = 0.0;
	if (ii == jj)
		delta = 1.0;
	return delta;
}

/**************************************************************************
   FEMLib-Method: CSolidProperties::Voigt_to_Kelvin_Strain
   Task: Maps a strain vector in Voigt notation into on in Kelvin notation
   This is an auxilliary routine that will not be needed when the entire FE
   code is set up in Kelvin notation
   Programing:
   06/2015 TN Implementation
**************************************************************************/
KVec CSolidProperties::Voigt_to_Kelvin_Strain(const double *voigt_strain)
{
	KVec kelvin_strain;
	for (size_t i=0; i<3; i++)
	{
		//Normal components
		kelvin_strain(i) = voigt_strain[i];
		//Shear components
		kelvin_strain(i+3) = voigt_strain[i+3]/sqrt(2.);
	}
	return kelvin_strain;
}

/**************************************************************************
   FEMLib-Method: CSolidProperties::Voigt_to_Kelvin_Stress
   Task: Maps a stress vector in Voigt notation into on in Kelvin notation
   This is an auxilliary routine that will not be needed when the entire FE
   code is set up in Kelvin notation
   Programing:
   06/2015 TN Implementation
**************************************************************************/
KVec CSolidProperties::Voigt_to_Kelvin_Stress(const double *voigt_stress)
{
	KVec kelvin_stress;
	for (size_t i=0; i<3; i++)
	{
		//Normal components
		kelvin_stress(i) = voigt_stress[i];
		//Shear components
		kelvin_stress(i+3) = voigt_stress[i+3]*sqrt(2.);
	}
	return kelvin_stress;
}

/**************************************************************************
   FEMLib-Method: CSolidProperties::Kelvin_to_Voigt_Strain()
   Task: Maps a strain vector in Kelvin notation into on in Voigt notation
   This is an auxilliary routine that will not be needed when the entire FE
   code is set up in Kelvin notation
   Programing:
   06/2015 TN Implementation
**************************************************************************/
void CSolidProperties::Kelvin_to_Voigt_Strain(const KVec &kelvin_strain, double *voigt_strain)
{
	for (size_t i=0; i<3; i++)
	{
		//Normal components
		voigt_strain[i] = kelvin_strain(i);
		//Shear components
		voigt_strain[i+3] = kelvin_strain(i+3)*sqrt(2.);
	}
	return;
}

/**************************************************************************
   FEMLib-Method: CSolidProperties::Kelvin_to_Voigt_Stress()
   Task: Maps a stress vector in Kelvin notation into on in Voigt notation
   This is an auxilliary routine that will not be needed when the entire FE
   code is set up in Kelvin notation
   Programing:
   06/2015 TN Implementation
**************************************************************************/
void CSolidProperties::Kelvin_to_Voigt_Stress(const KVec &kelvin_stress, double *voigt_stress)
{
	for (size_t i=0; i<3; i++)
	{
	    //Normal components
	    voigt_stress[i] = kelvin_stress(i);
	    //Shear components
	    voigt_stress[i+3] = kelvin_stress(i+3)/sqrt(2.);
	}
	return;
}


/**************************************************************************
   FEMLib-Method: CSolidProperties::ExtractConsistentTangent()
   Task: general routine to get 6x6 consistent tangent from local Newton iteration of material functionals
   Programing:
   06/2014 TN Implementation
**************************************************************************/
void CSolidProperties::ExtractConsistentTangent(const Eigen::MatrixXd &Jac, const Eigen::MatrixXd &dGdE, const bool pivoting, Eigen::Matrix<double,6,6> &dsigdE)
{
	const unsigned int local_dim(Jac.cols());
	//Check Dimensions
	assert(local_dim == dGdE.rows() && dGdE.cols() == 6);

	Eigen::MatrixXd dzdE(local_dim,6);
	//solve linear system
	if (pivoting)
		dzdE = Jac.fullPivLu().solve(-1.0*dGdE); //Could consider moving to different Eigen solver.
	else
		dzdE = Jac.householderQr().solve(-1.0*dGdE); //Could consider moving to different Eigen solver.
	//in-built Gauss elimination solver was at least 4 OoM more inaccurate.

	//Extract matrix part relevant for global tangent
	dsigdE = dzdE.block<6,6>(0,0);
}


/**************************************************************************
   FEMLib-Method: CSolidProperties::LocalNewtonBurgers()
   Task: general local Newton routine to integrate inelastic material models
   Programing:
   06/2014 TN Implementation
   03/2015 NB Modified
**************************************************************************/
void CSolidProperties::LocalNewtonBurgers(const double dt, double* strain_curr,
		double* stress_curr, double* strain_K_curr, double* strain_M_curr,
        Matrix* Consistent_Tangent, bool Output, double Temperature)
{
	//stress, strain, internal variable
	KVec sig_j, eps_K_j, eps_M_j;
	//deviatoric stress, strain
	KVec sigd_j;
	//local residual vector and Jacobian
	Eigen::Matrix<double,18,1> res_loc, inc_loc;
	Eigen::Matrix<double,18,18> K_loc;
	double sig_eff(1.);

	//initialisation of Kelvin vectors
	//Note: Can be done in one loop instead of 5 if done right here.
	const KVec eps_i(Voigt_to_Kelvin_Strain(strain_curr));
	const KVec eps_K_t(Voigt_to_Kelvin_Strain(strain_K_curr));
	const KVec eps_M_t(Voigt_to_Kelvin_Strain(strain_M_curr));
	eps_M_j = eps_M_t;
	eps_K_j = eps_K_t;

	//calculation of deviatoric and spherical parts
	const double e_i = eps_i(0) + eps_i(1) + eps_i(2);
	const KVec epsd_i(SolidMath::P_dev*eps_i);

	//dimensionless stresses
	sigd_j = 2.0 * (epsd_i - eps_M_t - eps_K_t); //initial guess as elastic predictor

	//Calculate effective stress and update material properties
	sig_eff = SolidMath::CalEffectiveStress(sigd_j);

	if (!T_Process)
	{
		Temperature= material_burgers->T_ref;
		material_burgers->B = 1;
		material_burgers->Q = 0; // for cutting off Arrhenius term
	}

	material_burgers->UpdateBurgersProperties(sig_eff*material_burgers->GM, Temperature);

	//initial evaluation of residual and Jacobian
	material_burgers->CalResidualBurgers(dt,epsd_i,sigd_j,eps_K_j,eps_K_t,eps_M_j,eps_M_t,res_loc);
	//initial evaluation of Jacobian
	material_burgers->CalJacobianBurgers(dt,K_loc,sig_eff,sigd_j,eps_K_j); //Note - With constant properties a single evaluation would be sufficient.

	//Loop variables
	int counter = 0;
	const int counter_max(20);
	const double local_tolerance(1.e-16);

	if (Output) //Need not be performant;
	{
		CRFProcess* m_pcs  = PCSGet("DEFORMATION");
		ofstream Dum("local.txt", ios::app);
		Dum << aktueller_zeitschritt << " " << m_pcs->GetIteSteps() << " " << counter+1 << " " << res_loc.norm() <<  " initial" << std::endl;
		Dum.close();
	};
	//    for (int counter(0); counter<counter_max && res_loc.norm() > local_tolerance; ++counter)
	while (res_loc.norm() > local_tolerance && counter < counter_max)
	{
		counter++;
		//Get Jacobian
		material_burgers->CalJacobianBurgers(dt,K_loc,sig_eff,sigd_j,eps_K_j);//for solution dependent Jacobians
		//Solve linear system
		//inc_loc = K_loc.fullPivHouseholderQr().solve(res_loc); //other linear solvers (faster but less accurate) can be considered.
		//inc_loc = K_loc.fullPivLu().solve(res_loc); //other linear solvers (faster but less accurate) can be considered
		//inc_loc = K_loc.colPivHouseholderQr().solve(res_loc);
		inc_loc = K_loc.householderQr().solve(-res_loc);
		//increment solution vectors
		sigd_j += inc_loc.block<6,1>(0,0);
		eps_K_j += inc_loc.block<6,1>(6,0);
		eps_M_j += inc_loc.block<6,1>(12,0);
		//Calculate effective stress and update material properties
		sig_eff = SolidMath::CalEffectiveStress(sigd_j);
		material_burgers->UpdateBurgersProperties(sig_eff*material_burgers->GM, Temperature);
		//evaluation of new residual
		material_burgers->CalResidualBurgers(dt,epsd_i,sigd_j,eps_K_j,eps_K_t,eps_M_j,eps_M_t,res_loc);
		if (Output)
		{
			CRFProcess* m_pcs  = PCSGet("DEFORMATION");
			ofstream Dum("local.txt", ios::app);
			Dum << aktueller_zeitschritt << " " << m_pcs->GetIteSteps() << " " << counter+1 << " " << res_loc.norm() <<  " visco" << std::endl;
			Dum.close();
		};
	}
	if (counter == counter_max)
		std::cout << "WARNING: Maximum iteration number needed in LocalNewtonBurgers. Convergence not guaranteed." << std::endl;

	//dGdE matrix and dsigdE matrix
	Eigen::Matrix<double,18,6> dGdE;
	Eigen::Matrix<double,6,6> dsigdE;
	
	//Calculate dGdE for time step
	material_burgers->CaldGdEBurgers(dGdE);
	//get dsigdE matrix
	ExtractConsistentTangent(K_loc,dGdE,false,dsigdE);

	//add hydrostatic part to stress and tangent
	sig_j = material_burgers->GM * sigd_j + material_burgers->KM * e_i * SolidMath::ivec;
	dsigdE = material_burgers->GM * dsigdE*SolidMath::P_dev + 3.*material_burgers->KM*SolidMath::P_sph;

	//Sort into Consistent Tangent matrix for global Newton iteration and into standard OGS arrays
	Kelvin_to_Voigt_Stress(sig_j,stress_curr);
	Kelvin_to_Voigt_Strain(eps_K_j,strain_K_curr);
	Kelvin_to_Voigt_Strain(eps_M_j,strain_M_curr);
	for (size_t i=0; i<3; i++)
	{
		for (size_t j=0; j<3; j++){
			(*Consistent_Tangent)(i,j) = dsigdE(i,j);
			(*Consistent_Tangent)(i,j+3) = dsigdE(i,j+3)/sqrt(2.);//from local to global shear components (Kelvin to Voigt)
			(*Consistent_Tangent)(i+3,j) = dsigdE(i+3,j)/sqrt(2.);//from local to global shear components (Kelvin to Voigt)
			(*Consistent_Tangent)(i+3,j+3) = dsigdE(i+3,j+3)/2.;//from local to global shear components (Kelvin to Voigt)
		}
	}
}

/**************************************************************************
   FEMLib-Method: CSolidProperties::LocalNewtonMinkley()
   Task: general local Newton routine to integrate inelastic material models
   Programing:
   06/2015 TN Implementation
**************************************************************************/
void CSolidProperties::LocalNewtonMinkley(const double dt, double* strain_curr, double* stress_curr, double* eps_K_curr,
										  double* eps_M_curr, double* eps_pl_curr, double& e_pl_v, double& e_pl_eff,
										  double& lam, Matrix* Consistent_Tangent,bool Output, double Temperature)
{
	//stress, strain, internal variable
	KVec sig_j, eps_K_j, eps_M_j, eps_pl_j;
	//deviatoric stress, strain
	KVec sigd_j;
	const double e_pl_v_t = e_pl_v, e_pl_eff_t = e_pl_eff;
	//local residual vector and Jacobian
	Eigen::Matrix<double,18,1> res_loc, inc_loc;
	Eigen::Matrix<double,18,18> K_loc;
	double sig_eff;
	Eigen::Matrix<double,6,6> dsigdE;

	//initialisation of Kelvin vectors
	//Note: Can be done in one loop instead of 5 if done right here.
	const KVec eps_i(Voigt_to_Kelvin_Strain(strain_curr));
	const KVec eps_K_t(Voigt_to_Kelvin_Strain(eps_K_curr));
	const KVec eps_M_t(Voigt_to_Kelvin_Strain(eps_M_curr));
	const KVec eps_pl_t(Voigt_to_Kelvin_Strain(eps_pl_curr));
	eps_M_j = eps_M_t;
	eps_K_j = eps_K_t;
	eps_pl_j = eps_pl_t;

	//calculation of deviatoric and spherical parts
	const double e_i = eps_i(0) + eps_i(1) + eps_i(2);
	const KVec epsd_i(SolidMath::P_dev*eps_i);

	//dimensionless stresses
	sigd_j = 2.0 * (epsd_i - eps_M_j - eps_K_j - eps_pl_j); //initial guess as elastic predictor
	sig_j = sigd_j + material_minkley->KM0/material_minkley->GM0 * (e_i - e_pl_v) * SolidMath::ivec;
	//std::cout << "Stress guesstimate_zz " << sig_j(2)*material_minkley->GM0 << std::endl;

	//Calculate effective stress and update material properties
	sig_eff = SolidMath::CalEffectiveStress(sigd_j);

	material_minkley->UpdateMinkleyProperties(sig_eff*material_minkley->GM0, e_pl_eff, Temperature);

	//initial evaluation of residual and Jacobian
	material_minkley->CalViscoelasticResidual(dt,epsd_i,e_i,e_pl_v,sig_j,eps_K_j,eps_K_t,eps_M_j,eps_M_t,eps_pl_j,res_loc);
	//initial evaluation of Jacobian
	material_minkley->CalViscoelasticJacobian(dt,sig_j,sig_eff,K_loc);

	//Loop variables
	int counter = 0;
	const int counter_max(20);
	const double local_tolerance(1.e-16);

	if (Output)
	{
		CRFProcess* m_pcs  = PCSGet("DEFORMATION");
		ofstream Dum("local.txt", ios::app);
		Dum << aktueller_zeitschritt << " " << m_pcs->GetIteSteps() << " " << counter+1 << " " << res_loc.norm() << " initial" << std::endl;
		Dum.close();
	}

	while (res_loc.norm() > local_tolerance && counter < counter_max)
	{
		counter++;
		//Get Jacobian
		material_minkley->CalViscoelasticJacobian(dt,sig_j,sig_eff,K_loc);
		//Solve linear system
		inc_loc = K_loc.fullPivHouseholderQr().solve(-res_loc);
		//increment solution vectors
		sig_j += inc_loc.block<6,1>(0,0);
		eps_K_j += inc_loc.block<6,1>(6,0);
		eps_M_j += inc_loc.block<6,1>(12,0);
		//Calculate effective stress and update material properties
		sig_eff = SolidMath::CalEffectiveStress(SolidMath::P_dev*sig_j);
		material_minkley->UpdateMinkleyProperties(sig_eff*material_minkley->GM0, e_pl_eff, Temperature);
		//evaluation of new residual
		material_minkley->CalViscoelasticResidual(dt,epsd_i,e_i,e_pl_v,sig_j,eps_K_j,eps_K_t,eps_M_j,eps_M_t,eps_pl_j,res_loc);
		if (Output)
		{
			CRFProcess* m_pcs  = PCSGet("DEFORMATION");
			ofstream Dum("local.txt", ios::app);
			Dum << aktueller_zeitschritt << " " << m_pcs->GetIteSteps() << " " << counter+1 << " " << res_loc.norm() <<  " visco" << std::endl;
			Dum.close();
		};
	}
	if (!(material_minkley->YieldMohrCoulomb(sig_j*material_minkley->GM0) < local_tolerance))
	{
		Eigen::Matrix<double,27,1> res_loc_p, inc_loc_p;
		Eigen::Matrix<double,27,27> K_loc_p;

		material_minkley->CalViscoplasticResidual(dt,epsd_i,e_i,sig_j,eps_K_j,eps_K_t,eps_M_j,eps_M_t,eps_pl_j,eps_pl_t, \
												  e_pl_v,e_pl_v_t,e_pl_eff,e_pl_eff_t,lam,res_loc_p);
		material_minkley->CalViscoplasticJacobian(dt,sig_j,sig_eff,lam,K_loc_p);
		while (res_loc_p.norm() > local_tolerance && counter < counter_max)
		{
			counter++;
			//std::cout << "iter " << counter << std::endl;
			//Get Jacobian
			material_minkley->CalViscoplasticJacobian(dt,sig_j,sig_eff,lam,K_loc_p);
			//Solve linear system
			inc_loc_p = K_loc_p.fullPivHouseholderQr().solve(-res_loc_p);
			//increment solution vectors
			sig_j += inc_loc_p.block<6,1>(0,0);
			eps_K_j += inc_loc_p.block<6,1>(6,0);
			eps_M_j += inc_loc.block<6,1>(12,0);
			eps_pl_j += inc_loc_p.block<6,1>(18,0);
			e_pl_v += inc_loc_p.block<1,1>(24,0)(0);
			e_pl_eff += inc_loc_p.block<1,1>(25,0)(0);
			lam += inc_loc_p.block<1,1>(26,0)(0);
			//Calculate effective stress and update material properties
			sig_eff = SolidMath::CalEffectiveStress(SolidMath::P_dev*sig_j);
			material_minkley->UpdateMinkleyProperties(sig_eff*material_minkley->GM0, e_pl_eff, Temperature);
			//evaluation of new residual
			material_minkley->CalViscoplasticResidual(dt,epsd_i,e_i,sig_j,eps_K_j,eps_K_t,eps_M_j,eps_M_t,eps_pl_j,eps_pl_t, \
													  e_pl_v,e_pl_v_t,e_pl_eff,e_pl_eff_t,lam,res_loc_p);
			if (Output)
			{
				CRFProcess* m_pcs  = PCSGet("DEFORMATION");
				ofstream Dum("local.txt", ios::app);
				Dum << aktueller_zeitschritt << " " << m_pcs->GetIteSteps() << " " << counter+1 << " " << res_loc_p.norm() <<  " plastic" << std::endl;
				Dum.close();
			}
		}
		//dGdE matrix and dsigdE matrix
		Eigen::Matrix<double,27,6> dGdE;
		Kelvin_to_Voigt_Strain(eps_pl_j,eps_pl_curr);
		//Calculate dGdE for time step
		material_minkley->CalEPdGdE(dGdE);
		//get dsigdE matrix
		ExtractConsistentTangent(K_loc_p,dGdE,true,dsigdE);
	}
	else
	{
		//dGdE matrix and dsigdE matrix
		Eigen::Matrix<double,18,6> dGdE;
		//Calculate dGdE for time step
		material_minkley->CaldGdE(dGdE);
		//get dsigdE matrix
		ExtractConsistentTangent(K_loc,dGdE,true,dsigdE);
	}
	if (counter == counter_max)
		std::cout << "WARNING: Maximum iteration number needed in LocalNewtonMinkley. Convergence not guaranteed." << std::endl;

	//add hydrostatic part to stress and tangent
	sig_j *= material_minkley->GM0;
	dsigdE *= material_minkley->GM0;

	//Sort into Consistent Tangent matrix for global Newton iteration and into standard OGS arrays
	Kelvin_to_Voigt_Stress(sig_j,stress_curr);
	Kelvin_to_Voigt_Strain(eps_K_j,eps_K_curr);
	Kelvin_to_Voigt_Strain(eps_M_j,eps_M_curr);
	//plastic strain dealt with further up
	for (size_t i=0; i<3; i++)
	{
		for (size_t j=0; j<3; j++){
			(*Consistent_Tangent)(i,j) = dsigdE(i,j);
			(*Consistent_Tangent)(i,j+3) = dsigdE(i,j+3)/sqrt(2.);//from local to global shear components (Kelvin to Voigt)
			(*Consistent_Tangent)(i+3,j) = dsigdE(i+3,j)/sqrt(2.);//from local to global shear components (Kelvin to Voigt)
			(*Consistent_Tangent)(i+3,j+3) = dsigdE(i+3,j+3)/2.;//from local to global shear components (Kelvin to Voigt)
		}
	}
}
/*************************************************************************
   FEMLib-Method: CSolidProperties::Calculate_Lame_Constant()
   Task: Get density
   Programing:
   08/2004 WW Implementation
**************************************************************************/
void CSolidProperties::Calculate_Lame_Constant()
{
	double nv = Poisson_Ratio();
	E = Youngs_Modulus(); // Constant at present
	// WX:1.2013. time dependet
	if (Time_Dependent_E_nv_mode > 0)
	{
		int valid = 1;
		switch (Time_Dependent_E_nv_mode)
		{
			case 1:
				nv *= GetCurveValue(Time_Dependent_E_nv_value[1], 0, aktuelle_zeit, &valid);
				E *= GetCurveValue(Time_Dependent_E_nv_value[0], 0, aktuelle_zeit, &valid);
				break;
			case 2:
				cout << "WARNING: for isotropic elasticity, please use TIME_DEPENDENT_YOUNG_POISSON mode=1" << endl;
				break;
			default:
				cout << "ERROR: not valid TIME_DEPENDENT_YOUNG_POISSON mode." << endl;
				break;
		}
	}
	Lambda = E * nv / ((1. + nv) * (1. - 2. * nv));
	G = 0.5 * E / (1. + nv);
	K = (3.0 * Lambda + 2.0 * G) / 3.0;
}

/*************************************************************************
   ROCKFLOW - Funktion: CSolidProperties::ElasticConsitutive

   Aufgabe:
   Fill the lastic constitutive

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
       const int Dimension  :  Dimension of the real space
       double *D_e          :  Elastic constitutive
   Ergebnis:
   - void -

   Programmaenderungen:
   11/2003   WW   Set plastic parameter

*************************************************************************/
void CSolidProperties::ElasticConsitutive(const int Dimension, Matrix* D_e) const
{
	(*D_e) = 0.0;
	(*D_e)(0, 0) = Lambda + 2 * G;
	(*D_e)(0, 1) = Lambda;
	(*D_e)(0, 2) = Lambda;

	(*D_e)(1, 0) = Lambda;
	(*D_e)(1, 1) = Lambda + 2 * G;
	(*D_e)(1, 2) = Lambda;

	(*D_e)(2, 0) = Lambda;
	(*D_e)(2, 1) = Lambda;
	(*D_e)(2, 2) = Lambda + 2 * G;

	(*D_e)(3, 3) = G;
	// Plane stress
	// plane stress, only for test
	//(*D_e)(0,0) = (1.0-Mu)*Lambda + 2 * G;
	//(*D_e)(0,1) = Lambda;

	//(*D_e)(1,0) = Lambda;
	// (*D_e)(1,1) = (1.0-Mu)*Lambda + 2 * G;
	// (*D_e)(3,3) = G;

	if (Dimension == 3)
	{
		(*D_e)(4, 4) = G;
		(*D_e)(5, 5) = G;
	}
}
/*************************************************************************
   ROCKFLOW - Funktion: CSolidProperties::ElasticConstitutiveTransverseIsotropic

   Aufgabe:
   Generate material matrix in case of transverse isotropic linear elasticity

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
       const int Dimension  :  Dimension of the real space
       double *D_e          :  Elastic constitutive
   Ergebnis:
   - void -

   Programmaenderungen:
   24.11.2009  UJG  first implementation

*************************************************************************/
void CSolidProperties::ElasticConstitutiveTransverseIsotropic(const int Dimension)
{
	double aii, aai, bii, bai, cii, cai;

	double ni = Poisson_Ratio();
	double Ei = (*data_Youngs)(0);
	double Ea = (*data_Youngs)(1);
	double nia = (*data_Youngs)(2);
	double Ga = (*data_Youngs)(3);

	// WX:01.2013, time dependet
	if (Time_Dependent_E_nv_mode > 0)
	{
		int valid = 1;
		switch (Time_Dependent_E_nv_mode)
		{
			case 1:
				cout << "WARNING: for transversely isotropic elasticity, please use TIME_DEPENDENT_YOUNG_POISSON mode=2"
				     << endl;
				break;
			case 2:
				Ei *= GetCurveValue(Time_Dependent_E_nv_value[0], 0, aktuelle_zeit, &valid);
				Ea *= GetCurveValue(Time_Dependent_E_nv_value[1], 0, aktuelle_zeit, &valid);
				ni *= GetCurveValue(Time_Dependent_E_nv_value[2], 0, aktuelle_zeit, &valid);
				nia *= GetCurveValue(Time_Dependent_E_nv_value[3], 0, aktuelle_zeit, &valid);
				Ga *= GetCurveValue(Time_Dependent_E_nv_value[4], 0, aktuelle_zeit, &valid);
				break;
			default:
				cout << "ERROR: not valid TIME_DEPENDENT_YOUNG_POISSON mode." << endl;
				break;
		}
	}

	double nai = nia * (Ea / Ei);
	double Disk = 1.0 + ni;

	Matrix* D_e = NULL;

	Disk *= (1.0 - ni - 2.0 * (nia * nai));
	Disk /= (Ei * Ei * Ea);

	aii = (1.0 - (nai * nia)) / (Ei * Ea * Disk);
	aai = (1.0 - (ni * ni)) / (Ei * Ei * Disk);
	bii = (ni + (nia * nai)) / (Ei * Ea * Disk);
	bai = (nia * (1.0 + ni)) / (Ei * Ei * Disk);
	cii = Ei / (2.0 * (1 + ni));
	cai = Ga;

	int size;
	size = Dimension * 2;

	switch (Youngs_mode)
	{
		case 10:
			int i, j, l, m;
			D_e = new Matrix(size, size);
			(*D_e) = 0.0;

			if (Dimension == 3)
			{
				(*D_e)(0, 0) = aii;
				(*D_e)(0, 1) = bii;
				(*D_e)(0, 2) = bai;

				(*D_e)(1, 0) = bii;
				(*D_e)(1, 1) = aii;
				(*D_e)(1, 2) = bai;

				(*D_e)(2, 0) = bai;
				(*D_e)(2, 1) = bai;
				(*D_e)(2, 2) = aai;

				(*D_e)(3, 3) = cii;
				(*D_e)(4, 4) = cai;
				(*D_e)(5, 5) = cai;
			}
			else
			{
				(*D_e)(0, 0) = aii;
				(*D_e)(0, 1) = bai;
				(*D_e)(0, 2) = bii;

				(*D_e)(1, 0) = bai;
				(*D_e)(1, 1) = aai;
				(*D_e)(1, 2) = bai;

				(*D_e)(2, 0) = bii;
				(*D_e)(2, 1) = bai;
				(*D_e)(2, 2) = aii;

				(*D_e)(3, 3) = cai;
			}

			D_tran = new Matrix(size, size);
			(*D_tran) = 0.;
			for (i = 0; i < size; i++)
				for (j = 0; j < size; j++)
					for (l = 0; l < size; l++)
						for (m = 0; m < size; m++)
							(*D_tran)(i, j) += (*Crotm)(l, i) * (*D_e)(l, m) * (*Crotm)(m, j);

			delete D_e;
			D_e = NULL;
			break;
		case 11:
			D_tran = new Matrix(size, size);
			(*D_tran)(0, 0) = aai;
			(*D_tran)(0, 1) = bai;
			(*D_tran)(0, 2) = bai;

			(*D_tran)(1, 0) = bai;
			(*D_tran)(1, 1) = aii;
			(*D_tran)(1, 2) = bii;

			(*D_tran)(2, 0) = bai;
			(*D_tran)(2, 1) = bii;
			(*D_tran)(2, 2) = aii;

			(*D_tran)(3, 3) = cai;

			if (Dimension == 3)
			{
				(*D_tran)(4, 4) = cai;
				(*D_tran)(5, 5) = cii;
			}
			break;
		case 12:
			D_tran = new Matrix(size, size);
			(*D_tran)(0, 0) = aii;
			(*D_tran)(0, 1) = bai;
			(*D_tran)(0, 2) = bii;

			(*D_tran)(1, 0) = bai;
			(*D_tran)(1, 1) = aai;
			(*D_tran)(1, 2) = bai;

			(*D_tran)(2, 0) = bii;
			(*D_tran)(2, 1) = bai;
			(*D_tran)(2, 2) = aii;

			(*D_tran)(3, 3) = cai;

			if (Dimension == 3)
			{
				(*D_tran)(4, 4) = cii;
				(*D_tran)(5, 5) = cai;
			}
			break;
		case 13:
			D_tran = new Matrix(size, size);
			(*D_tran)(0, 0) = aii;
			(*D_tran)(0, 1) = bii;
			(*D_tran)(0, 2) = bai;

			(*D_tran)(1, 0) = bii;
			(*D_tran)(1, 1) = aii;
			(*D_tran)(1, 2) = bai;

			(*D_tran)(2, 0) = bai;
			(*D_tran)(2, 1) = bai;
			(*D_tran)(2, 2) = aai;

			(*D_tran)(3, 3) = cii;

			if (Dimension == 3)
			{
				(*D_tran)(4, 4) = cai;
				(*D_tran)(5, 5) = cai;
			}
			break;
	}
}
/*************************************************************************
   ROCKFLOW - Funktion: CSolidProperties::CalculateTransformMatrixFromNormalVector

   Aufgabe:
   Generate rotation matrices for vectors and matrices based on a given normal vector

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
       const int Dimension  :  Dimension of the real space
       double *Crotv       :  Rotation matrix for vectors
       double *Crotm       :  Rotation matrix for matrices
   Ergebnis:
   - void -

   Programmaenderungen:
   25.11.2009  UJG  first implementation

*************************************************************************/
void CSolidProperties::CalculateTransformMatrixFromNormalVector(const int Dimension)
{
	if (!(Youngs_mode > 9 && Youngs_mode < 14)) // WW
		return;

	if (Youngs_mode != 10)
	{
		ElasticConstitutiveTransverseIsotropic(Dimension);
		return;
	}

	double e1[3] = {0.0}, e2[3] = {0.0}, e3[3] = {0.0};

	double nx = (*data_Youngs)(4);
	double ny = (*data_Youngs)(5);
	double nz = (*data_Youngs)(6);
	double Disk(1.0), t1(0.0), t2(0.0), t3(0.0), ax(1.0), ay(0.0), az(0.0);

	// rotation matrix for vectors: UJG 25.11.2009
	Matrix* Crotv = new Matrix(Dimension, Dimension);
	Crotm = new Matrix(Dimension * 2, Dimension * 2); // rotation matrix for matrices: UJG 25.11.2009

	if (Dimension == 3)
	{
		if (nz < -1.0)
		{
			t1 = 0.01745329252 * nx;
			t2 = 0.01745329252 * ny;
			e3[0] = sin(t1) * cos(t2);
			e3[1] = -sin(t2);
			e3[2] = cos(t1) * cos(t2);
		}
		else if (nz > 1.0)
		{
			t1 = 0.01745329252 * nx;
			t2 = 0.01745329252 * ny;
			e3[0] = sin(t2) * sin(t1);
			e3[1] = sin(t2) * cos(t1);
			e3[2] = cos(t2);
		}
		else
		{
			e3[0] = nx;
			e3[1] = ny;
			e3[2] = nz;
		}

		ax = e3[0];
		ay = e3[1];
		az = e3[2];

		if (ax == 0.0)
		{
			e1[0] = 1.0;
			e1[1] = 0.0;
			e1[2] = 0.0;

			e2[0] = 0.0;
			e2[1] = az;
			e2[2] = -ay;

			e3[0] = 0.0;
			e3[1] = ay;
			e3[2] = az;
		}
		else
		{
			Disk = 1.0 / sqrt(1.0 + ((az * az) / (ax * ax)));

			e1[0] = (az / ax) * Disk;
			e1[1] = 0.0;
			e1[2] = -Disk;

			t1 = ay * e1[2];
			t2 = az * e1[0] - ax * e1[2];
			t3 = -ay * e1[0];

			Disk = t1 * t1;
			Disk += t2 * t2;
			Disk += t3 * t3;

			Disk = 1.0 / sqrt(Disk);

			e2[0] = Disk * t1;
			e2[1] = Disk * t2;
			e2[2] = Disk * t3;
		}

		(*Crotv)(0, 0) = e1[0];
		(*Crotv)(0, 1) = e1[1];
		(*Crotv)(0, 2) = e1[2];

		(*Crotv)(1, 0) = e2[0];
		(*Crotv)(1, 1) = e2[1];
		(*Crotv)(1, 2) = e2[2];

		(*Crotv)(2, 0) = e3[0];
		(*Crotv)(2, 1) = e3[1];
		(*Crotv)(2, 2) = e3[2];

		(*Crotm)(0, 0) = (*Crotv)(0, 0) * (*Crotv)(0, 0);
		(*Crotm)(0, 1) = (*Crotv)(0, 1) * (*Crotv)(0, 1);
		(*Crotm)(0, 2) = (*Crotv)(0, 2) * (*Crotv)(0, 2);
		(*Crotm)(0, 3) = (*Crotv)(0, 0) * (*Crotv)(0, 1);
		(*Crotm)(0, 4) = (*Crotv)(0, 0) * (*Crotv)(0, 2);
		(*Crotm)(0, 5) = (*Crotv)(0, 1) * (*Crotv)(0, 2);

		(*Crotm)(1, 0) = (*Crotv)(1, 0) * (*Crotv)(1, 0);
		(*Crotm)(1, 1) = (*Crotv)(1, 1) * (*Crotv)(1, 1);
		(*Crotm)(1, 2) = (*Crotv)(1, 2) * (*Crotv)(1, 2);
		(*Crotm)(1, 3) = (*Crotv)(1, 0) * (*Crotv)(1, 1);
		(*Crotm)(1, 4) = (*Crotv)(1, 0) * (*Crotv)(1, 2);
		(*Crotm)(1, 5) = (*Crotv)(1, 1) * (*Crotv)(1, 2);

		(*Crotm)(2, 0) = (*Crotv)(2, 0) * (*Crotv)(2, 0);
		(*Crotm)(2, 1) = (*Crotv)(2, 1) * (*Crotv)(2, 1);
		(*Crotm)(2, 2) = (*Crotv)(2, 2) * (*Crotv)(2, 2);
		(*Crotm)(2, 3) = (*Crotv)(2, 0) * (*Crotv)(2, 1);
		(*Crotm)(2, 4) = (*Crotv)(2, 0) * (*Crotv)(2, 2);
		(*Crotm)(2, 5) = (*Crotv)(2, 1) * (*Crotv)(2, 2);

		(*Crotm)(3, 0) = 2.0 * (*Crotv)(0, 0) * (*Crotv)(1, 0);
		(*Crotm)(3, 1) = 2.0 * (*Crotv)(0, 1) * (*Crotv)(1, 1);
		(*Crotm)(3, 2) = 2.0 * (*Crotv)(0, 2) * (*Crotv)(1, 2);
		(*Crotm)(3, 3) = (*Crotv)(0, 0) * (*Crotv)(1, 1) + (*Crotv)(1, 0) * (*Crotv)(0, 1);
		(*Crotm)(3, 4) = (*Crotv)(0, 0) * (*Crotv)(1, 2) + (*Crotv)(0, 2) * (*Crotv)(1, 0);
		(*Crotm)(3, 5) = (*Crotv)(0, 1) * (*Crotv)(1, 2) + (*Crotv)(0, 2) * (*Crotv)(1, 1);

		(*Crotm)(4, 0) = 2.0 * (*Crotv)(0, 0) * (*Crotv)(2, 0);
		(*Crotm)(4, 1) = 2.0 * (*Crotv)(0, 1) * (*Crotv)(2, 1);
		(*Crotm)(4, 2) = 2.0 * (*Crotv)(0, 2) * (*Crotv)(2, 2);
		(*Crotm)(4, 3) = (*Crotv)(0, 0) * (*Crotv)(2, 1) + (*Crotv)(0, 1) * (*Crotv)(2, 0);
		(*Crotm)(4, 4) = (*Crotv)(0, 0) * (*Crotv)(2, 2) + (*Crotv)(2, 0) * (*Crotv)(0, 2);
		(*Crotm)(4, 5) = (*Crotv)(0, 1) * (*Crotv)(2, 2) + (*Crotv)(0, 2) * (*Crotv)(2, 1);

		(*Crotm)(5, 0) = 2.0 * (*Crotv)(1, 0) * (*Crotv)(2, 0);
		(*Crotm)(5, 1) = 2.0 * (*Crotv)(1, 1) * (*Crotv)(2, 1);
		(*Crotm)(5, 2) = 2.0 * (*Crotv)(1, 2) * (*Crotv)(2, 2);
		(*Crotm)(5, 3) = (*Crotv)(1, 0) * (*Crotv)(2, 1) + (*Crotv)(1, 1) * (*Crotv)(2, 0);
		(*Crotm)(5, 4) = (*Crotv)(1, 0) * (*Crotv)(2, 2) + (*Crotv)(1, 2) * (*Crotv)(2, 0);
		(*Crotm)(5, 5) = (*Crotv)(1, 1) * (*Crotv)(2, 2) + (*Crotv)(2, 1) * (*Crotv)(1, 2);
	}
	else
	{
		if (ny > 1.0)
		{
			e1[0] = cos(0.01745329252 * nx);
			e1[1] = sin(0.01745329252 * nx);

			e2[0] = -e1[1];
			e2[1] = e1[0];
		}
		else
		{
			e2[0] = nx;
			e2[1] = ny;

			e1[0] = ny;
			e1[1] = -nx;
		}

		(*Crotv)(0, 0) = e1[0];
		(*Crotv)(0, 1) = e1[1];

		(*Crotv)(1, 0) = e2[0];
		(*Crotv)(1, 1) = e2[1];

		(*Crotm)(0, 0) = (*Crotv)(0, 0) * (*Crotv)(0, 0);
		(*Crotm)(0, 1) = (*Crotv)(1, 0) * (*Crotv)(1, 0);
		(*Crotm)(0, 2) = 0.0;
		(*Crotm)(0, 3) = -(*Crotv)(0, 0) * (*Crotv)(1, 0);

		(*Crotm)(1, 0) = (*Crotv)(0, 1) * (*Crotv)(0, 1);
		(*Crotm)(1, 1) = (*Crotv)(1, 1) * (*Crotv)(1, 1);
		(*Crotm)(1, 2) = 0.0;
		(*Crotm)(1, 3) = -(*Crotv)(0, 1) * (*Crotv)(1, 1);

		(*Crotm)(2, 0) = 0.0;
		(*Crotm)(2, 1) = 0.0;
		(*Crotm)(2, 2) = 1.0;
		(*Crotm)(2, 3) = 0.0;

		(*Crotm)(3, 0) = -2.0 * (*Crotv)(0, 0) * (*Crotv)(0, 1);
		(*Crotm)(3, 1) = -2.0 * (*Crotv)(1, 0) * (*Crotv)(1, 1);
		(*Crotm)(3, 2) = 0.0;
		(*Crotm)(3, 3) = (*Crotv)(0, 0) * (*Crotv)(1, 1) + (*Crotv)(1, 0) * (*Crotv)(0, 1);
	}
	delete Crotv;
	Crotv = NULL;

	ElasticConstitutiveTransverseIsotropic(Dimension);
}
//
// WW. 09/02. Compute dilatancy for yield function
// For plane strain
double CSolidProperties::GetAngleCoefficent_DP(const double Angle)
{
	double val = 0.0;
	// Input as a coefficent
	if (Angle < 0.0)
		val = fabs(Angle);
	else // Input as a dialatant angle
	{
		double D_Angle = Angle * PI / 180.0;
		double sinA = sin(D_Angle);
		//     val = sinA/sqrt(9.0+4.0*sinA*sinA);
		//     val = 2.0*MSqrt2Over3*sinA/(3.0+sinA);
		val = 2.0 * MSqrt2Over3 * sinA / (3.0 + sinA); //(3.0-sinA)
		//     val = 2.0*MSqrt2Over3*sinA/(3.0-sinA);
	}
	return val;
}
// WW. 09/02. Cumpute yield coefficient, beta
// For plane strain
// al = 6.0*c*cos(a)/sqrt(3)/(3-sin(a))
double CSolidProperties::GetYieldCoefficent_DP(const double Angle)
{
	double val = 0.0;
	// Input as a coefficent
	if (Angle < 0.0 || Angle < MKleinsteZahl)
		val = 1.0;
	else // Input as a dialatant angle
	{
		double D_Angle = Angle * PI / 180.0;
		double sinA = sin(D_Angle);
		val = 6.0 * cos(D_Angle) / (3.0 + sinA); // (3.0-sinA)
		//     val = 2.0*sqrt(3.0)*cos(D_Angle)/sqrt(9.0+4.0*sinA*sinA);
	}
	return val;
}

void CSolidProperties::CalulateCoefficent_DP()
{
	Y0 = (*data_Plasticity)(0);
	Hard = (*data_Plasticity)(1);
	Al = GetAngleCoefficent_DP((*data_Plasticity)(2));
	Xi = GetAngleCoefficent_DP((*data_Plasticity)(3));
	BetaN = GetYieldCoefficent_DP((*data_Plasticity)(2));
	if (fabs(Al) < MKleinsteZahl && fabs(Xi) < MKleinsteZahl)
		BetaN = 1.0;
	BetaN *= sqrt(2.0 / 3.0);
	Hard_Loc = (*data_Plasticity)(4);
	tension = (*data_Plasticity)(5); // WX:
}

void CSolidProperties::CalculateCoefficent_MOHR(double ep, double scalar_comp, double scalar_tens) // WX:11.2010,
// 09.2011
{
	int valid = 1;
	double scalar_comp_tmp = scalar_comp;
	double scalar_tens_tmp = scalar_tens;
	double theta = (*data_Plasticity)(1) * PI / 180;
	double phi = (*data_Plasticity)(2) * PI / 180;
	Y0 = (*data_Plasticity)(0);
	tension = (*data_Plasticity)(3);

	if ((*data_Plasticity)(5) > 0 && (*data_Plasticity)(5) < 100)
		theta = GetCurveValue((int)(*data_Plasticity)(5), 0, ep, &valid) * PI / 180;
	if ((*data_Plasticity)(4) > 0 && (*data_Plasticity)(4) < 100)
		Y0 = GetCurveValue((int)(*data_Plasticity)(4), 0, ep, &valid);

	if (Plasticity_Bedding && ((scalar_comp_tmp != 0) || (scalar_tens_tmp != 0)))
	{
		// csn=GetCurveValue(bedding_uc_curve,0,scalar_comp_tmp,&valid);
		// theta = GetCurveValue(bedding_fric_curve,0,scalar_comp_tmp,&valid)*hPI/180;
		// tension=GetCurveValue(bedding_tens_curve,0,scalar_tens_tmp,&valid);
		csn = 0;
		tension = 0;
		int i;
		for (i = 0; i < bedding_uc_curve_order + 1; i++)
			csn += comp_para[i] * pow(scalar_comp_tmp, i);
		for (i = 0; i < bedding_tens_curve_order + 1; i++)
			tension += tens_para[i] * pow(scalar_tens_tmp, i);
		Ntheta = (1 + sin(theta)) / (1 - sin(theta));
		Nphi = (1 + sin(phi)) / (1 - sin(phi));
		if (tension < 0 || tension > Y0 / tan(theta))
			tension = Y0 / tan(theta);
	}
	else
	{
		Ntheta = (1 + sin(theta)) / (1 - sin(theta));
		Nphi = (1 + sin(phi)) / (1 - sin(phi));
		csn = 2 * Y0 * sqrt(Ntheta);
		if (tension < 0 || tension > Y0 / tan(theta))
			tension = Y0 / tan(theta);
	}
}

void CSolidProperties::CalculateCoefficent_HOEKBROWN() // WX: 02.2011
{
	HoekB_a = (*data_Plasticity)(0);
	HoekB_s = (*data_Plasticity)(1);
	HoekB_mb = (*data_Plasticity)(2);
	HoekB_sigci = (*data_Plasticity)(3);
	HoekB_tens = HoekB_s * HoekB_sigci / HoekB_mb;
	HoekB_cohe = HoekB_sigci * pow(HoekB_s, HoekB_a);
}

/**************************************************************************
   ROCKFLOW - Funktion: StressIntegrationDP

   Aufgabe:
   Computing the stresses at a point and return the plastical status of this
   point (Return mapping method)

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E :
     long index: Elementnummer
   double *TryStress                   :   Try incremental stresses as input
   New stresses as output
   const int GPi, const int GPj        :   Local indeces of Gauss Points

   const double G, double K            :   Shear modulus, (2*lambda+2*G)/3
   const double Al, const double Xi,
   const double Y0, const double BetaN :   Coefficient for Drucker-Prager model
   double* dPhi                        :   Plastic multiplier.
   double &jt                          :   enhanced parameters

   const int Update                    :   Indicator to store the stress or not.

   Ergebnis:
   - double* - Deviatoric effetive stresses, s11,s22,s12,s33

   Programmaenderungen:
   09/2002   WW  Erste Version
   02/2004   WW  Modify for the 3D case
   08/2004   WW  Set As a member of material class
   03/2007   WW  Multi-yield surface
**************************************************************************/
bool CSolidProperties::StressIntegrationDP(const int GPiGPj, const ElementValue_DM* ele_val, double* TryStress,
                                           double& dPhi, const int Update)
{
	int i = 0;
	double I1 = 0.0;
	double ep, ep0;
	double F = 0.0, F0 = 0.0; //, yl;
	//  double RF0 = 0.0;
	double sqrtJ2 = 0.0;
	double Beta = 0.0;
	double p3, normXi, Jac, err_corner = 0.0, fac = 0.0;
	int max_ite = 10;
	int ite = 0;
	// double dstrs[6];

	bool isLoop = true; // Used only to avoid warnings with .net
	bool ploading = false;

	const int Size = ele_val->Stress->Rows();
	// static double DevStress[6];

	int Dim = 2;
	if (Size > 4)
		Dim = 3;

	// Get the total effective plastic strain
	ep = (*ele_val->pStrain)(GPiGPj);

	for (i = 0; i < Size; i++)
	{
		//     dstrs[i] = TryStress[i];  // d_stress
		TryStress[i] += (*ele_val->Stress)(i, GPiGPj);
		devS[i] = TryStress[i];
	}
	// I_tr
	I1 = DeviatoricStress(devS);
	// s_tr
	sqrtJ2 = sqrt(TensorMutiplication2(devS, devS, Dim));
	//
	normXi = sqrtJ2;
	p3 = I1;
	/* If yield, compute plastic multiplier dPhi */
	dPhi = 0.0;
	F0 = sqrtJ2 + Al * I1;
	F = F0 - BetaN * (Y0 + Hard * ep);
	// if(pcs_deformation==1) F=-1.0;
	// yl = sqrt(TensorMutiplication2(devS, dstrs, Dim))/sqrtJ2;
	// yl += (dstrs[0]+dstrs[1]+dstrs[2])*Al;
	// if(yl<=0.0)
	//   F = -1.0;
	//
	if (F0 <= (*ele_val->y_surface)(GPiGPj)) // unloading
		F = -1.0;
	//
	if (F > 0.0 && (!PreLoad)) // in yield status
	{
		ploading = true;
		// err = 1.0e+5;
		ep0 = ep;

		// TEST
		// Check the corner region
		err_corner = 4.5 * Xi * K * normXi / G
		             + BetaN * (0.5 * Hard * sqrt(1.0 + 3.0 * Xi * Xi) * normXi / G + Y0 + Hard * ep) / Al;
		// Multi-surface
		if (p3 > err_corner)
		{
			// RF0 = F;
			dl2 = 0.0;
			while (isLoop)
			{
				ite++;
				// dl1
				dPhi = 0.5 * normXi / G;
				fac = sqrt(dPhi * dPhi + 3.0 * Xi * Xi * (dPhi + dl2) * (dPhi + dl2));
				Jac = 9.0 * Xi * K + 3.0 * Xi * Xi * BetaN * Hard * (dPhi + dl2) / (fac * Al);
				F = 9.0 * Xi * K * (dPhi + dl2) + BetaN * (Y0 + Hard * (ep0 + fac)) / Al - p3;
				dl2 -= F / Jac;
				if (fabs(F) < 1000.0 * Tolerance_Local_Newton)
					break;
				if (ite > max_ite)
					break;
			}
			ep = ep0 + fac;
			for (i = 0; i < 3; i++)
				TryStress[i] = I1 / 3.0 - 3.0 * (dPhi + dl2) * K * Xi;
			for (i = 3; i < Size; i++)
				TryStress[4] = 0.0;
		}
		else
		{
			// Local Newton-Raphson procedure
			// If non-perfect plasticity, the below line has to be change
			Jac = -2.0 * G - 9.0 * K * Al * Xi - BetaN * Hard * sqrt(1.0 + 3.0 * Xi * Xi);
			// RF0 = F;
			while (isLoop)
			{
				ite++;
				if (ite > max_ite)
					break;
				if (F < 0.0 || fabs(F) < 10.0 * Tolerance_Local_Newton)
					break;
				// if(err<TolLocalNewT) break;
				dPhi -= F / Jac;
				//
				p3 = I1 - 9.0 * dPhi * Xi * K;
				normXi = sqrtJ2 - 2.0 * G * dPhi;
				ep = ep0 + dPhi * sqrt(1.0 + 3.0 * Xi * Xi);
				F0 = normXi + Al * p3;
				F = F0 - BetaN * (Y0 + Hard * ep);
				/* Jac = fun(); if non-linear hardening is involved */
				// err = fabs(F)/RF0;
			}
			// update stress
			Beta = 1.0 - 2.0 * dPhi * G / sqrtJ2;
			for (i = 0; i < Size; i++)
				TryStress[i] = Beta * devS[i];
			for (i = 0; i < 3; i++)
				TryStress[i] += I1 / 3.0 - 3.0 * dPhi * K * Xi;
		}
	}
	else
	{
		//
		for (i = 0; i < Size; i++)
			TryStress[i] = devS[i];
		//
		for (i = 0; i < 3; i++)
			TryStress[i] += I1 / 3.0;
	}
	// Save the current stresses
	if (Update > 0)
	{
		if (dPhi > 0.0)
			(*ele_val->pStrain)(GPiGPj) = ep;
		(*ele_val->y_surface)(GPiGPj) = F0;
	}
	return ploading;
}

/**************************************************************************
   ROCKFLOW - Funktion: DirectStressIntegrationDP
   Computing the stresses at a point and return the plastical status of this
   point by direct integration.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E :
     long index: Elementnummer
     double *TryStress                   :   Try incremental stresses as input
                                             New stresses as output
     const int GPi, const int GPj        :   Local indeces of Gauss Points

   const double G, double K            :   Shear modulus, (2*lambda+2*G)/3
   const double Al, const double Xi,
   const double Y0, const double BetaN :   Coefficient for Drucker-Prager model
   double* dPhi                        :   Plastic multiplier.
   double &jt                          :   enhanced parameters

   const int Update                    :   Indicator to store the stress or not.

   return:
   - double* - Deviatoric effetive stresses, s11,s22,s12,s33

   Programmaenderungen:
   02/2006   WW  Erste Version
**************************************************************************/
bool CSolidProperties::DirectStressIntegrationDP(const int GPiGPj,
                                                 const ElementValue_DM* ele_val,
                                                 double* TryStress,
                                                 const int Update)
{
	int i = 0, m = 0, m_max = 100;
	double I1 = 0.0;
	double sy0, sy, ep, dlambda = 0.0;
	double F = 0.0, yy; //, yl;
	double R = 1.0;
	double sqrtJ2 = 0.0;
	double A_H = 0.0, domA = 0.0;
	double dstrs[6];
	bool ploading = false;
	const int Size = ele_val->Stress->Rows();
	// static double DevStress[6];

	int Dim = 2;
	if (Size > 4)
		Dim = 3;

	// Get the total effective plastic strain
	ep = (*ele_val->pStrain)(GPiGPj);

	for (i = 0; i < Size; i++)
	{
		dstrs[i] = TryStress[i]; // d_stress
		// stress_0
		TryStress[i] = (*ele_val->Stress)(i, GPiGPj);
		devS[i] = TryStress[i] + dstrs[i];
	}

	I1 = DeviatoricStress(devS);
	sqrtJ2 = sqrt(TensorMutiplication2(devS, devS, Dim));
	//
	sy = sqrtJ2 + Al * I1;
	yy = BetaN * (Y0 + Hard * ep);
	F = sy - yy;
	sy0 = (*ele_val->y_surface)(GPiGPj);
	// yl = sqrt(TensorMutiplication2(devS, dstrs, Dim))/sqrtJ2;
	// yl += (dstrs[0]+dstrs[1]+dstrs[2])*Al;
	// if(yl<=0.0)
	//  F = -1.0;
	if (sy <= sy0) // unloading
		F = -1.0;
	if (F > 0.0 && (!PreLoad)) // in yield status
	{
		if (ep < MKleinsteZahl) // Elastic in previous load step
			R = F / (sy - sy0);
		m = (int)(8.0 * F / yy) + 1;
		for (i = 0; i < Size; i++)
		{
			TryStress[i] += (1.0 - R) * dstrs[i];
			dstrs[i] *= R / (double)m;
		}
		if (m > m_max)
			m = m_max;
		// sub-inrement
		while (m > 0)
		{
			// Compute dlamda
			A_H = BetaN * Hard * sqrt(1 + 3.0 * Xi * Xi); // Hard: if it is not constant....
			for (i = 0; i < Size; i++)
				devS[i] = TryStress[i];
			I1 = DeviatoricStress(devS);
			sqrtJ2 = sqrt(TensorMutiplication2(devS, devS, Dim));
			for (i = 0; i < Size; i++)
			{
				devS[i] /= sqrtJ2;
				dFds[i] = devS[i];
			}
			for (i = 0; i < 3; i++)
				dFds[i] += Al;
			// dlambda
			dlambda = 0.0;
			domA = A_H + 2.0 * G + 9.0 * Al * Xi * K;
			for (i = 0; i < Size; i++)
				dlambda += dFds[i] * dstrs[i];
			dlambda /= domA;
			if (dlambda < 0.0)
				dlambda = 0.0;
			ep += dlambda * sqrt(1.0 + 3.0 * Xi * Xi);
			// Update stress
			for (i = 0; i < Size; i++)
				TryStress[i] += dstrs[i] - 2.0 * dlambda * G * devS[i];
			dlambda *= 3.0 * Xi * K;
			for (i = 0; i < 3; i++)
				TryStress[i] -= dlambda;
			m--;
		}
		for (i = 0; i < Size; i++)
			devS[i] = TryStress[i];
		I1 = DeviatoricStress(devS);
		sqrtJ2 = sqrt(TensorMutiplication2(devS, devS, Dim));
		sy = sqrtJ2 + Al * I1;
		yy = BetaN * (Y0 + Hard * ep);
		R = 1.0;
		if (sy > yy)
			R = yy / sy;
		for (i = 0; i < Size; i++)
		{
			TryStress[i] *= R;
			devS[i] *= R / sqrtJ2;
		}
		sy *= R;
		ploading = true;
	}
	else
		for (i = 0; i < Size; i++)
			TryStress[i] += dstrs[i];
	// Save the current stresses
	if (Update > 0)
	{
		(*ele_val->pStrain)(GPiGPj) = ep;
		(*ele_val->y_surface)(GPiGPj) = sy;
	}
	return ploading;
}

// WX
int CSolidProperties::DirectStressIntegrationDPwithTension(const int GPiGPj, Matrix* De, const ElementValue_DM* ele_val,
                                                           double* TryStress, const int Update, double& mm)
{
	int i = 0, j = 0; //, m=0, m_max=100;
	double I1 = 0.0;
	// WW double sy0,
	// WW double sy, ep, dlambda=0.0, tmpvalue=0., tmpvalue2=0., sqrtJ2I1;	//dTaun, dSign,
	double sy, ep, tmpvalue = 0., sqrtJ2I1; // dTaun, dSign,
	double F = 0.0, Ft = 0, yy; //, yl; Ft: tension failure WX: 11.08.2010
	// WW double R=1.0;
	double sqrtJ2 = 0.0;
	// WW double A_H = 0.0, domA=0.0;
	double dstrs[6];
	double dTempStr[6] = {0.}, dTempStr2[6] = {0.}, TmpStress0[6] = {0}, dStrainP[6] = {0.};
	// WW bool ploading = false;
	const int Size = ele_val->Stress->Rows();
	// WW double H=0;	//
	int failurestate = 0;

	// static double DevStress[6];

	int Dim = 2;
	if (Size > 4)
		Dim = 3;

	// Get the total effective plastic strain
	ep = (*ele_val->pStrain)(GPiGPj);

	for (i = 0; i < Size; i++)
	{
		dstrs[i] = TryStress[i]; // d_stress
		TryStress[i] = (*ele_val->Stress)(i, GPiGPj); // stress_0
		devS[i] = TryStress[i] + dstrs[i];
		tmpvalue += fabs(dstrs[i]);
	}

	I1 = DeviatoricStress(devS);
	sqrtJ2 = sqrt(TensorMutiplication2(devS, devS, Dim));
	Hard = 0; // WX:20.09.2010. no hardning in this model

	//
	sy = sqrtJ2 + Al * I1;
	yy = BetaN * (Y0 + Hard * ep);

	sqrtJ2I1 = yy - Al * 3 * tension;
	if (tension < 0)
		tension = 0;
	if (tension > (yy / Al / 3.0))
		tension = yy / Al / 3.0; // WX: tension strength must be positiv and has a max limit.

	F = sy - yy;
	Ft = I1 / 3 - tension; // WX: 11.08.2010  tension strength.
	// WW double Tau_P = (BetaN*Y0 - 3*Al*tension)/sqrt(2.0);
	// WW double Al_P = sqrt(1+4.5*Al*Al) - (3 * Al / (sqrt(2.0)));
	// WW  H = sqrtJ2/sqrt(2.0) - Tau_P - (Al_P * (I1/3. - tension));	//WX:16.09.2010.correct H.
	// WW sy0 = (*ele_val->y_surface)(GPiGPj);

	if (tmpvalue == 0)
		Ft = F = -1;

	if (F > 0.0 && (!PreLoad))
	{
		// return to Fs
		Matrix* tmpMatrix = new Matrix(Size, Size);
		for (i = 0; i < Size; i++)
		{
			D_dGds[i] = 0.0; // initialisation
			dFds[i] = devS[i] / (sqrtJ2 * sqrt(2.0));
			dGds[i] = dFds[i];
			if (i < 3)
			{
				dFds[i] += Al / sqrt(2.0);
				dGds[i] += Xi / sqrt(2.0);
			}
		}

		De->multi(dGds, D_dGds);
		tmpvalue = 0.;
		for (i = 0; i < Size; i++)
			tmpvalue += dFds[i] * D_dGds[i];

		for (i = 0; i < Size; i++)
			for (j = 0; j < Size; j++)
				(*tmpMatrix)(i, j) = D_dGds[i] * dFds[j];

		for (i = 0; i < Size; i++)
		{
			if (i < 3)
				dTempStr2[i] = TryStress[i] + dstrs[i] - yy / Al / 3.0;
			else
				dTempStr2[i] = TryStress[i] + dstrs[i];
		}
		tmpMatrix->multi(dTempStr2, dTempStr);

		for (i = 0; i < Size; i++)
		{
			dTempStr[i] /= tmpvalue;
			TryStress[i] += dstrs[i] - dTempStr[i];
		}

		tmpvalue = (TryStress[0] + TryStress[1] + TryStress[2]) / 3.0;
		if ((tmpvalue - tension) < MKleinsteZahl)
		{
			failurestate = 1; // shear
			delete tmpMatrix;
		}
		else
		{
			// return to Ft
			for (i = 0; i < Size; i++)
			{
				dFtds[i] = dGtds[i] = 1.0 / 3.0;
				if (i > 2)
					dFtds[i] = dGtds[i] = 0.0;
			}
			for (i = 0; i < 3; i++)
				TryStress[i] = (*ele_val->Stress)(i, GPiGPj) + dstrs[i];
			I1 = TryStress[0] + TryStress[1] + TryStress[2];
			dTempStr[0] = dTempStr[1] = dTempStr[2] = I1 / 3.0 - tension;
			for (i = 3; i < Size; i++)
				dTempStr[i] = 0.0;
			for (i = 0; i < Size; i++)
				TryStress[i] = (*ele_val->Stress)(i, GPiGPj) + dstrs[i] - dTempStr[i];

			for (i = 0; i < Size; i++)
				devS[i] = TryStress[i];

			I1 = DeviatoricStress(devS);
			sqrtJ2 = sqrt(TensorMutiplication2(devS, devS, Dim));

			if (sqrtJ2 <= sqrtJ2I1)
			{
				if ((dstrs[0] + dstrs[1] + dstrs[2]) != 0.0)
					mm = (3 * tension - ((*ele_val->Stress)(0, GPiGPj) + (*ele_val->Stress)(1, GPiGPj)
					                     + (*ele_val->Stress)(2, GPiGPj)))
					     / (dstrs[0] + dstrs[1] + dstrs[2]);
				else
					mm = 1.;
				// mm:  Ft(stress0 + mm*dstrs)=0	for Dep calculation
				failurestate = 2; // tensile failure
			}
			else
			{
				// return to corner between Fs and Ft
				double A, B, C;
				Matrix* tmpMatrix = new Matrix(Size, Size);
				for (i = 0; i < Size; i++) // initialize stress devstress
				{
					TryStress[i] = (*ele_val->Stress)(i, GPiGPj);
					devS[i] = TryStress[i] + dstrs[i];
					//  cout<<devS[i]<<"\n";
				}
				// cout<<"sigmaB_End"<<"\n";

				I1 = DeviatoricStress(devS);

				for (i = 0; i < Size; i++)
				{
					D_dGds[i] = 0.0; // initialization
					dFds[i] = devS[i] / (sqrtJ2 * sqrt(2.0));
					dGds[i] = dFds[i];
					if (i < 3)
					{
						dFds[i] += Al / sqrt(2.0);
						dGds[i] += Xi / sqrt(2.0);
					}
				}

				De->multi(dGds, D_dGds);
				tmpvalue = 0.;
				for (i = 0; i < Size; i++)
					tmpvalue += dFds[i] * D_dGds[i];

				for (i = 0; i < Size; i++)
					for (j = 0; j < Size; j++)
						(*tmpMatrix)(i, j) = D_dGds[i] * dFds[j];

				for (i = 0; i < Size; i++)
				{
					dTempStr[i] = 0.0;
					if (i < 3)
						dTempStr2[i] = TryStress[i] + dstrs[i] - yy / Al / 3.0;
					else
						dTempStr2[i] = TryStress[i] + dstrs[i];
				}
				tmpMatrix->multi(dTempStr2, dTempStr);
				for (i = 0; i < Size; i++)
				{
					dTempStr[i] /= tmpvalue;
					TryStress[i] += dstrs[i];
					// cout<<dTempStr[i]<<"||"<<TryStress[i]<<"\n";
				}
				//  cout<<"delta sigmaP__SigmaB__End"<<"\n";

				double dI1;
				I1 = DeviatoricStress(TryStress);
				dI1 = DeviatoricStress(dTempStr);

				A = 0.;
				for (i = 0; i < Size; i++)
				{
					if (i < 3)
						A += dTempStr[i] * dTempStr[i];
					else
						A += 2 * dTempStr[i] * dTempStr[i];
				}

				B = 0.;
				for (i = 0; i < Size; i++)
				{
					if (i < 3)
						B += -2 * TryStress[i] * dTempStr[i];
					else
						B += -4 * TryStress[i] * dTempStr[i];
				}

				C = 0.;
				for (i = 0; i < Size; i++)
				{
					if (i < 3)
						C += TryStress[i] * TryStress[i];
					else
						C += 2 * TryStress[i] * TryStress[i];
				}
				C -= sqrtJ2I1 * sqrtJ2I1;

				if ((B * B - 4 * A * C) < 0.0 || (A == 0))
					tmpvalue = 1.0;
				else
				{
					tmpvalue = (-B - sqrt(B * B - 4 * A * C)) / (2 * A);
					if (tmpvalue < 0)
					{
						tmpvalue = (-B + sqrt(B * B - 4 * A * C)) / (2 * A);
						if (tmpvalue > 1)
							tmpvalue = 1;
					}
				}

				double tmpStr[6] = {0.};
				for (i = 0; i < Size; i++)
				{
					if (i < 3)
						tmpStr[i] = (TryStress[i] + I1 / 3.0) - tmpvalue * (dTempStr[i] + dI1 / 3.0);
					else
						tmpStr[i] = TryStress[i] - tmpvalue * dTempStr[i];
				}
				double tmpI1 = 0.;
				tmpI1 = DeviatoricStress(tmpStr);

				for (i = 0; i < Size; i++)
				{
					if (i < 3)
						TryStress[i] = tmpStr[i] + tension;
					else
						TryStress[i] = tmpStr[i];
					// test
					// cout<<TryStress[i]<<"\n";
				}
				failurestate = 3; // corner

				if (Update < 1) // calculate mm for Dep
				{
					double m1, m2, n1, n2, dI2;

					m1 = tmpvalue;

					for (i = 0; i < Size; i++) // initialize stress devstress
					{
						TmpStress0[i] = (*ele_val->Stress)(i, GPiGPj);
						dTempStr[i] = dstrs[i];
						dTempStr2[i] = TmpStress0[i] + dstrs[i]; // stress0+dstrs
						//  cout<<devS[i]<<"\n";
					}
					I1 = DeviatoricStress(TmpStress0);
					dI1 = DeviatoricStress(dTempStr);
					dI2 = DeviatoricStress(dTempStr2);

					if (tmpI1 > 3 * tension)
						n1 = (tmpI1 - 3 * tension)
						     / (dI2 - 3 * tension); // n1 =  (I1(sigtmp)-I1(tension))/ (I1(SigB)-I1(tension))
					else
						n1 = 0;
					if (dI1 != 0.0)
					{
						n2 = (3 * tension - I1) / dI1;
						if (n2 < 0)
							n2 = 0.;
					}
					else
						n2 = 1.;

					A = 0.;
					for (i = 0; i < Size; i++)
					{
						if (i < 3)
							A -= dTempStr[i] * dTempStr[i];
						else
							A -= 2 * dTempStr[i] * dTempStr[i];
					}
					A += (Al * dI1) * (Al * dI1);

					B = 0.;
					for (i = 0; i < Size; i++)
					{
						if (i < 3)
							B -= 2 * TmpStress0[i] * dTempStr[i];
						else
							B -= 4 * TmpStress0[i] * dTempStr[i];
					}
					B += 2 * (Al * I1 - yy) * Al * dI1;

					C = 0.;
					for (i = 0; i < Size; i++)
					{
						if (i < 3)
							C -= TmpStress0[i] * TmpStress0[i];
						else
							C -= 2 * TmpStress0[i] * TmpStress0[i];
					}
					C += (Al * I1 - yy) * (Al * I1 - yy);

					if (((B * B - 4 * A * C) < 0) || (A == 0))
						m2 = 1;
					else
					{
						m2 = (-B - sqrt(B * B - 4 * A * C)) / (2 * A);
						if (m2 < 0)
						{
							m2 = (-B + sqrt(B * B - 4 * A * C)) / (2 * A);
							if (m2 > 1)
								m2 = 0;
						}
					}

					Matrix* Dp_shear = new Matrix(Size, Size);
					Matrix* Dp_tension = new Matrix(Size, Size);
					Matrix* M_ds = new Matrix(Size, Size); // dGsds*dFsdsT

					tmpvalue = 0.; // dFsdsT*De*dGsds
					// WW			  tmpvalue2 =0.;		//dFtdsT*De*dGtds

					for (i = 0; i < 3; i++)
						for (j = 0; j < 3; j++)
							(*Dp_tension)(i, j) = K;

					for (i = 0; i < Size; i++)
						tmpvalue += dFds[i] * D_dGds[i];

					for (i = 0; i < Size; i++)
						for (j = 0; j < Size; j++)
							(*M_ds)(i, j) = dGds[i] * dFds[j];
					De->multi(*M_ds, *De, *Dp_shear);

					*Dp_shear /= tmpvalue;

					ConstitutiveMatrix->resize(Size, Size);

					for (i = 0; i < Size; i++)
						for (j = 0; j < Size; j++)
							(*ConstitutiveMatrix)(i, j)
							    = (*De)(i, j) - (1 - m2) * m1 * (*Dp_shear)(i, j) - (1 - n2) * n1 * (*Dp_tension)(i, j);

					delete Dp_shear;
					delete Dp_tension;
					delete M_ds;
				}

				delete tmpMatrix;
			}
		}
	}
	else if (Ft > 0) // muss be tensile
	{
		for (i = 0; i < Size; i++)
		{
			dFtds[i] = dGtds[i] = 1.0 / 3.0;
			if (i > 2)
				dFtds[i] = dGtds[i] = 0.0;
		}

		for (i = 0; i < Size; i++)
			TryStress[i] = (*ele_val->Stress)(i, GPiGPj) + dstrs[i];
		I1 = TryStress[0] + TryStress[1] + TryStress[2];

		dTempStr[0] = dTempStr[1] = dTempStr[2] = I1 / 3.0 - tension;
		for (i = 3; i < Size; i++)
			dTempStr[i] = 0.0;

		for (i = 0; i < Size; i++)
			TryStress[i] -= dTempStr[i];
		if ((dstrs[0] + dstrs[1] + dstrs[2]) != 0.0)
			mm = (3 * tension
			      - ((*ele_val->Stress)(0, GPiGPj) + (*ele_val->Stress)(1, GPiGPj) + (*ele_val->Stress)(2, GPiGPj)))
			     / (dstrs[0] + dstrs[1] + dstrs[2]);
		else
			mm = 1.;
		failurestate = 2; // tensile failure
	}

	if (F <= 0 && Ft <= 0)
	{
		for (i = 0; i < Size; i++)
			TryStress[i] += dstrs[i];
		failurestate = 0;
	}

	// update ep
	if (failurestate != 0)
	{
		for (i = 0; i < Size; i++)
		{
			dTempStr[i] = (*ele_val->Stress)(i, GPiGPj) + dstrs[i] - TryStress[i]; // d plas stress
			dStrainP[i] = 0.;
		}

		Matrix* invDe = new Matrix(Size, Size);
		Cal_Inv_Matrix(Size, De, invDe);

		invDe->multi(dTempStr, dStrainP);

		ep += sqrt(TensorMutiplication2(dStrainP, dStrainP, Dim)) * sqrt(2.0 / 3.0);
		delete invDe;
	}
	// Save the current stresses
	if (Update > 0)
	{
		(*ele_val->pStrain)(GPiGPj) = ep;
		(*ele_val->y_surface)(GPiGPj) = sy;
	}
	// return ploading;
	return failurestate;
}

/**************************************************************************
   ROCKFLOW - Funktion: CSolidProperties::ConsistentTangentialDP

   Local assembly of elasto-plastic tangential matrix C^ep
   (Drucker-Prager model)

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E
   double *Dep            : Consistent tangential matrix
   const Dim              : Space dimension

   Programmaenderungen:
   10/2006   WW  Erste Version
   02/2006   WW  programmed
**************************************************************************/
void CSolidProperties::TangentialDP(Matrix* Dep)
{
	int i, j, Size;
	double domA;
	//
	Size = Dep->Rows();
	domA = BetaN * Hard * sqrt(1 + 3.0 * Xi * Xi); // Hard: if it is not constant....
	//
	for (i = 0; i < Size; i++)
	{
		D_dFds[i] = 2.0 * G * devS[i];
		D_dGds[i] = 2.0 * G * devS[i];
	}
	for (i = 0; i < 3; i++)
	{
		D_dFds[i] += 3.0 * Al * K;
		D_dGds[i] += 3.0 * Xi * K;
	}
	//
	domA += 2.0 * G + 9.0 * Al * Xi * K;
	//
	for (i = 0; i < Size; i++)
		for (j = 0; j < Size; j++)
			(*Dep)(i, j) -= D_dGds[i] * D_dFds[j] / domA;

	// Dep->Write();
}

// WX: return to shear
void CSolidProperties::TangentialDP2(Matrix* Dep)
{
	int i, j, Size;
	double sqrtJ2;
	//
	Size = Dep->Rows();

	double dTemp2 = 0;
	double dTempStr2[6] = {0};
	Matrix D_temp(Size, Size);
	Matrix D_temp2(Size, Size);
	//
	int Dim = 2;
	if (Size > 4)
		Dim = 3;
	sqrtJ2 = sqrt(TensorMutiplication2(devS, devS, Dim));

	for (i = 0; i < Size; i++)
	{
		D_dFds[i] = 0.5 * devS[i] / sqrtJ2;
		D_dGds[i] = 0.5 * devS[i] / sqrtJ2;
	}
	for (i = 0; i < 3; i++)
	{
		D_dFds[i] += Al;
		D_dGds[i] += Xi;
	}

	Dep->multi(D_dGds, dTempStr2);
	for (i = 0; i < Size; i++)
		dTemp2 += D_dFds[i] * dTempStr2[i];

	for (i = 0; i < Size; i++)
		for (j = 0; j < Size; j++)
			D_temp(i, j) = D_dGds[i] * D_dFds[j];

	Dep->multi(D_temp, *Dep, D_temp2);

	//
	for (i = 0; i < Size; i++)
		for (j = 0; j < Size; j++)
			(*Dep)(i, j) -= D_temp2(i, j) / dTemp2;
}

// WX: return to tension
void CSolidProperties::TangentialDPwithTension(Matrix* Dep, double mm)
{
	int i, j, Size;
	//  double domA;
	Size = Dep->Rows();
	// double dTemp2 = 0;
	// double dTempStr3[6] = {0};
	Matrix* D_temp = new Matrix(Size, Size);
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			(*D_temp)(i, j) = K;
	// Matrix D_temp2(Size, Size);
	// mm=0.;
	/*
	   Dep->multi(dGtds, dTempStr3);	// dTemp = De * dGt/ds

	   for(i=0; i<Size; i++)
	   dTemp2 += dFtds[i] * dTempStr3[i];	// dFt/ds^T * De * dGt/ds

	   double dTemp3 = 1.0/dTemp2;	// 1/ dFt/ds^T * De * dGt/ds

	   for(int i=0; i<3; i++)
	   for(int j=0; j<3; j++)
	       D_temp(i,j)=1.0/9.0;

	   Dep->multi(D_temp, *Dep, D_temp2);
	 */

	for (i = 0; i < Size; i++)
		for (j = 0; j < Size; j++)
			(*Dep)(i, j) = mm * (*Dep)(i, j) + (1 - mm) * ((*Dep)(i, j) - (*D_temp)(i, j));
	delete D_temp;
}

// WX: return to corner
void CSolidProperties::TangentialDPwithTensionCorner(Matrix* Dep, double /*mm*/)
{
	// return to corner
	*Dep = *ConstitutiveMatrix;
	// Dep->Write();
}

/*******************************************************
WX: Mohr coulomb,
directe stress integration, also for aniso.
*******************************************************/
int CSolidProperties::StressIntegrationMOHR_Aniso(const int GPiGPj, const ElementValue_DM* ele_val, double* TryStress,
                                                  const int Update, Matrix* Dep)
{
	int i, j, counter;
	int yield = 0;
	int Dim = 2;
	int Size = std::min(std::size_t(6ul), ele_val->Stress->Rows());
	double normdstr = 0., TmpValue1, TmpValue2;
	double LodeAngle, I1, J2, J3, sqrtJ2;
	double AnisoParaComp, AnisoParaTens = 0.0;
	double shearsurf, tensionsurf, ep, dlamda = 0.0 /*, dlamda_0, dlamda_1*/; //, ddlamda, Jacob;

	// initialize all vectors
	double dstrs[6] = {0.}, TryStress_0[6] = {0.}, TryStr_buff[6] = {0.}, TmpStress[6] = {0.};
	double tmp_prin_str[6] = {0.}, tmp_prin_dir[9] = {0.}, TmpPrinStrComp[6] = {0.}, TmpPrinStrTens[6] = {0.};
	double TmpStrComp[6] = {0.}, TmpStrTens[6] = {0.}, devStr[6] = {0.}, dstrainP[6] = {0.}, dstressP[6] = {0.};
	double sqrt3 = sqrt(3.0);
	bool lode_out_range = false;
	double TryStress_local[6]={0,}, dstrs_local[6]={0.};

	*TmpDe = (0.);
	*Inv_De = (0.);
	*TMatrix = (0.);
	*TmpMatrix = (0.);
	*TmpMatrix2 = (0.);

	*TmpDe = *Dep;

	// ConstitutiveMatrix->resize(Size,Size);		//in head already defined, and is used for later as global variable

	*ConstitutiveMatrix = (0.);

	ep = (*ele_val->pStrain)(GPiGPj); // get eff plas strain

	if (Size > 4)
		Dim = 3;

	for (i = 0; i < Size; i++)
	{
		dstrs[i] = TryStress[i]; // d_stress
		TryStress[i] = (*ele_val->Stress)(i, GPiGPj); // stress_0
		TryStr_buff[i] = TryStress[i] + dstrs[i];
		TryStress_0[i] = TryStr_buff[i];
		normdstr += fabs(dstrs[i]);
	}
	if (normdstr < MKleinsteZahl)
		return 0;
	// if(Size==4)
	//	devStr[4]=devStr[5]=0.;

	if (!Plasticity_Bedding)
	{
		CalculateCoefficent_MOHR(ep, 0, 0);
		for (i = 0; i < Size; i++)
		{
			TmpStress[i] = TryStr_buff[i];
			dstrs_local[i] = dstrs[i];
			TryStress_local[i] = TryStress[i];
		}
	}
	else
	{
		TransMicroStru_TInv->multi(TryStr_buff, TmpStress);
		TransMicroStru_TInv->multi(dstrs, dstrs_local);
		TransMicroStru_TInv->multi(TryStress, TryStress_local);
		for (i = 0; i < Size; i++)
		{
			TryStr_buff[i] = TmpStress[i]; // all stress in local cordinate sys.
		}

		// TEST
		// TmpStress[0] = -3099268.72638144;
		// TmpStress[1] = -2125514.62848589;
		// TmpStress[2] = -9114414.6324919;
		// TmpStress[3] = 26098.3556044748;
		// TmpStress[4] = -4635944.67744482;
		// TmpStress[5] = 46718.3740402566;
		// for(i=0;i<Size;i++)
		//{
		//	TryStr_buff[i]=TmpStress[i];
		//}
		//

		*TmpMatrix = (0.);
		*TmpMatrix2 = (0.);
		CalPrinStrDir(TmpStress, tmp_prin_str, tmp_prin_dir, Dim);
		CalTransMatrixA(tmp_prin_dir, TmpMatrix, Size);
		TmpMatrix->GetTranspose(*TmpMatrix2);

		for (i = 0; i < Size; i++)
		{
			TmpStrComp[i] = 0.;
			TmpStrTens[i] = 0.;
			if (tmp_prin_str[i] < MKleinsteZahl)
				TmpPrinStrComp[i] = tmp_prin_str[i];
			else
				TmpPrinStrTens[i] = tmp_prin_str[i];
		}
		TmpMatrix2->multi(TmpPrinStrComp, TmpStrComp);
		TmpMatrix2->multi(TmpPrinStrTens, TmpStrTens);
		AnisoParaComp = CalAnisoPara(TmpStrComp, MicroStruTensor);
		AnisoParaTens = CalAnisoPara(TmpStrTens, MicroStruTensor);
		CalculateCoefficent_MOHR(ep, AnisoParaComp, AnisoParaTens);
	}

	double Kronecker[6] = {1, 1, 1, 0, 0, 0};

	I1 = TmpStress[0] + TmpStress[1] + TmpStress[2];
	for (i = 0; i < Size; i++)
		devStr[i] = TmpStress[i] - I1 / 3.0 * Kronecker[i];
	J2 = 0.5 * TensorMutiplication2(devStr, devStr, Dim);
	sqrtJ2 = sqrt(J2);
	J3 = TensorMutiplication3(devStr, devStr, devStr, Dim);

	if (J2 == 0)
	{
		LodeAngle = 0; // avoid error
		lode_out_range = true;
	}
	else
	{
		if ((-3 * sqrt3 / 2. * J3 / pow(sqrtJ2, 3)) > 1 - 1e-6)
		{
			LodeAngle = PI / 6.;
			lode_out_range = true;
		}
		else if ((-3 * sqrt3 / 2. * J3 / pow(sqrtJ2, 3)) < -1 + 1e-6)
		{
			LodeAngle = -PI / 6.;
			lode_out_range = true;
		}
		else
		{
			LodeAngle = asin(-3 * sqrt3 / 2. * J3 / pow(sqrtJ2, 3)) / 3.;
			lode_out_range = false;
		}
	}

	shearsurf = Ntheta * (I1 / 3. + 2 / sqrt3 * sqrtJ2 * sin(LodeAngle + 2 / 3. * PI))
	            - (I1 / 3. + 2 / sqrt3 * sqrtJ2 * sin(LodeAngle - 2 / 3. * PI)) - csn;
	tensionsurf = (I1 / 3. + 2 / sqrt3 * sqrtJ2 * sin(LodeAngle + 2 / 3. * PI)) - tension;

	if (shearsurf > 1e-5 || tensionsurf > 1e-5)
	{
		yield = 1;
		double dsqrtJ2_dsig[6] = {0.}, dlode_dsig[6] = {0.};
		double dAniso_dsig_tens[6] = {0.}, dAniso_dsig_comp[6] = {0.};
		double dcsn_dsig[6] = {0.}, dtens_dsig[6] = {0.};
		double dsig_dlamda[6] = {0.}, /*dsig_dlamda_k1[6]={0.},*/ dfs_dsig[6] = {0.}, dft_dsig[6] = {0.},
		       dgs_dsig[6] = {0.}, dgt_dsig[6] = {0.};
		bool first_step;
		double shearsurf_k1, tensionsurf_k1 /*, local_damp, tmp_pos = 1., tmp_neg=1.*/; // Newton downhill
		double lamda_pos = 0., lamda_neg = 0., shearsurf_pos, shearsurf_neg, tensionsurf_pos, tensionsurf_neg;
		if (tensionsurf > MKleinsteZahl) // if tension
		{
			// dlamda = 0, dlamda_1=1;
			counter = 0;
			first_step = true;
			tensionsurf_k1 = tensionsurf;
			for (i = 0; i < Size; i++)
				TryStr_buff[i] = TmpStress[i];
			while (1)
			{
				for (i = 0; i < Size; i++)
					dsqrtJ2_dsig[i] = devStr[i] / 2. / sqrtJ2;

				if (lode_out_range || first_step)
				{
					for (i = 0; i < Size; i++)
						dlode_dsig[i] = 0;
				}
				else
				{
					dlode_dsig[0] = sqrt3 / (2. * MathLib::fastpow(sqrtJ2, 3) * cos(3. * LodeAngle))
					                * (3. * J3 / 2. / J2 * devStr[0]
					                   - (devStr[0] * devStr[0] + devStr[3] * devStr[3] + devStr[4] * devStr[4])
					                   + 2. * J2 * Kronecker[0] / 3.0);

					dlode_dsig[1] = sqrt3 / (2. * MathLib::fastpow(sqrtJ2, 3) * cos(3. * LodeAngle))
					                * (3. * J3 / 2. / J2 * devStr[1]
					                   - (devStr[3] * devStr[3] + devStr[1] * devStr[1] + devStr[5] * devStr[5])
					                   + 2. * J2 * Kronecker[1] / 3.0);

					dlode_dsig[2] = sqrt3 / (2. * MathLib::fastpow(sqrtJ2, 3) * cos(3. * LodeAngle))
					                * (3. * J3 / 2. / J2 * devStr[2]
					                   - (devStr[4] * devStr[4] + devStr[5] * devStr[5] + devStr[2] * devStr[2])
					                   + 2. * J2 * Kronecker[2] / 3.0);

					dlode_dsig[3] = sqrt3 / (2. * MathLib::fastpow(sqrtJ2, 3) * cos(3. * LodeAngle))
					                * (3. * J3 / 2. / J2 * devStr[3]
					                   - (devStr[0] * devStr[3] + devStr[3] * devStr[1] + devStr[4] * devStr[5])
					                   + 2. * J2 * Kronecker[3] / 3.0);

					if (Dim == 3)
					{
						dlode_dsig[4] = sqrt3 / (2. * MathLib::fastpow(sqrtJ2, 3) * cos(3. * LodeAngle))
						                * (3. * J3 / 2. / J2 * devStr[4]
						                   - (devStr[0] * devStr[4] + devStr[3] * devStr[5] + devStr[4] * devStr[2])
						                   + 2. * J2 * Kronecker[4] / 3.0);

						dlode_dsig[5] = sqrt3 / (2. * MathLib::fastpow(sqrtJ2, 3) * cos(3. * LodeAngle))
						                * (3. * J3 / 2. / J2 * devStr[5]
						                   - (devStr[3] * devStr[4] + devStr[1] * devStr[5] + devStr[5] * devStr[2])
						                   + 2. * J2 * Kronecker[5] / 3.0);
					}
				}

				if (Plasticity_Bedding)
				{
					double a_pk_sig_pq_sig_kq = 0., sig_pq_sig_pq = 0., dTens_daPara = 0.;
					double a_ki_sig_kj[6] = {0.};
					int p, k, q;
					*TmpMatrix = (0.);
					*TmpMatrix2 = (0.);
					for (i = 0; i < 3; i++)
					{
						for (j = 0; j < 3; j++)
						{
							if (i == j)
							{
								(*TmpMatrix)(i, j) = MicroStruTensor[i];
								//(*TmpMatrix2)(i,j) = TmpStress[i];
								(*TmpMatrix2)(i, j) = TmpStrTens[i];
							}
							else
								//(*TmpMatrix2)(i,j) = TmpStress[i+j+2];
								(*TmpMatrix2)(i, j) = TmpStrTens[i + j + 2];
						}
					}

					for (p = 0; p < 3; p++) // aki_sigkj, sigpq_sigpa = sigmn_sigmn, apk_sigpq_sigkq
					{
						a_ki_sig_kj[0] += (*TmpMatrix)(p, 0) * (*TmpMatrix2)(p, 0);
						a_ki_sig_kj[1] += (*TmpMatrix)(p, 1) * (*TmpMatrix2)(p, 1);
						a_ki_sig_kj[2] += (*TmpMatrix)(p, 2) * (*TmpMatrix2)(p, 2);
						a_ki_sig_kj[3] += (*TmpMatrix)(p, 0) * (*TmpMatrix2)(p, 1);
						a_ki_sig_kj[4] += (*TmpMatrix)(p, 0) * (*TmpMatrix2)(p, 2);
						a_ki_sig_kj[5] += (*TmpMatrix)(p, 1) * (*TmpMatrix2)(p, 2);
						for (q = 0; q < 3; q++)
						{
							sig_pq_sig_pq += (*TmpMatrix2)(p, q) * (*TmpMatrix2)(p, q);
							for (k = 0; k < 3; k++)
								a_pk_sig_pq_sig_kq += (*TmpMatrix)(p, k) * (*TmpMatrix2)(p, q) * (*TmpMatrix2)(k, q);
						}
					}
					for (i = 0; i < Size; i++)
						dAniso_dsig_tens[i] = 2 * (a_ki_sig_kj[i] * sig_pq_sig_pq - a_pk_sig_pq_sig_kq * TmpStrTens[i])
						                      / (sig_pq_sig_pq * sig_pq_sig_pq); // instead TmpStress[i] with TmpStrTens

					for (i = 1; i < bedding_tens_curve_order + 1; i++)
						dTens_daPara += i * tens_para[i] * MathLib::fastpow(AnisoParaTens, (i - 1)); // dcsn/deta
					for (i = 0; i < Size; i++)
						dtens_dsig[i] = dTens_daPara * dAniso_dsig_tens[i]; // dcsn/dsig=dcsn/deta * deta/dsig
				} // end if plasticity bedding
				if (first_step)
				{
					for (i = 0; i < Size; i++)
					{
						dft_dsig[i] = 1 / 3. * Kronecker[i] + 2 / sqrt3 * sin(LodeAngle + 2 / 3. * PI) * dsqrtJ2_dsig[i]
						              + 2 / sqrt3 * sqrtJ2 * dlode_dsig[i] * cos(LodeAngle + 2 / 3. * PI)
						              - dtens_dsig[i];
						dgt_dsig[i] = 1 / 3. * Kronecker[i] + 2 / sqrt3 * sin(LodeAngle + 2 / 3. * PI) * dsqrtJ2_dsig[i]
						              + 2 / sqrt3 * sqrtJ2 * dlode_dsig[i] * cos(LodeAngle + 2 / 3. * PI);
						//-dcsn_dsig[i];
					}
				}

				//
				if (counter == 0)
				{
					if (tensionsurf_k1 > 0)
					{
						lamda_pos = dlamda;
						TmpValue1 = 0.;
						if (first_step)
						{
							for (i = 0; i < Size; i++)
								dsig_dlamda[i] = 0.;
							Dep->multi(dgt_dsig, dsig_dlamda, -1);
							lamda_pos = 0.;
							tensionsurf_pos = tensionsurf;
							first_step = false;
						}
						for (i = 0; i < Size; i++)
						{
							TmpValue1 += dft_dsig[i] * (-1) * dsig_dlamda[i]; //(dF/dstr)T De dG/dstr
							//	TmpValue2 += dft_dsig[i]*dstrs[i];
						}
						if (TmpValue1 == 0)
							dlamda = 0.;
						else
							dlamda = tensionsurf / TmpValue1;
						tensionsurf_pos = tensionsurf_k1;
						dlamda += lamda_pos;
					}
					else
					{
						tensionsurf_neg = tensionsurf_k1;
						lamda_neg = dlamda;
						counter++;
					}
				}
				else
				{
					counter++;
					if (fabs(tensionsurf_k1) < 1e-3)
					{
						// for(i=0;i<Size;i++)
						//	dsig_dlamda[i]=0.;
						// Dep->multi(dgt_dsig,dsig_dlamda,-1);
						break;
					}
					else if (tensionsurf_k1 < 0)
					{
						// for(i=0;i<Size;i++)
						//	dsig_dlamda[i]=0.;
						// Dep->multi(dgt_dsig,dsig_dlamda,-1);
						shearsurf_neg = tensionsurf_k1;
						lamda_neg = dlamda;
						dlamda = lamda_pos
						         - tensionsurf_pos * (lamda_neg - lamda_pos) / (tensionsurf_neg - tensionsurf_pos);
						if (fabs(lamda_neg - lamda_pos) < 1e-8)
							break;
						// dlamda = (lamda_pos+lamda_neg)/2.;
					}
					else
					{
						// for(i=0;i<Size;i++)
						//	dsig_dlamda[i]=0.;
						// Dep->multi(dgt_dsig,dsig_dlamda,-1);
						lamda_pos = dlamda;
						tensionsurf_pos = tensionsurf_k1;
						dlamda = lamda_pos
						         - tensionsurf_pos * (lamda_neg - lamda_pos) / (tensionsurf_neg - tensionsurf_pos);
						if (fabs(lamda_neg - lamda_pos) < 1e-8)
							break;
						// dlamda = (lamda_pos+lamda_neg)/2.;
					}
				}

				for (i = 0; i < Size; i++)
					TmpStress[i] = TryStr_buff[i] + dlamda * dsig_dlamda[i]; // sig(n+1)=sig(try)-dlambda*De*dG/dsig

				if (!Plasticity_Bedding)
				{
					CalculateCoefficent_MOHR(ep, 0, 0);
				}
				else
				{
					*TmpMatrix = (0.);
					*TmpMatrix2 = (0.);
					CalPrinStrDir(TmpStress, tmp_prin_str, tmp_prin_dir, Dim);
					CalTransMatrixA(tmp_prin_dir, TmpMatrix, Size);
					TmpMatrix->GetTranspose(*TmpMatrix2);

					for (i = 0; i < Size; i++)
					{
						TmpStrComp[i] = 0.;
						TmpStrTens[i] = 0.;
						if (tmp_prin_str[i] < MKleinsteZahl)
							TmpPrinStrComp[i] = tmp_prin_str[i];
						else
							TmpPrinStrTens[i] = tmp_prin_str[i];
					}
					TmpMatrix2->multi(TmpPrinStrComp, TmpStrComp);
					TmpMatrix2->multi(TmpPrinStrTens, TmpStrTens);
					AnisoParaComp = CalAnisoPara(TmpStrComp, MicroStruTensor);
					AnisoParaTens = CalAnisoPara(TmpStrTens, MicroStruTensor);
					CalculateCoefficent_MOHR(ep, AnisoParaComp, AnisoParaTens);
				}

				I1 = TmpStress[0] + TmpStress[1] + TmpStress[2];
				for (i = 0; i < Size; i++)
					devStr[i] = TmpStress[i] - I1 / 3.0 * Kronecker[i];
				J2 = 0.5 * TensorMutiplication2(devStr, devStr, Dim);
				sqrtJ2 = sqrt(J2);
				J3 = TensorMutiplication3(devStr, devStr, devStr, Dim);

				if (J2 == 0)
				{
					LodeAngle = 0; // avoid error
					lode_out_range = true;
				}
				else
				{
					if ((-3 * sqrt3 / 2. * J3 / pow(sqrtJ2, 3)) > 1 - 1e-6)
					{
						LodeAngle = PI / 6.;
						lode_out_range = true;
					}
					else if ((-3 * sqrt3 / 2. * J3 / pow(sqrtJ2, 3)) < -1 + 1e-6)
					{
						LodeAngle = -PI / 6.;
						lode_out_range = true;
					}
					else
					{
						LodeAngle = asin(-3 * sqrt3 / 2. * J3 / pow(sqrtJ2, 3)) / 3.;
						lode_out_range = false;
					}
				}
				shearsurf = Ntheta * (I1 / 3 + 2 / sqrt3 * sqrtJ2 * sin(LodeAngle + 2 / 3. * PI))
				            - (I1 / 3. + 2 / sqrt3 * sqrtJ2 * sin(LodeAngle - 2 / 3. * PI)) - csn;
				tensionsurf_k1 = (I1 / 3. + 2 / sqrt3 * sqrtJ2 * sin(LodeAngle + 2 / 3. * PI)) - tension;
			} // end while for dlamda
			// calculate Dep
			// dsig_dlamda = -De dG/dstr
			TmpValue1 = 0.;
			TmpValue2 = 0.;
			for (i = 0; i < Size; i++)
			{
				TmpValue1 += dft_dsig[i] * (-1) * dsig_dlamda[i]; //(dF/dstr)T De dG/dstr
				TmpValue2 += dft_dsig[i] * dstrs_local[i];
			}
			*TmpMatrix = (0.);
			*TmpMatrix2 = (0.);
			// TEST
			// TmpMatrix->Write();
			// TmpMatrix2->Write();
			//
			for (i = 0; i < Size; i++)
				for (j = 0; j < Size; j++)
					(*TmpMatrix)(i, j) = dgt_dsig[i] * dft_dsig[j];
			Dep->multi(*TmpMatrix, *Dep, *TmpMatrix2); // De (dG/dstr)T dF/dstr De
			//
			// TmpDe->Write();
			// TmpMatrix2->Write();
			//
			for (i = 0; i < Size; i++)
			{
				for (j = 0; j < Size; j++)
					(*TmpDe)(i, j) -= (*TmpMatrix2)(i, j) / TmpValue1; // Dep = De - De dG/ds dF/ds De / dF/ds de dG/ds
			}
			// TEST
			// TmpDe->Write();
			//
		} // end if tension

		if (shearsurf > 1e-5)
		{
			dlamda = 0; //, dlamda_1=1, dlamda_0=0.;
			counter = 0;
			first_step = true;
			// local_damp = 1.;
			shearsurf_k1 = shearsurf;
			for (i = 0; i < Size; i++)
			{
				dstrs_local[i] = TmpStress[i] - TryStress_local[i];
				TryStr_buff[i] = TmpStress[i];
			}
			while (1)
			{
				for (i = 0; i < Size; i++)
					dsqrtJ2_dsig[i] = devStr[i] / 2. / sqrtJ2;

				if (lode_out_range || first_step)
				{
					for (i = 0; i < Size; i++)
						dlode_dsig[i] = 0;
				}
				else
				{
					dlode_dsig[0] = sqrt3 / (2. * MathLib::fastpow(sqrtJ2, 3) * cos(3. * LodeAngle))
					                * (3. * J3 / 2. / J2 * devStr[0]
					                   - (devStr[0] * devStr[0] + devStr[3] * devStr[3] + devStr[4] * devStr[4])
					                   + 2. * J2 * Kronecker[0] / 3.0);

					dlode_dsig[1] = sqrt3 / (2. * MathLib::fastpow(sqrtJ2, 3) * cos(3. * LodeAngle))
					                * (3. * J3 / 2. / J2 * devStr[1]
					                   - (devStr[3] * devStr[3] + devStr[1] * devStr[1] + devStr[5] * devStr[5])
					                   + 2. * J2 * Kronecker[1] / 3.0);

					dlode_dsig[2] = sqrt3 / (2. * MathLib::fastpow(sqrtJ2, 3) * cos(3. * LodeAngle))
					                * (3. * J3 / 2. / J2 * devStr[2]
					                   - (devStr[4] * devStr[4] + devStr[5] * devStr[5] + devStr[2] * devStr[2])
					                   + 2. * J2 * Kronecker[2] / 3.0);

					dlode_dsig[3] = sqrt3 / (2. * MathLib::fastpow(sqrtJ2, 3) * cos(3. * LodeAngle))
					                * (3. * J3 / 2. / J2 * devStr[3]
					                   - (devStr[0] * devStr[3] + devStr[3] * devStr[1] + devStr[4] * devStr[5])
					                   + 2. * J2 * Kronecker[3] / 3.0);

					if (Dim == 3)
					{
						dlode_dsig[4] = sqrt3 / (2. * MathLib::fastpow(sqrtJ2, 3) * cos(3. * LodeAngle))
						                * (3. * J3 / 2. / J2 * devStr[4]
						                   - (devStr[0] * devStr[4] + devStr[3] * devStr[5] + devStr[4] * devStr[2])
						                   + 2. * J2 * Kronecker[4] / 3.0);

						dlode_dsig[5] = sqrt3 / (2. * MathLib::fastpow(sqrtJ2, 3) * cos(3. * LodeAngle))
						                * (3. * J3 / 2. / J2 * devStr[5]
						                   - (devStr[3] * devStr[4] + devStr[1] * devStr[5] + devStr[5] * devStr[2])
						                   + 2. * J2 * Kronecker[5] / 3.0);
					}
				}

				if (Plasticity_Bedding)
				{
					double a_pk_sig_pq_sig_kq = 0., sig_pq_sig_pq = 0., dComp_daPara = 0.;
					double a_ki_sig_kj[6] = {0.};
					int p, k, q;
					*TmpMatrix = (0.);
					*TmpMatrix2 = (0.);
					for (i = 0; i < 3; i++)
					{
						for (j = 0; j < 3; j++)
						{
							if (i == j)
							{
								(*TmpMatrix)(i, j) = MicroStruTensor[i];
								//(*TmpMatrix2)(i,j) = TmpStress[i];
								(*TmpMatrix2)(i, j) = TmpStrComp[i];
							}
							else
								//(*TmpMatrix2)(i,j) = TmpStress[i+j+2];
								(*TmpMatrix2)(i, j) = TmpStrComp[i + j + 2];
						}
					}

					for (p = 0; p < 3; p++) // aki_sigkj, sigpq_sigpq = sigmn_sigmn, apk_sigpq_sigkq
					{
						a_ki_sig_kj[0] += (*TmpMatrix)(p, 0) * (*TmpMatrix2)(p, 0);
						a_ki_sig_kj[1] += (*TmpMatrix)(p, 1) * (*TmpMatrix2)(p, 1);
						a_ki_sig_kj[2] += (*TmpMatrix)(p, 2) * (*TmpMatrix2)(p, 2);
						a_ki_sig_kj[3] += (*TmpMatrix)(p, 0) * (*TmpMatrix2)(p, 1);
						a_ki_sig_kj[4] += (*TmpMatrix)(p, 0) * (*TmpMatrix2)(p, 2);
						a_ki_sig_kj[5] += (*TmpMatrix)(p, 1) * (*TmpMatrix2)(p, 2);
						for (q = 0; q < 3; q++)
						{
							sig_pq_sig_pq += (*TmpMatrix2)(p, q) * (*TmpMatrix2)(p, q);
							for (k = 0; k < 3; k++)
								a_pk_sig_pq_sig_kq += (*TmpMatrix)(p, k) * (*TmpMatrix2)(p, q) * (*TmpMatrix2)(k, q);
						}
					}
					for (i = 0; i < Size; i++)
						dAniso_dsig_comp[i] = 2 * (a_ki_sig_kj[i] * sig_pq_sig_pq - a_pk_sig_pq_sig_kq * TmpStrComp[i])
						                      / (sig_pq_sig_pq * sig_pq_sig_pq); // instead TmpStress with TmpStrComp

					for (i = 0; i < bedding_uc_curve_order + 1; i++)
						dComp_daPara += i * comp_para[i] * pow(AnisoParaComp, (i - 1)); // dcsn/deta
					for (i = 0; i < Size; i++)
						dcsn_dsig[i] = dComp_daPara * dAniso_dsig_comp[i]; // dcsn/dsig=dcsn/deta * deta/dsig
				} // end if plasticity bedding
				// counter++;
				if (first_step)
				{
					for (i = 0; i < Size; i++)
					{
						dfs_dsig[i] = (Ntheta - 1) / 3. * Kronecker[i]
						              + 2 / sqrt3
						                    * (Ntheta * sin(LodeAngle + 2 / 3. * PI) - sin(LodeAngle - 2 / 3. * PI))
						                    * dsqrtJ2_dsig[i]
						              + 2 / sqrt3 * sqrtJ2 * dlode_dsig[i]
						                    * (Ntheta * cos(LodeAngle + 2 / 3. * PI) - cos(LodeAngle - 2 / 3. * PI))
						              - dcsn_dsig[i];
						dgs_dsig[i] = (Nphi - 1) / 3. * Kronecker[i]
						              + 2 / sqrt3 * (Nphi * sin(LodeAngle + 2 / 3. * PI) - sin(LodeAngle - 2 / 3. * PI))
						                    * dsqrtJ2_dsig[i]
						              + 2 / sqrt3 * sqrtJ2 * dlode_dsig[i]
						                    * (Nphi * cos(LodeAngle + 2 / 3. * PI) - cos(LodeAngle - 2 / 3. * PI));
						//-dcsn_dsig[i];
					}
				}
				//
				if (counter == 0)
				{
					if (shearsurf_k1 > 0)
					{
						lamda_pos = dlamda;
						TmpValue1 = 0.;
						if (first_step)
						{
							for (i = 0; i < Size; i++)
								dsig_dlamda[i] = 0.;
							Dep->multi(dgs_dsig, dsig_dlamda, -1);
							lamda_pos = 0.;
							shearsurf_pos = shearsurf;
							// first_step = false;
						}
						for (i = 0; i < Size; i++)
						{
							TmpValue1 += dfs_dsig[i] * (-1) * dsig_dlamda[i]; //(dF/dstr)T De dG/dstr
							//	TmpValue2 += dfs_dsig[i]*dstrs[i];
						}
						if (TmpValue1 == 0)
							dlamda = 0.;
						else
							dlamda = shearsurf / TmpValue1;
						shearsurf_pos = shearsurf_k1;
						if (first_step)
						{
							dlamda = lamda_pos + dlamda;
							first_step = false;
						}
						else
							dlamda = lamda_pos + dlamda / 10.;
					}
					else
					{
						shearsurf_neg = shearsurf_k1;
						lamda_neg = dlamda;
						counter++;
					}
				}
				else
				{
					counter++;
					if (fabs(shearsurf_k1) < 1e-3)
					{
						// for(i=0;i<Size;i++)
						//	dsig_dlamda[i]=0.;
						// Dep->multi(dgs_dsig,dsig_dlamda,-1);//update De dG/dsig
						break;
					}
					else if (shearsurf_k1 < 0)
					{
						// for(i=0;i<Size;i++)
						//	dsig_dlamda[i]=0.;
						// Dep->multi(dgs_dsig,dsig_dlamda,-1);//update De dG/dsig
						shearsurf_neg = shearsurf_k1;
						lamda_neg = dlamda;
						dlamda = lamda_pos - shearsurf_pos * (lamda_neg - lamda_pos) / (shearsurf_neg - shearsurf_pos);
						if (fabs(lamda_neg - lamda_pos) < 1e-12)
							break;
						// dlamda = (lamda_pos+lamda_neg)/2.;
					}
					else
					{
						// for(i=0;i<Size;i++)
						//	dsig_dlamda[i]=0.;
						// Dep->multi(dgs_dsig,dsig_dlamda,-1);//update De dG/dsig
						lamda_pos = dlamda;
						shearsurf_pos = shearsurf_k1;
						dlamda = lamda_pos - shearsurf_pos * (lamda_neg - lamda_pos) / (shearsurf_neg - shearsurf_pos);
						if (fabs(lamda_neg - lamda_pos) < 1e-12)
							break;
						// dlamda = (lamda_pos+lamda_neg)/2.;
					}
				}

				for (i = 0; i < Size; i++) // local_damp for Newton downhill
					TmpStress[i] = TryStr_buff[i] + dlamda * dsig_dlamda[i]; // sig(n+1)=sig(try)-dlambda*De*dG/dsig

				if (!Plasticity_Bedding)
				{
					CalculateCoefficent_MOHR(ep, 0, 0);
				}
				else
				{
					*TmpMatrix = (0.);
					*TmpMatrix2 = (0.);
					CalPrinStrDir(TmpStress, tmp_prin_str, tmp_prin_dir, Dim);
					CalTransMatrixA(tmp_prin_dir, TmpMatrix, Size);
					TmpMatrix->GetTranspose(*TmpMatrix2);

					for (i = 0; i < Size; i++)
					{
						TmpStrComp[i] = 0.;
						TmpStrTens[i] = 0.;
						if (tmp_prin_str[i] < MKleinsteZahl)
							TmpPrinStrComp[i] = tmp_prin_str[i];
						else
							TmpPrinStrTens[i] = tmp_prin_str[i];
					}
					TmpMatrix2->multi(TmpPrinStrComp, TmpStrComp);
					TmpMatrix2->multi(TmpPrinStrTens, TmpStrTens);
					AnisoParaComp = CalAnisoPara(TmpStrComp, MicroStruTensor);
					AnisoParaTens = CalAnisoPara(TmpStrTens, MicroStruTensor);
					CalculateCoefficent_MOHR(ep, AnisoParaComp, AnisoParaTens);
				}

				I1 = TmpStress[0] + TmpStress[1] + TmpStress[2];
				for (i = 0; i < Size; i++)
					devStr[i] = TmpStress[i] - I1 / 3.0 * Kronecker[i];
				J2 = 0.5 * TensorMutiplication2(devStr, devStr, Dim);
				sqrtJ2 = sqrt(J2);
				J3 = TensorMutiplication3(devStr, devStr, devStr, Dim);

				if (J2 == 0)
				{
					LodeAngle = 0; // avoid error
					lode_out_range = true;
				}
				else
				{
					if ((-3 * sqrt3 / 2. * J3 / pow(sqrtJ2, 3)) > 1 - 1e-6)
					{
						LodeAngle = PI / 6.;
						lode_out_range = true;
					}
					else if ((-3 * sqrt3 / 2. * J3 / pow(sqrtJ2, 3)) < -1 + 1e-6)
					{
						LodeAngle = -PI / 6.;
						lode_out_range = true;
					}
					else
					{
						LodeAngle = asin(-3 * sqrt3 / 2. * J3 / pow(sqrtJ2, 3)) / 3.;
						lode_out_range = false;
					}
				}
				// shearsurf=Ntheta*(I1/3+2/sqrt3*sqrtJ2*sin(LodeAngle+2/3.*PI))
				//	-(I1/3.+2/sqrt3*sqrtJ2*sin(LodeAngle-2/3.*PI))-csn;
				shearsurf_k1 = Ntheta * (I1 / 3 + 2 / sqrt3 * sqrtJ2 * sin(LodeAngle + 2 / 3. * PI))
				               - (I1 / 3. + 2 / sqrt3 * sqrtJ2 * sin(LodeAngle - 2 / 3. * PI)) - csn;
				tensionsurf = (I1 / 3. + 2 / sqrt3 * sqrtJ2 * sin(LodeAngle + 2 / 3. * PI)) - tension;
			} // end while for dlamda
			// calculate Dep
			// dsig_dlamda = -De dG/dstr
			TmpValue1 = 0.;
			TmpValue2 = 0.;
			for (i = 0; i < Size; i++)
			{
				TmpValue1 += dfs_dsig[i] * (-1) * dsig_dlamda[i]; //(dF/dstr)T De dG/dstr
				TmpValue2 += dfs_dsig[i] * dstrs_local[i]; //(df/dsig)T dstress
			}
			*TmpMatrix = (0.);
			*TmpMatrix2 = (0.);
			// TEST
			// TmpMatrix->Write();
			// TmpMatrix2->Write();
			//
			for (i = 0; i < Size; i++)
				for (j = 0; j < Size; j++)
					(*TmpMatrix)(i, j) = dgs_dsig[i] * dfs_dsig[j];
			Dep->multi(*TmpMatrix, *Dep, *TmpMatrix2); // De (dG/dstr)T dF/dstr De
			//
			// TmpDe->Write();
			// TmpMatrix2->Write();
			//
			for (i = 0; i < Size; i++)
			{
				for (j = 0; j < Size; j++)
					(*TmpDe)(i, j) -= (*TmpMatrix2)(i, j) / TmpValue1; // Dep = De - De dG/ds dF/ds De / dF/ds de dG/ds
			}
			// TEST
			// TransMicroStru->Write();
			// TransMicroStru_T->Write();
			// TransMicroStru_TInv->Write();
			// TmpDe->Write();
			// Dep->Write();
			//
		} // end if shear

		for (i = 0; i < Size; i++)
			TryStress[i] = 0;
		TransMicroStru_T->multi(TmpStress, TryStress);

		*ConstitutiveMatrix = (0.);
		TransMicroStru_T->multi(*TmpDe, *TransMicroStru, *ConstitutiveMatrix);
		// TEST
		// ConstitutiveMatrix->Write();
		//
	} // end if shear if tension
	else
	{
		for (i = 0; i < Size; i++)
			TryStress[i] += dstrs[i];
	}
	if (Update)
	{
		Cal_Inv_Matrix(Size, Dep, Inv_De);
		for (i = 0; i < Size; i++)
			dstressP[i] = TryStress[i] - TryStress_0[i];
		Inv_De->multi(dstressP, dstrainP);
		ep += sqrt(2.0 / 3.0 * (TensorMutiplication2(dstrainP, dstrainP, Dim)));
		(*ele_val->pStrain)(GPiGPj) = ep;
	}
	return yield;
}
/*******************************************************
WX: calculate aniso. plas. parameters
*******************************************************/
double CSolidProperties::CalAnisoPara(double* Stress, double* MicroStruTensor)
{
	double L1, L2, L3, LkLk, AnisoPara;
	double l[3] = {0.};
	L1 = sqrt(Stress[0] * Stress[0] + Stress[3] * Stress[3] + Stress[4] * Stress[4]);
	L2 = sqrt(Stress[3] * Stress[3] + Stress[1] * Stress[1] + Stress[5] * Stress[5]);
	L3 = sqrt(Stress[4] * Stress[4] + Stress[5] * Stress[5] + Stress[2] * Stress[2]);

	LkLk = sqrt(L1 * L1 + L2 * L2 + L3 * L3);
	if (LkLk < MKleinsteZahl)
		return 1;

	l[0] = L1 / LkLk;
	l[1] = L2 / LkLk;
	l[2] = L3 / LkLk;

	AnisoPara = 0;
	for (int i = 0; i < 3; i++)
		AnisoPara += MicroStruTensor[i] * l[i] * l[i];

	return AnisoPara;
}
//*/
/*******************************************************
WX: Mohr coulomb,
return mapping (not direct stress integ.) also for aniso.
*******************************************************/
int CSolidProperties::DirectStressIntegrationMOHR(const int GPiGPj, ElementValue_DM* ele_val, double* TryStress,
                                                  const int Update, Matrix* Dep, int itesteps)
{
	int i, j;
	int yield = 0;
	int Dim = 2;
	int Size = ele_val->Stress->Rows();
	double TmpValue1, TmpValue2, dstrNorm = 0;
	// double LodeAngle, I1, J2, J3;
	double shearsurf, tensionsurf, ep;

	// initialize all vectors
	double dstrs[6] = {0.}, TmpStress[6] = {0.}, prin_str[6] = {0.}, prin_str0[6] = {0.}, prin_dir[9] = {0.};

	*TransMatrixA = (0.);
	*TransMatrixA_T = (0.);
	*Inv_TransA = (0.);
	*Inv_TransA_T = (0.);
	*TmpDe = (0.);
	*Inv_De = (0.);
	*TMatrix = (0.);
	*TmpMatrix = (0.);
	*TmpMatrix2 = (0.);

	*TmpDe = *Dep;

	*ConstitutiveMatrix = (0.); // in head already defined, and is used for later as global variable

	ep = (*ele_val->pStrain)(GPiGPj); // get equ plas strain

	if (Size > 4)
		Dim = 3;

	for (i = 0; i < Size; i++)
	{
		dstrs[i] = TryStress[i]; // d_stress
		TryStress[i] = (*ele_val->Stress)(i, GPiGPj); // stress_0
		devS[i] = TryStress[i] + dstrs[i];
		TmpStress[i] = devS[i];
		dstrNorm += dstrs[i] * dstrs[i];
	}
	if (Size == 4)
		devS[4] = devS[5] = 0.;

	////////////
	/*/test
	   devS[0]=TmpStress[0]=-1.5356885185724423e-006;
	   devS[1]=TmpStress[1]=5250.0000001923481;
	   devS[2]=TmpStress[2]=2.1768577340708362e-6;
	   devS[3]=TmpStress[3]=-5.7126080120270688e-007;
	   devS[4]=TmpStress[4]=-2.8646737121460246e-007;
	   devS[5]=TmpStress[5]=7.9612875121132550e-007;

	 */ ///////////
	// test
	CalPrinStrDir(devS, prin_str, prin_dir, Dim);

	// CalPrinStrs(devS, prin_str, Size);  //prin. stresses guess
	// CalPrinDir(prin_str, TmpStress, prin_dir, Size);
	CalTransMatrixA(prin_dir, TransMatrixA, Size);
	TransMatrixA->GetTranspose(*TransMatrixA_T);
	Cal_Inv_Matrix(Size, TransMatrixA, Inv_TransA);
	Inv_TransA->GetTranspose(*Inv_TransA_T);
	// Inv_TransA_T->multi(*Dep, *Inv_TransA, *PrinDe);  //De in prin. (A^-T De A^-1)

	*TmpDe = (0.);
	Inv_TransA_T->multi(*Dep, *Inv_TransA, *TmpDe); // De in prin. coord.

	double TmpStress0[6] = {0.}, tmp_prin_str[6] = {0.}, tmp_prin_dir[9] = {0.};
	for (i = 0; i < Size; i++)
	{
		TmpStress0[i] = TryStress[i];
	}
	if (Size == 4)
		TmpStress0[4] = TmpStress[5] = 0.;

	// CalPrinStrs(TmpStress0, tmp_prin_str, Size);
	CalPrinStrDir(TmpStress0, tmp_prin_str, tmp_prin_dir, Dim); // prin. stresses t0

	/*/////////
	   double tmpresult[6]={0.};
	   TransMatrixA_T->multi(prin_str, tmpresult);
	   for(i=0; i<Size; i++)
		cout<<tmpresult[i]<<endl;
	 */ /////////

	// if(m_pcs->GetIteSteps()==1)

	if (Plasticity_Bedding) // WX:09.2011
	{
		if (itesteps == 2 || Update > 0)
		{
			/*double tmp_load_dir_major = 0;
			double tmp_load_dir_minor = 0;
			for(i=0; i<3; i++)
			{
			    tmp_load_dir_major += Bedding_Norm[i]*prin_dir[3*i+2];
			    tmp_load_dir_minor += Bedding_Norm[i]*prin_dir[3*i];
			}*/

			double Stress_Bedding_comp[6] = {0.}, Stress_Bedding_tens[6] = {0.}; // Stress_Bedding[6]={0.}
			double prin_str_comp[6] = {0.}, prin_str_tens[6] = {0.};
			double TmpStress_comp[6] = {0.}, TmpStress_tens[6] = {0.};
			double L_sqr_comp[3] = {0.}, L_sqr_tens[3] = {0.};
			double tr_sig_comp = 0., tr_sig_tens = 0., tr_a_sig_comp = 0., tr_a_sig_tens = 0.;
			double scalar_comp = 0, scalar_tens = 0;

			for (i = 0; i < 3; i++)
			{
				if (prin_str[i] <= 0)
				{
					prin_str_comp[i] = prin_str[i];
					prin_str_tens[i] = 0;
				}
				else
				{
					prin_str_comp[i] = 0;
					prin_str_tens[i] = prin_str[i];
				}
			}
			TransMatrixA_T->multi(prin_str_comp, TmpStress_comp);
			TransMatrixA_T->multi(prin_str_tens, TmpStress_tens);
			TransMicroStru_TInv->multi(TmpStress_comp, Stress_Bedding_comp);
			TransMicroStru_TInv->multi(TmpStress_tens, Stress_Bedding_tens);
			/*for(i=0;i<6;i++)
			{
			    if(Stress_Bedding[i]<=0)
			    {
			        Stress_Bedding_comp[i]=Stress_Bedding[i];
			        Stress_Bedding_tens[i]=0;
			    }
			    else
			    {
			        Stress_Bedding_comp[i]=0;
			        Stress_Bedding_tens[i]=Stress_Bedding[i];
			    }
			}*/
			// for(i=0;i<6;i++)
			//{
			//	 tr_sig_comp += Stress_Bedding_comp[i]*Stress_Bedding_comp[i];
			//	 tr_sig_tens += Stress_Bedding_tens[i]*Stress_Bedding_tens[i];
			//}
			L_sqr_comp[0] = Stress_Bedding_comp[0] * Stress_Bedding_comp[0]
			                + Stress_Bedding_comp[3] * Stress_Bedding_comp[3]
			                + Stress_Bedding_comp[4] * Stress_Bedding_comp[4];
			L_sqr_comp[1] = Stress_Bedding_comp[1] * Stress_Bedding_comp[1]
			                + Stress_Bedding_comp[3] * Stress_Bedding_comp[3]
			                + Stress_Bedding_comp[5] * Stress_Bedding_comp[5];
			L_sqr_comp[2] = Stress_Bedding_comp[2] * Stress_Bedding_comp[2]
			                + Stress_Bedding_comp[5] * Stress_Bedding_comp[5]
			                + Stress_Bedding_comp[4] * Stress_Bedding_comp[4];
			tr_sig_comp = L_sqr_comp[0] + L_sqr_comp[1] + L_sqr_comp[2];

			L_sqr_tens[0] = Stress_Bedding_tens[0] * Stress_Bedding_tens[0]
			                + Stress_Bedding_tens[3] * Stress_Bedding_tens[3]
			                + Stress_Bedding_tens[4] * Stress_Bedding_tens[4];
			L_sqr_tens[1] = Stress_Bedding_tens[1] * Stress_Bedding_tens[1]
			                + Stress_Bedding_tens[3] * Stress_Bedding_tens[3]
			                + Stress_Bedding_tens[5] * Stress_Bedding_tens[5];
			L_sqr_tens[2] = Stress_Bedding_tens[2] * Stress_Bedding_tens[2]
			                + Stress_Bedding_tens[5] * Stress_Bedding_tens[5]
			                + Stress_Bedding_tens[4] * Stress_Bedding_tens[4];
			tr_sig_tens = L_sqr_tens[0] + L_sqr_tens[1] + L_sqr_tens[2];
			for (i = 0; i < 3; i++)
			{
				tr_a_sig_comp += MicroStruTensor[i] * L_sqr_comp[i];
				tr_a_sig_tens += MicroStruTensor[i] * L_sqr_tens[i];
			}

			//
			// tr_a_sig_comp += MicroStruTensor[i]*(Stress_Bedding_comp[i]*Stress_Bedding_comp[i]
			//+);
			// tr_a_sig_tens += MicroStruTensor[i]*Stress_Bedding_tens[i]*Stress_Bedding_tens[i];
			if (tr_sig_comp == 0)
				scalar_comp = 0;
			else
				scalar_comp = tr_a_sig_comp / tr_sig_comp;
			if (tr_sig_tens == 0)
				scalar_tens = 0;
			else
				scalar_tens = tr_a_sig_tens / tr_sig_tens;
			if (scalar_tens < MKleinsteZahl)
				scalar_tens = 1;
			// if((fabs(tmp_load_dir_major)>(1+1e-6))||(fabs(tmp_load_dir_minor)>(1+1e-6)))
			//{
			//	cout<<"ERROR: Please check input value for BEDDING_NORM in .msp"<<endl;
			//	exit(1);
			//}
			// if(fabs(tmp_load_dir_major)>1)
			//	tmp_load_dir_major = 1;
			// if (fabs(tmp_load_dir_minor)>1)
			//	tmp_load_dir_minor = 1;
			// tmp_load_dir_major = acos(fabs(tmp_load_dir_major))*180/PI; //because ||bedding_norm|| and ||prin_dir||
			// is 1
			// tmp_load_dir_minor = acos(fabs(tmp_load_dir_minor))*180/PI;
			(*ele_val->scalar_aniso_comp)(GPiGPj) = scalar_comp;
			(*ele_val->scalar_aniso_tens)(GPiGPj) = scalar_tens;
		}

		CalculateCoefficent_MOHR(ep, (*ele_val->scalar_aniso_comp)(GPiGPj), (*ele_val->scalar_aniso_tens)(GPiGPj));
	}
	else
		CalculateCoefficent_MOHR(ep, 0, 0);

	shearsurf = Ntheta * prin_str[0] - prin_str[2] - csn;
	tensionsurf = prin_str[0] - tension;
	if (dstrNorm == 0)
	{
		shearsurf = -1;
		tensionsurf = -1;
	}
	if ((shearsurf > 0) || (tensionsurf > 0))
	{
		double P_12, P_31, P_41, P_63, P_64, P_52, P_85, P_74, P_78, P_98, P_45 /*, P_X7*/;
		double t1, t2 /*, t1ra*/, t1r1 /*, t2ra, t2r2*/, t3r1 /*, t3r2*/;
		double fkt1, fkt2;
		double mm = 0.;

		double l1[6] = {0.}, l2[6] = {0.}, l1g[6] = {0.}, l2g[6] = {0.}, l1R[6] = {0.}, l2R[6] = {0.}, l3R[6] = {0.};
		double rsp[6] = {0.}, rtp[6] = {0.};
		double sigA[6] = {0.}, sig1R[6] = {0.}, sig2R[6] = {0.}, sigaR[6] = {0.};
		double dFsdprin_s[6] = {0.}, dFtdprin_s[6] = {0.}, dGsdprin_s[6] = {0.}, dGtdprin_s[6] = {0.};
		double De_dGsdprin_s[6] = {0.}, De_dGtdprin_s[6] = {0.};
		double dStressP[6] = {0.}, dStrainP[6] = {0.};
		double cVec[6] = {0};

		yield = 1;
		fkt1 = 0.001;
		fkt2 = 0.01;
		sigA[0] = sigA[1] = sigA[2] = csn / (Ntheta - 1);
		sig1R[0] = sig1R[1] = sig2R[0] = tension;
		sig1R[2] = sig2R[1] = sig2R[2] = Ntheta * tension - csn;
		sigaR[0] = sigaR[1] = sigaR[2] = tension;

		l1[0] = l1[1] = l1g[0] = l1g[1] = l2[0] = l2g[0] = 1.;
		l1[2] = l2[1] = l2[2] = Ntheta;
		l1g[2] = l2g[1] = l2g[2] = Nphi;
		l1R[2] = l2R[1] = l2R[2] = l3R[1] = 1.;

		dFsdprin_s[0] = Ntheta;
		dGsdprin_s[0] = Nphi;
		dFsdprin_s[2] = dGsdprin_s[2] = -1.;
		TmpDe->multi(dGsdprin_s, De_dGsdprin_s);
		// PrinDe->multi(dGsdprin_s, De_dGsdprin_s);//
		TmpValue1 = 0.;
		for (i = 0; i < Size; i++)
			TmpValue1 += dFsdprin_s[i] * De_dGsdprin_s[i];
		for (i = 0; i < Size; i++)
			rsp[i] = De_dGsdprin_s[i] / TmpValue1;

		dFtdprin_s[0] = dGtdprin_s[0] = 1.;
		TmpDe->multi(dGtdprin_s, De_dGtdprin_s);
		// PrinDe->multi(dGtdprin_s, De_dGtdprin_s);//
		TmpValue2 = 0.;
		for (i = 0; i < Size; i++)
			TmpValue2 += dFtdprin_s[i] * De_dGtdprin_s[i];
		for (i = 0; i < Size; i++)
			rtp[i] = De_dGtdprin_s[i] / TmpValue2;

		double tmp_shearsurf = Ntheta * tmp_prin_str[0] - tmp_prin_str[2] - csn;
		double tmp_tensionsurf = tmp_prin_str[0] - tension;
		if (((tmp_tensionsurf) == 0 && (tmp_shearsurf) <= 0) || ((tmp_tensionsurf) <= 0 && (tmp_shearsurf) == 0))
			mm = 0.;
		else if (prin_str[0] != tmp_prin_str[0])
		{
			mm = (tension - tmp_prin_str[0]) / (prin_str[0] - tmp_prin_str[0]);
			if (mm >= 0 && mm <= 1)
			{
				double tmp_prin_str_3 = tmp_prin_str[2] + mm * (prin_str[2] - tmp_prin_str[2]);
				if (tmp_prin_str_3 < (Ntheta * tension - csn) || tmp_prin_str_3 > tension)
					mm = (csn + tmp_prin_str[2] - Ntheta * tmp_prin_str[0])
					     / (Ntheta * (prin_str[0] - tmp_prin_str[0]) - (prin_str[2] - tmp_prin_str[2]));
			}
			else
				mm = (csn + tmp_prin_str[2] - Ntheta * tmp_prin_str[0])
				     / (Ntheta * (prin_str[0] - tmp_prin_str[0]) - (prin_str[2] - tmp_prin_str[2]));
		}
		else
			mm = (csn + tmp_prin_str[2] - Ntheta * tmp_prin_str[0])
			     / (Ntheta * (prin_str[0] - tmp_prin_str[0]) - (prin_str[2] - tmp_prin_str[2]));

		Cal_Inv_Matrix(Size, TmpDe, Inv_De);
		// Cal_Inv_Matrix(Size, PrinDe, Inv_De);

		t1 = CalVar_t(l1, l1g, Inv_De, prin_str, sig1R, Size);
		t2 = CalVar_t(l2, l2g, Inv_De, prin_str, sig2R, Size);
		t1r1 = CalVar_t(l1R, l1R, Inv_De, prin_str, sig1R, Size);
		// t1ra = CalVar_t(l1R, l1R, Inv_De, prin_str, sigaR, Size);
		// t2r2 = CalVar_t(l2R, l2R, Inv_De, prin_str, sig2R, Size);
		// t2ra = CalVar_t(l2R, l2R, Inv_De, prin_str, sigaR, Size);
		t3r1 = CalVar_t(l3R, l3R, Inv_De, prin_str, sig1R, Size);
		// t3r2 = CalVar_t(l3R, l3R, Inv_De, prin_str, sig2R, Size);

		P_12 = CalVarP(rsp, l1, prin_str, sig1R);
		P_31 = CalVarP(rsp, l2, prin_str, sig2R);
		P_41 = CalVarP(rsp, l3R, prin_str, sig1R);
		// P_63 = 0;
		// P_63 = CalVarP(rsp, l2R, prin_str, sig2R);
		P_63 = P_41;
		P_45 = CalVarP(rsp, rtp, prin_str, sig1R);
		// P_64 = CalVarP(rsp, rtp, prin_str, sig2R);
		P_64 = CalVarP(rsp, rtp, prin_str, sig2R);
		// P_52 = 0;
		// P_52 = CalVarP(rsp, l1R, prin_str, sig1R);
		P_52 = P_41;
		P_74 = CalVarP(rtp, l3R, prin_str, sig1R);
		// P_85 = 0;
		// P_85 = CalVarP(rtp, l1, prin_str, sig1R);
		P_85 = P_74;
		P_78 = CalVarP(rtp, l1R, prin_str, sig1R);
		// P_98 = 0;
		P_98 = CalVarP(rtp, l2R, prin_str, sig2R);
		// P_X7 = CalVarP(rtp, l2R, prin_str, sig2R);

		if (P_12 >= 0 && P_31 <= 0 && P_41 <= 0) // return to fmc
		{
			double tmpvalue = 0.;
			// Matrix *tmpMatrix2 = new Matrix (Size,Size);
			// Matrix *dGds_dFds = new Matrix (Size,Size);
			*TmpMatrix2 = (0.); // WX:08.2011
			*dGds_dFds = (0.); // WX:08.2011
			for (i = 0; i < 3; i++)
			{
				tmpvalue += dFsdprin_s[i] * (prin_str[i] - sigA[i]);
			}
			for (i = 0; i < 3; i++)
			{
				prin_str0[i] = prin_str[i];
				prin_str[i] -= tmpvalue * rsp[i];
				dStressP[i] = prin_str[i] - prin_str0[i];
			}
			Inv_De->multi(dStressP, dStrainP);
			for (i = 0; i < 3; i++)
				for (j = 0; j < 3; j++)
					(*dGds_dFds)(i, j) = rsp[i] * dFsdprin_s[j]; //(D*dGds/(dFds*D*dGds))*dFds
			dGds_dFds->multi(*TmpDe, *TmpMatrix2);
			// dGds_dFds->multi(*PrinDe, *tmpMatrix2);
			*TmpDe -= *TmpMatrix2;
			//*PrinDe -= *tmpMatrix2;
			// for(i=3; i<Size; i++)
			//	(*TmpDe)(i,i) = (*Dep)(i,i);
			//*TmpDe = *PrinDe;

			// delete tmpMatrix2;
			// delete dGds_dFds;
		}
		// else if( P_12<0 && t1<0 )		//return to l1
		else if (P_12 < 0 && P_52 < 0) // return to l1
		{
			for (i = 0; i < 3; i++)
			{
				prin_str0[i] = prin_str[i];
				prin_str[i] = t1 * l1[i] + sig1R[i];
				dStressP[i] = prin_str[i] - prin_str0[i];
			}
			Inv_De->multi(dStressP, dStrainP); // dstrainP = D-1 * dstressp
			VecCrossProduct(dStrainP, l1g, cVec); // cVec = dstrainP X l1g

			CalDep_l(l1, l1g, Inv_De, Dep_l, 1.0); // Dep_l = l*lgT/(lT*D-1*lg)
			CalDep_l(cVec, cVec, Inv_De, dDep_l, fkt2); // dDep_l = fkt2*cVec*cVecT/(cVecT*D-1*cVec)
			//*TmpDe=(0.);
			for (i = 0; i < 3; i++)
				for (j = 0; j < 3; j++)
					(*TmpDe)(i, j) = (*Dep_l)(i, j) + (*dDep_l)(i, j);
			// for(i=3; i<Size; i++)
			//	(*TmpDe)(i,i) = (*Dep)(i,i);
		}
		// else if( P_31>0 && t2<0 )		//return to l2
		else if (P_31 > 0 && P_63 < 0) // return to l2
		{
			for (i = 0; i < 3; i++)
			{
				prin_str0[i] = prin_str[i];
				prin_str[i] = t2 * l2[i] + sig2R[i];
				dStressP[i] = prin_str[i] - prin_str0[i];
			}
			Inv_De->multi(dStressP, dStrainP);
			VecCrossProduct(dStrainP, l2g, cVec);
			CalDep_l(l2, l2g, Inv_De, Dep_l, 1.0);
			CalDep_l(cVec, cVec, Inv_De, dDep_l, fkt2);
			//*TmpDe=(0.);
			for (i = 0; i < 3; i++)
				for (j = 0; j < 3; j++)
					(*TmpDe)(i, j) = (*Dep_l)(i, j) + (*dDep_l)(i, j);
			// for(i=3; i<Size; i++)
			//	(*TmpDe)(i,i) = (*Dep)(i,i);
		}
		// else if( P_41>0 && P_74<0 && t3r2>0 && t3r1<0 )		//return to l3R
		else if (P_41 > 0 && P_45 > 0 && P_64 < 0 && P_74 < 0) // return to l3R
		{
			for (i = 0; i < 3; i++)
			{
				prin_str0[i] = prin_str[i];
				prin_str[i] = t3r1 * l3R[i] + sig1R[i];
				dStressP[i] = prin_str[i] - prin_str0[i];
			}
			Inv_De->multi(dStressP, dStrainP);
			VecCrossProduct(dStrainP, l3R, cVec);
			CalDep_l(l3R, l3R, Inv_De, Dep_l, 1.0);
			CalDep_l(cVec, cVec, Inv_De, dDep_l, fkt2);
			//*TmpDe=(0.);
			for (i = 0; i < 3; i++)
				for (j = 0; j < 3; j++)
					(*TmpDe)(i, j) = (*Dep_l)(i, j) + (*dDep_l)(i, j);
			// for(i=3; i<Size; i++)
			//	(*TmpDe)(i,i) = (*Dep)(i,i);
		}
		else if (P_74 >= 0 && P_78 >= 0) // return to ft, && P_X7<=0
		{
			double tmpvalue = 0.;
			// Matrix *tmpMatrix2 = new Matrix (Size,Size);
			// Matrix *dGds_dFds = new Matrix (Size,Size);
			*TmpMatrix2 = (0.); // WX:08.2011
			*dGds_dFds = (0.); // WX:08.2011
			for (i = 0; i < 3; i++)
				tmpvalue += dFtdprin_s[i] * (prin_str[i] - sigaR[i]);
			for (i = 0; i < 3; i++)
			{
				prin_str0[i] = prin_str[i];
				prin_str[i] -= tmpvalue * rtp[i];
				dStressP[i] = prin_str[i] - prin_str0[i];
			}
			Inv_De->multi(dStressP, dStrainP);

			for (i = 0; i < 3; i++)
				for (j = 0; j < 3; j++)
					(*dGds_dFds)(i, j) = rtp[i] * dFtdprin_s[j];
			dGds_dFds->multi(*TmpDe, *TmpMatrix2);
			// dGds_dFds->multi(*PrinDe, *tmpMatrix2);
			*TmpDe -= *TmpMatrix2;
			//*PrinDe -= *tmpMatrix2;
			for (i = 0; i < 3; i++)
				for (j = 0; j < 3; j++)
					if (abs((*TmpDe)(i, j)) < MKleinsteZahl)
						(*TmpDe)(i, j) = 1.; //(*Dep)(i,j)*1e-9;//avoid "0" in Dep
			// for(i=3; i<Size; i++)
			//	(*TmpDe)(i,i) = (*Dep)(i,i);
			//*TmpDe = *PrinDe;
			// TmpDe->Write();
			//*TmpDe = *Dep;				//to be improved
			// delete tmpMatrix2;
			// delete dGds_dFds;
		}
		// else if(P_X7>0 && t2r2>0 && t2ra<0)			//return to l2R
		//{
		//	for(i=0; i<3; i++)
		//	{
		//		prin_str0[i] = prin_str[i];
		//		prin_str[i] = t2r2*l2R[i] + sig2R[i];
		//		dStressP[i] = prin_str[i] - prin_str0[i];
		//	}
		//	Inv_De->multi(dStressP,dStrainP);
		//	VecCrossProduct(dStrainP,l2R,cVec);
		//	CalDep_l(l2R, l2R, Inv_De, Dep_l, fkt2);
		//	CalDep_l(cVec, cVec, Inv_De, dDep_l, fkt2);
		//	*TmpDe=(0.);
		//	for(i=0; i<3; i++)
		//		for(j=0; j<3; j++)
		//			(*TmpDe)(i,j) = (*Dep_l)(i,j)+(*dDep_l)(i,j);
		//	for(i=3; i<Size; i++)
		//		(*TmpDe)(i,i) = G;
		//}
		// else if( P_78<0 && t1r1>0 && t1ra<0 )		//return to l1R
		// else if( P_52>=0 && P_45<=0 && P_85<=0 )		//return to l1R
		else if (P_85 >= 0 && P_78 <= 0 && P_98 <= 0) // return to l1R
		{
			for (i = 0; i < 3; i++)
			{
				prin_str0[i] = prin_str[i];
				prin_str[i] = t1r1 * l1R[i] + sig1R[i];
				dStressP[i] = prin_str[i] - prin_str0[i];
			}
			Inv_De->multi(dStressP, dStrainP);
			VecCrossProduct(dStrainP, l1R, cVec);
			CalDep_l(l1R, l1R, Inv_De, Dep_l, 1.0);
			CalDep_l(cVec, cVec, Inv_De, dDep_l, fkt2);
			//*TmpDe=(0.);
			for (i = 0; i < 3; i++)
				for (j = 0; j < 3; j++)
					(*TmpDe)(i, j) = (*Dep_l)(i, j) + (*dDep_l)(i, j);
			// for(i=3; i<Size; i++)
			//	(*TmpDe)(i,i) = (*Dep)(i,i);
		}
		// else if( t2ra>=0 && t1ra>=0 )		//return to sigRa
		else if (P_98 >= 0) // return to sigRa
		{
			double tmpvalue = 0.;
			// Matrix *dGds_dFds = new Matrix (Size,Size);
			*dGds_dFds = (0.); // WX:08.2011

			for (i = 0; i < 3; i++)
			{
				prin_str0[i] = prin_str[i];
				prin_str[i] = sigaR[i];
				dStressP[i] = prin_str[i] - prin_str0[i];
			}
			Inv_De->multi(dStressP, dStrainP);
			for (i = 0; i < 3; i++)
				for (j = 0; j < 3; j++)
					(*dGds_dFds)(i, j) = dStressP[i] * dStressP[j];
			tmpvalue = 0.;
			for (i = 0; i < 3; i++)
				tmpvalue += dStrainP[i] * dStressP[i];
			*dGds_dFds /= tmpvalue;
			*TmpDe -= *dGds_dFds;
			//*PrinDe -= *dGds_dFds;
			//*TmpDe *= fkt1;
			for (i = 0; i < Size; i++)
				for (j = 0; j < Size; j++)
					(*TmpDe)(i, j) *= fkt1;
			//*PrinDe *= fkt1;

			// for(i=3; i<Size; i++)
			//	(*TmpDe)(i,i) = (*Dep)(i,i);
			//*TmpDe = *PrinDe;

			// delete dGds_dFds;
		}
		// else if(t3r1>=0)				//return to sig1R
		else if (P_52 >= 0 && P_45 <= 0 && P_85 <= 0) // return to sig1R
		{
			double tmpvalue = 0.;
			// Matrix *dGds_dFds = new Matrix (Size,Size);
			*dGds_dFds = (0.); // WX:08.2011

			for (i = 0; i < 3; i++)
			{
				prin_str0[i] = prin_str[i];
				prin_str[i] = sig1R[i];
				dStressP[i] = prin_str[i] - prin_str0[i];
			}
			Inv_De->multi(dStressP, dStrainP);
			for (i = 0; i < 3; i++)
				for (j = 0; j < 3; j++)
					(*dGds_dFds)(i, j) = dStressP[i] * dStressP[j];
			tmpvalue = 0.;
			for (i = 0; i < 3; i++)
				tmpvalue += dStrainP[i] * dStressP[i];
			*dGds_dFds /= tmpvalue;
			*TmpDe -= *dGds_dFds;
			//*PrinDe -= *dGds_dFds;
			//*TmpDe *= fkt1;
			for (i = 0; i < Size; i++)
				for (j = 0; j < Size; j++)
					(*TmpDe)(i, j) *= fkt1;
			//*PrinDe *= fkt1;
			//*TmpDe = *PrinDe;
			// for(i=3; i<Size; i++)
			//	(*TmpDe)(i,i) = (*Dep)(i,i);
			// delete dGds_dFds;
		}
		// else if(t3r2<=0)		//return to sig2R
		else if (P_63 >= 0 && P_64 >= 0) // return to sig2R
		{
			double tmpvalue = 0.;
			// Matrix *dGds_dFds = new Matrix (Size,Size);
			*dGds_dFds = (0.); // WX:08.2011

			for (i = 0; i < 3; i++)
			{
				prin_str0[i] = prin_str[i];
				prin_str[i] = sig2R[i];
				dStressP[i] = prin_str[i] - prin_str0[i];
			}
			Inv_De->multi(dStressP, dStrainP);
			for (i = 0; i < 3; i++)
				for (j = 0; j < 3; j++)
					(*dGds_dFds)(i, j) = dStressP[i] * dStressP[j];
			tmpvalue = 0.;
			for (i = 0; i < 3; i++)
				tmpvalue += dStrainP[i] * dStressP[i];
			*dGds_dFds /= tmpvalue;
			*TmpDe -= *dGds_dFds;
			//*PrinDe -= *dGds_dFds;
			//*TmpDe *= fkt1;
			for (i = 0; i < Size; i++)
				for (j = 0; j < Size; j++)
					(*TmpDe)(i, j) *= fkt1;
			//*PrinDe *= fkt1;
			//*TmpDe = *PrinDe;
			// for(i=3; i<Size; i++)
			//	(*TmpDe)(i,i) = (*Dep)(i,i);
			// delete dGds_dFds;
		}
		else
		{
			cout << "error" << prin_str[0] << "||" << prin_str[1] << "||" << prin_str[2] << endl;
			cout << "error"
			     << "csn=" << csn << "||"
			     << "Ntheta=" << Ntheta << "||"
			     << "tension=" << tension << endl;
		}
		// update prin str to normal coordinate
		for (i = 0; i < Size; i++)
			TryStress[i] = 0.;
		*ConstitutiveMatrix = (0.);

		TransMatrixA_T->multi(prin_str, TryStress); // updated stress

		/*if(fabs(prin_str[0]-prin_str[1])<MKleinsteZahl)
		    (*TMatrix)(3,3) = 0.;
		else
		    (*TMatrix)(3,3)=(dStressP[0]-dStressP[1])/(prin_str[0]-prin_str[1]);
		if(Size==6)
		{
		    if(fabs(prin_str[0]-prin_str[2])<MKleinsteZahl)
		        (*TMatrix)(4,4) = 0.;
		    else
		        (*TMatrix)(4,4)=(dStressP[0]-dStressP[2])/(prin_str[0]-prin_str[2]);
		    if(fabs(prin_str[1]-prin_str[2])<MKleinsteZahl)
		        (*TMatrix)(5,5) = 0.;
		    else
		        (*TMatrix)(5,5)=(dStressP[1]-dStressP[2])/(prin_str[1]-prin_str[2]);
		}
		for(i=3; i<Size; i++)
		    (*TmpDe)(i,i)=(1-(*TMatrix)(i,i))*(*TmpDe)(i,i);*/

		TransMatrixA_T->multi(*TmpDe, *TransMatrixA, *ConstitutiveMatrix); // updated Depc

		for (i = 3; i < Size; i++)
			(*ConstitutiveMatrix)(i, i) = (*Dep)(i, i);

		for (i = 0; i < Size; i++)
			for (j = 3; j < Size; j++)
			{
				if (i < j)
				{
					(*ConstitutiveMatrix)(i, j) = 0.;
					(*ConstitutiveMatrix)(j, i) = 0.;
				}
			}

		// mm = 0;//WX:12.2011 corecte Dep

		for (i = 0; i < Size; i++)
			for (j = 0; j < Size; j++)
				(*ConstitutiveMatrix)(i, j) = mm * (*Dep)(i, j) + (1 - mm) * (*ConstitutiveMatrix)(i, j);
		ep += sqrt(2.0 / 3.0 * (TensorMutiplication2(dStrainP, dStrainP, Dim))); // updated eff plas strain
		// ConstitutiveMatrix->Write();
		// TransMatrixA->Write();
		// TransMatrixA_T->Write();
		// if (mm<ele_val->HM_load_Factor)
		//	ele_val->HM_load_Factor = mm;//WX subincrement. not sure
	}
	else
	{
		for (i = 0; i < Size; i++)
			TryStress[i] += dstrs[i];
	}

	if (Update > 0)
	{
		(*ele_val->pStrain)(GPiGPj) = ep;
	}

	return yield;
}

// WX: calculte prin. stress and direc.
void CSolidProperties::CalPrinStrDir(double* stress, double* prin_str, double* prin_dir, int Dim)
{
	int i, j, p, q, u, w, t, s, l;
	double fm, cn, sn, omega, x, y, d;
	double eps = 1e-12, TmpValue1;
	int jt = 100;
	int n = 3;
	double a[9] = {0.}, v[9] = {0.};
	int Tmp[3] = {0}, TmpValue2;
	l = 1;
	p = 0;
	q = 0;
	for (i = 0; i < 3; i++)
		a[i + i * n] = stress[i];
	a[1] = stress[3];
	a[3] = a[1];

	if (Dim == 2)
	{
		a[2] = 0;
		a[6] = 0;
		a[5] = 0;
		a[7] = 0;
	}
	else
	{
		a[2] = stress[4];
		a[6] = a[2];
		a[5] = stress[5];
		a[7] = a[5];
	}

	for (i = 0; i <= n - 1; i++)
	{
		v[i * n + i] = 1.0;
		for (j = 0; j <= n - 1; j++)
			if (i != j)
				v[i * n + j] = 0.0;
	}
	while (1 == 1)
	{
		fm = 0.0;
		for (i = 0; i <= n - 1; i++)
			for (j = 0; j <= n - 1; j++)
			{
				d = fabs(a[i * n + j]);
				if ((i != j) && (d > fm))
				{
					fm = d;
					p = i;
					q = j;
				}
			}
		if (fm < eps)
			break;
		if (l > jt)
			break;
		l = l + 1;
		u = p * n + q;
		w = p * n + p;
		t = q * n + p;
		s = q * n + q;
		x = -a[u];
		y = (a[s] - a[w]) / 2.0;
		omega = x / sqrt(x * x + y * y);
		if (y < 0.0)
			omega = -omega;
		sn = 1.0 + sqrt(1.0 - omega * omega);
		sn = omega / sqrt(2.0 * sn);
		cn = sqrt(1.0 - sn * sn);
		fm = a[w];
		a[w] = fm * cn * cn + a[s] * sn * sn + a[u] * omega;
		a[s] = fm * sn * sn + a[s] * cn * cn - a[u] * omega;
		a[u] = 0.0;
		a[t] = 0.0;
		for (j = 0; j <= n - 1; j++)
			if ((j != p) && (j != q))
			{
				u = p * n + j;
				w = q * n + j;
				fm = a[u];
				a[u] = fm * cn + a[w] * sn;
				a[w] = -fm * sn + a[w] * cn;
			}
		for (i = 0; i <= n - 1; i++)
			if ((i != p) && (i != q))
			{
				u = i * n + p;
				w = i * n + q;
				fm = a[u];
				a[u] = fm * cn + a[w] * sn;
				a[w] = -fm * sn + a[w] * cn;
			}
		for (i = 0; i <= n - 1; i++)
		{
			u = i * n + p;
			w = i * n + q;
			fm = v[u];
			v[u] = fm * cn + v[w] * sn;
			v[w] = -fm * sn + v[w] * cn;
		}
	}

	for (i = 0; i < 3; i++)
	{
		prin_str[i] = a[i * n + i];
		Tmp[i] = i;
	}

	for (i = 0; i < 3; i++)
	{
		for (j = i + 1; j < 3; j++)
			if (prin_str[i] < prin_str[j])
			{
				TmpValue1 = prin_str[i];
				prin_str[i] = prin_str[j];
				prin_str[j] = TmpValue1;
				TmpValue2 = Tmp[i];
				Tmp[i] = Tmp[j];
				Tmp[j] = TmpValue2;
			}
	}

	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			prin_dir[j * 3 + i] = v[j * 3 + Tmp[i]];
}

// WX: calculate transform matrix between normal stress and principal stress
void CSolidProperties::CalTransMatrixA(double* v, Matrix* A, int Size)
{
	if (Size == 4)
	{
		// row 1
		(*A)(0, 0) = v[0] * v[0];
		(*A)(0, 1) = v[3] * v[3];
		(*A)(0, 2) = v[6] * v[6];
		(*A)(0, 3) = v[0] * v[3];
		// row 2
		(*A)(1, 0) = v[1] * v[1];
		(*A)(1, 1) = v[4] * v[4];
		(*A)(1, 2) = v[7] * v[7];
		(*A)(1, 3) = v[1] * v[4];
		// row 3
		(*A)(2, 0) = v[2] * v[2];
		(*A)(2, 1) = v[5] * v[5];
		(*A)(2, 2) = v[8] * v[8];
		(*A)(2, 3) = v[2] * v[5];
		// row 4
		(*A)(3, 0) = 2 * v[0] * v[1];
		(*A)(3, 1) = 2 * v[3] * v[4];
		(*A)(3, 2) = 2 * v[6] * v[7];
		(*A)(3, 3) = v[0] * v[4] + v[1] * v[3];
	}
	else
	{
		// row 1
		(*A)(0, 0) = v[0] * v[0];
		(*A)(0, 1) = v[3] * v[3];
		(*A)(0, 2) = v[6] * v[6];
		(*A)(0, 3) = v[0] * v[3];
		(*A)(0, 4) = v[0] * v[6];
		(*A)(0, 5) = v[3] * v[6];
		// row 2
		(*A)(1, 0) = v[1] * v[1];
		(*A)(1, 1) = v[4] * v[4];
		(*A)(1, 2) = v[7] * v[7];
		(*A)(1, 3) = v[1] * v[4];
		(*A)(1, 4) = v[7] * v[1];
		(*A)(1, 5) = v[4] * v[7];
		// row 3
		(*A)(2, 0) = v[2] * v[2];
		(*A)(2, 1) = v[5] * v[5];
		(*A)(2, 2) = v[8] * v[8];
		(*A)(2, 3) = v[2] * v[5];
		(*A)(2, 4) = v[2] * v[8];
		(*A)(2, 5) = v[5] * v[8];
		// row 4
		(*A)(3, 0) = 2 * v[0] * v[1];
		(*A)(3, 1) = 2 * v[3] * v[4];
		(*A)(3, 2) = 2 * v[6] * v[7];
		(*A)(3, 3) = v[0] * v[4] + v[1] * v[3];
		(*A)(3, 4) = v[6] * v[1] + v[7] * v[0];
		(*A)(3, 5) = v[3] * v[7] + v[4] * v[6];
		// row 5
		(*A)(4, 0) = 2 * v[2] * v[0];
		(*A)(4, 1) = 2 * v[5] * v[3];
		(*A)(4, 2) = 2 * v[8] * v[6];
		(*A)(4, 3) = v[2] * v[3] + v[0] * v[5];
		(*A)(4, 4) = v[8] * v[0] + v[6] * v[2];
		(*A)(4, 5) = v[5] * v[6] + v[3] * v[8];
		// row 6
		(*A)(5, 0) = 2 * v[1] * v[2];
		(*A)(5, 1) = 2 * v[4] * v[5];
		(*A)(5, 2) = 2 * v[7] * v[8];
		(*A)(5, 3) = v[1] * v[5] + v[2] * v[4];
		(*A)(5, 4) = v[7] * v[2] + v[8] * v[1];
		(*A)(5, 5) = v[4] * v[8] + v[7] * v[5];
	}
}
// calculate P for mohr coulmob return mapping
double CSolidProperties::CalVarP(double* vec1, double* vec2, double* sigma_B, double* sigma_l)
{
	return (vec1[1] * vec2[2] - vec1[2] * vec2[1]) * (sigma_B[0] - sigma_l[0])
	       + (vec1[2] * vec2[0] - vec1[0] * vec2[2]) * (sigma_B[1] - sigma_l[1])
	       + (vec1[0] * vec2[1] - vec1[1] * vec2[0]) * (sigma_B[2] - sigma_l[2]);
}
// calculate factor t for mohr coulmob
double CSolidProperties::CalVar_t(double* vecl, double* veclg, Matrix* D, double* sigma_B, double* sigma_l, int Size)
{
	double TmpVec[6] = {0.};
	double deltaStr[6] = {0.};
	double TmpVal = 0.;
	double TmpVal2 = 0.;
	D->multi(vecl, TmpVec);
	for (int i = 0; i < Size; i++)
	{
		TmpVal += veclg[i] * TmpVec[i];
		TmpVec[i] = 0.;
	}

	for (int i = 0; i < Size; i++)
		deltaStr[i] = sigma_B[i] - sigma_l[i];
	D->multi(deltaStr, TmpVec);
	for (int i = 0; i < Size; i++)
		TmpVal2 += veclg[i] * TmpVec[i];
	return TmpVal2 / TmpVal;
}
// calculate the cross product of two vectors
void CSolidProperties::VecCrossProduct(double* vec1, double* vec2, double* result_vec)
{
	result_vec[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
	result_vec[1] = -vec1[0] * vec2[2] + vec1[2] * vec2[0];
	result_vec[2] = vec1[0] * vec2[1] + vec1[1] * vec1[0];
}
// calculate constitutive matrix when return to a line
void CSolidProperties::CalDep_l(double* vecl, double* veclg, Matrix* D, Matrix* Dep_l, double /*fkt*/)
{
	double TmpVec[6] = {0.};
	double TmpVal = 0.;
	D->multi(veclg, TmpVec);
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
			(*Dep_l)(i, j) = vecl[i] * veclg[j];

		TmpVal += vecl[i] * TmpVec[i];
	}
	*Dep_l /= TmpVal;
}

void CSolidProperties::TangentialMohrShear(Matrix* Dep)
{
	*Dep = *ConstitutiveMatrix;
}

void CSolidProperties::TangentialMohrTension(Matrix* Dep)
{
	*Dep = *ConstitutiveMatrix;
}
// WX: calculate inverse matrix
void CSolidProperties::Cal_Inv_Matrix(int Size, Matrix* MatrixA, Matrix* xx)
{
	int i, j, k, jj, lk, j_col = 0; // i_row,
	double var, R;
	int L[6];
	double rhs[6];
	Matrix AA(Size, Size);
	AA = *MatrixA;

	for (i = 0; i < Size; i++)
	{
		L[i] = i;
		var = 0.0;
		for (j = 0; j < Size; j++)
		{
			if (fabs(AA(i, j)) > var)
				var = fabs(AA(i, j));
			L[i] = i;
		}

		for (j_col = 0; j_col < Size; j_col++)
			(*xx)(i, j_col) = var;
	}

	for (k = 0; k < Size - 1; k++)
	{
		var = 0.0;
		jj = 0;
		for (i = k; i < Size; i++)
		{
			R = fabs(AA(L[i], k) / (*xx)(L[i], j_col));

			if (R > var)
			{
				jj = i;
				var = R;
			}
		}
		lk = L[jj];
		L[jj] = L[k];
		L[k] = lk;

		for (i = k + 1; i < Size; i++)
		{
			var = AA(L[i], k) / AA(lk, k);

			for (j = k + 1; j < Size; j++)
				AA(L[i], j) -= var * AA(lk, j);
			AA(L[i], k) = var;
		}
	}

	for (j_col = 0; j_col < Size; j_col++)
	{
		for (i = 0; i < std::min(6, Size); i++)
		{
			rhs[i] = 0.;
			if (i == j_col)
				rhs[i] = 1.;
		}

		/* Back substituting */
		for (k = 0; k < Size - 1; k++)
			for (i = k + 1; i < Size; i++)
				rhs[L[i]] -= AA(L[i], k) * rhs[L[k]];

		(*xx)(Size - 1, j_col) = rhs[L[Size - 1]] / AA(L[Size - 1], Size - 1);
		for (i = Size - 2; i >= 0; i--)
		{
			var = rhs[L[i]];
			for (j = i + 1; j < Size; j++)
				var -= AA(L[i], j) * (*xx)(j, j_col);
			(*xx)(i, j_col) = var / AA(L[i], i);
		}
	}
}

/**************************************************************************
   ROCKFLOW - Funktion: CSolidProperties::ConsistentTangentialDP

   Aufgabe:
   Local assembly of elasto-plastic consistent tangential matrix C^ep
   (Drucker-Prager model)

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E
   double *Dep            : Consistent tangential matrix
   const double dPhi      : Plastic multiplier
   const Dim              : Space dimension
   Ergebnis:
   - double* - Matrix of 4x4, stored in an array

   Programmaenderungen:
   10/2002   WW  Erste Version
   02/2004   WW  Modification for the 3D case
   03/2007   WW  Multi-phase yield surface
**************************************************************************/
void CSolidProperties::ConsistentTangentialDP(Matrix* Dep, const double dPhi, const int Dim)
{
	double s11, s22, s12, s33, s13, s23;
	double NormX = 0.0;
	double d = 0.0;
	double c1, c2, c3, c4;
	//
	s11 = devS[0];
	s22 = devS[1];
	s33 = devS[2];
	s12 = devS[3];
	s13 = 0.0;
	s23 = 0.0;
	if (Dim == 3)
	{
		s13 = devS[4];
		s23 = devS[5];
		NormX = sqrt(s11 * s11 + s22 * s22 + s33 * s33 + 2.0 * s12 * s12 + 2.0 * s13 * s13 + 2.0 * s23 * s23);
	}
	else
		NormX = sqrt(s11 * s11 + s22 * s22 + s33 * s33 + 2.0 * s12 * s12);
	//
	if (dl2 > 0.0) // Multi-surface
	{
		c1 = Xi * BetaN * Hard * (dPhi + dl2);
		c3 = K * c1 / (c1 + 3.0 * Al * K * sqrt(dPhi * dPhi + 3.0 * Xi * Xi * (dPhi + dl2) * (dPhi + dl2)));
		c2 = 0.5 * c3 / (Xi * G * (dPhi + dl2));
		(*Dep) = 0.0;
		//
		// Row 1
		(*Dep)(0, 0) = c3 + c2 * s11;
		(*Dep)(0, 1) = c3 + c2 * s22;
		(*Dep)(0, 3) = c2 * s12;
		// Row 2
		(*Dep)(1, 0) = c3 + c2 * s11;
		(*Dep)(1, 1) = c3 + c2 * s22;
		(*Dep)(1, 3) = c2 * s12;
		// Row 3
		(*Dep)(2, 0) = c3 + c2 * s11;
		(*Dep)(2, 1) = c3 + c2 * s22;
		(*Dep)(2, 3) = c2 * s12;
		// Row 4
		if (axisymmetry || Dim == 3)
		{
			(*Dep)(0, 2) = c3 + c2 * s33;
			(*Dep)(1, 2) = c3 + c2 * s33;
			(*Dep)(2, 2) = c3 + c2 * s33;
		}
		if (Dim == 3)
		{
			(*Dep)(0, 4) = c2 * s13;
			(*Dep)(0, 5) = c2 * s23;
			(*Dep)(1, 4) = c2 * s13;
			(*Dep)(1, 5) = c2 * s23;
			(*Dep)(2, 4) = c2 * s13;
			(*Dep)(2, 5) = c2 * s23;
		}
	}
	else
	{
		//
		d = 9.0 * K * Al * Xi + 2.0 * G + BetaN * Hard * sqrt(1.0 + 3.0 * Xi * Xi);

		c1 = 2.0 * G * (1.0 - 2.0 * G * dPhi / NormX);
		c2 = (1.0 - 9.0 * Al * Xi * K / d) * K;
		c3 = -6.0 * G * K / (NormX * d);
		c4 = -4.0 * G * G * (1.0 / d - dPhi / NormX) / (NormX * NormX);

		switch (Dim)
		{
			case 2: // 2D
				// Row 1
				(*Dep)(0, 0) = 2.0 * c1 / 3.0 + c2 + c3 * (Al + Xi) * s11 + c4 * s11 * s11;
				(*Dep)(0, 1) = -c1 / 3.0 + c2 + c3 * (Xi * s22 + Al * s11) + c4 * s11 * s22;
				(*Dep)(0, 2) = 0.0; //-c1/3.0+c2+c3*(Xi*s33+Al*s11)+c4*s11*s33;
				(*Dep)(0, 3) = c3 * Xi * s12 + c4 * s11 * s12;

				// Row 2
				(*Dep)(1, 0) = -c1 / 3.0 + c2 + c3 * (Xi * s11 + Al * s22) + c4 * s11 * s22;
				(*Dep)(1, 1) = 2.0 * c1 / 3.0 + c2 + c3 * (Al + Xi) * s22 + c4 * s22 * s22;
				(*Dep)(1, 2) = 0.0; //-c1/3.0+c2+c3*(Xi*s33+Al*s22)+c4*s33*s22;
				(*Dep)(1, 3) = c3 * Xi * s12 + c4 * s22 * s12;

				// Row 3
				(*Dep)(2, 0) = 0.0; //-c1/3.0+c2+c3*(Xi*s11+Al*s33)+c4*s11*s33;
				(*Dep)(2, 1) = 0.0; // -c1/3.0+c2+c3*(Xi*s22+Al*s33)+c4*s22*s33;
				(*Dep)(2, 2) = 0.0; // 2.0*c1/3.0+c2+c3*(Al+Xi)*s33+c4*s33*s33;
				(*Dep)(2, 3) = 0.0; // c3*Xi*s12+c4*s33*s12;

				// Row 4
				(*Dep)(3, 0) = c3 * Al * s12 + c4 * s12 * s11;
				(*Dep)(3, 1) = c3 * Al * s12 + c4 * s12 * s22;
				(*Dep)(3, 2) = 0.0; // c3*Al*s12+c4*s12*s33;
				(*Dep)(3, 3) = c1 / 2.0 + c4 * s12 * s12;

				if (axisymmetry)
				{
					(*Dep)(0, 2) = -c1 / 3.0 + c2 + c3 * (Xi * s33 + Al * s11) + c4 * s11 * s33;
					// Row 2
					(*Dep)(1, 2) = -c1 / 3.0 + c2 + c3 * (Xi * s33 + Al * s22) + c4 * s33 * s22;
					// Row 3
					(*Dep)(2, 0) = -c1 / 3.0 + c2 + c3 * (Xi * s11 + Al * s33) + c4 * s11 * s33;
					(*Dep)(2, 1) = -c1 / 3.0 + c2 + c3 * (Xi * s22 + Al * s33) + c4 * s22 * s33;
					(*Dep)(2, 2) = 2.0 * c1 / 3.0 + c2 + c3 * (Al + Xi) * s33 + c4 * s33 * s33;
					(*Dep)(2, 3) = c3 * Xi * s12 + c4 * s33 * s12;
					// Row 4
					(*Dep)(3, 2) = c3 * Al * s12 + c4 * s12 * s33;
				}
				break;
			case 3: // 3D
				// Row 1
				(*Dep)(0, 0) = 2.0 * c1 / 3.0 + c2 + c3 * (Al + Xi) * s11 + c4 * s11 * s11;
				(*Dep)(0, 1) = -c1 / 3.0 + c2 + c3 * (Xi * s22 + Al * s11) + c4 * s11 * s22;
				(*Dep)(0, 2) = -c1 / 3.0 + c2 + c3 * (Xi * s33 + Al * s11) + c4 * s11 * s33;
				(*Dep)(0, 3) = c3 * Xi * s12 + c4 * s11 * s12;
				(*Dep)(0, 4) = c3 * Xi * s13 + c4 * s11 * s13;
				(*Dep)(0, 5) = c3 * Xi * s23 + c4 * s11 * s23;
				// Row 2
				(*Dep)(1, 0) = -c1 / 3.0 + c2 + c3 * (Xi * s11 + Al * s22) + c4 * s11 * s22;
				(*Dep)(1, 1) = 2.0 * c1 / 3.0 + c2 + c3 * (Al + Xi) * s22 + c4 * s22 * s22;
				(*Dep)(1, 2) = -c1 / 3.0 + c2 + c3 * (Xi * s33 + Al * s22) + c4 * s33 * s22;
				(*Dep)(1, 3) = c3 * Xi * s12 + c4 * s22 * s12;
				(*Dep)(1, 4) = c3 * Xi * s13 + c4 * s22 * s13;
				(*Dep)(1, 5) = c3 * Xi * s23 + c4 * s22 * s23;
				// Row 3
				(*Dep)(2, 0) = -c1 / 3.0 + c2 + c3 * (Xi * s11 + Al * s33) + c4 * s11 * s33;
				(*Dep)(2, 1) = -c1 / 3.0 + c2 + c3 * (Xi * s22 + Al * s33) + c4 * s22 * s33;
				(*Dep)(2, 2) = 2.0 * c1 / 3.0 + c2 + c3 * (Al + Xi) * s33 + c4 * s33 * s33;
				(*Dep)(2, 3) = c3 * Xi * s12 + c4 * s33 * s12;
				(*Dep)(2, 4) = c3 * Xi * s13 + c4 * s33 * s23;
				(*Dep)(2, 5) = c3 * Xi * s23 + c4 * s33 * s23;
				// Row 4
				(*Dep)(3, 0) = c3 * Al * s12 + c4 * s12 * s11;
				(*Dep)(3, 1) = c3 * Al * s12 + c4 * s12 * s22;
				(*Dep)(3, 2) = c3 * Al * s12 + c4 * s12 * s33;
				(*Dep)(3, 3) = c1 / 2.0 + c4 * s12 * s12;
				(*Dep)(3, 4) = c4 * s12 * s13;
				(*Dep)(3, 5) = c4 * s12 * s23;
				// Row 5
				(*Dep)(4, 0) = c3 * Al * s13 + c4 * s13 * s11;
				(*Dep)(4, 1) = c3 * Al * s13 + c4 * s13 * s22;
				(*Dep)(4, 2) = c3 * Al * s13 + c4 * s13 * s33;
				(*Dep)(4, 3) = c4 * s13 * s12;
				(*Dep)(4, 4) = c1 / 2.0 + c4 * s13 * s13;
				(*Dep)(4, 5) = c4 * s13 * s23;
				// Row 6
				(*Dep)(5, 0) = c3 * Al * s23 + c4 * s23 * s11;
				(*Dep)(5, 1) = c3 * Al * s23 + c4 * s23 * s22;
				(*Dep)(5, 2) = c3 * Al * s23 + c4 * s23 * s33;
				(*Dep)(5, 3) = c4 * s23 * s12;
				(*Dep)(5, 4) = c4 * s23 * s13;
				(*Dep)(5, 5) = c1 / 2.0 + c4 * s23 * s23;
				break;
		}
	}
}
//-------------------------------------------------------------------------
void CSolidProperties::AllocateMemoryforSYS()
{
	d2G_dSdS = new Matrix(6, 6);
	d2G_dSdM = new Matrix(6, 4);
	// Stresses, 0-5, xi, 6-10, mat, 11-17, plastci multiplier, 18.
	LocalJacobi = new Matrix(19, 19);
	inv_Jac = new Matrix(18, 18);
	sumA_Matrix = new Matrix(18, 6);
	rhs_l = new double[18];
	x_l = new double[18];
	Li = new int[18];
}

void CSolidProperties::ResizeMatricesSYS(const int Dim)
{
	if (Dim == 2)
	{
		d2G_dSdS->LimitSize(4, 4);
		d2G_dSdM->LimitSize(4, 4);
		LocalJacobi->LimitSize(15, 15);
		inv_Jac->LimitSize(14, 14);
		sumA_Matrix->LimitSize(14, 4);
	}
	else
	{
		d2G_dSdS->LimitSize(6, 6);
		d2G_dSdM->LimitSize(6, 4);
		LocalJacobi->LimitSize(19, 19);
		inv_Jac->LimitSize(18, 18);
		sumA_Matrix->LimitSize(18, 6);
	}
}

/**************************************************************************
   ROCKFLOW - Funktion: CSolidProperties::CalStress_and_TangentialMatrix_SYS

   Aufgabe:
     Using the substteping techniques to compute the integral of stresses and
   the consistent tangential matrix.
   (Weimar's model)

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E :
   long index: Elementnummer
   const int GPi, const int GPj        :   Local indeces of Gauss Points
   double *D_ep                        :   the consistent tangential matrix
   const double *De                    :   Elastic constitutive matrix
   double *dStress                     :   Incremental stresses as input
   New stress as output
   const int Update                    :   0: Do not save Gauss points values
   1: Save Gauss points values

   Ergebnis:
   - double* - Deviatoric effetive stresses, s11,s22,s12,s33

   Programmaenderungen:
   04/2003   WW  Erste Version
   02/2004   WW  Modify for the 3D case
   08/2004   WW  Set as a member of CSolidProperties
**************************************************************************/
int CSolidProperties::CalStress_and_TangentialMatrix_SYS(const int GPiGPj, const ElementValue_DM* ele_val,
                                                         const Matrix* De, Matrix* D_ep, double* dStress,
                                                         const int Update)
{
	const int LengthMat = 7;
	const int LengthStrs = std::min(std::size_t(6ul), De->Cols());
	int dim = 2;
	const int minSub = 4;
	const int maxSub = 1000;
	const double den = 0.01; // 0.01
	const int MaxIter = 1000;

	const double TolF = Tolerance_Local_Newton;

	double ErrLoc, ErrLoc0;
	int PLASTIC, subPLASTIC;

	int NPStep;
	int preSub = 1;

	int i, j, k, l, iSub, elasticSub;
	double factor, dSub;

	double I_n = 0.0;
	double I_n1 = 0.0;

	double I_p2 = 0.0;
	double I_p4 = 0.0;

	double ep = 0.0;
	double PHI = 0.0;
	double PSI = 0.0;
	double ang = 0.0;
	double dlmd = 0.0;

	const double psi1 = (*data_Plasticity)(14);
	const double Ch = (*data_Plasticity)(16);
	const double Cd = (*data_Plasticity)(17);
	const double br = (*data_Plasticity)(18);
	const double mr = (*data_Plasticity)(19);

	double Co = 0.0;
	double t1, t2, c1, c2, dW;
	double damping;
	double NormR0, NormR1; //, NormU0, NormU1;

	double F = 0.0;

	static double Stress_Inv[7];

	static double Stress_n[6];
	static double Stress_n1[6];
	static double nStress[6];
	static double xi_n[6];
	static double xi_n1[6];
	static double Mat_n[7];
	static double Mat_n1[7];
	static double supMat[7];

	static double dF_dS[6];
	static double dF_dM[7];
	static double dG_dS[6];
	const int LocDim = LocalJacobi->Cols();

	// Limits of material parameters
	// Get material parameters of the previous step
	for (i = 0; i < LengthMat; i++)
	{
		supMat[i] = (*data_Plasticity)(i + LengthMat);
		Mat_n[i] = (*ele_val->MatP)(i, GPiGPj);
	}

	// Get the converged stresses of the previous step
	for (i = 0; i < LengthStrs; i++)
		Stress_n[i] = (*ele_val->Stress)(i, GPiGPj);

	// Get the total effective plastic strain
	ep = (*ele_val->pStrain)(GPiGPj);

	// Get the rotational hardening varriables of the previous step
	for (i = 0; i < LengthStrs - 1; i++)
		xi_n[i] = (*ele_val->xi)(i, GPiGPj);
	xi_n[2] = -xi_n[0] - xi_n[1]; // x_ii*delta_ii = 0

	l = 0;

	// Total elastic try
	for (i = 0; i < LengthStrs; i++)
	{
		Stress_n1[i] = Stress_n[i] + dStress[i];
		nStress[i] = Stress_n1[i];
	}

	I_n1 = DeviatoricStress(nStress);
	for (i = 0; i < LengthStrs; i++)
	{
		nStress[i] -= I_n1 * xi_n[i] / 3.0;
		xi_n1[i] = xi_n[i];
	}

	for (i = 0; i < LengthMat; i++)
		Mat_n1[i] = Mat_n[i];

	Stress_Inv[1] = 0.5 * TensorMutiplication2(nStress, nStress, dim);
	Stress_Inv[2] = TensorMutiplication3(nStress, nStress, nStress, dim);

	I_p2 = I_n1 * I_n1;
	I_p4 = I_p2 * I_p2;

	if (Stress_Inv[1] > 0.0)
	{
		ang = Stress_Inv[2] / MathLib::fastpow(sqrt(Stress_Inv[1]), 3);

		PHI = Stress_Inv[1] * pow(1 + Mat_n[5] * ang, Mat_n[6]) + 0.5 * Mat_n[0] * I_p2 + Mat_n[2] * Mat_n[2] * I_p4;

		/* Compute yield surface */
		F = sqrt(PHI) + Mat_n[1] * I_n1 + Mat_n[3] * I_p2 - Mat_n[4];
	}
	else
		F = -1.0;

	if (pcs_deformation == 1)
		F = -1.0;

	PLASTIC = 0;
	if (F > TolF && !PreLoad) /* In Yield Status */
	{
		PLASTIC = 1;
		subPLASTIC = 0;

		/* Compute the ratio of elastic increment */
		I_n = DeviatoricStress(Stress_n);
		for (i = 0; i < LengthStrs; i++)
			Stress_n[i] -= I_n * xi_n[i] / 3.0;

		/* Compute the number of sub step */
		Co = sqrt(2.0 * Stress_Inv[1] * TensorMutiplication2(Stress_n, Stress_n, dim));
		if (Co > 0.0)
			preSub = 1 + (int)(acos(TensorMutiplication2(Stress_n, nStress, dim) / Co) / den);
		else
			preSub = minSub;
		if (preSub > maxSub)
			preSub = maxSub;
		if (preSub < minSub)
			preSub = minSub;

		// Recover stress at n_th step
		for (i = 0; i < LengthStrs; i++)
			Stress_n[i] += I_n * xi_n[i] / 3.0;
		for (i = 0; i < 3; i++) // principle stress
			Stress_n[i] += I_n / 3.0;

		dSub = 1.0 / (float)preSub;

		// Approximate the yield surface
		elasticSub = 0;
		// The previous solution as the initial approximation
		for (i = 0; i < LengthMat; i++)
			Mat_n1[i] = Mat_n[i];
		for (iSub = 0; iSub < preSub; iSub++)
		{
			// Elastic try in each sub step
			for (i = 0; i < LengthStrs; i++)
			{
				// Elastic try in each sub step as the initial approximation
				Stress_n1[i] = Stress_n[i] + (iSub + 1.0) * dSub * dStress[i];
				// The previous solution as the initial approximation
				xi_n1[i] = xi_n[i];
				nStress[i] = Stress_n1[i];
			}
			I_n1 = DeviatoricStress(nStress);
			for (i = 0; i < LengthStrs; i++)
				nStress[i] -= I_n1 * xi_n1[i] / 3.0;

			Stress_Inv[0] = I_n1;
			Stress_Inv[1] = 0.5 * TensorMutiplication2(nStress, nStress, dim);
			Stress_Inv[2] = TensorMutiplication3(nStress, nStress, nStress, dim);

			I_p2 = I_n1 * I_n1;
			I_p4 = I_p2 * I_p2;
			ang = Stress_Inv[2] / MathLib::fastpow(sqrt(Stress_Inv[1]), 3);
			PHI = sqrt(Stress_Inv[1] * pow(1 + Mat_n1[5] * ang, Mat_n1[6]) + 0.5 * Mat_n1[0] * I_p2
			           + Mat_n1[2] * Mat_n1[2] * I_p4);
			// The yield function
			F = PHI + Mat_n1[1] * I_n1 + Mat_n1[3] * I_p2 - Mat_n1[4];
			if (F > TolF)
			{
				elasticSub = iSub;
				Co = dSub * (float)elasticSub;

				/*
				   for(i=0; i<LengthStrs; i++)
				   {
				   Stress_n[i] += Co*dStress[i];
				   dStress[i] *= 1-Co;
				   }
				 */
				break;
			}
		}

		subPLASTIC = 0;
		/* Compute stresses by substepping  */
		// for(iSub=0; iSub<preSub; iSub++)
		for (iSub = elasticSub; iSub < preSub; iSub++)
		{
			factor = dSub;
			if (iSub == elasticSub)
				factor += Co;
			/* Elastic try in each sub step */
			for (i = 0; i < LengthStrs; i++)
			{
				/* Elastic try in each sub step as the initial approximation*/
				Stress_n1[i] = Stress_n[i] + factor * dStress[i];
				/* The previous solution as the initial approximation */
				xi_n1[i] = xi_n[i];

				nStress[i] = Stress_n1[i];
			}

			I_n1 = DeviatoricStress(nStress);
			for (i = 0; i < LengthStrs; i++)
				nStress[i] -= I_n1 * xi_n1[i] / 3.0;

			/* The previous solution as the initial approximation */
			for (i = 0; i < LengthMat; i++)
				Mat_n1[i] = Mat_n[i];

			/* Local Newton-Raphson method */
			ErrLoc0 = 1.0e8;
			ErrLoc = 1.0e8;
			// if(iSub==0)
			dlmd = 0.0;
			NPStep = 0;
			// NormU0 = 1.0;
			NormR0 = 1.0;
			NormR1 = 1.0;
			while (ErrLoc > Tolerance_Local_Newton)
			{
				NPStep++;

				Stress_Inv[0] = I_n1;
				Stress_Inv[1] = 0.5 * TensorMutiplication2(nStress, nStress, dim);
				Stress_Inv[2] = TensorMutiplication3(nStress, nStress, nStress, dim);

				I_p2 = I_n1 * I_n1;
				I_p4 = I_p2 * I_p2;
				ang = Stress_Inv[2] / MathLib::fastpow(sqrt(Stress_Inv[1]), 3);
				PHI = sqrt(Stress_Inv[1] * pow(1 + Mat_n1[5] * ang, Mat_n1[6]) + 0.5 * Mat_n1[0] * I_p2
				           + Mat_n1[2] * Mat_n1[2] * I_p4);
				PSI = sqrt(Stress_Inv[1] / psi1 + 0.5 * Mat_n1[0] * I_p2 + Mat_n1[2] * Mat_n1[2] * I_p4);

				/* The yield function */
				F = PHI + Mat_n1[1] * I_n1 + Mat_n1[3] * I_p2 - Mat_n1[4];

				if (F <= TolF && NPStep == 1 && subPLASTIC == 0)
					break;
				else
				{
					if (NPStep == 1)
						subPLASTIC++;
					Stress_Inv[3] = ang;
					Stress_Inv[4] = PHI;
					Stress_Inv[5] = PSI;
					Stress_Inv[6] = PSI * PSI * PSI;

					dG_dNStress(dG_dS, nStress, Stress_Inv, Mat_n1, LengthStrs);
					dG__dNStress_dNStress(nStress, Stress_Inv, Mat_n1, LengthStrs);
					dG_dSTress_dMat(nStress, Stress_Inv, Mat_n1, LengthStrs);

					//--------------- Local Jacibin ------------------
					(*LocalJacobi) = 0.0;

					// dr_stress/d... --------------------------
					for (i = 0; i < LengthStrs; i++)
					{
						l = i;
						// dr_stress/dxi and dr_stress/dM  -------
						for (j = 0; j < LengthStrs; j++)
						{
							c1 = 1.0;
							if (j > 2)
								c1 = 2.0;
							t1 = 0.0;
							t2 = 0.0;
							for (k = 0; k < LengthStrs; k++)
							{
								c2 = 1.0;

								if (k > 2)
									c2 = 2.0;
								t1 += c1 * c2 * (*De)(i, k) * (*d2G_dSdS)(k, j);
								if (j >= 4)
									continue;
								t2 += c2 * (*De)(i, k) * (*d2G_dSdM)(k, j);
							}

							// dG_dStress_dXi
							if (br > 0.0)
							{
								if (j == 2) // xi_33=-(xi_11+xi_22+x_33)
								{
									(*LocalJacobi)(l, LengthStrs) += dlmd * t1 * I_n1 / 3.0;
									(*LocalJacobi)(l, LengthStrs + 1) += dlmd * t1 * I_n1 / 3.0;
								}
								else if (j < 2)
									(*LocalJacobi)(l, LengthStrs + j) += -dlmd * t1 * I_n1 / 3.0;
								else if (j > 2)
									(*LocalJacobi)(l, LengthStrs + j - 1) += -dlmd * t1 * I_n1 / 3.0;
							}

							// dG_dStress_dMat
							if ((Ch > 0.0) && (Cd > 0.0))
								(*LocalJacobi)(l, 2 * LengthStrs - 1 + j) += dlmd * t2;
						}

						// dG_dStress_dLambda
						t1 = 0.0;
						for (k = 0; k < LengthStrs; k++)
						{
							c2 = 1.0;
							if (k > 2)
								c2 = 2.0;
							t1 += c2 * (*De)(i, k) * dG_dS[k];
						}
						(*LocalJacobi)(l, LocDim - 1) += t1;
					}

					dW = 0.0;
					for (k = 0; k < LengthStrs; k++)
					{
						c2 = 1.0;
						if (k > 2)
							c2 = 2.0;
						dW += c2 * Stress_n1[k] * dG_dS[k];
					}

					// dr_M/d... --------------------------
					for (i = 0; i < LengthMat; i++)
					{
						l = i + 2 * LengthStrs - 1;
						Co = Ch;
						if (i > 4)
							Co = Cd;
						for (j = 0; j < LengthMat; j++)
							if (i == j)
								(*LocalJacobi)(l, 2 * LengthStrs - 1 + j) += 1.0 + dlmd * Co * dW;

						// dr_m/dnStress --------------------------
						for (j = 0; j < LengthStrs; j++)
						{
							c1 = 1.0;
							if (j > 2)
								c1 = 2.0;
							t1 = 0.0;
							for (k = 0; k < LengthStrs; k++)
							{
								c2 = 1.0;
								if (k > 2)
									c2 = 2.0;
								t1 += c1 * c2 * Stress_n1[k] * (*d2G_dSdS)(k, j);
							}
							t2 = -dlmd * Co * (supMat[i] - Mat_n1[i]) * t1;

							// d(.)/dxi
							if (br > 0.0)
							{
								if (j == 2) // xi_33 = -(xi_11+xi_22)
								{
									(*LocalJacobi)(l, LengthStrs) += I_n1 * t2 / 3.0;
									(*LocalJacobi)(l, LengthStrs + 1) += I_n1 * t2 / 3.0;
								}
								else if (j < 2)
									(*LocalJacobi)(l, LengthStrs + j) += -I_n1 * t2 / 3.0;
								else if (j > 2)
									(*LocalJacobi)(l, LengthStrs + j - 1) += -I_n1 * t2 / 3.0;
							}
						}

						// dr_m/dMat --------------------------
						for (j = 0; j < LengthStrs; j++)
						{
							t1 = 0.0; // Stress times d2G_dSdM
							for (k = 0; k < 4; k++)
							{
								c1 = 1.0;
								if (k > 2)
									c1 = 2.0;
								t1 += c1 * Stress_n1[k] * (*d2G_dSdM)(k, j);
							}
							t2 = (supMat[i] - Mat_n1[i]) * t1;
							(*LocalJacobi)(l, 2 * LengthStrs - 1 + j) += -dlmd * Co * t2;
						}

						// dr_M/dLambda --------------------------
						(*LocalJacobi)(l, LocDim - 1) += -Co * (supMat[i] - Mat_n1[i]) * dW;
					}

					// d(r_stress)__dstress --------------------------
					dG__dStress_dStress(nStress, xi_n1, Stress_Inv, Mat_n1, LengthStrs);
					for (i = 0; i < LengthStrs; i++)
					{
						l = i;
						for (j = 0; j < LengthStrs; j++)
						{
							if (i == j)
								(*LocalJacobi)(l, j) += 1.0;
							t1 = 0.0;
							t2 = 0.0;
							for (k = 0; k < LengthStrs; k++)
							{
								c2 = 1.0;
								if (k > 2)
								{
									if (j > 2)
										c2 = 4.0;
									else
										c2 = 2.0;
								}
								t1 += c2 * (*De)(i, k) * (*d2G_dSdS)(k, j);
							}
							// dG_dStress_dStress
							(*LocalJacobi)(i, j) += dlmd * t1;
						}
					}

					// d(r_M)__dstress --------------------------
					for (i = 0; i < LengthMat; i++)
					{
						Co = Ch;
						if (i > 4)
							Co = Cd;
						l = i + 2 * LengthStrs - 1;
						// df3_dStress --------------------------
						for (j = 0; j < LengthStrs; j++)
						{
							c1 = 1.0;
							if (j > 2)
								c1 = 2.0;
							t1 = 0.0;
							for (k = 0; k < LengthStrs; k++)
							{
								c2 = 1.0;
								if (k > 2)
									c2 = 2.0;
								t1 += c1 * c2 * Stress_n1[k] * (*d2G_dSdS)(k, j);
							}
							t2 = -dlmd * Co * (supMat[i] - Mat_n1[i]) * (c1 * dG_dS[j] + t1);
							(*LocalJacobi)(l, j) += t2;
						}
					}

					// d(r_xi)/d... --------------------------
					// d2G_dSdS: df2_dS. d2G_dSdM: df2_dXi.
					dfun2(nStress, xi_n1, Stress_Inv, Mat_n1, LengthStrs);
					for (i = 0; i < LengthStrs; i++)
					{
						if (i < 2)
							l = i + LengthStrs;
						if (i == 2)
							continue;
						if (i > 2)
							l = i - 1 + LengthStrs;

						for (j = 0; j < LengthStrs; j++)
						{
							c1 = 1.0;
							if (j > 2)
								c1 = 2.0;
							(*LocalJacobi)(l, j) += c1 * dlmd * br * (*d2G_dSdS)(i, j) / psi1;
							if (j == 2) // x_ii*delta_ii = 0;
							{
								(*LocalJacobi)(l, LengthStrs) -= c1 * dlmd * br * (*d2G_dSdM)(i, j) / psi1;
								(*LocalJacobi)(l, LengthStrs + 1) -= c1 * dlmd * br * (*d2G_dSdM)(i, j) / psi1;
							}
							if (j < 2)
								(*LocalJacobi)(l, LengthStrs + j) += c1 * dlmd * br * (*d2G_dSdM)(i, j) / psi1;
							if (j > 2)
								(*LocalJacobi)(l, LengthStrs + j) += c1 * dlmd * br * (*d2G_dSdM)(i, j - 1) / psi1;

							if (i == j) // dxi/dxi
							{
								if (j < 2)
									(*LocalJacobi)(l, LengthStrs + j) += 1.0;
								if (j == 2)
									continue;
								if (j > 2)
									(*LocalJacobi)(l, LengthStrs + j - 1) += 1.0;
							}
						}
						// d(r_xi)/dMat
						t1 = sqrt(Stress_Inv[1] / 3.0);
						(*LocalJacobi)(l, 2 * LengthStrs - 1) -= 0.25 * dlmd * br * t1 * I_n1
						                                         * (I_n1 * xi_n1[i] - mr * nStress[i])
						                                         / (Stress_Inv[6] * psi1 * I_n1);
						(*LocalJacobi)(l, 2 * LengthStrs + 1) -= dlmd * br * t1 * Mat_n1[2] * I_p4
						                                         * (I_n1 * xi_n1[i] - mr * nStress[i])
						                                         / (Stress_Inv[6] * psi1 * I_n1);
						// d(r_xi)/dLambda
						(*LocalJacobi)(l, LocDim - 1)
						    -= br * t1 * (I_n1 * xi_n1[i] - mr * nStress[i]) / (PSI * psi1 * I_n1);
					}

					// d(r_F)/d...  --------------------------
					dF_dNStress(dF_dS, nStress, Stress_Inv, Mat_n1, LengthStrs);
					dF_dMat(dF_dM, Stress_Inv, Mat_n1);
					l = LocDim - 1;
					for (j = 0; j < LengthStrs; j++)
					{
						c1 = 1.0;
						if (j > 2)
							c1 = 2.0;
						// dr_F/dxi
						if (br > 0.0)
						{
							if (j == 2)
							{
								(*LocalJacobi)(l, LengthStrs) += c1 * I_n1 * dF_dS[j] / 3.0;
								(*LocalJacobi)(l, LengthStrs + 1) += c1 * I_n1 * dF_dS[j] / 3.0;
							}
							else if (j < 2)
								(*LocalJacobi)(l, LengthStrs + j) += -c1 * I_n1 * dF_dS[j] / 3.0;
							else if (j > 2)
								(*LocalJacobi)(l, LengthStrs + j - 1) += -c1 * I_n1 * dF_dS[j] / 3.0;
						}
					}

					// dr_F/dStress
					dF_dStress(dF_dS, xi_n1, Stress_Inv, Mat_n1, LengthStrs);
					for (j = 0; j < LengthStrs; j++)
					{
						c1 = 1.0;
						if (j > 2)
							c1 = 2.0;
						(*LocalJacobi)(l, j) += c1 * dF_dS[j];
					}

					// dr_F/dM
					if ((Ch > 0.0) && (Cd > 0.0))
						for (j = 0; j < LengthMat; j++)
							(*LocalJacobi)(l, 2 * LengthStrs - 1 + j) += dF_dM[j];

					//-------- End Local Jacibin ----------

					//-------- RHS ------------
					for (i = 0; i < LengthStrs; i++)
					{
						rhs_l[i] = Stress_n1[i];
						for (k = 0; k < LengthStrs; k++)
						{
							c2 = 1.0;
							if (k > 2)
								c2 = 2.0;
							rhs_l[i] += c2 * dlmd * (*De)(i, k) * dG_dS[k];
						}
						rhs_l[i] -= Stress_n[i] + factor * dStress[i];

						// For Xi
						if (i == 2)
							continue;
						if (i < 2)
							l = LengthStrs + i;
						else if (i > 2)
							l = LengthStrs + i - 1;
						rhs_l[l] = xi_n1[i] - xi_n[i]
						           - dlmd * br * sqrt(Stress_Inv[1] / 3.0) * (I_n1 * xi_n1[i] - mr * nStress[i])
						                 / (PSI * psi1 * I_n1);
					}
					for (i = 0; i < LengthMat; i++)
					{
						Co = Ch;
						if (i > 4)
							Co = Cd;
						rhs_l[2 * LengthStrs - 1 + i] = Mat_n1[i] - Mat_n[i] - Co * dlmd * (supMat[i] - Mat_n1[i]) * dW;
					}
					// The yield function
					rhs_l[LocDim - 1] = F;
					//-------- End RHS -----------------

					//--------- Compute the error of the residual  --------
					if (NPStep == 1)
					{
						NormR0 = 0.0;
						for (i = 0; i < LocDim; i++)
							NormR0 += rhs_l[i] * rhs_l[i];
					}

					if (sqrt(NormR0) < TolF)
						ErrLoc = 0.01 * TolF;
					else
					{
						NormR1 = 0.0;
						for (i = 0; i < LocDim; i++)
							NormR1 += rhs_l[i] * rhs_l[i];
						ErrLoc = sqrt(NormR1 / NormR0);
					}
					if (fabs(F) < TolF || F < 0.0)
						ErrLoc = 0.01 * TolF;
					if (NormR1 < TolF)
						ErrLoc = 0.01 * TolF;
					//------ End: Compute the error of the residual

					damping = 1.0;
					if (NPStep > 1)
						if ((ErrLoc / ErrLoc0) > 0.01)
							damping = 0.5;
					if (NPStep > MaxIter - 1)
						damping = 0.2;

					//------  Solve the linear equation
					Gauss_Elimination(LocDim, *LocalJacobi, Li, x_l);
					Gauss_Back(LocDim, *LocalJacobi, rhs_l, Li, x_l);
					//------  End Solve the linear equation

					//------ Compute the error of the solution
					/*
					   if(NPStep==1)
					   {
					   NormU0 = 0.0;
					   for(i=0; i<LocDim; i++)
					   {
					      NormU0 += x_l[i]*x_l[i];
					   }
					   }

					   if(sqrt(NormU0)<TolLocalNewT) ErrLoc = 0.1*TolLocalNewT;
					   else
					   {
					   NormU1 = 0.0;
					   for(i=0; i<LocDim; i++)
					   {
					   NormU1 += x_l[i]*x_l[i];
					   }
					   ErrLoc = min(ErrLoc, sqrt(NormU1/NormU0));
					   }
					   if(fabs(F)<TolLocalNewT) ErrLoc = 0.1*TolLocalNewT;
					 */
					//----- End: Compute the error of the solution
					ErrLoc0 = ErrLoc;

					//----- Update the Newton-Raphson step
					for (i = 0; i < LocDim; i++)
						x_l[i] *= damping;

					for (i = 0; i < LengthStrs; i++)
					{
						Stress_n1[i] -= x_l[i];
						nStress[i] = Stress_n1[i];
						if (i == 2)
							continue;
						else if (i < 2)
							xi_n1[i] -= x_l[i + LengthStrs];
						else if (i > 2)
							xi_n1[i] -= x_l[i - 1 + LengthStrs];
					}
					I_n1 = DeviatoricStress(nStress);
					xi_n1[2] = -xi_n1[0] - xi_n1[1];

					for (i = 0; i < LengthStrs; i++)
						nStress[i] -= I_n1 * xi_n1[i] / 3.0;

					for (i = 0; i < LengthMat; i++)
						Mat_n1[i] -= x_l[i + 2 * LengthStrs - 1];

					dlmd -= x_l[LocDim - 1];
				} // End if (F>0.0)

				if (NPStep > MaxIter)
				{
					printf("\n ~O~ ~O~  Local Newton-Raphson has problem in convergence in MatCalcStressWM().. \n");
					printf("\n F = %g\n", F);
					break; // abort();
				}
			} // End of local Newton-Raphson

			// Accumulated plastic strain
			if (subPLASTIC > 0)
				ep += dlmd * sqrt(2.0 * TensorMutiplication2(dG_dS, dG_dS, dim) / 3.0);

			//-------  Compute the consistent tangential matrix
			if (Update <= 0 && subPLASTIC > 0)
			{
				//--- 1.  Compute the inverse of the Jacobian -
				for (j = 0; j < LocDim - 1; j++)
				{
					for (i = 0; i < LocDim; i++)
					{
						rhs_l[i] = 0.0;
						if (i == j)
							rhs_l[i] = 1.0;
					}
					// the i_th column of the invJac matrix
					Gauss_Back(LocDim, *LocalJacobi, rhs_l, Li, x_l);
					for (i = 0; i < LocDim - 1; i++)
						(*inv_Jac)(i, j) = x_l[i];
				}

				//- 2.  A*A*A*... -
				if (subPLASTIC == 1) // First substep, which has plasticity

					for (i = 0; i < LocDim - 1; i++)
						for (j = 0; j < LengthStrs; j++)
							(*sumA_Matrix)(i, j) = (*inv_Jac)(i, j) * factor;

				else
				{
					for (i = 0; i < LocDim - 1; i++)
						for (j = 0; j < LengthStrs; j++)
							if (i == j)
								(*sumA_Matrix)(i, j) += factor;

					LocalJacobi->LimitSize(LocDim - 1, LengthStrs);
					(*LocalJacobi) = 0.0;
					inv_Jac->multi(*sumA_Matrix, *LocalJacobi);
					for (i = 0; i < LocDim - 1; i++)
						for (j = 0; j < LengthStrs; j++)
							(*sumA_Matrix)(i, j) = (*LocalJacobi)(i, j);

					LocalJacobi->LimitSize(LocDim, LocDim);
				}
				//- 3.  D_ep -
				if (iSub == preSub - 1)
				{
					for (i = 0; i < LengthStrs; i++)
						for (j = 0; j < LengthStrs; j++)
						{
							(*D_ep)(i, j) = 0.0;
							for (k = 0; k < LengthStrs; k++)
								(*D_ep)(i, j) += (*sumA_Matrix)(i, k) * (*De)(k, j);
						}
				}
			}
			///  End Compute the consistent tangential matrix

			// Update the substep
			for (i = 0; i < LengthStrs; i++)
			{
				Stress_n[i] = Stress_n1[i];
				xi_n[i] = xi_n1[i];
			}

			for (i = 0; i < LengthStrs; i++)
				Mat_n[i] = Mat_n1[i];
		} // End of Compute stresses by substepping
	} // If F>0.0

	// Save the current stresses
	if (Update > 0)
	{
		(*ele_val->pStrain)(GPiGPj) = ep;
		// for(i=0; i<LengthStrs; i++)
		//    (*ele_val->Stress)(i, GPiGPj) = Stress_n1[i];
		for (i = 0; i < LengthStrs - 1; i++)
			(*ele_val->xi)(i, GPiGPj) = xi_n1[i];
		for (i = 0; i < LengthMat; i++)
			(*ele_val->MatP)(i, GPiGPj) = Mat_n1[i];
	}
	//   else
	//   {  // New stresses passed through dStress for the residual computation
	for (i = 0; i < LengthStrs; i++)
		dStress[i] = Stress_n1[i];
	//   }

	/*
	   // Already considered in element_dm, ele_val_dm
	   // For the case of the initial stress being given, the contribution of
	   // the initial stress to the right hand side should be taken into account
	   dStress[0] -= (*data_Plasticity)(20); // Initial stress_xx
	   dStress[1] -= (*data_Plasticity)(21); // Initial stress_yy
	   dStress[2] -= (*data_Plasticity)(22); // Initial stress_zz
	 */

	return PLASTIC;
}

/**************************************************************************
   ROCKFLOW - Funktion: CSolidProperties::dF_dNStress

   Aufgabe:
   Computing the derivatives of the yield function with respect to
   the normalized stresses (Weimar's model)

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E :
       double *dFdS               : The derivatives of the yield function
   with respect to stresses
   const long EleIndex        : Element index
   const double *DevS         : The deviatoric stresses
   const double *S_Invariants : Stress invariants
   const double *MatN1        : Material parameters

   Ergebnis:

   Programmaenderungen:
   08/2003   WW  Erste Version (for 2D)
   02/2004   WW  Modification for the 3D case
   08/2004   WW  Set as a member of CSolidProperties
 ***************************************************************************/
void CSolidProperties::dF_dNStress(double* dFdS, const double* DevS, const double* S_Invariants, const double* MatN1,
                                   const int LengthStrs)
{
	int i;
	double In1 = S_Invariants[0];
	double I_p3 = In1 * In1 * In1;
	double ang = S_Invariants[3];
	double PHI = S_Invariants[4];

	/* 1. Compute the derivatives of the yield function with respect to stresses */
	/* 1.1 s_{jk}s_{ki} */
	if (LengthStrs == 4) // 2D
	{
		// s11*s11+s12*s21
		dFdS[0] = DevS[0] * DevS[0] + DevS[3] * DevS[3];
		// s21*s12+s22*s22
		dFdS[1] = DevS[3] * DevS[3] + DevS[1] * DevS[1];
		dFdS[2] = DevS[2] * DevS[2]; // s33*s33
		// s11*s12+s12*s22
		dFdS[3] = DevS[0] * DevS[3] + DevS[3] * DevS[1];
	}
	else
	{
		// s11*s11+s12*s21+s13*s31
		dFdS[0] = DevS[0] * DevS[0] + DevS[3] * DevS[3] + DevS[4] * DevS[4];
		// s21*s12+s22*s22+s23*s32
		dFdS[1] = DevS[3] * DevS[3] + DevS[1] * DevS[1] + DevS[5] * DevS[5];
		// s31*s13+s32*s23+s33*s33
		dFdS[2] = DevS[4] * DevS[4] + DevS[5] * DevS[5] + DevS[2] * DevS[2];
		// s11*s12+s12*s22+s13*s32
		dFdS[3] = DevS[0] * DevS[3] + DevS[3] * DevS[1] + DevS[4] * DevS[5];
		// s11*s13+s12*s23+s13*s33
		dFdS[4] = DevS[0] * DevS[4] + DevS[3] * DevS[5] + DevS[4] * DevS[2];
		// s21*s13+s22*s23+s23*s33
		dFdS[5] = DevS[3] * DevS[4] + DevS[1] * DevS[5] + DevS[5] * DevS[2];
	}
	/* 1.2 dtheta/ds*/
	for (i = 0; i < LengthStrs; i++)
	{
		dFdS[i] -= 1.5 * S_Invariants[2] * DevS[i] / S_Invariants[1];
		if (i < 3)
			dFdS[i] -= 2.0 * S_Invariants[1] / 3.0;
		dFdS[i] /= MathLib::fastpow(sqrt(S_Invariants[1]), 3);
	}
	/* 1.3 dF/ds   ..  */
	for (i = 0; i < LengthStrs; i++)
	{
		dFdS[i] *= MatN1[5] * MatN1[6] * S_Invariants[1] * pow(1 + MatN1[5] * ang, MatN1[6] - 1.0);
		dFdS[i] += pow(1 + MatN1[5] * ang, MatN1[6]) * DevS[i];
		dFdS[i] /= (2 * PHI);
		if (i < 3)
		{
			dFdS[i] += (MatN1[0] * In1 + 4.0 * MatN1[2] * MatN1[2] * I_p3) / (2 * PHI);
			dFdS[i] += MatN1[1] + 2.0 * MatN1[3] * In1;
		}
	}
}

/**************************************************************************
   ROCKFLOW - Funktion: CSolidProperties::dF_dStress

   Aufgabe:
   Computing the derivatives of the yield function with respect to
   stresses (Weimar's model)
   = dF_dNStress- dF_dStress;

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E :
   double *dFdS               : The derivatives of the yield function
   with respect to stresses
   const double *RotV         : Rotational variables
   const double *S_Invariants : Stress invariants
   const double *MatN1        : Material parameters

   Ergebnis:

   Programmaenderungen:
   08/2003   WW  Erste Version (for 2D)
   02/2004   WW  Modification for the 3D case
   08/2004   WW  Set as a member of CSolidProperties
 ***************************************************************************/
void CSolidProperties::dF_dStress(double* dFdS, const double* RotV, const double* S_Invariants, const double* MatN1,
                                  const int LengthStrs)
{
	int i;
	double In1 = S_Invariants[0];
	double I_p3 = In1 * In1 * In1;
	double PHI = S_Invariants[4];

	for (i = 0; i < LengthStrs; i++)
		dFdS[i] -= ((MatN1[0] * In1 + 4.0 * MatN1[2] * MatN1[2] * I_p3) / (2 * PHI) + MatN1[1] + 2.0 * MatN1[3] * In1)
		           * RotV[i];
}

/**************************************************************************
   ROCKFLOW - Funktion: CSolidProperties::dF_dMat

   Aufgabe:
   Computing the derivatives of the yield function with respect to
   material parameters (Weimar's model)

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E :
       double *dFdM               : The derivatives of the yield function
   with respect to material parameters
   const long EleIndex        : Element index
   const double *DevS         : The deviatoric stresses
   const double *S_Invariants : Stress invariants
   const double *MatN1        : Material parameters

   Ergebnis:

   Programmaenderungen:
   08/2003   WW  Erste Version (for 2D)
   08/2004   WW  Set as a member of CSolidProperties
 ***************************************************************************/
void CSolidProperties::dF_dMat(double* dFdM, const double* S_Invariants, const double* MatN1)
{
	double In1 = S_Invariants[0];
	double I_p2 = In1 * In1;
	double I_p4 = I_p2 * I_p2;
	double ang = S_Invariants[3];
	double PHI = S_Invariants[4];
	/*dF_dAlpha*/
	dFdM[0] = 0.25 * I_p2 / PHI;
	/*dF_dBeta*/
	dFdM[1] = In1;
	/*dF_dDelta*/
	dFdM[2] = MatN1[2] * I_p4 / PHI;
	/*dF_dEpsilon*/
	dFdM[3] = I_p2;
	/*dF_dKappa*/
	dFdM[4] = -1.0;
	/*dF_dGamma*/
	dFdM[5] = 0.5 * S_Invariants[1] * MatN1[6] * ang * pow(1.0 + MatN1[5] * ang, MatN1[6] - 1.0) / PHI;
	/*dF_dM*/
	dFdM[6] = 0.5 * S_Invariants[1] * pow(1.0 + MatN1[5] * ang, MatN1[6]) * log(1.0 + MatN1[5] * ang) / PHI;
}

/**************************************************************************
   ROCKFLOW - Funktion: CSolidProperties::dG_dNStress

   Aufgabe:
   Computing the derivatives of the plastic potential with respect to
   normal stresses (Weimar's model)

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E :
       double *dGdS               : The derivatives of the plastic potential
   with respect to stresses
   const long EleIndex        : Element index
   const double *DevS         : The deviatoric stresses
   const double *S_Invariants : Stress invariants
   const double *MatN1        : Material parameters

   Ergebnis:

   Programmaenderungen:
   08/2003   WW  Erste Version (for 2D)
   02/2004   WW  Modification for the 3D case
   08/2004   WW  Set as a member of CSolidProperties
 ***************************************************************************/
void CSolidProperties::dG_dNStress(double* dGdS, const double* DevS, const double* S_Invariants, const double* MatN1,
                                   const int LengthStrs)
{
	int i;
	const double psi1 = (*data_Plasticity)(14);
	const double psi2 = (*data_Plasticity)(15);
	double In1 = S_Invariants[0];
	double I_p3 = In1 * In1 * In1;
	double PSI = S_Invariants[5];

	for (i = 0; i < LengthStrs; i++)
	{
		dGdS[i] = 0.5 * DevS[i] / PSI / psi1;
		if (i < 3)
			dGdS[i] += 0.5 * (MatN1[0] * In1 + 4.0 * MatN1[2] * MatN1[2] * I_p3) / PSI + (1.0 + psi2) * MatN1[1]
			           + 2.0 * MatN1[3] * In1;
	}
}

/**************************************************************************
   ROCKFLOW - Funktion: CSolidProperties::dG_dNStress_dNStress

   Aufgabe:
   Computing the second derivatives of the plastic potential with respect to
   normal stresses (Weimar's model)

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E :
       double *dG_dSdS            : The second derivatives of the plastic potential
   with respect to stresses
   const long EleIndex        : Element index
   const double *DevS         : The deviatoric stresses
   const double *S_Invariants : Stress invariants
   const double *MatN1        : Material parameters

   Ergebnis:

   Programmaenderungen:
   08/2003   WW  Erste Version (for 2D)
   02/2004   WW  Extended to 3D
   08/2004   WW  Set as a member of CSolidProperties
 ***************************************************************************/
void CSolidProperties::dG__dNStress_dNStress(const double* DevS,
                                             const double* S_Invariants,
                                             const double* MatN1,
                                             const int LengthStrs)
{
	int i, j;
	int ii, jj, kk, ll;
	double delta_ij_kl;
	const double psi1 = (*data_Plasticity)(14);
	double In1 = S_Invariants[0];
	double I_p2 = In1 * In1;
	double I_p3 = In1 * In1 * In1;
	double PSI = S_Invariants[5];
	double PSI_p3 = S_Invariants[6];

	for (i = 0; i < LengthStrs; i++)
	{
		ii = i;
		jj = i;
		if (i == 3) // s_12
		{
			ii = 1;
			jj = 2;
		}
		// For 3D
		else if (i == 4) // s_13
		{
			ii = 1;
			jj = 3;
		}
		else if (i == 5) // s_23
		{
			ii = 2;
			jj = 3;
		}

		for (j = 0; j < LengthStrs; j++)
		{
			kk = j;
			ll = j;
			if (j == 3)
			{
				kk = 1;
				ll = 2;
			}
			// For 3D
			else if (j == 4) // s_13
			{
				kk = 1;
				ll = 3;
			}
			else if (j == 5) // s_23
			{
				kk = 2;
				ll = 3;
			}

			delta_ij_kl = Kronecker(ii, jj) * Kronecker(kk, ll);
			// dG_dSdS[i*LengthStrs+j] =
			(*d2G_dSdS)(i, j)
			    = 0.5 * ((Kronecker(ii, kk) * Kronecker(jj, ll) - delta_ij_kl / 3.0) / psi1
			             + (MatN1[0] + 12.0 * MatN1[2] * MatN1[2] * I_p2) * delta_ij_kl)
			          / PSI
			      - 0.25 * (DevS[i] / psi1 + (MatN1[0] * In1 + 4.0 * MatN1[2] * MatN1[2] * I_p3) * Kronecker(ii, jj))
			            * (DevS[j] / psi1 + (MatN1[0] * In1 + 4.0 * MatN1[2] * MatN1[2] * I_p3) * Kronecker(kk, ll))
			            / PSI_p3
			      + 2.0 * MatN1[3] * delta_ij_kl;
		}
	}
}

/**************************************************************************
   ROCKFLOW - Funktion: CSolidProperties::dG_dStress_dStress

   Aufgabe:
   Computing the second derivatives of the plastic potential with respect to
   stresses (Weimar's model)

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E :
       double *dG_dSdS            : The second derivatives of the plastic potential
   with respect to stresses
   const long EleIndex        : Element index
   const double *DevS         : The deviatoric stresses
   const double *RotV         : Rotational variables
   const double *S_Invariants : Stress invariants
   const double *MatN1        : Material parameters

   Ergebnis:

   Programmaenderungen:
   08/2003   WW  Erste Version (for 2D)
   02/2004   WW  Extended to 3D
   08/2004   WW  Set as a member of CSolidProperties
 ***************************************************************************/
void CSolidProperties::dG__dStress_dStress(const double* DevS, const double* RotV, const double* S_Invariants,
                                           const double* MatN1, const int LengthStrs)
{
	int i, j;
	int ii, jj;
	const double psi1 = (*data_Plasticity)(14);
	double In1 = S_Invariants[0];
	double I_p2 = In1 * In1;
	double I_p3 = In1 * In1 * In1;
	double PSI = S_Invariants[5];
	double PSI_p3 = S_Invariants[6];

	for (i = 0; i < LengthStrs; i++)
	{
		ii = i;
		jj = i;
		if (i == 3)
		{
			ii = 1;
			jj = 2;
		}
		// For 3D
		else if (i == 4) // s_13
		{
			ii = 1;
			jj = 3;
		}
		else if (i == 5) // s_23
		{
			ii = 2;
			jj = 3;
		}

		for (j = 0; j < LengthStrs; j++)
			(*d2G_dSdS)(i, j)
			    -= 0.5 * (MatN1[0] + 12.0 * MatN1[2] * MatN1[2] * I_p2) * Kronecker(ii, jj) * RotV[j] / PSI
			       - 0.25 * (DevS[i] / psi1 + (MatN1[0] * In1 + 4.0 * MatN1[2] * MatN1[2] * I_p3) * Kronecker(ii, jj))
			             * (MatN1[0] * In1 + 4.0 * MatN1[2] * MatN1[2] * I_p3) * RotV[j] / PSI_p3
			       + 2.0 * MatN1[3] * Kronecker(ii, jj) * RotV[j];
	}
}

/**************************************************************************
   ROCKFLOW - Funktion: CSolidProperties::dG_dSTress_dMat

   Aufgabe:
   Computing the second derivatives of the plastic potential with respect to
   stresses and with recpect to material parameters (Weimar's model)

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E :
       double *dG_dSdM            : The second derivatives of the plastic potential
   with respect to stresses and material parameters

   2D:   dG_dS11dM1 dG_dS11dM2 dG_dS11dM3 dG_dS11dM4
   dG_dS22dM1 dG_dS22dM2 dG_dS22dM3 dG_dS22dM4
   dG_dS12dM1 dG_dS12dM2 dG_dS12dM3 dG_dS12dM4
   dG_dS33dM1 dG_dS33dM2 dG_dS33dM3 dG_dS33dM4

   const long EleIndex        : Element index
   const double *DevS         : The deviatoric stresses
   const double *S_Invariants : Stress invariants
   const double *MatN1        : Material parameters

   Ergebnis:

   Programmaenderungen:
   08/2003   WW  Erste Version (for 2D)
   02/2004   WW  Extended to 3D
   08/2004   WW  Set as a member of CSolidProperties
 ***************************************************************************/
void CSolidProperties::dG_dSTress_dMat(const double* DevS,
                                       const double* S_Invariants,
                                       const double* MatN1,
                                       const int LengthStrs)
{
	int i;
	const double psi1 = (*data_Plasticity)(14);
	const double psi2 = (*data_Plasticity)(15);
	double In1 = S_Invariants[0];
	double I_p2 = In1 * In1;
	double I_p3 = In1 * I_p2;
	double I_p4 = In1 * I_p3;
	double PSI = S_Invariants[5];
	double PSI_p3 = S_Invariants[6];
	int ii, jj;
	double delta_ij;

	for (i = 0; i < LengthStrs; i++)
	{
		ii = i;
		jj = i;
		if (i == 3)
		{
			ii = 1;
			jj = 2;
		}
		// For 3D
		else if (i == 4) // s_13
		{
			ii = 1;
			jj = 3;
		}
		else if (i == 5) // s_23
		{
			ii = 2;
			jj = 3;
		}

		delta_ij = Kronecker(ii, jj);
		// dG_dSdM[i*LengthStrs]   //dG_dS_dAlpha
		(*d2G_dSdM)(i, 0) // dG_dS_dAlpha
		    = 0.5 * In1 * delta_ij / PSI
		      - 0.125 * (DevS[i] / psi1 + (MatN1[0] * In1 + 4.0 * MatN1[2] * MatN1[2] * I_p3) * delta_ij) * I_p2
		            / PSI_p3;
		// dG_dSdM[i*LengthStrs+1]   //dG_dS_dBeta
		(*d2G_dSdM)(i, 1) // dG_dS_dBeta
		    = (1 + psi2) * delta_ij;
		// dG_dSdM[i*LengthStrs+2]   // dG_dS_dDelta
		(*d2G_dSdM)(i, 2) // dG_dS_dDelta
		    = 4.0 * MatN1[2] * I_p3 * delta_ij / PSI
		      - 0.5 * (DevS[i] / psi1 + (MatN1[0] * In1 + 4.0 * MatN1[2] * MatN1[2] * I_p3) * delta_ij) * MatN1[2]
		            * I_p4 / PSI_p3;
		// dG_dSdM[i*LengthStrs+3]   // dG_dS_dEpsilon
		(*d2G_dSdM)(i, 3) // dG_dS_dEpsilon
		    = 2.0 * In1 * delta_ij;
	}
}

/**************************************************************************
   ROCKFLOW - Funktion: CSolidProperties::dfun2

   Aufgabe:
   Computing (Weimar's model)
    1:  the derivatives of the rotational variable equation with
         respect to stresses
    2:  the derivatives of the rotational variable equation with
         respect to rotational variables

   Note : xi_11+xi_22+xi_33=0.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E :
   const long EleIndex        : Element index
   const double *DevS         : The deviatoric stresses
   const double *RotV         : The hardening rotational variables
   const double *S_Invariants : Stress invariants
   const double *MatN1        : Material parameters

   Ergebnis:
   double *d2G_dSdS              : Results 1
   double *d2G_dSdM              : Results 2

   Programmaenderungen:
   08/2003   WW  Erste Version (for 2D)
   02/2004   WW  Extended to 3D
   08/2004   WW  Set as a member of CSolidProperties
 ***************************************************************************/
void CSolidProperties::dfun2(const double* DevS, const double* RotV, const double* S_Invariants, const double* MatN1,
                             const int LengthStrs)
{
	int i, j, l;
	int ii, jj, kk, ll;
	const double mr = (*data_Plasticity)(19);
	const double psi1 = (*data_Plasticity)(14);

	double In1 = S_Invariants[0];
	double I_p3 = In1 * In1 * In1;
	double PSI = S_Invariants[5];
	double PSI_p3 = S_Invariants[6];
	double var = sqrt(S_Invariants[1] / 3.0);
	l = 0;
	for (i = 0; i < LengthStrs; i++)
	{
		// if(i<2) l = i*LengthStrs;
		if (i == 2)
			continue; //
		// if(i>2) l = (i-1)*LengthStrs;

		ii = i;
		jj = i;
		if (i == 3)
		{
			ii = 1;
			jj = 2;
		}
		// For 3D
		else if (i == 4) // s_13
		{
			ii = 1;
			jj = 3;
		}
		else if (i == 5) // s_23
		{
			ii = 2;
			jj = 3;
		}

		for (j = 0; j < LengthStrs; j++)
		{
			kk = j;
			ll = j;
			if (j == 3)
			{
				kk = 1;
				ll = 2;
			}
			// For 3D
			else if (j == 4) // s_13
			{
				kk = 1;
				ll = 3;
			}
			else if (j == 5) // s_23
			{
				kk = 2;
				ll = 3;
			}

			// Derivative with respect to normal stresses
			(*d2G_dSdS)(l, j) = -(
			    (In1 * RotV[i] - mr * DevS[i]) * DevS[j] / (6.0 * var * PSI)
			    - 0.5 * var * (In1 * RotV[i] - mr * DevS[i])
			          * (DevS[j] / psi1 + (MatN1[0] * In1 + 4.0 * MatN1[2] * MatN1[2] * I_p3) * Kronecker(kk, ll))
			          / PSI_p3
			    + var * (RotV[i] * Kronecker(kk, ll)
			             - mr * (Kronecker(ii, kk) * Kronecker(jj, ll) - Kronecker(ii, jj) * Kronecker(kk, ll) / 3.0))
			          / PSI);

			(*d2G_dSdS)(l, j) /= In1;
			(*d2G_dSdS)(l, j) -= 0.5 * (In1 * RotV[i] - mr * DevS[i]) * Kronecker(kk, ll) / (PSI * In1 * In1);

			if (j < 4)
				(*d2G_dSdM)(l, j)
				    = -In1 * (*d2G_dSdS)(i, j) / 3.0 - var * In1 * Kronecker(ii, kk) * Kronecker(jj, ll) / (PSI * In1);

			// Derivative with respect to stresses
			(*d2G_dSdS)(l, j) += -0.5 * var * (In1 * RotV[i] - mr * DevS[i])
			                         * (MatN1[0] * In1 + 4.0 * MatN1[2] * MatN1[2] * I_p3) * RotV[j] / PSI_p3
			                     + var * RotV[i] * RotV[j] / PSI
			                     + 0.5 * (In1 * RotV[i] - mr * DevS[i]) * RotV[j] / (PSI * In1 * In1);
		}
	}
}

/**************************************************************************
   ROCKFLOW - Funktion: CSolidProperties::Gauss_Elimination and Gauss_Back

   Aufgabe: Mini linear equation solver by Gasssian elemination with
           scaled partial pivoting.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E :
       const int DimE             : Dimension of the equations
       double *AA                 : Matrix
   double * rhs               : Right hand side
   int *L                     : Teporarily used array
   double *xx                 : Solution
   Ergebnis:

   Programmaenderungen:
   08/2003   WW  Erste Version
   08/2004   WW  Set as a member of CSolidProperties
 ***************************************************************************/
void CSolidProperties::Gauss_Elimination(const int DimE, Matrix& AA, int* L, double* xx)
{
	int i, j, k, jj, lk;
	double var, R;

	for (i = 0; i < DimE; i++)
	{
		L[i] = i;
		var = 0.0;
		for (j = 0; j < DimE; j++)
		{
			if (fabs(AA(i, j)) > var)
				var = fabs(AA(i, j));
			L[i] = i;
		}
		xx[i] = var;
	}

	for (k = 0; k < DimE - 1; k++)
	{
		var = 0.0;
		jj = 0;
		for (i = k; i < DimE; i++)
		{
			R = fabs(AA(L[i], k) / xx[L[i]]);

			if (R > var)
			{
				jj = i;
				var = R;
			}
		}
		lk = L[jj];
		L[jj] = L[k];
		L[k] = lk;

		for (i = k + 1; i < DimE; i++)
		{
			var = AA(L[i], k) / AA(lk, k);

			for (j = k + 1; j < DimE; j++)
				AA(L[i], j) -= var * AA(lk, j);
			AA(L[i], k) = var;
		}
	}
}

void CSolidProperties::Gauss_Back(const int DimE, Matrix& AA, double* rhs, int* L, double* xx)
{
	int i, j, k;
	double var;

	/* Back substituting */
	for (k = 0; k < DimE - 1; k++)
		for (i = k + 1; i < DimE; i++)
			rhs[L[i]] -= AA(L[i], k) * rhs[L[k]];

	xx[DimE - 1] = rhs[L[DimE - 1]] / AA(L[DimE - 1], DimE - 1);
	for (i = DimE - 2; i >= 0; i--)
	{
		var = rhs[L[i]];
		for (j = i + 1; j < DimE; j++)
			var -= AA(L[i], j) * xx[j];
		xx[i] = var / AA(L[i], i);
	}
}

/*=========================================================================

                    CAM-CLAY Model
   =========================================================================*/
/**************************************************************************
   ROCKFLOW - Funktion: CSolidProperties::CalStress_and_TangentialMatrix_CC

   Aufgabe:
   Computing the stresses at Gauss point and return the plastical status of this
   point. Plastic model: Cam-Clay.
   (Explicit integration of elastic model)

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E :
   double *TryStress                   :   Incremental straines as input
   New stresses as output
   double *Dep                         :   The consist tangential matrix
   const int GPiGPj                    :   Local indeces of Gauss Points
   const int Update                    :   Indicator. If true, update the gauss point
   values

   Ergebnis:
   - double* - Deviatoric effetive stresses, s11,s22,s12,s33

   Programmaenderungen:
   11/2003   WW  Erste Version

**************************************************************************/
//#define New
#define associative
void CSolidProperties::CalStress_and_TangentialMatrix_CC(const int GPiGPj, const ElementValue_DM* ele_val,
                                                         double* dStrain, Matrix* Dep, const int Update)
{
	int i, ns;
	double p, q, p_tr, q_tr, p_c, p_cn;
	double ep;
	double e_0, dev;
	double F0, F, vep, vep0, TolP;
	double var1, alpha1, alpha2, beta1, beta2;
	double gamma1, gamma2, gamma3, gamma4, gamma5;
	double dfdp, dfdq, dampFac;
	const double fac = sqrt(2.0 / 3.0);
#ifdef New
	double J11, J12, J21, J22;
#else
	double alpha3, Jac;
#endif
	static double DevStress[6], TryStress[6];

	const int MaxI = 400;
	int NPStep;

	const double SwellIndex = (*data_Plasticity)(2);
	const double CompressIndex = (*data_Plasticity)(1);
	const double M_s = (*data_Plasticity)(0);
	const double M2 = M_s * M_s;
	const double pmin = (*data_Plasticity)(9);

	double vartheta = 0.0;

	TolP = Tolerance_Local_Newton * 1.0e4;

	int dim = 2;
	ns = 4;
	if (Dep->Cols() > 4)
	{
		dim = 3;
		ns = 6;
	}
	bool isLoop = true; // Used only to avoid warnings with .net

	// Get the total effective plastic strain
	ep = (*ele_val->pStrain)(GPiGPj);

	// Get the preconsolidation pressure of the previous time step
	p_cn = (*ele_val->prep0)(GPiGPj);

	// Get the void ratio of the previous time step
	e_0 = (*ele_val->e_i)(GPiGPj);

	p = -((*ele_val->Stress)(0, GPiGPj) + (*ele_val->Stress)(1, GPiGPj) + (*ele_val->Stress)(2, GPiGPj)) / 3.0;

	if (fabs(p) < pmin)
		p = pmin;

	// Volume strain increment
	dev = -(dStrain[0] + dStrain[1] + dStrain[2]);
	// TEST. Sign?
	vartheta = (1.0 + e_0) / (CompressIndex - SwellIndex);

	// if(fabs(dev)<MKleinsteZahl)
	if (SwellingPressureType == 3)
	{
		K = (1.0 + e_0) * fabs(p) / (*data_Youngs)(6);
		if (K < 10.e6)
			K = 10.e6;
	}
	else
		K = (1.0 + e_0) * fabs(p) / SwellIndex;
	// else
	//    K = (1.0+e_0)*fabs(p)*(exp(PoissonRatio*dev/SwellIndex)-1.0)/dev;

	G = 1.5 * K * (1 - 2.0 * PoissonRatio) / (1 + PoissonRatio);

	Lambda = K - 2.0 * G / 3.0;

	ElasticConsitutive(dim, Dep);

	for (i = 0; i < ns; i++)
		TryStress[i] = 0.0;
	Dep->multi(dStrain, TryStress);

	for (i = 0; i < ns; i++)
	{
		TryStress[i] += (*ele_val->Stress)(i, GPiGPj);
		DevStress[i] = TryStress[i];
	}

	p_tr = -DeviatoricStress(DevStress) / 3.0;

	q_tr = sqrt(1.5 * TensorMutiplication2(DevStress, DevStress, dim));

	// If yield, integrate the stress
	F = q_tr * q_tr / M2 + p_tr * (p_tr - p_cn);

	p = p_tr;
	q = q_tr;
	p_c = p_cn;
	vep = 0.0;
	vep0 = 0.0;

	dampFac = 1.0;
	NPStep = 0;

	F0 = F;
	if (pcs_deformation == 1)
		F = -1.0;
	if ((*data_Plasticity)(3) < MKleinsteZahl) // p_c0=0
		F = -1.0;
	// TEST CAM-CLAY
	if (p_tr < 0)
		F = -1.0;
	if (F > 0.0 && !PreLoad) // in yield status
	{
		// Local Newton-Raphson procedure to compute the volume plastic strain
		vep = 0.0;

#ifdef associative
#ifdef New
		// Associative flow rule
		while (isLoop) // Newton step for the plastic multiplier
		{
			NPStep++;
			if (NPStep > MaxI)
				// printf("\n Too much iteration in Newton step in the integration of Cam-Clay \n");
				break;
			// abort();
			dfdp = 2.0 * p - p_c;
			dfdq = 2.0 * q / M2;
			gamma1 = 1.0 + 2.0 * vep * K;
			// J11: df/ddphi. J12: df/dpc
			// dp/dphi
			alpha1 = -K * dfdp / gamma1;
			// dq/dphi
			alpha2 = -q / (vep + M2 / (G * 6.0));
			//
			J11 = alpha1 * dfdp + dfdq * alpha2;
			J12 = dfdp * vep * K / gamma1;
			// J21: dg/ddphi. J22: dg/dpc
			gamma2 = p_cn * exp(vartheta * dev * dfdp);
			J21 = gamma2 * vartheta * (dfdp + 2.0 * vep * alpha1);
			J22 = gamma2 * vartheta * vep * (2.0 * dev * K / gamma1 - 1.0) - 1.0;

			F = q * q / M2 + p * (p - p_c);
			beta1 = gamma2 - p_c; // g;
			if (fabs(beta1) < TolP)
				//         if(fabs(F)<TolF) break;
				if (fabs(F / F0) < TolP)
					break;
			beta2 = J11 * J22 - J21 * J12;

			p_c -= (beta1 * J11 - F * J21) / beta2;
			vep -= (F * J22 - beta1 * J12) / beta2;

			p = (p_tr + vep * K * p_c) / gamma1;
			K = q = q_tr / (1.0 + 6.0 * G * vep / M2);
		}
#else // ifdef New
		// Associative flow rule
		while (isLoop) // Newton step for the plastic multiplier
		{
			NPStep++;

			alpha1 = (2.0 * p - p_c) / (1.0 + (2.0 * K + vartheta * p_c) * vep);
			alpha2 = -q / (vep + M2 / (G * 6.0));
			alpha3 = vartheta * p_c * alpha1;
			alpha1 *= -K;

			dfdp = 2.0 * p - p_c;
			dfdq = 2.0 * q / M2;

			Jac = alpha1 * dfdp + dfdq * alpha2 - p * alpha3;

			while (isLoop) // Damp
			{
				vep = vep0 - dampFac * F / Jac;
				p_c = 0.0;

				dampFac = 1.0;
				while (isLoop) // Newton step for p_c
				{
					NPStep++;
					if (NPStep > MaxI)
						//    printf("\n Too much iteration in Newton step the integration of Cam-Clay \n");
						// TEST	abort();
						// test
						break;

					//
					alpha1 = vartheta * vep / (1.0 + 2.0 * vep * K);
					alpha2 = p_cn * exp(alpha1 * (2.0 * p_tr - p_c));
					beta1 = -alpha1 * alpha2 - 1.0; // dG(p_c)
					alpha2 -= p_c; // G(p_c)
					p_c -= alpha2 / beta1;
					if (fabs(alpha2) < TolP)
						break;
					if (p_c < 0.0)
						break;
				}

				p = (p_tr + vep * K * p_c) / (1.0 + 2.0 * vep * K);
				q = q_tr / (1.0 + 6.0 * G * vep / M2);

				F = q * q / M2 + p * (p - p_c);
				if (F > 0.0)
					break;
				if (fabs(F / F0) < TolP)
					break;
				dampFac = 0.8;
			}
			vep0 = vep;
			if (fabs(F / F0) < TolP)
				break;
		}
#endif
		// Plastic strain
		//    alpha1 = 6.0*q*q/(M2*M2*(2.0*p-p_c)*(2.0*p-p_c));
		//    ep += fabs(vep)*sqrt(2.0*(1.0/9.0+alpha1)/3.0);
		ep += 3.0 * fabs(vep) * q / M2;

#else // ifdef associative
		double var2 = 0.0;
		while (1)
		{
			NPStep++;
			if (NPStep > MaxI)
			{
				printf("\n Too much iteration in Newton step in the integration of Cam-Clay \n");
				abort();
			}

			var1 = (p_tr - p) / (p_tr - 0.5 * p_c);

			dfdp = 2.0 * p - p_c;
			dfdq = 2.0 * q / M2;
			Jac = -K * dfdp // dF/dp *  dp/dv
			      - dfdq * q_tr * (K + 0.5 * p_c * vartheta * var1) / (p_tr - 0.5 * p_c) // dF/dq  * dq/dv
			      - p * vartheta * p_c; // df/dp_c  * dp_c/dv

			// Update
			dampFac = 1.0;

			while (1)
			{
				vep = vep0 - dampFac * F / Jac;

				p_c = p_cn * exp(vep * vartheta);
				p = p_tr - K * vep;
				q = q_tr * (p - 0.5 * p_c) / (p_tr - 0.5 * p_c);

				F = q * q / M2 + p * (p - p_c);
				if (F > 0.0)
					break;
				// if(fabs(F)/RF0<TolLocalNewT*1.0e-3) break;
				if (fabs(F) < Tolerance_Local_Newton)
					break;
				dampFac = 0.5;
			}
			vep0 = vep;
			// if(fabs(F)/RF0<TolLocalNewT*1.0e-3) break;
			if (fabs(F) < Tolerance_Local_Newton)
				break;
		}

		ep += 3.0 * fabs(vep) * q / M2;
#endif

		//-------------------------------------------------------------
		// Consistent tangential matrix

		if (Update < 1)
		{
#ifdef associative
			// a
			alpha1 = 1.0 + 2.0 * K * vep + p_c * vartheta * vep;
			// a1
			double a1 = (1.0 + p_c * vartheta * vep) / alpha1;
			double a2 = -(2.0 * p - p_c) / alpha1; // a2
			// a3
			double a3 = 2.0 * p_c * vartheta * vep / alpha1;
			// a4
			double a4 = vartheta * p_c * (2.0 * p - p_c) / (K * alpha1);
			// a5
			double a5 = sqrt(1.5) / (1.0 + 6.0 * G * vep / M2);
			// a6
			double a6 = -3.0 * q / (1.0 + 6.0 * G * vep / M2) / M2;

			alpha1 = -4.0 * G * q * a6 / M2 - K * ((2.0 * a2 - a4) * p - a2 * p_c);
			double b1 = -K * ((a3 - 2.0 * a1) * p + a1 * p_c) / alpha1;
			double b2 = 4.0 * G * a5 * q / (alpha1 * M2);

			gamma1 = 2.0 * G * q / q_tr;
			gamma2 = K * (a1 + a2 * b1) - gamma1 / 3.0;
			gamma3 = -K * a2 * b2;
			gamma4 = -2.0 * fac * G * a6 * b1;
			gamma5 = 2.0 * G * fac * (a5 + a6 * b2) - gamma1;
#else // ifdef associative
			dfdp = 2.0 * p - p_c;
			dfdq = 2.0 * q / M2;
			var1 = (p_tr - p) / (p_tr - 0.5 * p_c);
			var2 = (p - 0.5 * p_c) / (p_tr - 0.5 * p_c);

			alpha1 = q_tr * var1 / (p_tr - 0.5 * p_c);
			alpha2 = sqrt(3.0 / 2.0) * var2;
			alpha3 = -q_tr * (1.0 + 0.5 * vartheta * p_c * var1 / K) / (p_tr - 0.5 * p_c);

			double beta0 = K * dfdp - K * dfdq * alpha3 + vartheta * p * p_c;
			beta1 = K * (dfdq * alpha1 + dfdp) / beta0;
			beta2 = 2.0 * G * dfdq * alpha2 / beta0;

			gamma1 = 2.0 * G * q / q_tr;
			gamma2 = K * (1.0 - beta1) - gamma1 / 3.0;
			gamma3 = K * beta2;
			gamma4 = -fac * K * (alpha1 + alpha3 * beta1);
			gamma5 = 2.0 * G * fac * alpha2 - gamma1 + fac * alpha3 * beta2 * K;
#endif

			// Normalize the stress
			for (i = 0; i < ns; i++)
			{
				TryStress[i] = DevStress[i];
				TryStress[i] /= q_tr * fac;
			}

			// Row 1
			beta1 = gamma2 + gamma4 * TryStress[0];
			beta2 = gamma3 + gamma5 * TryStress[0];
			(*Dep)(0, 0) = gamma1 + beta1 + beta2 * TryStress[0];
			(*Dep)(0, 1) = beta1 + beta2 * TryStress[1];
			(*Dep)(0, 2) = beta1 + beta2 * TryStress[2];
			(*Dep)(0, 3) = beta2 * TryStress[3];
			if (dim == 3)
			{
				(*Dep)(0, 4) = beta2 * TryStress[4];
				(*Dep)(0, 5) = beta2 * TryStress[5];
			}

			// Row 2
			beta1 = gamma2 + gamma4 * TryStress[1];
			beta2 = gamma3 + gamma5 * TryStress[1];
			(*Dep)(1, 0) = beta1 + beta2 * TryStress[0];
			(*Dep)(1, 1) = gamma1 + beta1 + beta2 * TryStress[1];
			(*Dep)(1, 2) = beta1 + beta2 * TryStress[2];
			(*Dep)(1, 3) = beta2 * TryStress[3];
			if (dim == 3)
			{
				(*Dep)(1, 4) = beta2 * TryStress[4];
				(*Dep)(1, 5) = beta2 * TryStress[5];
			}

			// Row 3
			beta1 = gamma2 + gamma4 * TryStress[2];
			beta2 = gamma3 + gamma5 * TryStress[2];
			(*Dep)(2, 0) = beta1 + beta2 * TryStress[0];
			(*Dep)(2, 1) = beta1 + beta2 * TryStress[1];
			(*Dep)(2, 2) = gamma1 + beta1 + beta2 * TryStress[2];
			(*Dep)(2, 3) = beta2 * TryStress[3];
			if (dim == 3)
			{
				(*Dep)(2, 4) = beta2 * TryStress[4];
				(*Dep)(2, 5) = beta2 * TryStress[5];
			}

			// Row 4
			beta1 = gamma4 * TryStress[3];
			beta2 = gamma5 * TryStress[3];
			(*Dep)(3, 0) = beta1 + beta2 * TryStress[0];
			(*Dep)(3, 1) = beta1 + beta2 * TryStress[1];
			(*Dep)(3, 2) = beta1 + beta2 * TryStress[2];
			(*Dep)(3, 3) = gamma1 + beta2 * TryStress[3];
			if (dim == 3)
			{
				(*Dep)(3, 4) = beta2 * TryStress[4];
				(*Dep)(3, 5) = beta2 * TryStress[5];
				// End row 4

				// Row 5
				beta1 = gamma4 * TryStress[4];
				beta2 = gamma5 * TryStress[4];
				(*Dep)(4, 0) = beta1 + beta2 * TryStress[0];
				(*Dep)(4, 1) = beta1 + beta2 * TryStress[1];
				(*Dep)(4, 2) = beta1 + beta2 * TryStress[2];
				(*Dep)(4, 3) = beta2 * TryStress[3];
				(*Dep)(4, 4) = gamma1 + beta2 * TryStress[4];
				(*Dep)(4, 5) = beta2 * TryStress[5];

				// Row 6
				beta1 = gamma4 * TryStress[5];
				beta2 = gamma5 * TryStress[5];
				(*Dep)(5, 0) = beta1 + beta2 * TryStress[0];
				(*Dep)(5, 1) = beta1 + beta2 * TryStress[1];
				(*Dep)(5, 2) = beta1 + beta2 * TryStress[2];
				(*Dep)(5, 3) = beta2 * TryStress[3];
				(*Dep)(5, 4) = beta2 * TryStress[4];
				(*Dep)(5, 5) = gamma1 + beta2 * TryStress[5];
			}
		}

//-------------------------------------------------------------

// Update stresses
#ifdef associative
		for (i = 0; i < ns; i++)
		{
			DevStress[i] /= 1.0 + 6.0 * G * vep / M2;
			TryStress[i] = DevStress[i];
		}
#else
		for (i = 0; i < ns; i++)
		{
			DevStress[i] /= 1.0 + 2.0 * K * vep / (2.0 * p + p_c);
			TryStress[i] = DevStress[i];
		}
#endif

		// True stress
		for (i = 0; i < 3; i++)
			TryStress[i] -= p;
	}
	else if (Update < 1)
		ElasticConsitutive(dim, Dep);

	for (i = 0; i < ns; i++)
		dStrain[i] = TryStress[i];

	// Save the current stresses
	if (Update > 0)
	{
		if ((*data_Plasticity)(3) < MKleinsteZahl) // p_c0=0
		{
			p_c = p + q * q / (M2 * p);
			(*ele_val->prep0)(GPiGPj) = p_c;
		}
		if (ep > 0.0)
		{
			var1 = dev * (1. + e_0);
			e_0 -= var1;
			(*ele_val->pStrain)(GPiGPj) = ep;
			(*ele_val->e_i)(GPiGPj) = e_0;
			(*ele_val->prep0)(GPiGPj) = p_c;
		}
	}

	// For the case of the initial stress being given, the contribution of
	// the initial stress to the right hand side should be taken into account
	// If the initial stress is not accounted into the final stress
	//
	//  for(i=0; i<3; i++)
	//    dStrain[i] -= (*data_Plasticity)(6+i); // Initial stress
	//
}

/**************************************************************************
   GeoSys - Function: Integration with substep for CAM-Clay like model
   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E :
     double *TryStress                   :   Incremental straines as input
                                             New stresses as output
     double *Dep                         :   The consist tangential matrix
     const int GPiGPj                    :   Local indeces of Gauss Points
     const int Update                    :   Indicator. If true, update the gauss point
                                             values

   Ergebnis:
   - double* - Deviatoric effetive stresses, s11,s22,s12,s33

   Programmaenderungen:
   06/2008   WW  Programming
**************************************************************************/
void CSolidProperties::CalStress_and_TangentialMatrix_CC_SubStep(const int GPiGPj, const ElementValue_DM* ele_val,
                                                                 double* dStrain, Matrix* Dep, const int Update)
{
	int i, ns;
	double p, q, p_tr, q_tr, p_c, p_cn;
	double ep;
	double e_0, dev;
	double F0, F, vep, vep0;
	double var1, alpha1, alpha2, beta1, beta2;
	double gamma1, gamma2, gamma3, gamma4, gamma5;
	double dfdp, dfdq, dampFac;
	const double fac = sqrt(2.0 / 3.0);
	double alpha3, Jac;
	static double DevStress[6], TryStress[6];
	static double dStress0[6], dStress1[6];

	const int MaxI = 400;
	int NPStep;
	double sub_step = 0.5;
	double sub_step_sum = 0.0;
	//
	const double SwellIndex = (*data_Plasticity)(2);
	const double CompressIndex = (*data_Plasticity)(1);
	const double M_s = (*data_Plasticity)(0);
	const double M2 = M_s * M_s;
	const double pmin = (*data_Plasticity)(9);
	//
	double vartheta = 0.0;
	double suc = 0.0;
	double dsuc = 0.0;

	p = q = q_tr = p_c = vep = 0.;
	NPStep = 1;
	//
	int dim = 2;
	ns = 4;
	if (Dep->Cols() > 4)
	{
		dim = 3;
		ns = 6;
	}
	bool isLoop = true; // Used only to avoid warnings with .net

	// Get the total effective plastic strain
	ep = (*ele_val->pStrain)(GPiGPj);

	// Get the preconsolidation pressure of the previous time step
	p_cn = (*ele_val->prep0)(GPiGPj);

	// Get the void ratio of the previous time step
	e_0 = (*ele_val->e_i)(GPiGPj);
	// Stress of previous load step
	dev = -(dStrain[0] + dStrain[1] + dStrain[2]);
	for (i = 0; i < ns; i++)
	{
		TryStress[i] = (*ele_val->Stress)(i, GPiGPj);
		dStress0[i] = dStress1[i] = 0.;
	}
	// TEST
	int test_c = 0;
	// Begin substeps
	// WW bool OK = true; //OK411
	for (;;) // WW while(OK)
	{
		/*
		   if(test_c>2)
		   test_c = test_c;
		 */
		// TEST
		test_c++;
		// if(test_c>20)
		//  test_c = test_c;
		if (test_c > 100)
			break;

		for (int kk = 0; kk < 2; kk++)
		{
			double* dsig = NULL;
			double de_vsw = 0.;
			if (kk == 0)
			{
				dsig = dStress0;
				for (i = 0; i < ns; i++)
					dsig[i] = TryStress[i];
			}
			else
			{
				dsig = dStress1;
				for (i = 0; i < ns; i++)
					dsig[i] = TryStress[i] + dStress0[0];
			}
			p = -(dsig[0] + dsig[1] + dsig[2]) / 3.0;
			if (fabs(p) < pmin)
				p = pmin;
			// Volume strain increment
			vartheta = (1.0 + e_0) / (CompressIndex - SwellIndex);
			// if(fabs(dev)<MKleinsteZahl)
			if (SwellingPressureType == 3)
			{
				suc = (*data_Youngs)(7);
				dsuc = (*data_Youngs)(8);
				// 0: sign is nagive -
				de_vsw = TEPSwellingParameter(p) * dsuc / (suc + 1.0e5);
				for (i = 0; i < 3; i++)
					dStrain[i] += de_vsw;
				K = (1.0 + e_0) * fabs(p) / (*data_Youngs)(6);
				if (K < 10.e6)
					K = 10.e6;
				vartheta = (1.0 + e_0) / (CompressIndex - (*data_Youngs)(6));
			}
			else
				K = (1.0 + e_0) * fabs(p) / SwellIndex;
			// else
			//    K = (1.0+e_0)*fabs(p)*(exp(PoissonRatio*dev/SwellIndex)-1.0)/dev;

			G = 1.5 * K * (1 - 2.0 * PoissonRatio) / (1 + PoissonRatio);
			Lambda = K - 2.0 * G / 3.0;
			ElasticConsitutive(dim, Dep);
			//
			for (i = 0; i < ns; i++)
				dsig[i] = 0.0;
			Dep->multi(dStrain, dsig, sub_step);
			// Recover strain increment
			if (SwellingPressureType == 3)
				for (i = 0; i < 3; i++)
					dStrain[i] -= de_vsw;
		}
		// Stress estimation
		// double norm_ds = 0.;
		// double norm_s = 0.;
		for (i = 0; i < ns; i++)
		{
			TryStress[i] += 0.5 * (dStress0[i] + dStress1[i]);
			// DevStress as temporary buffer
			DevStress[i] = dStress0[i] - dStress1[i];
		}
		double norms = StressNorm(TryStress, dim);
		double q_f = 1.0;
		double R_n = 1.0;
		if (norms > DBL_EPSILON)
		{
			R_n = 0.5 * StressNorm(DevStress, dim) / norms;
			if (R_n > DBL_EPSILON)
				q_f = 0.95 * sqrt(s_tol / R_n);
		}
		else
			R_n = 0.0;
		if (q_f < 0.2)
			q_f = 0.2;
		if (q_f > 1.2)
			q_f = 1.2;
		sub_step *= q_f;
		if ((sub_step_sum + sub_step) > 1.0)
			sub_step = 1.0 - sub_step_sum;
		// Ckeck convergence
		if (R_n < s_tol) // accepted
			sub_step_sum += sub_step;
		else
		{
			for (i = 0; i < ns; i++)
				TryStress[i] -= 0.5 * (dStress0[i] + dStress1[i]);
			continue;
		}
		//

		for (i = 0; i < ns; i++)
			DevStress[i] = TryStress[i];

		p_tr = -DeviatoricStress(DevStress) / 3.0;

		q_tr = sqrt(1.5 * TensorMutiplication2(DevStress, DevStress, dim));

		// If yield, integrate the stress
		F = q_tr * q_tr / M2 + p_tr * (p_tr - p_cn);

		p = p_tr;
		q = q_tr;
		p_c = p_cn;
		vep = 0.0;
		vep0 = 0.0;

		dampFac = 1.0;
		NPStep = 0;

		F0 = F;
		if (pcs_deformation == 1)
			F = -1.0;
		if ((*data_Plasticity)(3) < MKleinsteZahl) // p_c0=0
			F = -1.0;
		// TEST CAM-CLAY
		if (p_tr < 0)
			F = -1.0;

		if (F > f_tol && !PreLoad) // in yield status
		{
			// Local Newton-Raphson procedure to compute the volume plastic strain
			vep = 0.0;

			// Associative flow rule
			while (isLoop) // Newton step for the plastic multiplier
			{
				NPStep++;

				alpha1 = (2.0 * p - p_c) / (1.0 + (2.0 * K + vartheta * p_c) * vep);
				alpha2 = -q / (vep + M2 / (G * 6.0));
				alpha3 = vartheta * p_c * alpha1;
				alpha1 *= -K;

				dfdp = 2.0 * p - p_c;
				dfdq = 2.0 * q / M2;

				Jac = alpha1 * dfdp + dfdq * alpha2 - p * alpha3;

				while (isLoop) // Damp
				{
					vep = vep0 - dampFac * F / Jac;
					p_c = 0.0;

					dampFac = 1.0;
					while (isLoop) // Newton step for p_c
					{
						NPStep++;
						if (NPStep > MaxI)
							//    printf("\n Too much iteration in Newton step the integration of Cam-Clay \n");
							// TEST	abort();
							// test
							break;

						//
						alpha1 = vartheta * vep / (1.0 + 2.0 * vep * K);
						alpha2 = p_cn * exp(alpha1 * (2.0 * p_tr - p_c));
						beta1 = -alpha1 * alpha2 - 1.0; // dG(p_c)
						alpha2 -= p_c; // G(p_c)
						p_c -= alpha2 / beta1;
						if (fabs(alpha2) < s_tol)
							break;
						if (p_c < 0.0)
							break;
					}

					p = (p_tr + vep * K * p_c) / (1.0 + 2.0 * vep * K);
					q = q_tr / (1.0 + 6.0 * G * vep / M2);

					F = q * q / M2 + p * (p - p_c);
					if (F > 0.0)
						break;
					if (fabs(F / F0) < s_tol)
						break;
					dampFac = 0.8;
				}
				vep0 = vep;
				if (fabs(F / F0) < s_tol)
					break;
			}
		}
		// End substep
		// Plastic strain
		//    alpha1 = 6.0*q*q/(M2*M2*(2.0*p-p_c)*(2.0*p-p_c));
		//    ep += fabs(vep)*sqrt(2.0*(1.0/9.0+alpha1)/3.0);
		ep += 3.0 * fabs(vep) * q / M2;
		if (fabs(sub_step_sum - 1.0) < DBL_EPSILON)
			break;
	}
	//-------------------------------------------------------------
	// Consistent tangential matrix

	if (Update < 1 && NPStep > 0)
	{
		alpha1 = 1.0 + 2.0 * K * vep + p_c * vartheta * vep; // a
		// a1
		double a1 = (1.0 + p_c * vartheta * vep) / alpha1;
		double a2 = -(2.0 * p - p_c) / alpha1; // a2
		double a3 = 2.0 * p_c * vartheta * vep / alpha1; // a3
		// a4
		double a4 = vartheta * p_c * (2.0 * p - p_c) / (K * alpha1);
		double a5 = sqrt(1.5) / (1.0 + 6.0 * G * vep / M2); // a5
		// a6
		double a6 = -3.0 * q / (1.0 + 6.0 * G * vep / M2) / M2;

		alpha1 = -4.0 * G * q * a6 / M2 - K * ((2.0 * a2 - a4) * p - a2 * p_c);
		double b1 = -K * ((a3 - 2.0 * a1) * p + a1 * p_c) / alpha1;
		double b2 = 4.0 * G * a5 * q / (alpha1 * M2);

		gamma1 = 2.0 * G * q / q_tr;
		gamma2 = K * (a1 + a2 * b1) - gamma1 / 3.0;
		gamma3 = -K * a2 * b2;
		gamma4 = -2.0 * fac * G * a6 * b1;
		gamma5 = 2.0 * G * fac * (a5 + a6 * b2) - gamma1;
		//
		// Normalize the stress
		for (i = 0; i < ns; i++)
		{
			TryStress[i] = DevStress[i];
			TryStress[i] /= q_tr * fac;
		}

		// Row 1
		beta1 = gamma2 + gamma4 * TryStress[0];
		beta2 = gamma3 + gamma5 * TryStress[0];
		(*Dep)(0, 0) = gamma1 + beta1 + beta2 * TryStress[0];
		(*Dep)(0, 1) = beta1 + beta2 * TryStress[1];
		(*Dep)(0, 2) = beta1 + beta2 * TryStress[2];
		(*Dep)(0, 3) = beta2 * TryStress[3];
		if (dim == 3)
		{
			(*Dep)(0, 4) = beta2 * TryStress[4];
			(*Dep)(0, 5) = beta2 * TryStress[5];
		}

		// Row 2
		beta1 = gamma2 + gamma4 * TryStress[1];
		beta2 = gamma3 + gamma5 * TryStress[1];
		(*Dep)(1, 0) = beta1 + beta2 * TryStress[0];
		(*Dep)(1, 1) = gamma1 + beta1 + beta2 * TryStress[1];
		(*Dep)(1, 2) = beta1 + beta2 * TryStress[2];
		(*Dep)(1, 3) = beta2 * TryStress[3];
		if (dim == 3)
		{
			(*Dep)(1, 4) = beta2 * TryStress[4];
			(*Dep)(1, 5) = beta2 * TryStress[5];
		}

		// Row 3
		beta1 = gamma2 + gamma4 * TryStress[2];
		beta2 = gamma3 + gamma5 * TryStress[2];
		(*Dep)(2, 0) = beta1 + beta2 * TryStress[0];
		(*Dep)(2, 1) = beta1 + beta2 * TryStress[1];
		(*Dep)(2, 2) = gamma1 + beta1 + beta2 * TryStress[2];
		(*Dep)(2, 3) = beta2 * TryStress[3];
		if (dim == 3)
		{
			(*Dep)(2, 4) = beta2 * TryStress[4];
			(*Dep)(2, 5) = beta2 * TryStress[5];
		}

		// Row 4
		beta1 = gamma4 * TryStress[3];
		beta2 = gamma5 * TryStress[3];
		(*Dep)(3, 0) = beta1 + beta2 * TryStress[0];
		(*Dep)(3, 1) = beta1 + beta2 * TryStress[1];
		(*Dep)(3, 2) = beta1 + beta2 * TryStress[2];
		(*Dep)(3, 3) = gamma1 + beta2 * TryStress[3];
		if (dim == 3)
		{
			(*Dep)(3, 4) = beta2 * TryStress[4];
			(*Dep)(3, 5) = beta2 * TryStress[5];
			// End row 4

			// Row 5
			beta1 = gamma4 * TryStress[4];
			beta2 = gamma5 * TryStress[4];
			(*Dep)(4, 0) = beta1 + beta2 * TryStress[0];
			(*Dep)(4, 1) = beta1 + beta2 * TryStress[1];
			(*Dep)(4, 2) = beta1 + beta2 * TryStress[2];
			(*Dep)(4, 3) = beta2 * TryStress[3];
			(*Dep)(4, 4) = gamma1 + beta2 * TryStress[4];
			(*Dep)(4, 5) = beta2 * TryStress[5];

			// Row 6
			beta1 = gamma4 * TryStress[5];
			beta2 = gamma5 * TryStress[5];
			(*Dep)(5, 0) = beta1 + beta2 * TryStress[0];
			(*Dep)(5, 1) = beta1 + beta2 * TryStress[1];
			(*Dep)(5, 2) = beta1 + beta2 * TryStress[2];
			(*Dep)(5, 3) = beta2 * TryStress[3];
			(*Dep)(5, 4) = beta2 * TryStress[4];
			(*Dep)(5, 5) = gamma1 + beta2 * TryStress[5];
		}
		//-------------------------------------------------------------

		// Update stresses
		for (i = 0; i < ns; i++)
		{
			DevStress[i] /= 1.0 + 6.0 * G * vep / M2;
			TryStress[i] = DevStress[i];
		}

		// True stress
		for (i = 0; i < 3; i++)
			TryStress[i] -= p;
	}
	else if (Update < 1)
		ElasticConsitutive(dim, Dep);

	for (i = 0; i < ns; i++)
		dStrain[i] = TryStress[i];

	// Save the current stresses
	if (Update > 0)
	{
		if ((*data_Plasticity)(3) < MKleinsteZahl) // p_c0=0
		{
			p_c = p + q * q / (M2 * p);
			(*ele_val->prep0)(GPiGPj) = p_c;
		}
		if (ep > 0.0)
		{
			var1 = dev * (1. + e_0);
			e_0 -= var1;
			(*ele_val->pStrain)(GPiGPj) = ep;
			(*ele_val->e_i)(GPiGPj) = e_0;
			(*ele_val->prep0)(GPiGPj) = p_c;
		}
	}

	// For the case of the initial stress being given, the contribution of
	// the initial stress to the right hand side should be taken into account
	// If the initial stress is not accounted into the final stress
	//
	//  for(i=0; i<3; i++)
	//    dStrain[i] -= (*data_Plasticity)(6+i); // Initial stress
	//
}
/**************************************************************************
   FEMLib-Method:
   Task: Caculate increment of strain deduced by creep
   Programing:
   12/2005 WW
   last modified:
**************************************************************************/
void CSolidProperties::AddStain_by_Creep(const int ns, double* stress_n, double* dstrain, double temperature)
{
	int i, dim;
	double norn_S, fac = 0.0;
	DeviatoricStress(stress_n);
	dim = 2;
	if (ns > 4)
		dim = 3;
	norn_S = sqrt(3.0 * TensorMutiplication2(stress_n, stress_n, dim) / 2.0);
	if (norn_S < DBL_MIN)
		return;
	switch (Creep_mode)
	{
		case 1:
			//  fac = pow(1.5, (*data_Creep)(1)+1.0)*(*data_Creep)(0)*pow(norn_S, (*data_Creep)(1)-1.0)*dt;
			// fac = pow(2.0/3.0, (*data_Creep)(1))*(*data_Creep)(0)*pow(norn_S, (*data_Creep)(1))*dt;
			fac = (*data_Creep)(0) * pow(norn_S, (*data_Creep)(1)) * dt;
			break;
		case 2:
			// gas constant = R = 8.314472(15) J ?K-1 ?mol-1
			// ec= A*exp(-G/RT)s^n
			fac = 1.5 * dt * (*data_Creep)(0) * exp(-(*data_Creep)(2) / (8.314472 * (temperature + 273.15)))
			      * pow(norn_S, (*data_Creep)(1));
			break;
		// TN: BGRb
		case 3:
			fac = 1.5 * dt * ((*data_Creep)(0) * exp(-(*data_Creep)(2) / (8.314472 * (temperature + 273.15)))
			                      * pow(norn_S, (*data_Creep)(1))
			                  + (*data_Creep)(3) * exp(-(*data_Creep)(5) / (8.314472 * (temperature + 273.15)))
			                        * pow(norn_S, (*data_Creep)(4)));
			break;
		// TN: BGRsf
		case 4:
			fac = 1.5 * dt * ((*data_Creep)(0) * exp(-(*data_Creep)(2) / (8.314472 * (temperature + 273.15)))
			                      * pow(norn_S, (*data_Creep)(1))
			                  + (*data_Creep)(4) * pow(norn_S, 2));
			break;
	}
	for (i = 0; i < ns; i++)
		dstrain[i] -= fac * stress_n[i] / norn_S;
}
/**************************************************************************
   FEMLib-Method:
   Task: Caculate increment of strain induced by HL creep model
   Programing:
   10/2008 UJG/WW
   last modified:
**************************************************************************/
void CSolidProperties::CleanTrBuffer_HL_ODS()
{
	for (int i = 0; i < 6; i++)
		(*data_Creep)(i, 1) = 0.0;
}
/**************************************************************************
   FEMLib-Method:
   Task: Caculate increment of strain induced by HL creep model
   Programing:
   10/2008 UJG/WW
   last modified:
**************************************************************************/
void CSolidProperties::AccumulateEtr_HL_ODS(const ElementValue_DM* ele_val, const int nGS)
{
	int i, ns;
	ns = ele_val->xi->Size();
	//
	for (i = 0; i < ns; i++)
		(*ele_val->xi)(i) += (*data_Creep)(i, 1) / (double)nGS;
}
/**************************************************************************
   FEMLib-Method:
   Task: Caculate increment of strain induced by HL creep model
   Programing:
   01/2009 UJG/WW
   last modified:
**************************************************************************/
void CSolidProperties::CalcYoungs_SVV(const double strain_v)
{
	double nv = Poisson_Ratio();
	E = (*data_Youngs)(0) / (1. + (*data_Youngs)(1) * pow(strain_v, (*data_Youngs)(2)));
	Lambda = E * nv / ((1. + nv) * (1. - 2. * nv));
	G = 0.5 * E / (1. + nv);
	K = (3.0 * Lambda + 2.0 * G) / 3.0;
}
/**************************************************************************
   FEMLib-Method:
   Task: Caculate increment of strain induced by HL creep model
   Programing:
   10/2008 UJG/WW
   last modified:
**************************************************************************/
void CSolidProperties::AddStain_by_HL_ODS(const ElementValue_DM* ele_val, double* stress_n, double* dstrain,
                                          double temperature)
{
	int i, ns, dim;
	double norn_S, norm_str;
	static double epsilon_tr[6];
	ns = ele_val->xi->Size();
	dim = 2;
	if (ns > 4)
		dim = 3;
	DeviatoricStress(stress_n);
	norn_S = sqrt(1.5 * TensorMutiplication2(stress_n, stress_n, dim));
	// WX:12.2012 threshold dev str for lubby2
	if (norn_S <= threshold_dev_str)
		return;
	//
	for (i = 0; i < ns; i++)
		epsilon_tr[i] = (*ele_val->xi)(i);
	//
	norm_str = sqrt(2.0 * TensorMutiplication2(epsilon_tr, epsilon_tr, dim) / 3.0);
	double max_etr = norn_S / ((*data_Creep)(6, 0) * exp((*data_Creep)(4, 0) * norn_S)); // WX:12.2012 bug fixed
	double eta_k = (*data_Creep)(3, 0) * exp((*data_Creep)(5, 0) * norn_S);
	double eta_m
	    = (*data_Creep)(0, 0) * exp((*data_Creep)(1, 0) * norn_S) * exp((temperature + 273.16) * (*data_Creep)(2, 0));
	if (max_etr < DBL_EPSILON)
		return;
	if (threshold_dev_str >= 0)
		norm_str = min(norm_str, max_etr);
	// WX:12.2012, change "norm_str" to "min(norm_str,max_etr)"
	// the norm_str should not be higher than max_etr
	for (i = 0; i < ns; i++)
	{
		(*data_Creep)(i, 1) += 1.5 * dt * (1 - norm_str / max_etr) * stress_n[i] / eta_k;
		dstrain[i] -= 1.5 * dt * ((1 - norm_str / max_etr) / eta_k + 1 / eta_m) * stress_n[i];
	}
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   11/2005 CMCD From fluid properties
   last modified:
**************************************************************************/
void CSolidProperties::CalPrimaryVariable(vector<string>& pcs_name_vector)
{
	CRFProcess* m_pcs = NULL;

	int nidx0, nidx1;
	if (!Fem_Ele_Std) // OK
		return;

	for (size_t i = 0; i < pcs_name_vector.size(); i++)
	{
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
	}
}

/**************************************************************************
   FEMLib-Method:
   01/2006 OK Implementation
   05/2009 OK DENSITY
   09/2009 OK Bugfix, write only existing data
**************************************************************************/
void CSolidProperties::Write(std::fstream* msp_file)
{
	//----------------------------------------------------------------------
	// KEYWORD
	*msp_file << "#SOLID_PROPERTIES"
	          << "\n";
	//-----------------------------------------------------------------------
	// NAME
	if (name.length() > 0)
	{
		*msp_file << " $NAME"
		          << "\n";
		*msp_file << "  ";
		*msp_file << name << "\n";
	}
	//-----------------------------------------------------------------------
	// GEO_TYPE
	//-----------------------------------------------------------------------
	// DIMENSION
	//-----------------------------------------------------------------------
	// PROPERTIES
	//.......................................................................
	*msp_file << " $DENSITY"
	          << "\n";
	*msp_file << "  " << Density_mode;
	if (data_Density) // OK410
		*msp_file << " " << (*data_Density)(0) << "\n";
	else
		*msp_file << " Warning: no density data"
		          << "\n";
	//.......................................................................
	// Elasticity properties
	*msp_file << " $ELASTICITY"
	          << "\n";
	if (Poisson_Ratio()) // OK410
		*msp_file << "  POISSION " << Poisson_Ratio() << "\n";
	*msp_file << "  YOUNGS_MODULUS"
	          << "\n";
	*msp_file << "  " << Youngs_mode;
	if (data_Youngs) // OK410
		*msp_file << " " << (*data_Youngs)(0) << "\n";
	//.......................................................................
	// Thermal properties
	*msp_file << " $THERMAL"
	          << "\n";
	if (ThermalExpansion >= 0)
	{
		*msp_file << "  EXPANSION"
		          << "\n";
		*msp_file << "  " << ThermalExpansion << "\n";
	}
	if (Capacity_mode > 0)
	{
		*msp_file << "  CAPACITY"
		          << "\n";
		*msp_file << "  " << Capacity_mode;
		*msp_file << " " << (*data_Capacity)(0) << "\n";
	}
	if (this->Conductivity_mode > 0) // OK410
	{
		*msp_file << "  CONDUCTIVITY"
		          << "\n";
		*msp_file << "  " << Conductivity_mode;
		*msp_file << " " << (*data_Conductivity)(0) << "\n";
	}
	//-----------------------------------------------------------------------
}
/**************************************************************************
   FEMLib-Method:
   03/2008 WW Implementation
**************************************************************************/
void CSolidProperties::TEPSwellingParameter_kis(const double suction)
{
	double val = 0.;
	double alf_i = (*data_Youngs)(1);
	// k_(s)
	// if(suction<=-0.999*1e-6/alf_i)
	val = 1.0 + alf_i * suction;
	// else
	//   val = 0.001;
	if (val < 0.0)
		val = 0.001;
	// Swelling index. Kappa
	(*data_Youngs)(6) = val * (*data_Youngs)(0);
}

/**************************************************************************
   FEMLib-Method:
   03/2008 WW Implementation
**************************************************************************/
double CSolidProperties::TEPSwellingParameter(const double mean_stress)
{
	double val = 0.;
	double alf_sp = (*data_Youngs)(3);
	double pref = (*data_Youngs)(5);
	double e0 = (*data_Plasticity)(4);
	double suction = (*data_Youngs)(7);
	// k_(s)
	TEPSwellingParameter_kis(suction);
	//
	if (mean_stress < 1.0e-20)
		val = 1.0 + alf_sp * log(1.0e-20 / pref);
	else if (mean_stress >= pref * exp(-1.0 / alf_sp))
	{
		val = 0.;
		return val;
	}
	else
		val = 1. + alf_sp * log(mean_stress / pref);
	return val * (*data_Youngs)(2) * exp((*data_Youngs)(4) * suction) / (3. + 3. * e0);
}

/************************************************************************
11.2011 WX calculate the transform matrix for micro structure tensor
*************************************************************************/
void CSolidProperties::CalTransMatrixMicroStru(Matrix* Trans, double* bedding_norm)
{
	double tmp_dir[9] = {0.};
	double i, j, k, tmp_sqrt;
	i = bedding_norm[0];
	j = bedding_norm[1];
	k = bedding_norm[2];
	tmp_sqrt = sqrt(1 - k * k);
	if (fabs(i) == 1)
		tmp_dir[2] = tmp_dir[3] = tmp_dir[7] = 1.;
	else if (fabs(j) == 1)
		tmp_dir[1] = tmp_dir[5] = tmp_dir[6] = 1.;
	else if (fabs(k) == 1)
		tmp_dir[0] = tmp_dir[4] = tmp_dir[8] = 1.;
	else
	{
		tmp_dir[0] = j / tmp_sqrt;
		tmp_dir[1] = k * i / tmp_sqrt;
		tmp_dir[2] = i;
		tmp_dir[3] = -i / tmp_sqrt;
		tmp_dir[4] = k * j / tmp_sqrt;
		tmp_dir[5] = j;
		tmp_dir[6] = 0;
		tmp_dir[7] = -(i * i + j * j) / tmp_sqrt;
		tmp_dir[8] = k;
	}
	for (int ii = 0; ii < 9; ii++)
	{
		if (fabs(tmp_dir[ii]) > 1)
			cout << "WARNING: Please check the input value for $BEDDING_NORM !" << endl;
		break;
	}
	CalTransMatrixA(tmp_dir, Trans, 6);
}
// WX:06.2012 E dependents on stress, strain ...
double CSolidProperties::E_Function(int dim, const ElementValue_DM* ele_val, int ngp)
{
	int nGauss = ngp, ele_dim = dim, valid = 1, size = 6;
	double tmp_stress = 0., return_value = 1., I1 = 0.;
	double stress[6] = {0.};
	if (ele_dim == 2)
		size = 4;
	switch (E_Function_Model)
	{
		case 1:
			for (int i = 0; i < 3; i++)
			{
				tmp_stress = 0.;
				for (int j = 0; j < nGauss; j++)
					tmp_stress += (*ele_val->Stress)(i, j) / nGauss;
				stress[i] = tmp_stress;
			}
			I1 = (stress[0] + stress[1] + stress[2]) / 3.;
			return_value = GetCurveValue((int)E_Function_Model_Value[0], 0, I1, &valid);
			break;
		case 2:
		{
			double prin_str[6] = {0.}, prin_dir[6] = {0.};
			for (int i = 0; i < size; i++)
			{
				tmp_stress = 0.;
				for (int j = 0; j < nGauss; j++)
					tmp_stress += (*ele_val->Stress)(i, j) / nGauss;
				stress[i] = tmp_stress;
			}
			CalPrinStrDir(stress, prin_str, prin_dir, size);
			return_value = GetCurveValue((int)E_Function_Model_Value[0], 0, prin_str[0], &valid);
		}
		default:
			return_value = 1.;
			break;
	}
	return return_value;
}

// TN - added for TNEQ/TEQ process
void CSolidProperties::setSolidReactiveSystem(FiniteElement::SolidReactiveSystem reactive_system)
{
	_reactive_system = reactive_system;
}

FiniteElement::SolidReactiveSystem CSolidProperties::getSolidReactiveSystem() const
{
	return _reactive_system;
}

} // end namespace

/////////////////////////////////////////////////////////////////////////////

/**************************************************************************
   FEMLib-Method:
   Task: Master read function
   Programing:
   08/2004 OK Implementation for fluid properties
   08/2004 WW Modification for solid properties
   01/2005 OK Boolean type
   01/2005 OK Destruct before read
**************************************************************************/
bool MSPRead(const std::string& given_file_base_name)
{
	//----------------------------------------------------------------------
	// OK  MSPDelete();
	//----------------------------------------------------------------------
	SolidProp::CSolidProperties* m_msp = NULL;
	char line[MAX_ZEILE];
	std::string sub_line;
	std::string line_string;
	std::ios::pos_type position;
	//========================================================================
	// File handling
	std::string msp_file_name = given_file_base_name + MSP_FILE_EXTENSION;
	std::ifstream msp_file(msp_file_name.data(), std::ios::in);
	if (!msp_file.good())
		return false;
	msp_file.seekg(0L, std::ios::beg);
	//========================================================================
	// Keyword loop
	std::cout << "MSPRead"
	          << "\n";
	while (!msp_file.eof())
	{
		msp_file.getline(line, MAX_ZEILE);
		line_string = line;
		if (line_string.find("#STOP") != string::npos)
			return true;
		//----------------------------------------------------------------------
		// keyword found
		if (line_string.find("#SOLID_PROPERTIES") != std::string::npos)
		{
			m_msp = new SolidProp::CSolidProperties();
			m_msp->file_base_name = given_file_base_name;
			position = m_msp->Read(&msp_file);
			msp_vector.push_back(m_msp);
			msp_file.seekg(position, std::ios::beg);
		} // keyword found
	} // eof
	return true;
	//========================================================================
}

/**************************************************************************
   ROCKFLOW - Funktion: TensorMutiplication2

   Aufgabe:
   Calculate tensor mutiplication: a_{ij}b_{ij}

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E
   const double *s1: Stress tensor in an array (s1_11, s1_22, s1_12, s1_33)
   const double *s2: Stress tensor in an array (s2_11, s2_22, s2_12, s2_33)

   Ergebnis:
   - double - The value of third stress invariant

   Programmaenderungen:
   04/2003   WW  Erste Version

**************************************************************************/
double TensorMutiplication2(const double* s1, const double* s2, const int Dim)
{
	switch (Dim)
	{
		case 2:
			return s1[0] * s2[0] + s1[1] * s2[1] + s1[2] * s2[2] + 2.0 * s1[3] * s2[3];
			break;
		case 3:
			return s1[0] * s2[0] + s1[1] * s2[1] + s1[2] * s2[2] + 2.0 * s1[3] * s2[3] + 2.0 * s1[4] * s2[4]
			       + 2.0 * s1[5] * s2[5];
			break;
	}
	return 0.0; // To avoid warnings
}

/**************************************************************************
   GeoSys: Norm of stress
   06/2008   WW  Programming

**************************************************************************/
double StressNorm(const double* s, const int Dim)
{
	double val = 0.0;
	double mean_s = (s[0] + s[1] + s[2]) / 3.0;
	double s1 = s[0] - mean_s;
	double s2 = s[1] - mean_s;
	double s3 = s[2] - mean_s;
	val = s1 * s1 + s2 * s2 + s3 * s3 + 2.0 * s[3] * s[3];
	if (Dim == 3)
		val += +2.0 * s[4] * s[4] + 2.0 * s[5] * s[5];
	return val;
}

/**************************************************************************
   ROCKFLOW - Funktion: TensorMutiplication3

   Aufgabe:
   Calculate tensor mutiplication: a_{ij}b_{jk}c_{ki} for plane strain problem

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E
   const double *s1: Stress tensor in an array (s1_11, s1_22, s1_12, s1_33)
   const double *s2: Stress tensor in an array (s2_11, s2_22, s2_12, s2_33)
   const double *s3: Stress tensor in an array (s3_11, s3_22, s3_12, s3_33)

   Ergebnis:
   - double - The value of third stress invariant

   Programmaenderungen:
   04/2003   WW  Erste Version

**************************************************************************/
double TensorMutiplication3(const double* s1, const double* s2, const double* s3, const int Dim)
{
	switch (Dim)
	{
		case 2:
			return (s1[0] * (s2[0] * s3[0] + s2[3] * s3[3]) // s1_11*(s2_11*s3_11+s2_12*s3_21)
			        + s1[3] * (s2[0] * s3[3] + s2[3] * s3[1]) // s1_12*(s2_11*s3_12+s2_12*s3_22)
			        + s1[3] * (s2[3] * s3[0] + s2[1] * s3[3]) // s1_21*(s2_21*s3_11+s2_22*s3_21)
			        + s1[1] * (s2[3] * s3[3] + s2[1] * s3[1]) // s1_22*(s2_21*s3_12+s2_22*s3_22)
			        + s1[2] * s2[2] * s3[2])
			       / 3.0; // s33*s33*s33
			break;
		case 3:
			return (
			           // s1_11*(s2_11*s3_11+s2_12*s3_21+s2_13*s3_31)
			           s1[0] * (s2[0] * s3[0] + s2[3] * s3[3] + s2[4] * s3[4])
			           // s1_12*(s2_11*s3_12+s2_12*s3_22+s2_13*s3_32)
			           + s1[3] * (s2[0] * s3[3] + s2[3] * s3[1] + s2[4] * s3[5])
			           // s1_13*(s2_11*s3_13+s2_12*s3_23+s2_13*s3_33)
			           + s1[4] * (s2[0] * s3[4] + s2[3] * s3[5] + s2[4] * s3[2])
			           // s1_21*(s2_21*s3_11+s2_22*s3_21+s2_23*s3_31)
			           + s1[3] * (s2[3] * s3[0] + s2[1] * s3[3] + s2[5] * s3[4])
			           // s1_22*(s2_21*s3_12+s2_22*s3_22+s2_23*s3_32)
			           + s1[1] * (s2[3] * s3[3] + s2[1] * s3[1] + s2[5] * s3[5])
			           // s1_23*(s2_21*s3_13+s2_22*s3_23+s2_23*s3_33)
			           + s1[5] * (s2[3] * s3[4] + s2[1] * s3[5] + s2[5] * s3[2])
			           // s1_31*(s2_31*s3_11+s2_32*s3_21+s2_33*s3_31)
			           + s1[4] * (s2[4] * s3[0] + s2[5] * s3[3] + s2[2] * s3[4]) // WX:bug fixed s3_31 is s3[4]
			           // s1_32*(s2_31*s3_12+s2_32*s3_22+s2_33*s3_32)
			           + s1[5] * (s2[4] * s3[3] + s2[5] * s3[1] + s2[2] * s3[5])
			           // s1_33*(s2_31*s3_13+s2_32*s3_23+s2_33*s3_33)
			           + s1[2] * (s2[4] * s3[4] + s2[5] * s3[5] + s2[2] * s3[2]))
			       / 3.0;
			break;
	}
	return 0.0; // To avoid warnings
}

/**************************************************************************
   ROCKFLOW - Funktion: DeviatoricStress

   Aufgabe:
   Computing the deviatoric stresses
   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E :
     double *Stress                   :   Stresses

   Ergebnis:
   - double* - The first stress invariant

   Programmaenderungen:
   08/2003   WW  Erste Version (for 2D)

**************************************************************************/
double DeviatoricStress(double* Stress)
{
	int i;
	double I1 = Stress[0] + Stress[1] + Stress[2];
	for (i = 0; i < 3; i++)
		Stress[i] -= I1 / 3.0;

	return I1;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   01/2005 OK Implementation
   last modified:
**************************************************************************/
void MSPDelete()
{
	long i;
	int no_msp = (int)msp_vector.size();
	for (i = 0; i < no_msp; i++)
		delete msp_vector[i];
	msp_vector.clear();
}

/**************************************************************************
   FEMLib-Method:
   01/2006 OK Implementation
**************************************************************************/
void MSPWrite(const std::string& base_file_name)
{
	SolidProp::CSolidProperties* m_msp = NULL;
	//----------------------------------------------------------------------
	// File handling
	std::fstream msp_file;
	std::string msp_file_name = base_file_name + MSP_FILE_EXTENSION;
	msp_file.open(msp_file_name.data(), ios::trunc | ios::out);
	msp_file.setf(ios::scientific, ios::floatfield);
	msp_file.precision(12);
	if (!msp_file.good())
		return;
	//----------------------------------------------------------------------
	msp_file << "GeoSys-MSP: Material Solid Properties -------------"
	         << "\n";
	//----------------------------------------------------------------------
	for (int i = 0; i < (int)msp_vector.size(); i++)
	{
		m_msp = msp_vector[i];
		m_msp->Write(&msp_file);
	}
	msp_file << "#STOP";
	msp_file.close();
	//----------------------------------------------------------------------
}
