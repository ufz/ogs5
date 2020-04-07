/**
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   FEMLib-Object: MAT-MP
   Task: MediumProperties
   Programing:
   01/2004 OK Implementation
**************************************************************************/

//#include "makros.h"
// C++ STL
//#include <iostream>
#include <cfloat>
#include "display.h"

// FEMLib
#include "tools.h"
//#include "rf_pcs.h"
//#include "femlib.h"
extern double* GEOGetELEJacobianMatrix(long number, double* detjac);
#include "mathlib.h"
//#include "rf_mfp_new.h"
#include "rf_msp_new.h"
//#include "material.h"
#include "rf_tim_new.h"
#include "rfmat_cp.h"
extern double gravity_constant;
// using SolidProp::CSolidProperties;
// LIB
#include "files0.h"
// this
#include "rf_mmp_new.h"
//#include "rf_react.h"
// Gauss point veclocity
#include "fem_ele_std.h"
#include "fem_ele_vec.h"
// MSHLib
//#include "msh_lib.h"
#include "pcs_dm.h"  //WX

#include "PhysicalConstant.h"

// MAT-MP data base lists
list<string> keywd_list;
list<string> mat_name_list;
list<char*> mat_name_list_char;
list<CMediumProperties*> db_mat_mp_list;
// MAT-MP list
vector<CMediumProperties*> mmp_vector;
list<CMediumPropertiesGroup*> mmp_group_list;

using namespace std;
using namespace PhysicalConstant;
using FiniteElement::CElement;
using FiniteElement::CFiniteElementStd;
using FiniteElement::ElementValue;
using FiniteElement::ElementValue_DM;
using namespace Display;

/**************************************************************************
   FEMLib-Method: CMediumProperties
   Task: constructor
   Programing:
   02/2004 OK Implementation
   last modification:
**************************************************************************/
CMediumProperties::CMediumProperties()
    : geo_dimension(0), _mesh(NULL), _geo_type(GEOLIB::GEODOMAIN)
{
    name = "DEFAULT";
    mode = 0;
    selected = false;
    m_pcs = NULL;  // OK
    // GEO
    geo_area = 1.0;  // OK
    porosity_model = -1;
    porosity_model_values[0] = 0.1;
    tortuosity_model = -1;
    tortuosity_model_values[0] = 1.;
    // flow
    storage_model = -1;
    storage_model_values[0] = 0.;
    permeability_model = -1;
    permeability_tensor_type = 0;
    tortuosity_tensor_type = 0;
    permeability_tensor[0] = 1.e-13;
    local_permeability = 1e-13;
    residual_saturation[0] = 0.0;  // sgr: residual saturation, this phase
    maximum_saturation[0] = 1.0;   // sgm: maximum saturation, this phase
    saturation_exponent[0] =
        1.0;  // (set exponent = 1 results in a linear k_rel function)
    conductivity_model = -1;
    flowlinearity_model = 0;
    capillary_pressure_model = -1;
    capillary_pressure_values[4] = 1.0 / DBL_EPSILON;  // JT: max Pc
    entry_pressure_conversion = false;
    permeability_saturation_model[0] = -1;
    minimum_relative_permeability = 1.0e-9;  // JT: the default value
    unconfined_flow_group = -1;
    permeability_stress_mode = -1;  // WW
    c_coefficient = NULL;           // WW
    // surface flow
    friction_coefficient = -1;
    //  friction_model = -1;
    // mass transport
    // heat transport
    heat_dispersion_model = -1;  // WW
    heat_dispersion_longitudinal = 0.;
    heat_dispersion_transverse = 0.;
    lgpn = 0.0;
    mass_dispersion_transverse = 0.0;
    mass_dispersion_longitudinal = 0.0;
    heat_diffusion_model = -1;                  // WW
    base_heat_diffusion_coefficient = 2.16e-5;  // JM
    geo_area = 1.0;
    //	geo_type_name = "DOMAIN";             //OK
    vol_mat = 0.0;
    vol_bio = 0.0;
    vol_mat_model = 0;
    vol_bio_model = 0;
    foc = 0.0;
    alpha_t_model = -1;
    graindiameter = 0;  // CB Chiogna et al alpha-t model
    hydraulicrad = 0;
    betaexpo = 0;
    ElementVolumeMultiplyer = 1.0;  // SB / JOD 2014-11-10

    permeability_pressure_model = -1;  // 01.09.2011. WW
    permeability_strain_model = -1;    // 01.09.2011. WW
    forchheimer_cf = 0.0;              // NW
    forchheimer_De = .0;
    forchheimer_a1 = .0;
    forchheimer_a2 = .0;
    heat_transfer_model = 0;
    effective_heat_transfer_model = 0;
    heat_transfer_model_value = .0;
    particle_diameter_model = 0;
    particle_diameter_model_value = .0;

    PhaseHeatedByFriction = "SOLID";
    _fric_phase = FiniteElement::SOLID;
    storage_effstress_model = 0;
    permeability_effstress_model = 0;
    evaporation = -1;
}

/**************************************************************************
   FEMLib-Method: CMediumProperties
   Task: destructor
   Programing:
   02/2004 OK Implementation
   last modification:
**************************************************************************/
CMediumProperties::~CMediumProperties(void)
{
    if (c_coefficient)
        delete[] c_coefficient;  // WW
    geo_name_vector.clear();
}

////////////////////////////////////////////////////////////////////////////
// IO functions
////////////////////////////////////////////////////////////////////////////

/**************************************************************************
   FEMLib-Method: MMPRead
   Task: master read functionn
   Programing:
   02/2004 OK Implementation
   last modification:
**************************************************************************/
bool MMPRead(std::string base_file_name)
{
    //----------------------------------------------------------------------
    // OK  MMPDelete();
    //----------------------------------------------------------------------
    ScreenMessage("MMPRead ... ");;
    CMediumProperties* m_mat_mp = NULL;
    char line[MAX_ZEILE];
    std::string sub_line;
    std::string line_string;
    std::ios::pos_type position;
    //========================================================================
    // file handling
    std::string mp_file_name;
    mp_file_name = base_file_name + MMP_FILE_EXTENSION;
    std::ifstream mp_file(mp_file_name.data(), std::ios::in);
    if (!mp_file.good())
    {
        ScreenMessage("! Error in MMPRead: No material data !\n");
        return false;
    }
    mp_file.seekg(0L, std::ios::beg);
    //========================================================================
    // keyword loop
    while (!mp_file.eof())
    {
        mp_file.getline(line, MAX_ZEILE);
        line_string = line;
        if (line_string.find("#STOP") != string::npos)
        {
            ScreenMessage("done, read %d medium properties\n",
                                   mmp_vector.size());

            return true;
        }
        //----------------------------------------------------------------------
        // keyword found
        if (line_string.find("#MEDIUM_PROPERTIES") != string::npos)
        {
            m_mat_mp = new CMediumProperties();
            position = m_mat_mp->Read(&mp_file);
            // OK41
            m_mat_mp->number = (int)mmp_vector.size();
            mmp_vector.push_back(m_mat_mp);
            mp_file.seekg(position, std::ios::beg);
        }  // keyword found
    }      // eof
    return true;
    // Tests
}

/**************************************************************************
   FEMLib-Method: CMediumProperties::Read
   Task: read functionn
   Programing:
   02/2004 OK Template
   08/2004 CMCD Implementation
   10/2004 MX/OK Porosity model 3, swelling
   11/2004 CMCD String streaming
   07/2005 MB porosity_file, permeability_file, GEO_TYPE layer
   10/2005 OK GEO_TYPE geo_name_vector
   01/2006 YD PCS_TYPE
   05/2007 PCH Tortuosity tensor
   05/2007 WW Stress permeability coorector. Two models.
   last modification:
**************************************************************************/
// Order of Key Words
/*
         0. $NAME
            (i)    _BORDEN
         1. $GEOTYPE
            (i)		_CLAY
            (ii)	_SILT
            (iii)	_SAND
            (iv)	_GRAVEL
            (v)		_CRYSTALINE
         2. $GEOMETRY
            (i)		_DIMENSION
   (ii)	_AREA
   3. $POROSITY
   4. $TORTUOSITY
   5. $MOBILE_IMOBILE_MODEL
   6. $LITHOLOGY_GRAIN_CLASS
   7. $FLOWLINEARITY
   8. $SORPTION_MODEL
   9. $STORAGE
   11.$PERMEABILITY
   12.$PERMEABILITY_FUNCTION_
   (1)		DEFORMATION
   (2)  PRESSURE
   (3)	    SATURATION
   (4)	    STRESS
   (5)		VELOCITY
   (6)		POROSITY
   13.$CAPILLARY_PRESSURE
   14.$MASSDISPERSION
   (i)		_LONGITUDINAL
   (ii)	_TRANSVERSE
   15.$HEATDISPERSION
   (i)		_LONGITUDINAL
   (ii)	_TRANSVERSE
   19.$ELECTRIC_CONDUCTIVITY
   20.$UNCONFINED_FLOW_GROUP
   21.$FLUID_EXCHANGE_WITH_OTHER_CONTINUA
 */
std::ios::pos_type CMediumProperties::Read(std::ifstream* mmp_file)
{
    int i, j, k = 0;
    std::string line_string;
    std::stringstream in;
    std::ios::pos_type position;
    std::string dollar("$");
    std::string hash("#");
    // WW bool new_subkeyword = false;
    bool new_keyword = false;
    std::string m_string;
    // WW
    std::stringstream buff;
    std::vector<string> tokens;
    char* pch;
    char seps[] = "+\n";
    char seps1[] = "*";
    double f_buff;

    while (!new_keyword)
    {
        // WW new_subkeyword = false;
        position = mmp_file->tellg();
        line_string = GetLineFromFile1(mmp_file);
        if (line_string.size() < 1)
            break;
        if (line_string.find(hash) != std::string::npos)
        {
            new_keyword = true;
            break;
        }
        //--------------------------------------------------------------------
        // PCS                         //YD
        if (line_string.find("$PCS_TYPE") != std::string::npos)
        {
            in.str(GetLineFromFile1(mmp_file));
            in >> pcs_type_name;
            in.clear();
            continue;
        }

        // NAME
        // subkeyword found
        if (line_string.find("$NAME") != std::string::npos)
        {
            in.str(GetLineFromFile1(mmp_file));
            in >> name;  // sub_line
            in.clear();
            continue;
        }
        //--------------------------------------------------------------------
        // GEO_TYPE

        // subkeyword found
        if (line_string.find("$GEO_TYPE") != std::string::npos)
        {
            while (!(m_string.find("$") != std::string::npos) &&
                   (!(m_string.find("#") != std::string::npos)))
            {
                position = mmp_file->tellg();
                in.str(GetLineFromFile1(mmp_file));
                in >> m_string >> geo_name;
                in.clear();
                if (!(m_string.find("$") != std::string::npos) &&
                    (!(m_string.find("#") != std::string::npos)))
                {
                    //					geo_type_name = m_string;
                    _geo_type = GEOLIB::convertGeoType(m_string);
                    geo_name_vector.push_back(geo_name);
                }
            }
            mmp_file->seekg(position, std::ios::beg);
            continue;
        }
        //....................................................................
        // ToDo to GeoLib
        // 2i..GEOMETRY_DIMENSION
        //------------------------------------------------------------------------
        // subkeyword found
        if (line_string.find("$GEOMETRY_DIMENSION") != std::string::npos)
        {
            in.str(GetLineFromFile1(mmp_file));
            in >> geo_dimension;
            in.clear();
            continue;
        }
        //------------------------------------------------------------------------
        // subkeyword found
        if (line_string.find("$GEOMETRY_INCLINATION") != std::string::npos)
        {
            in.str(GetLineFromFile1(mmp_file));
            in >> geo_inclination;
            in.clear();
            continue;
        }
        //------------------------------------------------------------------------
        // ToDo to GeoLib
        // 2ii..GEOMETRY_AREA
        //------------------------------------------------------------------------
        // subkeyword found
        if (line_string.find("$GEOMETRY_AREA") != std::string::npos)
        {
            in.str(GetLineFromFile1(mmp_file));
            in >> line_string;
            if (line_string.find("FILE") != string::npos)
            {
                in >> geo_area_file;
                geo_area_file = FilePath + geo_area_file;
                // End of new lines
            }
            else
                geo_area = strtod(line_string.data(), NULL);
            in.clear();
            continue;
        }
        //------------------------------------------------------------------------
        // 3..POROSITY
        //------------------------------------------------------------------------
        // CB
        // subkeyword found
        if ((line_string.find("$POROSITY") != string::npos) &&
            (!(line_string.find("_DISTRIBUTION") != std::string::npos)))
        {
            in.str(GetLineFromFile1(mmp_file));
            in >> porosity_model;
            switch (porosity_model)
            {
                case 0:  // n=f(x)
                    in >> porosity_curve;
                    break;
                case 1:  // n=const
                    in >> porosity_model_values[0];
                    break;
                case 2:  // f(normal effective stress for fracture systems)
                    in >> porosity_model_values[0];
                    in >> porosity_model_values[1];
                    in >> porosity_model_values[2];
                    in >> porosity_model_values[3];
                    porosity_pcs_name_vector.push_back("PRESSURE1");
                    break;
                case 3:                              // Chemical swelling model
                    in >> porosity_model_values[0];  // Initial porosity
                    in >> porosity_model_values[1];  // Specific surface[m^2/g]
                    in >> porosity_model_values[2];  // Expansive min. fragtion
                    in >> porosity_model_values[3];  // m
                    in >> porosity_model_values[4];  // I
                    in >> porosity_model_values[5];  // S^l_0
                    in >> porosity_model_values[6];  // beta
                    porosity_pcs_name_vector.push_back("SATURATION2");
                    porosity_pcs_name_vector.push_back("TEMPERATURE1");
                    break;
                case 4:  // Chemical swelling model (constrained swelling,
                         // constant I)
                    in >> porosity_model_values[0];  // Initial porosity
                    in >> porosity_model_values[1];  // Specific surface[m^2/g]
                    in >> porosity_model_values[2];  // Expansive min. fragtion
                    in >> porosity_model_values[3];  // m
                    in >> porosity_model_values[4];  // I
                    in >> porosity_model_values[5];  // S^l_0
                    in >> porosity_model_values[6];  // beta
                    in >> porosity_model_values[7];  // n_min
                    // for richard flow only
                    porosity_pcs_name_vector.push_back("SATURATION1");
                    porosity_pcs_name_vector.push_back("TEMPERATURE1");
                    break;
                case 5:  // Chemical swelling model (free swelling, constant I)
                    in >> porosity_model_values[0];  // Initial porosity
                    in >> porosity_model_values[1];  // Specific surface[m^2/g]
                    in >> porosity_model_values[2];  // Expansive min. fragtion
                    in >> porosity_model_values[3];  // m
                    in >> porosity_model_values[4];  // I
                    in >> porosity_model_values[5];  // S^l_0
                    in >> porosity_model_values[6];  // beta
                    porosity_pcs_name_vector.push_back("SATURATION2");
                    porosity_pcs_name_vector.push_back("TEMPERATURE1");
                    break;
                case 6:  // Chemical swelling model (constrained swelling)
                    in >> porosity_model_values[0];  // Initial porosity
                    in >> porosity_model_values[1];  // Specific surface[m^2/g]
                    in >> porosity_model_values[2];  // Expansive min. fragtion
                    in >> porosity_model_values[3];  // m
                    in >> porosity_model_values[4];  // I
                    in >> porosity_model_values[5];  // S^l_0
                    in >> porosity_model_values[6];  // beta
                    in >> porosity_model_values[7];  // n_min
                    porosity_pcs_name_vector.push_back("SATURATION2");
                    porosity_pcs_name_vector.push_back("TEMPERATURE1");
                    break;
                case 7:  // n=f(stress_mean) WW
                    in >> porosity_curve;
                    break;
                case 10:  // Chemical swelling model (constrained swelling,
                          // constant I)
                {
                    int m;
                    in >> porosity_model_values[0];  // Initial porosity
                    in >> m;                         // m
                    if (m > 15)
                        std::cout << "Maximal number of solid phases is now "
                                     "limited to be 15!!!"
                                  << "\n";
                    for (int i = 0; i < m + 1; i++)
                        // molar volume [l/mol]
                        in >> porosity_model_values[i + 1];
                    break;
                }
                case 11:  // MB: read from file ToDo
                    // in >> porosity_file; // CB
                    in >> porosity_model_values[0];  // CB some dummy default
                                                     // value is read
                    // CB $POROSITY_DISTRIBUTION should be given as keyword in
                    // *.mmp file,
                    //     porosities then are to be read in from file by fct.
                    //     CMediumProperties::SetDistributedELEProperties
                    break;
                case 12:
                    in >> porosity_model_values[0];  // WX 03.2011, dependent on
                                                     // strain
                    break;
                case 13:  // mineral precipitation / dissolution by ABM
                    in >> porosity_model_values[0];  // Initial porosity
                    break;
#ifdef GEM_REACT
                case 15:
                    in >> porosity_model_values[0];  // set a default value for
                                                     // GEMS calculation
                    // save this seperately;
                    break;
#endif
#ifdef BRNS
                case 16:
                    in >> porosity_model_values[0];  // set a default value for
                                                     // BRNS calculation
                    break;
#endif
                default:
                    std::cerr << "Error in MMPRead: no valid porosity model"
                              << "\n";
                    break;
            }
            in.clear();
            continue;
        }
        // subkeyword found
        if (line_string.find("$VOL_MAT") != string::npos)
        {
            in.str(GetLineFromFile1(mmp_file));
            // in >> idummy >> this->vol_mat; CB
            in >> vol_mat_model >> this->vol_mat;
            switch (vol_mat_model)
            {
                case 1:  // do nothing
                    break;
                case 2:  // do nothing
                    break;
                default:
                    std::cout << "Error in MMPRead: no valid vol_mat_model"
                              << "\n";
                    break;
            }
            in.clear();
            continue;
        }
        // subkeyword found
        if (line_string.find("$VOL_BIO") != string::npos)
        {
            in.str(GetLineFromFile1(mmp_file));
            // in >> idummy >> this->vol_bio; CB
            in >> vol_bio_model >> this->vol_bio;
            switch (vol_bio_model)
            {
                case 1:  // do nothing
                    break;
                case 2:  // do nothing
                    break;
                default:
                    std::cout << "Error in MMPRead: no valid vol_bio_model"
                              << "\n";
                    break;
            }
            in.clear();
            continue;
        }
        //------------------------------------------------------------------------
        // 4..TORTUOSITY
        //------------------------------------------------------------------------
        // subkeyword found
        if (line_string.find("$TORTUOSITY") != std::string::npos)
        {
            in.str(GetLineFromFile1(mmp_file));
            in >> tortuosity_tensor_type_name;
            switch (tortuosity_tensor_type_name[0])
            {
                case '0':  // n=f(x) <- Case zero
                    break;
                case '1':  // n=const
                    tortuosity_model = 1;
                    in >> tortuosity_model_values[0];
                    break;
                case '2':  // t=n*factor
                    tortuosity_model = 2;
                    in >> tortuosity_model_values[0];
                    break;
                case 'I':  // isotropic
                    tortuosity_tensor_type = 0;
                    tortuosity_model = 1;
                    in >> tortuosity_model_values[0];
                    // CMCD to pick up 2D and 3D Isotropic case.
                    tortuosity_model_values[1] = tortuosity_model_values[2] =
                        tortuosity_model_values[0];
                    break;
                case 'O':  // orthotropic  <- Case Alphabet O
                    tortuosity_model = 1;
                    tortuosity_tensor_type = 1;
                    if (geo_dimension == 0)
                        std::cout << "Error in CMediumProperties::Read: no "
                                     "geometric dimension"
                                  << "\n";
                    if (geo_dimension == 2)
                    {
                        in >> tortuosity_model_values[0];
                        in >> tortuosity_model_values[1];
                    }
                    if (geo_dimension == 3)
                    {
                        in >> tortuosity_model_values[0];
                        in >> tortuosity_model_values[1];
                        in >> tortuosity_model_values[2];
                    }
                    break;
                case 'A':  // anisotropic
                    tortuosity_model = 1;
                    tortuosity_tensor_type = 2;
                    if (geo_dimension == 0)
                        std::cout << "Error in CMediumProperties::Read: no "
                                     "geometric dimension"
                                  << "\n";
                    if (geo_dimension == 2)
                    {
                        in >> tortuosity_model_values[0];
                        in >> tortuosity_model_values[1];
                        in >> tortuosity_model_values[2];
                        in >> tortuosity_model_values[3];
                    }
                    if (geo_dimension == 3)
                    {
                        in >> tortuosity_model_values[0];
                        in >> tortuosity_model_values[1];
                        in >> tortuosity_model_values[2];
                        in >> tortuosity_model_values[3];
                        in >> tortuosity_model_values[4];
                        in >> tortuosity_model_values[5];
                        in >> tortuosity_model_values[6];
                        in >> tortuosity_model_values[7];
                        in >> tortuosity_model_values[8];
                    }
                    break;
                case 'F':                  // SB: read from file
                    tortuosity_model = 2;  // OK
                    in >> permeability_file;
                    break;
                default:
                    std::cout
                        << "Error in MMPRead: no valid tortuosity tensor type"
                        << "\n";
                    break;
            }
            in.clear();
            continue;
        }

        //------------------------------------------------------------------------
        // 5..MOBILE_IMOBILE_MODEL
        //------------------------------------------------------------------------
        // To do as necessary

        //------------------------------------------------------------------------
        // 6..LITHOLOGY_GRAIN_CLASS
        //------------------------------------------------------------------------
        // To do as necessary

        //------------------------------------------------------------------------
        // 7..FLOWLINEARITY
        //------------------------------------------------------------------------
        // subkeyword found
        if (line_string.find("$FLOWLINEARITY") != std::string::npos)
        {
            in.str(GetLineFromFile1(mmp_file));
            in >> flowlinearity_model;
            switch (flowlinearity_model)
            {
                case 0:  // k=f(x)
                    break;
                case 1:  // Read in alpha
                    // Alpha
                    in >> flowlinearity_model_values[0];
                    break;
                case 2:  // For equivalent flow in trianglular elements
                    // Alpha
                    in >> flowlinearity_model_values[0];
                    // Number of Fractures in Equivalent Medium
                    in >> flowlinearity_model_values[1];
                    // Reynolds Number above which non linear flow occurs.
                    in >> flowlinearity_model_values[2];
                    pcs_name_vector.push_back("PRESSURE1");
                    break;
                case 3:  // EE formulation for Forchheimer
                    in >> forchheimer_cf;
                    std::cout << "->Forchheimer nonlinear flow for EE"
                              << std::endl;
                    break;
                case 4:  // GW formulation for Forchheimer
                    in >> forchheimer_a1 >> forchheimer_a2;
                    std::cout
                        << "->Forchheimer nonlinear flow with given a1, a2"
                        << std::endl;
                    break;
                case 5:  // GW formulation for Forchheimer if a1=1/k0
                    in >> forchheimer_a2;
                    forchheimer_a1 = .0;
                    std::cout << "->Forchheimer nonlinear flow assuming a1=1/K0"
                              << std::endl;
                    break;
                case 6:  // EE formulation 2 for Forchheimer
                    in >> forchheimer_De;
                    std::cout << "->Forchheimer nonlinear flow for EE" << '\n';
                    break;
                default:
                    std::cout
                        << "Error in MMPRead: no valid flow linearity model"
                        << "\n";
                    break;
            }
            in.clear();
            continue;
        }
        //------------------------------------------------------------------------
        // 8..SORPTION_MODEL
        //------------------------------------------------------------------------
        // subkeyword found
        if (line_string.find("$ORGANIC_CARBON") != std::string::npos)
        {
            in.str(GetLineFromFile1(mmp_file));
            in >> foc;
            in.clear();
            continue;
        }

        //------------------------------------------------------------------------
        // 9..STORAGE
        //------------------------------------------------------------------------

        // subkeyword found
        if (line_string.find("$STORAGE") != std::string::npos)
        {
            in.str(GetLineFromFile1(mmp_file));
            in >> storage_model;
            switch (storage_model)
            {
                case 0:                             // S=f(x)
                    in >> storage_model_values[0];  // Function of pressure
                                                    // defined by curve
                    pcs_name_vector.push_back("PRESSURE1");
                    break;
                case 1:                             // S=const
                    in >> storage_model_values[0];  // Constant value in Pa
                    break;
                case 2:
                    in >> storage_model_values[0];  // S0
                    in >> storage_model_values[1];  // increase
                    in >> storage_model_values[2];  // sigma(z0)
                    in >> storage_model_values[3];  // d_sigma/d_z
                    break;
                case 3:
                    in >> storage_model_values[0];  // curve number (as real
                                                    // number)
                    in >> storage_model_values[1];  // sigma(z0)
                    in >> storage_model_values[2];  //_sigma/d_z
                    break;
                case 4:
                    in >> storage_model_values[0];  // curve number (as real
                                                    // number)
                    in >> storage_model_values[1];  // time collation
                    in >> storage_model_values[2];  // solid density
                    in >> storage_model_values[3];  // curve fitting factor,
                                                    // default 1
                    pcs_name_vector.push_back("PRESSURE1");
                    break;
                case 5:  // Storativity is a function of normal effective stress
                         // defined by curve, set up for KTB.
                    in >> storage_model_values[0];  // curve number
                    in >> storage_model_values[1];  // time collation
                    in >> storage_model_values[2];  // Default storage value for
                                                    // material groups > 0
                    in >>
                        storage_model_values[3];  // Angular difference between
                                                  // Y direction and Sigma 1
                    pcs_name_vector.push_back("PRESSURE1");
                    break;
                case 6:  // Storativity is a function of normal effective stress
                         // defined by curve and distance from borehole, set up
                         // for KTB.
                    in >> storage_model_values[0];  // curve number
                    in >> storage_model_values[1];  // time collation
                    in >> storage_model_values[2];  // Default storage value for
                                                    // material groups > 0
                    in >>
                        storage_model_values[3];  // Angular difference between
                                                  // Y direction and Sigma 1
                    in >> storage_model_values[4];  // Borehole (x) coordinate
                    in >> storage_model_values[5];  // Borehole (y) coordinate
                    in >> storage_model_values[6];  // Borehole (z) coordinate
                    in >> storage_model_values[7];  // Maximum thickness of
                                                    // shear zone
                    in >> storage_model_values[8];  // Fracture density
                    pcs_name_vector.push_back("PRESSURE1");
                    break;
                case 10:  // S=const for permeability_saturation_model = 10
                    storage_model_values[0] = 1;
                    break;
                case 11:
                    storage_model_values[0] = 1;
                    break;
                default:
                    cout << "Error in MMPRead: no valid storativity model"
                         << "\n";
                    break;
                case 7:  // RW/WW
                {
                    if (!PCSGet("DEFORMATION"))
                    {
                        ScreenMessage(
                            "Error: Porosity model 7 must be combined with "
                            "deformation process");
                        exit(EXIT_FAILURE);
                    }
                }
                break;
            }
            in.clear();
            continue;
        }
        //------------------------------------------------------------------------
        // 10..CONDUCTIVITY_MODEL
        //------------------------------------------------------------------------
        // subkeyword found
        if (line_string.find("$CONDUCTIVITY_MODEL") != std::string::npos)
        {
            in.str(GetLineFromFile1(mmp_file));
            in >> conductivity_model;
            switch (conductivity_model)
            {
                case 0:  // K=f(x)
                    break;
                case 1:  // K=const
                    in >> conductivity;
                    break;
                case 2:  // Manning
                    break;
                case 3:  // Chezy
                    break;
                default:
                    std::cout << "Error in MMPRead: no valid conductivity model"
                              << "\n";
                    break;
            }
            in.clear();
            continue;
        }

        // subkeyword found
        if (line_string.find("$UNCONFINED") != string::npos)
        {
            // unconfined_flow_group = 1;
            in.str(GetLineFromFile1(mmp_file));
            in >> unconfined_flow_group;  // 5.3.07 JOD
            in.clear();
            continue;
        }

        //------------------------------------------------------------------------
        // 11..PERMEABILITY_TENSOR
        //------------------------------------------------------------------------
        // subkeyword found
        if (line_string.find("$PERMEABILITY_TENSOR") != std::string::npos)
        {
            in.str(GetLineFromFile1(mmp_file));
            in >> permeability_tensor_type_name;
            switch (permeability_tensor_type_name[0])
            {
                case 'I':  // isotropic
                    permeability_tensor_type = 0;
                    permeability_model = 1;
                    in >> permeability_tensor[0];
                    // CMCD to pick up 2D and 3D Isotropic case.
                    permeability_tensor[1] = permeability_tensor[2] =
                        permeability_tensor[0];
                    break;
                case 'O':  // orthotropic
                    permeability_tensor_type = 1;
                    if (geo_dimension == 0)
                        std::cout << "Error in CMediumProperties::Read: no "
                                     "geometric dimension"
                                  << "\n";
                    if (geo_dimension == 2)
                    {
                        in >> permeability_tensor[0];
                        in >> permeability_tensor[1];
                    }
                    if (geo_dimension == 3)
                    {
                        in >> permeability_tensor[0];
                        in >> permeability_tensor[1];
                        in >> permeability_tensor[2];
                    }
                    break;
                case 'A':  // anisotropic
                    permeability_tensor_type = 2;
                    if (geo_dimension == 0)
                        std::cout << "Error in CMediumProperties::Read: no "
                                     "geometric dimension"
                                  << "\n";
                    if (geo_dimension == 2)
                    {
                        in >> permeability_tensor[0];
                        in >> permeability_tensor[1];
                        in >> permeability_tensor[2];
                        in >> permeability_tensor[3];
                    }
                    if (geo_dimension == 3)
                    {
                        in >> permeability_tensor[0];
                        in >> permeability_tensor[1];
                        in >> permeability_tensor[2];
                        in >> permeability_tensor[3];
                        in >> permeability_tensor[4];
                        in >> permeability_tensor[5];
                        in >> permeability_tensor[6];
                        in >> permeability_tensor[7];
                        in >> permeability_tensor[8];
                    }
                    break;
                case 'F':                    // SB: read from file
                    permeability_model = 2;  // OK
                    in >> permeability_file;
                    break;
                default:
                    std::cout
                        << "Error in MMPRead: no valid permeability tensor type"
                        << "\n";
                    break;
            }
            in.clear();
            continue;
        }

        //------------------------------------------------------------------------
        // 12. $PERMEABILITY_FUNCTION
        //				(i)		_DEFORMATION
        //				(ii)	_PRESSURE
        //				(iii)	_SATURATION
        //				(iv)	_STRESS
        //				(v)		_VELOCITY
        //				(vi)    _POROSITY
        //------------------------------------------------------------------------
        //------------------------------------------------------------------------
        // 12.1 PERMEABILITY_FUNCTION_DEFORMATION
        //------------------------------------------------------------------------
        // subkeyword found
        if (line_string.find("$PERMEABILITY_FUNCTION_DEFORMATION") !=
            std::string::npos)
        {
            in.str(GetLineFromFile1(mmp_file));
            in >> permeability_model;
            switch (permeability_model)
            {
                case 0:  // k=f(x)
                    break;
                case 1:  // k=const
                    in >> permeability;
                    break;
                default:
                    std::cout << "Error in MMPRead: no valid permeability model"
                              << "\n";
                    break;
            }
            in.clear();
            continue;
        }
        // WX: 05.2010
        if (line_string.find("$PERMEABILITY_FUNCTION_STRAIN") != string::npos)
        {
            in.str(GetLineFromFile1(mmp_file));
            in >> permeability_strain_model;
            switch (permeability_strain_model)
            {
                case 0:  // strain_volume
                    break;
                case 1:  // strain_volumeeff plas strain
                    in >> permeability_strain_model_value[0];
                    break;
                case 2:  // eff plas strainif eff plas strain>0, f(strainp).
                         // else stain volume
                    in >> permeability_strain_model_value[0];
                    break;
                case 3:  // if eff plas strain>0, f(strainp). else stain volume
                    in >> permeability_strain_model_value[0];  // for strain
                                                               // volume
                    in >> permeability_strain_model_value[1];  // for eff plas
                                                               // strain
                    break;
                case 4:  // strain volume + eff plas strain
                    in >> permeability_strain_model_value[0];
                    in >> permeability_strain_model_value[1];
                    break;
                case 5:  // strain volume (threshold value, can be changed with
                         // plas strain)
                    in >> permeability_strain_model_value[0];  // threshold vol.
                                                               // strain
                    in >> permeability_strain_model_value
                              [1];  // d_fac/d_volStrain
                                    // when vol. strain
                                    // <= threshold
                    in >> permeability_strain_model_value
                              [2];  // d_fac/d_volStrain
                                    // when vol. strain
                                    // > threshold
                    in >> permeability_strain_model_value
                              [3];  // curve numer for dependenc between
                                    // threshold and plas strain
                    // if -1, threshold is constant
                    in >> permeability_strain_model_value[4];  // lower limit
                    in >> permeability_strain_model_value[5];  // uper limit
                    break;
                default:
                    cout << "Error in MMPRead: no valid permeability strain "
                            "model"
                         << "\n";
                    break;
            }
            in.clear();
            continue;
        }
        //------------------------------------------------------------------------
        // 12.2 PERMEABILITY_FUNCTION_PRESSURE
        //------------------------------------------------------------------------
        // subkeyword found
        if (line_string.find("$PERMEABILITY_FUNCTION_PRESSURE") !=
            std::string::npos)
        {
            in.str(GetLineFromFile1(mmp_file));
            in >> permeability_pressure_model;
            switch (permeability_pressure_model)
            {
                case 0:  // k=f(x)
                    break;
                case 1:  // k=const
                    in >> permeability_pressure_model_values[0];
                    pcs_name_vector.push_back("PRESSURE1");
                    break;
                case 2:  // Permebility is a function of effective stress
                    in >> permeability_pressure_model_values[0];
                    in >> permeability_pressure_model_values[1];
                    in >> permeability_pressure_model_values[2];
                    in >> permeability_pressure_model_values[3];
                    pcs_name_vector.push_back("PRESSURE1");
                    break;
                case 3:  // Permeability is a function of non linear flow
                    in >> permeability_pressure_model_values[0];
                    in >> permeability_pressure_model_values[1];
                    in >> permeability_pressure_model_values[2];
                    in >> permeability_pressure_model_values[3];
                    pcs_name_vector.push_back("PRESSURE1");
                    break;
                case 4:  // Function of effective stress from a curve
                    in >> permeability_pressure_model_values[0];
                    in >> permeability_pressure_model_values[1];
                    in >> permeability_pressure_model_values[2];
                    pcs_name_vector.push_back("PRESSURE1");
                    break;
                case 5:  // Function of overburden converted to effective stress
                         // and related to a curve.
                    in >> permeability_pressure_model_values[0];
                    in >> permeability_pressure_model_values[1];
                    in >> permeability_pressure_model_values[2];
                    pcs_name_vector.push_back("PRESSURE1");
                    break;
                case 6:  // Normal effective stress calculated according to
                         // fracture orientation, related over a curve : KTB
                         // site
                    in >> permeability_pressure_model_values[0];
                    in >> permeability_pressure_model_values[1];
                    in >> permeability_pressure_model_values[2];
                    in >> permeability_pressure_model_values[3];
                    pcs_name_vector.push_back("PRESSURE1");
                    break;
                case 7:  // Normal effective stress calculated according to
                         // fracture orientation, related over a curve : KTB
                         // site, and distance to the borehole.
                    in >> permeability_pressure_model_values[0];
                    in >> permeability_pressure_model_values[1];
                    in >> permeability_pressure_model_values[2];
                    in >> permeability_pressure_model_values[3];
                    in >> permeability_pressure_model_values[4];
                    in >> permeability_pressure_model_values[5];
                    in >> permeability_pressure_model_values[6];
                    in >> permeability_pressure_model_values[7];
                    in >> permeability_pressure_model_values[8];
                    pcs_name_vector.push_back("PRESSURE1");
                    break;
                case 8:  // Effective stress related to depth and curve, Urach
                    in >> permeability_pressure_model_values[0];
                    in >> permeability_pressure_model_values[1];
                    pcs_name_vector.push_back("PRESSURE1");
                    break;
                case 9:  // Effective stress related to depth and curve, and to
                         // temperature, Urach
                    in >> permeability_pressure_model_values[0];
                    in >> permeability_pressure_model_values[1];
                    in >> permeability_pressure_model_values[2];
                    in >> permeability_pressure_model_values[3];
                    in >> permeability_pressure_model_values[4];
                    pcs_name_vector.push_back("PRESSURE1");
                    pcs_name_vector.push_back("TEMPERATURE1");
                    break;
                case 10:  // WX:05.2010 directly from curve
                    in >> permeability_pressure_model_values
                              [0];  // WX: curve number 05.2010
                    break;
                default:
                    std::cout << "Error in MMPRead: no valid permeability model"
                              << "\n";
                    break;
            }
            in.clear();
            continue;
        }
        //------------------------------------------------------------------------
        // 12.3 PERMEABILITY_FUNCTION_SATURATION
        //------------------------------------------------------------------------
        if (line_string.find("$PERMEABILITY_SATURATION") != std::string::npos)
        {
            int num_phases = 0;
            if (H3_Process)
                num_phases =
                    3;  // the existence of additional fluids in the .mfp file
                        // is not an indicator of a 4-phase system.
            else if (H2_Process)
                num_phases =
                    2;  // the existence of additional fluids in the .mfp file
                        // is not an indicator of a 3-phase system.
            else if (H_Process)
                num_phases = 1;
            //
            for (k = 0; k < num_phases; k++)
            {
                in.str(GetLineFromFile1(mmp_file));
                in >> permeability_saturation_model[k];
                //
                switch (permeability_saturation_model[k])
                {
                    case 0:                              // Curve of [Sw,k_rel]
                        in >> perm_saturation_value[k];  // curve number
                        break;
                    //
                    case 1:                              // Constant
                        in >> perm_saturation_value[k];  // constant value of
                                                         // k_rel
                        break;
                    //
                    case 2:  // krg = 1.0 - krl (only for gas phase)
                        if (k == 0)
                        {
                            ScreenMessage(
                                "ERROR in MMPRead: Relative permeability model "
                                "2 is only valid for the gas phase.\n");
                            exit(0);
                        }
                        in >> minimum_relative_permeability;
                        break;
                    //
                    case 3:  // JT: Function (power or linear)      WETTING
                        // krw = (m*Se)^exp
                        // Se  = (sl - slr) / (slm - slr)
                        in >> residual_saturation[k];  // slr: residual
                                                       // saturation, this phase
                        in >> maximum_saturation[k];   // slm: maximum
                                                       // saturation, this phase
                        in >> saturation_exponent[k];  // (set exponent = 1
                                                       // results in a linear
                                                       // k_rel function)
                        in >> perm_saturation_value[k];  // multiplier "m"
                        in >>
                            minimum_relative_permeability;  // minimum relative
                                                            // permeability this
                                                            // phase
                        break;
                    //
                    case 33:  // JT: Function (power or linear)     NON-WETTING
                        // krg = (m*(1-Se))^exp
                        // Se  = (sl - slr) / (slm - slr) --> slr = 1 - sgm -->
                        // slm = 1 - sgr
                        in >> residual_saturation[k];  // sgr: residual
                                                       // saturation, this phase
                        in >> maximum_saturation[k];   // sgm: maximum
                                                       // saturation, this phase
                        in >> saturation_exponent[k];  // (set exponent = 1
                                                       // results in a linear
                                                       // k_rel function)
                        in >> perm_saturation_value[k];  // multiplier "m"
                        in >>
                            minimum_relative_permeability;  // minimum relative
                                                            // permeability this
                                                            // phase
                        break;
                    //
                    case 4:  // 2-phase Van Genuchten/Mualem Model  WETTING
                        // krw = pow(se,0.5) *
                        // pow(1.0-pow(1.0-pow(se,1.0/m),m),2) Se  = (sl - slr)
                        // / (slm - slr)
                        in >> residual_saturation[k];  // slr: residual
                                                       // saturation, this phase
                        in >> maximum_saturation[k];   // slm: maximum
                                                       // saturation, this phase
                        in >> saturation_exponent[k];  // exponent (always
                                                       // <= 1.0) --> (typical
                                                       // is 0.5) i.e. n = 1 /
                                                       // (1 - exponent) == 2.0
                        in >>
                            minimum_relative_permeability;  // minimum relative
                                                            // permeability this
                                                            // phase
                        break;
                    //
                    case 44:  // 2-phase Van Genuchten/Mualem Model NON-WETTING
                        // Water Resour. Res. VOL. 23, pp2197-2206 1987. Air
                        in >> residual_saturation[k];  // sgr: residual
                                                       // saturation, this phase
                        in >> maximum_saturation[k];   // sgm: maximum
                                                       // saturation, this phase
                        in >> saturation_exponent[k];  // exponent (always
                                                       // <= 1.0) --> (typical
                                                       // is 0.5) i.e. n = 1 /
                                                       // (1 - exponent) == 2.0
                        in >>
                            minimum_relative_permeability;  // minimum relative
                                                            // permeability this
                                                            // phase
                        break;
                    //
                    case 6:  // Brooks/Corey						   WETTING
                        // krw = pow(se,3.0+2.0/m)
                        // Se  = (sl - slr) / (slm - slr)
                        in >> residual_saturation[k];  // slr: residual
                                                       // saturation, this phase
                        in >> maximum_saturation[k];   // slm: maximum
                                                       // saturation, this phase
                        in >>
                            saturation_exponent[k];  // exponent (always >= 1.0)
                                                     // (typical might be 2.0)
                        in >>
                            minimum_relative_permeability;  // minimum relative
                                                            // permeability this
                                                            // phase
                        break;
                    //
                    case 10:  // unconfined GW 6/2012 JOD
                        break;
                    case 61:  // Brooks-Corey CB modified Seff
                        in >> residual_saturation[k];  // sgr: residual
                                                       // saturation, this phase
                        in >> maximum_saturation[k];   // sgm: maximum
                                                       // saturation, this phase
                        in >>
                            saturation_exponent[k];  // exponent (always >= 1.0)
                                                     // (typical might be 2.0)
                        break;
                    case 66:  // Brooks/Corey					   NON-WETTING
                        // krg = pow(1.0-se,2)*(1.0-pow(se,1.0+2.0/m))
                        // Se  = (sl - slr) / (slm - slr) --> slr = 1 - sgm -->
                        // slm = 1 - sgr
                        in >> residual_saturation[k];  // sgr: residual
                                                       // saturation, this phase
                        in >> maximum_saturation[k];   // sgm: maximum
                                                       // saturation, this phase
                        in >>
                            saturation_exponent[k];  // exponent (always >= 1.0)
                                                     // (typical might be 2.0)
                        in >>
                            minimum_relative_permeability;  // minimum relative
                                                            // permeability this
                                                            // phase
                        break;
                    //
                    case 7:  // Corey's curves					   WETTING
                        // krw = pow(se,4)
                        // Se  = (sl - slr) / (slm - slr)
                        in >> residual_saturation[k];  // slr: residual
                                                       // saturation, this phase
                        in >> maximum_saturation[k];   // slm: maximum
                                                       // saturation, this phase
                        in >>
                            saturation_exponent[k];  // exponent (always >= 1.0)
                        in >>
                            minimum_relative_permeability;  // minimum relative
                                                            // permeability this
                                                            // phase
                        break;
                    //
                    case 77:  // Corey's curves					   NON-WETTING
                        // krg = pow((1.0-se),2)*(1.0-pow(se,2))
                        // Se  = (sl - slr) / (slm - slr) --> slr = 1 - sgmax
                        // --> slm = 1 - sgr
                        in >> residual_saturation[k];  // sgr: residual
                                                       // saturation, this phase
                        in >> maximum_saturation[k];   // sgm: maximum
                                                       // saturation, this phase
                        in >>
                            saturation_exponent[k];  // exponent (always >= 1.0)
                        in >>
                            minimum_relative_permeability;  // minimum relative
                                                            // permeability this
                                                            // phase
                        break;
                    //
                    default:
                    {
                        ScreenMessage(
                            "Error in MMPRead: no valid permeability "
                            "saturation model.\n");
                        abort();
                    }
                    break;
                }
                in.clear();
            }
            continue;
        }
        //------------------------------------------------------------------------
        // 12.4 PERMEABILITY_FUNCTION_STRESS
        //------------------------------------------------------------------------
        // subkeyword found
        if (line_string.find("$PERMEABILITY_FUNCTION_STRESS") !=
            std::string::npos)
        {
            in.str(GetLineFromFile1(mmp_file));
            in >> permeability_stress_mode;
            switch (permeability_stress_mode)
            {
                case 0:  // k=f(x)
                    break;
                case 1:  // k=const
                    in >> permeability;
                    break;
                case 2:  // Modified LBNL model. WW
                    c_coefficient = new double[18];
                    in >> c_coefficient[0]          // b_0
                        >> c_coefficient[1]         // alpha
                        >> c_coefficient[2]         // br_1
                        >> c_coefficient[3]         // br_2
                        >> c_coefficient[4]         // br_3
                        >> c_coefficient[5] >> ws;  // fracture frequency
                    for (i = 6; i < 18; i++)
                        c_coefficient[i] = 0.0;
                    for (i = 0; i < 3; i++)
                    {
                        in.clear();
                        in.str(GetLineFromFile1(mmp_file));
                        in >> m_string;
                        if (m_string.find("_XX") != string::npos)
                            k = 0;
                        else if (m_string.find("_YY") != string::npos)
                            k = 1;
                        else if (m_string.find("_ZZ") != string::npos)
                            k = 2;
                        m_string.clear();
                        in >> m_string;
                        pch = strtok(const_cast<char*>(m_string.c_str()), seps);
                        buff << pch;
                        buff >> c_coefficient[6 + k * 4];
                        buff.clear();
                        while (pch != NULL)
                        {
                            pch = strtok(NULL, seps);
                            if (pch == NULL)
                                break;
                            string token = pch;
                            tokens.push_back(token);
                        }
                        for (j = 0; j < (int)tokens.size(); j++)
                        {
                            pch = strtok(const_cast<char*>(tokens[j].c_str()),
                                         seps1);
                            buff << pch;
                            buff >> f_buff;
                            buff.clear();
                            pch = strtok(NULL, seps1);
                            switch (pch[0])
                            {
                                case 'x':
                                    c_coefficient[k * 4 + 7] = f_buff;
                                    break;
                                case 'y':
                                    c_coefficient[k * 4 + 8] = f_buff;
                                    break;
                                case 'z':
                                    c_coefficient[k * 4 + 9] = f_buff;
                                    break;
                            }
                        }
                        tokens.clear();
                    }
                    break;
                case 3:  // Barton-Bandis  WW
                    c_coefficient = new double[24];
                    in >> c_coefficient[0]          // JRC
                        >> c_coefficient[1]         // JCS        //an0
                        >> c_coefficient[2]         // UCS        //Kn
                        >> c_coefficient[7]         // sig_h
                        >> c_coefficient[8] >> ws;  // sig_H
                    for (i = 3; i < 7; i++)
                        c_coefficient[i] = 0.0;
                    for (i = 9; i < 24; i++)
                        c_coefficient[i] = 0.0;
                    for (i = 0; i < 3; i++)
                    {
                        in.clear();
                        in.str(GetLineFromFile1(mmp_file));
                        in >> m_string;
                        if (m_string.find("_XX") != string::npos)
                            k = 0;
                        else if (m_string.find("_YY") != string::npos)
                            k = 1;
                        else if (m_string.find("_ZZ") != string::npos)
                            k = 2;
                        m_string.clear();
                        in >> m_string;
                        pch = strtok(const_cast<char*>(m_string.c_str()), seps);
                        buff << pch;
                        buff >> c_coefficient[9 + k * 4];
                        buff.clear();
                        while (pch != NULL)
                        {
                            pch = strtok(NULL, seps);
                            if (pch == NULL)
                                break;
                            string token = pch;
                            tokens.push_back(token);
                        }
                        for (j = 0; j < (int)tokens.size(); j++)
                        {
                            pch = strtok(const_cast<char*>(tokens[j].c_str()),
                                         seps1);
                            buff << pch;
                            buff >> f_buff;
                            buff.clear();
                            pch = strtok(NULL, seps1);
                            switch (pch[0])
                            {
                                case 'x':
                                    c_coefficient[k * 4 + 10] = f_buff;
                                    break;
                                case 'y':
                                    c_coefficient[k * 4 + 11] = f_buff;
                                    break;
                                case 'z':
                                    c_coefficient[k * 4 + 12] = f_buff;
                                    break;
                            }
                        }
                        tokens.clear();
                    }
                    //
                    CalStressPermeabilityFactor3_Coef();
                    //
                    break;

                default:
                    std::cout << "Error in MMPRead: no valid permeability model"
                              << "\n";
                    break;
            }
            in.clear();
            continue;
        }
        // AS,WX:08.2012
        if (line_string.find("$PERMEABILITY_FUNCTION_EFFSTRESS") !=
            std::string::npos)
        {
            in.str(GetLineFromFile1(mmp_file));
            in >> permeability_effstress_model;
            switch (permeability_effstress_model)
            {
                case 1:
                    in >> permeability_effstress_model_value[0];
                    break;
                default:
                    cout << "Error in MMPRead: no valid permeability stress "
                            "model"
                         << "\n";
                    break;
            }
            in.clear();
            continue;
        }
        //------------------------------------------------------------------------
        // 12.5 PERMEABILITY_FUNCTION_VELOCITY
        //------------------------------------------------------------------------
        // WW
        if (line_string.find("$PERMEABILITY_FUNCTION_VELOCITY") !=
            std::string::npos)
        {
            // WW
            // if(line_string.find("$PERMEABILITY_FUNCTION_STRESS")!=string::npos)
            // { //subkeyword found
            in.str(GetLineFromFile1(mmp_file));
            in >> permeability_model;
            switch (permeability_model)
            {
                case 0:  // k=f(x)
                    break;
                case 1:  // k=const
                    in >> permeability;
                    break;
                default:
                    std::cout << "Error in MMPRead: no valid permeability model"
                              << "\n";
                    break;
            }
            in.clear();
            continue;
        }
        //------------------------------------------------------------------------
        // 12.6 PERMEABILITY_FUNCTION_POROSITY
        //------------------------------------------------------------------------

        // subkeyword found
        if (line_string.find("$PERMEABILITY_FUNCTION_POROSITY") !=
            std::string::npos)
        {
            in.str(GetLineFromFile1(mmp_file));
            in >> permeability_porosity_model;
            switch (permeability_porosity_model)
            {
                case 0:  // k=f(x)
                    break;
                case 1:  // k=const
                    in >> permeability_porosity_model_values[0];
                    mmp_file->ignore(MAX_ZEILE, '\n');
                    break;
                case 2:  // Model from Ming Lian
                    in >> permeability_porosity_model_values[0];
                    in >> permeability_porosity_model_values[1];
                    break;
                case 3:  // HS: 11.2008, Kozeny-Carman relationship
                    // we set the tensor type first to isotropic and constant as
                    // initial value
                    permeability_tensor_type = 0;
                    permeability_model = 3;  // this means permeability depends
                                             // on K-C relationship
                    // initial values
                    in >> permeability_porosity_model_values[0];  // initial
                                                                  // porosity
                    in >>
                        permeability_porosity_model_values[1];  // initial
                                                                // permeability
                    KC_permeability_initial =
                        permeability_porosity_model_values[1];
                    KC_porosity_initial = permeability_porosity_model_values[0];

                    break;
                case 4:  // HS: 11.2008, Kozeny-Carman_normalized relationship
                    // we set the tensor type first to isotropic and constant as
                    // initial value
                    permeability_tensor_type = 0;
                    permeability_model = 4;  // this means permeability depends
                                             // on K-C_normalized relationship
                    // initial values
                    in >> permeability_porosity_model_values[0];  // initial
                                                                  // porosity
                    in >>
                        permeability_porosity_model_values[1];  // initial
                                                                // permeab ility
                    KC_permeability_initial =
                        permeability_porosity_model_values[1];
                    KC_porosity_initial = permeability_porosity_model_values[0];
                    break;
                case 5:  // HS: 01.2010, Clement 1996 model
                    // M. Thullner et al. 2004, J Contaminant Hydrology 70:
                    // 37-62, pp42
                    permeability_tensor_type = 0;
                    permeability_model = 5;  // Clement original model
                    // this is initial porosity
                    in >> permeability_porosity_model_values[0];
                    // this is initial permeability
                    in >> permeability_porosity_model_values[1];
                    break;
                case 6:  // HS: 01.2010, ,modified Clement, biomass colonies
                         // clogging
                    // M. Thullner et al. 2004, J Contaminant Hydrology 70:
                    // 37-62, pp42
                    permeability_tensor_type = 0;
                    permeability_model =
                        6;  // modified Clement, biomass growing in colonies
                    // this is initial porosity
                    in >> permeability_porosity_model_values[0];
                    // this is initial permeability
                    in >> permeability_porosity_model_values[1];
                    // this is parameter a
                    in >> permeability_porosity_model_values[2];
                    break;
                case 7:  // HS: 01.2010, ,modified Clement, biofilm clogging
                    // M. Thullner et al. 2004, J Contaminant Hydrology 70:
                    // 37-62, pp42
                    permeability_tensor_type = 0;
                    permeability_model =
                        7;  // modified Clement, biomass growing in biofilm
                    // this is initial porosity
                    in >> permeability_porosity_model_values[0];
                    // this is initial permeability
                    in >> permeability_porosity_model_values[1];
                    // this is parameter b
                    in >> permeability_porosity_model_values[2];
                    // this is prarameter k_fmin
                    in >> permeability_porosity_model_values[3];
                    break;
                case 8:  // AB: 2010, Updating permeability due to change in
                         // porosity
                    // Inital permeability provided in $PERMEABILITY_TENSOR will
                    // be read for initial condition
                    permeability_model = 8;
                    in >> permeability_porosity_updating_type_name;
                    switch (permeability_porosity_updating_type_name[0])
                    {
                        case 'K':  // Kozeny-Carman formulation
                            permeability_porosity_updating_type = 0;
                            break;
                        case 'V':  // Verma-Pruess formulation
                            permeability_porosity_updating_type = 1;
                            in >> permeability_porosity_updating_values
                                      [0];  // critical porosity
                            in >> permeability_porosity_updating_values
                                      [1];  // n: a power law exponent
                            break;
                    }
                    if (this->permeability_tensor_type == 2)
                        std::cout << "Warning in MMPRead: permeability_model 8 "
                                     "is not implemented for "
                                     "permeability_tensor_type 2"
                                  << "\n";
                    break;
                default:
                    std::cout << "Error in MMPRead: no valid permeability model"
                              << "\n";
                    break;
            }
            in.clear();
            continue;
        }

        //....................................................................
        // subkeyword found
        if (line_string.find("$CAPILLARY_PRESSURE") != std::string::npos)
        {
            in.str(GetLineFromFile1(mmp_file));
            in >> capillary_pressure_model;
            i = 0;
            // Value needed to check for old version formats.
            capillary_pressure_values[2] = -1.0;
            bool old_format = false;
            //
            switch (capillary_pressure_model)
            {
                case 0:                                  // k=f(Se)
                    in >> capillary_pressure_values[0];  // curve
                    in >> capillary_pressure_values[1];  // Slr
                    in >> capillary_pressure_values[2];  // Slmax
                    //
                    // JT: Check for old version format.
                    if (capillary_pressure_values[2] < 0.0)
                    {
                        capillary_pressure_values[1] =
                            residual_saturation[0];  // old version uses
                                                     // relative permeabilty
                                                     // values for this
                        capillary_pressure_values[2] =
                            maximum_saturation[0];  // old version uses relative
                                                    // permeabilty values for
                                                    // this
                        old_format = true;
                    }
                    break;
                case 1:  // const
                    in >>
                        capillary_pressure_values[0];  // the constant Pc value
                    //
                    // This next value is NOT REQUIRED. It is here because in
                    // the old version this entry was used to represent a
                    // constant value of saturation in PP models. SHOULD DELETE
                    // THIS CHECK EVENTUALLY!
                    in >> capillary_pressure_values[2];
                    if (capillary_pressure_values[2] >= 0.0)
                    {  // Then a constant saturation value has been entered.
                       // This is model #2.
                        ScreenMessage(
                            "WARNING in MMPRead. Capillary pressure model 1 "
                            "used for a constant saturation. THIS IS "
                            "NOW MODEL #2. PLEASE SWITCH TO MODEL #2.\n");
                        capillary_pressure_model = 2;
                        capillary_pressure_values[0] =
                            capillary_pressure_values[2];
                    }
                    //
                    // Assign bounds
                    capillary_pressure_values[1] = 0.0;  // Slr
                    capillary_pressure_values[2] = 1.0;  // Slmax
                    break;
                case 2:  // Constant saturation for pp models (for WX, from JT)
                    in >> capillary_pressure_values[0];  // The fixed saturation
                    // Assign bounds
                    capillary_pressure_values[1] = 0.0;  // Slr
                    capillary_pressure_values[2] = 1.0;  // Slmax
                    break;
                case 4:                                  // van Genuchten
                    in >> capillary_pressure_values[0];  // Pb (or "alpha" if
                                                         // [alpha_switch>0])
                    in >> capillary_pressure_values[1];  // Slr
                    in >> capillary_pressure_values[2];  // Slmax
                    in >> capillary_pressure_values
                              [3];  // exponent (always <= 1.0) --> (typical is
                                    // 0.5) i.e. n = 1 / (1
                                    // - exponent) == 2.0
                    in >> capillary_pressure_values[4];  // maximum Pc
                    in >> i;  // alpha_switch (default = 0)
                    if (i > 0)
                        entry_pressure_conversion = true;
                    //
                    // JT: Check for old version format.
                    if (capillary_pressure_values[2] < 0.0)
                    {
                        entry_pressure_conversion =
                            true;  // entry is alpha in old version
                        capillary_pressure_values[1] =
                            residual_saturation[0];  // old version uses
                                                     // relative permeabilty
                                                     // values for this
                        capillary_pressure_values[2] =
                            maximum_saturation[0];  // old version uses relative
                                                    // permeabilty values for
                                                    // this
                        capillary_pressure_values[3] =
                            saturation_exponent[0];  // old version uses
                                                     // relative permeabilty
                                                     // values for this
                        capillary_pressure_values[4] = 1.0e10;
                        old_format = true;
                    }
                    break;
                case 6:  //  Brook & Corey. 2.05.2008. WW
                    in >> capillary_pressure_values[0];  // Pb
                    in >> capillary_pressure_values[1];  // Slr
                    in >> capillary_pressure_values[2];  // Slmax
                    in >> capillary_pressure_values[3];  // exponent (always
                                                         // >= 1.0) (typical
                                                         // might be 2.0)
                    in >> capillary_pressure_values[4];  // maximum Pc
                    //
                    // JT: Check for old version format.
                    if (capillary_pressure_values[2] < 0.0)
                    {
                        capillary_pressure_values[1] =
                            residual_saturation[0];  // old version uses
                                                     // relative permeabilty
                                                     // values for this
                        capillary_pressure_values[2] =
                            maximum_saturation[0];  // old version uses relative
                                                    // permeabilty values for
                                                    // this
                        capillary_pressure_values[3] =
                            saturation_exponent[0];  // old version uses
                                                     // relative permeabilty
                                                     // values for this
                        capillary_pressure_values[4] = 1.0e10;
                        old_format = true;
                    }
                    break;
                case 10:  // unconfined 3D GW. 5.3.07 JOD
                    in >> capillary_pressure_values[1];  // Slr
                    in >> capillary_pressure_values[0];  // Pb
                    break;
                default:
                    ScreenMessage(
                        "Error in MMPRead: no valid capillary pressure "
                        "model.\n");
                    exit(1);
                    break;
            }
            if (old_format)
            {
                ScreenMessage(
                    "\n--\n Adopting capillary pressure saturation parameters "
                    "from the\n");
                ScreenMessage(
                    " relative permeability function for phase 0. "
                    "Alternatively, you\n");
                ScreenMessage(
                    " may enter capillary pressure specific parameters "
                    "directly.\n--/n");
            }
            in.clear();
            continue;
        }
        //....................................................................
        // Dual Richards
        if (line_string.find("$TRANSFER_COEFFICIENT") != std::string::npos)
        {
            in.str(GetLineFromFile1(mmp_file));
            in >> transfer_coefficient;  //(-)
            //   in >> unsaturated_hydraulic_conductivity;      //(L/T)
            in.clear();
            continue;
        }
        //....................................................................
        // Dual Richards
        if (line_string.find("$SPECIFIC_STORAGE") != std::string::npos)
        {
            in.str(GetLineFromFile1(mmp_file));
            in >> specific_storage;  //(Pa-1)
            in.clear();
            continue;
        }
        // AS:08.2012 storage function eff stress
        if (line_string.find("$STORAGE_FUNCTION_EFFSTRESS") !=
            std::string::npos)
        {
            in.str(GetLineFromFile1(mmp_file));
            in >> storage_effstress_model;
            switch (storage_effstress_model)
            {
                case 1:
                    in >> storage_effstress_model_value[0];
                    break;
                default:
                    cout << "Error in MMPRead: no valid storage stress model"
                         << "\n";
                    break;
            }
            in.clear();
            continue;
        }
        //------------------------------------------------------------------------
        // 14 MASSDISPERSION_
        //			(1) LONGITUDINAL
        //			(2) TRANSVERSE
        //------------------------------------------------------------------------
        // subkeyword found
        if (line_string.find("$MASS_DISPERSION") != std::string::npos)
        {
            in.str(GetLineFromFile1(mmp_file));
            in >> mass_dispersion_model;
            switch (mass_dispersion_model)
            {
                case 0:  // f(x)
                    break;
                case 1:  // Constant value
                    in >> mass_dispersion_longitudinal;
                    in >> mass_dispersion_transverse;
                    in >> lgpn;
                    if (lgpn > 0)
                        cout << "      Limiting Grid Peclet Numbers to " << lgpn
                             << "\n";
                    break;
                default:
                    std::cout << "Error in CMediumProperties::Read: no valid "
                                 "mass dispersion model"
                              << "\n";
                    break;
            }
            in.clear();
            continue;
        }
        // Chiogna et al ES&T model: Daq-dependent alpha_t
        if (line_string.find("$COMPOUND_DEPENDENT_DT") != std::string::npos)
        {  // subkeyword found
            in.str(GetLineFromFile1(mmp_file));
            in >> alpha_t_model >> graindiameter >> hydraulicrad >> betaexpo;
            in.clear();
            continue;
        }
        //------------------------------------------------------------------------
        // 15 HEATDISPERSION
        //			(1) LONGTIDUINAL
        //			(2) TRANSVERSE
        //------------------------------------------------------------------------
        // subkeyword found
        if (line_string.find("$HEAT_DISPERSION") != std::string::npos)
        {
            in.str(GetLineFromFile1(mmp_file));
            in >> heat_dispersion_model;
            switch (heat_dispersion_model)
            {
                case 0:  // f(x)
                    break;
                case 1:  // Constant value
                    in >> heat_dispersion_longitudinal;
                    in >> heat_dispersion_transverse;
                    break;
                default:
                    std::cout << "Error in CMediumProperties::Read: no valid "
                                 "heat dispersion model"
                              << "\n";
                    break;
            }
            in.clear();
            continue;
        }
        // subkeyword found
        if (line_string.find("$DIFFUSION") != std::string::npos)
        {
            in.str(GetLineFromFile1(mmp_file));
            in >> heat_diffusion_model;
            if (heat_diffusion_model == 1)
                in >> base_heat_diffusion_coefficient;
            else if (heat_diffusion_model == 273)  // JM allow old model number
                heat_diffusion_model = 1;
            in.clear();
            continue;
        }
        // subkeyword found
        if (line_string.find("$EVAPORATION") != string::npos)
        {
            in.str(GetLineFromFile1(mmp_file));
            in >> evaporation;
            in >> heatflux;
            in >> vaporfraction;
            in.clear();
            continue;
        }
        //------------------------------------------------------------------------
        // 16. Surface water
        //------------------------------------------------------------------------
        // subkeyword found
        if (line_string.find("$SURFACE_FRICTION") != std::string::npos)
        {
            in.str(GetLineFromFile1(mmp_file));
            in >> friction_coefficient >> friction_exp_slope >>
                friction_exp_depth;
            in.clear();
            continue;
        }

        // subkeyword found
        if (line_string.find("$WIDTH") != std::string::npos)
        {
            in.str(GetLineFromFile1(mmp_file));
            in >> overland_width;
            in.clear();
            continue;
        }

        // subkeyword found
        if (line_string.find("$RILL") != std::string::npos)
        {
            in.str(GetLineFromFile1(mmp_file));
            in >> rill_height >> rill_epsilon;
            in.clear();
            continue;
        }

        // subkeyword found
        if (line_string.find("$CHANNEL") != std::string::npos)
        {
            channel = 1;
            continue;
        }

        //------------------------------------------------------------------------
        // 19 ELECTRIC_CONDUCTIVITY
        //------------------------------------------------------------------------
        //------------------------------------------------------------------------
        // 20 UNCONFINED_FLOW_GROUP
        //------------------------------------------------------------------------
        //------------------------------------------------------------------------
        // 21 FLUID_EXCHANGE_WITH_OTHER_CONTINUA
        //------------------------------------------------------------------------

        //------------------------------------------------------------------------
        // 11..PERMEABILITY_DISTRIBUTION
        //------------------------------------------------------------------------
        // subkeyword found
        if (line_string.find("$PERMEABILITY_DISTRIBUTION") != std::string::npos)
        {
            in.str(GetLineFromFile1(mmp_file));
            in >> permeability_file;
            permeability_file = FilePath + permeability_file;
            //--------------------------------------
            // WW
            std::ifstream mmp_file(permeability_file.data(), std::ios::in);
            if (!mmp_file.good())
                std::cout << "Fatal error in MMPRead: no "
                             "PERMEABILITY_DISTRIBUTION file"
                          << "\n";
            mmp_file.close();
            permeability_model = 2;
            in.clear();
            continue;
        }

        //------------------------------------------------------------------------
        // 11..POROSITY_DISTRIBUTION
        //------------------------------------------------------------------------
        // subkeyword found
        if (line_string.find("$POROSITY_DISTRIBUTION") != std::string::npos)
        {
            in.str(GetLineFromFile1(mmp_file));
            in >> porosity_file;
            porosity_file = FilePath + porosity_file;
            // WW
            ifstream mmp_file(porosity_file.data(), ios::in);
            if (!mmp_file.good())
                std::cout
                    << "Fatal error in MMPRead: no POROSITY_DISTRIBUTION file"
                    << "\n";
            mmp_file.close();
            porosity_model = 11;
            in.clear();
            continue;
        }

        //------------------------------------------------------------------------
        // HEAT_TRANSFER
        //------------------------------------------------------------------------
        if (line_string.find("$HEAT_TRANSFER") != std::string::npos)  // NW
        {
            in.str(GetLineFromFile1(mmp_file));
            in >> heat_transfer_model;
            switch (heat_transfer_model)
            {
                case 0:  // f(x)
                    break;
                case 1:  // Constant value
                    in >> heat_transfer_model_value;
                    break;
                case 2:  // Effective heat transfer coefficient between fluid
                         // and particles
                    in >> effective_heat_transfer_model;
                    if (effective_heat_transfer_model == 1)
                    {
                        std::cout
                            << "-> Heat transfer model: Schaube11 with given h"
                            << '\n';
                        in >> heat_transfer_model_value;
                    }
                    else if (effective_heat_transfer_model == 2)
                    {
                        std::cout << "-> Heat transfer model: Schaube11 with h "
                                     "by Gnielinskl"
                                  << '\n';
                    }
                    else
                    {
                        std::cout << "Error in CMediumProperties::Read: no "
                                     "valid effective heat transfer model"
                                  << '\n';
                    }
                    break;
                default:
                    std::cout << "Error in CMediumProperties::Read: no valid "
                                 "heat transfer model"
                              << '\n';
                    break;
            }
            in.clear();
            continue;
        }

        //------------------------------------------------------------------------
        // Particle diameter
        //------------------------------------------------------------------------
        if (line_string.find("$PARTICLE_DIAMETER") != std::string::npos)  // NW
        {
            in.str(GetLineFromFile1(mmp_file));
            in >> particle_diameter_model;
            switch (particle_diameter_model)
            {
                case 0:  // f(x)
                    break;
                case 1:  // Constant value
                    std::cout << "-> Particle dimaeter is set" << '\n';
                    in >> particle_diameter_model_value;
                    break;
                default:
                    std::cout << "Error in CMediumProperties::Read: no valid "
                                 "heat transfer model"
                              << '\n';
                    break;
            }
            in.clear();
            continue;
        }

        if (line_string.find("$INTERPHASE_FRICTION") != std::string::npos)
        {
            in.str(GetLineFromFile1(mmp_file));
            in >> PhaseHeatedByFriction;
            this->setFrictionPhase(
                FiniteElement::convertFrictionPhase(PhaseHeatedByFriction));
            if (PhaseHeatedByFriction.compare("SOLID") != 0 &&
                PhaseHeatedByFriction.compare("FLUID") != 0 &&
                PhaseHeatedByFriction.compare("NONE") != 0)
                std::cout
                    << "Error in CMediumProperties::Read: $INTERPHASE_FRICTION "
                       "must be either SOLID or FLUID or NONE";
            in.clear();
            continue;
        }
        //------------------------------------------------------------------------
        // Element volme multiplyer
        //------------------------------------------------------------------------
        if (line_string.find("$ELEMENT_VOLUME_MULTIPLYER") != std::string::npos)
        {
            in.str(GetLineFromFile1(mmp_file));
            in >> ElementVolumeMultiplyer;
            std::cout << " Setting ElementVolumeMultiplyer to "
                      << ElementVolumeMultiplyer << "- times the grid value \n";
            in.clear();
            continue;
        }
    }
    return position;
}

/**************************************************************************
   FEMLib-Method: MMPWrite
   Task: master write functionn
   Programing:
   02/2004 OK Implementation
   last modification:
**************************************************************************/
void MMPWrite(std::string file_base_name)
// void MATWriteMediumProperties(fstream *mp_file)
{
    CMediumProperties* m_mat_mp = NULL;
    std::string sub_line;
    std::string line_string;
    //========================================================================
    // File handling
    std::string mp_file_name = file_base_name + MMP_FILE_EXTENSION;
    std::fstream mp_file(mp_file_name.c_str(), std::ios::trunc | ios::out);
    mp_file.setf(std::ios::scientific, std::ios::floatfield);
    mp_file.precision(12);
    if (!mp_file.good())
        return;
    mp_file << "GeoSys-MMP: Material Medium Properties "
               "------------------------------------------------\n";
    //========================================================================
    // MAT-MP list
    int no_mat = (int)mmp_vector.size();
    int i;
    for (i = 0; i < no_mat; i++)
    {
        m_mat_mp = mmp_vector[i];
        m_mat_mp->Write(&mp_file);
    }
    mp_file << "#STOP";
    mp_file.close();
}

/**************************************************************************
   FEMLib-Method: CMediumProperties::Write
   Task: read functionn
   Programing:
   02/2004 OK Template
   01/2005 OK new formats
   05/2005 OK CAPILLARY_PRESSURE,
   05/2005 OK Surface flow parameter
   06/2005 MB file names
   10/2005 OK/YD bugfix model 4
   10/2005 OK geo_name_vector
   ToDo: MB CONDUCTIVITY_MODEL, PERMEABILITY_SATURATION
**************************************************************************/
void CMediumProperties::Write(std::fstream* mmp_file)
{
    // KEYWORD
    *mmp_file << "#MEDIUM_PROPERTIES\n";
    //--------------------------------------------------------------------
    // NAME
    if (name.length() > 0)
    {
        *mmp_file << " $NAME"
                  << "\n";
        *mmp_file << "  ";
        *mmp_file << name << "\n";
    }
    //--------------------------------------------------------------------
    // GEO_TYPE
    //	if(geo_type_name.compare("DOMAIN") == 0) //OK
    if (_geo_type == GEOLIB::GEODOMAIN)  // OK / TF
    {
        *mmp_file << " $GEO_TYPE"
                  << "\n";
        *mmp_file << "  DOMAIN"
                  << "\n";
    }
    else if ((int)geo_name_vector.size() > 0)  // OK
    {
        *mmp_file << " $GEO_TYPE"
                  << "\n";
        for (int i = 0; i < (int)geo_name_vector.size(); i++)
        {
            *mmp_file << "  ";
            //			*mmp_file << geo_type_name;
            *mmp_file << GEOLIB::convertGeoTypeToString(_geo_type);
            *mmp_file << " ";
            *mmp_file << geo_name_vector[i] << "\n";
        }
    }
    //--------------------------------------------------------------------
    // DIMENSION
    *mmp_file << " $GEOMETRY_DIMENSION"
              << "\n";
    *mmp_file << "  ";
    *mmp_file << geo_dimension << "\n";
    *mmp_file << " $GEOMETRY_AREA"
              << "\n";
    *mmp_file << "  ";
    *mmp_file << geo_area << "\n";
    //--------------------------------------------------------------------
    // PROPERTIES
    //....................................................................
    // POROSITY
    if (porosity_model > -1)
    {
        *mmp_file << " $POROSITY"
                  << "\n";
        *mmp_file << "  ";
        *mmp_file << porosity_model << " ";
        switch (porosity_model)
        {
            case 1:
                *mmp_file << porosity_model_values[0] << "\n";
                break;
            case 11:                                 // MB ToDo
                *mmp_file << porosity_file << "\n";  // MB
                break;
        }
    }
    //....................................................................
    // TORTUOSITY //OK4104
    if (tortuosity_model > -1)
    {
        *mmp_file << " $TORTUOSITY"
                  << "\n";
        *mmp_file << "  ";
        *mmp_file << tortuosity_model << " ";
        switch (tortuosity_model)
        {
            case 1:
                *mmp_file << tortuosity_model_values[0] << "\n";
                break;
            case 2:
                *mmp_file << tortuosity_model_values[0] << "\n";
                break;
        }
    }
    //....................................................................
    // CONDUCTIVITY //MB ToDo
    if (conductivity_model > -1)
    {
        *mmp_file << " $CONDUCTIVITY_MODEL"
                  << "\n";
        *mmp_file << "  ";
        *mmp_file << conductivity_model << " ";
        switch (conductivity_model)
        {
            case 1:
                *mmp_file << conductivity;
                break;
        }
        *mmp_file << "\n";
    }
    //....................................................................
    // STORAGE
    if (storage_model > -1)
    {
        *mmp_file << " $STORAGE"
                  << "\n";
        switch (storage_model)
        {
            case 1:
                *mmp_file << "  " << storage_model << " "
                          << storage_model_values[0] << "\n";
                break;
        }
    }
    //....................................................................
    // PERMEABILITY_TENSOR
    if (permeability_tensor_type > -1)
    {
        *mmp_file << " $PERMEABILITY_TENSOR"
                  << "\n";
        switch (permeability_tensor_type)
        {
            case 0:
                *mmp_file << "  "
                          << "ISOTROPIC"
                          << " " << permeability_tensor[0] << "\n";
                break;
            case 3:  // MB for heterogeneous fields //OK
                *mmp_file << "  "
                          << "FILE"
                          << " " << permeability_file << "\n";
                break;
        }
    }
    //....................................................................
    // PERMEABILITY_DISTRIBUTION
    if (permeability_model == 2)
    {
        *mmp_file << " $PERMEABILITY_DISTRIBUTION"
                  << "\n";
        *mmp_file << "  " << permeability_file << "\n";
    }
    //....................................................................
    // PERMEABILITY
    if (permeability_saturation_model[0] > -1)
    {
        *mmp_file << " $PERMEABILITY_SATURATION"
                  << "\n";
        for (int i = 0; i < (int)mfp_vector.size(); i++)
        {
            *mmp_file << "  ";
            *mmp_file << permeability_saturation_model[0] << " ";
            switch (permeability_saturation_model[0])
            {
                case 0:
                    *mmp_file << (int)perm_saturation_value[0] << "\n";
                    break;
                case 4:  // VG
                    *mmp_file << perm_saturation_value[0] << " ";
                    *mmp_file << 1.0 - perm_saturation_value[0] << " ";
                    *mmp_file << perm_saturation_value[0] << "\n";
                    break;
            }
        }
    }
    //....................................................................
    // CAPILLARY PRESSURE
    if (capillary_pressure_model > -1)
    {
        *mmp_file << " $CAPILLARY_PRESSURE"
                  << "\n";
        *mmp_file << "  ";
        *mmp_file << capillary_pressure_model << " ";
        switch (capillary_pressure_model)
        {
            case 0:
                *mmp_file << (int)capillary_pressure_values[0] << "\n";
                break;
            case 4:  // VG
                *mmp_file << capillary_pressure_values[0] << "\n";
                break;
            case 9:  // HydroSphere
                *mmp_file << capillary_pressure_values[0] << " ";
                *mmp_file << capillary_pressure_values[1] << "\n";
                *mmp_file << capillary_pressure_values[2] << "\n";
                break;
        }
    }
    //....................................................................
    // HEAT DISPERSION
    if (heat_dispersion_model > -1)
    {
        *mmp_file << " $HEAT_DISPERSION"
                  << "\n";
        switch (heat_dispersion_model)
        {
            case 1:
                // CMCD permeability
                *mmp_file << "  " << heat_dispersion_model << " "
                          << heat_dispersion_longitudinal << " "
                          << heat_dispersion_transverse << "\n";
                break;
        }
    }
    //....................................................................
    // MASS DISPERSION
    if (mass_dispersion_model > -1)
    {
        *mmp_file << " $MASS_DISPERSION"
                  << "\n";
        switch (mass_dispersion_model)
        {
            case 1:
                // CMCD permeability
                *mmp_file << "  " << mass_dispersion_model << " "
                          << mass_dispersion_longitudinal << " "
                          << mass_dispersion_transverse << "\n";
                break;
        }
    }

    //----------------------------------------------------------------------
}

////////////////////////////////////////////////////////////////////////////
// Access functions
////////////////////////////////////////////////////////////////////////////

/**************************************************************************
   FEMLib-Method: CMediumProperties
   Task: get instance by name
   Programing:
   02/2004 OK Implementation
   last modification:
**************************************************************************/
CMediumProperties* MMPGet(const std::string& mat_name)
{
    CMediumProperties* m_mat = NULL;
    int no_mat = (int)mmp_vector.size();
    int i;
    for (i = 0; i < no_mat; i++)
    {
        m_mat = mmp_vector[i];
        if (mat_name.compare(m_mat->name) == 0)
            return m_mat;
    }
    return NULL;
}

/**************************************************************************
   FEMLib-Method: CMediumProperties
   Task: get instance by name
   Programing:
   02/2004 OK Implementation
   last modification:
**************************************************************************/
CMediumProperties* CMediumProperties::GetByGroupNumber(int group_number)
{
    CMediumProperties* m_mat = NULL;
    int no_mat = (int)mmp_vector.size();
    int i;
    for (i = 0; i < no_mat; i++)
    {
        m_mat = mmp_vector[i];
        if (m_mat->number == group_number)
            return m_mat;
    }
    return NULL;
}

////////////////////////////////////////////////////////////////////////////
// Properties functions
////////////////////////////////////////////////////////////////////////////

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2012 JT
**************************************************************************/
double CMediumProperties::GetEffectiveSaturationForPerm(
    const double wetting_saturation, int phase)
{
    double sl, se, slr, slm;
    sl = wetting_saturation;
    //
    switch (phase)
    {
        case 0:
            slr = residual_saturation[phase];
            slm = maximum_saturation[phase];
            sl = MRange(slr, sl, slm);
            se = (sl - slr) / (slm - slr);
            break;
        //
        case 1:
            slr = 1.0 - maximum_saturation[phase];   // slr = 1.0 - sgm
            slm = 1.0 - residual_saturation[phase];  // slm = 1.0 - sgr
            sl = MRange(slr, sl, slm);
            se = (sl - slr) / (slm - slr);
            break;
        //
        default:
            return sl;
    }
    //
    return se;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/1997 CT Erste Version
   09/1998 CT Kurven ueber Schluesselwort #CURVES
   08/1999 CT Erweiterung auf n-Phasen begonnen
   09/2002 OK case 13
   05/2003 MX case 14
   08/2004 OK Template for MMP implementation
   05/2005 MX Case 14, 15
   05/2005 OK MSH
   08/2005 WW Remove interpolation
   03/2007 WW Brooks/Corey:
   03/2012 JT All new
**************************************************************************/
double CMediumProperties::PermeabilitySaturationFunction(
    const double wetting_saturation, int phase)
{
    double kr = 0.0, sl, se, slr, slm, m, b;
    int model, gueltig;
    bool phase_shift = false;
    sl = wetting_saturation;
    //
    model = permeability_saturation_model[phase];
    if (model == 2)
    {  // krg = 1.0 - krl : get paramters for liquid phase calculation
        phase_shift = true;
        phase = 0;
        model = permeability_saturation_model[phase];
    }
    //
    switch (model)
    {
        default:
            ScreenMessage(
                "ERROR in PermeabilitySaturationFunction(): Unrecognized "
                "relative permeability method.\n");
            exit(0);
            break;
        //
        case 0:  // CURVE
            kr = GetCurveValue((int)perm_saturation_value[phase], 0, sl,
                               &gueltig);
            if (kr < minimum_relative_permeability)
                kr = minimum_relative_permeability;
            break;
        //
        case 1:  // CONSTANT VALUE
            kr = perm_saturation_value[phase];
            break;
        //
        case 2:  // krg = 1.0 - krl
            // No need to come here. Method will have been shifted to the liquid
            // phase.
            ScreenMessage(
                "ERROR in PermeabilitySaturationFunction(). Shouldn't be "
                "here.\n");
            break;
        //
        case 3:  // FUNCTION: LINEAR OR POWER --> WETTING: krw = (b*Se)^m
            slr = residual_saturation[phase];
            slm = maximum_saturation[phase];
            m = saturation_exponent[phase];
            b = perm_saturation_value[phase];
            sl = MRange(slr, sl, slm);
            se = (sl - slr) / (slm - slr);
            //
            kr = pow((b * se), m);
            if (kr < minimum_relative_permeability)
                kr = minimum_relative_permeability;
            break;
        //
        case 10:  // for unconfined 3D GW 5.3.07 JOD
            // MW correct function. did only get constants before
            /*
            b = capillary_pressure_values[0];
            slr = capillary_pressure_values[1];
            if(sl > 0 && sl < 1) {
              kr = max(0., 1 - (sl / b));
              kr = pow(kr, 2* (1 - kr));    //
              kr = slr  + (1 - slr) * kr;
            }
            else
              kr = 1;
            */
            kr = sl;
            break;
        case 33:  // FUNCTION: LINEAR OR POWER --> NON-WETTING krg =
                  // (b*(1-Se))^m
            slr = 1.0 - maximum_saturation[phase];   // slr = 1.0 - sgm
            slm = 1.0 - residual_saturation[phase];  // slm = 1.0 - sgr
            m = saturation_exponent[phase];
            b = perm_saturation_value[phase];
            sl = MRange(slr, sl, slm);
            se = (sl - slr) / (slm - slr);
            //
            kr = pow((b * (1.0 - se)), m);
            if (kr < minimum_relative_permeability)
                kr = minimum_relative_permeability;
            break;
        //
        case 4:  // 2-phase VAN GENUCHTEN/MUALEM --> WETTING
            slr = residual_saturation[phase];
            slm = maximum_saturation[phase];
            m = saturation_exponent[phase];  // always <= 1.0.  Input is
                                             // exponent = 1 / (1-lambda)
            sl = MRange(slr, sl, slm);
            se = (sl - slr) / (slm - slr);
            //
            kr = sqrt(se) * pow(1.0 - pow(1.0 - pow(se, 1.0 / m), m), 2);
            if (kr < minimum_relative_permeability)
                kr = minimum_relative_permeability;
            break;
        //
        case 44:  // 2-phase VAN GENUCHTEN/MUALEM --> NON-WETTING
            // Water Resour. Res. VOL. 23, pp2197-2206 1987. Air
            slr = 1.0 - maximum_saturation[phase];   // slr = 1.0 - sgm
            slm = 1.0 - residual_saturation[phase];  // slm = 1.0 - sgr
            m = saturation_exponent[phase];          // always <= 1.0.  Input is
                                             // exponent = 1 / (1-lambda)
            sl = MRange(slr, sl, slm);
            se = (sl - slr) / (slm - slr);
            //
            kr =
                pow(1.0 - se, 1.0 / 3.0) * pow(1.0 - pow(se, 1.0 / m), 2.0 * m);
            if (kr < minimum_relative_permeability)
                kr = minimum_relative_permeability;
            // TF: the following commented code is a new implementation of JOD
            // with this code we have different results in benchmark
            // H2/LabGasInjec/H2_Permeability_GasPressure
            //			if (sl > (maximum_saturation[phase] - MKleinsteZahl)) sl
            //= maximum_saturation[phase]
            //							- MKleinsteZahl;
            //			if (sl < (residual_saturation[phase] + MKleinsteZahl))
            // sl = residual_saturation[phase]
            //							+ MKleinsteZahl;
            //
            //			slr = residual_saturation[0];
            //			slr1 = residual_saturation[1];
            //
            //			m = saturation_exponent[phase];
            //
            //			se = (sl - slr1) / (1 - slr);
            //			//
            //			kr = pow(1.0 - se, 1.0 / 3.0) * pow(1.0 - pow(se, 1.0 /
            // m), 2.0 * m); 			if (kr < minimum_relative_permeability)
            // kr = minimum_relative_permeability;
            break;
        //
        case 6:  // 2-phase BROOKS/COREY --> WETTING
            slr = residual_saturation[phase];
            slm = maximum_saturation[phase];
            m = saturation_exponent[phase];
            sl = MRange(slr, sl, slm);
            se = (sl - slr) / (slm - slr);
            //
            kr = pow(se, 3.0 + 2.0 / m);
            if (kr < minimum_relative_permeability)
                kr = minimum_relative_permeability;
            break;
        //
        case 61:  // Brooks Corey:  CB with modified effective saturation
            slr = residual_saturation[phase];
            slm = maximum_saturation[phase];
            sl = MRange(slr, sl, slm);
            se = (sl - slr) / (1 - slr);
            m = saturation_exponent[phase];
            kr = pow(se, (2.0 + 3.0 * m) / m);
            kr = MRange(DBL_EPSILON, kr, 1.);
            break;

        case 66:  // 2-phase BROOKS/COREY --> NON-WETTING
            slr = 1.0 - maximum_saturation[phase];   // slr = 1.0 - sgm
            slm = 1.0 - residual_saturation[phase];  // slm = 1.0 - sgr
            m = saturation_exponent[phase];
            sl = MRange(slr, sl, slm);
            se = (sl - slr) / (slm - slr);
            //
            kr = pow(1.0 - se, 2) * (1.0 - pow(se, 1.0 + 2.0 / m));
            if (kr < minimum_relative_permeability)
                kr = minimum_relative_permeability;

            // TF: the following commented code is a new implementation of JOD
            // with this code we have different results in benchmark of
            // Liakopoulos type and Buckley Leverett
            //			if (sl > (maximum_saturation[phase] - MKleinsteZahl))
            //				sl = maximum_saturation[phase]- MKleinsteZahl;
            //			if (sl < (residual_saturation[phase] + MKleinsteZahl))
            //				sl = residual_saturation[phase]+ MKleinsteZahl;
            //			//
            //			se = (sl - residual_saturation[1]) / (1. -
            // residual_saturation[0] - residual_saturation[1]); 			kr =
            // pow(1.0
            // - se, 2) * (1.0 - pow(se, (2.0 + saturation_exponent[phase])
            // / saturation_exponent[phase])); 			kr =
            // MRange(minimum_relative_permeability,kr,1.);

            break;
        //
        case 7:  // 2-phase COREY'S CURVES (JT) --> WETTING
            slr = residual_saturation[phase];
            slm = maximum_saturation[phase];
            m = saturation_exponent[phase];
            sl = MRange(slr, sl, slm);
            se = (sl - slr) / (slm - slr);
            //
            kr = pow(se, 4);
            if (kr < minimum_relative_permeability)
                kr = minimum_relative_permeability;
            break;
        //
        case 77:  // 2-phase COREY'S CURVES (JT) --> NON-WETTING
            slr = 1.0 - maximum_saturation[phase];   // slr = 1.0 - sgm
            slm = 1.0 - residual_saturation[phase];  // slm = 1.0 - sgr
            m = saturation_exponent[phase];
            sl = MRange(slr, sl, slm);
            se = (sl - slr) / (slm - slr);
            //
            kr = pow((1.0 - se), 2) * (1.0 - pow(se, 2));
            if (kr < minimum_relative_permeability)
                kr = minimum_relative_permeability;
            break;
    }
    //
    if (phase_shift)
    {  // krg = 1.0 - krl : revert to gaseous phase
        kr = 1.0 - kr;
        if (kr < minimum_relative_permeability)
            kr = minimum_relative_permeability;
    }
    //
    return kr;
}

/**************************************************************************
   FEMLib-Method:
   Task: Calculate heat capacity of porous medium
   c rho = (c rho)^f + (c rho)^s
   Programing:
   04/2002 OK Implementation
   08/2004 OK MMP implementation
   03/2005 WW Case of no fluids
   10/2005 YD/OK: general concept for heat capacity
   10/2010 TF changed access to process type
**************************************************************************/
double CMediumProperties::HeatCapacity(long number, double theta,
                                       CFiniteElementStd* assem)
{
    SolidProp::CSolidProperties* m_msp = NULL;
    double heat_capacity_fluids, specific_heat_capacity_solid;
    double density_solid;
    double porosity, Sat, PG;
    int group;
    double T0, T1 = 0.0;
    //  double H0,H1;
    // ???
    bool FLOW = false;  // WW
    for (size_t ii = 0; ii < pcs_vector.size(); ii++)
        // if(pcs_vector[ii]->pcs_type_name.find("FLOW")!=string::npos) {
        if (isFlowProcess(pcs_vector[ii]->getProcessType()))
        {
            FLOW = true;
            break;
        }

    // WW     if(pcs_vector[ii]->pcs_type_name.find("RICHARDS")!=string::npos)
    // WW       nphase = 3;

    //----------------------------------------------------------------------
    switch (assem->SolidProp->GetCapacityModel())
    {
        //....................................................................
        case 0:  // f(x) user-defined function
            break;
        //....................................................................
        case 1:  // const
            // OK411
            group = m_pcs->m_msh->ele_vector[number]->GetPatchIndex();
            m_msp = msp_vector[group];
            specific_heat_capacity_solid = m_msp->Heat_Capacity();
            density_solid = fabs(m_msp->Density());
            if (FLOW)
            {
                porosity = assem->MediaProp->Porosity(number, theta);
                heat_capacity_fluids = MFPCalcFluidsHeatCapacity(assem);
            }
            else
            {
                heat_capacity_fluids = 0.0;
                porosity = 0.0;
            }
            heat_capacity =
                porosity * heat_capacity_fluids +
                (1.0 - porosity) * specific_heat_capacity_solid * density_solid;
            break;
        case 2:  // boiling model for YD
            // YD/OK: n c rho = n S^g c^g rho^g + n S^l c^l rho^l + (1-n) c^s
            // rho^s assem->GetNodalVal(1); WW
            T0 = assem->interpolate(assem->NodalVal0);
            // This following lines moved from fem_ele_std but wrong.. WW
            /*
               if(assem->FluidProp->heat_phase_change_curve>0){ //
               if(assem->FluidProp->heat_phase_change_curve>0||assem->heat_phase_change)
               { //
               if(fabs(assem->TG-T0)<1.0e-8) T1 +=1.0e-8;
               H0 = assem->interpolate(assem->NodalVal2);
               H1 = assem->interpolate(assem->NodalVal3);
               heat_capacity = (H1-H0)/(assem->TG-T0);
               }
               else //WW
               {
             */
            if (FLOW)
            {
                PG = assem->interpolate(assem->NodalValC1);
                if (assem->cpl_pcs->type == 1212)  // Multi-phase WW
                    PG *= -1.0;
                Sat = SaturationCapillaryPressureFunction(-PG);
            }
            else
                Sat = 1.0;
            T1 = assem->TG;
            if ((T1 - T0) < DBL_MIN)
                T1 *= -1;
            heat_capacity =
                assem->SolidProp->Heat_Capacity(T1, Porosity(assem), Sat);
            //  }
            break;
        case 3:  // D_THM1 - Richards model //WW
            T1 = assem->TG;
            heat_capacity = assem->SolidProp->Heat_Capacity(T1) *
                                fabs(assem->SolidProp->Density()) +
                            Porosity(assem) * MFPCalcFluidsHeatCapacity(assem);
            break;
        //....................................................................
        default:
            std::cout << "Error in CMediumProperties::HeatCapacity: no valid "
                         "material model"
                      << "\n";
            break;
            //....................................................................
    }
    return heat_capacity;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   04/2002 OK Implementation
   08/2004 MB case 5 tet, case 6 prism
   09/2004 OK MMP implementation
   03/2005 WW Case of no fluids and symmetry
   11/2005 CMCD
   03/2007 WW Conductivity for multi-phase flow
   last modification:
   ToDo:
   10/2010 TF changed access to process type
**************************************************************************/
double* CMediumProperties::HeatConductivityTensor(int number)
{
    const int group = m_pcs->m_msh->ele_vector[number]->GetPatchIndex();
    SolidProp::CSolidProperties* const m_msp = msp_vector[group];

    bool FLOW = false;  // WW
    for (size_t ii = 0; ii < pcs_vector.size(); ii++)
    {
        //		if (pcs_vector[ii]->pcs_type_name.find("FLOW") != string::npos)
        // TF
        if (isFlowProcess(pcs_vector[ii]->getProcessType()))
            FLOW = true;
    }

    double heat_conductivity_fluids = 0.0;
    double porosity = 0.0;
    if (FLOW)  // WW
    {
        CFluidProperties* m_mfp = Fem_Ele_Std->FluidProp;  // WW
        porosity = this->porosity_model_values[0];
        CRFProcess const* cpl_pcs = Fem_Ele_Std->cpl_pcs;
        if (cpl_pcs && cpl_pcs->type == 1212)  // Multi-phase WW
        {
            const double PG = Fem_Ele_Std->interpolate(
                Fem_Ele_Std->NodalValC1);  // Capillary pressure
            const double Sw =
                Fem_Ele_Std->MediaProp->SaturationCapillaryPressureFunction(PG);
            //
            m_mfp = mfp_vector[0];
            heat_conductivity_fluids = Sw * m_mfp->HeatConductivity();
            m_mfp = mfp_vector[1];
            heat_conductivity_fluids += (1.0 - Sw) * m_mfp->HeatConductivity();
        }
        else
        {
            if (Fem_Ele_Std->FluidProp->density_model == 14 &&
                Fem_Ele_Std->MediaProp->heat_diffusion_model == 1 &&
                Fem_Ele_Std->cpl_pcs)
            {
                double dens_arg[3];
                dens_arg[0] = Fem_Ele_Std->interpolate(
                    Fem_Ele_Std->NodalValC1);  // Pressure
                dens_arg[1] = Fem_Ele_Std->interpolate(
                    Fem_Ele_Std->NodalVal1);       // Temperature
                dens_arg[2] = Fem_Ele_Std->Index;  // ELE index
                heat_conductivity_fluids =
                    Fem_Ele_Std->FluidProp->HeatConductivity(dens_arg);
            }
            else
                heat_conductivity_fluids =
                    Fem_Ele_Std->FluidProp->HeatConductivity();

            if (cpl_pcs && cpl_pcs->type != 1)
            {
                double PG = Fem_Ele_Std->interpolate(
                    Fem_Ele_Std->NodalValC1);  // Capillary pressure

                if (PG < 0.0)
                {
                    const double Sw =
                        Fem_Ele_Std->MediaProp
                            ->SaturationCapillaryPressureFunction(-PG);
                    heat_conductivity_fluids *= Sw;
                    if (Fem_Ele_Std->GasProp != 0)
                        heat_conductivity_fluids +=
                            (1. - Sw) *
                            Fem_Ele_Std->GasProp->HeatConductivity();
                }
            }
        }
    }

    const int dimen = m_pcs->m_msh->GetCoordinateFlag() / 10;
    for (int i = 0; i < dimen * dimen; i++)  // MX
        heat_conductivity_tensor[i] = 0.0;

    m_msp->HeatConductivityTensor(dimen, heat_conductivity_tensor,
                                  group);  // MX

    for (int i = 0; i < dimen * dimen; i++)
        heat_conductivity_tensor[i] *= (1.0 - porosity);
    for (int i = 0; i < dimen; i++)
        heat_conductivity_tensor[i * dimen + i] +=
            porosity * heat_conductivity_fluids;

    if (evaporation == 647)
    {
        double a, b, rhow, rho_gw, rho_ga, rho_g, p_gw, mat_fac_w;
        double mat_fac_g, A, B, H_vap, dp_gw, dPc, dA, dB, dT, q;
        int GravityOn = 1;
        if ((Fem_Ele_Std->coordinate_system) % 10 != 2 &&
            (!Fem_Ele_Std->axisymmetry))
            GravityOn = 0;
        double PG2 = Fem_Ele_Std->interpolate(Fem_Ele_Std->NodalVal_p2);
        // Capillary pressure
        double PG = Fem_Ele_Std->interpolate(Fem_Ele_Std->NodalValC1);
        // Temperature
        double TG = Fem_Ele_Std->interpolate(Fem_Ele_Std->NodalVal1);
        double Sw =
            Fem_Ele_Std->MediaProp->SaturationCapillaryPressureFunction(PG);
        heat_conductivity_fluids = 0.0;
        H_vap = 2257000;  // pow((Tc - TG),0.38)*2.65E+5;
        a = 19.81;
        b = 4975.9;
        CFluidProperties* m_mfp = mfp_vector[0];
        rhow = m_mfp->Density();
        const double Rv = SpecificGasConstant::WaterVapour;
        const double expfactor = 1.0 / (rhow * Rv * TG);
        rho_gw = m_mfp->vaporDensity(TG) * exp(-PG * expfactor);
        p_gw = rho_gw * Rv * TG;
        double dens_arg[3];  // AKS
        dens_arg[0] = PG2 - p_gw;
        dens_arg[1] = TG;
        m_mfp = mfp_vector[1];
        rho_ga = m_mfp->Density(dens_arg);
        rho_g = rho_ga + rho_gw;
        m_mfp = mfp_vector[0];
        mat_fac_w = PermeabilitySaturationFunction(Sw, 0) / m_mfp->Viscosity();
        m_mfp = mfp_vector[1];
        mat_fac_g = PermeabilitySaturationFunction(Sw, 1) / m_mfp->Viscosity();
        A = b + PG / (rhow * Rv);
        B = a - log(p_gw / 30024.895431831395);
        q = heatflux;

        double Kx[3];
        for (int i = 0; i < dimen; i++)
            Kx[i] = 0.0;
        dPc = (q / (H_vap * 1.0e-13)) *
              ((1 / (rhow * mat_fac_w)) + (1 / (rho_g * mat_fac_g)));
        dA = dPc / (rhow * Rv);
        dp_gw = q / (H_vap * rho_gw * mat_fac_g * 1.0e-13);
        dB = -dp_gw / p_gw;
        dT = (B * dA - A * dB) / ((B * B) + (0 / TG));
        heat_conductivity_fluids = 2 * q / dT;
        Kx[0] = heat_conductivity_fluids;
        if (GravityOn)
        {
            dPc -= (rhow - rho_g) * gravity_constant;
            dA = dPc / (rhow * Rv);
            dp_gw -= rho_g * gravity_constant;
            dB = -dp_gw / p_gw;
            dT = (B * dA - A * dB) / ((B * B) + (0 / TG));
            heat_conductivity_fluids = 2 * q / dT;
            Kx[dimen - 1] = heat_conductivity_fluids;
        }
        for (int i = 0; i < dimen * dimen; i++)
            heat_conductivity_tensor[i] = 0.0;
        for (int i = 0; i < dimen; i++)
            heat_conductivity_tensor[i * dimen + i] = Kx[i];
    }
    return heat_conductivity_tensor;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   04/2002 OK Implementation
   09/2004 OK MMP implementation
   12/2005 CMCD
   last modification:
   ToDo:
**************************************************************************/
double* CMediumProperties::HeatDispersionTensorNew(int ip)
{
    static double heat_dispersion_tensor[9];
    double* heat_conductivity_porous_medium;
    double vg, D[9];
    double heat_capacity_fluids = 0.0;
    double fluid_density;
    double alpha_t, alpha_l;
    long index = Fem_Ele_Std->GetMeshElement()->GetIndex();
    CFluidProperties* m_mfp;
    int Dim = m_pcs->m_msh->GetCoordinateFlag() / 10;
    int i;
    ElementValue* gp_ele = ele_gp_value[index];

    // Materials
    // MX, add index
    heat_conductivity_porous_medium = HeatConductivityTensor(index);
    m_mfp = Fem_Ele_Std->FluidProp;
    fluid_density = m_mfp->Density();
    heat_capacity_fluids = m_mfp->getSpecificHeatCapacity();

    // Global Velocity
    double velocity[3] = {0., 0., 0.};
    gp_ele->getIPvalue_vec(ip, velocity);  // gp velocities
    vg = MBtrgVec(velocity, 3);

    // Dl in local coordinates
    alpha_l = heat_dispersion_longitudinal;
    alpha_t = heat_dispersion_transverse;

    if (abs(vg) > MKleinsteZahl  // For the case of diffusive transport only
                                 // WW
        && (alpha_l > MKleinsteZahl || alpha_t > MKleinsteZahl))
    {
        switch (Dim)
        {
            case 1:  // line elements
                heat_dispersion_tensor[0] =
                    heat_conductivity_porous_medium[0] +
                    alpha_l * heat_capacity_fluids * fluid_density * vg;
                break;
            case 2:
                D[0] = (alpha_t * vg) +
                       (alpha_l - alpha_t) * (velocity[0] * velocity[0]) / vg;
                D[1] = ((alpha_l - alpha_t) * (velocity[0] * velocity[1])) / vg;
                D[2] = ((alpha_l - alpha_t) * (velocity[1] * velocity[0])) / vg;
                D[3] = (alpha_t * vg) +
                       (alpha_l - alpha_t) * (velocity[1] * velocity[1]) / vg;
                for (i = 0; i < 4; i++)
                    heat_dispersion_tensor[i] =
                        heat_conductivity_porous_medium[i] +
                        (D[i] * heat_capacity_fluids * fluid_density);
                break;
            case 3:
                D[0] = (alpha_t * vg) +
                       (alpha_l - alpha_t) * (velocity[0] * velocity[0]) / vg;
                D[1] = ((alpha_l - alpha_t) * (velocity[0] * velocity[1])) / vg;
                D[2] = ((alpha_l - alpha_t) * (velocity[0] * velocity[2])) / vg;
                D[3] = ((alpha_l - alpha_t) * (velocity[1] * velocity[0])) / vg;
                D[4] = (alpha_t * vg) +
                       (alpha_l - alpha_t) * (velocity[1] * velocity[1]) / vg;
                D[5] = ((alpha_l - alpha_t) * (velocity[1] * velocity[2])) / vg;
                D[6] = ((alpha_l - alpha_t) * (velocity[2] * velocity[0])) / vg;
                D[7] = ((alpha_l - alpha_t) * (velocity[2] * velocity[1])) / vg;
                D[8] = (alpha_t * vg) +
                       (alpha_l - alpha_t) * (velocity[2] * velocity[2]) / vg;
                for (i = 0; i < 9; i++)
                    heat_dispersion_tensor[i] =
                        heat_conductivity_porous_medium[i] +
                        (D[i] * heat_capacity_fluids * fluid_density);
                break;
        }
    }
    else
        for (i = 0; i < Dim * Dim; i++)
            heat_dispersion_tensor[i] = heat_conductivity_porous_medium[i];
    return heat_dispersion_tensor;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   04/2002 OK Implementation
   09/2004 OK MMP implementation
   10/2004 SB Adapted to mass dispersion
   11/2005 CMCD
   05/2007 PCH Diffusion Coefficient corrected and Anisotropic diffusion
         coeefient computed with tortuosity added
   ToDo:
**************************************************************************/
double* CMediumProperties::MassDispersionTensorNew(int ip,
                                                   int tr_phase)  // SB + BG
{
    static double
        advection_dispersion_tensor[9];  // Name change due to static conflict
    int component = Fem_Ele_Std->pcs->pcs_component_number;
    int i;
    long index = Fem_Ele_Std->GetMeshElement()->GetIndex();
    double molecular_diffusion[9], molecular_diffusion_value;
    double vg;
    double D[9];
    double alpha_l, alpha_t;
    double theta = Fem_Ele_Std->pcs->m_num->ls_theta;
    double g[3] = {0., 0., 0.};
    double l_char = 0.0;  // OK411 volume=0.0;
    double saturation = 1.0;
    int set = 0;
    double arg, Daq, Pec;
    ElementValue* gp_ele = ele_gp_value[index];
    CompProperties* m_cp = cp_vec[component];
    // CFluidProperties* m_mfp;
    // m_mfp = Fem_Ele_Std->FluidProp;
    MshElemType::type eleType =
        m_pcs->m_msh->ele_vector[number]->GetElementType();
    int Dim = m_pcs->m_msh->GetCoordinateFlag() / 10;
    //----------------------------------------------------------------------
    // Materials
    Daq = m_cp->CalcDiffusionCoefficientCP(index, theta, m_pcs);
    molecular_diffusion_value = Daq * TortuosityFunction(index, g, theta);

    molecular_diffusion_value =
        m_cp->CalcDiffusionCoefficientCP(index, theta, m_pcs) *
        TortuosityFunction(index, g, theta);

    molecular_diffusion_value *= Porosity(index, theta);
    // CB, SB
    saturation = PCSGetEleMeanNodeSecondary_2(
        index, Fem_Ele_Std->pcs->flow_pcs_type, "SATURATION1", 1);
    if (tr_phase == 10)  // multi phase transport
        saturation = 1.0 - saturation;
    molecular_diffusion_value *= saturation;
    for (i = 0; i < Dim * Dim; i++)
        molecular_diffusion[i] = 0.0;
    for (i = 0; i < Dim; i++)
        molecular_diffusion[i * Dim + i] = molecular_diffusion_value;
    //----------------------------------------------------------------------

    // Anisotropic diffusion coefficient
    if (tortuosity_tensor_type_name[0] == 'A')
    {
        switch (Dim)
        {
            //--------------------------------------------------------------------
            case 1:  // line elements
                ;    // Do nothing
                break;
            //--------------------------------------------------------------------
            case 2:
                molecular_diffusion[0] =
                    molecular_diffusion_value * tortuosity_model_values[0];
                molecular_diffusion[1] =
                    molecular_diffusion_value * tortuosity_model_values[1];
                molecular_diffusion[2] =
                    molecular_diffusion_value * tortuosity_model_values[2];
                molecular_diffusion[3] =
                    molecular_diffusion_value * tortuosity_model_values[3];
                break;
            //--------------------------------------------------------------------
            case 3:
                molecular_diffusion[0] =
                    molecular_diffusion_value * tortuosity_model_values[0];
                molecular_diffusion[1] =
                    molecular_diffusion_value * tortuosity_model_values[1];
                molecular_diffusion[2] =
                    molecular_diffusion_value * tortuosity_model_values[2];
                molecular_diffusion[3] =
                    molecular_diffusion_value * tortuosity_model_values[3];
                molecular_diffusion[4] =
                    molecular_diffusion_value * tortuosity_model_values[4];
                molecular_diffusion[5] =
                    molecular_diffusion_value * tortuosity_model_values[5];
                molecular_diffusion[6] =
                    molecular_diffusion_value * tortuosity_model_values[6];
                molecular_diffusion[7] =
                    molecular_diffusion_value * tortuosity_model_values[7];
                molecular_diffusion[8] =
                    molecular_diffusion_value * tortuosity_model_values[8];
        }
    }

    // Global Velocity
    double velocity[3] = {0., 0., 0.};
    gp_ele->getIPvalue_vec_phase(ip, tr_phase,
                                 velocity);  // gp velocities // SB
    vg = MBtrgVec(velocity, 3);
    //  if(index < 10) cout <<" Velocity in MassDispersionTensorNew():
    //  "<<velocity[0]<<", "<<velocity[1]<<",
    //  "<<velocity[2]<<", "<<vg<<"\n";
    // test bubble velocity
    if (m_cp->bubble_velocity_model == 1)
    {
        if (index < 0)
            cout << " Velocity in MassDispersionTensorNew(): " << velocity[0]
                 << ", " << velocity[1] << ", " << velocity[2] << ", " << vg
                 << "\n";
        velocity[0] = m_cp->bubble_velocity[0];
        velocity[1] = m_cp->bubble_velocity[1];
        velocity[2] = m_cp->bubble_velocity[2];
        vg = MBtrgVec(velocity, 3);
        if (index == 100)
            std::cout << " Bubble velocity in MassDispersionTensorNew(): "
                      << velocity[0] << ", " << velocity[1] << ", "
                      << velocity[2] << ", " << vg << "\n";
    }
    // end bubble velocity
    // Dl in local coordinates
    alpha_l = mass_dispersion_longitudinal;
    alpha_t = mass_dispersion_transverse;

    // transverse Dispersion model by Chiogna et al. ES&T 2010
    // here, the dispersive part of Dmech = Dp + Ddisp is divided by vg to
    // obtain a Daq-dependent alpha_t which is used in the dispersion tensor in
    // the usual way
    if (alpha_t_model == 1)
    {
        Pec = vg * graindiameter / Daq;
        arg =
            pow(Pec, 2) / (Pec + 2 + 4 * pow(graindiameter / hydraulicrad, 2));
        alpha_t =
            Daq * pow(arg, betaexpo);  // D_disp_t (without the diffusion term)
        alpha_t /= vg;                 // aT = D_disp_t/vg
    }
    // hard stabilization
    if (this->lgpn > 0.0)
    {
        MeshLib::CElem* m_ele = NULL;
        m_ele = m_pcs->m_msh->ele_vector[index];
        if (eleType == 2)
            l_char = sqrt(m_ele->GetVolume());
        if (eleType == 4)
            l_char = sqrt(m_ele->GetVolume());
        // cout << " Element number: " << index << ", Volume: " <<
        // m_ele->GetVolume() << ", l_char: " << l_char << "\n";
        set = 0;
        if (alpha_l < l_char / lgpn)
        {
            set = 1;  // flag for output
            alpha_l = l_char / lgpn;
        }
        if (alpha_t < l_char / lgpn)
        {
            set = 1;
            alpha_t = l_char / lgpn;
        }

        // cout << " alpha_L = " << alpha_l << " < l_char/Pe; setting alpha_L =
        // " << l_char/lgpn << " for element " << index << "\n";
        if ((set > 0) & (aktueller_zeitschritt == 1) & (component < 1) &
            (ip < 1))
            std::cout << "element " << index << " " << l_char << " " << alpha_l
                      << " " << alpha_t << "\n";
    }
    //----------------------------------------------------------------------

    if (abs(vg) > MKleinsteZahl)  // For the case of diffusive transport only.
    {
        switch (Dim)
        {
            //--------------------------------------------------------------------
            case 1:  // line elements
                advection_dispersion_tensor[0] =
                    molecular_diffusion[0] + alpha_l * vg;
                break;
            //--------------------------------------------------------------------
            case 2:
                D[0] = (alpha_t * vg) +
                       (alpha_l - alpha_t) * (velocity[0] * velocity[0]) / vg;
                D[1] = ((alpha_l - alpha_t) * (velocity[0] * velocity[1])) / vg;
                D[2] = ((alpha_l - alpha_t) * (velocity[1] * velocity[0])) / vg;
                D[3] = (alpha_t * vg) +
                       (alpha_l - alpha_t) * (velocity[1] * velocity[1]) / vg;

                if (tortuosity_tensor_type_name[0] == 'A')
                    for (i = 0; i < Dim * Dim; i++)
                        advection_dispersion_tensor[i] =
                            molecular_diffusion[i] + D[i];
                else
                {
                    advection_dispersion_tensor[0] =
                        molecular_diffusion[0] + D[0];
                    // SB added - CHP, all parts of tensor required
                    advection_dispersion_tensor[1] = D[1];
                    // SB added - CHP, all parts of tensor required
                    advection_dispersion_tensor[2] = D[2];
                    advection_dispersion_tensor[3] =
                        molecular_diffusion[3] + D[3];
                }
                break;
            //--------------------------------------------------------------------
            case 3:
                D[0] = (alpha_t * vg) +
                       (alpha_l - alpha_t) * (velocity[0] * velocity[0]) / vg;
                D[1] = ((alpha_l - alpha_t) * (velocity[0] * velocity[1])) / vg;
                D[2] = ((alpha_l - alpha_t) * (velocity[0] * velocity[2])) / vg;
                D[3] = ((alpha_l - alpha_t) * (velocity[1] * velocity[0])) / vg;
                D[4] = (alpha_t * vg) +
                       (alpha_l - alpha_t) * (velocity[1] * velocity[1]) / vg;
                D[5] = ((alpha_l - alpha_t) * (velocity[1] * velocity[2])) / vg;
                D[6] = ((alpha_l - alpha_t) * (velocity[2] * velocity[0])) / vg;
                D[7] = ((alpha_l - alpha_t) * (velocity[2] * velocity[1])) / vg;
                D[8] = (alpha_t * vg) +
                       (alpha_l - alpha_t) * (velocity[2] * velocity[2]) / vg;
                if (tortuosity_tensor_type_name[0] == 'A')
                    for (i = 0; i < Dim * Dim; i++)
                        advection_dispersion_tensor[i] =
                            molecular_diffusion[i] + D[i];
                else
                {
                    /* SB changed - CHP, check - aktually the main coordinate
                       method should work advection_dispersion_tensor[0] =
                       molecular_diffusion[0] + D[0];
                       advection_dispersion_tensor[4] = molecular_diffusion[4] +
                       D[4]; advection_dispersion_tensor[8] =
                       molecular_diffusion[8] + D[8];
                     */
                    for (i = 0; i < Dim * Dim; i++)
                        advection_dispersion_tensor[i] = D[i];
                    advection_dispersion_tensor[0] += molecular_diffusion[0];
                    advection_dispersion_tensor[4] += molecular_diffusion[4];
                    advection_dispersion_tensor[8] += molecular_diffusion[8];
                }
        }
    }
    else
        for (i = 0; i < Dim * Dim; i++)
            advection_dispersion_tensor[i] = molecular_diffusion[i];
    return advection_dispersion_tensor;
}

/**************************************************************************
   FEMLib-Method:
   Task: Material tensor calculation for
   MULTI COMPONENTIAL FLOW Global Approach
   Implementaion:
   03/2011 AKS
**************************************************************************/
double* CMediumProperties::DispersionTensorMCF(int ip, int PCSIndex, int CIndex,
                                               double* variables)
{
    int k;
    double Material[9], D[9], multiplier = 1.0;
    double set, vg, fac = 0.0, alpha_l, alpha_t, g[3] = {0., 0., 0.},
                    l_char = 0.0, theta = Fem_Ele_Std->pcs->m_num->ls_theta;
    static double tensor[9];
    CFluidProperties* m_mfp;
    SolidProp::CSolidProperties* m_msp = NULL;
    int group = m_pcs->m_msh->ele_vector[number]->GetPatchIndex();
    m_msp = msp_vector[group];
    MshElemType::type eleType =
        m_pcs->m_msh->ele_vector[number]->GetElementType();
    long index = Fem_Ele_Std->GetMeshElement()->GetIndex();
    ElementValue* gp_ele = ele_gp_value[index];
    m_mfp = Fem_Ele_Std->FluidProp;
    m_mfp = mfp_vector[0];
    porosity = this->porosity_model_values[0];
    int Dim = m_pcs->m_msh->GetCoordinateFlag() / 10;
    for (k = 0; k < Dim * Dim; k++)
    {
        tensor[k] = 0.0;
        Material[k] = 0.0;
    }
    switch (PCSIndex)
    {
        case 0:  // FLOW
            fac = permeability_tensor[0] / m_mfp->Viscosity(variables);
            multiplier = 0.0;
            break;

        case 1:  // HEAT
            fac = porosity * m_mfp->HeatConductivity(variables) +
                  (1.0 - porosity) * m_msp->Heat_Conductivity(0);
            multiplier = m_mfp->Density(variables) *
                         m_mfp->SpecificHeatCapacity(variables);
            break;

        case 2:  // MASS
            fac = porosity * TortuosityFunction(index, g, theta) *
                  m_mfp->EffectiveDiffusionCoef(CIndex, variables);
            multiplier = 1.0;
            break;
    }

    for (k = 0; k < Dim; k++)
        Material[k * Dim + k] = fac;

    // Global Velocity
    double velocity[3] = {0., 0., 0.};
    gp_ele->getIPvalue_vec_phase(ip, 0, velocity);  // gp velocities // SB
    vg = MBtrgVec(velocity, 3);

    alpha_l = mass_dispersion_longitudinal;
    alpha_t = mass_dispersion_transverse;

    // hard stabilization
    if (this->lgpn > 0.0)
    {
        MeshLib::CElem* m_ele = NULL;
        m_ele = m_pcs->m_msh->ele_vector[index];
        if (eleType == 2)
            l_char = sqrt(m_ele->GetVolume());
        if (eleType == 4)
            l_char = sqrt(m_ele->GetVolume());
        // cout << " Element number: " << index << ", Volume: " <<
        // m_ele->GetVolume() << ", l_char: " << l_char << endl;
        set = 0;
        if (alpha_l < l_char / lgpn)
        {
            set = 1;  // flag for output
            alpha_l = l_char / lgpn;
        }
        if (alpha_t < l_char / lgpn)
        {
            set = 1;
            alpha_t = l_char / lgpn;
        }

        // cout << " alpha_L = " << alpha_l << " < l_char/Pe; setting alpha_L =
        // " << l_char/lgpn << " for element " << index << endl;
        if ((set > 0) & (aktueller_zeitschritt == 1) & (CIndex < 2) & (ip < 1))
            std::cout << "element " << index << " " << l_char << " " << alpha_l
                      << " " << alpha_t << std::endl;
    }
    //----------------------------------------------------------------------

    if (abs(vg) > MKleinsteZahl &&
        PCSIndex > 0)  // For the case of diffusive transport only.
    {
        switch (Dim)
        {
            //--------------------------------------------------------------------
            case 1:
                tensor[0] = Material[0] + alpha_l * vg * multiplier;
                break;
            //--------------------------------------------------------------------
            case 2:
                D[0] = (alpha_t * vg) +
                       (alpha_l - alpha_t) * (velocity[0] * velocity[0]) / vg;
                D[1] = ((alpha_l - alpha_t) * (velocity[0] * velocity[1])) / vg;
                D[2] = ((alpha_l - alpha_t) * (velocity[1] * velocity[0])) / vg;
                D[3] = (alpha_t * vg) +
                       (alpha_l - alpha_t) * (velocity[1] * velocity[1]) / vg;
                for (k = 0; k < Dim * Dim; k++)
                    tensor[k] = D[k] * multiplier;
                tensor[0] += Material[0];
                tensor[3] += Material[3];
                break;
            //--------------------------------------------------------------------
            case 3:
                D[0] = (alpha_t * vg) +
                       (alpha_l - alpha_t) * (velocity[0] * velocity[0]) / vg;
                D[1] = ((alpha_l - alpha_t) * (velocity[0] * velocity[1])) / vg;
                D[2] = ((alpha_l - alpha_t) * (velocity[0] * velocity[2])) / vg;
                D[3] = ((alpha_l - alpha_t) * (velocity[1] * velocity[0])) / vg;
                D[4] = (alpha_t * vg) +
                       (alpha_l - alpha_t) * (velocity[1] * velocity[1]) / vg;
                D[5] = ((alpha_l - alpha_t) * (velocity[1] * velocity[2])) / vg;
                D[6] = ((alpha_l - alpha_t) * (velocity[2] * velocity[0])) / vg;
                D[7] = ((alpha_l - alpha_t) * (velocity[2] * velocity[1])) / vg;
                D[8] = (alpha_t * vg) +
                       (alpha_l - alpha_t) * (velocity[2] * velocity[2]) / vg;
                for (k = 0; k < Dim * Dim; k++)
                    tensor[k] = D[k] * multiplier;
                tensor[0] += Material[0];
                tensor[4] += Material[4];
                tensor[8] += Material[8];
                break;
        }
    }
    else
    {
        for (k = 0; k < Dim * Dim; k++)
            tensor[k] = Material[k];
    }
    return tensor;
}

////////////////////////////////////////////////////////////////////////////
// DB functions
////////////////////////////////////////////////////////////////////////////

/**************************************************************************
   FEMLib-Method: CMediumProperties
   Task: get instance by name
   Programing:
   02/2004 OK Implementation
   last modification:
**************************************************************************/
CMediumProperties* CMediumProperties::GetDB(std::string mat_name)
{
    CMediumProperties* m_mat = NULL;
    std::list<CMediumProperties*>::const_iterator p_mat =
        db_mat_mp_list.begin();
    while (p_mat != db_mat_mp_list.end())
    {
        m_mat = *p_mat;
        if (mat_name.compare(m_mat->name) == 0)
            return m_mat;
        ++p_mat;
    }
    return NULL;
}

/**************************************************************************
   FEMLib-Method: CMediumProperties
   Task: set properties
   Programing:
   02/2004 OK Implementation
   last modification:
**************************************************************************/
void CMediumProperties::SetDB(std::string mat_name, std::string prop_name,
                              double value)
{
    CMediumProperties* m_mat = NULL;
    m_mat = GetDB(mat_name);
    switch (GetPropertyType(prop_name))
    {
        case 0:
            m_mat->conductivity = value;
            break;
        case 1:
            m_mat->permeability = value;
            break;
        case 2:
            m_mat->porosity = value;
            break;
    }
}

/**************************************************************************
   FEMLib-Method: CMediumProperties::GetPropertyType
   Task: get property type
   Programing:
   02/2004 OK Implementation
   last modification:
**************************************************************************/
int CMediumProperties::GetPropertyType(std::string prop_name)
{
    int counter = 0;
    string keyword_name;
    list<string>::const_iterator p_keywords = keywd_list.begin();
    while (p_keywords != keywd_list.end())
    {
        keyword_name = *p_keywords;
        if (prop_name.compare(keyword_name) == 0)
            return counter;
        ++p_keywords, counter++;
    }
    return -1;
}

////////////////////////////////////////////////////////////////////////////
// MAT-MP data base
////////////////////////////////////////////////////////////////////////////
std::string read_MAT_name(std::string in, std::string* z_rest_out)
{
    std::string mat_name;
    std::string z_rest;
    std::string delimiter(";");
    // wenn eine mg gefunden wird nach dem Schlüsselwort
    if (in.find_first_not_of(delimiter) != std::string::npos)
    {
        z_rest = in.substr(in.find_first_not_of(delimiter));
        // string für die erste (1) material group
        mat_name = z_rest.substr(0, z_rest.find_first_of(delimiter));
        *z_rest_out = z_rest.substr(mat_name.length());
        return mat_name;
    }
    else
        return "";
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   01/2004 OK/JG Implementation
   last modification:
**************************************************************************/
void MATLoadDB(std::string csv_file_name)
{
    std::string keyword("MATERIALS");
    std::string in;
    std::string line;
    std::string z_rest;
    std::string mat_name;
    std::string mat_name_tmp("MAT_NAME");
    //========================================================================
    // Read available MAT-MP keywords
    read_keywd_list();
    //========================================================================
    // File handling
    ifstream csv_file(csv_file_name.data(), ios::in);
    if (!csv_file.good())
        return;
    csv_file.seekg(0L, ios::beg);  // rewind
    //========================================================================
    // Read MATERIALS group names
    comp_keywd_list(csv_file_name);
    /*
       //iostream
       // 2.1 Create MAT-MP instances
       // search "MATERIALS"
       SOIL_PROPERTIES *sp = NULL;
       sp = MATAddSoilPropertiesObject(); //OK
       strcpy(sp->mat_type_name,"Brooklawn");
       // 2.2 Insert to db_mat_mp_list
       db_material_mp_list.push_back(sp);
       // 2.3 Read material properties
       string line;
       mat_read_hydraulic_conductivity_mean(line,sp->mat_type_name);
       // select corresponding MAT-MP
       }
     */
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   01/2004 JG Implementation
   last modification:
**************************************************************************/
void read_keywd_list(void)
{
    // File handling=============================
    std::ifstream eingabe("mat_mp_keywords.dat", ios::in);
    if (eingabe.good())
    {
        eingabe.seekg(0L, ios::beg);
        //==========================================
        std::string keyword("MATERIALS");
        std::string in;
        std::string line;
        std::string z_rest;
        std::string mat_name;
        std::string delimiter(";");
        std::string keywd_tmp("KEYWD");
        char line_char[MAX_ZEILE];
        // Read MATERIALS group names
        // 1 get string with keywords
        while (!eingabe.eof())
        {
            eingabe.getline(line_char, MAX_ZEILE);
            line = line_char;
            if (line.find_first_not_of(delimiter) != string::npos)
            {
                // schneidet delimiter ab
                in = line.substr(line.find_first_not_of(delimiter));
                keywd_tmp = in.substr(0, in.find_first_of(delimiter));
                keywd_list.push_back(keywd_tmp);
            }
            // keywd_list.remove(keyword);
        }  // eof
    }      // if eingabe.good
    else
        printf("No keyword file: mat_mp_keywords.dat");
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   01/2004 JG/OK Implementation
   last modification:
**************************************************************************/
void comp_keywd_list(std::string csv_file_name)
{
    std::string mat_names("MATERIALS");
    std::string in;
    std::string line;
    std::string z_rest;
    std::string mat_name;
    std::string mat_name_tmp("MAT_NAME");
    char line_char[MAX_ZEILE];
    std::string in1;  // zwischenstring zum abschneiden der einheit
    std::string delimiter(";");
    std::string keyword("MATERIALS");
    double kwvalue;

    // File handling------------------------------------
    std::string sp;
    std::ifstream eingabe(csv_file_name.data(), ios::in);
    if (eingabe.good())
    {
        eingabe.seekg(0L, ios::beg);  // rewind um materialgruppen auszulesen
        eingabe.getline(line_char, MAX_ZEILE);
        line = line_char;
        //-----------------------------------------------
        // 1 - read MAT names
        // ReadMATNames(csv_file_name);
        if (line.find(mat_names) != string::npos)
        {
            in = line.substr(mat_names.length() + 1);
            // 2 Read material group names
            while (!mat_name_tmp.empty())
            {
                mat_name_tmp = read_MAT_name(in, &z_rest);
                if (mat_name_tmp.empty())
                    break;
                else
                {
                    mat_name = mat_name_tmp;
                    mat_name_list.push_back(mat_name);
                    in = z_rest;
                }
            }  // while mat_name
        }      // keyword found
        //-----------------------------------------------
        // 2 - create MAT-SP instances
        CMediumProperties* m_mat_mp = NULL;
        list<string>::const_iterator p_mat_names = mat_name_list.begin();
        while (p_mat_names != mat_name_list.end())
        {
            mat_name = *p_mat_names;
            m_mat_mp = new CMediumProperties;
            m_mat_mp->name = mat_name;
            db_mat_mp_list.push_back(m_mat_mp);
            ++p_mat_names;
        }
        //-----------------------------------------------
        // 3 - read MAT properties

        // 1 get string where keyword hits
        CMediumProperties* m_mat_mp1 = NULL;
        while (!eingabe.eof())
        {
            eingabe.getline(line_char, MAX_ZEILE);
            line = line_char;
            string sp;
            list<string>::const_iterator pm = keywd_list.begin();
            while (pm != keywd_list.end())
            {
                sp = *pm;
                if (line.find(sp) != string::npos)
                {
                    // schneidet keyword ab
                    in1 = line.substr(line.find_first_not_of(delimiter) +
                                      sp.length() + 1);
                    // schneidet keyword ab
                    in = in1.substr(in1.find_first_of(delimiter));
                    // Schleife über alle MAT-Gruppen
                    list<string>::const_iterator p_mat_names =
                        mat_name_list.begin();
                    while (p_mat_names != mat_name_list.end())
                    {
                        mat_name = *p_mat_names;
                        m_mat_mp1 = m_mat_mp1->GetDB(mat_name);
                        mat_name_tmp = read_MAT_name(in, &z_rest);
                        kwvalue = atof(mat_name_tmp.data());
                        // Val = strtod(pDispInfo->item.strText, NULL);
                        m_mat_mp1->SetDB(mat_name, sp, kwvalue);
                        ++p_mat_names;
                    }
                }  //
                ++pm;
            }  // kwlist
        }      // eof
    }          // if eingabe.good
}

///*************************************************************************************************
////////////////////////////////////////////////////////////////////////////
// Properties functions
////////////////////////////////////////////////////////////////////////////

/*************************************************************************
   ROCKFLOW - Funktion: SetSoilPropertiesDefaultsClay

   Aufgabe:

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E: SOIL_PROPERTIES *sp: Zeiger auf SOIL_PROPERTIES

   Ergebnis:
   - void -

   Programmaenderungen:
   10/2003   CMCD   Erste Version

*************************************************************************/

void CMediumProperties::SetMediumPropertiesDefaultsClay(void)
{
    int i;
    geo_dimension = 3;
    geo_area = 1.0;
    porosity = 0.4;
    tortuosity = 1.0;
    storage = 1.5e-10;  // m/pa
    for (i = 0; i < 9; i++)
        permeability_tensor[i] = 1.e-17;
    heat_dispersion_longitudinal = 0.01;
    heat_dispersion_transverse = 0.01;
    mass_dispersion_longitudinal = 0.01;
    mass_dispersion_transverse = 0.01;
    for (i = 0; i < 9; i++)
        heat_conductivity_tensor[i] = 3.2;
}

/*************************************************************************
   ROCKFLOW - Funktion: SetSoilPropertiesDefaultsSilt

   Aufgabe:

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E: SOIL_PROPERTIES *sp: Zeiger auf SOIL_PROPERTIES

   Ergebnis:
   - void -

   Programmaenderungen:
   10/2003   CMCD   Erste Version

*************************************************************************/

void CMediumProperties::SetMediumPropertiesDefaultsSilt(void)
{
    int i;
    geo_dimension = 3;
    geo_area = 1.0;
    porosity = 0.4;
    tortuosity = 1.0;
    storage = 5.e-7;  // m/pa
    // m²
    for (i = 0; i < 9; i++)
        permeability_tensor[i] = 1.e-14;
    heat_dispersion_longitudinal = 0.01;
    heat_dispersion_transverse = 0.01;
    mass_dispersion_longitudinal = 0.01;
    mass_dispersion_transverse = 0.01;
    for (i = 0; i < 9; i++)
        heat_conductivity_tensor[i] = 3.2;
}

/*************************************************************************
   ROCKFLOW - Funktion: SetSoilPropertiesDefaultsSand

   Aufgabe:

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E: SOIL_PROPERTIES *sp: Zeiger auf SOIL_PROPERTIES

   Ergebnis:
   - void -

   Programmaenderungen:
   10/2003   CMCD   Erste Version

*************************************************************************/
void CMediumProperties::SetMediumPropertiesDefaultsSand(void)
{
    int i;
    geo_dimension = 3;
    geo_area = 1.0;
    porosity = 0.3;
    tortuosity = 1.0;
    storage = 5.e-8;  // m/pa
    // m²
    for (i = 0; i < 9; i++)
        permeability_tensor[i] = 1.e-11;
    heat_dispersion_longitudinal = 0.01;
    heat_dispersion_transverse = 0.01;
    mass_dispersion_longitudinal = 0.01;
    mass_dispersion_transverse = 0.01;
    for (i = 0; i < 9; i++)
        heat_conductivity_tensor[i] = 3.2;
}

/*************************************************************************
   ROCKFLOW - Funktion: SetSoilPropertiesDefaultsGravel

   Aufgabe:

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E: SOIL_PROPERTIES *sp: Zeiger auf SOIL_PROPERTIES

   Ergebnis:
   - void -

   Programmaenderungen:
   10/2003   CMCD   Erste Version

*************************************************************************/

void CMediumProperties::SetMediumPropertiesDefaultsGravel(void)
{
    int i;
    geo_dimension = 3;
    geo_area = 1.0;
    porosity = 0.32;
    tortuosity = 1.0;
    storage = 5.e-9;  // m/pa
    // m²
    for (i = 0; i < 9; i++)
        permeability_tensor[i] = 1.e-9;
    heat_dispersion_longitudinal = 0.01;
    heat_dispersion_transverse = 0.01;
    mass_dispersion_longitudinal = 0.01;
    mass_dispersion_transverse = 0.01;
    for (i = 0; i < 9; i++)
        heat_conductivity_tensor[i] = 3.2;
}

/*************************************************************************
   ROCKFLOW - Funktion: SetMediumPropertiesDefaultsCrystaline

   Aufgabe:

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E: SOIL_PROPERTIES *sp: Zeiger auf SOIL_PROPERTIES

   Ergebnis:
   - void -

   Programmaenderungen:
   10/2003   CMCD   Erste Version

*************************************************************************/
void CMediumProperties::SetMediumPropertiesDefaultsCrystalline(void)
{
    int i;
    geo_dimension = 3;
    geo_area = 1.0;
    porosity = 0.05;
    tortuosity = 1.0;
    storage = 5.e-11;  // m/pa
    // m²
    for (i = 0; i < 9; i++)
        permeability_tensor[i] = 1.e-14;
    heat_dispersion_longitudinal = 0.01;
    heat_dispersion_transverse = 0.01;
    mass_dispersion_longitudinal = 0.01;
    mass_dispersion_transverse = 0.01;
    for (i = 0; i < 9; i++)
        heat_conductivity_tensor[i] = 3.2;
}

/*************************************************************************
   ROCKFLOW - Funktion: SetMediumPropertiesDefaultsCrystaline

   Aufgabe:

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E: SOIL_PROPERTIES *sp: Zeiger auf SOIL_PROPERTIES

   Ergebnis:
   - void -

   Programmaenderungen:
   10/2003   CMCD   Erste Version

*************************************************************************/
void CMediumProperties::SetMediumPropertiesDefaultsBordenAquifer(void)
{
    int i;
    geo_dimension = 3;
    geo_area = 1.0;
    porosity = 0.10;
    tortuosity = 1.0;
    storage = 5.e-11;  // m/pa
    // m²
    for (i = 0; i < 9; i++)
        permeability_tensor[i] = 1.e-12;
    heat_dispersion_longitudinal = 0.01;
    heat_dispersion_transverse = 0.01;
    mass_dispersion_longitudinal = 0.01;
    mass_dispersion_transverse = 0.01;
    for (i = 0; i < 9; i++)
        heat_conductivity_tensor[i] = 3.2;
}

/*------------------------------------------------------------------------*/
/* MAT-MP geometric properties */
/* 3 porosity */
/*------------------------------------------------------------------------*/

/**************************************************************************
   FEMLib-Method:
   Task: Porosity calc function
   Case overview
   0 Curve function
   1 Constant Value
   2 Function of normal stress from Geomechnical model
   3 Free swelling
   4 Constraint swelling
   Programing:
   07/2004 OK C++ Implementation
   08/2004	CMCD Re-written based on MATCalcPorosity
   10/2004 MX Swelling processes
   04/2007 WW Porosity by gauss stress value
   last modification:
 *************************************************************************/
double CMediumProperties::Porosity(long number, double theta)
{
    int nidx0, nidx1;
    double primary_variable[PCS_NUMBER_MAX];
    int gueltig;
#ifdef GEM_REACT
    int idx;
#endif
    double porosity_sw;
    CFiniteElementStd* assem = m_pcs->GetAssember();
    string str;
    ///
    ElementValue_DM* gval = NULL;

    // CB Get idx of porosity in elements mat vector for het porosity
    size_t por_index(0);

    if (porosity_model == 11)
        for (por_index = 0; por_index < m_pcs->m_msh->mat_names_vector.size();
             por_index++)
            if (m_pcs->m_msh->mat_names_vector[por_index].compare("POROSITY") ==
                0)
                break;

    // Functional dependencies
    CRFProcess* pcs_temp;

    const size_t no_pcs_names(this->porosity_pcs_name_vector.size());
    for (size_t i = 0; i < no_pcs_names; i++)
    {
        str = porosity_pcs_name_vector[i];
        pcs_temp = PCSGet(str, true);  // MX
        nidx0 = pcs_temp->GetNodeValueIndex(porosity_pcs_name_vector[i]);
        nidx1 = nidx0 + 1;
        if (mode == 0)  // Gauss point values
        {
            // assem->ComputeShapefct(1);
            primary_variable[i] =
                (1. - theta) * assem->interpolate(nidx0, pcs_temp) +
                theta * assem->interpolate(nidx1, pcs_temp);
        }  // Node values
        else if (mode == 1)
            primary_variable[i] =
                (1. - theta) * pcs_temp->GetNodeValue(number, nidx0) +
                theta * pcs_temp->GetNodeValue(number, nidx1);
        // Element average value
        else if (mode == 2)
            primary_variable[i] =
                (1. - theta) * assem->elemnt_average(nidx0, pcs_temp) +
                theta * assem->elemnt_average(nidx1, pcs_temp);
        else if (mode == 1)  // Node values

            primary_variable[i] =
                (1. - theta) * pcs_temp->GetNodeValue(number, nidx0) +
                theta * pcs_temp->GetNodeValue(number, nidx1);
        else if (mode == 2)  // Element average value

            primary_variable[i] =
                (1. - theta) * assem->elemnt_average(nidx0, pcs_temp) +
                theta * assem->elemnt_average(nidx1, pcs_temp);
    }

    //----------------------------------------------------------------------
    // Material models
    /*
       if(permeability_stress_mode==3) //Barton-Bandis WW
       {
       int i;
       double w[3], TG = 0.0;
       if(assem->cpl_pcs)
          TG = assem->interpolate(assem->NodalValC1);
       else
          TG = 296.0;
       CalStressPermeabilityFactor(w, TG);
       porosity = 0.0;
       for(i=0; i<geo_dimension; i++)
       {
       if(i==0)
       porosity =
       pow(18.0e+6*permeability_tensor[i]*w[i]/(c_coefficient[21+i]*c_coefficient[21+i]),
       1.0/(float)geo_dimension);
       else
       porosity *=
       pow(18.0e+6*permeability_tensor[i]*w[i]/(c_coefficient[21+i]*c_coefficient[21+i]),
       1.0/(float)geo_dimension);
       }
       //
       return porosity;
       }
     */
    switch (porosity_model)
    {
        case 0:  // n = f(x)
            porosity =
                GetCurveValue(fct_number, 0, primary_variable[0], &gueltig);
            break;
        case 1:  // n = const
            porosity = porosity_model_values[0];
            break;
        case 2:  // n = f(sigma_eff), Stress dependance
            porosity = PorosityEffectiveStress(number, primary_variable[0]);
            break;
        case 3:  // n = f(S), Free chemical swelling
            porosity = PorosityVolumetricFreeSwellingConstantIonicstrength(
                number, primary_variable[0], primary_variable[1]);
            break;
        case 4:  // n = f(S), Constrained chemical swelling
            porosity =
                PorosityEffectiveConstrainedSwellingConstantIonicStrength(
                    number, primary_variable[0], primary_variable[1],
                    &porosity_sw);
            break;
        case 5:  // n = f(S), Free chemical swelling, I const
            porosity = PorosityVolumetricFreeSwelling(
                number, primary_variable[0], primary_variable[1]);
            break;
        case 6:  // n = f(S), Constrained chemical swelling, I const
            porosity = PorosityEffectiveConstrainedSwelling(
                number, primary_variable[0], primary_variable[1], &porosity_sw);
            break;
        case 7:  // n = f(mean stress) WW
            gval = ele_value_dm[number];
            primary_variable[0] = -gval->MeanStress(assem->gp) / 3.0;
            porosity =
                GetCurveValue(porosity_curve, 0, primary_variable[0], &gueltig);
            break;
        case 10:
            /* porosity change through dissolution/precipitation */
            porosity = PorosityVolumetricChemicalReaction(number);
            break;
        case 11:  // n = temp const, but spatially distributed CB
            // porosity = porosity_model_values[0];
            porosity = _mesh->ele_vector[number]->mat_vector[por_index];
            break;
        case 12:  // n = n0 + vol_strain, WX: 03.2011
            porosity =
                PorosityVolStrain(number, porosity_model_values[0], assem);
            break;
        case 13:
        {
            CRFProcess* m_pcs_flow = PCSGetFlow();
            const int idx_n = m_pcs_flow->GetElementValueIndex("POROSITY");

            // porosity change through dissolution/precipitation
            // Here, you should access porosity from the element value vector of
            // the flow process so you have to get the index of porosity above,
            // if porosity model = 13
            porosity = m_pcs_flow->GetElementValue(number, idx_n + 1);
            break;
        }
#ifdef GEM_REACT
        case 15:
            porosity = porosity_model_values[0];  // default value as backup

            for (size_t i = 0; i < pcs_vector.size(); i++)
                //		if ((pcs_vector[i]->pcs_type_name.find("FLOW") !=
                // string::npos)) {
                if (isFlowProcess(pcs_vector[i]->getProcessType()))
                {
                    idx = pcs_vector[i]->GetElementValueIndex("POROSITY");
                    porosity = pcs_vector[i]->GetElementValue(
                        number, idx + 1);  // always return new/actual value
                    if (porosity < 0.0 || porosity > 1.0)
                    {
                        cout << "Porosity: error getting porosity for model "
                                "15. porosity: "
                             << porosity << " at node " << number << "\n";
                        porosity = porosity_model_values[0];
                    }
                }

            break;
#endif
#ifdef BRNS
        case 16:
            porosity = porosity_model_values[0];  // default value as backup
            if (aktueller_zeitschritt > 1)
                for (size_t i = 0; i < pcs_vector.size(); i++)
                {
                    pcs_temp = pcs_vector[i];
                    //	            if (
                    // pcs_temp->pcs_type_name.compare("GROUNDWATER_FLOW") == 0
                    //||
                    //                     pcs_temp->pcs_type_name.compare("LIQUID_FLOW")
                    //                     == 0         ) {
                    if (pcs_temp->getProcessType() ==
                            FiniteElement::GROUNDWATER_FLOW ||
                        pcs_temp->getProcessType() ==
                            FiniteElement::LIQUID_FLOW)
                    {
                        int idx;
                        idx = pcs_temp->GetElementValueIndex("POROSITY");

                        porosity = pcs_temp->GetElementValue(number, idx);
                        if (porosity < 1.e-6)
                            cout << "error for porosity1 " << porosity
                                 << " node " << number << "\n";
                    }
                }
            break;
#endif
        default:
            cout << "Unknown porosity model!"
                 << "\n";
            break;
    }
    return porosity;
}

/*------------------------------------------------------------------------*/
/* MAT-MP geometric properties */
/* 3 porosity */
/*------------------------------------------------------------------------*/

/**************************************************************************
   FEMLib-Method:
   Task: Porosity calc function
   Case overview
   0 Curve function
   1 Constant Value
   2 Function of normal stress from Geomechnical model
   3 Free swelling
   4 Constraint swelling
   Programing:
   07/2004 OK C++ Implementation
   08/2004	CMCD Re-written based on MATCalcPorosity
   10/2004 MX Swelling processes
   04/2007 WW Porosity by gauss stress value
   last modification:
 *************************************************************************/
// WW
double CMediumProperties::Porosity(CElement* assem)
{
    static int nidx0, nidx1, idx_n;
    double primary_variable[PCS_NUMBER_MAX];
    int gueltig;
    double porosity_sw, theta;
    std::string str;
    ///
    ElementValue_DM* gval = NULL;
    CFiniteElementStd* assem_tmp =
        m_pcs->GetAssember();  // WX: for poro vol strain. 03.2011

    //----------------------------------------------------------------------
    // Functional dependencies
    number = assem->GetElementIndex();
    CRFProcess* pcs_temp;
    CRFProcess* m_pcs_flow;

    const size_t no_pcs_names(porosity_pcs_name_vector.size());
    for (size_t i = 0; i < no_pcs_names; i++)
    {
        str = porosity_pcs_name_vector[i];
        pcs_temp = PCSGet(str, true);       // MX
        theta = pcs_temp->m_num->ls_theta;  // WW
        nidx0 = pcs_temp->GetNodeValueIndex(porosity_pcs_name_vector[i]);
        nidx1 = nidx0 + 1;
        if (mode == 0)  // Gauss point values
        {
            // assem->ComputeShapefct(1);
            primary_variable[i] =
                (1. - theta) * assem->interpolate(nidx0, pcs_temp) +
                theta * assem->interpolate(nidx1, pcs_temp);
        }
        else if (mode == 1)  // Node values

            primary_variable[i] =
                (1. - theta) * pcs_temp->GetNodeValue(number, nidx0) +
                theta * pcs_temp->GetNodeValue(number, nidx1);
        else if (mode == 2)  // Element average value

            primary_variable[i] =
                (1. - theta) * assem->elemnt_average(nidx0, pcs_temp) +
                theta * assem->elemnt_average(nidx1, pcs_temp);
    }
    //----------------------------------------------------------------------
    // Material models
    switch (porosity_model)
    {
        case 0:  // n = f(x)
            porosity =
                GetCurveValue(fct_number, 0, primary_variable[0], &gueltig);
            break;
        case 1:  // n = const
            porosity = porosity_model_values[0];
            break;
        case 2:  // n = f(sigma_eff), Stress dependance
            porosity = PorosityEffectiveStress(number, primary_variable[0]);
            break;
        case 3:  // n = f(S), Free chemical swelling
            porosity = PorosityVolumetricFreeSwellingConstantIonicstrength(
                number, primary_variable[0], primary_variable[1]);
            break;
        case 4:  // n = f(S), Constrained chemical swelling
            porosity =
                PorosityEffectiveConstrainedSwellingConstantIonicStrength(
                    number, primary_variable[0], primary_variable[1],
                    &porosity_sw);
            break;
        case 5:  // n = f(S), Free chemical swelling, I const
            porosity = PorosityVolumetricFreeSwelling(
                number, primary_variable[0], primary_variable[1]);
            break;
        case 6:  // n = f(S), Constrained chemical swelling, I const
            porosity = PorosityEffectiveConstrainedSwelling(
                number, primary_variable[0], primary_variable[1], &porosity_sw);
            break;
        case 7:  // n = f(mean stress) WW
            gval = ele_value_dm[number];
            primary_variable[0] = -gval->MeanStress(assem->GetGPindex()) / 3.0;
            porosity =
                GetCurveValue(porosity_curve, 0, primary_variable[0], &gueltig);
            break;
        case 10:
            /* porosity change through dissolution/precipitation */
            porosity = PorosityVolumetricChemicalReaction(number);
            break;
        case 11:  // n = const, but spatially distributed CB
            porosity = porosity_model_values[0];
            break;
        case 12:  // n = n0 + vol_strain
            porosity = PorosityVolStrain(number, porosity_model_values[0],
                                         assem_tmp);  // WX:03.2011
            break;
        case 13:
            m_pcs_flow = PCSGetFlow();
            idx_n = m_pcs_flow->GetElementValueIndex("POROSITY");
            porosity = m_pcs_flow->GetElementValue(number, idx_n + 1);
            break;
#ifdef GEM_REACT
        case 15:

            porosity = porosity_model_values[0];  // default value as backup

            //                CRFProcess* mf_pcs = NULL;
            for (size_t i = 0; i < pcs_vector.size(); i++)
            {
                pcs_temp = pcs_vector[i];
                //			if
                //((pcs_temp->pcs_type_name.compare("GROUNDWATER_FLOW")
                //== 0) || (pcs_temp->pcs_type_name.compare("RICHARDS_FLOW") ==
                // 0)||(pcs_temp->pcs_type_name.compare("MULTI_PHASE_FLOW") ==
                // 0))
                if ((pcs_temp->getProcessType() ==
                     FiniteElement::GROUNDWATER_FLOW) ||
                    (pcs_temp->getProcessType() == FiniteElement::RICHARDS_FLOW)
                    // TF
                    || (pcs_temp->getProcessType() ==
                        FiniteElement::MULTI_PHASE_FLOW))
                {
                    int idx = pcs_temp->GetElementValueIndex("POROSITY");
                    porosity = pcs_temp->GetElementValue(number, idx);
                    if (porosity < 1.e-6 || porosity > 1.0)
                    {
                        std::cout << " error getting porosity model 15 "
                                  << porosity << " node " << number << "\n";
                        porosity = porosity_model_values[0];
                    }
                }
            }
            // KG44: TODO!!!!!!!!!!!!! check the above  ***************
            break;
#endif

        default:
            DisplayMsgLn("Unknown porosity model!");
            break;
    }
    //----------------------------------------------------------------------
    return porosity;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   05/2005 MX Implementation
   last modification:
**************************************************************************/
double
CMediumProperties::PorosityEffectiveConstrainedSwellingConstantIonicStrength(
    long index, double saturation, double temperature, double* porosity_sw)
{
    porosity_sw = porosity_sw;  // WW Remove this argument
    /*  Soil Properties */

    double mat_mp_m, beta;
    double porosity = 0.0;
    // WW static double theta;
    double porosity_n, porosity_IL, d_porosity, density_rock, fmon,
        porosity_max;
    double n_sw, det_n_sw;
    double old_value, new_value;
    static double porosity_IP0;
    static double S_0;
    /* Component Properties */
    static double satu, epsilon;
    static double satu_0 = 0.20;
    static double porosity_min = 0.05;
    static double ion_strength;
    static double F_const = 96484.6, epsilon_0 = 8.854e-12;
    static double R = 8.314510, psi = 1.0;
    /* Fluid properies */
    // WW int phase=1;

    //--------------------------------------------------------------------
    // MMP medium properties
    S_0 = porosity_model_values[1];       // Specific surface area [m^2/g]
    fmon = porosity_model_values[2];      // Anteil quelfaehige mineral [-]
    mat_mp_m = porosity_model_values[3];  // Schichtenanzahl eines quellfähigen
                                          // Partikels [-]
    beta = porosity_model_values[6];  // modifications coefficient (z.B. free
                                      // swelling Weimar beta=3.0)
    //--------------------------------------------------------------------
    // MSP solid properties
    SolidProp::CSolidProperties* m_msp = NULL;
    // long group = ElGetElementGroupNumber(index);
    long group = m_pcs->m_msh->ele_vector[index]->GetPatchIndex();
    m_msp = msp_vector[group];
    density_rock = fabs(m_msp->Density(1));
    //--------------------------------------------------------------------

    /* Component Properties */
    ion_strength = porosity_model_values[4]; /*Ionic strength [M]*/
    satu_0 = porosity_model_values[5]; /*Initial saturation, if the sample is
                                          homogenous [-]*/
    porosity_min =
        porosity_model_values[7]; /*minimal porosity after swelling compaction*/

    /* Field v0ariables */
    // WW theta = m_pcs->m_num->ls_theta;
    // WW phase=1;
    satu = saturation; /*only for fluid, phase=1*/

    /*-----------------------------------------------------------------------*/
    /* Interlayer Porositaet berechnen */
    if (temperature < 0.0 || temperature > 600.0)
        temperature = 298.0;  // ToDo, MX
    epsilon = 87.0 + exp(-0.00456 *
                         (temperature - PhysicalConstant::CelsiusZeroInKelvin));
    porosity_n = porosity_model_values[0];

    /* Maximal inter layer porosity */
    porosity_max = fmon * psi * mat_mp_m * S_0 * (density_rock * 1.0e3) *
                   sqrt(epsilon * epsilon_0 * R * temperature /
                        (2000.0 * F_const * F_const * ion_strength));
    d_porosity = porosity_max * (pow(satu, beta) - pow(satu_0, beta));

    /*-----------Interlayer porosity calculation------------------*/

    /*  porosity_IL = porosity_IL*satu; */
    porosity_IL = porosity_max * pow(satu, beta);
    m_pcs->SetElementValue(
        index, m_pcs->GetElementValueIndex("POROSITY_IL") + 1, porosity_IL);

    /*-----------Effective porosity calculation------------------*/

    porosity = porosity_n - d_porosity;

    if (porosity < porosity_min)
        porosity = porosity_min;

    m_pcs->SetElementValue(index, m_pcs->GetElementValueIndex("POROSITY") + 1,
                           porosity);

    /*-----------Swelling strain rate, xie, wang and Kolditz,
     * (23)------------------*/

    porosity_IP0 = porosity_n - porosity_max * pow(satu_0, beta);
    n_sw = d_porosity - (porosity_IP0 - porosity_min);

    if (n_sw <= 0.0)
        n_sw = 0.0;
    int indx_sw0 = m_pcs->GetElementValueIndex("n_sw");
    int indx_sw1 = indx_sw0 + 1;

    if (aktuelle_zeit == dt)
    {
        m_pcs->SetElementValue(index, indx_sw0, n_sw);
        m_pcs->SetElementValue(index, indx_sw1, n_sw);
    }
    if (aktuelle_zeit >= 2 * dt)
        m_pcs->SetElementValue(index, indx_sw1, n_sw);

    old_value = m_pcs->GetElementValue(index, indx_sw0);
    new_value = m_pcs->GetElementValue(index, indx_sw1);

    if (index == 399 && aktuelle_zeit == 10 * dt)

        index = index;  // only for debug

    // change rate of swelling strain
    det_n_sw = new_value - old_value;
    if (det_n_sw < -1.0e-8)
        det_n_sw = det_n_sw;
    if (aktuelle_zeit == dt)
        m_pcs->SetElementValue(index, m_pcs->GetElementValueIndex("n_sw_rate"),
                               det_n_sw);
    m_pcs->SetElementValue(index, m_pcs->GetElementValueIndex("n_sw_rate") + 1,
                           det_n_sw);

    return porosity;
}

//------------------------------------------------------------------------
// 4..TORTUOSITY
//------------------------------------------------------------------------

/* 4B.4.3 non-linear flow */

/**************************************************************************
   9 Storage

**************************************************************************
   11 Permeability
**************************************************************************
   FEMLib-Method:
   Task: Master calc function
   Programing:
   08/2004 OK MMP implementation
           based on GetPermeabilityTensor by OK
   10/2004 SB adapted to het_file
   last modification:
   10/2010 TF changed access to process type
**************************************************************************/
double* CMediumProperties::PermeabilityTensor(long index)
{
    static double tensor[9];
    unsigned int perm_index = 0;

    int idx_k, idx_n;
    double /*k_old, n_old,*/ k_new, n_new, k_rel, n_rel;
    CRFProcess* m_pcs_tmp(NULL);
    if ((permeability_model == 8) && (porosity_model == 13))
        m_pcs_tmp = PCSGetFlow();

    // HS: move the following loop into the "if ( permeability_tensor_type == 0
    // )" scope.---- this is not necessary for in-isotropic case;
    // if(permeability_model==2)
    //    for(perm_index=0;perm_index<(int)m_pcs->m_msh->mat_names_vector.size();perm_index++)
    //        if(m_pcs->m_msh->mat_names_vector[perm_index].compare("PERMEABILITY")==0)
    //              break;
    // end of comment out
    // codes--------------------------------------------------------------

    // -------------------------------------------------------------------------------------------------------
    // Start of K-C relationship. This will write a value into tensor[0] first,
    // the values depends on the value in permability_model. for 3 and 4, it
    // gets the values from K-C relationship. this will only influence then case
    // when permeability_tensor_type = 0
    // -------------------------------------------------------------------------------------------------------

    if (permeability_tensor_type == 0)
    {  // only when permeability is isotropic
        tensor[0] = permeability_tensor[0];

        if (permeability_model == 3)
        {  // HS: 11.2008, for K-C relationship
            k_new = tensor[0];
            CRFProcess* pcs_tmp(NULL);
            for (size_t i = 0; i < pcs_vector.size(); i++)
            {
                pcs_tmp = pcs_vector[i];
                if (pcs_tmp->getProcessType() ==
                        FiniteElement::GROUNDWATER_FLOW ||
                    pcs_tmp->getProcessType() == FiniteElement::LIQUID_FLOW ||
                    pcs_tmp->getProcessType() == FiniteElement::RICHARDS_FLOW)
                    break;
            }
            // get indexes
            idx_k = pcs_tmp->GetElementValueIndex("PERMEABILITY");
            idx_n = pcs_tmp->GetElementValueIndex("POROSITY");

            // get values of k0, n0, and n.
            k_new = pcs_tmp->GetElementValue(index, idx_k + 1);
            n_new = pcs_tmp->GetElementValue(index, idx_n + 1);

            // if first time step, get the k_new from material class
            if (aktueller_zeitschritt == 1)  // for the first time step.....
            {
                // get the permeability.
                //				KC_permeability_initial = k_new = tensor[0];
                //				KC_porosity_initial = n_new;
            }
            // save old permeability
            pcs_tmp->SetElementValue(index, idx_k, k_new);

            // calculate new permeability
            k_new = CMediumProperties::KozenyCarman(KC_permeability_initial,
                                                    KC_porosity_initial, n_new);

            // save new permeability
            pcs_tmp->SetElementValue(index, idx_k + 1, k_new);

            // now gives the newly calculated value to tensor[]
            tensor[0] = k_new;
        }
        else if (permeability_model == 4)
        {  // HS: 11.2008, for K-C_normalized relationship
            k_new = tensor[0];
            CRFProcess* pcs_tmp(NULL);
            for (size_t i = 0; i < pcs_vector.size(); i++)
            {
                pcs_tmp = pcs_vector[i];
                if (pcs_tmp->getProcessType() ==
                        FiniteElement::GROUNDWATER_FLOW ||
                    pcs_tmp->getProcessType() == FiniteElement::LIQUID_FLOW ||
                    pcs_tmp->getProcessType() == FiniteElement::RICHARDS_FLOW)
                    break;
            }
            // get indexes
            idx_k = pcs_tmp->GetElementValueIndex("PERMEABILITY");
            idx_n = pcs_tmp->GetElementValueIndex("POROSITY");

            // get values of k0, n0, and n.
            k_new = pcs_tmp->GetElementValue(index, idx_k + 1);
            n_new = pcs_tmp->GetElementValue(index, idx_n + 1);

            // if first time step, get the k_new from material class
            if (aktueller_zeitschritt == 1)  // for the first time step.....
            {
                // get the permeability.
                //				KC_permeability_initial = k_new = tensor[0];
                //				KC_porosity_initial = n_new;
            }
            // save old permeability
            pcs_tmp->SetElementValue(index, idx_k, k_new);

            // calculate new permeability
            k_new = CMediumProperties::KozenyCarman_normalized(
                KC_permeability_initial, KC_porosity_initial, n_new);

            // save new permeability
            pcs_tmp->SetElementValue(index, idx_k + 1, k_new);

            // now gives the newly calculated value to tensor[]
            tensor[0] = k_new;
        }
        else if (permeability_model == 5)
        {  // HS: 11.2008, for Clement clogging model
            // if first time step, do nothing. otherwise,
            if (aktueller_zeitschritt > 1)
            {
                CRFProcess* pcs_tmp(NULL);
                for (size_t i = 0; i < pcs_vector.size(); i++)
                {
                    pcs_tmp = pcs_vector[i];
                    if (pcs_tmp->getProcessType() ==
                            FiniteElement::GROUNDWATER_FLOW ||
                        pcs_tmp->getProcessType() == FiniteElement::LIQUID_FLOW)
                        break;
                }
                // get index
                idx_k = pcs_tmp->GetElementValueIndex("PERMEABILITY");
                idx_n = pcs_tmp->GetElementValueIndex("POROSITY");

                // get values of n.
                n_new = pcs_tmp->GetElementValue(index, idx_n + 1);

                // calculate new permeability
                // k_rel(n) = n_rel^{19/6}
                // first relative porosity change
                n_rel = n_new / permeability_porosity_model_values[0];
                // then relative permeability change
                k_rel = pow(n_rel, 19.0 / 6.0);
                // finially permeability
                k_new = k_rel * permeability_porosity_model_values[1];
                // save new permeability
                m_pcs->SetElementValue(index, idx_k + 1, k_new);

                // now gives the newly calculated value to tensor[]
                tensor[0] = k_new;
            }
        }
        else if (permeability_model == 6)
        {  // HS: 11.2008, for Clement biomass colony clogging
            // if first time step, do nothing. otherwise,
            if (aktueller_zeitschritt > 1)
            {
                CRFProcess* pcs_tmp(NULL);
                for (size_t i = 0; i < pcs_vector.size(); i++)
                {
                    pcs_tmp = pcs_vector[i];
                    //	                if (
                    // m_pcs_tmp->pcs_type_name.compare("GROUNDWATER_FLOW") == 0
                    //||
                    //                                   m_pcs_tmp->pcs_type_name.compare("LIQUID_FLOW")
                    //                                   == 0) TF
                    if (pcs_tmp->getProcessType() ==
                            FiniteElement::GROUNDWATER_FLOW ||
                        pcs_tmp->getProcessType() == FiniteElement::LIQUID_FLOW)
                        break;
                }
                // get index
                idx_k = pcs_tmp->GetElementValueIndex("PERMEABILITY");
                idx_n = pcs_tmp->GetElementValueIndex("POROSITY");

                // get values of n.
                n_new = pcs_tmp->GetElementValue(index, idx_n + 1);

                // calculate new permeability
                // k_rel(n) = n_rel^{19/6}
                // first relative porosity change
                n_rel = n_new / permeability_porosity_model_values[0];

                // then relative permeability change
                const double temp(
                    (n_rel - permeability_porosity_model_values[0]) /
                    (1 - permeability_porosity_model_values[0]));
                k_rel =
                    permeability_porosity_model_values[2] *
                        MathLib::fastpow(
                            (n_rel - permeability_porosity_model_values[0]) /
                                (1 - permeability_porosity_model_values[0]),
                            3) +
                    (1 - permeability_porosity_model_values[2]) * temp * temp;
                // finially permeability
                k_new = k_rel * permeability_porosity_model_values[1];
                // save new permeability
                m_pcs->SetElementValue(index, idx_k + 1, k_new);

                // now gives the newly calculated value to tensor[]
                tensor[0] = k_new;
            }
        }
        else if (permeability_model == 7)
        {  // HS: 11.2008, for Clement biofilm clogging
            // if first time step, do nothing. otherwise,
            if (aktueller_zeitschritt > 1)
            {
                CRFProcess* pcs_tmp(NULL);
                for (size_t i = 0; i < pcs_vector.size(); i++)
                {
                    pcs_tmp = pcs_vector[i];
                    //	                if (
                    // m_pcs_tmp->pcs_type_name.compare("GROUNDWATER_FLOW") == 0
                    //|| m_pcs_tmp->pcs_type_name.compare("LIQUID_FLOW") == 0)
                    // TF
                    if (pcs_tmp->getProcessType() ==
                            FiniteElement::GROUNDWATER_FLOW ||
                        pcs_tmp->getProcessType() == FiniteElement::LIQUID_FLOW)
                        break;
                }
                // get index
                idx_k = pcs_tmp->GetElementValueIndex("PERMEABILITY");
                idx_n = pcs_tmp->GetElementValueIndex("POROSITY");

                // get values of n.
                n_new = pcs_tmp->GetElementValue(index, idx_n + 1);

                // calculate new permeability
                // k_rel(n) = n_rel^{19/6}
                // first relative porosity change
                n_rel = n_new / permeability_porosity_model_values[0];

                // then relative permeability change
                k_rel = (pow((n_rel - permeability_porosity_model_values[0]) /
                                 (1 - permeability_porosity_model_values[0]),
                             permeability_porosity_model_values[2]) +
                         permeability_porosity_model_values[3]) /
                        (1 + permeability_porosity_model_values[3]);
                // finially permeability
                k_new = k_rel * permeability_porosity_model_values[1];
                // save new permeability
                m_pcs->SetElementValue(index, idx_k + 1, k_new);

                // now gives the newly calculated value to tensor[]
                tensor[0] = k_new;
            }
        }
        else if ((permeability_model == 8) && (porosity_model == 13))
        {  // ABM 11.2010
            idx_k = m_pcs_tmp->GetElementValueIndex("PERMEABILITY");
            tensor[0] = m_pcs_tmp->GetElementValue(index, idx_k + 1);
        }
        // end of K-C
        // relationship-----------------------------------------------------------------------------------
    }
    // For heterogeneous values
    if (permeability_model == 2)
    {
        for (perm_index = 0; perm_index < m_pcs->m_msh->mat_names_vector.size();
             perm_index++)
            if (m_pcs->m_msh->mat_names_vector[perm_index].compare(
                    "PERMEABILITY") == 0)
                break;
        // end of getting the
        // index---------------------------------------------------------

        tensor[0] = _mesh->ele_vector[index]->mat_vector[perm_index];
        // CMCD
        // 01.09.2011 WW.  int edx =
        // m_pcs->GetElementValueIndex("PERMEABILITY"); CMCD 01.09.2011 WW.
        // m_pcs->SetElementValue(index,edx,tensor[0]);
    }

    switch (geo_dimension)
    {
        case 1:  // 1-D
            // HS: tensor[0] already set, doing nothing here;
            break;
        case 2:  // 2-D
            if (permeability_tensor_type == 0)
            {
                // tensor[0] = permeability_tensor[0]; // HS: done already;
                tensor[1] = 0.0;
                tensor[2] = 0.0;
                // tensor[3] = permeability_tensor[0];
                tensor[3] = tensor[0];  // HS: use the existing value;

            }
            else if (permeability_tensor_type == 1)
            {
                if ((permeability_model == 8) && (porosity_model == 13))
                {  // AB 11.2010
                    idx_k = m_pcs_tmp->GetElementValueIndex("PERMEABILITY");
                    tensor[0] = m_pcs_tmp->GetElementValue(index, idx_k + 1);
                    tensor[1] = 0.0;
                    tensor[2] = 0.0;
                    idx_k = m_pcs_tmp->GetElementValueIndex("PERMEABILITY_YY");
                    tensor[3] = m_pcs_tmp->GetElementValue(index, idx_k + 1);
                }
                else if (permeability_model == 2)
                {
                    // tensor[0] = already set
                    tensor[1] = 0.0;
                    tensor[2] = 0.0;
                    tensor[3] =
                        _mesh->ele_vector[index]->mat_vector[perm_index + 1];
                }
                else
                {
                    tensor[0] = permeability_tensor[0];
                    tensor[1] = 0.0;
                    tensor[2] = 0.0;
                    tensor[3] = permeability_tensor[1];
                }
            }
            else if (permeability_tensor_type == 2)
            {
                if (permeability_model == 2)
                {
                    // tensor[0] = already se
                    tensor[1] =
                        _mesh->ele_vector[index]->mat_vector[perm_index + 1];
                    tensor[2] =
                        _mesh->ele_vector[index]->mat_vector[perm_index + 2];
                    tensor[3] =
                        _mesh->ele_vector[index]->mat_vector[perm_index + 3];
                }
                else
                {
                    tensor[0] = permeability_tensor[0];
                    tensor[1] = permeability_tensor[1];
                    tensor[2] = permeability_tensor[2];
                    tensor[3] = permeability_tensor[3];
                }
            }
            break;
        case 3:  // 3-D
            if (permeability_tensor_type == 0)
            {
                // tensor[0] = permeability_tensor[0]; // HS: not needed.
                // already done before
                tensor[1] = 0.0;
                tensor[2] = 0.0;
                tensor[3] = 0.0;
                // tensor[4] = permeability_tensor[0]; // HS: using the line
                // below instead;
                tensor[4] = tensor[0];  // using the existing value
                tensor[5] = 0.0;
                tensor[6] = 0.0;
                tensor[7] = 0.0;
                // tensor[8] = permeability_tensor[0]; // HS: using the line
                // below instead;
                tensor[8] = tensor[0];  // using the existing value
                                        // HS: this is not needed any
                                        // more--------------------------------
                                        // if(permeability_model==2) {
                // SB 4218	tensor[0] = GetHetValue(index,"permeability");
                //      tensor[0] =
                //      m_pcs->m_msh->ele_vector[index]->mat_vector[perm_index];
                //      tensor[4] = tensor[0];
                //      tensor[8] = tensor[0];
                // }
                // end of comment out
                // section-------------------------------------
            }
            else if (permeability_tensor_type == 1)
            {
                if ((permeability_model == 8) && (porosity_model == 13))
                {  // AB 11.2010
                    idx_k = m_pcs_tmp->GetElementValueIndex("PERMEABILITY");
                    tensor[0] = m_pcs_tmp->GetElementValue(index, idx_k + 1);
                    tensor[1] = 0.0;
                    tensor[2] = 0.0;
                    tensor[3] = 0.0;
                    idx_k = m_pcs_tmp->GetElementValueIndex("PERMEABILITY_YY");
                    tensor[4] = m_pcs_tmp->GetElementValue(index, idx_k + 1);
                    tensor[5] = 0.0;
                    tensor[6] = 0.0;
                    tensor[7] = 0.0;
                    idx_k = m_pcs_tmp->GetElementValueIndex("PERMEABILITY_ZZ");
                    tensor[8] = m_pcs_tmp->GetElementValue(index, idx_k + 1);
                }
                else if (permeability_model == 2)
                {
                    // tensor[0] = already set
                    tensor[1] = tensor[2] = tensor[3] = 0.0;
                    tensor[4] =
                        _mesh->ele_vector[index]->mat_vector[perm_index + 1];
                    tensor[5] = tensor[6] = tensor[7] = 0.0;
                    tensor[8] =
                        _mesh->ele_vector[index]->mat_vector[perm_index + 2];
                }
                else
                {
                    tensor[0] = permeability_tensor[0];
                    tensor[1] = 0.0;
                    tensor[2] = 0.0;
                    tensor[3] = 0.0;
                    tensor[4] = permeability_tensor[1];
                    tensor[5] = 0.0;
                    tensor[6] = 0.0;
                    tensor[7] = 0.0;
                    tensor[8] = permeability_tensor[2];
                }
            }
            else if (permeability_tensor_type == 2)
            {
                if (permeability_model == 2)
                {
                    // tensor[0] = already se
                    tensor[1] =
                        _mesh->ele_vector[index]->mat_vector[perm_index + 1];
                    tensor[2] =
                        _mesh->ele_vector[index]->mat_vector[perm_index + 2];
                    tensor[3] =
                        _mesh->ele_vector[index]->mat_vector[perm_index + 3];
                    tensor[4] =
                        _mesh->ele_vector[index]->mat_vector[perm_index + 4];
                    tensor[5] =
                        _mesh->ele_vector[index]->mat_vector[perm_index + 5];
                    tensor[6] =
                        _mesh->ele_vector[index]->mat_vector[perm_index + 6];
                    tensor[7] =
                        _mesh->ele_vector[index]->mat_vector[perm_index + 7];
                    tensor[8] =
                        _mesh->ele_vector[index]->mat_vector[perm_index + 8];
                }
                else
                {
                    tensor[0] = permeability_tensor[0];
                    tensor[1] = permeability_tensor[1];
                    tensor[2] = permeability_tensor[2];
                    tensor[3] = permeability_tensor[3];
                    tensor[4] = permeability_tensor[4];
                    tensor[5] = permeability_tensor[5];
                    tensor[6] = permeability_tensor[6];
                    tensor[7] = permeability_tensor[7];
                    tensor[8] = permeability_tensor[8];
                }
            }
            break;
    }
    return tensor;
}

//------------------------------------------------------------------------
// 12.(i) PERMEABILITY_FUNCTION_DEFORMATION
//------------------------------------------------------------------------
//
//------------------------------------------------------------------------
// 12.(ii) PERMEABILITY_FUNCTION_PRESSURE
//------------------------------------------------------------------------
// WX: implementation of ths permeability_function_pressure. 1. version only for
// multi_phase_flow. 05.2010
double CMediumProperties::PermeabilityFunctionPressure(long /*index*/,
                                                       double PG2)
{
    int gueltig;  // WX: for function GetCurveValue(). 11.05.2010
    double fac_perm_pressure = 1;
    // WX: permeability as function of gas pressure. 11.05.2010
    switch (permeability_pressure_model)
    {
        case 10:  // WX: case 10, factor directly calculated from curve. 05.2010
            if (m_pcs->getProcessType() == FiniteElement::MULTI_PHASE_FLOW)
                // WX: now it's only works for Multi_Phase_Flow. 05.2010
                fac_perm_pressure =
                    GetCurveValue((int)permeability_pressure_model_values[0], 0,
                                  PG2, &gueltig);
            // WX: PG2 = gas pressue. 11.05.2010
            break;
        default:
            break;
    }
    return fac_perm_pressure;
}
// WX: Permeability_function_stain: 05.2010
double CMediumProperties::PermeabilityFunctionStrain(
    long index, int nnodes, CFiniteElementStd* h_fem)  // WW:02.08.2010
{
    double fac_perm_strain = 1.;
    // WW int ele_index,
    int gueltig;
    // WW CRFProcessDeformation *dm_pcs = (CRFProcessDeformation *) this;

    // get plas strain and volume strain
    int idStrainP;
    double strainp_nodes[20] = {0.};
    double strainp = 0.;
    idStrainP = h_fem->dm_pcs->GetNodeValueIndex("STRAIN_PLS");
    for (int i = 0; i < nnodes; i++)
        strainp_nodes[i] = h_fem->dm_pcs->GetNodeValue(
            h_fem->dm_pcs->m_msh->ele_vector[index]->getNodeIndices()[i],
            idStrainP);
    strainp = h_fem->interpolate(strainp_nodes);

    double strain_temp[3] = {0}, vol_strain_temp = 0;
    int idx_temp[3];
    int dim = m_pcs->m_msh->GetCoordinateFlag() / 10;
    if (dim == 2)
        if (h_fem->axisymmetry)
            dim = 3;
    idx_temp[0] = h_fem->dm_pcs->GetNodeValueIndex("STRAIN_XX");
    idx_temp[1] = h_fem->dm_pcs->GetNodeValueIndex("STRAIN_YY");
    idx_temp[2] = h_fem->dm_pcs->GetNodeValueIndex("STRAIN_ZZ");
    double strain_nodes[20] = {0};
    for (int j = 0; j < dim; j++)
    {
        for (int i = 0; i < nnodes; i++)
            strain_nodes[i] = h_fem->dm_pcs->GetNodeValue(
                h_fem->dm_pcs->m_msh->ele_vector[index]->getNodeIndices()[i],
                idx_temp[j]);
        strain_temp[j] = h_fem->interpolate(strain_nodes);
    }
    for (int j = 0; j < dim; j++)
        vol_strain_temp += strain_temp[j];

    switch (permeability_strain_model)
    {
        case 1:
        {
            fac_perm_strain = GetCurveValue(permeability_strain_model_value[0],
                                            0, vol_strain_temp, &gueltig);
            if (fac_perm_strain <= 0.)
                fac_perm_strain = 1.;
            break;
        }
        case 2:  // equivalent plasical strain
        {
            fac_perm_strain = GetCurveValue(permeability_strain_model_value[0],
                                            0, strainp, &gueltig);
            if (fac_perm_strain <= 0.)
                fac_perm_strain = 1.;
            break;
        }
        case 3:  // if StrainP>0, factor=f(StrainP), else
                 // factor=f(strain_Volume)
        {
            if (strainp > 0)
                fac_perm_strain = GetCurveValue(
                    permeability_strain_model_value[1], 0, strainp, &gueltig);
            else
            {
                fac_perm_strain =
                    GetCurveValue(permeability_strain_model_value[0], 0,
                                  vol_strain_temp, &gueltig);
            }
            if (fac_perm_strain <= 0.)
                fac_perm_strain = 1.;
            break;
        }
        case 4:  // factor = f(strainP+strain_Volume)
        {
            double tmpfkt = 1.;
            if (strainp > 0.)
            {
                tmpfkt = GetCurveValue(permeability_strain_model_value[1], 0,
                                       strainp, &gueltig);
                if (tmpfkt < 1.)
                    tmpfkt = 1.;
            }
            fac_perm_strain = GetCurveValue(permeability_strain_model_value[0],
                                            0, vol_strain_temp, &gueltig);
            if (fac_perm_strain <= 0.)
                fac_perm_strain = 1.;
            fac_perm_strain *= tmpfkt;
            break;
        }
        case 5:
        {
            double threshold = 0.;
            threshold = permeability_strain_model_value[0];

            // TODO: Error index out of bounds
            if (permeability_strain_model_value[3] > MKleinsteZahl)
                threshold = GetCurveValue(permeability_strain_model_value[3], 0,
                                          vol_strain_temp, &gueltig);
            if (vol_strain_temp <= threshold)
                fac_perm_strain = 1 - permeability_strain_model_value[1] *
                                          (threshold - vol_strain_temp);
            else
                fac_perm_strain = 1 + permeability_strain_model_value[2] *
                                          (vol_strain_temp - threshold);
            fac_perm_strain =
                MRange(permeability_strain_model_value[4], fac_perm_strain,
                       permeability_strain_model_value[5]);
            break;
        }
        default:
            break;
    }
    return fac_perm_strain;
}
//------------------------------------------------------------------------
// 12.(iv) PERMEABILITY_FUNCTION_EFFSTRESS
//------------------------------------------------------------------------
// AS: permeability as function of effective stress for liquid flow. 08.2012
double CMediumProperties::PermeabilityFunctionEffStress(
    long index, int nnodes, CFiniteElementStd* h_fem)
{
    int i, j, size;
    double perm_stress = 1.0;
    // calculate principal effective stress
    double stress[6] = {0.}, prin_str[6] = {0.}, prin_dir[9] = {0.};
    int stress_index[6];
    // model dimension
    int dim = h_fem->dm_pcs->m_msh->GetCoordinateFlag() / 10;
    size = 6;
    if (dim == 2)
    {
        size = 4;
        if (h_fem->axisymmetry)
            dim = 3;
    }
    stress_index[0] = h_fem->dm_pcs->GetNodeValueIndex("STRESS_XX");
    stress_index[1] = h_fem->dm_pcs->GetNodeValueIndex("STRESS_YY");
    stress_index[2] = h_fem->dm_pcs->GetNodeValueIndex("STRESS_ZZ");
    stress_index[3] = h_fem->dm_pcs->GetNodeValueIndex("STRESS_XY");
    if (size == 6)
    {
        stress_index[4] = h_fem->dm_pcs->GetNodeValueIndex("STRESS_XZ");
        stress_index[5] = h_fem->dm_pcs->GetNodeValueIndex("STRESS_YZ");
    }
    double stress_nodes[20] = {0.};

    switch (permeability_effstress_model)
    {
        // AS: case 1, permeability directly calculated from curve. 16.08.2012
        case 1:
            for (j = 0; j < size; j++)
            {
                for (i = 0; i < nnodes; i++)
                    stress_nodes[i] = h_fem->dm_pcs->GetNodeValue(
                        h_fem->dm_pcs->m_msh->ele_vector[index]
                            ->getNodeIndices()[i],
                        stress_index[j]);
                stress[j] = h_fem->interpolate(stress_nodes);
            }
            h_fem->SolidProp->CalPrinStrDir(stress, prin_str, prin_dir, dim);
            // permeability from curve with minimum (absolute value) principal
            // effective stress as input
            perm_stress = GetCurveValue(
                (int)permeability_effstress_model_value[0], 0, prin_str[0], &i);
            break;
        default:
            break;
    }
    return perm_stress;
}
//------------------------------------------------------------------------
// 12.(i) PERMEABILITY_FUNCTION_SATURATION
//------------------------------------------------------------------------
//------------------------------------------------------------------------
// 12.(vi) PERMEABILITY_FUNCTION_POROSITY
//------------------------------------------------------------------------

//---------------------------------------------------------------------------------
// 12.(vii) PERMEABILITY_FUNCTION_FRAC_APERTURE
// RFW 07/2005
//---------------------------------------------------------------------------------
//------------------------------------------------------------------------
// 13. CAPILLARY_PRESSURE_FUNCTION
//------------------------------------------------------------------------

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/1997 CT Erste Version
   09/1998 CT Kurven ueber Schluesselwort #CURVES
   08/1999 CT Erweiterung auf n-Phasen begonnen
   08/2004 OK Template for MMP implementation
   01/2004 OK Mode 3, given saturation
   05/2008 WW Brook & Corey
   03/2012 JT All new
   last modification:
**************************************************************************/
double CMediumProperties::CapillaryPressureFunction(
    const double wetting_saturation)
{
    double pc, pb, sl, slr, slm, se, m;
    int gueltig;
    sl = wetting_saturation;
    //
    switch (capillary_pressure_model)
    {
        default:
            ScreenMessage(
                "Error in CFluidProperties::CapillaryPressure: no valid "
                "material model.\n");
            exit(0);
            break;
        //
        case 0:  // k=f(x)
            pc = GetCurveValue((int)capillary_pressure_values[0], 0, sl,
                               &gueltig);
            break;
        //
        case 1:  // Constant capillary pressure for ps models
            pc = capillary_pressure_values[0];
            break;
        //
        case 2:  // Constant saturation for pp models (for WX, from JT) (MUST BE
                 // A PP MODEL, SO WON'T COME HERE)
            ScreenMessage("ERROR: in CFluidProperties::CapillaryPressure:");
            ScreenMessage(
                "Constant saturation is not possible for a PS model "
                "(PwSnw).\n");
            exit(0);
            break;
        //
        case 4:  // van Genuchten
            pb = capillary_pressure_values[0];
            slr = capillary_pressure_values[1];
            slm = capillary_pressure_values[2];
            m = capillary_pressure_values[3];  // always <= 1.0.  Input is
                                               // exponent = 1 / (1-lambda)
            //
            // Convert alpha to entry pressure?
            if (entry_pressure_conversion)
                pb = (mfp_vector[0]->Density() * 9.81) / pb;
            //
            // sl  = MRange(slr, sl, slm);
            sl = MRange(slr + DBL_EPSILON, sl, slm - DBL_EPSILON);
            se = (sl - slr) / (slm - slr);
            pc = pb * pow(pow(se, (-1.0 / m)) - 1.0, 1.0 - m);
            pc = MRange(DBL_EPSILON, pc, capillary_pressure_values[4]);
            break;
        //
        case 6:  //  Brook & Corey
            pb = capillary_pressure_values[0];
            slr = capillary_pressure_values[1];
            slm = capillary_pressure_values[2];
            m = capillary_pressure_values[3];  // always >= 1.0
            //
            sl = MRange(slr, sl, slm);
            se = (sl - slr) / (slm - slr);
            pc = pb * pow(se, -1.0 / m);
            pc = MRange(DBL_EPSILON, pc, capillary_pressure_values[4]);
            break;
    }
    return pc;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2002 CT CECalcSatuFromCap
        SB Extensions case 4 and case 5
   02/2005 OK CMediumProperties function
   08/2005 WW Remove interploation
   03/2012 JT All new.
   Last modified:
**************************************************************************/
double CMediumProperties::SaturationCapillaryPressureFunction(
    const double capillary_pressure)
{
    double se, sl, slr, slm, m, pb, pc;
    int gueltig;
    pc = capillary_pressure;
    //
    // Get Se
    switch (capillary_pressure_model)
    {
        default:
            ScreenMessage(
                "Error in "
                "CFluidProperties::SaturationCapillaryPressureFunction: no "
                "valid material model.\n");
            exit(0);
            break;
        //
        case 0:  // k=f(x)
            // Note: Use Inverse curve value because .rfd file has x=S, y=Pc as
            // columns, and Pc is our input.
            sl = GetCurveValueInverse((int)capillary_pressure_values[0], 0, pc,
                                      &gueltig);
            sl = MRange(capillary_pressure_values[1] + DBL_EPSILON, sl,
                        capillary_pressure_values[2] - DBL_EPSILON);
            break;
        //
        case 1:  // Constant capillary pressure for ps models
            ScreenMessage(
                "ERROR: in "
                "CFluidProperties::SaturationCapillaryPressureFunction:");
            ScreenMessage(
                "Constant capillary pressure is not possible for a pressure "
                "model (PcPnw, PwPnw, or Richards).\n");
            exit(0);
            break;
        //
        case 2:  // Constant saturation for pp models (for WX, from JT)
            sl = capillary_pressure_values[0];
            break;
        //
        case 4:  // van Genuchten
            pb = capillary_pressure_values[0];
            slr = capillary_pressure_values[1];
            slm = capillary_pressure_values[2];
            m = capillary_pressure_values[3];  // always <= 1.0.  Input is
                                               // exponent = 1 / (1-lambda)
            //
            // Convert alpha to entry pressure?
            if (entry_pressure_conversion)
                pb = (mfp_vector[0]->Density() * 9.81) / pb;
            //
            if (pc < 0.0)
                pc = 0.0;
            se = pow(pc / pb, 1.0 / (1.0 - m)) + 1.0;
            se = pow(se, -m);
            sl = se * (slm - slr) + slr;
            sl = MRange(slr + DBL_EPSILON, sl, slm - DBL_EPSILON);
            break;
        //
        case 6:  //  Brook & Corey
            pb = capillary_pressure_values[0];
            slr = capillary_pressure_values[1];
            slm = capillary_pressure_values[2];
            m = capillary_pressure_values[3];  // always >= 1.0
            if (pc < pb)
                pc = pb;
            se = pow(pc / pb, -m);
            sl = se * (slm - slr) + slr;
            sl = MRange(slr + DBL_EPSILON, sl, slm - DBL_EPSILON);
            break;
        case 10:  //  unconfined 3D GW.  5.3.07 JOD
            // MW: remove comment to provide variables, not constants to
            // PermeabilitySaturationFunction
            pb = capillary_pressure_values[0];
            slr = capillary_pressure_values[1];
            if (pc > 0)
            {
                sl = max(0., 1 - (pc / pb));
                sl = pow(sl, 2 * (1 - sl));
                sl = slr + (1 - slr) * sl;
            }
            else
                sl = 1;
            break;
    }
    return sl;
}

/**************************************************************************
   FEMLib-Method:
   Task: returns dSw/dPc
   Programing:
   08/1999 CT Erste Version (MMTM0699GetSaturationPressureDependency)
   05/2001 OK Verallgemeinerung
   02/2005 OK CMediumProperties function
   08/2005 WW Remove interpolation
   03/2007 WW Analytical solution:
   02/2008 PCH Brooks-Corey dPc/dSw added
   03/2012 JT All new. Added van Genuchten, Brooks Corey, fixed iterative, etc.
           Plus, input should be Pc, not Sw.
   Last modified:
**************************************************************************/

/* JT: No longer used. But it is accurate. Use PressureSaturationDependency()
with "invert" = true

double CMediumProperties::SaturationPressureDependency(const double
capillary_pressure, bool allow_zero)
{
    double dsdp,v1,v2,pc,pb,slr,slm,m;
    int gueltig;
    pc = capillary_pressure;
    //
    if(allow_zero && pc < DBL_EPSILON) // If we allow dSw/dPc = 0.0 (i.e.
Richard's flow) return 0.0;
    //
    switch(capillary_pressure_model)
    {
        default:
            ScreenMessage("Error in
CFluidProperties::SaturationPressureDependency: no valid material model.\n");
            exit(0);
            break;

        case 0: // Curve value
            // Note: Use Inverse curve value because .rfd file has x=S, y=Pc as
columns, and Pc is our input. dsdp =
GetCurveInverseDerivative((int)capillary_pressure_values[0],1,pc,&gueltig);
            break;

        case 1: //  Pc = CONSTANT
            dsdp = 0.0;
            break;

        case 2: //  Sw = CONSTANT
            dsdp = 0.0;
            break;

        case 4: //  Van Genuchten: 01.3.2007 WW  // 05.2010 JT.
            pb  = capillary_pressure_values[0];
            slr = capillary_pressure_values[1];
            slm = capillary_pressure_values[2];
            m   = capillary_pressure_values[3];			// always <= 1.0. Input
is exponent = 1 / (1-lambda) pc  =
MRange(FLT_EPSILON,pc,capillary_pressure_values[4]);
            //
            // Convert alpha to entry pressure?
            if(entry_pressure_conversion)
                pb = (mfp_vector[0]->Density()*9.81)/pb;
            //
            v1 = pow((pc/pb),(1.0/(1.0-m)));
            v2 = pow((1.0+v1),(-1.0-m));
            dsdp = (m*v1*v2*(slm-slr)) / ((m-1.0)*pc);
            break;

        case 6: //  Brooks & Corey. 10/2010 JT
            pb = capillary_pressure_values[0];
            pc  = MRange(FLT_EPSILON,pc,capillary_pressure_values[4]);
            slr = capillary_pressure_values[1];
            slm = capillary_pressure_values[2];
            m   = capillary_pressure_values[3];		// always >= 1.0
            //
            v1 = pow((pc/pb),-m);
            dsdp = (m*v1*(slr - slm)) / pc;
            break;
    }
    //
    return dsdp;
}
*/

/**************************************************************************
   FEMLib-Method:
   Task: returns dPc/dSw
   "invert" = false. Return dPc/dSw
   "invert" = true.  Return dSw/dPc
   Programing:
   03/2009 PCH Brooks-Corey dPc/dSw added
   03/2012 JT: All new. + van Genuchten, use curve derivative, etc.
   Last modified:
**************************************************************************/
double CMediumProperties::PressureSaturationDependency(
    const double wetting_saturation, bool invert)
{
    double dpds, dsdp, v1, v2, pb, sl, slr, slm, m, lim, ds, dpc;
    int gueltig;
    sl = wetting_saturation;
    //
    switch (capillary_pressure_model)
    {
        default:
            ScreenMessage(
                "Error in CFluidProperties::PressureSaturationDependency: no "
                "valid material model.\n");
            exit(0);
            break;

        case 0:  // curve value
            sl = MRange(capillary_pressure_values[1] + DBL_EPSILON, sl,
                        capillary_pressure_values[2] - DBL_EPSILON);
            dpds = GetCurveDerivative((int)capillary_pressure_values[0], 1, sl,
                                      &gueltig);
            break;

        case 1:  //  Pc = CONSTANT
            return 0.0;

        case 2:  //  Sw = CONSTANT
            return 0.0;

        case 4:  //  Van Genuchten: 01.3.2007 WW  // 05.2010 JT.
            pb = capillary_pressure_values[0];
            slr = capillary_pressure_values[1];
            slm = capillary_pressure_values[2];
            m = capillary_pressure_values[3];  // always <= 1.0.  Input is
                                               // exponent = 1 / (1-lambda)
            sl = MRange(slr + DBL_EPSILON, sl,
                        slm - DBL_EPSILON);  // (infinity also occurs at
                                             // sl=slmax for van Genuchten)
            //
            // Convert alpha to entry pressure?
            if (entry_pressure_conversion)
                pb = (mfp_vector[0]->Density() * 9.81) / pb;
            //
            // Get dPc/dSw
            v1 = pow(((sl - slr) / (slm - slr)), (-1.0 / m));
            v2 = pow(v1 - 1.0, -m);
            dpds = (pb * (m - 1.0) * v1 * v2) / (m * (sl - slr));
            break;

        case 6:  //  Brooks & Corey. 10/2010 JT
            pb = capillary_pressure_values[0];
            slr = capillary_pressure_values[1];
            slm = capillary_pressure_values[2];
            m = capillary_pressure_values[3];  // always >= 1.0
            sl = MRange(slr + DBL_EPSILON, sl,
                        slm);  // No upper bound needed for B&C
            //
            // Get dPc/dSw
            v1 = pow(((sl - slr) / (slm - slr)), (-1.0 / m));
            dpds = (pb * v1) / (m * (slr - sl));
            break;

        case 10:       //  unconfined 3D GW  6/2012 JOD
            return 0;  // set phi/(rho * g) in storage
        case 99:  // The old iterative method. Should anyone need it (but it is
                  // somewhat inaccurate at low and high saturations)
            ds = 1.0e-2;
            do
            {
                ds /= 10.;
                v1 = CapillaryPressureFunction(sl - ds);
                v2 = CapillaryPressureFunction(sl + ds);
                dpc = v1 - v2;
            } while ((ds > DBL_EPSILON) && (v2 < DBL_EPSILON / 100.));
            if (((v1 > DBL_EPSILON) || (v2 > DBL_EPSILON)) &&
                (dpc > DBL_EPSILON))
                dsdp = (-2.0 * ds) / dpc;
            else
                dsdp = 0.0;
            //
            dpds = 1.0 / dsdp;
            break;
    }
    //
    // Note: Sw and Pc are inversly related.. i.e. dsdp and dpds are always <=
    // 0.0
    lim = (-1.0) / DBL_EPSILON;
    //
    if (invert)
    {  // Return dSw/dPc
        dsdp = 1.0 / dpds;
        if (dsdp < lim)
            dsdp = lim;
        return dsdp;
    }
    else
    {  // Return dPc/dSw
        if (dpds < lim)
            dpds = lim;
        return dpds;
    }
}

/*************************************************************************************************************
   14.$MASSDISPERSION
 ************************************************************************************************************/

/*************************************************************************************************************
   Density
 ************************************************************************************************************/
/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   11/2004 OK Implementation
**************************************************************************/
double CMediumProperties::Density(long element, double* gp, double theta)
{
    // OK411
    gp = gp;

    int no_phases = (int)mfp_vector.size();
    double density = 0.0;
    int i;
    CFluidProperties* m_mfp = NULL;
    // OK411 CSolidProperties* m_msp = NULL;
    char saturation_name[15];
    if (no_phases == 1)
    {
        m_mfp = mfp_vector[0];
        density = Porosity(element, theta) * m_mfp->Density();
    }
    else
        for (i = 0; i < no_phases; i++)
        {
            m_mfp = mfp_vector[i];
            sprintf(saturation_name, "SATURATION%i", i + 1);
            // OK411 density += Porosity(element,theta) * m_mfp->Density() *
            // PCSGetELEValue(element,gp,theta,saturation_name);
        }
    /*OK411
       long group = ElGetElementGroupNumber(element);
       m_msp = msp_vector[group];
       density += (1.-Porosity(element,theta))*fabs(m_msp->Density());
     */
    return density;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   01/2005 OK Implementation
   last modified:
**************************************************************************/
void MMPDelete()
{
    long i;
    int no_mmp = (int)mmp_vector.size();
    for (i = 0; i < no_mmp; i++)
        delete mmp_vector[i];
    mmp_vector.clear();
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   06/2005 OK Implementation
   07/2005 WW Change due to geometry element object applied
   last modified:
**************************************************************************/
void MMP2PCSRelation(CRFProcess* m_pcs)
{
    CMediumProperties* m_mmp = NULL;
    MeshLib::CElem* m_ele = NULL;
    if (m_pcs->m_msh)
        for (long i = 0; i < (long)m_pcs->m_msh->ele_vector.size(); i++)
        {
            m_ele = m_pcs->m_msh->ele_vector[i];
            if (m_ele->GetPatchIndex() < mmp_vector.size())
            {
                m_mmp = mmp_vector[m_ele->GetPatchIndex()];
                m_mmp->m_pcs = m_pcs;
            }
        }
    else
        for (int j = 0; j < (int)mmp_vector.size(); j++)
        {
            m_mmp = mmp_vector[j];
            m_mmp->m_pcs = m_pcs;
        }
}

/**************************************************************************
   PCSLib-Function
   Liest zu jedem Knoten einen Wert der Permeabilität ein.
   Identifikation über Koordinaten
   Programing:
   10/2003     SB  First Version
   01/2004     SB  2. Version
   01/2005 OK Check MAT groups //OK41
   06/2005 MB msh, loop over mmp groups
   09/2005 MB EleClass
   //SB/MB ? member function of CFEMesh / CMediumProperties
   11/2005 OK GEOMETRY_AREA
   04/2006 CMCD Constant area
**************************************************************************/
// MMPGetHeterogeneousFields
void GetHeterogeneousFields()
{
    // OK411 int ok=0;
    // OK411 char* name_file=NULL;
    CMediumProperties* m_mmp = NULL;
    //----------------------------------------------------------------------
    // File handling
    string file_path;
    string file_path_base_ext;

    //----------------------------------------------------------------------
    // Tests
    if (mmp_vector.size() == 0)
        return;
    //----------------------------------------------------------------------
    // Schleife über alle Gruppen
    for (int i = 0; i < (int)mmp_vector.size(); i++)
    {
        m_mmp = mmp_vector[i];
        //....................................................................
        // For Permeability
        if (m_mmp->permeability_file.size() > 0)
        {
            // OK name_file = (char *) m_mmp->permeability_file.data();
            // OK if(name_file != NULL)
            // OK ok = FctReadHeterogeneousFields(name_file,m_mmp);

            // WW file_path_base_ext = file_path + m_mmp->permeability_file;
            // WW
            m_mmp->SetDistributedELEProperties(m_mmp->permeability_file);
            m_mmp->WriteTecplotDistributedProperties();
        }
        //....................................................................
        // For Porosity
        if (m_mmp->porosity_file.size() > 0)
        {
            // CB name_file = (char *) m_mmp->porosity_file.data();
            // CB if(name_file != NULL)
            // CB  ok = FctReadHeterogeneousFields(name_file,m_mmp);
            // file_path_base_ext = file_path + m_mmp->porosity_file;
            // m_mmp->SetDistributedELEProperties(file_path_base_ext); // CB
            // Removed bugs in this function CB Removed bugs in this function
            m_mmp->SetDistributedELEProperties(m_mmp->porosity_file);
            m_mmp->WriteTecplotDistributedProperties();
        }
        //....................................................................
        // GEOMETRY_AREA
        if (m_mmp->geo_area_file.size() > 0)
        {
            file_path_base_ext = file_path + m_mmp->geo_area_file;
            m_mmp->SetDistributedELEProperties(file_path_base_ext);
            m_mmp->WriteTecplotDistributedProperties();
        }
        // NW    else m_mmp->SetConstantELEarea(m_mmp->geo_area,i);
        //....................................................................
    }
    //----------------------------------------------------------------------
}

/**************************************************************************
   PCSLib-Method:
   Programing:
   04/2006 CMCD Implementation
**************************************************************************/
void CMediumProperties::SetConstantELEarea(double area, int group)
{
    long i, j, ele_group;
    long no_ele;
    int no_processes = (int)pcs_vector.size();
    if (area != 1.0)
        for (i = 0; i < no_processes; i++)
        {
            _mesh = FEMGet(
                convertProcessTypeToString(pcs_vector[i]->getProcessType()));
            if (!_mesh)
                return;  // WW
            no_ele = (long)_mesh->ele_vector.size();
            for (j = 0; j < no_ele; j++)
            {
                ele_group = _mesh->ele_vector[j]->GetPatchIndex();
                if (ele_group == group)
                    _mesh->ele_vector[j]->SetFluxArea(area);
            }
        }
}

/**************************************************************************
   PCSLib-Method:
   Programing:
   11/2005 OK Implementation
**************************************************************************/
void CMediumProperties::SetDistributedELEProperties(string file_name)
{
    string line_string, line1;
    string mmp_property_name;
    string mmp_property_dis_type;
    string mmp_property_mesh;
    MeshLib::CElem* m_ele_geo = NULL;
    bool element_area = false;
    long i, ihet;
    double mmp_property_value;
    double ddummy, conversion_factor = 1.0;  // init WW
    vector<double> xvals, yvals, zvals, mmpvals;
    vector<double> temp_store;
    double x, y, z, mmpv;
    std::stringstream in;
    // CB
    vector<double> garage;
    int por_index = 0;
    int vol_bio_index = 0;
    string outfile;
    int k;

    // default: scalar property
    unsigned int n_components{1};

    cout << " SetDistributedELEProperties: ";
    //----------------------------------------------------------------------
    // File handling
    ifstream mmp_property_file(file_name.data(), ios::in);
    if (!mmp_property_file.good())
    {
        cout << "Warning in CMediumProperties::SetDistributedELEProperties: no "
                "MMP property data"
             << "\n";
        return;
    }
    mmp_property_file.clear();
    mmp_property_file.seekg(0, ios::beg);
    //----------------------------------------------------------------------
    line_string = GetLineFromFile1(&mmp_property_file);
    if (!(line_string.find("#MEDIUM_PROPERTIES_DISTRIBUTED") != string::npos))
    {
        cout << "Keyword #MEDIUM_PROPERTIES_DISTRIBUTED not found"
             << "\n";
        return;
    }
    //----------------------------------------------------------------------
    while (!mmp_property_file.eof())
    {
        line_string = GetLineFromFile1(&mmp_property_file);
        if (line_string.find("STOP") != string::npos)
            return;
        if (line_string.empty())
        {
            cout << "Error in CMediumProperties::SetDistributedELEProperties - "
                    "no enough data sets"
                 << "\n";
            return;
        }
        //....................................................................
        if (line_string.find("$MSH_TYPE") != string::npos)
        {
            line_string = GetLineFromFile1(&mmp_property_file);
            mmp_property_mesh = line_string;
            _mesh = FEMGet(line_string);
            if (!_mesh)
            {
                cout << "CMediumProperties::SetDistributedELEProperties: no "
                        "MSH data"
                     << "\n";
                return;
            }
            continue;
        }
        //....................................................................
        if (line_string.find("$MMP_TYPE") != string::npos)
        {
            element_area = false;
            mmp_property_file >> mmp_property_name;
            cout << mmp_property_name << "\n";
            _mesh->mat_names_vector.push_back(mmp_property_name);
            if (mmp_property_name == "GEOMETRY_AREA")
                element_area = true;
            continue;
        }
        //....................................................................
        // optional parameter (default is zero)
        //
        if (line_string.find("$COMPONENTS") != string::npos)
        {
            mmp_property_file >> n_components;
            if (!mmp_property_name.empty())
                std::cerr << "Error in CMediumProperties::"
                          << "SetDistributedELEProperties:\n"
                          << "$MMP_TYPE hast to be set before $COMPONENTS"
                          << "\n";
            for (unsigned int i = 0; i < n_components - 1; ++i)
            {
                _mesh->mat_names_vector.push_back(mmp_property_name +
                                                  std::to_string(i + 1));
            }

            // check consistencies
            if (mmp_property_name.find("PERMEABILITY") != string::npos)
                switch (this->permeability_tensor_type)
                {
                    case 0:  // isotropic
                        assert(n_components == 1);
                        break;
                    case 1:  // orthotropic?
                        assert(n_components == this->geo_dimension);
                        break;
                    case 2:  // anisotropic
                        assert(
                            (this->geo_dimension == 2 && n_components == 4) ||
                            (this->geo_dimension == 3 && n_components == 9));
                        break;
                    default:
                        assert(!"PERMEABILITY_TENSOR_TYPE not recognized!");
                }
        }
        //....................................................................
        if (line_string.find("$DIS_TYPE") != string::npos)
        {
            mmp_property_file >> mmp_property_dis_type;
            continue;
        }
        //....................................................................
        if (line_string.find("$CONVERSION_FACTOR") != string::npos)
        {
            mmp_property_file >> conversion_factor;
            continue;
        }
        //....................................................................
        if (line_string.find("$DATA") != string::npos)
        {
            switch (mmp_property_dis_type[0])
            {
                case 'N':  // Next neighbour
                case 'G':  // Geometric mean
                           // add checks since implementation for 'G' is
                           // broken&unused
                    if (n_components > 1)
                    {
                        std::cerr << "more than one component of heterogeneous"
                                  << "fields only available with"
                                  << "$DIS_TYPE ELEMENT!\n";
                        assert(n_components == 1);  // must fail!
                    }
                    // Read in all values given, store in vectors for x, y, z
                    // and value
                    i = 0;
                    while (i == 0)
                    {
                        line1 = GetLineFromFile1(&mmp_property_file);
                        if (line1.find("STOP") != string::npos)
                            break;
                        in.str((string)line1);
                        in >> x >> y >> z >> mmpv;
                        in.clear();
                        mmpv *= conversion_factor;  // convert values
                        xvals.push_back(x);
                        yvals.push_back(y);
                        zvals.push_back(z);
                        mmpvals.push_back(mmpv);
                    }
                    // sort values to mesh
                    for (i = 0; i < (long)_mesh->ele_vector.size(); i++)
                    {
                        m_ele_geo = _mesh->ele_vector[i];
                        if (mmp_property_dis_type[0] == 'N')
                        {
                            // Search for all elements of the mesh, which is the
                            // nearest given value in the input file Return
                            // value ihet is the index of the het. val in the
                            // mmpval-vector
                            ihet = GetNearestHetVal2(i, _mesh, xvals, yvals,
                                                     zvals, mmpvals);
                            m_ele_geo->mat_vector.push_back(mmpvals[ihet]);
                        }
                        if (mmp_property_dis_type[0] == 'G')
                        {
                            mmpv = GetAverageHetVal2(i, _mesh, xvals, yvals,
                                                     zvals, mmpvals);
                            m_ele_geo->mat_vector.push_back(mmpv);
                        }
                    }
                    break;
                case 'E':  // Element data
                    for (i = 0; i < (long)_mesh->ele_vector.size(); i++)
                    {
                        m_ele_geo = _mesh->ele_vector[i];
                        mmp_property_file >> ddummy;
                        if (ddummy != i)
                        {
                            std::cout << "SetDistributedELEProperties:"
                                      << "element indices not contiguous!\n";
                        }
                        for (unsigned int c = 0; c < n_components; ++c)
                        {
                            mmp_property_file >> mmp_property_value;
                            m_ele_geo->mat_vector.push_back(mmp_property_value);
                        }
                        if (element_area)
                            _mesh->ele_vector[i]->SetFluxArea(
                                mmp_property_value);
                        if (line_string.empty())
                        {
                            cout << "Error in "
                                    "CMediumProperties::"
                                    "SetDistributedELEProperties - not enough "
                                    "data sets"
                                 << "\n";
                            return;
                        }
                    }
                    break;
                default:
                    cout << " Unknown interpolation option for the values!"
                         << "\n";
                    break;
            }
            continue;
        }
        //....................................................................
    }
    // CB now set VOL_MAT & VOL_BIO as heterogeneous values, if defined as model
    // 2 and het Porosity
    if ((mmp_property_name == "POROSITY") && (this->vol_bio_model == 2))
    {
        _mesh->mat_names_vector.push_back("VOL_BIO");
        for (i = 0; i < (long)_mesh->ele_vector.size(); i++)
        {
            m_ele_geo = _mesh->ele_vector[i];  // Get the element
            // Set the VOL_BIO value from mmp file input
            m_ele_geo->mat_vector.push_back(this->vol_bio);
        }
    }
    if ((mmp_property_name == "POROSITY") && (this->vol_mat_model == 2))
    {
        _mesh->mat_names_vector.push_back("VOL_MAT");
        // Get the porosity index
        for (por_index = 0; por_index < (int)_mesh->mat_names_vector.size();
             por_index++)
            if (_mesh->mat_names_vector[por_index].compare("POROSITY") == 0)
                break;
        // Get the vol_bio index
        for (vol_bio_index = 0;
             vol_bio_index < (int)_mesh->mat_names_vector.size();
             vol_bio_index++)
            if (_mesh->mat_names_vector[vol_bio_index].compare("VOL_BIO") == 0)
                break;
        for (i = 0; i < (long)_mesh->ele_vector.size(); i++)
        {
            m_ele_geo = _mesh->ele_vector[i];  // Get the element
            m_ele_geo->mat_vector.push_back(
                1 - m_ele_geo->mat_vector[por_index] -
                m_ele_geo->mat_vector[vol_bio_index]);
        }
    }
    //----------------------------------------------------------------------
    // Write sorted output file
    //----------------------------------------------------------------------
    // File handling

    // CB
    for (k = 0; k < (int)_mesh->mat_names_vector.size(); k++)
    {
        // file_name +="_sorted";
        outfile = _mesh->mat_names_vector[k] + "_sorted";
        ofstream mmp_property_file_out(outfile.data());
        if (!mmp_property_file_out.good())
        {
            cout << "Warning in "
                    "CMediumProperties::WriteDistributedELEProperties: no MMP "
                    "property data file to write to"
                 << "\n";
            return;
        }
        mmp_property_file_out << "#MEDIUM_PROPERTIES_DISTRIBUTED"
                              << "\n";
        mmp_property_file_out << "$MSH_TYPE"
                              << "\n"
                              << "  " << mmp_property_mesh << "\n";
        // mmp_property_file_out << "$MSH_TYPE" << "\n" << "  " <<
        // mmp_property_mesh << "\n"; mmp_property_file_out << "$MMP_TYPE" <<
        // "\n" << "  " << "PERMEABILITY" << "\n";
        mmp_property_file_out << "$MMP_TYPE"
                              << "\n"
                              << "  " << _mesh->mat_names_vector[k] << "\n";
        mmp_property_file_out << "$DIS_TYPE"
                              << "\n"
                              << "  "
                              << "ELEMENT"
                              << "\n";
        mmp_property_file_out << "$DATA"
                              << "\n";
        for (i = 0; i < (long)_mesh->ele_vector.size(); i++)
        {
            m_ele_geo = _mesh->ele_vector[i];
            mmp_property_file_out << i << "  " << m_ele_geo->mat_vector[k]
                                  << "\n";
        }
        mmp_property_file_out << "#STOP"
                              << "\n";
        mmp_property_file_out.close();
        //----------------------------------------------------------------------
    }
}

/**************************************************************************
   PCSLib-Method:
   Programing:
   11/2005 OK Implementation
**************************************************************************/
void CMediumProperties::WriteTecplotDistributedProperties()
{
    int j, k;
    long i;
    string element_type;
    string m_string = "MAT";
    double m_mat_prop_nod;
    //----------------------------------------------------------------------
    // Path
    string path;
    //--------------------------------------------------------------------
    // MSH
    MeshLib::CNode* m_nod = NULL;
    MeshLib::CElem* m_ele = NULL;
    if (!_mesh)
        return;
    //--------------------------------------------------------------------
    // File handling
    string mat_file_name = path + name + "_" + _mesh->pcs_name + "_PROPERTIES" +
                           TEC_FILE_EXTENSION;
    fstream mat_file(mat_file_name.data(), ios::trunc | ios::out);
    mat_file.setf(ios::scientific, ios::floatfield);
    mat_file.precision(12);
    if (!mat_file.good())
        return;
    mat_file.seekg(0L, ios::beg);
    //--------------------------------------------------------------------
    if ((long)_mesh->ele_vector.size() > 0)
    {
        m_ele = _mesh->ele_vector[0];
        switch (m_ele->GetElementType())
        {
            case MshElemType::LINE:
                element_type = "QUADRILATERAL";
                break;
            case MshElemType::QUAD:
                element_type = "QUADRILATERAL";
                break;
            case MshElemType::HEXAHEDRON:
                element_type = "BRICK";
                break;
            case MshElemType::TRIANGLE:
                element_type = "TRIANGLE";
                break;
            case MshElemType::TETRAHEDRON:
                element_type = "TETRAHEDRON";
                break;
            case MshElemType::PRISM:
                element_type = "BRICK";
                break;
            default:
                std::cerr
                    << "CMediumProperties::WriteTecplotDistributedProperties "
                       "MshElemType not handled"
                    << "\n";
        }
    }
    //--------------------------------------------------------------------
    // Header
    mat_file << "VARIABLES = X,Y,Z";
    for (j = 0; j < (int)_mesh->mat_names_vector.size(); j++)
        mat_file << "," << _mesh->mat_names_vector[j];
    mat_file << "\n";
    mat_file << "ZONE T = " << name << ", "
             << "N = " << (long)_mesh->nod_vector.size() << ", "
             << "E = " << (long)_mesh->ele_vector.size() << ", "
             << "F = FEPOINT"
             << ", "
             << "ET = " << element_type << "\n";
    //--------------------------------------------------------------------
    // Nodes
    for (i = 0; i < (long)_mesh->nod_vector.size(); i++)
    {
        m_nod = _mesh->nod_vector[i];
        double const* const pnt(m_nod->getData());
        mat_file << pnt[0] << " " << pnt[1] << " " << pnt[2];
        for (size_t j = 0; j < _mesh->mat_names_vector.size(); j++)
        {
            m_mat_prop_nod = 0.0;
            for (k = 0; k < (int)m_nod->getConnectedElementIDs().size(); k++)
            {
                m_ele = _mesh->ele_vector[m_nod->getConnectedElementIDs()[k]];
                m_mat_prop_nod += m_ele->mat_vector[j];
            }
            m_mat_prop_nod /= (int)m_nod->getConnectedElementIDs().size();
            mat_file << " " << m_mat_prop_nod;
        }
        mat_file << "\n";
    }
    //--------------------------------------------------------------------
    // Elements
    for (i = 0; i < (long)_mesh->ele_vector.size(); i++)
    {
        m_ele = _mesh->ele_vector[i];
        // OK if(m_ele->GetPatchIndex()==number) {
        switch (m_ele->GetElementType())
        {
            case MshElemType::LINE:
                mat_file << m_ele->getNodeIndices()[0] + 1 << " "
                         << m_ele->getNodeIndices()[1] + 1 << " "
                         << m_ele->getNodeIndices()[1] + 1 << " "
                         << m_ele->getNodeIndices()[0] + 1 << "\n";
                element_type = "QUADRILATERAL";
                break;
            case MshElemType::QUAD:
                mat_file << m_ele->getNodeIndices()[0] + 1 << " "
                         << m_ele->getNodeIndices()[1] + 1 << " "
                         << m_ele->getNodeIndices()[2] + 1 << " "
                         << m_ele->getNodeIndices()[3] + 1 << "\n";
                element_type = "QUADRILATERAL";
                break;
            case MshElemType::HEXAHEDRON:
                mat_file << m_ele->getNodeIndices()[0] + 1 << " "
                         << m_ele->getNodeIndices()[1] + 1 << " "
                         << m_ele->getNodeIndices()[2] + 1 << " "
                         << m_ele->getNodeIndices()[3] + 1 << " "
                         << m_ele->getNodeIndices()[4] + 1 << " "
                         << m_ele->getNodeIndices()[5] + 1 << " "
                         << m_ele->getNodeIndices()[6] + 1 << " "
                         << m_ele->getNodeIndices()[7] + 1 << "\n";
                element_type = "BRICK";
                break;
            case MshElemType::TRIANGLE:
                mat_file << m_ele->getNodeIndices()[0] + 1 << " "
                         << m_ele->getNodeIndices()[1] + 1 << " "
                         << m_ele->getNodeIndices()[2] + 1 << "\n";
                element_type = "TRIANGLE";
                break;
            case MshElemType::TETRAHEDRON:
                mat_file << m_ele->getNodeIndices()[0] + 1 << " "
                         << m_ele->getNodeIndices()[1] + 1 << " "
                         << m_ele->getNodeIndices()[2] + 1 << " "
                         << m_ele->getNodeIndices()[3] + 1 << "\n";
                element_type = "TETRAHEDRON";
                break;
            case MshElemType::PRISM:
                mat_file << m_ele->getNodeIndices()[0] + 1 << " "
                         << m_ele->getNodeIndices()[0] + 1 << " "
                         << m_ele->getNodeIndices()[1] + 1 << " "
                         << m_ele->getNodeIndices()[2] + 1 << " "
                         << m_ele->getNodeIndices()[3] + 1 << " "
                         << m_ele->getNodeIndices()[3] + 1 << " "
                         << m_ele->getNodeIndices()[4] + 1 << " "
                         << m_ele->getNodeIndices()[5] + 1 << "\n";
                element_type = "BRICK";
                break;
            default:
                std::cerr
                    << "CMediumProperties::WriteTecplotDistributedProperties "
                       "MshElemType not handled"
                    << "\n";
        }
    }
}

/**************************************************************************
   MSHLib-Method: GetNearestHetVal2
   Task:
   Programing:
   0?/2004 SB Implementation
   09/2005 MB EleClass
   01/2006 SB ReImplementation with new concept by Olaf, moved here
**************************************************************************/
long GetNearestHetVal2(long EleIndex,
                       CFEMesh* m_msh,
                       vector<double>
                           xvals,
                       vector<double>
                           yvals,
                       vector<double>
                           zvals,
                       vector<double>
                           mmpvals)
{
    (void)mmpvals;
    long i, nextele, no_values;
    double ex, ey, ez, dist, dist1;  // WW , dist2;
    double x, y, z;
    MeshLib::CElem* m_ele = NULL;
    no_values = (long)xvals.size();

    x = 0.0;
    y = 0.0;
    z = 0.0;
    dist = 10000000.0;  // Startwert
    // WW dist2 = 0.01;                                  // Abstand zwischen
    // eingelesenen Knoten und Geometrieknoten-RF; Achtung, doppelbelegung
    // möglich bei kleinen Gitterabständen
    nextele = -1;

    // Get element data
    m_ele = m_msh->ele_vector[EleIndex];
    double const* center(m_ele->GetGravityCenter());
    x = center[0];
    y = center[1];
    z = center[2];

    // Calculate distances
    for (i = 0; i < no_values; i++)
    {
        ex = xvals[i];
        ey = yvals[i];
        ez = zvals[i];
        dist1 = (ex - x) * (ex - x) + (ey - y) * (ey - y) + (ez - z) * (ez - z);
        if (dist1 < dist)
        {
            dist = dist1;
            nextele = i;
        }
    }

    return nextele;
}

/**************************************************************************
   MSHLib-Method: GetAverageHetVal2
   Task:
   Programing:
   06/2005 MB Implementation
   01/2006 SB Adapted to new structure
**************************************************************************/
double GetAverageHetVal2(long EleIndex,
                         CFEMesh* m_msh,
                         vector<double>
                             xvals,
                         vector<double>
                             yvals,
                         vector<double>
                             zvals,
                         vector<double>
                             mmpvals)
{
    long i, j, ihet;
    double average;
    double xp[3], yp[3];
    double value;
    double NumberOfValues;
    // WW double InvNumberOfValues;
    CGLPoint* m_point = NULL;
    MeshLib::CElem* m_ele = NULL;
    long no_values = (long)xvals.size();

    j = 0;  // only for 1 value

    //-----------------------------------------------------------------------
    // Get element data
    m_ele = m_msh->ele_vector[EleIndex];
    for (j = 0; j < 3; j++)
    {
        double const* const pnt(m_ele->GetNode(j)->getData());
        xp[j] = pnt[0];
        yp[j] = pnt[1];
        // zp[j] = 0.0;
    }

    //-----------------------------------------------------------------------
    // Find data points in the element
    NumberOfValues = 0;
    // WW InvNumberOfValues = 0;
    m_point = new CGLPoint;

    average = -1;
    value = 0;

    for (i = 0; i < no_values; i++)
        if (mmpvals[i] != -999999.0)  // Data point not within an element yet
        {
            m_point->x = xvals[i];
            m_point->y = yvals[i];
            m_point->z = 0.0;

            //....................................................................
            // Calculate the product of values in element
            // CC 10/05
            if (m_point->IsInTriangleXYProjection(xp, yp))
            {
                value = value + zvals[i];
                NumberOfValues++;
                mmpvals[i] = -999999.0;  // used as marker
            }
        }
    // end for
    //........................................................................
    if (NumberOfValues ==
        0)  // if no data points in element --> get neares value
    {
        ihet = GetNearestHetVal2(EleIndex, m_msh, xvals, yvals, zvals, mmpvals);
        if (ihet < 0)
            DisplayMsgLn(" Error getting nearest het_value location");
        else
            average = mmpvals[ihet];
    }
    //........................................................................
    else  // if data points in element --> Calculate arithmetic mean

        average = value / NumberOfValues;
    delete m_point;
    return average;
}

/**************************************************************************
   FEMLib-Method:
   Task: set MMP group member
   Programing:
   01/2006 YD Implementation
**************************************************************************/
void CMediumPropertiesGroup::Set(CRFProcess* m_pcs)
{
    long j, k;
    CFEMesh* m_msh = m_pcs->m_msh;
    CMediumProperties* m_mmp = NULL;
    MeshLib::CElem* elem = NULL;
    //----------------------------------------------------------------------
    // Tests //
    if (!m_msh)
        cout << "Warning in CSourceTermGroup::Set - no MSH data"
             << "\n";
    // return;
    //----------------------------------------------------------------------
    long no_mmp = (long)mmp_vector.size();
    for (j = 0; j < no_mmp; j++)
    {
        m_mmp = mmp_vector[j];
        //====================================================================
        if (m_mmp->pcs_type_name.compare(pcs_type_name) == 0)
        {
            m_mmp = mmp_vector[j];
            for (k = 0; k < (long)m_msh->ele_vector.size(); k++)
            {
                elem = m_msh->ele_vector[k];
                elem->SetPatchIndex(j);
            }
        }
    }
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   05/2007 WW Implementation
**************************************************************************/
void CMediumProperties::CalStressPermeabilityFactor(double* kfac,
                                                    const double T)
{
    switch (permeability_stress_mode)
    {
        case 2:
            CalStressPermeabilityFactor2(kfac, T);
            break;
        case 3:
            CalStressPermeabilityFactor3(kfac);
            break;
    }
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   05/2007 WW Implementation
**************************************************************************/
void CMediumProperties::CalStressPermeabilityFactor2(double* kfac,
                                                     const double T)
{
    int i, ia, ib, ele_index;
    double xyz[3], sig[3], b[3], b0;
    ElementValue_DM* e_valDM = NULL;
    ele_index = Fem_Ele_Std->Index;
    Fem_Ele_Std->RealCoordinates(xyz);
    b0 = pow(6.0 * permeability_tensor[0] / c_coefficient[5], 1.0 / 3.0);
    if ((long)ele_value_dm.size() > 0)  // Deformation process coupled
    {
        e_valDM = ele_value_dm[ele_index];
        for (i = 0; i < 3; i++)
        {
            sig[i] = (*e_valDM->Stress)(i, Fem_Ele_Std->GetGPindex());
            b[i] = c_coefficient[2 + i] +
                   c_coefficient[0] *
                       exp(c_coefficient[1] * sig[i] /
                           (T * PhysicalConstant::IdealGasConstant));
            //
        }
    }
    else
        for (i = 0; i < 3; i++)
        {
            sig[i] = c_coefficient[6 + i * 4] +
                     c_coefficient[7 + i * 4] * xyz[0] +
                     c_coefficient[8 + i * 4] * xyz[1] +
                     c_coefficient[9 + i * 4] * xyz[2];
            b[i] = c_coefficient[2 + i] +
                   c_coefficient[0] *
                       exp(c_coefficient[1] * sig[i] /
                           (T * PhysicalConstant::IdealGasConstant));
        }
    for (i = 0; i < 3; i++)
    {
        ia = (i + 1) % 3;
        ib = (i + 2) % 3;
        kfac[i] = 0.5 *
                  (MathLib::fastpow(b[ia], 3) + MathLib::fastpow(b[ib], 3)) /
                  MathLib::fastpow(b0, 3);
    }
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   05/2007 WW Implementation
**************************************************************************/
void CMediumProperties::CalStressPermeabilityFactor3(double* kfac)
{
    int i, ele_index;
    double am, pG, fy;
    double xyz[3], sig[3], ah[3];
    double JRC = c_coefficient[0];
    double a01 = c_coefficient[1];
    double Kn = c_coefficient[2];
    //
    ElementValue_DM* e_valDM = NULL;
    ele_index = Fem_Ele_Std->Index;
    pG = Fem_Ele_Std->PG;
    if (pG < 0.0)
        pG = 0.0;
    Fem_Ele_Std->RealCoordinates(xyz);
    if ((long)ele_value_dm.size() > 0)  // Deformation process coupled
    {
        e_valDM = ele_value_dm[ele_index];
        for (i = 0; i < 3; i++)
            sig[i] = 1.e-6 * ((*e_valDM->Stress)(i, Fem_Ele_Std->GetGPindex()) -
                              max(pG, 0.0));
    }
    else
        for (i = 0; i < 3; i++)
            sig[i] =
                1.0e-6 *
                (c_coefficient[9 + i * 4] + c_coefficient[10 + i * 4] * xyz[0] +
                 c_coefficient[11 + i * 4] * xyz[1] +
                 c_coefficient[12 + i * 4] * xyz[2] - max(pG, 0.0));
    // am at 100
    double am0_h =
        a01 - (c_coefficient[7] / (-Kn + c_coefficient[7] / c_coefficient[3]) +
               c_coefficient[4] + c_coefficient[5] + c_coefficient[6]);
    double am0_H =
        a01 - (c_coefficient[8] / (-Kn + c_coefficient[8] / c_coefficient[3]) +
               c_coefficient[4] + c_coefficient[5] + c_coefficient[6]);
    double ah0_h = am0_h * am0_h;
    double ah0_H = am0_H * am0_H;
    if (ah0_h > am0_h)
        ah0_h = am0_h;
    if (ah0_H > am0_H)
        ah0_H = am0_h;
    for (i = 0; i < 3; i++)
    {
        am = a01 - (sig[i] / (-Kn + sig[i] / c_coefficient[3]) +
                    c_coefficient[4] + c_coefficient[5] + c_coefficient[6]);
        ah[i] = am * am;
        if (ah[i] > am)
            ah[i] = am;
        //
        c_coefficient[21 + i] = ah[i] / MathLib::fastpow(sqrt(JRC), 5);
    }
    kfac[0] = ah[0] * ah[0] / (ah0_h * ah0_h);
    if (geo_dimension == 2)
    {
        fy = ah[2] * ah[2] / (ah0_H * ah0_H);
        kfac[1] = 0.5 * (fy + kfac[0]);
    }
    else if (geo_dimension == 3)
    {
        fy = ah[1] * ah[1] / (ah0_H * ah0_H);
        kfac[1] = fy;
        kfac[2] = 0.5 * (kfac[0] + kfac[1]);
    }
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   06/2007 WW Implementation
**************************************************************************/
void CMediumProperties::CalStressPermeabilityFactor3_Coef()
{
    int i;
    double d_max = 0.0;
    double a[4], delta[4], Kni[4];
    double A[] = {-0.296, -0.1005, -0.1031, -0.1031};
    double B[] = {-0.0056, -0.0073, -0.0074, -0.0074};
    double C[] = {2.241, 1.0082, 1.135, 1.135};
    double D[] = {-0.245, -0.23, -0.251, -0.251};
    //
    double C1[] = {84.77, 44.37, 31.38, 20};
    double C2[] = {0.02, 0.01, 0.01, 0.01};
    //
    double JRC = c_coefficient[0];
    double JCS = c_coefficient[1];
    double UCS = c_coefficient[2];
    double a0 = 0.2 * JRC * (0.2 * UCS / JCS - 0.1);
    //
    for (i = 0; i < 4; i++)
    {
        a[i] = a0;
        d_max = A[i] + B[i] * JRC + C[i] * pow(JCS / a0, D[i]);
        Kni[i] = 0.0178 * JCS / a[i] + 1.748 * JRC - 7.155;
        delta[i] = 0.01 * (C1[i] - C2[i] * JCS / a0) * d_max;
        a0 -= delta[i];
    }
    //
    //
    c_coefficient[1] = a[0];
    c_coefficient[2] = Kni[3];
    c_coefficient[3] = d_max;
    for (i = 4; i < 7; i++)
        c_coefficient[i] = delta[i - 4];
    // Unit of stresses is MPa
    c_coefficient[7] *= 1.0e-6;
    c_coefficient[8] *= 1.0e-6;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   01/2006 YD Implementation
   last modification:
**************************************************************************/
CMediumPropertiesGroup* MMPGetGroup(const string& pcs_type_name)
{
    CMediumPropertiesGroup* m_mmp_group = NULL;
    list<CMediumPropertiesGroup*>::const_iterator p_mmp_group =
        mmp_group_list.begin();
    while (p_mmp_group != mmp_group_list.end())
    {
        m_mmp_group = *p_mmp_group;
        if (m_mmp_group->pcs_type_name.compare(pcs_type_name) == 0)
            return m_mmp_group;
        ++p_mmp_group;
    }
    return NULL;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   01/2006 YD Implementation
   last modified:
**************************************************************************/
void MMPGroupDelete(/*string pcs_type_name*/)
{
    CMediumPropertiesGroup* m_mmp_group = NULL;
    list<CMediumPropertiesGroup*>::const_iterator p = mmp_group_list.begin();
    while (p != mmp_group_list.end())
    {
        m_mmp_group = *p;
        // if(m_mmp_group->pcs_type_name.compare(pcs_type_name)==0){
        delete m_mmp_group;
        // mmp_group_list.remove(m_mmp_group);
        // return;
        ++p;
    }
    mmp_group_list.clear();
}

/**************************************************************************
   FEMLib-Method:
   08/2007 OK Implementation
**************************************************************************/
bool MMPExist(ifstream* mmp_file)
{
    string line_string;
    std::stringstream in;
    string name;
    //======================================================================
    while (!mmp_file->eof())
    {
        line_string = GetLineFromFile1(mmp_file);
        //--------------------------------------------------------------------
        // NAME
        if (line_string.find("$NAME") != string::npos)
        {
            in.str(GetLineFromFile1(mmp_file));
            in >> name;  // sub_line
            // Test whether MMP already exists
            for (int i = 0; i < (int)mmp_vector.size(); i++)
                if (name.compare(mmp_vector[i]->name) == 0)
                    return true;

            in.clear();
            continue;
        }
        //--------------------------------------------------------------------
    }
    //======================================================================
    return false;
}

/**************************************************************************
   FEMLib-Method:
   08/2007 OK Implementation
   06/2009 OK Bug fix
**************************************************************************/
bool MMPExist()
{
    for (int i = 0; i < (int)mmp_vector.size(); i++)
        for (int j = 0; j < (int)mmp_vector.size(); j++)
        {
            if (i == j)
                continue;
            if (mmp_vector[i]->name.compare(mmp_vector[j]->name) == 0)
                return true;
        }
    return false;
}

/**************************************************************************
   FEMLib-Method:
   Task: retrun the new permeability based on original permeability
      and old/new porosity
   Programing:
   11/2008 HS Implementation
   last modification:
**************************************************************************/
double CMediumProperties::KozenyCarman(double k0, double n0, double n)
{
    double rt = 0.0;

    // TODO: here it should be the min_porosity and max_porosity instead of 0
    // and 1
    if (k0 < 1.0e-40 || k0 > 1.0 || n0 <= 0 || n0 >= 1 || n <= 0 || n >= 1)
        return k0;
    else
    {
        rt = k0;

        rt *= MathLib::fastpow(n / n0, 3);
    }

    return rt;
}

/**************************************************************************
   FEMLib-Method:
   Task: retrun the new permeability based on original permeability
      and old/new porosity (Koseny-Carman normalized)
   Programing:
   11/2008 HS Implementation
   last modification:
**************************************************************************/
double CMediumProperties::KozenyCarman_normalized(double k0, double n0,
                                                  double n)
{
    double rt = 0.0;

    // TODO: here it should be the min_porosity and max_porosity instead of 0
    // and 1
    if (k0 < 1.0e-40 || k0 > 1.0 || n0 <= 0 || n0 >= 1 || n <= 0 || n >= 1)
        return k0;
    else
    {
        rt = k0;

        rt *= MathLib::fastpow(n / n0, 3);

        rt *= ((1 - n0) / (1 - n) * (1 - n0) / (1 - n));
    }

    return rt;
}

/**************************************************************************
FEMLib-Method:
Task: Returns the new permeability which will be calculated using
      Kozeny-Carman formulation. It takes intial permeability& porosity
      and also the new/updated porosity as an input.
Programing:
11.2010 AB Implementation
last modification:
**************************************************************************/
double CMediumProperties::KozenyCarmanNew(double k_init, double n_init,
                                          double n_t)
{
    double k_t = 0.0;

    // Limit minimum and maximum values.
    if (k_init <= 0.0 || n_init <= 0. || n_init >= 1. || n_t <= 0. || n_t >= 1.)
        return k_init;
    else
    {
        k_t = k_init * (pow((1. - n_init) / (1 - n_t), 2)) *
              (pow(n_t / n_init, 3));
    }

    return k_t;
}

/*************************************************************************************
FEMLib-Method:
Task: Returns the new permeability which will be calculated using
      Verma-Pruess formulation. It takes intial permeability& porosity,
      the new/updated porosity, critical porosity and a power law exponent as an
input. Programing: 11.2010 AB Implementation last modification:
*************************************************************************************/
double CMediumProperties::VermaPruess(double k_init, double n_init, double n_t)
{
    double k_t = 0.0;
    double n_crit, a;

    n_crit = permeability_porosity_updating_values[0];  // critical porosity
    a = permeability_porosity_updating_values[1];       // a power law exponent

    // Limit minimum and maximum values.
    if (k_init <= 0.0 || k_init >= 1.0 || n_init <= 0 || n_init >= 1 ||
        n_t <= 0 || n_t >= 1)
        return k_init;
    else
    {
        k_t = k_init * (pow((n_t - n_crit) / (n_init - n_crit), a));
    }

    return k_t;
}
///////////////////////////////////////////////////////////////////////////
// old data structure functions //OK411

/**************************************************************************
   //3.1 Subfunction of Porosity
   ROCKFLOW - Funktion: MATCalcSoilPorosityMethod1

   Aufgabe:
   Takes the porosity from a Geomechanical model, calulcates the effective
stress and then takes the value of porosity from a curve.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E SOIL_PROPERTIES *sp: Zeiger auf eine Instanz vom Typ SOIL_PROPERTIES.

   Ergebnis:
   - s.o. -
   Programming Example
   Programmaenderungen:
   03/2004  CMCD First version

**************************************************************************/
double CMediumProperties::PorosityEffectiveStress(long index,
                                                  double element_pressure)
{
    // OK411
    element_pressure = element_pressure;
    index = index;

    /*OK411
        int  nn, i, idummy, Type;
        long *element_nodes;
        double p;
       double znodes[8],ynodes[8],xnodes[8];
        double zelnodes[8],yelnodes[8],xelnodes[8];
       double coords[3];
       double Pie,angle;
       double sigma1,sigma2,sigma3;
       double a1,a2,a3,b1,b2,b3;
       double normx,normy,normz,normlen;
       double dircosl, dircosm, dircosn;
       double tot_norm_stress, eff_norm_stress;
       int material_group;
       double x_mid, y_mid, z_mid;

       // Normal stress calculated according to the orientation of the fracture
       element

       Type=ElGetElementType(index);
       material_group = ElGetElementGroupNumber(index);

       dircosl = dircosm = dircosn = 0.0;//Initialise variable

       if (Type == 2||Type == 3||Type == 4)  //Function defined for square,
       triangular and cubic elements
       {
       nn = ElNumberOfNodes[Type - 1];
       element_nodes = ElGetElementNodes(index);

       // Calculate directional cosins, note that this is currently set up
       // Sigma1 is in the y direction
       // Sigma2 is in the z direction
       // Sigma3 is in the x direction
       // This correspondes approximately to the KTB site conditions

       for (i=0;i<nn;i++)
       {
       zelnodes[i]=GetNodeZ(element_nodes[i]);
       yelnodes[i]=GetNodeY(element_nodes[i]);
       xelnodes[i]=GetNodeX(element_nodes[i]);
       }

       // Coordinate transformation to match sigma max direction
       // y direction matches the north direction
       Pie=3.141592654;
       angle=(porosity_model_values[3]*Pie)/180.;
       for (i=0;i<nn;i++)
       {
       znodes[i]=zelnodes[i];
       xnodes[i]=xelnodes[i]*cos(angle)-yelnodes[i]*sin(angle);
       ynodes[i]=xelnodes[i]*sin(angle)+yelnodes[i]*cos(angle);

       }

       if (Type == 2) //Square
       {
       a1=xnodes[2]-xnodes[0];
       a2=ynodes[2]-ynodes[0];
       a3=znodes[2]-znodes[0];
       b1=xnodes[3]-xnodes[1];
       b2=ynodes[3]-ynodes[1];
       b3=znodes[3]-znodes[1];
       normx=a2*b3-a3*b2;
       normy=a3*b1-a1*b3;
       normz=a1*b2-a2*b1;
       normlen=sqrt(pow(normx,2.0)+pow(normy,2.0)+pow(normz,2.0));
       dircosl= normy/normlen;
       dircosm= normz/normlen;
       dircosn= normx/normlen;
       }

       if (Type == 3) //Cube
       {
       a1=xnodes[2]-xnodes[0];
       a2=ynodes[2]-ynodes[0];
       a3=znodes[2]-znodes[0];
       b1=xnodes[3]-xnodes[1];
       b2=ynodes[3]-ynodes[1];
       b3=znodes[3]-znodes[1];
       normx=a2*b3-a3*b2;
       normy=a3*b1-a1*b3;
       normz=a1*b2-a2*b1;
       normlen=sqrt(pow(normx,2.0)+pow(normy,2.0)+pow(normz,2.0));
       dircosl= normy/normlen;
       dircosm= normz/normlen;
       dircosn= normx/normlen;
       }

       if (Type == 4) //Triangle
       {
       a1=xnodes[1]-xnodes[0];
       a2=ynodes[1]-ynodes[0];
       a3=znodes[1]-znodes[0];
       b1=xnodes[2]-xnodes[1];
       b2=ynodes[2]-ynodes[1];
       b3=znodes[2]-znodes[1];
       normx=a2*b3-a3*b2;
       normy=a3*b1-a1*b3;
       normz=a1*b2-a2*b1;
       normlen=sqrt(pow(normx,2.0)+pow(normy,2.0)+pow(normz,2.0));
       dircosl= normy/normlen;
       dircosm= normz/normlen;
       dircosn= normx/normlen;
       }

       // Calculation of average location of element
       CalculateSimpleMiddelPointElement(index,coords);
       x_mid = coords[0];
       y_mid = coords[1];
       z_mid = coords[2];

       //Calculate fluid pressure in element
       p = element_pressure;

       //Calcualtion of stress according to Ito & Zoback 2000 for KTB hole

       sigma1=z_mid*0.045*-1e6;
       sigma2=z_mid*0.028*-1e6;
       sigma3=z_mid*0.02*-1e6;

       //Calculate total normal stress on element
       //Note in this case sigma2 corresponds to the vertical stress
       tot_norm_stress=sigma1*dircosl*dircosl+sigma2*dircosm*dircosm+sigma3*dircosn*dircosn;

       //Calculate normal effective stress
       eff_norm_stress=tot_norm_stress-p;

       //Take value of storage from a curve
       porosity=GetCurveValue((int) porosity_model_values[0], 0,
       eff_norm_stress, &idummy);

       }
       // default value if element type is not included in the method for
       calculating the normal else porosity=porosity_model_values[0];
     */
    return porosity;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   05/2005 MX Implementation
   last modification:
**************************************************************************/
double CMediumProperties::PorosityVolumetricFreeSwellingConstantIonicstrength(
    long index, double saturation, double temperature)
{
    // OK411
    index = index;
    saturation = saturation;
    temperature = temperature;

    /*OK411
       double mat_mp_m, beta;
       double porosity = 0.0;
       static double theta;
       static double porosity_n, porosity_IL, d_porosity, \
                    density_rock, fmon;
       static double S_0;
       static double epsilon;
       static double satu_0=0.20;
       //  static double porosity_min=0.05;
       static double ion_strength;
       static double F_const=96484.6, epsilon_0=8.854e-12;
       static double R=8.314510, psi=1.0;
       theta = 1.0;
       //--------------------------------------------------------------------
       // MMP medium properties
       S_0 =      porosity_model_values[1];      // Specific surface area
     [m^2/g] fmon =     porosity_model_values[2];     // Anteil quelfaehige
     mineral [-] mat_mp_m = porosity_model_values[3]; // Schichtenanzahl eines
     quellfähigen Partikels [-] beta =     porosity_model_values[6];     //
     modifications coefficient (z.B. free swelling Weimar beta=3.0)
       //--------------------------------------------------------------------
       // MSP solid properties
       CSolidProperties *m_msp = NULL;
       long group = ElGetElementGroupNumber(index);
       m_msp = msp_vector[group];
       density_rock  = m_msp->Density(1);
       //--------------------------------------------------------------------
       // State properties
       ion_strength = porosity_model_values[4]; // Ionic strength [M]
       satu_0 =       porosity_model_values[5];       // Initial saturation, if
     the sample is homogenous [-]
       //--------------------------------------------------------------------
       // Interlayer porosity calculation
       if (abs(temperature)>1.0e10) temperature=298.0;  //TODO MX
       epsilon = 87.0 +
     exp(-0.00456*(temperature-PhysicalConstant::CelsiusZeroInKelvin));
       porosity_n = porosity_model_values[0];
       // Maximal inter layer porosity
       porosity_IL = fmon * psi * mat_mp_m * S_0 * (density_rock * 1.0e3) \
     * sqrt(epsilon * epsilon_0 * R * temperature / (2000.0 * F_const * F_const
     * ion_strength ));
       d_porosity=porosity_IL*(pow(saturation,beta)-pow(satu_0,beta));
       porosity_IL *=pow(saturation,beta);
       ElSetElementVal(index,PCSGetELEValueIndex("POROSITY_IL"),porosity_IL);
       // Total porosity calculation
       porosity = porosity_n+d_porosity;
       ElSetElementVal(index,PCSGetELEValueIndex("POROSITY"),porosity);
     */
    return porosity;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   05/2005 MX Implementation
   last modification:
**************************************************************************/
double CMediumProperties::PorosityEffectiveConstrainedSwelling(
    long index, double saturation, double temperature, double* porosity_sw)
{
    // OK411
    index = index;
    saturation = saturation;
    temperature = temperature;
    porosity_sw = porosity_sw;

    /*OK411
       // Soil Properties

       double mat_mp_m, beta;
       double porosity = 0.0;
       static double theta;
       double porosity_n, porosity_IL, d_porosity, \
                    n_total, density_rock, fmon;
       //  double n_max, n_min;
       static double S_0;
       // Component Properties
       static double  satu, epsilon;
       static double satu_0=0.20;
       static double porosity_min=0.05;
       static double ion_strength;
       static double F_const=96484.6, epsilon_0=8.854e-12;
       static double R=8.314510, psi=1.0;
       // Fluid properies
       int phase=1;

       //--------------------------------------------------------------------
       // MMP medium properties
       S_0 =      porosity_model_values[1];      // Specific surface area
     [m^2/g] fmon =     porosity_model_values[2];     // Anteil quelfaehige
     mineral [-] mat_mp_m = porosity_model_values[3]; // Schichtenanzahl eines
     quellfähigen Partikels [-] beta =     porosity_model_values[6];     //
     modifications coefficient (z.B. free swelling Weimar beta=3.0)
       //--------------------------------------------------------------------
       // MSP solid properties
       CSolidProperties *m_msp = NULL;
       //long group = ElGetElementGroupNumber(index);

       long group = m_pcs->m_msh->ele_vector[index]->GetPatchIndex();
       m_msp = msp_vector[group];
       density_rock  = m_msp->Density(1);
       //--------------------------------------------------------------------

       // Component Properties
       ion_strength = MATCalcIonicStrengthNew(index);
       if (ion_strength == 0.0){
       ion_strength = porosity_model_values[4]; //Ionic strength [M]
       }
       satu_0 = porosity_model_values[5];       //Initial saturation, if the
     sample is homogenous [-] porosity_min = porosity_model_values[7]; //minimal
     porosity after swelling compaction

       //  ion_strength = MATGetSoilIonicStrength(index);

       // Field variables
       //theta = GetNumericalTimeCollocation("TRANSPORT");
       theta = 1.0;
       // T = MATGetTemperatureGP (index,0.0,0.0,0.0,theta);
       //  T=298.0;
       phase=1;
       satu = saturation; //MATGetSaturationGP(phase, index, 0.0, 0.0, 0.0,
     theta); //only for fluid, phase=1

       //-----------------------------------------------------------------------
       // Interlayer Porositaet berechnen
       if (abs(temperature)>1.0e10) temperature=298.0;  //TODO MX
       epsilon
     =87.0+exp(-0.00456*(temperature-PhysicalConstant::CelsiusZeroInKelvin));
       porosity_n = porosity_model_values[0];

       // Maximal inter layer porosity
       porosity_IL=fmon * psi * mat_mp_m * S_0 * (density_rock * 1.0e3) \
     * sqrt(epsilon * epsilon_0 * R * temperature / (2000.0 * F_const * F_const
     * ion_strength )); d_porosity=porosity_IL*(pow(satu, beta)-pow(satu_0,
     beta));

       //-----------Interlayer porosity calculation------------------

       //  porosity_IL = porosity_IL*satu;
       porosity_IL  *=pow(satu,beta);
       ElSetElementVal(index,PCSGetELEValueIndex("POROSITY_IL"),porosity_IL);
       // constrained swelling

       //-----------Effective porosity calculation------------------
       porosity = porosity_n-d_porosity;

       if (porosity<porosity_min)
       porosity =porosity_min;

       ElSetElementVal(index,PCSGetELEValueIndex("POROSITY"),porosity);
       //-----------Void ratio for swelling pressure
     calculation------------------
       //  e = porosity/(1.-porosity);
       //  ElSetElementVal(index,PCSGetELEValueIndex("VoidRatio"),e);
       //-----------Swelling potential calculation------------------
       // constrained swelling
       //  n_total=porosity_n - d_porosity;
       n_total=porosity_n + d_porosity-porosity_min;

       //  if(n_total>porosity_min)
       if(n_total>=0)
       {
     ***porosity_sw = n_total;
       }
       else
       //    *porosity_sw=-porosity_IL*(satu-satu_0)+(porosity_n-porosity_min);
     ***porosity_sw=d_porosity+(porosity_n-porosity_min);
       ElSetElementVal(index,PCSGetELEValueIndex("POROSITY_SW"),*porosity_sw);

       //-----------Swelling pressure calculation------------------
       // MATCalSoilSwellPressMethod0(index);
     */
    return porosity;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   05/2005 MX Implementation
   last modification:
**************************************************************************/
double CMediumProperties::PorosityVolumetricFreeSwelling(long index,
                                                         double saturation,
                                                         double temperature)
{
    // OK411
    index = index;
    saturation = saturation;
    temperature = temperature;

    /*OK411
       double mat_mp_m, beta;
       double porosity = 0.0;
       static double theta;
       static double porosity_n, porosity_IL, d_porosity, \
                    density_rock, fmon;
       static double S_0;
       static double epsilon;
       static double satu_0=0.20;
       //  static double porosity_min=0.05;
       static double ion_strength;
       static double F_const=96484.6, epsilon_0=8.854e-12;
       static double R=8.314510, psi=1.0;
       theta = 1.0;
       //--------------------------------------------------------------------
       // MMP medium properties
       S_0 =      porosity_model_values[1];      // Specific surface area
     [m^2/g] fmon =     porosity_model_values[2];     // Anteil quelfaehige
     mineral [-] mat_mp_m = porosity_model_values[3]; // Schichtenanzahl eines
     quellfähigen Partikels [-] beta =     porosity_model_values[6];     //
     modifications coefficient (z.B. free swelling Weimar beta=3.0)
       //--------------------------------------------------------------------
       // MSP solid properties
       CSolidProperties *m_msp = NULL;
       long group = ElGetElementGroupNumber(index);
       m_msp = msp_vector[group];
       density_rock  = m_msp->Density(1);
       //--------------------------------------------------------------------
       // State properties
       ion_strength = MATCalcIonicStrengthNew(index);
       if (ion_strength == 0.0){
       ion_strength = porosity_model_values[4]; // Ionic strength [M]
       }
       satu_0 =       porosity_model_values[5];       // Initial saturation, if
     the sample is homogenous [-]
       //--------------------------------------------------------------------
       // Interlayer porosity calculation
       epsilon = 87.0 +
     exp(-0.00456*(temperature-PhysicalConstant::CelsiusZeroInKelvin));
       porosity_n = porosity_model_values[0];
       // Maximal inter layer porosity
       porosity_IL = fmon * psi * mat_mp_m * S_0 * (density_rock * 1.0e3) \
     * sqrt(epsilon * epsilon_0 * R * temperature / (2000.0 * F_const * F_const
     * ion_strength ));
       d_porosity=porosity_IL*(pow(saturation,beta)-pow(satu_0,beta));
       porosity_IL *=pow(saturation,beta);
       ElSetElementVal(index,PCSGetELEValueIndex("POROSITY_IL"),porosity_IL);
       // Total porosity calculation
       porosity = porosity_n+d_porosity;
       ElSetElementVal(index,PCSGetELEValueIndex("POROSITY"),porosity);
     */
    return porosity;
}

/**************************************************************************
   FEMLib-Method: PorosityVolumetricChemicalReaction
   Task: Porosity variation owing to chemical reactions
   Programing:
   05/2005 MX Implementation
   last modification:
**************************************************************************/
double CMediumProperties::PorosityVolumetricChemicalReaction(long index)
{
    // OK411
    index = index;

    /*OK411
       int n=0, i, timelevel, nn=0, m=0;
       static long *element_nodes;
       //  long component;
       double porosity=0.0, tot_mineral_volume=0.0, tot_mineral_volume0=0.0;
       double conc[100];
       double mineral_volume[100],molar_volume[100];
       //  double *ergebnis=NULL;
       // REACTION_MODEL *rcml=NULL;

       //  rcml=GETReactionModel(index);
       REACT *rc = NULL; //SB
       rc->GetREACT();

       // MMP Medium Properties
       porosity = porosity_model_values[0];
       if (!rc) return porosity;

       //  m=rcml->number_of_equi_phases;
       m = rc->rcml_number_of_equi_phases;
       if (m ==0 || ElGetElement(index)==NULL) // wenn solid phases and Element
       existiert return porosity;

       // mineral molar volume abholen
       for (i=0; i<m; i++){
       molar_volume[i]=porosity_model_values[i+1];
       }

       tot_mineral_volume0=porosity_model_values[m+1]; //initial total volume of
       the minerals

       //  if (!rcml) {return porosity;}

       // calculate the concentrations of each solid phases (in mol/l soil) as
       element value from adjecent node values
       //  n=rcml->number_of_master_species;
       n = rc->rcml_number_of_master_species;
       //  n = get_rcml_number_of_master_species(rcml);
       for (int component=0; component<m+2+n+rc->rcml_number_of_ion_exchanges;
       component++) { conc[component] =0.0;
       // Not used: CompProperties *m_cp = cp_vec[component];
       //    int z_i = m_cp->valence;
       //	 m_cp->compname; //What's this for
       if (component>=n+2 && component<n+2+m){

       if (ElGetElementActiveState(index)){
       nn = ElGetElementNodesNumber(index);
       element_nodes = ElGetElementNodes(index);
       for (int j=0;j<nn;j++) {
       timelevel=1;
       conc[component] +=
       PCSGetNODConcentration(element_nodes[j],component,timelevel);
       }
       conc[component] /= (double)nn;
       element_nodes = NULL;
       }

       // calculate the solid phase: volume =v_mi*Ci
       timelevel=0;
       //      conc[i]=CalcElementMeanConcentration (index, i, timelevel,
       ergebnis); mineral_volume[component-n-2] =
       conc[component]*molar_volume[component-n-2]; tot_mineral_volume +=
       mineral_volume[component-n-2];
       }
       } //for

       porosity += tot_mineral_volume0 - tot_mineral_volume;
       //  ElSetElementVal(index,"POROSITY",porosity);
       ElSetElementVal(index,PCSGetELEValueIndex("POROSITY"),porosity);
     */
    return porosity;
}

// WX: 03.2011. Porosity = n0 + vol_strain
double CMediumProperties::PorosityVolStrain(long index, double val0,
                                            CFiniteElementStd* assem)
{
    double val = val0, vol_strain_temp = 0., strain_temp[3] = {0.},
           strain_nodes[20] = {0.};
    // WW int idx_temp[3]={0}, ele_index, nnodes = assem->nnodes;
    int idx_temp[3] = {0}, nnodes = assem->nnodes;
    // WW CRFProcessDeformation *dm_pcs = (CRFProcessDeformation *) this;
    int dim = m_pcs->m_msh->GetCoordinateFlag() / 10;
    if (dim == 2)
        if (assem->axisymmetry)
            dim = 3;
    idx_temp[0] = assem->dm_pcs->GetNodeValueIndex("STRAIN_XX");
    idx_temp[1] = assem->dm_pcs->GetNodeValueIndex("STRAIN_YY");
    idx_temp[2] = assem->dm_pcs->GetNodeValueIndex("STRAIN_ZZ");
    // WW ele_index = index;
    for (int j = 0; j < dim; j++)
    {
        for (int i = 0; i < nnodes; i++)
            strain_nodes[i] = assem->dm_pcs->GetNodeValue(
                assem->dm_pcs->m_msh->ele_vector[index]->getNodeIndices()[i],
                idx_temp[j]);
        strain_temp[j] = assem->interpolate(strain_nodes);
    }
    for (int j = 0; j < dim; j++)
        vol_strain_temp += strain_temp[j];
    val += vol_strain_temp;
    if (val < MKleinsteZahl)
        val = 1e-6;  // lower limit of porostity
    return val;
}

/**************************************************************************
   FEMLib-Method: TortuosityFunction
   Task:
   Programing:
   05/2007 PCH Diffusion tensor is handled by tortuosity tensor as in
      permeability. Agreed with GK
   last modification:
**************************************************************************/
double CMediumProperties::TortuosityFunction(long number, double* gp,
                                             double theta,
                                             CFiniteElementStd* assem)
{
    // OK411
    theta = theta;
    gp = gp;
    number = number;
    assem = assem;

    switch (tortuosity_model)
    {
        case 0:  // Tortuosity is read from a curve
            // To do
            break;
        case 1:  // Constant Tortuosity
            tortuosity = tortuosity_model_values[0];
            break;
        case 2: /* Tortuosity is lin. fct. of porosity*/
            tortuosity = Porosity(number, theta) * tortuosity_model_values[0];
            break;
        default:
            DisplayMsgLn("Unknown tortuosisty model!");
            break;
    }
    return tortuosity;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: NonlinearFlowFunction
 */
/* Aufgabe:
   Berechnung der relativen Permeabilitaet
   fuer nichtlineare Fliessgesetze
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
 */
/* Ergebnis:
   0 bei Fehler, sonst 1 (dient nur zum Abbrechen der Funktion)
 */
/* Programmaenderungen:
   06/1999   OK   Implementierung
   09/2004   CMCD In GeoSys 4
 */
/**************************************************************************/
double CMediumProperties::NonlinearFlowFunction(long index, int gp,
                                                double /*theta*/,
                                                CFiniteElementStd* assem)
{
    double k_rel = 1.0;
    // OK411
    /*
    if (flowlinearity_model==3 || flowlinearity_model==6) // Forchheimer for DLR
    {
       assert (assem->PcsType==EnumProcessType::A ||
    assem->PcsType==EnumProcessType::S); ElementValue* gp_ele =
    ele_gp_value[index];
       //Velocity
       double vel[3] = {};
       vel[0] = gp_ele->Velocity(0, gp);
       vel[1] = gp_ele->Velocity(1, gp);
       vel[2] = gp_ele->Velocity(2, gp);
       double v_mag = sqrt(vel[0]*vel[0]+vel[1]*vel[1]+vel[2]*vel[2]);
       if (v_mag != 0.0) {
           double* k = PermeabilityTensor(index); //temporary assuming isotropic
           double dens_arg[3] = {};
           dens_arg[0] = (1.0-theta)*assem->interpolate(assem->NodalVal0) +
    theta*assem->interpolate(assem->NodalVal1); dens_arg[1] =
    (1.0-theta)*assem->interpolate(assem->NodalVal_t0) +
    theta*assem->interpolate(assem->NodalVal_t1);
           //dens_arg[1] = assem->interpolate(assem->NodalValC1)+T_KILVIN_ZERO;
           dens_arg[2] = (1.0-theta)*assem->interpolate(assem->NodalVal_X0) +
    theta*assem->interpolate(assem->NodalVal_X1);
           //dens_arg[2] = index;
           double vis = assem->FluidProp->Viscosity(dens_arg);
           double rhog = assem->FluidProp->Density(dens_arg); // which model?
           if (flowlinearity_model==6) {
               if (ParticleDiameter()==.0) {
                   std::cout << "***Error: dp = .0" << std::endl;
                   exit(0);
               }
               forchheimer_cf = 0.55*(1.-5.5*ParticleDiameter()/forchheimer_De);
           }
           k_rel = 1.0/(1.0+forchheimer_cf*rhog*sqrt(k[0])/vis*v_mag);
       }
    }
    */
    if (flowlinearity_model ==
        4)  // general Forchheimer, -dp/dx = a1 q + a2 q^2
    {
        assert(gp != -1);
        ElementValue* gp_ele = ele_gp_value[index];
        // Velocity
        double vel[3] = {};
        vel[0] = gp_ele->Velocity(0, gp);
        vel[1] = gp_ele->Velocity(1, gp);
        vel[2] = gp_ele->Velocity(2, gp);
        // double v_mag = sqrt(vel[0]*vel[0]+vel[1]*vel[1]+vel[2]*vel[2]);
        double v_mag = MBtrgVec(vel, 3);
#ifdef STATIC_K1
        k_rel =
            1.0 / (1.0 + forchheimer_a2 / forchheimer_a1 * v_mag);  // k0=1/a1
#else
        k_rel = 1.0 / (forchheimer_a1 + forchheimer_a2 * v_mag);  // k0=1.0
#endif
        // if (k_rel!=1.0)
        //    std::cout << index << "-" << gp << ": k_rel=" << k_rel << endl;
    }
    else if (flowlinearity_model ==
             5)  // general Forchheimer, -dp/dx = a1 q + a2 q^2
    {
        assert(gp != -1);
        ElementValue* gp_ele = ele_gp_value[index];
        // Velocity
        double vel[3] = {};
        vel[0] = gp_ele->Velocity(0, gp);
        vel[1] = gp_ele->Velocity(1, gp);
        vel[2] = gp_ele->Velocity(2, gp);
        // double v_mag = sqrt(vel[0]*vel[0]+vel[1]*vel[1]+vel[2]*vel[2]);
        double v_mag = MBtrgVec(vel, 3);
        if (v_mag != 0.0)
        {
            double* k_tensor = assem->MediaProp->PermeabilityTensor(index);
            double k0 = k_tensor[0];
            k_rel = 1.0 / (1.0 + forchheimer_a2 * k0 * v_mag);
        }
    }
    else
    {
        std::cout << "***ERROR: not supported flow linearity model "
                  << flowlinearity_model << std::endl;
    }

    /*OK411
       //Pressure variable name (PRESSURE 1) not used as yet in this function
       int i;
       //Pointer to fluid properties
       //	int no_phases =(int)mfp_vector.size();
       CFluidProperties *m_mfp = NULL;
       int phase = (int)mfp_vector.size()-1;
       m_mfp = mfp_vector[phase];
        // Knotendaten
        int nn,Type;
        long *element_nodes;
       double p_element_node[4], h_element_node[4], z_element_node[4];
       // Materialdaten
       double alpha;
       double g, rho, mu;
       // Elementdaten
       double grad_h[2],grad_x,grad_y;
       double mult[2];
       double detjac, *invjac, invjac2d[4];
       double grad_omega[8];
       double grad_h_min = MKleinsteZahl;
       double xgt[3],ygt[3],zgt[3];				//CMCD Global x,y,z coordinates
       of traingular element
       double xt[3],yt[3];							//CMCD Local x,y coordinates
       of traingular element double pt_element_node[4], ht_element_node[4],
       zt_element_node[4]; //CMCD Pressure, depth head traingular element double
       dN_dx[3],dN_dy[3], area;				//CMCD Shape function derivates for
       triangles double isotropicgradient,porosity, Re, Lambda, apperture,
       hyd_radius, perm; double linear_q,turbulent_q,linear_grad,
       turbulent_grad, flow_rate, temp; double dircos[6];
       //CMCD 04 2004
       gp[0]= gp[1] = gp[2] = 0.0;					//Gaus points, triangular
       interpretation value not relevant

       element_nodes = ElGetElementNodes(index);
       g = gravity_constant;
       //OK4104 theta = GetNumericalTimeCollocation("PRESSURE1");
       Type = ElGetElementType(index);
       //--------------------------------------------------------------------
       // MMP medium properties
       CMediumProperties *m_mmp = NULL;
       long group = ElGetElementGroupNumber(index);
       m_mmp = mmp_vector[group];
       //--------------------------------------------------------------------
       switch (flowlinearity_model){
       case 1://Element Type Dependent
       alpha = flowlinearity_model_values[0];
       rho = m_mfp->Density();

       switch (Type){
       case 1:
       nn = ElNumberOfNodes[0];
       for (i = 0; i < nn; i++) {
       p_element_node[i] = GetNodeVal(element_nodes[i],1);
       z_element_node[i] = GetNode(element_nodes[i])->z;
       h_element_node[i] = (p_element_node[i]) / (g * rho) + z_element_node[i];
       }
       invjac = GEOGetELEJacobianMatrix(index, &detjac);
       //          if(fabs(h_element_node[1]-h_element_node[0])>MKleinsteZahl)
       if (fabs(h_element_node[1] - h_element_node[0]) > grad_h_min) {
       k_rel = \
       pow(fabs(0.5 * sqrt(MSkalarprodukt(invjac, invjac, 3))), (alpha - 1.0)) *
       \ pow(fabs(h_element_node[1] - h_element_node[0]), (alpha - 1.0)); if
       (k_rel > 1) k_rel = 1.; } else k_rel = 1.0;

       break;

       case 2:
       nn = ElNumberOfNodes[1];        // Knotenanzahl nn muss 4 sein !
       for (i = 0; i < nn; i++) {
       p_element_node[i] = GetNodeVal(element_nodes[i],1);
       z_element_node[i] = GetNode(element_nodes[i])->z;
       h_element_node[i] = p_element_node[i] / (g * rho) + z_element_node[i];
       }
       Calc2DElementJacobiMatrix(index, 0.0, 0.0, invjac2d, &detjac);
       MGradOmega2D(grad_omega, 0, 0);         // Gradientenmatrix
       MMultMatVec(grad_omega, 2, 4, h_element_node, 4, mult, 2);
       MMultVecMat(mult, 2, invjac2d, 2, 2, grad_h, 2);
       //          if(
       (fabs(grad_h[0])>MKleinsteZahl)||(fabs(grad_h[1])>MKleinsteZahl) ) ) if
       ((fabs(grad_h[0]) > grad_h_min) || (fabs(grad_h[1]) > grad_h_min)) {
       k_rel = \
       pow(fabs(sqrt((grad_h[0]) * (grad_h[0]) + (grad_h[1]) * (grad_h[1]))), \
       (alpha - 1.0));
       } else {
       k_rel = 1.0;
       }
       break;
       case 3:
       k_rel = 1.;
       break;
       }

       case 2://Equivalent Fractured Media represented by triangles   CMCD April
       2004

       //Geometry
       nn = ElNumberOfNodes[Type - 1];
       element_nodes = ElGetElementNodes(index);

       for (i = 0; i < nn; i++)
       {
       xgt[i] = GetNodeX(element_nodes[i]);
       ygt[i] = GetNodeY(element_nodes[i]);
       zgt[i] = GetNodeZ(element_nodes[i]);
       }
       //Input parameters
       porosity = CMediumProperties::Porosity(index,theta);
       alpha = flowlinearity_model_values[0];
       apperture = porosity / flowlinearity_model_values[1]; //Change equivalent
       porosity to individual fracture porosity Re =
       flowlinearity_model_values[2];

       //Fluid properties
       rho = m_mfp->Density();
       mu  = m_mfp->Viscosity();
       Re = 0.0; //Reynolds number for turbulent flow CMCD 04. 2004
       Lambda = 0.0; //Frictional Resistence

       //Flow status
       hyd_radius = 2*apperture;
       perm = (apperture * apperture)/12.;
       linear_q = (Re*mu)/(hyd_radius*rho); //max linear q

       //Fluid pressure at each node
       for (i = 0; i < nn; i++) {
       pt_element_node[i] = GetNodeVal(element_nodes[i],1);
       zt_element_node[i] = GetNode(element_nodes[i])->z;
       ht_element_node[i] = pt_element_node[i] + ((g * rho) *
       zt_element_node[i]); //In Pascal
       }
       Calc2DElementCoordinatesTriangle(index,xt,yt,dircos); //CMCD included
       03/2004 area = ElGetElementVolume(index)/m_mmp->geo_area; //CMCD March
       2004 removed wrong area in  code
       //Shape function derivatives
       dN_dx[0] = (yt[1] - yt[2]) / (2. * area);
       dN_dx[1] = (yt[2] - yt[0]) / (2. * area);
       dN_dx[2] = (yt[0] - yt[1]) / (2. * area);
       dN_dy[0] = (xt[2] - xt[1]) / (2. * area);
       dN_dy[1] = (xt[0] - xt[2]) / (2. * area);
       dN_dy[2] = (xt[1] - xt[0]) / (2. * area);
       grad_x = MSkalarprodukt(dN_dx, ht_element_node, nn);
       grad_y = MSkalarprodukt(dN_dy, ht_element_node, nn);
       //v2[0] = MSkalarprodukt(dN_dx, zg, nn);
       //v2[1] = MSkalarprodukt(dN_dy, zg, nn);
       //Assume isotropic nonlinear flow (p268 Kolditz 2001)
       linear_q = linear_q/3.0; // Here the whole element is considered hence 4*
       to remove avereaging effects isotropicgradient =
       pow((grad_x*grad_x+grad_y*grad_y),0.5); flow_rate = (perm *
       isotropicgradient)/mu; if (flow_rate > linear_q){ turbulent_q =
       flow_rate-linear_q; linear_grad = (linear_q *apperture * mu)/perm;
       turbulent_grad = isotropicgradient - linear_grad;
       temp = pow((turbulent_grad/(rho*g)),1-alpha)/(turbulent_grad/(rho*g));
       k_rel = ((linear_grad*1.0) +(turbulent_grad*temp))/isotropicgradient;
       }
       else {
       k_rel = 1.0;
       }

       //velovec[0] = (-k_x * k_rel_grad_p * k_rel_S*k_rel / mu)* (v1[0] + (rho
       * g * v2[0]));
       //velovec[1] = (-k_y * k_rel_grad_p * k_rel_S*k_rel / mu) * (v1[1] + (rho
       * g * v2[1]));

       // special stop CMCD

       //		if (index == 7021){
       //			printf("\n");
       //			printf("Element 4516 k_rel = %g\n",k_rel);
       //			}
       //
       //		break;
       }
     */
    return k_rel;
}

/**************************************************************************
   ROCKFLOW - Funktion: Storage Function

   Aufgabe:
   Berechnet Speicherkoeffizienten

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E long index: Elementnummer

   Ergebnis:
   rel. Permeabilitaet

   Programmaenderungen:
   1/2001   C.Thorenz   Erste Version
   7/2007   C.Thorenz   Div. Druckabhaengigkeiten
   06/2003 OK/CM case 4: Storage as function of effective stress read from curve
   11/2003   CMCD		Geothermal methods added (20)

   Case overview
   0		Curve
   1		Constant
   2		Funktion der effektiven Spannung und des Drucks in Elementmitte
   3		Funktion der effektiven Spannung und des Drucks in Elementmitte
ueber Kurve 4		Storage as function of effective stress read from curve 5
Storage as normal stress in element in stress field defined by KTB stress field.
6 Storage as normal stress in element in stress field defined by KTB stress
field, function to increase storage with distance from borehole.
**************************************************************************/
double CMediumProperties::StorageFunction(long index, double* gp, double theta)
{
    // OK411
    theta = theta;
    gp = gp;
    index = index;

    switch (storage_model)
    {
        case 0:

            // OK411 storage = GetCurveValue((int) storage_model_values[0], 0,
            // primary_variable[0], &idummy);

            break;

        case 1:
            // Konstanter Wert
            storage = storage_model_values[0];
            break;

        case 2:
// Funktion der effektiven Spannung und des Drucks in Elementmitte
#ifdef obsolete  // WW. 06.11.2008
            // Den Druck holen
            p = primary_variable[0];

            /* Mittlere Tiefe */
            nn = ElNumberOfNodes[ElGetElementType(index) - 1];
            element_nodes = ElGetElementNodes(index);

            for (i = 0; i < nn; i++)
                z[i] = GetNodeZ(element_nodes[i]);

            /* Spannung = sigma(z0) + d_sigma/d_z*z */
            // OKsigma = storage_model_values[2] + storage_model_values[3] *
            // InterpolValueVector(index, z, 0., 0., 0.);
            sigma =
                storage_model_values[2] +
                storage_model_values[3] *
                    InterpolValueVector(ElGetElementType(index), z, 0., 0., 0.);

            /* Auf effektive Spannung umrechnen */
            sigma -= p;

            storage = exp(storage_model_values[0] -
                          storage_model_values[1] * log(sigma));
#endif  //#ifdef obsolete //WW. 06.11.2008
            break;

        case 3:
/* Funktion der effektiven Spannung und des Drucks in Elementmitte
           ueber Kurve */
#ifdef obsolete  // WW. 06.11.2008
            /* Den Druck holen */
            p = primary_variable[0];

            /* Mittlere Tiefe */
            nn = ElNumberOfNodes[ElGetElementType(index) - 1];
            element_nodes = ElGetElementNodes(index);

            for (i = 0; i < nn; i++)
                z[i] = GetNodeZ(element_nodes[i]);

            /* Spannung = sigma(z0) + d_sigma/d_z*z */
            sigma =
                storage_model_values[1] +
                storage_model_values[2] *
                    InterpolValueVector(ElGetElementType(index), z, 0., 0., 0.);

            /* Auf effektive Spannung umrechnen */
            sigma -= p;

            storage = GetCurveValue((int)storage_model_values[0], 0, sigma, &i);
#endif
            break;

        case 4: /* McD Storage as function of effective stress read from curve
                 */
        /*OK411
              CalculateSimpleMiddelPointElement(index, coords);
                p = primary_variable[0];
               density_solid = storage_model_values[2];
              stress_eff = (fabs(coords[2])*gravity_constant*density_solid) - p;
                storage =GetCurveValue((int) storage_model_values[0], 0,
           stress_eff, &idummy); break;
         */
        case 5: /* Stroage : Normal stress calculated according to the
                   orientation of the fracture element*/
        /*OK411
              Type=ElGetElementType(index);
              material_group = ElGetElementGroupNumber(index);

                 if (Type == 2||Type == 3||Type == 4)  //Function defined for
           square, triangular and cubic elements
                 {
                    nn = ElNumberOfNodes[Type - 1];
                    element_nodes = ElGetElementNodes(index);

                 for (i=0;i<nn;i++)
                 {
           zelnodes[i]=GetNodeZ(element_nodes[i]);
           yelnodes[i]=GetNodeY(element_nodes[i]);
           xelnodes[i]=GetNodeX(element_nodes[i]);
           }

           Pie=3.141592654;
           angle=(storage_model_values[3]*Pie)/180.;
           for (i=0;i<nn;i++)
           {
           znodes[i]=zelnodes[i];
           xnodes[i]=xelnodes[i]*cos(angle)-yelnodes[i]*sin(angle);
           ynodes[i]=xelnodes[i]*sin(angle)+yelnodes[i]*cos(angle);

           }

           if (Type == 2) //Square
           {
           a1=xnodes[2]-xnodes[0];
           a2=ynodes[2]-ynodes[0];
           a3=znodes[2]-znodes[0];
           b1=xnodes[3]-xnodes[1];
           b2=ynodes[3]-ynodes[1];
           b3=znodes[3]-znodes[1];
           normx=a2*b3-a3*b2;
           normy=a3*b1-a1*b3;
           normz=a1*b2-a2*b1;
           normlen=sqrt(pow(normx,2.0)+pow(normy,2.0)+pow(normz,2.0));
           dircosl= normy/normlen;
           dircosm= normz/normlen;
           dircosn= normx/normlen;
           }

           if (Type == 3) //Cube
           {
           a1=xnodes[2]-xnodes[0];
           a2=ynodes[2]-ynodes[0];
           a3=znodes[2]-znodes[0];
           b1=xnodes[3]-xnodes[1];
           b2=ynodes[3]-ynodes[1];
           b3=znodes[3]-znodes[1];
           normx=a2*b3-a3*b2;
           normy=a3*b1-a1*b3;
           normz=a1*b2-a2*b1;
           normlen=sqrt(pow(normx,2.0)+pow(normy,2.0)+pow(normz,2.0));
           dircosl= normy/normlen;
           dircosm= normz/normlen;
           dircosn= normx/normlen;
           }

           if (Type == 4) //Triangle
           {
           a1=xnodes[1]-xnodes[0];
           a2=ynodes[1]-ynodes[0];
           a3=znodes[1]-znodes[0];
           b1=xnodes[2]-xnodes[1];
           b2=ynodes[2]-ynodes[1];
           b3=znodes[2]-znodes[1];
           normx=a2*b3-a3*b2;
           normy=a3*b1-a1*b3;
           normz=a1*b2-a2*b1;
           normlen=sqrt(pow(normx,2.0)+pow(normy,2.0)+pow(normz,2.0));
           dircosl= normy/normlen;
           dircosm= normz/normlen;
           dircosn= normx/normlen;
           }
           // Calculate average location of the element
           CalculateSimpleMiddelPointElement(index, coords);
           x_mid = coords[0];
           y_mid = coords[1];
           z_mid = coords[2];

           p = primary_variable[0];

           //Calcualtion of stress according to Ito & Zoback 2000 for KTB hole

           sigma1=z_mid*0.045*-1e6;
           sigma2=z_mid*0.028*-1e6;
           sigma3=z_mid*0.02*-1e6;

           ///Calculate total normal stress on element
           //Note in this case sigma2 corresponds to the vertical stress
           tot_norm_stress=sigma1*dircosl*dircosl+sigma2*dircosm*dircosm+sigma3*dircosn*dircosn;

           //Calculate normal effective stress
           eff_norm_stress=tot_norm_stress-p;

           // special stop CMCD
           if (eff_norm_stress>220000000.){
           phase=0;
           }

           //Take value of storage from a curve
           S=GetCurveValue((int) storage_model_values[0], 0, eff_norm_stress,
           &idummy);
           }

           // default value if element type is not included in the method for
           calculating the normal else S=storage_model_values[2]; storage = S;
           break;
         */
        case 6: /* Normal stress calculated according to the orientation of the
                   fracture element*/
            /*OK411
                  Type=ElGetElementType(index);
                  material_group = ElGetElementGroupNumber(index);

                  if (material_group == 0)
                  {
                     if (Type == 2||Type == 3||Type == 4)
                     {
                        nn = ElNumberOfNodes[Type - 1];
                        element_nodes = ElGetElementNodes(index);

               for (i=0;i<nn;i++)
               {
               zelnodes[i]=GetNodeZ(element_nodes[i]);
               yelnodes[i]=GetNodeY(element_nodes[i]);
               xelnodes[i]=GetNodeX(element_nodes[i]);
               }

               Pie=3.141592654;
               angle=(storage_model_values[3]*Pie)/180.;
               for (i=0;i<nn;i++)
               {
               znodes[i]=zelnodes[i];
               xnodes[i]=xelnodes[i]*cos(angle)-yelnodes[i]*sin(angle);
               ynodes[i]=xelnodes[i]*sin(angle)+yelnodes[i]*cos(angle);

               }

               if (Type == 2)
               {
               a1=xnodes[2]-xnodes[0];
               a2=ynodes[2]-ynodes[0];
               a3=znodes[2]-znodes[0];
               b1=xnodes[3]-xnodes[1];
               b2=ynodes[3]-ynodes[1];
               b3=znodes[3]-znodes[1];
               normx=a2*b3-a3*b2;
               normy=a3*b1-a1*b3;
               normz=a1*b2-a2*b1;
               normlen=sqrt(pow(normx,2.0)+pow(normy,2.0)+pow(normz,2.0));
               dircosl= normy/normlen;
               dircosm= normz/normlen;
               dircosn= normx/normlen;
               }

               if (Type == 3)
               {
               a1=xnodes[2]-xnodes[0];
               a2=ynodes[2]-ynodes[0];
               a3=znodes[2]-znodes[0];
               b1=xnodes[3]-xnodes[1];
               b2=ynodes[3]-ynodes[1];
               b3=znodes[3]-znodes[1];
               normx=a2*b3-a3*b2;
               normy=a3*b1-a1*b3;
               normz=a1*b2-a2*b1;
               normlen=sqrt(pow(normx,2.0)+pow(normy,2.0)+pow(normz,2.0));
               dircosl= normy/normlen;
               dircosm= normz/normlen;
               dircosn= normx/normlen;
               }

               if (Type == 4)
               {
               a1=xnodes[1]-xnodes[0];
               a2=ynodes[1]-ynodes[0];
               a3=znodes[1]-znodes[0];
               b1=xnodes[2]-xnodes[1];
               b2=ynodes[2]-ynodes[1];
               b3=znodes[2]-znodes[1];
               normx=a2*b3-a3*b2;
               normy=a3*b1-a1*b3;
               normz=a1*b2-a2*b1;
               normlen=sqrt(pow(normx,2.0)+pow(normy,2.0)+pow(normz,2.0));
               dircosl= normy/normlen;
               dircosm= normz/normlen;
               dircosn= normx/normlen;
               }

               CalculateSimpleMiddelPointElement(index,coords);
               x_mid = coords[0];
               y_mid = coords[1];
               z_mid = coords[2];

               p = primary_variable[0];

               x_bore = storage_model_values[4];
               y_bore = storage_model_values[5];
               z_bore = storage_model_values[6];

               distance = pow(pow((x_mid-x_bore),2.0) + pow((y_mid-y_bore),2.0)
               + pow((z_mid-z_bore),2.0),0.5);

               sigma1=z_mid*0.045*-1e6;
               sigma2=z_mid*0.028*-1e6;
               sigma3=z_mid*0.02*-1e6;

               ///Calculate total normal stress on element
               tot_norm_stress=sigma1*dircosl*dircosl+sigma2*dircosm*dircosm+sigma3*dircosn*dircosn;

               eff_norm_stress=tot_norm_stress-p;

               S=GetCurveValue((int) storage_model_values[0], 0,
               eff_norm_stress, &idummy); if (distance >
               storage_model_values[7]){ distance = storage_model_values[7];
               }
               if (distance > 2) S = S + (S * (distance-2) *
               storage_model_values[8]);
               }

               else S=storage_model_values[2];
               }

               else S=storage_model_values[2];
               storage = S;
             */
            break;
        case 7:  // poroelasticity RW
        {
            // Moved the following comment from double
            // CFiniteElementStd::CalCoefMass() to here by WW JT 2010, needed
            // storage term and fluid compressibility... We derive here the
            // storage at constant strain, or the inverse of Biot's "M"
            // coefficient Assumptions are the most general possible::
            // Invarience under "pi" (Detournay & Cheng) loading. Se = 1/M =
            // poro/Kf + (alpha-poro)/Ks    ::    Cf = 1/Kf = 1/rho * drho/dp ::
            // alpha = 1 - K/Ks Second term (of Se) below vanishes for
            // incompressible grains

            SolidProp::CSolidProperties* solid_prop =
                msp_vector[Fem_Ele_Std->GetMeshElement()->GetPatchIndex()];
            const double biots_constant = solid_prop->getBiotsConstant();
            const double porosity = Porosity(index, theta);
            return (biots_constant - porosity) * (1.0 - biots_constant) /
                   solid_prop->getBulkModulus();
        }
        case 10:
            if (permeability_saturation_model[0] == 10)  // MW
                storage = porosity_model_values[0] /
                          (gravity_constant * gravity_constant *
                           mfp_vector[0]->Density());
            // MW I have no idea, why I need 1/(g^2*rho) here; it should only be
            // 1/(g*rho) maybe, the mass term in richards flow has been
            // normalized on g ???
            else
                std::cout << "Wrong PERMEABILITY_SATURATION_MODEL for STORAGE "
                             "model 10."
                          << std::endl;
            break;
        case 11:
            if (m_pcs->getProcessType() == FiniteElement::LIQUID_FLOW)
                storage = porosity_model_values[0] /
                          (gravity_constant * mfp_vector[0]->Density());
            else
                std::cout << "Wrong process type for STORAGE model 11 (only "
                             "intended for LIQUID_FLOW)."
                          << std::endl;
            break;
        default:
            storage = 0.0;  // OK DisplayMsgLn("The requested storativity model
                            // is unknown!!!");
            break;
    }
    return storage;
}

// AS:08.2012 storage function of effective stress
double CMediumProperties::StorageFunctionEffStress(long index, int nnodes,
                                                   CFiniteElementStd* h_fem)
{
    int i, j, size;
    double storage_stress = 1.0;
    // calculate principal effective stress
    double stress[6] = {0.}, prin_str[6] = {0.}, prin_dir[9] = {0.};
    int stress_index[6];
    // model dimension
    int dim = h_fem->dm_pcs->m_msh->GetCoordinateFlag() / 10;
    size = 6;
    if (dim == 2)
    {
        size = 4;
        if (h_fem->axisymmetry)
            dim = 3;
    }
    stress_index[0] = h_fem->dm_pcs->GetNodeValueIndex("STRESS_XX");
    stress_index[1] = h_fem->dm_pcs->GetNodeValueIndex("STRESS_YY");
    stress_index[2] = h_fem->dm_pcs->GetNodeValueIndex("STRESS_ZZ");
    stress_index[3] = h_fem->dm_pcs->GetNodeValueIndex("STRESS_XY");
    if (size == 6)
    {
        stress_index[4] = h_fem->dm_pcs->GetNodeValueIndex("STRESS_XZ");
        stress_index[5] = h_fem->dm_pcs->GetNodeValueIndex("STRESS_YZ");
    }
    double stress_nodes[20] = {0.};

    switch (storage_effstress_model)
    {
        // AS: case 1, permeability directly calculated from curve. 16.08.2012
        case 1:
            for (j = 0; j < size; j++)
            {
                for (i = 0; i < nnodes; i++)
                    stress_nodes[i] = h_fem->dm_pcs->GetNodeValue(
                        h_fem->dm_pcs->m_msh->ele_vector[index]
                            ->getNodeIndices()[i],
                        stress_index[j]);
                stress[j] = h_fem->interpolate(stress_nodes);
            }
            h_fem->SolidProp->CalPrinStrDir(stress, prin_str, prin_dir, dim);
            // permeability from curve with minimum (absolute value) principal
            // effective stress as input
            storage_stress = GetCurveValue(
                (int)storage_effstress_model_value[0], 0, prin_str[0], &i);
            break;
        default:
            break;
    }
    return storage_stress;
}

double CMediumProperties::PermeabilityPressureFunction(long index, double* gp,
                                                       double theta)
{
    double k_rel = 0.0;
    // OK411
    theta = theta;
    gp = gp;
    index = index;

#ifdef obsolete  // OK411
    int nn, i, idummy, p_idx1;
    long* element_nodes;
    static double gh, p, eins_durch_rho_g, sigma, z[8], h[8], grad_h[3];
    static double invjac[8], detjac, grad_omega[8];
    static double mult[3];
    double x_mid, y_mid, z_mid, coords[3];
    double density_solid, stress_eff;
    //	double porosity, factora, factorb,valuelogk;
    CFluidProperties* m_mfp = NULL;

    // Collect primary variables
    static int nidx0, nidx1;
    double primary_variable[10];  // OK To Do
    int count_nodes;
    int no_pcs_names = (int)pcs_name_vector.size();
    for (i = 0; i < no_pcs_names; i++)
    {
        nidx0 = PCSGetNODValueIndex(pcs_name_vector[i], 0);
        nidx1 = PCSGetNODValueIndex(pcs_name_vector[i], 1);
        if (mode == 0)  // Gauss point values

            primary_variable[i] =
                (1. - theta) *
                    InterpolValue(number, nidx0, gp[0], gp[1], gp[2]) +
                theta * InterpolValue(number, nidx1, gp[0], gp[1], gp[2]);
        else if (mode == 1)  // Node values

            primary_variable[i] = (1. - theta) * GetNodeVal(number, nidx0) +
                                  theta * GetNodeVal(number, nidx1);
        else if (mode == 2)  // Element average value
        {
            count_nodes = ElNumberOfNodes[ElGetElementType(number) - 1];
            element_nodes = ElGetElementNodes(number);
            for (i = 0; i < count_nodes; i++)
                primary_variable[i] += GetNodeVal(element_nodes[i], nidx1);
            primary_variable[i] /= count_nodes;
        }
    }
    switch (permeability_pressure_model)
    {
        case 0:  // Curve function
            k_rel = 1.0;
            break;
        case 1:  // No functional dependence
            k_rel = 1.0;
            break;
        /* Funktion der effektiven Spannung */
        case 2:
#ifdef obsolete  // WW. 06.11.2008
            /* Den Druck holen */
            p = primary_variable[0];
            /* Mittlere Tiefe */
            nn = ElNumberOfNodes[ElGetElementType(index) - 1];
            element_nodes = ElGetElementNodes(index);
            for (i = 0; i < nn; i++)
                z[i] = GetNodeZ(element_nodes[i]);
            /* Spannung = sigma(z0) + d_sigma/d_z*z */
            sigma =
                permeability_pressure_model_values[2] +
                permeability_pressure_model_values[3] *
                    InterpolValueVector(ElGetElementType(index), z, 0., 0., 0.);
            /* Auf effektive Spannung umrechnen */
            sigma -= p;
            k_rel = exp(permeability_pressure_model_values[0] -
                        permeability_pressure_model_values[1] * log(sigma));
#endif
            break;
        case 3: /* Turbulentes Fliessen */
            /* k_rel = max(min((grad(h)*alpha1)^(alpha2-1), alpha4),alpha3) */
            m_mfp = MFPGet("LIQUID");
            // YDGetFluidDensity(0, index, 0., 0., 0., 1.);
            eins_durch_rho_g = 1. / gravity_constant / m_mfp->Density();
            nn = ElNumberOfNodes[ElGetElementType(index) - 1];
            element_nodes = ElGetElementNodes(index);
            p_idx1 = PCSGetNODValueIndex("PRESSURE1", 1);
            for (i = 0; i < nn; i++)
                h[i] = GetNodeVal(element_nodes[i], p_idx1) * eins_durch_rho_g +
                       GetNodeZ(element_nodes[i]);
            switch (ElGetElementType(index))
            {
                default:
                    DisplayMsgLn("Error in GetSoilRelPermPress!");
                    DisplayMsgLn("  Nonlinear permeability not available!");
                    abort();
                case 2:
                    Calc2DElementJacobiMatrix(index, 0.0, 0.0, invjac, &detjac);
                    MGradOmega2D(grad_omega, 0, 0); /* Gradientenmatrix */
                    MMultMatVec(grad_omega, 2, 4, h, 4, mult, 2);
                    MMultVecMat(mult, 2, invjac, 2, 2, grad_h, 2);
                    gh = sqrt(grad_h[0] * grad_h[0] + grad_h[1] * grad_h[1]) *
                         permeability_pressure_model_values[0];
                    if (gh < MKleinsteZahl)
                        k_rel = 1.0;
                    else
                        k_rel = max(
                            min(pow(gh, permeability_pressure_model_values[1] -
                                            1.0),
                                permeability_pressure_model_values[3]),
                            permeability_pressure_model_values[2]);
            }
        case 4:
#ifdef obsolete  // WW. 06.11.2008
            /* Funktion der effektiven Spannung ueber Kurve */
            /* Den Druck holen */
            p = primary_variable[0];
            /* Mittlere Tiefe */
            nn = ElNumberOfNodes[ElGetElementType(index) - 1];
            element_nodes = ElGetElementNodes(index);
            for (i = 0; i < nn; i++)
                z[i] = GetNodeZ(element_nodes[i]);
            /* Spannung = sigma(z0) + d_sigma/d_z*z */
            sigma =
                permeability_pressure_model_values[1] +
                permeability_pressure_model_values[2] *
                    InterpolValueVector(ElGetElementType(index), z, 0., 0., 0.);
            /* Auf effektive Spannung umrechnen */
            sigma -= p;
            k_rel = GetCurveValue((int)permeability_pressure_model_values[0], 0,
                                  sigma, &i);
#endif
            break;
        case 5:
            /* Funktion der effektiven Spannung ueber Kurve CMCD 26.06.2003*/
            /* Average depth */
            CalculateSimpleMiddelPointElement(index, coords);
            x_mid = coords[0];
            y_mid = coords[1];
            z_mid = coords[2];
            p = primary_variable[0];
            density_solid = permeability_pressure_model_values[2];
            stress_eff = (fabs(z_mid) * gravity_constant * density_solid) - p;
            k_rel = GetCurveValue((int)permeability_pressure_model_values[0], 0,
                                  stress_eff, &idummy);
            break;
        case 6:
            k_rel =
                PermeabilityPressureFunctionMethod1(index, primary_variable[0]);
            break;
        case 7:
            k_rel =
                PermeabilityPressureFunctionMethod2(index, primary_variable[0]);
            break;
        case 8:
            k_rel =
                PermeabilityPressureFunctionMethod3(index, primary_variable[0]);
            break;
        case 9:
            k_rel = PermeabilityPressureFunctionMethod4(
                index, primary_variable[0], primary_variable[1]);
            break;
        default:          // CMCD Einbau
            k_rel = 1.0;  // CMCD Einbau
            break;        // CMCD Einbau
    }
#endif
    return k_rel;
}

/**************************************************************************
   12.(ii)a Subfunction of Permeability_Function_Pressure
   ROCKFLOW - Funktion: MATCalcPressurePermeabilityMethod1
   Application to KTB
   Aufgabe:
   Calculates relative permeability from
   the normal stress according to orientation of the fractures in the KTB site
system converts the normal stress to effective stress, reads from a curve what
the permeability of a fracture under the given effective stress is.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E: long index: Elementnummer, double primary_variable pressure
   R: relative permeability

   Ergebnis:
   rel. Permeabilitaet

   Programmaenderungen:
   09/2004   CMCD Inclusion in GeoSys 4
**************************************************************************/
double CMediumProperties::PermeabilityPressureFunctionMethod1(long index,
                                                              double pressure)
{
    // OK411
    index = index;
    pressure = pressure;

#ifdef obsolete  // OK411
    double R, Pie, p;
    double znodes[8], ynodes[8], xnodes[8], angle;
    double x_mid, y_mid, z_mid, coords[3];
    double zelnodes[8], yelnodes[8], xelnodes[8];
    double sigma1, sigma2, sigma3;
    double a1, a2, a3, b1, b2, b3;
    double normx, normy, normz, normlen;
    double dircosl, dircosm, dircosn;
    double tot_norm_stress, eff_norm_stress;
    int material_group, Type;
    int phase, nn, i, idummy;
    long* element_nodes;
    dircosl = dircosm = dircosn = 0.0;
    Type = ElGetElementType(index);
    material_group = ElGetElementGroupNumber(index);

    /* special stop CMCD*/
    if (index == 6590)
        phase = 0;
    /* Stop here*/

    if (Type == 2 || Type == 3 ||
        Type ==
            4) /*Function defined for square, triangular and cubic elements*/
    {
        nn = ElNumberOfNodes[Type - 1];
        element_nodes = ElGetElementNodes(index);

        /* Calculate directional cosins, note that this is currently set up*/
        /* Sigma1 is in the y direction*/
        /* Sigma2 is in the z direction*/
        /* Sigma3 is in the x direction*/
        /* This correspondes approximately to the KTB site conditions*/

        for (i = 0; i < nn; i++)
        {
            zelnodes[i] = GetNodeZ(element_nodes[i]);
            yelnodes[i] = GetNodeY(element_nodes[i]);
            xelnodes[i] = GetNodeX(element_nodes[i]);
        }

        /*Coordinate transformation to match sigma max direction*/
        /* y direction matches the north direction*/

        Pie = 3.141592654;
        for (i = 0; i < nn; i++)
        {
            znodes[i] = zelnodes[i];
            angle = (permeability_pressure_model_values[3] * Pie) / 180.;
            xnodes[i] = xelnodes[i] * cos(angle) - yelnodes[i] * sin(angle);
            ynodes[i] = xelnodes[i] * sin(angle) + yelnodes[i] * cos(angle);
        }

        if (Type == 2) /*Square*/
        {
            a1 = xnodes[2] - xnodes[0];
            a2 = ynodes[2] - ynodes[0];
            a3 = znodes[2] - znodes[0];
            b1 = xnodes[3] - xnodes[1];
            b2 = ynodes[3] - ynodes[1];
            b3 = znodes[3] - znodes[1];
            normx = a2 * b3 - a3 * b2;
            normy = a3 * b1 - a1 * b3;
            normz = a1 * b2 - a2 * b1;
            normlen = sqrt(normx * normx + normy * normy + normz * normz);
            dircosl = normy / normlen;
            dircosm = normz / normlen;
            dircosn = normx / normlen;
        }

        if (Type == 3) /*Cube*/
        {
            a1 = xnodes[2] - xnodes[0];
            a2 = ynodes[2] - ynodes[0];
            a3 = znodes[2] - znodes[0];
            b1 = xnodes[3] - xnodes[1];
            b2 = ynodes[3] - ynodes[1];
            b3 = znodes[3] - znodes[1];
            normx = a2 * b3 - a3 * b2;
            normy = a3 * b1 - a1 * b3;
            normz = a1 * b2 - a2 * b1;
            normlen = sqrt(normx * normx + normy * normy + normz * normz);
            dircosl = normy / normlen;
            dircosm = normz / normlen;
            dircosn = normx / normlen;
        }

        if (Type == 4) /*Triangle*/
        {
            a1 = xnodes[1] - xnodes[0];
            a2 = ynodes[1] - ynodes[0];
            a3 = znodes[1] - znodes[0];
            b1 = xnodes[2] - xnodes[1];
            b2 = ynodes[2] - ynodes[1];
            b3 = znodes[2] - znodes[1];
            normx = a2 * b3 - a3 * b2;
            normy = a3 * b1 - a1 * b3;
            normz = a1 * b2 - a2 * b1;
            normlen = sqrt(normx * normx + normy * normy + normz * normz);
            dircosl = normy / normlen;
            dircosm = normz / normlen;
            dircosn = normx / normlen;
        }

        /* Calculation of average location of element*/
        CalculateSimpleMiddelPointElement(index, coords);
        x_mid = coords[0];
        y_mid = coords[1];
        z_mid = coords[2];

        /*Calculate fluid pressure in element*/
        p = pressure;
        /*Calcualtion of stress according to Ito & Zoback 2000 for KTB hole*/

        sigma1 = z_mid * 0.045 * -1e6;
        sigma2 = z_mid * 0.028 * -1e6;
        sigma3 = z_mid * 0.02 * -1e6;

        /*Calculate total normal stress on element
           Note in this case sigma2 corresponds to the vertical stress*/
        tot_norm_stress = sigma1 * dircosl * dircosl +
                          sigma2 * dircosm * dircosm +
                          sigma3 * dircosn * dircosn;

        /*Calculate normal effective stress*/
        eff_norm_stress = tot_norm_stress - p;

        /*Take value of storage from a curve*/
        R = GetCurveValue((int)permeability_pressure_model_values[0], 0,
                          eff_norm_stress, &idummy);

        return R;
    }
    /* default value if element type is not included in the method for
     * calculating the normal */
    else
        R = permeability_pressure_model_values[2];
#endif
    return 0.0;
}

/**************************************************************************
   12.(ii)b Subfunction of Permeability_Function_Pressure
   Function: MATCalcPressurePermeabilityMethod2
   Application to KTB
   Aufgabe:
   Calculates relative permeability from
   the normal stress according to orientation of the fractures in the KTB site
system converts the normal stress to effective stress, reads from a curve what
the permeability of a fracture under the given effective stress is. Permeability
is then related to the distance from the borehole.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E: long index: Elementnummer, double primary_variable pressure
   R: relative permeability

   Ergebnis:
   rel. Permeabilitaet

   Programmaenderungen:
   09/2004   CMCD Inclusion in GeoSys 4
**************************************************************************/
double CMediumProperties::PermeabilityPressureFunctionMethod2(long index,
                                                              double pressure)
{
    // OK411
    index = index;
    pressure = pressure;

#ifdef obsolete  // OK411
    double R, Pie, p;
    double znodes[8], ynodes[8], xnodes[8], angle;
    double x_mid, y_mid, z_mid, coords[3];
    double zelnodes[8], yelnodes[8], xelnodes[8];
    double sigma1, sigma2, sigma3;
    double a1, a2, a3, b1, b2, b3;
    double normx, normy, normz, normlen;
    double dircosl, dircosm, dircosn;
    double tot_norm_stress, eff_norm_stress;
    int material_group, Type;
    double x_bore, y_bore, z_bore, distance;
    int nn, i, idummy;
    long* element_nodes;
    dircosl = dircosm = dircosn = 0.0;
    Type = ElGetElementType(index);
    material_group = ElGetElementGroupNumber(index);

    /* special stop CMCD*/
    if (index == 4516)
    {
        /* Stop here */
    }

    if (Type == 2 || Type == 3 ||
        Type ==
            4) /*Function defined for square, triangular and cubic elements*/
    {
        nn = ElNumberOfNodes[Type - 1];
        element_nodes = ElGetElementNodes(index);

        /* Calculate directional cosins, note that this is currently set up*/
        /* Sigma1 is in the y direction*/
        /* Sigma2 is in the z direction*/
        /* Sigma3 is in the x direction*/
        /* This correspondes approximately to the KTB site conditions*/

        for (i = 0; i < nn; i++)
        {
            zelnodes[i] = GetNodeZ(element_nodes[i]);
            yelnodes[i] = GetNodeY(element_nodes[i]);
            xelnodes[i] = GetNodeX(element_nodes[i]);
        }

        Pie = 3.141592654;
        angle = (permeability_pressure_model_values[3] * Pie) / 180.;
        for (i = 0; i < nn; i++)
        {
            znodes[i] = zelnodes[i];
            xnodes[i] = xelnodes[i] * cos(angle) - yelnodes[i] * sin(angle);
            ynodes[i] = xelnodes[i] * sin(angle) + yelnodes[i] * cos(angle);
        }

        if (Type == 2) /*Square*/
        {
            a1 = xnodes[2] - xnodes[0];
            a2 = ynodes[2] - ynodes[0];
            a3 = znodes[2] - znodes[0];
            b1 = xnodes[3] - xnodes[1];
            b2 = ynodes[3] - ynodes[1];
            b3 = znodes[3] - znodes[1];
            normx = a2 * b3 - a3 * b2;
            normy = a3 * b1 - a1 * b3;
            normz = a1 * b2 - a2 * b1;
            normlen = sqrt(normx * normx + normy * normy + normz * normz);
            dircosl = normy / normlen;
            dircosm = normz / normlen;
            dircosn = normx / normlen;
        }

        if (Type == 3) /*Cube*/
        {
            a1 = xnodes[2] - xnodes[0];
            a2 = ynodes[2] - ynodes[0];
            a3 = znodes[2] - znodes[0];
            b1 = xnodes[3] - xnodes[1];
            b2 = ynodes[3] - ynodes[1];
            b3 = znodes[3] - znodes[1];
            normx = a2 * b3 - a3 * b2;
            normy = a3 * b1 - a1 * b3;
            normz = a1 * b2 - a2 * b1;
            normlen = sqrt(normx * normx + normy * normy + normz * normz);
            dircosl = normy / normlen;
            dircosm = normz / normlen;
            dircosn = normx / normlen;
        }

        if (Type == 4) /*Triangle*/
        {
            a1 = xnodes[1] - xnodes[0];
            a2 = ynodes[1] - ynodes[0];
            a3 = znodes[1] - znodes[0];
            b1 = xnodes[2] - xnodes[1];
            b2 = ynodes[2] - ynodes[1];
            b3 = znodes[2] - znodes[1];
            normx = a2 * b3 - a3 * b2;
            normy = a3 * b1 - a1 * b3;
            normz = a1 * b2 - a2 * b1;
            normlen = sqrt(normx * normx + normy * normy + normz * normz);
            dircosl = normy / normlen;
            dircosm = normz / normlen;
            dircosn = normx / normlen;
        }

        CalculateSimpleMiddelPointElement(index, coords);
        x_mid = coords[0];
        y_mid = coords[1];
        z_mid = coords[2];

        /*Calculate fluid pressure in element*/
        p = pressure;

        /*Calculate distance from borehole*/
        x_bore = permeability_pressure_model_values[4];
        y_bore = permeability_pressure_model_values[5];
        z_bore = permeability_pressure_model_values[6];

        distance = sqrt((x_mid - x_bore) * (x_mid - x_bore) +
                        (y_mid - y_bore) * (y_mid - y_bore) +
                        (z_mid - z_bore) * (z_mid - z_bore));

        /*Calcualtion of stress according to Ito & Zoback 2000 for KTB hole*/

        sigma1 = z_mid * 0.045 * -1e6;
        sigma2 = z_mid * 0.028 * -1e6;
        sigma3 = z_mid * 0.02 * -1e6;

        /// Calculate total normal stress on element
        /*Note in this case sigma2 corresponds to the vertical stress*/
        tot_norm_stress = sigma1 * dircosl * dircosl +
                          sigma2 * dircosm * dircosm +
                          sigma3 * dircosn * dircosn;

        /*Calculate normal effective stress*/
        eff_norm_stress = tot_norm_stress - p;

        /*Take value of storage from a curve*/
        R = GetCurveValue((int)permeability_pressure_model_values[0], 0,
                          eff_norm_stress, &idummy);
        if (distance > permeability_pressure_model_values[7])
            distance = permeability_pressure_model_values[7];

        if (distance > 2)
            R = R +
                (R * (distance - 2) * permeability_pressure_model_values[8]);

        return R;
    }
    /* default value if element type is not included in the method for
     * calculating the normal */
    else
        R = permeability_pressure_model_values[2];
#endif
    return 0.0;
}

/**************************************************************************
   12.(ii)c Subfunction of Permeability_Function_Pressure
   Function: MATCalcPressurePermeabilityMethod3
   Application to Urach

   Aufgabe:
   The normal stress across the fractures is calculated according to an
approximate formulation from the relationship of stress with depth. This normal
stress is then converted into a permeabilty by reference to effective stress and
a curve.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E: long index: Elementnummer, double primary_variable pressure
   R: relative permeability

   Ergebnis:
   rel. Permeabilitaet

   Programmaenderungen:
   09/2004   CMCD Inclusion in GeoSys 4
**************************************************************************/
double CMediumProperties::PermeabilityPressureFunctionMethod3(long index,
                                                              double pressure)
{
    // OK411
    index = index;
    pressure = pressure;

    double perm = 0.0;
#ifdef obsolete  // OK411
    /* Funktion der effektiven Spannung ueber Kurve CMCD 26.06.2003*/
    double x_mid, y_mid, z_mid, coords[3];
    double stress_eff, p, perm;
    int idummy;
    /* Average depth */
    CalculateSimpleMiddelPointElement(index, coords);
    x_mid = coords[0];
    y_mid = coords[1];
    z_mid = coords[2];

    /*Calculate fluid pressure in element*/
    p = pressure;
    stress_eff =
        (fabs(z_mid) * permeability_pressure_model_values[1] * 1e6) - p;
    perm = GetCurveValue((int)permeability_pressure_model_values[0], 0,
                         stress_eff, &idummy);
#endif
    return perm;
}

/**************************************************************************
   12.(ii)d Subfunction of Permeability_Function_Pressure
   Function: MATCalcPressurePermeabilityMethod4
   Application to Urach

   Aufgabe:
   The normal stress across the fractures is calculated according to an
approximate formulation from the relationship of stress with depth. This normal
stress is then adjusted to take account of thermal cooling, and the resulting
effective stress across the fracture isconverted into a permeabilty. by
reference to a curve.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E: long index: Elementnummer, double primary_variable pressure
   R: relative permeability

   Ergebnis:
   rel. Permeabilitaet

   Programmaenderungen:
   09/2004   CMCD Inclusion in GeoSys 4
**************************************************************************/
double CMediumProperties::PermeabilityPressureFunctionMethod4(
    long index, double pressure, double temperature)
{
    // OK411
    index = index;
    pressure = pressure;
    temperature = temperature;
    double perm = 0.0;

#ifdef obsolete  // OK411
    /* Funktion der effektiven Spannung ueber Kurve CMCD 26.06.2003*/
    double x_mid, y_mid, z_mid, coords[3];
    double stress_eff, p, perm, thermal_stress;
    int idummy;
    /* Average depth */
    CalculateSimpleMiddelPointElement(index, coords);
    x_mid = coords[0];
    y_mid = coords[1];
    z_mid = coords[2];

    /*Calculate effective stress in the element*/
    p = pressure;
    stress_eff =
        (fabs(z_mid) * permeability_pressure_model_values[1] * 1e6) - p;

    /*Impact of thermal stress*/
    thermal_stress = (permeability_pressure_model_values[2] - temperature) *
                     permeability_pressure_model_values[3] *
                     permeability_pressure_model_values[4];

    stress_eff = stress_eff - thermal_stress;
    if (stress_eff < 0.0)
        stress_eff = 0.0;
    /*Read value effective stress against curve*/
    perm = GetCurveValue((int)permeability_pressure_model_values[0], 0,
                         stress_eff, &idummy);
#endif
    return perm;
}

double CMediumProperties::PermeabilityPorosityFunction(long number, double* gp,
                                                       double theta)
{
    double k_rel = 0.0;
    // OK411
    theta = theta;
    gp = gp;
    number = number;
#ifdef obsolete  // OK411
    int i;
    long* element_nodes;
    double factora, factorb, valuelogk;

    // Collect primary variables
    static int nidx0, nidx1;
    double primary_variable[10];  // OK To Do
    int count_nodes;
    int no_pcs_names = (int)pcs_name_vector.size();
    for (i = 0; i < no_pcs_names; i++)
    {
        nidx0 = PCSGetNODValueIndex(pcs_name_vector[i], 0);
        nidx1 = PCSGetNODValueIndex(pcs_name_vector[i], 1);
        if (mode == 0)  // Gauss point values

            primary_variable[i] =
                (1. - theta) *
                    InterpolValue(number, nidx0, gp[0], gp[1], gp[2]) +
                theta * InterpolValue(number, nidx1, gp[0], gp[1], gp[2]);
        else if (mode == 1)  // Node values

            primary_variable[i] = (1. - theta) * GetNodeVal(number, nidx0) +
                                  theta * GetNodeVal(number, nidx1);
        else if (mode == 2)  // Element average value
        {
            count_nodes = ElNumberOfNodes[ElGetElementType(number) - 1];
            element_nodes = ElGetElementNodes(number);
            for (i = 0; i < count_nodes; i++)
                primary_variable[i] += GetNodeVal(element_nodes[i], nidx1);
            primary_variable[i] /= count_nodes;
        }
    }
    switch (permeability_porosity_model)
    {
        case 0:  // Reserved for a curve function
            k_rel = 1.0;
            break;
        case 1:  // Constant value function
            k_rel = 1.0;
            break;
        case 2:  // Ming Lian function
            factora = permeability_porosity_model_values[0];
            factorb = permeability_porosity_model_values[1];
            valuelogk = factora + (factorb * porosity);
            k_rel = 0.987e-12 * pow(10., valuelogk);
            CMediumProperties* m_mmp = NULL;
            long group = ElGetElementGroupNumber(number);
            m_mmp = mmp_vector[group];
            double* k_ij;
            k_ij = m_mmp->PermeabilityTensor(number);  // permeability;
            k_rel /= k_ij[0];
            break;
    }
#endif
    return k_rel;
}
double CMediumProperties::ParticleDiameter()
{
    double val = .0;
    switch (particle_diameter_model)
    {
        case 0:
            break;
        case 1:
            val = particle_diameter_model_value;
            break;
        default:
            break;
    }

    return val;
}

double CMediumProperties::HeatTransferCoefficient(long number, double theta,
                                                  CFiniteElementStd* assem)
{
    double val = .0;
    switch (heat_transfer_model)
    {
        case 0:
            break;
        case 1:
            val = heat_transfer_model_value;
            break;
        case 2:
        {
            const double dp = ParticleDiameter();
            assert(dp > .0);
            double h = .0;
            if (effective_heat_transfer_model == 1)
            {
                h = heat_transfer_model_value;
            }
            else if (effective_heat_transfer_model == 2)
            {
                h = .0;
            }
            else
            {
                std::cout << "***Error: Effective heat transfer model is not "
                             "supported. "
                          << effective_heat_transfer_model << '\n';
            }
            const double& n = Porosity(number, theta);
            const int dimen = assem->pcs->m_msh->GetCoordinateFlag() / 10;
            for (int i = 0; i < dimen * dimen; i++)
                heat_conductivity_tensor[i] = 0.0;
            assem->SolidProp->HeatConductivityTensor(
                dimen, heat_conductivity_tensor,
                assem->MeshElement->GetPatchIndex());
            const double lamda_s =
                heat_conductivity_tensor[0];  // assume isotropic
            val = 6. * (1 - n) / dp * 1. /
                  (1. / h + dp / 10.0 / lamda_s);  // confirmed with DLR people.
                                                   // typo in Schaube2011
        }
        break;
        default:
            break;
    }

    return val;
}
// TN - added for TNEQ/TEQ process
void CMediumProperties::setFrictionPhase(
    FiniteElement::FrictionPhase fric_phase)
{
    _fric_phase = fric_phase;
}

FiniteElement::FrictionPhase CMediumProperties::getFrictionPhase() const
{
    return _fric_phase;
}
