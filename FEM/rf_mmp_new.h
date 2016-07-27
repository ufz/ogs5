/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
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

#ifndef rfmat_mp_new_INC
#define rfmat_mp_new_INC
/* Schutz gegen mehrfaches Einfuegen */

// C++ STL
//#include <list>
//#include <string>
//#include <vector>
//#include <fstream>

// GeoLib
#include "GeoType.h"
#include "makros.h" // JT

// PCSLib
#include "rf_pcs.h"

#include "BHEAbstract.h"
#include "BHE_Net_ELE_HeatPump.h"
using namespace BHE; 

namespace FiniteElement
{
class CFiniteElementStd;
}
using FiniteElement::CFiniteElementStd;
class CMediumProperties
{
public:
	CFiniteElementStd* Fem_Ele_Std;

private:
	// WW
	friend class FiniteElement::CFiniteElementStd;
	// Data base
	void SetMediumPropertiesDefaultsClay(void); // CMCD 9/2004 GeoSys 4
	void SetMediumPropertiesDefaultsSilt(void); // CMCD 9/2004 GeoSys 4
	void SetMediumPropertiesDefaultsSand(void); // CMCD 9/2004 GeoSys 4
	// CMCD 9/2004 GeoSys 4
	void SetMediumPropertiesDefaultsGravel(void);
	// CMCD 9/2004 GeoSys 4
	void SetMediumPropertiesDefaultsCrystalline(void);
	// CMCD 9/2004 GeoSys 4
	void SetMediumPropertiesDefaultsBordenAquifer(void);
	// Porosity
	// CMCD 9/2004 GeoSys 4
	double PorosityEffectiveStress(long, double);
	double PorosityVolumetricFreeSwellingConstantIonicstrength(long, double, double);
	// MX 1/2005
	double PorosityEffectiveConstrainedSwelling(long, double, double, double*);
	// MX 1/2005
	double PorosityVolumetricFreeSwelling(long, double, double);
	// MX 1/2005
	double PorosityEffectiveConstrainedSwellingConstantIonicStrength(long, double, double, double*);
	// Permeability
	// Permeabilty stress corrector WW
	double* c_coefficient;
	unsigned geo_dimension;
	unsigned geo_inclination; // inclination of domain/sub domain from horizontal anticloclwise: AKS
	int permeability_stress_mode;
	//
	// CMCD 9/2004 GeoSys 4
	double PermeabilityPressureFunctionMethod1(long, double);
	// CMCD 9/2004 GeoSys 4
	double PermeabilityPressureFunctionMethod2(long, double);
	// CMCD 9/2004 GeoSys 4
	double PermeabilityPressureFunctionMethod3(long, double);
	// CMCD 9/2004 GeoSys 4
	double PermeabilityPressureFunctionMethod4(long, double, double);
	friend class CMediumPropertiesGroup;

public:
	//-------------------------------------------
	// Methods
	CMediumProperties(void); // constructor
	~CMediumProperties(void); // destructor
	CMediumProperties* Get(std::string);
	CMediumProperties* GetDB(std::string);
	CMediumProperties* GetByGroupNumber(int);
	void Set(std::string, std::string, double);
	void SetDB(std::string, std::string, double);
	int GetPropertyType(std::string);
	std::ios::pos_type Read(std::ifstream*);
	void Write(std::fstream*);
	void WriteTecplot(std::string);
	double* PermeabilityTensor(long index);
	// CMCD 9/2004 GeoSys 4
	double Porosity(FiniteElement::CElement* assem = NULL);
	// CMCD 9/2004 GeoSys 4
	double TortuosityFunction(long number, double* gp, double theta, CFiniteElementStd* assem = NULL);
	// CMCD 9/2004 GeoSys 4
	double NonlinearFlowFunction(long number, int gp, double theta, CFiniteElementStd* assem = NULL);
	// CMCD 9/2004 GeoSys 4
	double PermeabilityPressureFunction(long index, double* gp, double theta);
	// CMCD 9/2004 GeoSys 4
	double PermeabilityPorosityFunction(long index, double* gp, double theta);
	double PermeabilityFunctionPressure(long index, double PG2); // WX:  05.2010
	double PermeabilityFunctionStrain(long index, int nnodes, CFiniteElementStd* h_fem); // WX
	double PermeabilityFunctionEffStress(long index, int nnodes, CFiniteElementStd* h_fem); // AS:08.2012
	double StorageFunctionEffStress(long index, int nnodes, CFiniteElementStd* h_fem); // AS 08.2012
	double PorosityVolStrain(long index, double val0, CFiniteElementStd* assem); // WX: 03.2011
	double KozenyCarman(double k0 /*old permeability*/,
	                    double n0 /*old porosity*/,
	                    double n /*new porosity*/); // HS: 11.2008
	double KozenyCarman_normalized(double k0 /*old permeability*/,
	                               double n0 /*old porosity*/,
	                               double n /*new porosity*/); // HS: 11.2008
	// WW
	double KozenyCarmanNew(double k_init, double n_init, double n_t); // AB
	double VermaPruess(double k_init, double n_init, double n_t); // AB
	void CalStressPermeabilityFactor(double* kfac, const double T = 273.0);
	// WW
	void CalStressPermeabilityFactor2(double* kfac, const double T = 273.0);
	// WW
	void CalStressPermeabilityFactor3(double* kfac);
	void CalStressPermeabilityFactor3_Coef(); // WW
	void CalStressPermeabilityFactor4(double* kfac, double);
	// CMCD 9/2004 GeoSys 4
	double StorageFunction(long number, double* gp, double theta);
	double HeatCapacity(long number, double theta, CFiniteElementStd* assem = NULL);
    double CalcIceVolFrac(double T_in_dC, double freezing_sigmoid_coeff, double porosity); // TYZ
	double Calcsigmoidderive(double freezing_sigmoid_coeff, double phi_i, double porosity); // TYZ
	double* HeatConductivityTensor(int number); // MX
	double* HeatDispersionTensorNew(int ip); // CMCD
	double* MassDispersionTensorNew(int ip, int phase); // CMCD, SB, BG
	double* DispersionTensorMCF(int ip, int PCSIndex, int CIndex, double* variables = NULL);
	// OK
	double Density(long number, double* gp, double theta);
	// Capillary pressure functions
	double CapillaryPressureFunction(const double wetting_saturation);
	double PressureSaturationDependency(double wetting_saturation, bool invert);
	// JT: No longer used // double SaturationPressureDependency(const double capillary_pressure, bool allow_zero =
	// false);
	double SaturationCapillaryPressureFunction(const double capillary_pressure);
	// WW
	double PermeabilitySaturationFunction(const double wetting_saturation, int phase);
	double GetEffectiveSaturationForPerm(const double wetting_saturation, int phase); // JT
	// MX 1/2005
	double PorosityVolumetricChemicalReaction(long);
	double Porosity(long number, double theta); // CMCD 9/2004 GeoSys 4
	// CMCD 4/2005
	void SetConstantELEarea(double area, int group);
	// OK
	void SetDistributedELEProperties(std::string);

	void WriteTecplotDistributedProperties(); // OK
	double HeatTransferCoefficient(long number, double theta, CFiniteElementStd* assem); // NW
	double ParticleDiameter();
	unsigned GetGeoDimension(void) { return geo_dimension; }
	/**
	 * the type of the geometric entity, the material property is assigned to
	 * @return a value of the enum GEOLIB::GEOTYPE
	 */
	GEOLIB::GEOTYPE getGeoType() const { return _geo_type; }
	CFEMesh* getMesh(void) { return _mesh; }
	// Properties
private:
	// PCS
	std::string pcs_type_name; // YD
public:
	CRFProcess* m_pcs; // OK
	std::vector<std::string> pcs_name_vector;

private:
	std::vector<std::string> porosity_pcs_name_vector;
	CFEMesh* _mesh; // OK

	/**
	 * attribute describes the type of the geometric entity the
	 * material property is assigned to
	 */
	GEOLIB::GEOTYPE _geo_type;
	FiniteElement::FrictionPhase _fric_phase;

public:
	// GEO
	std::string geo_name;
	std::vector<std::string> geo_name_vector; // OK
	double geo_area;
	std::string geo_area_file; // OK

	double density;
	std::string name;
	int number;
	int porosity_model; // porosity
	int porosity_curve;
	double porosity_model_values[15];
	double porosity;
	double KC_porosity_initial; // HS 11.2008
	double KC_permeability_initial; // HS 11.2008
	std::string porosity_file; // OK/MB
	int tortuosity_model;
	double tortuosity_model_values[10];
	double tortuosity;
	int flowlinearity_model;
	double flowlinearity_model_values[10];
	int storage_model; // storativity
	double storage_model_values[10];
	double storage;
	int conductivity_model;
	double conductivity;
	int unconfined_flow_group;
	int permeability_model; // permeability
	double permeability;
	double local_permeability; // CB
	double permeability_tensor[9];
	double permeability_porosity_updating_values[5]; // ABM: Maximum of 2 values in Verma-Pruess case
	std::string permeability_tensor_type_name;
	std::string permeability_porosity_updating_type_name; // ABM
	std::string tortuosity_tensor_type_name;
	int permeability_tensor_type;
	int permeability_porosity_updating_type; // ABM
	int tortuosity_tensor_type;

	double ElementVolumeMultiplyer; // Multiplyer of element volume

	std::string
	    PhaseHeatedByFriction; // In TNEQ/TES models: dissipated heat due to friction into solid or fluid energy balance

	int permeability_pressure_model;
	double permeability_pressure_model_values[10];
	double permeability_pressure_rel;
	int permeability_strain_model; // WX: permeability function strain model. 05.2010
	int permeability_strain_model_value[6]; // WX:permeability fuction strain model value. 05.2010
	int permeability_effstress_model_value[3]; // AS:perlmeability function eff stress 08.2012
	int permeability_effstress_model;
	int storage_effstress_model_value[3]; // AS:storage function eff stress 08.2012
	int storage_effstress_model;
	//
	// Relative permeability (JT)
	int permeability_saturation_model[MAX_FLUID_PHASES];
	double minimum_relative_permeability;
	double residual_saturation[MAX_FLUID_PHASES];
	double maximum_saturation[MAX_FLUID_PHASES];
	double saturation_exponent[MAX_FLUID_PHASES];
	double perm_saturation_value[MAX_FLUID_PHASES];
	//
	std::string permeability_file; // SB //OK/MB string permeability_dis_type_file;
	std::string tortuosity_file; // PCH
	bool entry_pressure_conversion; // JT
	int capillary_pressure_model;
	int permeability_porosity_model;
	double permeability_porosity_model_values[10];
	double storativity;
	double capillary_pressure_values[5]; // JT2012
	double heat_capacity; // thermal properties
	int mass_dispersion_model;
	double mass_dispersion_longitudinal;
	double mass_dispersion_transverse;
	double lgpn; // local grid peclet number
	int heat_dispersion_model;
	double heat_dispersion_longitudinal;
	double heat_dispersion_transverse;
	double heat_conductivity_tensor[9];
	int fct_number; // functions
	int heat_diffusion_model;
	double base_heat_diffusion_coefficient;
	int evaporation; // if it is 647 then evaporation ON, else OFF: and corresponding heat loss will compensated by heat
	// ST
	double heatflux;
	double vaporfraction;
	// aux
	int m_color[3];
	bool selected;
	int mode;
	// surface water
	// JOD
	double friction_coefficient, friction_exp_slope, friction_exp_depth;
	// JOD
	double overland_width, rill_height, rill_epsilon;
	bool channel;
	double argument; // OK
	// Dual Richards transfer coefficient  YD
	double transfer_coefficient;
	// double unsaturated_hydraulic_conductivity;
	double specific_storage;
	int vol_mat_model; // CB
	double vol_mat; // SB
	int vol_bio_model; // CB
	double vol_bio; // SB
	double foc; // organic carbon content
	double forchheimer_cf; // NW geometric constant for Forchheimer equation
	double forchheimer_De; // NW equivalent diameter of the bed
	double forchheimer_a1; // NW
	double forchheimer_a2; // NW

	int heat_transfer_model; // NW
	int effective_heat_transfer_model; // NW
	double heat_transfer_model_value; // NW

	int particle_diameter_model;
	double particle_diameter_model_value;
	// Get friction phase - TN
	FiniteElement::FrictionPhase getFrictionPhase() const;

	// Set value for friction phase - TN
	void setFrictionPhase(FiniteElement::FrictionPhase fric_phase);

	// CB Chiogna et al alpha-t model
	int alpha_t_model;
	double graindiameter;
	double hydraulicrad;
	double betaexpo;
    // HS: borehole heat exchanger related parameters
    bool is_BHE;  // indicate whether this MMP is a borehole heat exchanger. 
    BHE::BHE_TYPE bhe_type;
    BHE::BHE_BOUNDARY_TYPE bhe_bound_type;
    BHE::BHE_DISCHARGE_TYPE bhe_2u_discharge_type;
    bool bhe_use_ext_therm_resis; 
	bool bhe_user_defined_therm_resis;
	bool bhe_use_flowrate_curve;
    double bhe_power_in_watt_val; 
    double bhe_delta_T_val; 
    double bhe_length, bhe_diameter, bhe_refrigerant_flow_rate, bhe_inner_radius_pipe;
    double bhe_outer_radius_pipe, bhe_pipe_in_wall_thickness, bhe_pipe_out_wall_thickness;
    double bhe_therm_conductivity_pipe_wall, bhe_therm_conductivity_grout, bhe_pipe_distance;
    std::size_t bhe_fluid_type_idx;
    std::size_t bhe_power_in_watt_curve_idx; 
    int bhe_heating_cop_curve_idx; // Heating COP curve index
	int bhe_cooling_cop_curve_idx; // Cooling COP curve index
	int bhe_flowrate_curve_idx;
    double bhe_refrigerant_viscosity;
    double bhe_refrigerant_density;
    double bhe_refrigerant_heat_capacity;
	double bhe_refrigerant_alpha_L; 
    double bhe_grout_density;
	double bhe_grout_porosity;
    double bhe_grout_heat_capacity;
    double bhe_regrigerant_heat_conductivity;
    double bhe_therm_resistance; 
    double bhe_intern_resistance; 
	double bhe_R_fig, bhe_R_fog, bhe_R_gg1, bhe_R_gg2, bhe_R_gs;
	double bhe_switch_off_threshold;
	// BHE Net heat pump parameters
	bool is_heat_pump;
	std::string heat_pump_name;
	BHE::HEAT_PUMP_BOUNDARY_TYPE heat_pump_boundary_type;
	double heat_pump_power_val;
	double heat_pump_delta_T_val;
	double heat_pump_flowrate;
	int heat_pump_power_curve_idx;
	int heat_pump_COP_curve_idx;
	int heat_pump_fluid_idx;
	double heat_pump_refrigerant_density;
	double heat_pump_refrigerant_heat_capacity;
	// BHE Net distributor parameters
	bool is_distributor;
	std::string distributor_name;
	int distributor_n_in;
	int distributor_n_out;
	Eigen::VectorXd distributor_in_ratios;
	Eigen::VectorXd distributor_out_ratios;
	// BHE Net pipe parameters
	bool is_pipe;
	std::string pipe_name, pipe_from, pipe_to;
	int pipe_from_port, pipe_to_port;
	double pipe_loss_coeff;
};

class CMediumPropertiesGroup // YD
{
public:
	CMediumPropertiesGroup() : m_msh(NULL), OrigSize(0) {}
	void Set(CRFProcess* m_pcs);
	std::string pcs_name;
	std::string pcs_type_name;
	CFEMesh* m_msh;
	std::vector<CMediumProperties*> mmp_group_vector;

private:
	int OrigSize; // For excavation simulation.
};

// YD
extern CMediumPropertiesGroup* MMPGetGroup(const std::string& pcs_type_name);
// YD
extern std::list<CMediumPropertiesGroup*> mmp_group_list;
extern void MMPGroupDelete(/*string pcs_type_name*/);

extern std::vector<CMediumProperties*> mmp_vector;
extern void MATLoadDB(std::string);
extern std::list<std::string> keywd_list; // keyword-referenzliste "kw"
extern void comp_keywd_list(std::string);
extern void read_keywd_list(void);
extern std::list<std::string> mat_name_list;

extern void MMPWrite(std::string);
extern bool MMPRead(std::string);
extern void MMPDelete();
extern CMediumProperties* MMPGet(const std::string&);
extern void MMP2PCSRelation(CRFProcess*);
extern void GetHeterogeneousFields(); // SB
extern long GetNearestHetVal2(long EleIndex,
                              CFEMesh* m_msh,
                              std::vector<double>
                                  xvals,
                              std::vector<double>
                                  yvals,
                              std::vector<double>
                                  zvals,
                              std::vector<double>
                                  mmpvals);
double GetAverageHetVal2(long EleIndex,
                         CFEMesh* m_msh,
                         std::vector<double>
                             xvals,
                         std::vector<double>
                             yvals,
                         std::vector<double>
                             zvals,
                         std::vector<double>
                             mmpvals);
extern bool MMPExist(std::ifstream* mmp_file); // OK
extern bool MMPExist(); // OK

#define MMP_FILE_EXTENSION ".mmp"

#endif
