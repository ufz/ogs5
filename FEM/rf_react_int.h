/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

 /* reaction interface header file */


//CB: What's the use of this?
//#ifndef reaction_INC
//#define reactions_INC
#include <string>
#include <vector>
//Water moles per kg of water
#define MOLH2OPERKG 55.50843506
#define REI_FILE_EXTENSION ".rei"

/* Structure for exchange of reaction related data */
class REACTINT{
  private:
	public:
    REACTINT(void);
	~REACTINT(void);
// Data
    int flowtype;
    bool unitconversion ;
    bool constantdensity ;
    bool constanttemperature;
    bool constantpressure;
	bool residual_back;
    double c_TT, c_PP;  // CB check for variable type c_PP
    int ssp_outstep;
    int pcs_outstep;
    bool dump_min_moles;
    bool dump_all_pcs ;
    bool dump_mass_integrals;
    bool s_water_limit ;
    double WaterSatLimit;
    std::string WaterSpeciesName;
    std::string NeutralCO2name;
    std::string SodiumSpeciesName;
    bool icOutput;
    bool icSolidUpdate ;
    bool readNodePoro;
    std::string porofile;
    bool vle_flag, vle_p_flag ;
    bool pcs_rename_init_flag ;
    bool pcs_rename_pre_flag ;
    bool pcs_rename_post_flag ;
    bool poroupdate_flag;
    bool heatpump_2DhTO2Dv;
    double heatpump_Z;
	int t_step;

    std::vector <double> Temp_store;
    std::vector <long> Temp_GHP_mapidx;

    std::vector <double> node_porosity;
    std::vector <double> node_ini_porosity;
    std::vector <double> water_conc;
    std::vector <std::string> water_conc_species;
    std::vector <std::vector <int> > ElementInSpecies;
    //std::vector <int> sp_pcsind; // store pcs_idx of all species
	std::vector <int> sp_varind; // store node_value_idx+timelevel of all species
    std::vector <bool> dried_out_nodes; // for eclipse coupling

	int nodenumber;
	std::vector<std::string>               pcs_rename0_init,      pcs_rename0_pre,       pcs_rename0_post;
	std::vector<int>                       pcs_rename0_idx_init,  pcs_rename0_idx_pre,   pcs_rename0_idx_post;
	std::vector<std::vector<std::string> > pcs_rename1_init,      pcs_rename1_pre,       pcs_rename1_post;
	std::vector<std::vector<int> >         pcs_rename1_idx_init,  pcs_rename1_idx_pre,   pcs_rename1_idx_post;
	std::vector<std::vector<double> >      pcs_rename1_stoi_init, pcs_rename1_stoi_pre,  pcs_rename1_stoi_post;
	std::vector<double>                    pow_stoi_init, pow_stoi_pre, pow_stoi_post;

	typedef struct
	{
		int    aq_idx,   vp_idx;        // check input species list and get the species idx
		double aq_value, vp_value; // get the value from processes
		std::string aq_name,  vp_name;   // from input file, name of process
		int idx_aq_species;
		double delta;
		double eq_ac;
	}VLE_type;
	std::vector<VLE_type> VLE_conditions, VLE_pressure;


// Member functions

    REACTINT* GetREACTINT(void);
    bool Read(std::ifstream*);
	void InitREACTINT(void);
    void CalcWaterConc(void);
    void CalcUnitConversionFactors(long index, double *fl, double *fs, bool molal);
    double GetCO2SolubilityDuan(long node);
    double GetPressure(long );
    double GetTemperature(long);
    double CalcDensityOfWater(void);
    std::vector<int> formula2index(std::string formula);
    static std::vector<std::string> string2vector(std::string line);	//split string line to pieces, and store in a vector
    //std::vector<std::string> string2vector(std::string line);
    double GetWaterSaturation( long );
    void ReactionPreProcessing(void);
    void ReactionPostProcessing(bool);
    void SetInitialPorosityAndPermToElementValue(void);
    void PorosityVolumetricReactionUpdate(void);
    void PorosityVolumetricReactionUpdate_2(void);
    void PermeabilityPorosityUpdate(void);
    void CopySymmetricConcentrationsInRadialModel(void);
    void DumpSolidSpeciesMoles(void);
    void DumpAllVariables(void);
    void DumpMassIntegrals(void);
    void CheckForDriedOutNodes(void);
    void CopyAllConcentrationsToOtherTimeLevel(bool);
    void ReadRestartNodePoro(long);
    double LiquidViscosity_Yaws_1976(double );
    double LiquidDensity_Busch(double);
    void Heatpump_2DhTO2Dv_Mapping(bool);


};
extern std::vector <REACTINT*> REACTINT_vec;
extern double GetNodePhaseVolume(long , double , int );
extern bool REACINTRead(std::string);

// CB moved here from rf_pcs.h file
typedef struct // CB DL CO2 phase transition
{
	std::string name;                     //fluid name
	double temperature;
	double pressure;
	double density;                       //density g/cm^3
	double viscosity;                     //viscosity mPa.s
	double volume;                        //volume cm^3
	double mass;                          //weight g
	double CO2;                           //mole of CO2
	double H2O;                           //mole of H2O
	double NaCl;                          //mole of NaCl
	double C_GAS;	  //mole of CO2 in Gas
	double H2;        //mole of H2
}Phase_Properties2;

extern void VLE_CalcNewPressure(double T, double &P, double &V_gas, double &V_liquid, double CO2_gas_old, double CO2_gas_new, double CO2_liquid_old, double CO2_liquid_new, double H2O_liquid_old, double H2O_liquid_new, double &rho);
extern void VLE_isobaric_fixphase(double T, double P, Phase_Properties2 &vapor, Phase_Properties2 &liquid, Phase_Properties2 &solid, int f);
extern void VLE_isochoric_fixphase(double T, double &P, Phase_Properties2 &vapor, Phase_Properties2 &liquid, Phase_Properties2 &solid, int f);

extern void VLE_isobaric_transphase_H2(double T, double P, Phase_Properties2 &vapor, Phase_Properties2 &liquid);
extern void VLE_isochoric_transphase_H2(double T, double P, Phase_Properties2 &vapor, Phase_Properties2 &liquid);

extern void VLE_isobaric_fixphase_H2(double T, double P, Phase_Properties2 &vapor, Phase_Properties2 &liquid);
extern void VLE_isochoric_fixphase_H2(double T, double P, Phase_Properties2 &vapor, Phase_Properties2 &liquid);




//#endif
