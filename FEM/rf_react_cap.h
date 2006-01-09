/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

 /* reaction package header file */

//#ifndef reaction_INC
//#define reactions_INC

#include <vector>
//using namespace std;

//CB2406 #ifdef OGS_FEM_CAP // CAP_REACT

/* Structure for exchange of reaction rates */
class REACT_CAP{
	private:
	public:
		REACT_CAP(void); // Initialize
		~REACT_CAP(void); // Clean up

// Data
	int *rateflag;     /* flag used for determining if reaction are calculated */
	long nodenumber;   /* number of nodes, on which reactions are calculated */
	bool flag_cap;     /* flag if *.cap file exists   DL 28,10,08*/
    std::string data_file;
	std::string data_format;
	int    mass_num;
	std::string mass_type;
	std::string species_define;
    std::vector<int>      pcs_mass_idx;      //index of mass_trans order-number in processes list, (some other process are not mass transport)
	std::vector<int>      pcs_ospecies_idx;  //index of output species order-number in processes list     21.01.2009
	std::vector<int>      pcs_nspecies_idx;  //index of nlog-output species order-number in processes list 02.02.2009
	std::vector<int>          nspecies_idx;  //index of nlog-output species order-number in NLOG_AC list

	std::vector<int>      phase_list;
	std::vector<char*>    phase_post;

	std::vector<char*>    species_list_name; //From .cap list section, species name
	std::vector<int>      species_list_phase;//From .cap list section, phase flag for species 0-solid, 1-gas, 2-aqueous
	std::vector<long int> species_list_idx;  //From .cap list section, index of species in ChemApp DATA file

	std::vector<char*>    species_name;   // From Process, list of species name for chemical reaction
	std::vector<int>      species_phase;  // From Process, phase flag for species 0-solid, 1-gas, 2-aqueous
	std::vector<long int> species_idx;    // From Process, index of species in ChemApp DATA file
	std::vector<int>      species_mobil;  // From Process, process mobile
	std::vector<int>      species_dormant;// From KINETIC_SPECIES, 0-equilibrium 1-kinetic

	std::vector<char*>    species_kin_name; //DL 03.2010
	std::vector<int>      species_kin_phase;//
	std::vector<long int> species_kin_idx;  //

	std::vector<char*>    species_relative_activity_name; //DL 03.2013
	std::vector<int>      species_relative_activity_idx;
	std::vector<long int> species_relative_activity_idx_pcs;
	std::vector<double>   species_relative_activity, species_relative_activity_ia,species_relative_activity_calc;
	std::vector<int>      species_relative_activity_state;

	std::vector<int>      Node_All_SI;
	bool show_data_file;

	std::vector<std::string>          pcs_rename0, pcs_rename0_pre;
	std::vector<std::vector<std::string> > pcs_rename1, pcs_rename1_pre;
	std::vector<std::vector<double> >      pcs_rename_stoi, pcs_rename_stoi_pre;
	std::vector<int>                  pcs_rename_idx0, pcs_rename_idx0_pre;
	std::vector<std::vector<int> >         pcs_rename_idx1, pcs_rename_idx1_pre;

	int current_node; // for Kin Reaction calc
	double current_TT, current_PP; // for Kin Reaction calc
	typedef struct
	{
		std::string  type;
		double  logK;
		double  rate;
		std::vector<char*>    species_name;
		std::vector<int>      species_phase;
		std::vector<long int> species_idx;
		std::vector<double>   species_stoi;
		//for Kinetic reaction
		bool   is_Kin;
		int    logK_type;
		double logK_kin;
		double theta;
		double eta;
		double SA;
		double PF;
        std::vector<double>   EA;
		std::vector<double>   K25;
		std::vector<std::string>   ref_species;
		std::vector<int>      ref_idx; // position in mass trans list --> node_ac matrix
		std::vector<double>   ref_expo;
	}reaction;
	reaction React;
	std::vector<reaction> Kin_Reactions, Kin_HKF_Reactions;
	std::vector<int>  Kin_HKF_index; //look up table between Reactions and HKF_Species

	typedef struct
	{
		//input
		std::string name0; //name in Helgeson datafile
		//from datafile
		std::string name;
		std::string abbrev; //abbreviation
		std::string scform; //structural chemical formula
		std::string ecform; //elemental chemical formula
		std::string ref[2];
		double charge;  //charge
		int    type  ;  //phase flag
						//0-minerals that do not undergo phase transitions
						//1-minerals that undergo one phase transition
						//2-minerals that undergo two phase transitions
						//3-minerals that undergo three phase transitions
						//4-gases
						//5-aqueous species
		double param[9][4]; //parameters
		//calculating result
		double S; //entropy
		double H; //enthalpy
		double G; //Gibbs free energy
		double V; //volume
		double Cp;//heat capacity
	}SpeciesData;


	std::vector<std::vector<double> > node_logK; // node, reaction
	std::vector<std::vector<double> > node_ac;   // node, spec
	std::vector<SpeciesData> Kin_HKF_species;
	std::vector<std::vector<double> > node_HKF_logK;


	std::vector<char*>            nlog_name; //correspond name with pcs and nlog
	std::vector<char*>    species_nlog_name; //use .cap to define -log(ac) value
	std::vector<int>      species_nlog_phase;
	std::vector<long int> species_nlog_idx;

	int pcs_redox;
	std::vector<char*>    species_redox_name;
	std::vector<int>      species_redox_phase;
	std::vector<long int> species_redox_idx;
	std::vector<long>             redox_stoi;
	bool check_no_reaction_nodes; /* flag if CheckNoReactionNodes has been performed */
//	double temperature; /* temperature of water at node as specified in input file */
//	int heatflag;      /* if 1, heat transport is active and passed to phreeqc */
	std::vector <std::string> ic_2_bc_geometry_type;  // Convert eq concentrations at this gemetry to BC there
	std::vector <std::string> ic_2_bc_geometry_name;
  std::vector<size_t> ic_2_bc_GeoID;

	std::vector <std::string> ic_2_bc_species;
	std::vector <double> ic_2_bc_species_value;
  /**
   * @param geo_obj object of class GEOObjects managing the geometric entities
   * @param unique_name unique name to access the geometric entities in geo_obj
   */
	void ConvertIC2BC(const GEOLIB::GEOObjects& geo_obj, const std::string& unique_name);

	std::ofstream warning_out;
// Member functions

/**
 * @param geo_obj object of class GEOObjects managing the geometric entities
 * @param unique_name unique name to access the geometric entities in geo_obj
 */
	std::ios::pos_type Read(std::ifstream*, const GEOLIB::GEOObjects& geo_obj, const std::string& unique_name);

	void CreateREACT(void);
	void InitChemApp(void);
	void ResetElement(void);
	void DefSpeciesInChemApp(void);
	void DefSpeciesInChemAppList(void);
	void DefSpeciesInChemAppAll(void);

	void UndefOutSpecies(void);

	void RecoverChemSystem(void);

	bool CmpSpeciesList(void);
	bool CmpSpeciesListDAT(void);
	bool CmpSpeciesListPCS(void);
	bool CmpSpeciesListNLOG(void);
	bool CmpSpeciesListREDOX(void);
	bool CmpSpeciesListKinReact(void); //DL 03.2010
	bool CmpSpeciesListKIN(void); //DL 03.2010

	bool CmpSpeciesListRelativeActivity(void); // DL 04.2013

	bool SetKinHKFspecies(void);  //DL 11.2010
	bool LoadHKFparam(void);
	bool HKFcalc(double T, double P);

	void LoadMassPCS(void);

	void CompNamePhase(std::string, char*&, int&);

	int  CheckNoReactionNodes(void);
	void LoopNodeReact(int, int);
		std::vector<double> KineticReact(std::vector<double>); //DL 05.10.10
		void KinInit(void);
		void KinRate(void);
		std::vector<double> derivs(std::vector<double>); //DL 05.10.10
		std::vector<double> ODE(std::vector<double>, std::vector<double>, int, double);
		void KinParamUpdata(int, int, std::vector <double> & spvc);
		void KinParamUpdataHKF(int,double T, double P, int err);
	void LoopNodeReact_Liquid_Vapor(int, int);
	void LoopNodeReact_Liquid_Solid(int, int);

	void ExecuteReactionsChemApp(int, int); //DL 28,10,08

	void ShowDataFile(void);
	void SetAllSolidAsDormant(void);

};
extern std::vector <REACT_CAP*> REACT_CAP_vec;
extern bool REACT_CAP_Read(std::string, const GEOLIB::GEOObjects& geo_obj, const std::string& unique_name);

#ifdef OGS_FEM_CAP // CAP_REACT
#endif
