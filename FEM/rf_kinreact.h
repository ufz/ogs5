/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   rf_kinreact.cpp

                                    KINETIC REACTIONS

   FEMLib-Object KinReact

   Programming:
   01/2004    Dirk Sch�fer       Original IMplementation
   02/2006    Sebastian Bauer    Adaption to C++ Class structure, new FEM concept

 ***************************************************************************/
#ifndef rf_kinreact_INC
#define rf_kinreact_INC
// C++ STL
#include <fstream>
#include <string>
#include <vector>
#include <time.h>
// GEOLIB
#include "GEOObjects.h"

/* residual for linearisation of critical functions */
#define maxMonod 5
#define maxInhibition 5
#define maxBioReactions 30
#define maxNumber_of_Components 30

#define KRC_FILE_EXTENSION ".krc"

/* New class KinReaction: contains the kinetic reactions and all necessary data structures for them */
// C++ Class Monodsubstruct
class MonodSubstruct
{
private:
public:
	std::string species; // Name of species
	int speciesnumber; // number of species;
	double concentration; // Monod concentration
	double order; // Order of monod term
	int isotopecouplenumber; // CB isotope fractionation : specis number of isotope partner
	// CB for Threshhold terms
	bool threshhold;
	double threshConc;
	double threshOrder;
	MonodSubstruct(void);
	~MonodSubstruct(void);
};

// CB _drmc_

/* New class KinReaction: contains the kinetic reactions and all necessary data structures for them */
// C++ Class Monodsubstruct
class MicrobeData
{
private:
public:
	std::string MicrobeName;
	std::string MonodReactionName;
	std::vector<double> Gibbs;
	std::vector<double> _drmc_level;
	double steepness;
	double k_incr;
	double k_decr;
	double G0;
	double dt;
	double decayrate;
	bool _drmc_;
	int MonodReactionIdx;

	double GetGibbsEnergy(long node);
	bool Read(std::ifstream* in);

	MicrobeData(void);
	~MicrobeData(void);
};

// C++ Class MinKinMech
class MinkinMech
{
private:
public:
	int no_mechSpec;
	std::vector<int> mechSpeciesIdx;
	std::vector<std::string> mechSpeciesNames;
	std::vector<double> mechSpeciesExpo;
	double Eact;
	double k25;

	MinkinMech(void);
	~MinkinMech(void);
	double Mech(double* c, int node);
};

// C++ Class CKinReact
class CKinReact
{
private:
	/**
	 * name of reaction
	 */
	std::string name;
	/**
	 * type of reaction: monod, exchange, NAPLdissolution, ...
	 */
	std::string type;

	friend bool KRRead(std::string const& file_base_name, GEOLIB::GEOObjects const& geo_obj,
	                   std::string const& unique_name);
	friend void KRConfig(const GEOLIB::GEOObjects& geo_obj, const std::string& unique_name);

public:
	CKinReact(void); // Constructor
	~CKinReact(void); // Destructor

	std::string const& getType() const { return type; }
	std::string const& getName() const { return name; }
	int number; /* counter */
	int number_reactionpartner; /* Number of chemical species involved in reaction */
	std::vector<std::string> reactionpartner; /* all names of reaction partners stored here */
	std::vector<double> stochmet; /* stochiometric coefficients for each reactionpartner stored here */
	double rateconstant; /* rateconstant */
	double rateorder; /* order of reaction */
	int number_monod; /* Number of Monod terms */
	int number_inhibit; /* Number of inhibition terms */
	int number_production; /* number of production terms */
	int number_isotope_couples; /* number of production terms */
	std::vector<MonodSubstruct*> monod; /* saves monod concentrations and names of species */
	std::vector<MonodSubstruct*> inhibit; /* saves inhibit concentrations and names of species */
	std::vector<MonodSubstruct*> production; /* saves production concentrations, orders and names of species */
	int grow; /* growth or no growth */
	std::string bacteria_name;
	int bacteria_number;
	std::vector<double> ProductionStoch; // stochiometry of reaction
	//    vector <double>	ProductionStoch2; // stochiometry of reaction - short version
	std::vector<MonodSubstruct*> ProdStochhelp; // store input values
	double specif_cap;
	// CB Isotope fractionation
	std::string Isotope_light;
	std::string Isotope_heavy;
	std::string degType;
	double isoenfac;
	bool T_dependence;
	int T_model;
	std::vector<double> T_params;

	// CB _drmc_
	bool _drmc_;
	std::string Microbename;
	std::string DormType;
	int DormTypeIdx;
	int MicrobeData_idx;
	double dG0; // standard Gibbs energy of reaction
	double Yieldcoefficient;

	// CB Not this particular reaction on specified GEO-Objects; Data structures
	std::vector<std::string> NotThisReactGeoName;
	std::vector<size_t> NotThisReactGeoID; // 06/2010 TF
	std::vector<std::string> NotThisReactGeoType;
	std::vector<bool> switched_off_node;
	long currentnode; // CB 19/10/09 This is eclusively for Brand model to allow porosity in Inhibition constant
	// calculation

	// exchange data
	std::vector<std::string> ex_species_names;
	std::vector<int> ex_species;
	std::vector<double> ex_param;
	int exSurfaceID;
	std::string exType; /* sorption type: linear, langmuir, exchange */

	//#ds NAPLdissolution data
	std::string blob_name; /* name of blob-class */
	int blob_ID; /* id number of blobs where the NAPL phase resides */
	double Csat_pure; /* maximum solubility of the pure NAPL phase */
	double current_Csat; /* current solubility after considering Roult's law, interally calculated */
	double Density_NAPL; /* density of the pure NAPL phase */
	//	double  ConversionFactor;       /* factor to convert concentrations to mol/kg */
	// SB speed-up flags
	int typeflag_monod; /* set to 1 if reaction is monod type */
	int typeflag_exchange; /* set to 1 if reaction is exchange type */
	int typeflag_exchange_linear; /* set to 1 if reaction is linear exchange type */
	int typeflag_exchange_langmuir; /* set to 1 if reaction is langmuir exchange type */
	int typeflag_exchange_freundlich; /* set to 1 if reaction is freundlich exchange type */
	int typeflag_napldissolution; /* set to 1 if reaction is NAPL dissolution */
	int typeflag_iso_fract; /* set to 1 if reaction is isotope fractionation */
	int typeflag_mineralkinetics; /* set to 1 if reaction is mineral kinetics */
	int typeflag_gasdissolution;

	// CB new parallel computation data structure, replaces double current_Csat;
	std::vector<double> Current_Csat; /* initial interfacial area, interfacial area of last iteration */
	// Mineral kinetics data
	std::vector<MinkinMech*> mechvec; // stores individual mechanisms of this reaction
	std::vector<int> minSpeciesIdx; // stores cp vector index of mineral species
	std::vector<double>
	    Am; // initial specific reactive mineral surface area m²/m³ aquifer = Am * rho_b = m²/kg * kg/m³aquifer
	std::vector<double> Km; // Equilibrium constant
	std::vector<double> Cminini; // initial mineral concentration
	int number_Mech; // Number of Mechanism terms */
	bool Km_uniform;
	bool Km_CHEMAPP;
	bool Km_HKF;
	bool scale_rate;
	bool lagneau;

	double Km_default;
	bool Am_constant;
	int Am_model;
	double Am_ini; // initial mineral Surface area
	double Theta;
	double Eta;
	bool precip_baseterm_only;
	double precipfactor;
	double precipexponent;
	std::string mineral_name;
	std::string chemapp_name;
	int mineral_number; // cp_vec idx of mineral
	int water_number; // cp_vec idx of mineral
	double Temperature;
	double Pressure;
	double OmegaThreshhold;
	double MinKinRateCoeff;
	bool ExplicitMinKinRateCoeff;

	/* Methods */
	/**
	 * read data from stream
	 * @param in input stream from file
	 * @param geo_obj object of class GEOObjects that manages the geometric entities
	 * @param unique_name the name of the project to access the right geometric entities
	 * @return true (in every case) ToDo
	 */
	bool Read(std::ifstream* in, const GEOLIB::GEOObjects& geo_obj, const std::string& unique_name);

	void Write(std::ofstream*); /* Class Write Function */
	void ReadReactionEquation(std::string); /* Read function for chemical equations */
	int CheckReactionDataConsistency(void); /* check data set */
	void TestWrite(void); // test output function

	// CB isotope fractionation + higher order terms
	double Monod(double, double, double, double);
	double Inhibition(double, double);
	double BacteriaGrowth(int r, double* c, double sumX, int exclude, long node);
	int GetPhase(int);
	//   double GetPorosity( int comp, long index );
	// CB replaced by
	double GetReferenceVolume(int comp, long index);
	double GetDensity(int comp, long index);
	double GetNodePoreVelocity(long node);
	double GetNodeDarcyVelocity(long node);
	double GetPhaseVolumeAtNode(long node, double theta, int phase);
	double GetSaturation(long node, int theta, int phase);
	double GetMaxSolubility(long node, double density);
	double EffectiveYield(long node);
	// MinKin
	double Omega(double* c, int node);
	// void IonicStrength(int nnode, int ncomp);
	// void ActvityCoefficients(int nnode, int ncomp);
	// double Arrhenius(double, double, double);
	double MinRate(double* c, int node, double dt, double Om);
	double MinRateConstant(double* c, int node, double disso);
};

//#ds Class for blob properties
class CKinBlob
{
public:
	CKinBlob();
	~CKinBlob();
	/**
	 * read values from stream to initialize the CKinBlob object
	 * @param in input file stream
	 * @param geo_obj object of class GEOObjects that manages the geometric entities
	 * @param unique_name the name of the project to access the right geometric entities
	 * @return true (always) TODO
	 */
	bool Read(std::ifstream* in, const GEOLIB::GEOObjects& geo_obj, const std::string& unique_name);
	void Write(std::ostream& out) const;
	const std::vector<size_t>& getBlobGeoID() { return BlobGeoID; }
	std::string name;
	double d50;
	double Sh_factor;
	double Re_expo;
	double Sc_expo;
	double Geometry_expo;
	double Mass;
	double Volume;
	double Masstransfer_k;
	double current_Interfacial_area;
	std::vector<double> Masstransfer_K; /* CB: current Mass transfer coefficient, now saved in node vector*/
	std::vector<double> Area_Value;
	std::vector<double> Interfacial_area;
	std::vector<std::string> BlobGeoType;
	std::vector<std::string> BlobGeoName;
	// New-Sherwood-Number
	int shidx;
	double NContent_expo;
	double WContent_expo;
	double Pe_expo;
	double D50_expo;
	double Poro_expo;
	double Delta_expo;
	double UI_expo;
	double NContent_ini_expo;
	double NContent_res_expo;
	double GSR_expo;
	double Length_expo;
	double Tort_expo;
	double Pfannkuch_constant;
	double beta_4;
	double one_three, two_three, four_nine, fife_nine, elev_nine, fife_three;
	double grain_var;
	double grain_expo;
	double dM;
	double dS;
	double UI;
	double Ui_expo;
	double NCont_ini;
	double NCont_res;
	double NCont_Sh1;
	double NCont_Sh2;
	double GSR;
	double Tort;
	double Length;
	bool modSherwood;
	bool Sherwood_model;

	bool gas_dissolution_flag;

private:
	std::vector<size_t> BlobGeoID;
};

class CKinReactData
{
public:
	/* Data */
	int SolverType;
	double relErrorTolerance;
	double minTimestep;
	double initialTimestep;
	double usedt;
	int NumberReactions;
	int NumberLinear;
	int NumberLangmuir;
	int NumberFreundlich;
	int NumberMonod;
	int NumberNAPLdissolution;
	int NumberMineralkinetics;

	// CB _drmc_
	int NumberMicrobeData;

	clock_t cpu_time_krc;

	// biodeg data
	double maxBacteriaCapacity;
	std::vector<int> is_a_bacterium;
	double minBactConcentration;
	int minConcentrationMode;
	double minConcentrationThreshhold;
	double minConcentrationSet;
	//	vector <int> is_a_bacterium2; // short version
	// exchange data
	int maxSurfaces;
	std::vector<double> exSurface;
	// output flag
	bool testoutput;
	// index vector for shortening vectors c in kinetic calculations (omitting nonreacting species)
	//	vector <int> sp_index;
	//	int kr_active_species;
	// std::vector<int> sp_pcsind; // HS
	std::vector<int> sp_varind;

	// No reactions on specified GEO-Objects; Data structures
	std::vector<std::string> NoReactGeoName;
	std::vector<size_t> NoReactGeoID;
	std::vector<std::string> NoReactGeoType;
	std::vector<bool> is_a_CCBC;
	std::vector<std::string> AllowReactGeoName;
	std::vector<size_t> AllowReactGeoID;
	std::vector<std::string> AllowReactGeoType;

	// CB ReactDeact no reaction switch
	bool ReactDeactFlag; // method flag
	int ReactDeactPlotFlag; // flag for tecplot plots of flags each timestep
	double ReactDeactEpsilon; // treshhold
	double ReactDeactCThresh;
	std::vector<bool> ReactDeact; // flags for individual nodes
	std::vector<double> React_dCdT; // Sum of reaction rates for individual nodes
	// node indices of local neighborhood around individual nodes
	std::vector<std::vector<int> > ReactNeighborhood;
	int ReactDeactMode;
	bool ReactDeactRelative;

	bool copy_concentrations; // SB 01.2011 Flag for copying concentrations in radial models to save reaction simulation
	// time
	std::vector<long> copy_nodes; // SB 01.2011
	bool radial; // default case
	bool batch;
	bool TwoDinThreeD;

	bool scale_dcdt; // stabilization for ODE solver in derivs, sclae dcdt vector
	bool lagneau; // special configuration of minrate and porosity update for Lagneaubenchmark
	bool sortnodes;
	long* substeps;
	long* noderanks;
	double OmegaThresh;
	long noototnodes;
	long noocalcnodes;

	bool debugoutflag;
	std::string debugoutfilename;
	std::string c_dumpfilename;
	std::ofstream debugoutstr;
	std::ofstream c_dump;

	std::vector<double> node_foc;
	std::vector<double> IonicStrengths;
	std::vector<std::vector<double> > ActivityCoefficients; //  For Mineral Dissolution, store gammas in 1 data object
	int activity_model;

	/* Methods */
	CKinReactData(void);
	~CKinReactData(void);

	/* Class Read Function */
	/**
	 * reading input data for kinetic reactions
	 * @param in input file stream
	 * @param geo_obj object of class GEOObjects that manages the geometric entities
	 * @param unique_name the name of the project to access the right geometric entities
	 * @return
	 */
	bool Read(std::ifstream* in, const GEOLIB::GEOObjects& geo_obj, const std::string& unique_name);
	void Write(std::ofstream*); /* Class Write Function */
	void TestWrite(void);
	void ExecuteKinReact(void);
	void Biodegradation(double* m_Conc, long node, double eps, double hmin, double* usedtneu, int* nok, int* nbad,
	                    int tt, int tsteps);

	// CB ReactDeact
	void ReactionDeactivation(long); // Sets nodes active / inactive
	void ReactDeactPlotFlagsToTec();
	void ReactDeactSetOldReactionTerms(long nonodes);
	void SortIterations(long* iterations, long* indexes, long len);

	void CopyConcentrations(void); // SB 01.2011
	double** concentrationmatrix;
	void Aromaticum(long nonodes);

	// CB 12/09 new data structures for parallel computation
	void NAPLDissolutionPreprocessing();
	void PreprocessMinKin();
	void PostprocessMinKin();
	void Calc_linearized_rates(double* m_Conc, long Number_of_Components, double deltaT, long node);

	// CB _drmc_
	void PreprocessMicrobe_drmc_(double steplength);
};

extern std::vector<CKinReact*> KinReact_vector; // declare extern instance of class CKinReact
// declare extern instance of class CKinReact
extern std::vector<CKinReactData*> KinReactData_vector;
extern std::vector<CKinBlob*> KinBlob_vector; // declare extern instance of class Blob

// CB _drmc_
extern std::vector<MicrobeData*> MicrobeData_vector; // declare extern instance of class Blob

/**
 * read file for kinetic reaction
 * @param file_base_name base file name (without extension) containing the data
 * @param geo_obj object of class GEOObjects managing the geometric entities
 * @param unique_name unique name to access the geometric entities in geo_obj
 * @return false if file can not be opened, else true
 */
bool KRRead(const std::string& file_base_name, const GEOLIB::GEOObjects& geo_obj, const std::string& unique_name);

extern bool KRWrite(std::string const&);
extern void KRCDelete(void);

/**
 * configure kinetic reaction
 * @param geo_obj object of class GEOObjects managing the geometric entities
 * @param unique_name unique name to access the geometric entities in geo_obj
 */
void KRConfig(const GEOLIB::GEOObjects& geo_obj, const std::string& unique_name);

/**
 * configure Blob-Object
 * @param geo_obj object of class GEOObjects managing the geometric entities
 * @param unique_name unique name to access the geometric entities in geo_obj
 */
void KBlobConfig(const GEOLIB::GEOObjects& geo_obj, const std::string& unique_name);

extern void KBlobCheck(void); /* check Blob-Object for input errors */
extern bool KNaplDissCheck(void); /* CB check if NAPL dissolution is modeled */
extern void KNaplCalcDensity(void); /* CB check if NAPL dissolution is modeled */
extern double CalcNAPLDens(int node); // CB 01/08
extern void CalcNewNAPLSat(); // CB 01/08
extern void CalcNAPLCompMasses();
extern void CalcNewPhasePressure();
extern bool KMinKinCheck(void);

// CB _drmc_
void MicrobeConfig(void);

/* Externe Subroutine-Deklarationen fuer Bulirsch-Stoer Gleichungsl�ser */

extern bool odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1, double hmin, double* nexth,
                   int* nok, int* nbad, void (*derivs)(double, double[], double[], int, long, double),
                   bool (*stifbs)(double[], double[], int, double*, double, double, double[], double*, double*,
                                  void (*)(double, double[], double[], int, long, double), long),
                   bool (*rkqs)(double[], double[], int, double*, double, double, double[], double*, double*,
                                void (*)(double, double[], double[], int, long, double), long),
                   long, int);

extern bool stifbs(double y[], double dydx[], int nv, double* xx, double htry, double eps, double yscal[], double* hdid,
                   double* hnext, void (*derivs)(double, double[], double[], int, long, double), long);

extern bool rkqs(double y[], double dydx[], int n, double* x, double htry, double eps, double yscal[], double* hdid,
                 double* hnext, void (*derivs)(double, double[], double[], int, long, double), long);

extern double* dvector(long nl, long nh);
extern void free_dvector(double* v, long nl, long nh);

/* interne Deklarationen */
extern void ExecuteKineticReactions();
extern double TBCGetExchange(long node, int sp);
extern void derivs(double x, double y[], double dydx[], int n, long node, double);
extern void jacobn(double x, double y[], double dfdx[], double** dfdy, int n, long node);
#endif
