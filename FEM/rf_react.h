/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/* reaction package header file */

#ifndef reaction_INC
#define reactions_INC

#include <vector>

/* Structure for exchange of reaction rates */
class REACT
{
private:
public:
	REACT(void);
	~REACT(void);
	// Data
	const char** name;                    /* names of reactants */
	double** val_in;                      /* input concentration values to reaction module */
	double** val_out;                     /* output concentration values of reaction module */
	double** rate;                        /* reaction rates for MTM2 */
	int* rateflag;                        /* flag used for determining if reaction are calculated with PHREEQC */
	int countsteps;                       /* number of timesteps, after which reactions are recalculated in the whole model domain */
	int number_of_comp;                   /* Number of components in above data structures: size = number_of_comp * NodeListLength */
	int heatflag;                         /* if 1, heat transport is active and passed to phreeqc */
	long nodenumber;                      /* number of nodes, on which reactions are calculated */
	bool flag_pqc;                        /* flag if *.pqc file exists */
	bool check_no_reaction_nodes;         /* flag if CheckNoReactionNodes has been performed */
	double temperature;                   /* temperature of water at node as specified in input file */
	long elenumber;                       //number of elements
	/* hier später arrays of reactions reinhängen ?*/

	// rcml moved here
	int rcml_number_of_master_species;    /* number of master species (inorgan. equilibrium) */
	int rcml_number_of_equi_phases;       /* number of phases (in equilibrium) */
	int rcml_number_of_kinetics;          /* number of kinetic reactions  */
	int rcml_number_of_ion_exchanges;     /* number of phases (in equilibrium) */
	int rcml_number_of_gas_species;       /* number of species in gas phase */
	size_t rcml_number_of_secondary_species; /* number of species in gas phase */
	/* hier später reaction models reinhängen ?*/
	int rcml_pH_flag;                     /* =0, pH constant; =1 (default), pH will change  */
	int rcml_pe_flag;                     /* =0, pe constant; =1 (default), pe will change  */
	int rcml_heat_flag;                   /* =0, temp constant (default); =1 , temp will change  */
	int rcml_number_of_pqcsteps;          /* Anzahl der Reaktionsschritte in PHREEQC aus dem Befehl: -steps "time" in "pqcsteps" steps */
	int rcml_pH_charge;                   /* =0, no charge balance for pH; =1, used for charge balance (keyword charge in line with pH*/
	char* outfile;                        /* Ausgabefile von PHREEQC */
	std::string file_name_pqc;            // Name of pqc file in GeoSys project (*.pqc)
      std::string file_name_database;             // Name of pqc database file in GeoSys project (*.pqc)
	std::string outfile_name;
	std::string results_file_name;
	std::vector < std::string > pqc_names; // species names in *-pqc input file
	std::vector < int > pqc_index;        // index in process array
	std::vector < int > pqc_process;      // process number in pcs_vector
	double gamma_Hplus;                   //activity coefficent of H+ ion
      std::vector < std::string > additional_punches;

	// Member functions
	REACT* GetREACT(void);
	void CreateREACT(void);
	void InitREACT(void);
	void ExecuteReactionsPHREEQC(void);
	void ExecuteReactionsPHREEQCNew(void);
	void TestPHREEQC(std::string);
	int  Call_Phreeqc(void);
	void GetTransportResults(void);
	int  ReadReactionModel(FILE* File);
	int  ReadReactionModelNew(std::ifstream*);
	//fsout removed 3912
	int  ReadInputPhreeqc( long index, FILE* fpqc, FILE* Fphinp);
	int  WriteInputPhreeqc(long, /*ifstream*,*/ std::ofstream*);
	int  ReadOutputPhreeqc(char* fout);
	int  ReadOutputPhreeqcNew(void);
	void ResetpHpe(void);
	void CalculateReactionRates(void);
	void SetConcentrationResults(void);
	void CalculateReactionRateFlag(void);
	void SetNeighborNodesActive(long startnode, long level, int* help);
	int  CheckNoReactionNodes(void);
      int Teststeps(long nodes);
	// Reaction at elements //MX
	void InitREACT0(void);
	void ExecuteReactionsPHREEQC0(void);
	void SetConcentrationResultsEle(void);
	void GetTransportResults2Element(void);
#ifdef LIBPHREEQC
	// MDL: libphreeqc
	void ExecuteReactionsPHREEQCNewLib(void); // MDL:
	                                          // MDL:
	int  WriteInputPhreeqcLib(long, std::stringstream*, int*);
	int  ReadOutputPhreeqcNewLib(double*); // MDL:
	int  Call_PhreeqcLib(int, int, int, std::stringstream*, double*);
#endif                                         // LIBPHREEQC
};
extern std::vector <REACT*> REACT_vec;

extern void DestroyREACT(void);
extern void RCRead(std::string);
extern double MATCalcIonicStrengthNew(long index);
extern void REACTInit();                          //OK

//Water moles per kg of water
//#define MOLH2OPERKG 55.50843506

#endif

#ifdef LIBPHREEQC
// MDL: libphreeqc
//#pragma comment (lib, "libphreeqc.lib")	 // LB Commented out, old stuff???
#endif
