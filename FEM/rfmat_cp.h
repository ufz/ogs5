/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************/
/* ROCKFLOW - Modul: rfmat_cp.h
 */
/* Task:
   Methods for ComponentProperties
 */
/* Programming:
   10/2004   SB  First Implemented
 */
/**************************************************************************/

/* Schutz gegen mehrfaches Einfuegen */
#ifndef rfmat_cp_INC
#define rfmat_cp_INC

#include "ProcessInfo.h"
#define CP_FILE_EXTENSION ".mcp" /* File extension for component properties input file */
#include <fstream>
#include <map>
#include <string>
#include <vector>

/*************************************************************************

   Class ComponentProperties

 **************************************************************************/

class CompProperties : public ProcessInfo
{
private:
public:
	/* constructor */
	CompProperties();
	/* destructor */
	~CompProperties(void);

	size_t idx; /* the unique index of this component. not effective, saved for future*/
	std::string compname; /* component name */
	std::string iupac_formula; // chemical formula
	int mobil; /* flag mobil */
	int transport_phase; /* number of phase, in which component is transported */
	int fluid_phase;
	int valence; // valence of ionic elements /*MX*/
	double molar_mass;
	double pc; // Critical pressure [Pa]
	double omega; // accentric factor [-]
	double Tc; // Critical temperature [K]
	int fluid_id; // Requred to detect a particular fluid from *.mcp files: 0 for CO2; 1 for H2O; 2 for CH4; 3 for N2
	double Vm; // Vm [m3/kmol]: Molar volume used in interation parameter calculation
	double Vd; //[cm3/mol] diffusion volume
	int OutputMassOfComponentInModel; // 05/2012 BG

	/* Diffusionsmodelle und zugehoerige Beschreibungswerte */
	int diffusion_model; /* Zerfallsmodell in geloester Phase */
	int count_of_diffusion_model_values; /* Anzahl der Parameter zur Spezifikation des Diffusionsmodells */
	double diffusion_model_values[10]; /* Parameter fuer das Diffusionsmodell */
	double diffusion_anisotropy_ratio[3]; /* Ratio of three diffusion coefficients */ // PCH
	int diffusion_function_name;
	/* Zugriff auf Number of Parameters */
	int GetNumberDiffusionValuesCompProperties(int);

	/* Zerfallsmodelle und zugehoerige Beschreibungswerte in der geloesten Phase */
	int decay_model; /* Zerfallsmodell in geloester Phase */
	int count_of_decay_model_values; /* Anzahl und Werte zur Spezifikation der */
	double decay_model_values[10]; /* Parameter fuer Zerfallsprozess wie z.B. Zerfallsrate */
	int decay_function_name;
	int GetNumberDecayValuesCompProperties(int); /* Zugriff auf Number of Parameters */

	/* Sorption */
	int isotherm_model; /* Isothermen-Typ */
	int count_of_isotherm_model_values; /* Anzahl der Isothermen-Koeffizienten */
	double isotherm_model_values[10]; /* Isothermen-Koeffizienten */
	int isotherm_function_name;
	/* Zugriff auf Number of Parameters */
	int GetNumberIsothermValuesCompProperties(int);
	/* bubble velocity */
	int bubble_velocity_model;
	double bubble_velocity[3]; /* velocity of rising bubbles */

	/* parameters for NAPL dissolution CB140708 */
	double molar_density;
	double molar_weight;
	double max_solubility;
	/* parameters for mineral kinetics */
	double a_zero;

	double mineral_density;

//// CB _ctx_
//    bool _ctx_;
//    std::string ct_substratename;
//    double _ctx_Coefficient;

#ifdef GEM_REACT
	// kg44 25.11.2008 kinetics...for coupling with GEMS
	//
	int kinetic_model; // only 1 = GEMS implemented right now
	int n_activities; // number of species for activities
	std::string active_species[10]; // name for species ...maximum 10 names
	double kinetic_parameters[41];
	//	0,1,2  double E_acid,E_neutral,E_base; // activation energies
	//      3-5  double k_acid, k_neutral,k_base; // dissolution/precipitation rate constants
	//      6-11  double p1,q1,p2,q2,p2,q2; // exponents for omega
	//      12,13, 14  double n_1, n_2, n_3; // exponents for acidic neutral and base cases for species one
	//      append for each species another set of n_1, n_2, n_3 (up to 10 sets -> up to ten species)
	int surface_model; // currently only 1 implemented
	double surface_area[10];
	double GetNodePorosityValue_MT(long node_Index, int timelevel);
#endif

	std::ios::pos_type Read(std::ifstream*); /* Lesefunktion f? eine Instanz von CompProperties */
	void Write(std::ofstream*); /* Schreibfunktion f? eine Instanz von CompProperties */

	/* Member - Functions */
	double CalcDiffusionCoefficientCP(long index, double theta, CRFProcess* m_pcs);
	double CalcDiffusionCoefficientCP_Method1(long index, double T, double P, double eta);
	// OK411 double CalcElementRetardationFactor( long index, double*gp, double theta );
	double CalcElementRetardationFactorNew(long index, double* gp, CRFProcess* m_pcs);
	double CalcElementMeanConc(long index);
	double CalcElementMeanConcNew(long index, CRFProcess* m_pcs);
	double CalcElementDecayRate(long index);
	double CalcElementDecayRateNew(long index, CRFProcess* m_pcs);
	// IO
	std::string file_base_name;
};

/**
 * HS 01.2011 instead of the old cp_vec:
 * extern vector <CompProperties*> cp_vec;
 * a new cp_vec is introduced with map structure.
 * cp_vec[i] still points the i-th component.
 */
extern std::map<int, CompProperties*> cp_vec;
/**
 * also introduce a second map structure, that store the relationship
 * of component name and its index value.
 * this will give access to the user who would like to get access
 * to a particular component by name.
 */
extern std::map<std::string, int> cp_name_2_idx;
/**
 * Here is a straight forward example:
 * cp_vec[cp_name_2_idx["O2"]];
 */

/* ----------------------------------------------------------------------- */
/* Read all Component Properties instance by instance from input file *.cp */
extern bool CPRead(std::string);
/* Write all Component Properties instance by instance to output file */
extern void CPWrite(std::string, int);
extern int CPGetMobil(long comp);
extern void MCPDelete();
#endif
