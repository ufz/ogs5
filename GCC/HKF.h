/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

//#pragma once

#include <string>
#include <vector>
// using namespace std;

class HKF
{
private:
public:
	HKF(void);
	~HKF(void);

	/* Data */
	typedef struct
	{
		// input
		std::string name0; // name in Helgeson datafile
		// from datafile
		std::string name;
		std::string abbrev; // abbreviation
		std::string scform; // structural chemical formula
		std::string ecform; // elemental chemical formula
		std::string ref[2];
		double charge; // charge
		int type; // phase flag
		// 0-minerals that do not undergo phase transitions
		// 1-minerals that undergo one phase transition
		// 2-minerals that undergo two phase transitions
		// 3-minerals that undergo three phase transitions
		// 4-gases
		// 5-aqueous species
		double param[9][4]; // parameters
		// calculating result
		double S; // entropy
		double H; // enthalpy
		double G; // Gibbs free energy
		double V; // volume
		double Cp; // heat capacity
	} SpeciesData;

	static std::vector<SpeciesData> HKFspecies;

	/* Methods */
	static int load_param(std::vector<SpeciesData>& spec, std::vector<int> name_type, std::vector<int> phase_type);
	// species storage
	// name type: 0-name, 1-scform, 2-abbrev, 3-ecform, others-all
	// phase type: 0-5 SpeciesData.type, others- last found in database
	static int HelgesonEquation(double T, double P, std::vector<SpeciesData>& spec);

	static double HKFcalc(double T, double P, int ghs, std::string name, int sub, int type);
	// T K, P bar,
	// ghs 0-gibbs energy 1-enthalpy 2-entropy,
	// name, sub 0-name 1-scform 2-abbrev 3-ecform other-all,
	// type 0-min 1-min 2-min 3-min 4-gases 5-aq others-last found in list

	static double HKFcalcw(double T, double P, int ghs);

	static int reaction_convert(void);

	static double diel_water(double T, double P);
	static double Qfunc(double T, double P);
	static double Ufunc(double T, double P);
	static double Nfunc(double T, double P);
	static double Xfunc(double T, double P);
	static double Yfunc(double T, double P);
	static double Ydiff(double T, double P);
	static double gfunc(double T, double P);

	static int HKF_OGS_loadparam(std::string species_name0, int& type, double& charge, double param[][4]);
	static int HKF_OGS_calc(double T, double P, double& G, double& H, double& S, int type, double charge,
	                        double param[][4]);

	static void entrance(void);
};
