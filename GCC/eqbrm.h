/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <iostream>
#include <iomanip>
#include <math.h>
#include <cmath>
#include <limits>
#include <stdio.h>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <ctype.h>
// using namespace std;

class EQ
{
private:
public:
	EQ(void);
	~EQ(void);

	static std::vector<std::string> species; // name of species in the last line of input file (eqbrm.in)
	static double T, P;

	/* Data */
	typedef struct
	{
		std::vector<std::string> species_name;
		std::vector<std::string> species_formula;
		std::vector<double> species_stoi;
		double logK;
		double stoi_water;
	} ReactionData;
	static std::vector<ReactionData> eq_reactions;

	typedef struct
	{
		std::string name;
		std::string formula;
		double charge;
		int idx_phase;
		std::vector<int> idx_element;
		std::vector<double> stoi_element;
	} SpeciesData;
	static std::vector<SpeciesData> eq_species;

	// set chemical system for eqbrm input value
	int NN, MM;
	double** VV; // mass action submatrix
	double** GG; // mass balance submatrix --> chemical element
	double* ZZ; // charge vector
	std::vector<std::string> chemical_species, chemical_element; // name and index of species and element

	// pass parameters
	double* LOGKK;
	double* BB; // mass amount --> chemical element
	int* xidx_phase; // phase index vector
	double *AA, *GAMMAA; // concentration and activity coefficent
	double AAW;
	double* xstoi_water; // check "H2O" in reactions, 0 no water, 1 or -1 water stoi number

	int activity_model; // 0-9

	/* Method */
	void eqbrm(int N, int M, double** V, double* B, double** G, double* LOGK, double* Z, double* A, double* GAMMA,
	           double& AW);
	void eqbrm_multiphases(int N, int M, double** V, double* B, double** G, double* LOGK, double* VW, double* Z, int* P,
	                       double* A, double* GAMMA, double& AW);
	void actcof(double N, double* Z, double* A, double* GAMMA, double& AW);
	void actcofWATEQ(double N, double* Z, double* A, double* GAMMA, double& AW);

	void ludcmp(double** a, int n, int* indx, double* d);
	void lubksb(double** a, int n, int* indx, double b[]);
	double** matrix(int row, int col);
	static double** matrix_new(int row, int col);
	void matrixdel(double** data);

	void entrance(void);
	void load_chemical_system(void);
	void show_chemical_system(void);
	void calc_chemical_system(void);

	void reset_estimated_value(void);

	void calc_activity(double* A, double* GAMMA, double& AW);

	void batch_calc(void);
	void batch_calc_pH(void);
	void single_calc(void);

	static void mainentrance(void);
};
