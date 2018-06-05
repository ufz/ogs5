/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <string>
#include <vector>

// using namespace std;

class DiffusionAL
{
private:
public:
	DiffusionAL(void);
	~DiffusionAL(void);

	/* Data */
	typedef struct
	{
		string N; // name
		int Z; // charge
		double A; // parameter A or C
		double B; // parameter B
		double R; // Radii
		double M; // molecular weight
		double D0; // Self-Diffusion Coefficient
	} AqData;
	vector<AqData> Ions;

	typedef struct
	{
		string N; // name of solute
		double C; // molality of solute
		double D; // diffusion coefficient(**final result**)
	} Solute;
	vector<Solute> ListS; // initial input condition

	typedef struct
	{
		string Nc; // name of cation
		string Na; // name of anion
		double Lc; // number of cation
		double La; // number of anion
		double S; // number of solvent
		double M; // molality
	} speciesE;

	typedef struct
	{
		string N; // name of neutral species
		double L; // number of neutral species
		double S; // number of solvent
		double M; // molality
	} speciesN;

	/* Methods */
	void LoadParameters(void);
	std::vector<string> string2vector(string); // split string line to pieces, and store in a vector

	void SelfDiffCoef(double, double);
	double SingleElectrolytesDiffCoef(double, double, double, string, string, string);
	// temperature(K),pressure(MPa),molality,name of cation,name of anion,
	// name of diffusive species(destination)(return its diffusion coefficient
	int GreatestCommonDivisor(int, int);
	void DivideSolution(std::vector<Solute>, double&, std::vector<speciesE>&, std::vector<speciesN>&);
	void MultiDiffCoef(double, double);

	void CreateListS(void);

	// void testDS(double, double);
	void Interface(void);

	void entrance(void);
};
