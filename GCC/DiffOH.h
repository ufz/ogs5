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

class Diffusion
{
private:
public:
	Diffusion(void);
	~Diffusion(void);

	/* Data */
	int index1;
	double t0Na, t0Cl; // transference numbers of Na and Cl ions Eq(52) Ref.(2)
	double AeNaCl, AeCaCl2, AeHCl, AeNaOH, AeNaBr, AeNaI, AeMgCl2,
	    AeHBr; // conductance of references electrolyte Eq(46) Ref.(2)
	typedef struct
	{
		std::string N; // name
		int Z; // charge
		double A; // conductance
		double D; // diffusion coefficient
		double S; // entropy
		double H; // enthalpy
		double G; // Gibbs free energy
		std::string M; // name in Helgeson datafile
		double C[11]; // parameters
	} AqData;
	std::vector<AqData> Ions;

	/* Methods */
	void IonsA(std::string, double); // write Conductance value
	void IonsD(std::string, double); // write Diff Coef value
	double IonsA(std::string); // return Conductance value
	double IonsD(std::string); // return Diff Coef value
	double IonsS(std::string); // return entropy value

	void InputIons(void); // add Ions into vector
	void TransferenceNumberNaCl(double); // calculate Transference Number NaCl
	void ConductancesRef(double, double); // calculate Conductances of Ref electrolytes

	void LoadParameters(void); // load in parameters for every species from database (slop98.dat)

	void HelgesonEquation(double, double);

	void DiffusionCoefficientIon(double);
	void PrintOut(void);

	void DiffCoefOH(double, double);
};

double gfun(double, double);
double diel_water(double, double);

double Qfun(double, double);
double Ufun(double, double);
double Nfun(double, double);
double Xfun(double, double);
double Yfun(double, double);
