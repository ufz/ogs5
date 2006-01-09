/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <vector>

//using namespace std;

//density of electrolyte
class density
{
private:

public:
	density(void);
	~density(void);

/* Data */
	static double Ds, Vs; //storage of the return density (g/cm^3) & infinite dilution apparent molar volume (cm^3/mol)

/* Methods */
	static void MaoModel(double T, double P, double m, int f);
	//T (K), P (bar), m (mol/kg)
	//flag 0-LiCl, 1-NaCl, 2-KCl, 3-MgCl2, 4-CaCl2, 5-SrCl2, 6-BaCl2

	static double viscosity(double T, double P, double m, int f);
	//T (K), P (bar), m (mol/kg)  flag 0-LiCl, 1-NaCl, 2-KCl

	static double MultiDensity(double T, double P, std::vector<double> mv, std::vector<int> fv); //mv-molality fv-flag
	void Interface(void);

	static double CO2brine(double T, double P, double mNaCl, double mCO2);



	//0-Li, 1-Na, 2-K, 3-Mg, 4-Ca, 5-Cl, 6-SO4, 7-CO3
	//static void PartialVolume(double T, double P, double &Vs[8]);
	static double ExcessM(double mv[8]);
	static double Multi_density(double T, double P, double mv[8]);

	static double CO2_MultiBrine_density(double T, double P, double mv[8], double mCO2);
	static double concentration_water(double T, double P, double cv[8], double cv_CO2);

	static void entrance(void);

};
