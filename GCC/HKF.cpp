/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/* References
(1).Shock, E. L., Oelkers, E. H., Johnson, J. W., Sverjensky, D. A., and Helgeson, H. C., 1992.
    Calculation of the Thermodynamic Properties of Aqueous Species at High-Pressures and Temperatures-
    Effective Electrostatic Radii, Dissociation-Constants and Standard Partial Molal Properties to
    1000-Degrees-C and 5-Kbar.
    Journal of the Chemical Society-Faraday Transactions 88, 803-826.

(2).Oelkers, E. H. and Helgeson, H. C., 1988. Calculation of the Thermodynamic and Transport-
    Properties of Aqueous Species at High-Pressures and Temperatures - Aqueous Tracer Diffusion-
    Coefficients of Ions to 1000-Degrees-C and 5-Kb.
    Geochimica Et Cosmochimica Acta 52, 63-85.

(3).Duan, Z. H., Moller, N., and Weare, J. H., 1992.
    An Equation of State for the Ch4-Co2-H2o System .1. Pure Systems from 0-Degrees-C to 1000-Degrees-C
    and 0 to 8000 Bar.
    Geochimica Et Cosmochimica Acta 56, 2605-2617.

(4).Helgeson, H. C., 1992. Effects of complex formation in flowing fluids on the hydrothermal
    solubilities of minerals as a function of fluid pressure and temperature in the
    critical and supercritical regions of the system H2O.
    Geochimica Et Cosmochimica Acta 56, 3191-3207.
*/

#include <limits>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <vector>
#include <stdio.h>
#include <ctype.h>
#include "VLE.h"
#include "NR.h"
#include "HKF.h"
#include "IO.h"
#include "IAPWS-IF97.h"
//#include "slop98.h"

using namespace std;

HKF::HKF(void)
{
}
HKF::~HKF(void)
{
}

int HKF::HKF_OGS_loadparam(std::string species_name0, int& type, double& charge, double param[][4])
{
	int i, j;
	HKF::SpeciesData spec;
	vector<HKF::SpeciesData> species_vector;
	vector<int> name_type, phase_type;

	spec.name0 = species_name0;
	species_vector.clear();
	species_vector.push_back(spec);
	name_type.clear();
	name_type.push_back(-1);
	phase_type.clear();
	phase_type.push_back(-1);

	if (HKF::load_param(species_vector, name_type, phase_type) != 1)
	{
		// cout << " warning... " << endl;
		// type=-1;
		// charge=-99;
		// for(i=0;i<9;i++)
		//	for(j=0;j<4;j++)
		//		param[i][j]=-99;
		return 0;
	}
	else
	{
		spec = species_vector.back();
		type = spec.type;
		charge = spec.charge;
		for (i = 0; i < 9; i++)
			for (j = 0; j < 4; j++)
				param[i][j] = spec.param[i][j];
		return 1;
	}
}

int HKF::HKF_OGS_calc(double T, double P, double& G, double& H, double& S, int type, double charge, double param[][4])
{
	int i, j, res = 1;
	HKF::SpeciesData spec;
	vector<HKF::SpeciesData> species_vector;
	vector<int> name_type, phase_type;

	spec.type = type;
	spec.charge = charge;
	for (i = 0; i < 9; i++)
		for (j = 0; j < 4; j++)
			spec.param[i][j] = param[i][j];

	species_vector.clear();
	species_vector.push_back(spec);
	if (HelgesonEquation(T, P, species_vector) == 1)
	{
		spec = species_vector.back();
		G = spec.G;
		H = spec.H;
		S = spec.S;
	}
	return res;
}

double HKF::HKFcalcw(double T, double P, int ghs)
{
	double GG, HH, SS, res;
	if (ghs > 2 || ghs < 0)
		exit(0);
	HH = (IF97::H(T, P / 10.0) / 55.51 - 287.45532) * 1000.0 / 4.18;
	SS = (IF97::S(T, P / 10.0) / 55.51 + 0.0632944) * 1000.0 / 4.18;
	GG = -125005.0 - HH;
	if (ghs == 0)
		res = GG;
	else if (ghs == 1)
		res = HH;
	else if (ghs == 2)
		res = SS;
	return res;
}

double HKF::HKFcalc(double T, double P, int ghs, std::string name, int sub, int type)
{
	if (ghs > 2 || ghs < 0)
		exit(0);
	SpeciesData ion;
	vector<HKF::SpeciesData> Ions;
	vector<int> ntv, ptv;
	int is_OK;
	double res;
	ion.name0 = name;
	Ions.clear();
	ntv.clear();
	ptv.clear();
	Ions.push_back(ion);
	ntv.push_back(sub);
	ptv.push_back(type);
	is_OK = HKF::load_param(Ions, ntv, ptv);
	if (is_OK != 1)
		exit(0);
	HelgesonEquation(T, P, Ions);
	if (ghs == 0)
		res = Ions[0].G;
	else if (ghs == 1)
		res = Ions[0].H;
	else if (ghs == 2)
		res = Ions[0].S;
	return res;
}

int HKF::load_param(vector<SpeciesData>& spec, vector<int> name_type, vector<int> phase_type)
{
	if ((int)spec.size() == (int)name_type.size() && (int)spec.size() == (int)phase_type.size())
	{
		int i, j, js, jx, jz, type = -1, is_spec, is_OK = 0;
		vector<string> pies, piesx, piesz, slop98;
		vector<int> id_spec; // storage index of name0, which can be found in the species list.
		int cek, cekz; // check species
		string first_char;

		id_spec.clear();
		is_OK = IO::file2vector("slop98.dat", slop98);
		// for(i=0;i<(int)slop98.size();i++){
		//	cin.get();
		//	cout << slop98[i] <<endl;
		//}

		if (is_OK)
		{
			for (i = 0; i < (int)slop98.size(); i++)
			{
				first_char = slop98[i].substr(0, 3);
				if (first_char == "min" || first_char == "gas" || first_char == "aqu")
				{
					type++;
					i++;
				}

				piesx = IO::string2vector(slop98[i]);
				piesz = IO::string2vector(slop98[i + 1]);

				for (j = 0; j < (int)spec.size(); j++)
				{
					// cout << i << " " << slop98[i] << " OK " << endl;
					is_spec = 0;

					// searching depend on name scform ecform and abbrev
					if (name_type[j] == 0 && (int)piesx.size() > 0)
					{
						if (spec[j].name0 == piesx[0])
							is_spec = 1;
					}
					else if (name_type[j] == 1 && (int)piesx.size() > 1)
					{
						if (spec[j].name0 == piesx[1])
							is_spec = 1;
					}
					else if (name_type[j] == 2 && (int)piesz.size() > 0)
					{
						if (spec[j].name0 == piesz[0])
							is_spec = 1;
					}
					else if (name_type[j] == 3 && (int)piesz.size() > 1)
					{
						if (spec[j].name0 == piesz[1])
							is_spec = 1;
					}
					else
					{
						// cout << " OK load " << endl;
						for (js = 0; js < (int)piesx.size(); js++)
							if (spec[j].name0 == piesx[js])
								is_spec = 1;
						for (js = 0; js < (int)piesz.size(); js++)
							if (spec[j].name0 == piesz[js])
								is_spec = 1;
					}

					if (is_spec == 1)
					{
						// searching depend on phase type

						is_spec = 0;
						if (phase_type[j] == type || phase_type[j] < 0 || phase_type[j] > 5)
							is_spec = 1;
					}

					if (is_spec == 1)
					{
						id_spec.push_back(j);
						spec[j].name = piesx[0];
						spec[j].scform = piesx[1];
						if ((int)piesz.size() >= 2)
						{
							spec[j].abbrev = piesz[0];
							spec[j].ecform = piesz[1];
						}
						else
						{
							spec[j].abbrev = "N/A";
							spec[j].ecform = piesz[0];
						}
						pies = IO::string2vector(slop98[i + 2]);
						spec[j].ref[0] = pies[0];
						spec[j].ref[1] = pies[1];

						spec[j].type = type;
						switch (type)
						{
							case 0:
							case 4:
							case 5:
								for (jx = 0; jx < 3; jx++)
								{
									pies = IO::string2vector(slop98[i + jx + 3]);
									for (jz = 0; jz < 4; jz++)
										if (jz < (int)pies.size())
											spec[j].param[jx][jz] = atof(pies[jz].c_str());
										else
											spec[j].param[jx][jz] = 0.0;
								}
								for (jx = 3; jx < 9; jx++)
									for (jz = 0; jz < 4; jz++)
										spec[j].param[jx][jz] = 0.0;
								spec[j].charge = spec[j].param[2][3];
								break;

							case 1:
								pies = IO::string2vector(slop98[i + 3]);
								for (jz = 0; jz < 4; jz++)
									spec[j].param[0][jz] = atof(pies[jz].c_str());

								pies = IO::string2vector(slop98[i + 4]);
								for (jz = 0; jz < 3; jz++)
									spec[j].param[1][jz] = atof(pies[jz].c_str());
								spec[j].param[1][3] = 0.0;

								for (jz = 0; jz < 4; jz++)
									spec[j].param[2][jz] = atof(pies[jz + 3].c_str());

								pies = IO::string2vector(slop98[i + 5]);
								for (jz = 0; jz < 3; jz++)
									spec[j].param[3][jz] = atof(pies[jz].c_str());
								spec[j].param[3][3] = 0.0;

								spec[j].param[4][0] = atof(slop98[i + 6].c_str());
								for (jz = 1; jz < 4; jz++)
									spec[j].param[4][jz] = 0.0;

								for (jx = 5; jx < 9; jx++)
									for (jz = 0; jz < 4; jz++)
										spec[j].param[jx][jz] = 0.0;
								spec[j].charge = 0.0;
								break;

							case 2:
								pies = IO::string2vector(slop98[i + 3]);
								for (jz = 0; jz < 4; jz++)
									spec[j].param[0][jz] = atof(pies[jz].c_str());

								pies = IO::string2vector(slop98[i + 4]);
								for (jz = 0; jz < 3; jz++)
									spec[j].param[1][jz] = atof(pies[jz].c_str());
								spec[j].param[1][3] = 0.0;
								for (jz = 0; jz < 4; jz++)
									spec[j].param[2][jz] = atof(pies[jz + 3].c_str());

								pies = IO::string2vector(slop98[i + 5]);
								for (jz = 0; jz < 3; jz++)
									spec[j].param[3][jz] = atof(pies[jz].c_str());
								spec[j].param[3][3] = 0.0;
								for (jz = 0; jz < 4; jz++)
									spec[j].param[4][jz] = atof(pies[jz + 3].c_str());

								pies = IO::string2vector(slop98[i + 6]);
								for (jz = 0; jz < 3; jz++)
									spec[j].param[5][jz] = atof(pies[jz].c_str());
								spec[j].param[5][3] = 0.0;

								spec[j].param[6][0] = atof(slop98[i + 7].c_str());
								for (jz = 1; jz < 4; jz++)
									spec[j].param[6][jz] = 0.0;

								for (jx = 7; jx < 9; jx++)
									for (jz = 0; jz < 4; jz++)
										spec[j].param[jx][jz] = 0.0;
								spec[j].charge = 0.0;
								break;

							case 3:
								pies = IO::string2vector(slop98[i + 3]);
								for (jz = 0; jz < 4; jz++)
									spec[j].param[0][jz] = atof(pies[jz].c_str());

								pies = IO::string2vector(slop98[i + 4]);
								for (jz = 0; jz < 3; jz++)
									spec[j].param[1][jz] = atof(pies[jz].c_str());
								spec[j].param[1][3] = 0.0;
								for (jz = 0; jz < 4; jz++)
									spec[j].param[2][jz] = atof(pies[jz + 3].c_str());

								pies = IO::string2vector(slop98[i + 5]);
								for (jz = 0; jz < 3; jz++)
									spec[j].param[3][jz] = atof(pies[jz].c_str());
								spec[j].param[3][3] = 0.0;
								for (jz = 0; jz < 4; jz++)
									spec[j].param[4][jz] = atof(pies[jz + 3].c_str());

								pies = IO::string2vector(slop98[i + 6]);
								for (jz = 0; jz < 3; jz++)
									spec[j].param[5][jz] = atof(pies[jz].c_str());
								spec[j].param[5][3] = 0.0;
								for (jz = 0; jz < 4; jz++)
									spec[j].param[6][jz] = atof(pies[jz + 3].c_str());

								pies = IO::string2vector(slop98[i + 7]);
								for (jz = 0; jz < 3; jz++)
									spec[j].param[7][jz] = atof(pies[jz].c_str());
								spec[j].param[7][3] = 0.0;

								spec[j].param[8][0] = atof(slop98[i + 8].c_str());
								for (jz = 1; jz < 4; jz++)
									spec[j].param[8][jz] = 0.0;
								spec[j].charge = 0.0;
								break;

							default:
								break;
						} // end switch
					} // end if
				} // end for

				switch (type)
				{
					case 0:
					case 4:
					case 5:
						i += 5;
						break;
					case 1:
						i += 6;
						break;
					case 2:
						i += 7;
						break;
					case 3:
						i += 8;
						break;
					default:
						break;
				}
			}
			if ((int)id_spec.size() > 0)
			{
				cekz = 1;
				for (i = 0; i < (int)spec.size(); i++)
				{
					cek = 0;
					for (j = 0; j < (int)id_spec.size(); j++)
						if (id_spec[j] == i)
							cek = 1;
					if (cek == 0)
					{
						cekz = 0;
						cout << endl
						     << " Warning! " << spec[i].name0 << " (" << name_type[i] << ") is NOT in the species list!"
						     << endl;
					}
				}
				if (cekz == 1)
					return 1; // OK
				else
					return -1; // species name is incorrect
			}
			else
				return -2; // did nothing
		}
		else
			return 0; // can NOT open species list vector
	}
	else
		exit(0);
}

int HKF::HelgesonEquation(double T, double P, vector<SpeciesData>& spec)
{
	double c1, c2, a1, a2, a3, a4, Z, dG, dH, dS, dV; // aqueous species parameters from database
	double Sc1, Sc2, Sa3, Sa4, Sw;
	double Hc1, Hc2, Ha1, Ha2, Ha3, Ha4, Hw;
	double Gc1, Gc2, Ga1, Ga2, Ga3, Ga4, Gw;
	double a0, b0, c0 /*, Tmax*/; // gases and minerals parameters from database

	double theta, psi, eta, Tr, Pr, e, er, Yr, Y, g; // independent on species
	double rejr, rej, wr, w, X1; // dependent on species
	int i, is_aq = 0;

	Tr = 298.15; // p819 Ref.(1)
	Pr = 1.0; // p819 Ref.(1)

	for (i = 0; i < (int)spec.size(); i++)
		if (spec[i].type == 5)
			is_aq = 1;

	if (is_aq == 1)
	{
		theta = 228.0; // p819 Ref.(1)
		psi = 2600.0; // p819 Ref.(1)
		eta = 1.66027E5; // p824 Ref.(1)

		e = diel_water(T, P); //(C1) p821 Ref.(1)
		er = diel_water(Tr, Pr);
		Yr = Yfunc(Tr, Pr); // Yr = -5.8E-5; p824 Ref.(1)
		Y = Yfunc(T, P); //(A25) Ref.(1)
		g = gfunc(T, P); // Eq(24) p807 Ref.(1)

		Sc1 = log(T / Tr);
		Sc2 = -(1 / (T - theta) - 1 / (Tr - theta) + 1 / theta * log(Tr * (T - theta) / T / (Tr - theta))) / theta;
		Sa3 = 1 / (T - theta) / (T - theta) * (P - Pr);
		Sa4 = 1 / (T - theta) / (T - theta) * log((psi + P) / (psi + Pr));
		//(A20) & (B21) Ref.(1)

		Hc1 = T - Tr;
		Hc2 = -(1.0 / (T - theta) - 1.0 / (Tr - theta));
		Ha1 = P - Pr;
		Ha2 = log((psi + P) / (psi + Pr));
		Ha3 = (2 * T - theta) / pow((T - theta), 2.0) * (P - Pr);
		Ha4 = (2 * T - theta) / pow((T - theta), 2.0) * log((psi + P) / (psi + Pr));
		// A(22) & B(23) Ref.(1)

		Gc1 = -(T * log(T / Tr) - T + Tr);
		Gc2 = -(((1.0 / (T - theta)) - (1.0 / (Tr - theta))) * (theta - T) / theta
		        - T / theta / theta * log(Tr * (T - theta) / T / (Tr - theta)));
		Ga1 = P - Pr;
		Ga2 = log((psi + P) / (psi + Pr));
		Ga3 = 1.0 / (T - theta) * (P - Pr);
		Ga4 = 1.0 / (T - theta) * log((psi + P) / (psi + Pr));
		// A(24) & B(25) Ref.(1)
	}
	else
	{
		Yr = g = e = er = Y = 0.0;
		Sc1 = Sc2 = Sa3 = Sa4 = 0.0;
		Hc1 = Hc2 = Ha1 = Ha2 = Ha3 = Ha4 = 0.0;
		Gc1 = Gc2 = Ga1 = Ga2 = Ga3 = Ga4 = 0.0;
	}

	for (i = 0; i < (int)spec.size(); i++)
	{
		if (spec[i].type == 5)
		{
			dG = spec[i].param[0][0];
			dH = spec[i].param[0][1];
			dS = spec[i].param[0][2];
			a1 = spec[i].param[1][0] * 1.0E-1;
			a2 = spec[i].param[1][1] * 1.0E+2;
			a3 = spec[i].param[1][2];
			a4 = spec[i].param[1][3] * 1.0E+4;
			c1 = spec[i].param[2][0];
			c2 = spec[i].param[2][1] * 1.0E+4;
			wr = spec[i].param[2][2] * 1.0E+5;
			Z = spec[i].charge;
			if (Z == 0)
			{
				Sw = wr * (Y - Yr); //(A20) Ref.(1)
				Gw = wr * (Yr * (T - Tr) + 1.0 / e - 1.0 / er); //(A24) Ref.(1)
				Hw = wr * (T * Y - Tr * Yr + 1.0 / e - 1.0 / er); //(A22) Ref.(1)
			}
			else
			{
				rejr = Z * Z * (eta * Yr - 100.0) / (dS - 71.5 * abs(Z)); //(D1) Ref.(1) & Eq(29) Ref.(4)
				rej = rejr + abs(Z) * g; // Eq(6) Ref.(1)
				wr = eta * Z * (Z / rejr - 1.0 / 3.082); //(D5) Ref.(1)
				w = eta * Z * (Z / (rejr + abs(Z) * g) - 1.0 / (3.082 + g)); //(D4) Ref.(1)
				X1 = -eta * (pow(abs(Z), 3.0) * pow(rej, -2.0) - Z * pow((3.082 + g), -2.0)); //(B19) Ref.(1)
				Sw = w * Y - (1.0 / e - 1.0) * X1 * NR::dfridrX(gfunc, T, P) - wr * Yr; //(B21) & (B14) Ref.(1)
				Hw = w * (1.0 / e - 1.0) + w * T * Y - T * (1.0 / e - 1.0) * X1 * NR::dfridrX(gfunc, T, P)
				     - wr * (1.0 / er - 1.0) - wr * Tr * Yr;
				Gw = w * (1.0 / e - 1.0) - wr * (1.0 / er - 1.0) + wr * Yr * (T - Tr);
			}
			spec[i].S = dS + c1 * Sc1 + c2 * Sc2 + a3 * Sa3 + a4 * Sa4 + Sw; //(A20) (B21) Ref.(1)
			spec[i].G = dG - dS * (T - Tr) + c1 * Gc1 + c2 * Gc2 + a1 * Ga1 + a2 * Ga2 + a3 * Ga3 + a4 * Ga4
			            + Gw; //(A24) (B25) Ref.(1)
			spec[i].H = dH + c1 * Hc1 + c2 * Hc2 + a1 * Ha1 + a2 * Ha2 + a3 * Ha3 + a4 * Ha4 + Hw; //(A22) (B23) Ref.(1)
		}
		else
		{
			dG = spec[i].param[0][0];
			dH = spec[i].param[0][1];
			dS = spec[i].param[0][2];
			dV = spec[i].param[0][3];
			a0 = spec[i].param[1][0];
			b0 = spec[i].param[1][1] * 1.0e-3;
			c0 = spec[i].param[1][2] * 1.0e5;

			switch (spec[i].type)
			{
				case 0:
				case 4:
					// Tmax= spec[i].param[2][0];
					spec[i].G = dG - dS * (T - Tr) + a0 * (T - Tr - T * log(T / Tr))
					            + (-1.0 * c0 - b0 * T * Tr * Tr) * (T - Tr) * (T - Tr) / (2.0 * T * Tr * Tr)
					            + dV * (P - Pr) * 0.1 / 4.18; // volume, pressure and energy unit convert factor
					spec[i].H = dH + a0 * (T - Tr) + b0 / 2.0 * (T * T - Tr * Tr) - c0 * (1.0 / T - 1.0 / Tr)
					            + dV * (P - Pr) * 0.1 / 4.18;
					spec[i].S = dS + a0 * log(T / Tr) + b0 * (T - Tr) - c0 / 2.0 * (1.0 / T / T - 1.0 / Tr / Tr);
					break;
				case 1:
				case 2:
				case 3:
				default:
					spec[i].G = 0.0;
					spec[i].H = 0.0;
					spec[i].S = 0.0;
					break;
			}
		}
	}

	return 1;
}

// calculting the dielectric constant (relative permittivity) of pure water from T,P 1000C 5000bar Ref.(1) Appendix C
double HKF::diel_water(double T, double P)
{ // T K, P bar
	double k0, k1, k2, k3, k4, Tr, z, dens;
	double a1, a2, a3, a4, a5, a6, a7, a8, a9, a10;
	Tr = 298.15;
	z = T / Tr;
	dens = VLE::density_H2O(T, P); // Ref.(3)
	a1 = 14.70333593; // Table C1 Ref.(1)
	a2 = 212.8462733;
	a3 = -115.4445173;
	a4 = 19.55210915;
	a5 = -83.30347980;
	a6 = 32.13240048;
	a7 = -6.694098645;
	a8 = -37.86202045;
	a9 = 68.87359646;
	a10 = -27.29401652;
	k0 = 1.0; //(C2)
	k1 = a1 / z; //(C3)
	k2 = a2 / z + a3 + a4 * z; //(C4)
	k3 = a5 / z + a6 * z + a7 * z * z; //(C5)
	k4 = a8 / z / z + a9 / z + a10; //(C6)
	return k0 + k1 * dens + k2 * dens * dens + k3 * pow(dens, 3.0) + k4 * pow(dens, 4.0); //(C1)
}

double HKF::Ydiff(double T, double P)
{
	return 1.0 / diel_water(T, P);
}

double HKF::Yfunc(double T, double P)
{
	return -NR::dfridrX(Ydiff, T, P);
}

double HKF::gfunc(double T, double P)
{
	double Dw, ag, bg, a1, a2, a3, Tf, Pf;
	Dw = VLE::density_H2O(T, P); // Ref.(3)
	if (Dw >= 1.0)
		return 0.0; // Fig.6 & p811 Ref.(1)
	else
	{
		T = T - 273.15;
		ag = -2.037662 + 5.747000E-3 * T - 6.557892E-6 * T * T; // Eq(25) Table 3  p807 Ref.(1)
		bg = 6.107361 - 1.074377E-2 * T + 1.268348E-5 * T * T; // Eq(26)
		if (T < 155.0 || P > 1000.0)
			return ag * pow((1.0 - Dw), bg); // Eq(24) Ref.(1) Fig.6
		else
		{
			a1 = 3.66666E-16; // Table.4 Ref.(1)
			a2 = -1.504956E-10;
			a3 = 5.01799E-14;
			Tf = (T - 155.0) / 300.0;
			Pf = 1000.0 - P;
			return ag * pow((1.0 - Dw), bg)
			       - (pow(Tf, 4.8) + a1 * pow(Tf, 16))
			             * (a2 * pow(Pf, 3.0) + a3 * pow(Pf, 4.0)); // Eq(32) & Eq(33) Ref.(1)
		}
	}
}
