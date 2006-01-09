/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

 //using namespace std;
class LOGK
{
private:
public:
	static double KsNaCl(double T, double P);//NaCl -> Na + Cl
	static double logKw(double T, double P);//H2O -> H + OH
	static double logK1(double T, double P);//H2O + CO2 -> HCO3 + H
	static double logK2(double T, double P);//HCO3 -> CO3 + H
	static double logK_CO2(double T, double P); //CO2.gas -> CO2,aq

	static double lnKs1(double T);
	static double dV1(double T);
	static double dk1(double T);
	static double lnKs2(double T);
	static double dV2(double T);
	static double dk2(double T);
	static double lnK(double T, double P, double dV, double dk, double lnKs);


};
