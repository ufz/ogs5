/**
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

class VLE
{
private:
public:
    VLE(void);
    ~VLE(void);
    /* Data */
    static double TT, PP;

    /* Methods */
    static double Psat(double T);  // unit MPa

    // unit T K, P bar, V cm^3, D g/cm^3
    static double density_CO2(double T, double P);  //~1300K, 8000bar
    static double LnPHI_CO2(double T, double P);
    static double Z_CO2(double T, double P, double V);
    static double dZ_CO2(double V);
    static double density_CH4(double T, double P);  //~1300K, 8000bar
    static double LnPHI_CH4(double T, double P);
    static double Z_CH4(double T, double P, double V);
    static double dZ_CH4(double V);
    static double density_H2O(double T, double P);  //~1300K, 8000bar
    static double LnPHI_H2O(double T, double P);
    static double Z_H2O(double T, double P, double V);
    static double dZ_H2O(double V);

    static double density_H2O_CO2(double T, double P, double x);
    static double LnPHI_H2O_H2O_CO2(double T, double P, double x);
    static double LnPHI_CO2_H2O_CO2(double T, double P, double x);

    static double PF(double T, double P);
    static double fraction_H2O(double T, double P, double AW);
    static double pressure_CO2(double T, double D);

    // unit CO2 mol/kgw T K, P bar, mNaCl mol/kg
    static double solubility_CO2(double T, double P, double mNaCl);
    static double u0_CO2(double T, double P);  // unit u/RT J/mol
    static double lnrCO2(double T, double P, double mNaCl);

    static double solubilityNEW_CO2(double T, double P, double mNaCl);
    static double uNEW_CO2(double T, double P);  // unit u/RT J/mol
    static double LGAMMA_CO2(double T, double P, double mCO2, double mNaCl);

    static double solubility_CH4(double T, double P, double mNaCl);
    static double u0_CH4(double T, double P);
    static double lnrCH4(double T, double P, double mNaCl);
    static double lamda_CH4_NaCl(double T, double P);

    static double Henry_const_CO2(double T);
    static double Henry_const_H2(double T);

    static void EoS_PR_H2(double T, double P, double& V, double& Z,
                          double& lnphi);
    static double solubility_H2_PR(double T, double P, double a);

    static double density_H2(double T, double P);
    static double pressure_H2(double T, double dens);

    static void entrance(void);
};
