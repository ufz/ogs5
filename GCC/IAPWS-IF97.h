/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// IAPWS-IF97 for density and viscosity
class IF97
{
private:
public:
    IF97(void);
    ~IF97(void);

    /* Data */
    static double R, Rm, M, Tc, Pc, Dc, Tt, Pt, Tb;
    static double TT, PP;

    /* Methods */
    static void ReferenceConstants(void);
    static double Psat(double);
    static double Tsat(double);
    static double Pb23(double);
    static double Tb23(double);
    static int region(double, double);

    static double G(double, double);  //(T K, P Mpa)
    static double H(double, double);  //(T K, P Mpa)
    static double S(double, double);  //(T K, P Mpa)

    static double density(double, double);     // density (kg m^-3) (T K, P Mpa)
    static double viscosity(double, double);   // viscosity (Pa s)  (T K, P Mpa)
    static double dielectric(double, double);  // dielectric constant ()

    static double g1PT(double, double);
    static double g2PT(double, double);
    static double f3DT(double, double);
    static double dpressure(double);
    static void entrance(void);
};
