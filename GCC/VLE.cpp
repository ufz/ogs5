/**
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <iostream>
#include <math.h>
#include <cmath>
#include <limits>
#include "VLE.h"
#include "PITZdata.h"
#include "NR.h"
#include "IAPWS-IF97.h"
#include "Brent/brent.hpp"

using namespace std;

double VLE::TT;
double VLE::PP;

VLE::VLE(void) {}
VLE::~VLE(void) {}

// unit T K, P MPa
double VLE::Psat(double T)
{
    double A, B, C, theta, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10;
    n1 = 0.11670521452767e4;
    n2 = -0.72421316703206e6;
    n3 = -0.17073846940092e2;
    n4 = 0.12020824702470e5;
    n5 = -0.32325550322333e7;
    n6 = 0.14915108613530e2;
    n7 = -0.48232657361591e4;
    n8 = 0.40511340542057e6;
    n9 = -0.23855557567849;
    n10 = 0.65017534844798e3;
    theta = T + n9 / (T - n10);
    A = pow(theta, 2.0) + n1 * theta + n2;
    B = n3 * pow(theta, 2.0) + n4 * theta + n5;
    C = n6 * pow(theta, 2.0) + n7 * theta + n8;
    return pow(2.0 * C / (-B + pow(B * B - 4.0 * A * C, 0.5)), 4.0);
}

//====> solubility section <====
double VLE::solubility_CO2(double T, double P, double mNaCl)
{
    return (P - Psat(T) * 10.0) /
           exp(-LnPHI_CO2(T, P) + u0_CO2(T, P) + lnrCO2(T, P, mNaCl));
}

double VLE::solubilityNEW_CO2(double T, double P, double mNaCl)
{
    double a, b, mCO2, dev, err = 1.0e-6;
    int i, iter_max = 100;
    mCO2 = VLE::solubility_CO2(T, P, mNaCl);
    a = mCO2 * 0.5;
    b = mCO2 * 1.5;
    for (i = 0; i < iter_max; i++)
    {
        mCO2 = (a + b) * 0.5;
        dev = mCO2 -
              (P - Psat(T) * 10.0) / exp(-LnPHI_CO2(T, P) + uNEW_CO2(T, P) +
                                         LGAMMA_CO2(T, P, mCO2, mNaCl));
        if (abs(dev) < err)
            break;
        else if (dev < 0.0)
            a = mCO2;
        else if (dev > 0.0)
            b = mCO2;
    }
    return mCO2;
}

double VLE::u0_CO2(double T, double P)
{
    double c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11;
    c1 = 2.89447706e+01;
    c2 = -3.54581768e-02;
    c3 = -4.77067077e+03;
    c4 = 1.02782768e-05;
    c5 = 3.38126098e+01;
    c6 = 9.04037140e-03;
    c7 = -1.14934031e-03;
    c8 = -3.07405726e-01;
    c9 = -9.07301486e-02;
    c10 = 9.32713393e-04;
    c11 = 0.0;
    return c1 + c2 * T + c3 / T + c4 * pow(T, 2) + c5 / (630 - T) + c6 * P +
           c7 * P * log(T) + c8 * P / T + c9 * P / (630 - T) +
           c10 * pow(P, 2) / pow((630 - T), 2) + c11 * T * log(P);
}
double VLE::lnrCO2(double T, double P, double mNaCl)
{
    double c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, d1, d2, d3, d4, d5, d6,
        d7, d8, d9, d10, d11, res;
    c1 = -4.11370585e-01;
    c2 = 6.07632013e-04;
    c3 = 97.5347708;
    c4 = 0.0;
    c5 = 0.0;
    c6 = 0.0;
    c7 = 0.0;
    c8 = -2.37622469e-02;
    c9 = 1.70656236e-02;
    c10 = 0.0;
    c11 = 1.41335834e-05;
    d1 = 3.36389723e-04;
    d2 = -1.98298980e-05;
    d3 = 0.0;
    d4 = 0.0;
    d5 = 0.0;
    d6 = 0.0;
    d7 = 0.0;
    d8 = 2.12220830e-03;
    d9 = -5.24873303e-03;
    d10 = 0.0;
    d11 = 0.0;
    res = 2.0 * mNaCl *
          (c1 + c2 * T + c3 / T + c4 * pow(T, 2) + c5 / (630 - T) + c6 * P +
           c7 * P * log(T) + c8 * P / T + c9 * P / (630 - T) +
           c10 * pow(P, 2) / pow((630 - T), 2) + c11 * T * log(P));
    res += mNaCl * mNaCl *
           (d1 + d2 * T + d3 / T + d4 * pow(T, 2) + d5 / (630 - T) + d6 * P +
            d7 * P * log(T) + d8 * P / T + d9 * P / (630 - T) +
            d10 * pow(P, 2) / pow((630 - T), 2) + d11 * T * log(P));
    return res;
}
double VLE::LGAMMA_CO2(double T, double P, double mCO2, double mNaCl)
{
    double LAMN, LAM, ZETA;
    PITZdata M;
    LAMN = M.pitzer_parameters(T, P, "LAMN_CO2_CO2");
    LAM = M.pitzer_parameters(T, P, "LAM_Na_CO2");
    ZETA = M.pitzer_parameters(T, P, "ZETA_NaCl_CO2");
    return 2.0 * LAMN * mCO2 + 2.0 * LAM * mNaCl + ZETA * mNaCl * mNaCl;
}

double VLE::solubility_CH4(double T, double P, double mNaCl)
{
    return (P - Psat(T) * 10.0) /
           exp(u0_CH4(T, P) - LnPHI_CH4(T, P) + lnrCH4(T, P, mNaCl));
}
double VLE::u0_CH4(double T, double P)
{
    double c1, c2, c3, c4, c5, c6, c7, c8, c9;
    c1 = 0.83143711e+01;
    c2 = -0.72772168e-03;
    c3 = 0.21489858e+04;
    c4 = -0.14019672e-04;
    c5 = -0.66743449e+06;
    c6 = 0.76985890e-02;
    c7 = -0.50253331e-05;
    c8 = -0.30092013e+01;
    c9 = 0.48468502e+03;
    return c1 + c2 * T + c3 / T + c4 * T * T + c5 / T / T + c6 * P +
           c7 * P * T + c8 * P / T + c9 * P / T / T;
}
double VLE::lnrCH4(double T, double P, double mNaCl)
{
    double c1, c2, c3, c4, c5;
    c1 = -0.81222036e+00;
    c2 = 0.10635172e-02;
    c3 = 0.18894036e+03;
    c4 = 0.44105635e-04;
    c5 = -0.46797718e-10;
    return 2.0 * (c1 + c2 * T + c3 / T + c4 * P + c5 * P * P * T) * mNaCl -
           0.29903571e-02 * mNaCl * mNaCl;
}

double VLE::uNEW_CO2(double T, double P)
{
    double lnkH, Vphi;
    lnkH = 1.3999520898e+01 - 1.3340507000e-02 * T - 5.5897820016e+02 / T -
           4.2257732219e+05 / T / T;
    Vphi = 35.663 - 5.960e-2 * (T - 273.15) + 6.308e-4 * pow((T - 273.15), 2.0);
    return lnkH + Vphi * (P - Psat(T) * 10.0) / (83.14471 * T);
}

//==========================================
//>       EoS section  D unit: g/cm^3      <
//==========================================
double VLE::density_CO2(double T, double P)
{
    double x1, x2;
    x1 = 0.8e1 * T / P;
    x2 = 2.0e2 * T / P;
    if (P > 20.0)
    {
        x1 = 2.0e1 * T / P;
        x2 = 2.0e2 * T / P;
    }

    if (P > 40.0)
    {
        x1 = 0.75e1 * T / P;
        x2 = 9.0e1 * T / P;
    }

    if (P > 520.0)
    {  // if(P>600.0){ //DL 2012.03.22
        x1 = 2.5e1 * T / P;
        x2 = 5.0e2 * T / P;
    }
    if (P > 1000.0)
    {
        x1 = 5.0e1 * T / P;
        x2 = 5.0e2 * T / P;
    }
    if (P > 3500.0)
    {
        x1 = 1.5e2 * T / P;
        x2 = 2.0e3 * T / P;
    }
    TT = T;
    PP = P;
    double result;
    brent::local_min(x1, x2, 1.0e-8, dZ_CO2, result);
    return 44.01 / result;
}
double VLE::density_CH4(double T, double P)
{  // g/cm^3
    double x1, x2;
    x1 = 2.0e1 * T / P;
    x2 = 2.0e2 * T / P;
    if (P > 1000.0)
    {
        x1 = 5.0e1 * T / P;
        x2 = 5.0e2 * T / P;
    }
    if (P > 3500.0)
    {
        x1 = 1.5e2 * T / P;
        x2 = 2.0e3 * T / P;
    }
    TT = T;
    PP = P;
    double result;
    brent::local_min(x1, x2, 1.0e-8, dZ_CH4, result);
    return 16.04 / result;
}
double VLE::density_H2O(double T, double P)
{  // g/cm^3
    double x1, x2, Ps;
    Ps = Psat(T) * 10.0;
    x1 = 0.8e1 * T / P;  // 0.5e1*T/P;
    x2 = 2.0e2 * T / P;
    if (P > 400)
    {
        x1 = 1.5e1 * T / P;
        x2 = 2.0e2 * T / P;
    }
    if (P > 1000.0)
    {
        x1 = 5.0e1 * T / P;
        x2 = 5.0e2 * T / P;
    }
    if (P > 3500.0)
    {
        x1 = 1.1e2 * T / P;
        x2 = 2.0e3 * T / P;
    }
    if (T < 423.0 && P > Ps)
    {
        x1 = 1.0e1;
        x2 = 2.5e1;
    }
    else if (T < 630.0 && P > Ps)
    {
        x1 = 1.5e1;
        x2 = 3.6e1;
    }
    TT = T;
    PP = P;
    double result;
    brent::local_min(x1, x2, 1.0e-8, dZ_H2O, result);
    return 18.015 / result;
}

double VLE::dZ_CO2(double V)
{
    double R, Tc, Pc, Vc, Tr, Pr, Vr;
    R = 83.14467;
    Tc = 304.2;
    Pc = 73.825;
    Vc = R * Tc / Pc;
    Tr = TT / Tc;
    Pr = PP / Pc;
    Vr = V / Vc;
    return Z_CO2(TT, PP, V) - Pr * Vr / Tr;
}
double VLE::dZ_CH4(double V)
{
    double R, Tc, Pc, Vc, Tr, Pr, Vr;
    R = 83.14467;
    Tc = 190.6;
    Pc = 46.00;
    Vc = R * Tc / Pc;
    Tr = TT / Tc;
    Pr = PP / Pc;
    Vr = V / Vc;
    return Z_CH4(TT, PP, V) - Pr * Vr / Tr;
}
double VLE::dZ_H2O(double V)
{
    double R, Tc, Pc, Vc, Tr, Pr, Vr;
    R = 83.14467;
    Tc = 647.25;
    Pc = 221.19;
    Vc = R * Tc / Pc;
    Tr = TT / Tc;
    Pr = PP / Pc;
    Vr = V / Vc;
    return Z_H2O(TT, PP, V) - Pr * Vr / Tr;
}

double VLE::Z_CO2(double T, double /*P*/, double V)
{  // K, bar, cm^3
    double a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15;
    double R, Tc, Pc, Vc, Tr, Vr, B, C, D, E, F;
    a1 = 8.99288497E-2;
    a2 = -4.94783127E-1;
    a3 = 4.77922245E-2;
    a4 = 1.03808883E-2;
    a5 = -2.82516861E-2;
    a6 = 9.49887563E-2;
    a7 = 5.20600880E-4;
    a8 = -2.93540971E-4;
    a9 = -1.77265112E-3;
    a10 = -2.51101973E-5;
    a11 = 8.93353441E-5;
    a12 = 7.88998563E-5;
    a13 = -1.66727022E-2;
    a14 = 1.39800000E-0;
    a15 = 2.96000000E-2;
    R = 83.14467;
    Tc = 304.2;
    Pc = 73.825;
    Vc = R * Tc / Pc;
    Tr = T / Tc;
    Vr = V / Vc;
    B = a1 + a2 / pow(Tr, 2) + a3 / pow(Tr, 3);
    C = a4 + a5 / pow(Tr, 2) + a6 / pow(Tr, 3);
    D = a7 + a8 / pow(Tr, 2) + a9 / pow(Tr, 3);
    E = a10 + a11 / pow(Tr, 2) + a12 / pow(Tr, 3);
    F = a13 / pow(Tr, 3);
    return 1.0 + B / Vr + C / pow(Vr, 2) + D / pow(Vr, 4) + E / pow(Vr, 5) +
           F / pow(Vr, 2) * (a14 + a15 / pow(Vr, 2)) * exp(-a15 / pow(Vr, 2));
}
double VLE::LnPHI_CO2(double T, double P)
{
    double a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15;
    double R, V, Tc, Pc, Vc, Tr, Vr, B, C, D, E, Z;
    V = 44.01 / density_CO2(T, P);
    Z = Z_CO2(T, P, V);
    a1 = 8.99288497E-2;
    a2 = -4.94783127E-1;
    a3 = 4.77922245E-2;
    a4 = 1.03808883E-2;
    a5 = -2.82516861E-2;
    a6 = 9.49887563E-2;
    a7 = 5.20600880E-4;
    a8 = -2.93540971E-4;
    a9 = -1.77265112E-3;
    a10 = -2.51101973E-5;
    a11 = 8.93353441E-5;
    a12 = 7.88998563E-5;
    a13 = -1.66727022E-2;
    a14 = 1.39800000E-0;
    a15 = 2.96000000E-2;
    R = 83.14467;
    Tc = 304.2;
    Pc = 73.825;
    Vc = R * Tc / Pc;
    Tr = T / Tc;
    Vr = V / Vc;
    B = a1 + a2 / pow(Tr, 2) + a3 / pow(Tr, 3);
    C = a4 + a5 / pow(Tr, 2) + a6 / pow(Tr, 3);
    D = a7 + a8 / pow(Tr, 2) + a9 / pow(Tr, 3);
    E = a10 + a11 / pow(Tr, 2) + a12 / pow(Tr, 3);
    return Z - 1 - log(Z) + B / Vr + C / 2 / pow(Vr, 2) + D / 4 / pow(Vr, 4) +
           E / 5 / pow(Vr, 5) +
           a13 / (2 * pow(Tr, 3) * a15) *
               (a14 + 1 -
                (a14 + 1 + a15 / pow(Vr, 2)) * exp(-a15 / pow(Vr, 2)));
}

double VLE::Z_CH4(double T, double /*P*/, double V)
{
    double a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15;
    double R, Tc, Pc, Vc, Tr, Vr, B, C, D, E, F;
    a1 = 8.72553928E-2;
    a2 = -7.52599476E-1;
    a3 = 3.75419887E-1;
    a4 = 1.07291342E-2;
    a5 = 5.49626360E-3;
    a6 = -1.84772802E-2;
    a7 = 3.18993183E-4;
    a8 = 2.11079375E-4;
    a9 = 2.01682801E-5;
    a10 = -1.65606189E-5;
    a11 = 1.19614546E-4;
    a12 = -1.08087289E-4;
    a13 = 4.48262295E-2;
    a14 = 7.53970000E-1;
    a15 = 7.71670000E-2;
    R = 83.14467;
    Tc = 190.6;
    Pc = 46.00;
    Vc = R * Tc / Pc;
    Tr = T / Tc;
    Vr = V / Vc;
    B = a1 + a2 / pow(Tr, 2) + a3 / pow(Tr, 3);
    C = a4 + a5 / pow(Tr, 2) + a6 / pow(Tr, 3);
    D = a7 + a8 / pow(Tr, 2) + a9 / pow(Tr, 3);
    E = a10 + a11 / pow(Tr, 2) + a12 / pow(Tr, 3);
    F = a13 / pow(Tr, 3);
    return 1.0 + B / Vr + C / pow(Vr, 2) + D / pow(Vr, 4) + E / pow(Vr, 5) +
           F / pow(Vr, 2) * (a14 + a15 / pow(Vr, 2)) * exp(-a15 / pow(Vr, 2));
}
double VLE::LnPHI_CH4(double T, double P)
{
    double a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15;
    double R, V, Tc, Pc, Vc, Tr, Vr, B, C, D, E, Z;
    V = 16.04 / density_CH4(T, P);
    Z = Z_CH4(T, P, V);
    a1 = 8.72553928E-2;
    a2 = -7.52599476E-1;
    a3 = 3.75419887E-1;
    a4 = 1.07291342E-2;
    a5 = 5.49626360E-3;
    a6 = -1.84772802E-2;
    a7 = 3.18993183E-4;
    a8 = 2.11079375E-4;
    a9 = 2.01682801E-5;
    a10 = -1.65606189E-5;
    a11 = 1.19614546E-4;
    a12 = -1.08087289E-4;
    a13 = 4.48262295E-2;
    a14 = 7.53970000E-1;
    a15 = 7.71670000E-2;
    R = 83.14467;
    Tc = 190.6;
    Pc = 46.00;
    Vc = R * Tc / Pc;
    Tr = T / Tc;
    Vr = V / Vc;
    B = a1 + a2 / pow(Tr, 2) + a3 / pow(Tr, 3);
    C = a4 + a5 / pow(Tr, 2) + a6 / pow(Tr, 3);
    D = a7 + a8 / pow(Tr, 2) + a9 / pow(Tr, 3);
    E = a10 + a11 / pow(Tr, 2) + a12 / pow(Tr, 3);
    return Z - 1 - log(Z) + B / Vr + C / 2 / pow(Vr, 2) + D / 4 / pow(Vr, 4) +
           E / 5 / pow(Vr, 5) +
           a13 / (2 * pow(Tr, 3) * a15) *
               (a14 + 1 -
                (a14 + 1 + a15 / pow(Vr, 2)) * exp(-a15 / pow(Vr, 2)));
}

double VLE::Z_H2O(double T, double /*P*/, double V)
{
    double a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15;
    double R, Tc, Pc, Vc, Tr, Vr, B, C, D, E, F;
    a1 = 8.64449220E-2;
    a2 = -3.96918955E-1;
    a3 = -5.73334886E-2;
    a4 = -2.93893000E-4;
    a5 = -4.15775512E-3;
    a6 = 1.99496791E-2;
    a7 = 1.18901426E-4;
    a8 = 1.55212063E-4;
    a9 = -1.06855859E-4;
    a10 = -4.93197687E-6;
    a11 = -2.73739155E-6;
    a12 = 2.65571238E-6;
    a13 = 8.96079018E-3;
    a14 = 4.02000000;
    a15 = 2.57000000E-2;
    R = 83.14467;
    Tc = 647.25;
    Pc = 221.19;
    Vc = R * Tc / Pc;
    Tr = T / Tc;
    Vr = V / Vc;
    B = a1 + a2 / pow(Tr, 2) + a3 / pow(Tr, 3);
    C = a4 + a5 / pow(Tr, 2) + a6 / pow(Tr, 3);
    D = a7 + a8 / pow(Tr, 2) + a9 / pow(Tr, 3);
    E = a10 + a11 / pow(Tr, 2) + a12 / pow(Tr, 3);
    F = a13 / pow(Tr, 3);
    return 1.0 + B / Vr + C / pow(Vr, 2) + D / pow(Vr, 4) + E / pow(Vr, 5) +
           F / pow(Vr, 2) * (a14 + a15 / pow(Vr, 2)) * exp(-a15 / pow(Vr, 2));
}
double VLE::LnPHI_H2O(double T, double P)
{
    double a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15;
    double R, V, Tc, Pc, Vc, Tr, Vr, B, C, D, E, Z;
    V = 18.015 / density_H2O(T, P);
    Z = Z_H2O(T, P, V);
    a1 = 8.64449220E-2;
    a2 = -3.96918955E-1;
    a3 = -5.73334886E-2;
    a4 = -2.93893000E-4;
    a5 = -4.15775512E-3;
    a6 = 1.99496791E-2;
    a7 = 1.18901426E-4;
    a8 = 1.55212063E-4;
    a9 = -1.06855859E-4;
    a10 = -4.93197687E-6;
    a11 = -2.73739155E-6;
    a12 = 2.65571238E-6;
    a13 = 8.96079018E-3;
    a14 = 4.02000000;
    a15 = 2.57000000E-2;
    R = 83.14467;
    Tc = 647.25;
    Pc = 221.19;
    Vc = R * Tc / Pc;
    Tr = T / Tc;
    Vr = V / Vc;
    B = a1 + a2 / pow(Tr, 2) + a3 / pow(Tr, 3);
    C = a4 + a5 / pow(Tr, 2) + a6 / pow(Tr, 3);
    D = a7 + a8 / pow(Tr, 2) + a9 / pow(Tr, 3);
    E = a10 + a11 / pow(Tr, 2) + a12 / pow(Tr, 3);
    return Z - 1 - log(Z) + B / Vr + C / 2 / pow(Vr, 2) + D / 4 / pow(Vr, 4) +
           E / 5 / pow(Vr, 5) +
           a13 / (2 * pow(Tr, 3) * a15) *
               (a14 + 1 -
                (a14 + 1 + a15 / pow(Vr, 2)) * exp(-a15 / pow(Vr, 2)));
}

//==MIX of H2O & CO2 fagucity
// ln_Phi_H2O
double VLE::PF(double TT, double PP)
{
    double cPFphi = 0.0, lnP;
    double d1, d2, d3, d4;
    double M, N, P, Q, A, B, C, D, E;

    A = -20.7600738043549 + 0.0897803300811 * TT - 0.0000954680712 * TT * TT;
    B = 3.08284 - 0.00956 * TT + 0.0000092069 * TT * TT;
    C = -0.38986 + 0.00228 * TT - 0.0000026679 * TT * TT;
    D = 6.6623328249993 - 0.057517079835 * TT + 0.0000643503862 * TT * TT;
    E = -0.445272747863 + 0.0130486543901 * TT - 0.000015177905 * TT * TT;

    M = A + 4.0 * B + 16.0 * C;
    N = D + 5.0 * E;
    P = B + 8.0 * C;
    Q = E;

    d4 = -(2.0 * M - 2.0 * N + P + Q) / 0.0025;
    d3 = (Q - P - 0.0225 * d4) / 2.0;
    d2 = P - 8.0 * d3 + d4 / 16.0;
    d1 = M - 4.0 * d2 - 16.0 * d3 - d4 / 4.0;

    lnP = log(PP);
    if (lnP >= 5.0)
        cPFphi = D + E * lnP;
    if (lnP <= 4.0)
        cPFphi = A + B * lnP + C * lnP * lnP;
    if (lnP > 4.0 && lnP < 5.0)
        cPFphi = d1 + d2 * lnP + d3 * lnP * lnP + d4 / lnP;

    return cPFphi;
}

double VLE::fraction_H2O(double T, double P, double AW)
{
    double Z01, Ps_H2O, lnS_phi_H2O;
    Ps_H2O = VLE::Psat(T) * 10.0;
    lnS_phi_H2O = VLE::LnPHI_H2O(T, Ps_H2O * 1.00001);
    Z01 = VLE::PF(T, P);
    return exp(Z01) * AW * exp(lnS_phi_H2O) * Ps_H2O / P / P;
}

double VLE::pressure_CO2(double T, double D)
{
    double a, b, P, dev, err = 1.0e-4;
    int i, iter_max = 100;
    a = 5.0e-1;
    b = 5.0e3;

    for (i = 0; i < iter_max; i++)
    {
        P = (a + b) * 0.5;
        dev = D - VLE::density_CO2(T, P);
        if (abs(dev) < err)
            break;
        if (dev > 0)
            a = P;
        if (dev < 0)
            b = P;
    }
    return P;
}

// m=1000x/((1-x)Mw), if x-->0, ==> m=1000x/Mw, ==> kH=f/x, kM=f/m, ==>
// kM=kH*x/m, ==> kM=kH*x/(1000x/Mw), ==> kM=kH*Mw/1000
double VLE::Henry_const_CO2(double T)
{  // eq(15) return exp value, ln(kH)
    double A, B, C;
    double Tr, Tc, Tau, Ps;
    Tc = 647.096;
    Tr = T / Tc;
    Tau = 1.0 - Tr;
    A = -8.55445;
    B = 4.01195;
    C = 9.52345;
    Ps = IF97::Psat(T) * 10.0;
    return A / Tr + B / Tr * pow(Tau, 0.355) + C * pow(Tr, -0.41) * exp(Tau) +
           log(Ps) - log(55.51);  //"-log(55.51)" for unit convert from mole
                                  // fraction to mol/kg water
}

double VLE::Henry_const_H2(double T)
{  // eq(15) return exp value, ln(kH)
    double A, B, C;
    double Tr, Tc, Tau, Ps;
    Tc = 647.096;
    Tr = T / Tc;
    Tau = 1.0 - Tr;
    A = -4.73284;
    B = 6.08954;
    C = 6.06066;
    Ps = IF97::Psat(T) * 10.0;
    return A / Tr + B / Tr * pow(Tau, 0.355) + C * pow(Tr, -0.41) * exp(Tau) +
           log(Ps) - log(55.51);  //"-log(55.51)" for unit convert from mole
                                  // fraction to mol/kg water
}

void VLE::EoS_PR_H2(double T, double P, double& V, double& Z, double& lnphi)
{
    double alpha, beta, q, Psi, Omega, sigma, epsilon, w;
    double Tc, Tr, Pc, Pr;
    double Znew, Zerr = 1.0e-8;
    int i, iter_max = 20;
    double I;

    w = -0.215;
    Tc = 33.3;   // K
    Pc = 12.97;  // bar
    Tr = T / Tc;
    Pr = P / Pc;

    Omega = 0.07780;
    Psi = 0.45742;
    sigma = 1.0 + pow(2.0, 0.5);
    epsilon = 1.0 - pow(2.0, 0.5);
    alpha = pow((1.0 + (0.37464 + 1.54226 * w - 0.26992 * w * w) *
                           (1.0 - pow(Tr, 0.5))),
                2.0);

    beta = Omega * Pr / Tr;
    q = Psi * alpha / Omega / Tr;

    Z = 1.0;
    Znew = 0.0;
    for (i = 0; i < iter_max; i++)
    {
        Znew =
            1.0 + beta -
            q * beta * (Z - beta) / (Z + epsilon * beta) / (Z + sigma * beta);
        if (abs(abs(Z) - abs(Znew)) < Zerr)
            break;
        else
            Z = Znew;
    }
    V = Z * 83.14 * T / P;  // unit V-cm^3 , P-bar, T-K
    I = 1.0 / (sigma - epsilon) *
        log((Z + sigma * beta) / (Z + epsilon * beta));
    lnphi = Z - 1.0 - log(Z - beta) - q * I;
}

double VLE::solubility_H2_PR(double T, double P, double a)
{
    double Z, V, lnphi, kH, v0, R, Ps = 1.0;
    VLE::EoS_PR_H2(T, P, V, Z, lnphi);
    kH = VLE::Henry_const_H2(T);  // return ln(kH), unit bar/(mol/kg)
    R = 83.14;                    // unit V-cm^3, P-bar, T-K
    // v0 partial moler volume of H2 in water, unit cm^3/mol
    v0 = 22.79545 - 0.12791 * (T - 273.15) +
         7.16364e-4 * (T - 273.15) *
             (T - 273.15);  // basd on the Figure 10-2, in J.M. Prausnitz, MTFPE
                            // 3rd edition.
    if (T > 373.15)
        Ps = IF97::Psat(T) * 10.0;
    kH += v0 * (P - Ps) / R / T;
    return exp(lnphi) * P / exp(kH) / a;
}

double VLE::density_H2(double T, double P)
{
    double Z, V, lnphi;  // V -cm^3
    double H2_molecular_weight = 2.016;
    VLE::EoS_PR_H2(T, P, V, Z, lnphi);
    return H2_molecular_weight / V;  // unit - g/cm^3
}

double VLE::pressure_H2(double T, double dens)
{
    double a, b, P, dev, err = 1.0e-4;
    int i, iter_max = 100;
    a = 5.0e-1;
    b = 5.0e3;

    for (i = 0; i < iter_max; i++)
    {
        P = (a + b) * 0.5;
        dev = dens - VLE::density_H2(T, P);
        if (abs(dev) < err)
            break;
        if (dev > 0)
            a = P;
        if (dev < 0)
            b = P;
    }
    return P;
}
