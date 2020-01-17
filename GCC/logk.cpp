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
#include <limits>
#include <cmath>
#include "logk.h"
#include "IAPWS-IF97.h"
#include "VLE.h"
using namespace std;

// calculating ion product of NaCl
double LOGK::KsNaCl(double T, double P)
{
    return 1.188919E+03 - 2.636581E+00 * T - 1.561700E+05 / T +
           1.795657E-03 * T * T + 3.237887E-03 * P + 5.409459E-05 * T * P +
           7.200665E-06 * P * P - 2.992790E-08 * T * P * P;
    // Knacl=1188.919-2.636581*B1-156170/B1+0.001795657*B1*B1+0.003237887*C1+0.00005409459*B1*C1+0.000007200665*C1*C1-0.0000000299279*B1*C1*C1
}

// calculating ion product of pure water
double LOGK::logKw(double T, double P)
{
    double a1, a2, a3, a4, a5, a6, a7, A, B;
    a1 = -4.098;
    a2 = -3245.2;
    a3 = 2.2362e5;
    a4 = -3.984e7;
    a5 = 13.957;
    a6 = -1262.3;
    a7 = 8.5641e5;
    A = a1 + a2 / T + a3 / (T * T) + a4 / (T * T * T);
    B = a5 + a6 / T + a7 / (T * T);
    return A + B * log10(IF97::density(T, P / 10.0) / 1.0e3);
}

// calculating logK1 logK2
double LOGK::logK1(double T, double P)
{
    double dV1, dk1, lnKs1;
    dV1 = LOGK::dV1(T);
    dk1 = LOGK::dk1(T);
    lnKs1 = LOGK::lnKs1(T);
    return LOGK::lnK(T, P, dV1, dk1, lnKs1) / 2.302585;
}

double LOGK::logK2(double T, double P)
{
    double dV2, dk2, lnKs2;
    dV2 = LOGK::dV2(T);
    dk2 = LOGK::dk2(T);
    lnKs2 = LOGK::lnKs2(T);
    return LOGK::lnK(T, P, dV2, dk2, lnKs2) / 2.302585;
}

double LOGK::lnK(double T, double P, double dV, double dk, double lnKs)
{
    double Ps;
    if (T >= 373.15)
        Ps = IF97::Psat(T) * 10.0;
    else
        Ps = 1.0;
    return -1.0 / 8.314 / T *
               (-1.0 * dV / 10.0 * (P - Ps) +
                0.5 * dk / 10000.0 * (P - Ps) * (P - Ps)) +
           lnKs;
}

double LOGK::lnKs1(double T)
{
    double a1, a2, a3, a4, a5;
    a1 = 233.5159304;
    a2 = 0.0;
    a3 = -11974.38348;
    a4 = 0.0;
    a5 = -36.50633536;
    return a1 + a2 * T + a3 / T + a4 / T / T + a5 * log(T);
}

double LOGK::dV1(double T)
{
    double a1, a2, a3, a4, a5;
    a1 = -3748.1678;
    a2 = 0.0;
    a3 = 177207.90885;
    a4 = 0.0;
    a5 = 558.25496;
    return a1 + a2 * T + a3 / T + a4 / T / T + a5 * log(T);
}

double LOGK::dk1(double T)
{
    double a1, a2, a3, a4, a5;
    a1 = -1395.81946;
    a2 = 0.0;
    a3 = 66772.55024;
    a4 = 0.0;
    a5 = 206.23006;
    return a1 + a2 * T + a3 / T + a4 / T / T + a5 * log(T);
}

double LOGK::lnKs2(double T)
{
    double a1, a2, a3, a4, a5;
    a1 = -151.1815202;
    a2 = -0.088695577;
    a3 = -1362.259146;
    a4 = 0.0;
    a5 = 27.79798156;
    return a1 + a2 * T + a3 / T + a4 / T / T + a5 * log(T);
}

double LOGK::dV2(double T)
{
    double a1, a2, a3, a4, a5;
    a1 = -2453.97326;
    a2 = 0.0;
    a3 = 115489.29266;
    a4 = 0.0;
    a5 = 367.46855;
    return a1 + a2 * T + a3 / T + a4 / T / T + a5 * log(T);
}

double LOGK::dk2(double T)
{
    double a1, a2, a3, a4, a5;
    a1 = -535.45092;
    a2 = 0.0;
    a3 = 27345.82051;
    a4 = 0.0;
    a5 = 78.76586;
    return a1 + a2 * T + a3 / T + a4 / T / T + a5 * log(T);
}

double LOGK::logK_CO2(double T, double P)
{
    double lnkH, Vphi, R, lnphi, Ps = 1.0;
    if (T > 373.15)
        Ps = IF97::Psat(T) * 10.0;
    R = 83.14471;  // unit V-cm^3, P-bar, T-K
    lnphi = VLE::LnPHI_CO2(T, P);
    lnkH = 1.3999520898E+01 - 1.3340507000E-02 * T - 5.5897820016E+02 / T -
           4.2257732219E+05 / T / T;
    Vphi = 35.663 - 5.960e-2 * (T - 273.15) +
           6.308e-4 * (T - 273.15) * (T - 273.15);
    return log10(exp(lnphi) * P /
                 exp(lnkH + Vphi * (P - Ps) / R /
                                T));  // return eq CO2 activity in aqueous phase
}
