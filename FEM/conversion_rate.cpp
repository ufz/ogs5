/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "conversion_rate.h"
#include <cmath>

#include "PhysicalConstant.h"

//#define SIMPLE_KINETICS //wenn definiert, dann einfache Kinetik, sonst Schaube

conversion_rate::conversion_rate(double T_solid,
                                 double T_gas,
                                 double p_gas,
                                 double x_reactive,
                                 double rho_s_initial,
                                 double phi_S,
                                 double delta_t,
                                 FiniteElement::SolidReactiveSystem system)
    : R(PhysicalConstant::IdealGasConstant),
      rho_s_0(rho_s_initial),
      p_eq(1.0),
      tol_l(1.0e-4),
      tol_u(1.0 - tol_l),
      tol_rho(0.1),
      x(Eigen::VectorXd(1))
{
    update_param(T_solid,
                 T_gas,
                 p_gas,
                 x_reactive,
                 rho_s_initial,
                 phi_S,
                 delta_t,
                 system);

    if (system == FiniteElement::CaOH2)
    {  // Definition auch in void
       // CSolidProperties::SetSolidReactiveSystemProperties()
        rho_low = 1665.1;
        rho_up = 2200.0;
        // only for equilibrium calculation. Do not touch. Reaction heat value
        // set in rf_msp
        reaction_enthalpy = -1.12e+05;  // in J/mol; negative for exothermic
                                        // composition reaction
        reaction_entropy = -143.5;      // in J/mol K
        M_carrier = PhysicalConstant::MolarMass::N2;
        M_react = PhysicalConstant::MolarMass::Water;
    }
    else if (system == FiniteElement::Mn3O4)
    {  // Definition auch in void
       // CSolidProperties::SetSolidReactiveSystemProperties()
        rho_low = 4500.0;
        rho_up = 4860.0;
        reaction_enthalpy = -1.376e+05;  // in J/mol; negative for exothermic
                                         // composition reaction
        reaction_entropy = -114.1;       // in J/mol K
        M_carrier = PhysicalConstant::MolarMass::N2;
        M_react = PhysicalConstant::MolarMass::O2;
    }
    else if (system == FiniteElement::Z13XBF)
    {  // Definition auch in void
       // CSolidProperties::SetSolidReactiveSystemProperties()
        // TODO [CL] read those values from some input file
        rho_low = 1150.0;
        rho_up = -1.0;                                // not needed
        reaction_enthalpy = 0.0;                      // see CalcEnthalpy13XBF()
        reaction_entropy = 0.0;                       // see CalcEntropy13XBF()
        M_carrier = PhysicalConstant::MolarMass::N2;  // consider switch to air
        M_react = PhysicalConstant::MolarMass::Water;
        W0 = 0.291 / 1.e3;  // in m^3/kg
        p_min = 0.0;        // in Pa
    }
}

conversion_rate::~conversion_rate(void) {}

void conversion_rate::update_param(double T_solid,
                                   double T_gas,
                                   double p_gas,
                                   double x_reactive,
                                   double rho_s_initial,
                                   double phi_S,
                                   double delta_t,
                                   FiniteElement::SolidReactiveSystem system)
{
    conversion_rate::T_s = T_solid;
    conversion_rate::T = T_gas;
    conversion_rate::p_gas = p_gas;  // should be in unit bar
    conversion_rate::x_react = x_reactive;
    conversion_rate::rho_s = rho_s_initial;
    x(0) = rho_s_initial;
    conversion_rate::phi_solid = phi_S;
    conversion_rate::dt = delta_t;
    conversion_rate::reaction_system = system;
}

// determine equilibrium temperature and pressure according to van't Hoff
void conversion_rate::set_chemical_equilibrium()
{
    X_D = (rho_s - rho_up - tol_rho) / (rho_low - rho_up - 2.0 * tol_rho);
    X_D = (X_D < 0.5)
              ? std::max(tol_l, X_D)
              : std::min(X_D, tol_u);  // constrain to interval [tol_l;tol_u]

    X_H = 1.0 - X_D;

    // calculate equilibrium
    // using the p_eq to calculate the T_eq - Clausius-Clapeyron
    // T_eq = (reaction_enthalpy/R) / ((reaction_entropy/R) + log(p_r_g)); //
    // unit of p in bar Alternative: Use T_s as T_eq and calculate p_eq - for
    // Schaube kinetics p_eq = exp((reaction_enthalpy/R)/T_s -
    // (reaction_entropy/R)); Schaube / Ostermeier
    const double A(-12845.), B(16.508);
    T_eq = A / (log(p_r_g) - B);
    p_eq = exp(A / T_s + B);
}

// determine equilibrium loading according to Dubinin
void conversion_rate::set_sorption_equilibrium()
{
    // determine adsorption potential
    const double A = get_potential(T_s, p_r_g);
    // determine adsorbed volume
    const double W = characteristic_curve(A);
    // determine equilibrium loading
    C_eq = W * get_adsorbate_density(T_s);  // kg/kg
}

double conversion_rate::get_mole_fraction(double xm)
{
    return M_carrier * xm / (M_carrier * xm + M_react * (1.0 - xm));
}

void conversion_rate::calculate_qR()
{
    // Convert mass fraction into mole fraction
    const double mol_frac_react = get_mole_fraction(x_react);

    switch (reaction_system)
    {
        case FiniteElement::
            CaOH2:  // Definition auch in void
                    // CSolidProperties::SetSolidReactiveSystemProperties()
        {
            p_r_g = std::max(mol_frac_react * p_gas,
                             1.0e-3);  // avoid illdefined log
            set_chemical_equilibrium();
            const double dXdt = Ca_hydration();
            qR = (rho_up - rho_low) * dXdt;
        }
        break;

        case FiniteElement::
            Mn3O4:  // Definition auch in void
                    // CSolidProperties::SetSolidReactiveSystemProperties()
        {
            p_r_g = std::max(mol_frac_react * p_gas,
                             1.0e-3);  // avoid illdefined log
            set_chemical_equilibrium();
            const double dXdt = Mn_redox();
            qR = (rho_up - rho_low) * dXdt;
        }
        break;

        case FiniteElement::Z13XBF:
        {
            // partial pressure
            p_r_g =
                std::max(mol_frac_react * p_gas * 1.0e5,
                         p_min);  // avoid illdefined log, gas pressure in Pa
            set_sorption_equilibrium();
            const double dCdt = Z13XBF_adsorption();
            qR = rho_low * dCdt;
        }
        break;

        default:
            qR = 0.;
            break;
    }
}

void conversion_rate::set_rho_s(double new_rho_s)
{
    rho_s = new_rho_s;
}

double conversion_rate::get_qR()
{
    return qR;
}

void conversion_rate::get_x(Eigen::VectorXd& output_x)
{
    output_x = x;
}

void conversion_rate::eval(double /*t*/,
                           Eigen::VectorXd const& y,
                           Eigen::VectorXd& dydx)
{
    assert(y.size() == dydx.size());

    set_rho_s(y(0));
    calculate_qR();
    dydx(0) = get_qR();
}

double conversion_rate::Ca_hydration()
{
    double dXdt;
// step 3, calculate dX/dt
#ifdef SIMPLE_KINETICS
    if (T_s < T_eq)  // hydration - simple model
#else
    if (p_r_g > p_eq)  // hydration - Schaube model
#endif
    {
// X_H = max(tol_l,X_H); //lower tolerance to avoid oscillations at onset of
// hydration reaction. Set here so that no residual reaction rate occurs at end
// of hydration.
#ifdef SIMPLE_KINETICS  // this is from P. Schmidt
        dXdt = -1.0 * (1.0 - X_H) * (T_s - T_eq) / T_eq * 0.2 *
               conversion_rate::x_react;
#else  // this is from Schaube
        if (X_H == tol_u || rho_s == rho_up)
            dXdt = 0.0;
        else if ((T_eq - T_s) >= 50.0)
            dXdt = 13945.0 * exp(-89486.0 / R / T_s) *
                   std::pow(p_r_g / p_eq - 1.0, 0.83) * 3.0 *
                   (X_D)*std::pow(-1.0 * log(X_D), 0.666);
        else
            dXdt =
                1.0004e-34 * exp(5.3332e4 / T_s) * std::pow(p_r_g, 6.0) * (X_D);
#endif
    }
    else  // dehydration
    {
// X_D = max(tol_l,X_D); //lower tolerance to avoid oscillations at onset of
// dehydration reaction. Set here so that no residual reaction rate occurs at
// end of dehydration.
#ifdef SIMPLE_KINETICS  // this is from P. Schmidt
        dXdt = -1.0 * (1.0 - X_D) * (T_s - T_eq) / T_eq * 0.05;
#else
        if (X_D == tol_u || rho_s == rho_low)
            dXdt = 0.0;
        else if (X_D < 0.2)
            dXdt = -1.9425e12 * exp(-1.8788e5 / R / T_s) *
                   std::pow(1.0 - p_r_g / p_eq, 3.0) * (X_H);
        else
            dXdt = -8.9588e9 * exp(-1.6262e5 / R / T_s) *
                   std::pow(1.0 - p_r_g / p_eq, 3.0) * 2.0 * std::pow(X_H, 0.5);
#endif
    }
    return dXdt;
}

double conversion_rate::Mn_redox()
{
    double dXdt, Avrami;

    if (T_s < T_eq)  // oxidation
    {
        if (X_H == tol_u || rho_s == rho_up)
            dXdt = 0.0;
        else
        {
            Avrami = 282.277 - 0.912 * T_s + 9.949e-4 * T_s * T_s -
                     3.620e-7 * T_s * T_s * T_s;
            dXdt = 55271.0 * exp(-95.493e3 / (R * T_s)) * Avrami *
                   (X_D)*std::pow(-1.0 * log(X_D), (Avrami - 1.0) / Avrami) *
                   std::pow(1.0 - T_s / T_eq, 0.86);
        }
        if (p_r_g <= 1.0e-3)
            dXdt = 0.0;
    }
    else  // reduction
    {
        dXdt = 0.0;
    }
    return dXdt;
}

double conversion_rate::Z13XBF_adsorption()
{
    const double k_rate = 6.0e-3;             // to be specified
    const double C = rho_s / rho_low - 1.;    // current degree of loading
    const double dCdt = k_rate * (C_eq - C);  // scaled with mass fraction
    if (dCdt > 0. && p_r_g < p_min * 1.01)
        return 0.;
    // automatic time stepping should be used instead of the following.
    // else if (dCdt > 0.) {
    //	double dens_guess = p_r_g*COMP_MOL_MASS_WATER/(R*T); //vapor density
    // guess
    //	//const double max_rate = dens_guess/(rho_low*dt);
    //	//dCdt = min(dCdt,max_rate);
    //}
    return dCdt;
}

// Density model for water found in the works of Hauer
double conversion_rate::get_adsorbate_density(const double Tads)
{
    // set reference state for adsorbate EOS in Hauer
    const double T0 = 293.15, rho0 = 998.084,
                 alpha0 = 2.06508e-4;  // K; kg/m^3; 1/K
    // double test = (rho0 * (1. - alpha0 * (Tads-T0)));
    return (rho0 * (1. - alpha0 * (Tads - T0)));  // in kg/m^3
}

// Thermal expansivity model for water found in the works of Hauer
double conversion_rate::get_alphaT(const double Tads)
{
    // set reference state for adsorbate EOS in Hauer
    const double T0 = 293.15,
                 /*rho0 = 998.084,*/ alpha0 = 2.06508e-4;  // K; kg/m^3; 1/K
    return (alpha0 / (1. - alpha0 * (Tads - T0)));         // in 1/K
}

// Saturation pressure for water used in Nunez
double conversion_rate::get_ps(const double Tads)
{
    // critical T and p
    const double Tc = 647.3;    // K
    const double pc = 221.2e5;  // Pa
    // dimensionless T
    const double Tr = Tads / Tc;
    const double theta = 1. - Tr;
    // empirical constants
    const double c[] = {-7.69123, -26.08023, -168.17065, 64.23285, -118.96462,
                        4.16717,  20.97506,  1.0e9,      6.0};
    const double K[] = {c[0] * theta + c[1] * pow(theta, 2) +
                            c[2] * pow(theta, 3) + c[3] * pow(theta, 4) +
                            c[4] * pow(theta, 5),
                        1. + c[5] * theta + c[6] * pow(theta, 2)};

    const double exponent =
        K[0] / (K[1] * Tr) - theta / (c[7] * pow(theta, 2) + c[8]);
    return pc * exp(exponent);  // in Pa
}

// Evaporation enthalpy of water from Nunez
double conversion_rate::get_hv(double Tads)  // in kJ/kg
{
    Tads -= PhysicalConstant::CelsiusZeroInKelvin;
    if (Tads <= 10.)
    {
        const double c[] = {2.50052e3,   -2.1068,     -3.57500e-1,
                            1.905843e-1, -5.11041e-2, 7.52511e-3,
                            -6.14313e-4, 2.59674e-5,  -4.421e-7};
        double hv = 0.;
        for (size_t i = 0; i < sizeof(c) / sizeof(c[0]); i++)
            hv += c[i] * std::pow(Tads, static_cast<int>(i));
        return hv;
    }
    else if (Tads <= 300.)
    {
        const double c[] = {2.50043e3,  -2.35209,    1.91685e-4,  -1.94824e-5,
                            2.89539e-7, -3.51199e-9, 2.06926e-11, -6.4067e-14,
                            8.518e-17,  1.558e-20,   -1.122e-22};
        double hv = 0.;
        for (size_t i = 0; i < sizeof(c) / sizeof(c[0]); i++)
            hv += c[i] * std::pow(Tads, static_cast<int>(i));
        return hv;
    }
    else
    {
        const double c[] = {2.99866e3, -3.1837e-3,  -1.566964e1,
                            -2.514e-6, 2.045933e-2, 1.0389e-8};
        return ((c[0] + c[2] * Tads + c[4] * pow(Tads, 2)) /
                (1. + c[1] * Tads + c[3] * pow(Tads, 2) + c[5] * pow(Tads, 3)));
    }
}

// evaluate specific heat capacity of adsorbate follwing Nunez
double conversion_rate::get_specific_heat_capacity(const double Tads)
{
    const double c[] = {4.224,         -3.716e-3,   9.351e-5,      -7.1786e-7,
                        -9.1266e-9,    2.69247e-10, -2.773104e-12, 1.553177e-14,
                        -4.982795e-17, 8.578e-20,   -6.12423e-23};
    double cp = 0.;
    for (unsigned i = 0; i < sizeof(c) / sizeof(c[0]); i++)
        cp += c[i] * std::pow(Tads, static_cast<int>(i));
    return cp;  // kJ/(kg*K)
}

// Evaluate adsorbtion potential A
double conversion_rate::get_potential(const double Tads, double pads)
{
    pads = std::max(pads, p_min);
    double A = R * Tads * log(get_ps(Tads) / pads) /
               (M_react * 1.e3);  // in kJ/kg = J/g
    if (A < 0.0)
    {
        // vapour partial pressure > saturation pressure
        // A = 0.0; // TODO [CL] debug output
    }
    return A;
}

// Characteristic curve. Return W (A)
double conversion_rate::characteristic_curve(const double A)
{
    // parameters from least squares fit (experimental data)
    const double c[] = {0.34102920966608297,     -0.0013106032830951296,
                        -0.00060754147575378876, 3.7843404172683339e-07,
                        4.0107503869519016e-07,  3.1274595098338057e-10,
                        -7.610441241719489e-11};
    double W =
        (c[0] + c[2] * A + c[4] * std::pow(A, 2) + c[6] * std::pow(A, 3)) /
        (1.0 + c[1] * A + c[3] * std::pow(A, 2) +
         c[5] * std::pow(A, 3));  // cm^3/g
    if (W < 0.0)
    {
        W = 0.0;  // TODO [CL] debug output
    }
    return W / 1.e3;  // m^3/kg
}

// Calculate sorption entropy
double conversion_rate::get_entropy(const double Tads, const double A)
{
    const double epsilon = 1.0e-8;
    const double dAdlnW = 2.0 * epsilon /
                          (log(characteristic_curve(A + epsilon)) -
                           log(characteristic_curve(A - epsilon)));
    return dAdlnW * get_alphaT(Tads);
}

// Calculate sorption enthalpy
double conversion_rate::get_enthalpy(const double Tads, const double pads)
{
    const double A = get_potential(Tads, pads);
    return (get_hv(Tads) + A - Tads * get_entropy(Tads, A)) *
           1000.0;  // in J/kg
}
