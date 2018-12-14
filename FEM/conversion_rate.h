/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef _CONVERSION_RATE_H
#define _CONVERSION_RATE_H

#include "Eigen/Eigen"
#include "FEMEnums.h"

class conversion_rate
{
public:
    conversion_rate(double T_solid,
                    double T_gas,
                    double p_gas,
                    double w_water,
                    double rho_s_initial,
                    double phi_S,
                    double delta_t,
                    FiniteElement::SolidReactiveSystem system);

    ~conversion_rate(void);

    void update_param(double T_solid,
                      double T_gas,
                      double p_gas,
                      double w_water,
                      double rho_s_initial,
                      double phi_S,
                      double delta_t,
                      FiniteElement::SolidReactiveSystem system);

    double get_mole_fraction(double xm);

    void eval(double t, Eigen::VectorXd const& y, Eigen::VectorXd& dydx);

    // Adsorbate EOS
    double get_adsorbate_density(const double Tads);
    // Adsorbate Dubinin stuff
    double get_enthalpy(const double Tads, const double pads);

private:
    void calculate_qR();
    double get_qR();
    void get_x(Eigen::VectorXd& output_x);
    void set_rho_s(double new_rho_s);

    // Adsorbate EOS
    double get_alphaT(const double Tads);
    double get_ps(const double Tads);
    double get_hv(const double Tads);

    double get_specific_heat_capacity(
        const double Tads);  // TODO [CL] why unused?

    // Adsorbate Dubinin stuff
    double get_potential(const double Tads, double pads);
    void set_sorption_equilibrium();
    double Z13XBF_adsorption();
    double characteristic_curve(const double A);
    double get_entropy(const double Tads, const double A);

    double Ca_hydration();
    double Mn_redox();
    void set_chemical_equilibrium();

    // const double rho_caoh2, rho_cao; // density for Ca(OH)2 and CaO
    const double R;        // [J/mol/K]
    double rho_s;          // solid phase density
    const double rho_s_0;  // initial solid phase density
    double phi_solid;      // solid volume fraction
    double p_gas;          // gas phase pressure;
    double p_r_g;          // pressure of H2O on gas phase;
    double p_eq;           // equilibrium pressure; // unit in bar
    double T_eq;           // equilibrium temperature;
    double T_s;            // solid phase temperature;
    double T;              // gas phase temperature;
    double qR;             // rate of solid density change;
    double x_react;        // mass fraction of water in gas phase;
    double X_D;      // mass fraction of dehydration (CaO) in the solid phase;
    double X_H;      // mass fraction of hydration in the solid phase;
    double dt;       // time step size;
    double rho_low;  // lower density limit
    double rho_up;   // upper density limit
    double reaction_enthalpy;
    double reaction_entropy;
    double M_carrier;  // inert component molar mass
    double M_react;    // reactive component molar mass
    double W0;         // maximum specific adsorbed volume
    double C_eq;       // equilibrium loading of sorbens
    double p_min;      // minimum pressure for which we find a valid adsorption
                       // potential
    FiniteElement::SolidReactiveSystem reaction_system;

    double tol_l;
    double tol_u;
    double tol_rho;

    Eigen::VectorXd x;
    size_t i;  // iterator
};

#endif
