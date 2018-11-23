/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * zero()-function: Copyright Richard Brent, John Burkhardt | GNU LGPL license
 *
 */

/**********************************************************************
   Module: Equation of State

   Task: This file includes coefficients and functions for calculating the
   thermal properties of liquids and gases in relation to density, pressure
   and temperature.

   Programming: NB
          Aug 2008
**********************************************************************/

#include "rf_mfp_new.h"
#include "rf_mmp_new.h"
#include "fem_ele_std.h"
//#include "rf_num_new.h"
#include "eos.h"
#include "tools.h"
#include "PhysicalConstant.h"

#include "Material/Fluid/Viscosity/WaterViscosityIAPWS.h"

using namespace std;

/**********************************************************************
   Some functions and derivations as shown in [Span&Wagner 1994]
 ***********************************************************************/
double theta_fn(double tau, double A, double delta, double beta)
{
	return (1 - tau) + A * pow((delta - 1) * (delta - 1), 1 / (2 * beta));
	// TF  return (1-tau)+A*pow(pow(delta-1,2),1/(2*beta));
}

double phi_fn(double C, double delta, double D, double tau)
{
	return exp(-C * (delta - 1) * (delta - 1) - D * (tau - 1) * (tau - 1));
	// TF  return exp (-C*pow(delta-1,2)-D*pow(tau-1,2));
}

double delta_fn(double theta_fn, double B, double delta, double alpha)
{
	return theta_fn * theta_fn + B * pow((delta - 1) * (delta - 1), alpha);
	// TF  return pow(theta_fn,2)+B*pow(pow(delta-1,2),alpha);
}

double dDELTApowb_ddelta(double b, double delta_fn, double dDELTA_deriv)
{
	return b * pow(delta_fn, b - 1) * dDELTA_deriv;
}

double dDELTA_ddelta(double delta, double A, double theta_fn, double beta, double B, double a)
{
	return (delta - 1) * (A * theta_fn * 2 / beta * pow((delta - 1) * (delta - 1), (1 / (2 * beta) - 1))
	                      + 2 * B * a * pow((delta - 1) * (delta - 1), (a - 1)));
	//   return ((delta-1)*(A*theta_fn*2/beta*pow(pow((delta-1),2),(1/(2*beta)-1))+2*B*a*pow(pow((delta-1),2),(a-1))));
}

double d2DELTA_ddelta2(double delta, double dDELTA_deriv, double A, double theta_fn, double beta, double B, double a)
{
	const double delta_m1(delta - 1);
	const double delta_m1_quad((delta_m1) * (delta_m1));
	const double t0(pow(delta_m1_quad, (1 / (2 * beta) - 1)));
	return 1 / (delta - 1) * dDELTA_deriv
	       + delta_m1_quad
	             * (4 * B * a * (a - 1) * pow(delta_m1_quad, (a - 2)) + 2 * A * A * pow((1 / beta), 2) * t0 * t0
	                + A * theta_fn * 4 / beta * (1 / (2 * beta) - 1) * pow(delta_m1_quad, (1 / (2 * beta) - 2)));

	//   return 1/(delta-1)*dDELTA_deriv+pow((delta-1),2)*
	//      (4*B*a*(a-1)*pow(pow((delta-1),2),(a-2))+2*pow(A,2)*pow((1/beta),2)*pow(pow(pow((delta-1),2),(1/(2*beta)-1)),2)
	//      +A*theta_fn*4/beta*(1/(2*beta)-1)*pow(pow((delta-1),2),(1/(2*beta)-2)));
}

double d2DELTApowb_ddelta2(double b, double delta_fn, double d2DELTA_deriv, double dDELTA_deriv)
{
	return b
	       * (pow(delta_fn, (b - 1)) * d2DELTA_deriv + (b - 1) * pow(delta_fn, (b - 2)) * dDELTA_deriv * dDELTA_deriv);
	//   return b*(pow(delta_fn,(b-1))*d2DELTA_deriv+(b-1)*pow(delta_fn,(b-2))*pow(dDELTA_deriv,2));
}

double dphi_ddelta(double C, double delta, double phi)
{
	return -2 * C * (delta - 1) * phi;
}

double dphi_dtau(double D, double tau, double phi_fn)
{
	return -2 * D * (tau - 1) * phi_fn;
}

double d2phi_dtau2(double D, double tau, double phi_fn)
{
	return (2 * D * pow((tau - 1), 2) - 1) * 2 * D * phi_fn;
}

double d2phi_ddelta2(double C, double delta, double phi_fn)
{
	return (2 * C * (delta - 1) * (delta - 1) - 1) * 2 * C * phi_fn;
	//   return (2*C*pow((delta-1),2)-1)*2*C*phi_fn;
}

double dDELTA_dtau(double theta_fn, double b, double delta_fn)
{
	return -2 * theta_fn * b * pow(delta_fn, (b - 1));
}

double d2DELTA_dtau2(double b, double delta_fn, double theta_fn)
{
	return 2 * b * pow(delta_fn, (b - 1)) + 4 * theta_fn * theta_fn * b * (b - 1) * pow(delta_fn, (b - 2));
	//   return 2*b*pow(delta_fn,(b-1))+4*pow(theta_fn,2)*b*(b-1)*pow(delta_fn,(b-2));
}

double d2DELTApowb_ddeltadtau(
    double A, double b, double beta, double delta_fn, double delta, double theta_fn, double dDELTA_deriv)
{
	return -A * b * 2 / beta * pow(delta_fn, (b - 1)) * (delta - 1)
	           * pow((delta - 1) * (delta - 1), (1 / (2 * beta) - 1))
	       - 2 * theta_fn * b * (b - 1) * pow(delta_fn, (b - 2)) * dDELTA_deriv;

	//   return -A*b*2/beta*pow(delta_fn,(b-1))*(delta-1)*pow(pow((delta-1),2),(1/(2*beta)-1))
	//      -2*theta_fn*b*(b-1)*pow(delta_fn,(b-2))*dDELTA_deriv;
}

double d2phi_ddeltadtau(double C, double D, double delta, double tau, double phi_fn)
{
	return 4 * C * D * (delta - 1) * (tau - 1) * phi_fn;
}

/**********************************************************************
   A derivation of the free energy function phi
   last change: NB JUN 09
 ***********************************************************************/
double CFluidProperties::phi_r_d(double rho, double T) const
{
	double phi_a = 0, phi_b = 0, phi_c = 0, phi_d = 0;
	double delta, tau, DELTA, THETA, PHI, DPHI, dDELTA_deriv, dDELTApowbddelta;
	int i;

	tau = Tc / T;
	delta = rho / rhoc;

	for (i = 0; i < limit[3]; i++)
	{
		if (i < limit[0])
			phi_a = phi_a + (K[0][i] * K[1][i] * pow(delta, (K[1][i] - 1)) * pow(tau, K[2][i]));
		else
		{
			if (i < limit[1])
				phi_b = phi_b + (K[0][i] * exp(-pow(delta, K[3][i])) * (pow(delta, K[1][i] - 1) * pow(tau, K[2][i])
				                                                        * (K[1][i] - K[3][i] * pow(delta, K[3][i]))));
			else if (i < limit[2])
				phi_c = phi_c + (K[0][i] * pow(delta, K[1][i]) * pow(tau, K[2][i])
				                 * exp(-K[10][i] * (delta - K[13][i]) * (delta - K[13][i])
				                       - K[11][i] * (tau - K[12][i]) * (tau - K[12][i]))
				                 * (K[1][i] / delta - 2 * K[10][i] * (delta - K[13][i])));
			else if (i < limit[3])
			{
				THETA = theta_fn(tau, K[6][i], delta, K[11][i]);
				DELTA = delta_fn(THETA, K[7][i], delta, K[4][i]);
				PHI = phi_fn(K[8][i], delta, K[9][i], tau);
				dDELTA_deriv = dDELTA_ddelta(delta, K[6][i], THETA, K[11][i], K[7][i], K[4][i]);
				dDELTApowbddelta = dDELTApowb_ddelta(K[5][i], DELTA, dDELTA_deriv);
				DPHI = dphi_ddelta(K[8][i], delta, PHI);

				phi_d = phi_d + K[0][i] * (pow(DELTA, K[5][i]) * (PHI + delta * DPHI) + dDELTApowbddelta * delta * PHI);
			}
		}
	}
	return phi_a + phi_b + phi_c + phi_d;
}

/**********************************************************************
   A derivation of the free energy function phi
   last change: NB JUN 09
 ***********************************************************************/
double CFluidProperties::phi_r_tt(double rho, double T) const
{
	// CFluidProperties *FP;
	double phi_a = 0, phi_b = 0, phi_c = 0, phi_d = 0;
	double delta, tau, THETA, PHI, DELTA, DDELTA, D2DELTA, DPHI, D2PHI;
	int i = 0;

	// FP=MFPGet(c);

	tau = Tc / T;
	delta = rho / rhoc;

	for (i = 0; i < limit[3]; i++)
	{
		if (i < limit[0])
			phi_a = phi_a + (K[0][i] * K[2][i] * (K[2][i] - 1) * pow(delta, K[1][i]) * pow(tau, (K[2][i] - 2)));
		else if (i < limit[1])
			phi_b = phi_b + (K[0][i] * K[2][i] * (K[2][i] - 1) * pow(delta, K[1][i]) * pow(tau, (K[2][i] - 2))
			                 * exp(-pow(delta, K[3][i])));
		else if (i < limit[2])
			phi_c = phi_c + (K[0][i] * pow(delta, K[1][i]) * pow(tau, K[2][i])
			                 * exp(-K[10][i] * ((delta - K[13][i]) * (delta - K[13][i]))
			                       - K[11][i] * ((tau - K[12][i]) * (tau - K[12][i])))
			                 * (pow((K[2][i] / tau - 2 * K[11][i] * (tau - K[12][i])), 2) - K[2][i] / (tau * tau)
			                    - 2 * K[11][i]));
		else if (i < limit[3])
		{
			THETA = theta_fn(tau, K[6][i], delta, K[11][i]);
			DELTA = delta_fn(THETA, K[7][i], delta, K[4][i]);
			PHI = phi_fn(K[8][i], delta, K[9][i], tau);

			D2DELTA = d2DELTA_dtau2(K[5][i], DELTA, THETA);
			DDELTA = dDELTA_dtau(THETA, K[5][i], DELTA);
			DPHI = dphi_dtau(K[9][i], tau, PHI);
			D2PHI = d2phi_dtau2(K[9][i], tau, PHI);

			phi_d = phi_d + (K[0][i] * delta * (D2DELTA * PHI + 2 * DDELTA * DPHI + pow(DELTA, K[5][i]) * D2PHI));
		}
	}
	return phi_a + phi_b + phi_c + phi_d;
}

/**********************************************************************
   A derivation of the free energy function phi
   last change: NB JUN 09
 ***********************************************************************/
double CFluidProperties::phi_0_t(double T) const
{
	double phi_c = 0, phi_d = 0, phi_e = 0;
	double tau;
	int i;

	tau = Tc / T;

	phi_c = k[0][1];
	phi_d = k[0][2] / tau;

	for (i = 3; i < 8; i++)
		phi_e = phi_e + (k[0][i] * k[1][i] * (1 / (1 - exp(-k[1][i] * tau)) - 1));

	return phi_c + phi_d + phi_e;
}

/**********************************************************************
   A derivation of the free energy function phi
   last change: NB JUN 09
 ***********************************************************************/
double CFluidProperties::phi_0_tt(double T) const
{
	if (fluid_id == 3) // N2
	{
		double a[6];
		const double tau(Tc / T);

		a[0] = 2.5;
		a[1] = -1.934819e-4;
		a[2] = -1.247742e-5;
		a[3] = 6.678326e-8;
		a[4] = 1.012941;
		a[5] = 26.65788;

		const double exp_a8(exp(a[5] * tau));
		const double phi_zero_tt = (-a[0] + 2.0 * a[1] / tau + 6.0 * a[2] / tau / tau + 12 * a[3] / tau / tau / tau
		                            - a[4] * a[5] * a[5] * tau * tau * exp_a8 / ((exp_a8 - 1.0) * (exp_a8 - 1.0)))
		                           / (tau * tau);
		return phi_zero_tt;
	}
	else
	{
		double phi_d = 0, phi_e = 0;
		double tau;
		int i;

		tau = Tc / T;
		phi_d = k[0][2] / (tau * tau);
		for (i = 3; i < 8; i++)
			phi_e = phi_e + (k[0][i] * (k[1][i] * k[1][i]) * exp(-k[1][i] * tau) * pow(1 - exp(-k[1][i] * tau), -2));

		return 0 - phi_d - phi_e;
	}
}

/**********************************************************************
   A derivation of the free energy function phi
   last change: NB JUN 09
 ***********************************************************************/
double CFluidProperties::phi_r_t(double rho, double T) const
{
	double phi_a = 0, phi_b = 0, phi_c = 0, phi_d = 0, h;
	int i;
	double delta, tau;
	double thetafn, deltafn, ddeltatau, phifn, dphitau;

	tau = Tc / T;
	delta = rho / rhoc;

	for (i = 0; i < limit[0]; i++)
		phi_a = phi_a + (K[0][i] * K[2][i] * pow(delta, K[1][i]) * pow(tau, (K[2][i] - 1)));

	for (i = limit[0]; i < limit[1]; i++)
		phi_b = phi_b + (K[0][i] * K[2][i] * pow(delta, K[1][i]) * pow(tau, (K[2][i] - 1)) * exp(-pow(delta, K[3][i])));

	for (i = limit[1]; i < limit[2]; i++)
		phi_c = phi_c + (K[0][i] * pow(delta, K[1][i]) * pow(tau, (K[2][i]))
		                 * exp(-K[10][i] * ((delta - K[13][i]) * (delta - K[13][i]))
		                       - K[11][i] * (tau - K[12][i]) * (tau - K[12][i]))
		                 * (K[2][i] / tau - 2 * K[11][i] * (tau - K[12][i])));

	for (i = limit[2]; i < limit[3]; i++)
	{
		thetafn = theta_fn(tau, K[6][i], delta, K[11][i]);
		deltafn = delta_fn(thetafn, K[7][i], delta, K[4][i]);
		ddeltatau = dDELTA_dtau(thetafn, K[5][i], deltafn);
		phifn = phi_fn(K[8][i], delta, K[9][i], tau);
		dphitau = dphi_dtau(K[9][i], tau, phifn);

		phi_d = phi_d + (K[0][i] * delta * (ddeltatau * phifn + pow(deltafn, K[5][i]) * dphitau));
	}
	h = phi_a + phi_b + phi_c + phi_d;
	return h;
}

/**********************************************************************
   A derivation of the free energy function phi
   last change: NB 4.9.05
 ***********************************************************************/
double CFluidProperties::phi_r_dt(double rho, double T) const
{
	double phi_a = 0, phi_b = 0, phi_c = 0, phi_d = 0;
	int i;
	double delta, tau;
	double phifn, thetafn, deltafn, d2phideriv, dDELTAderiv, dDELTApowbddelta, dphidtau, dDELTApowbdtau, dphiddelta,
	    d2DELTApowbddeltadtau;

	tau = Tc / T;
	delta = rho / rhoc;

	for (i = 0; i < limit[0]; i++)
		phi_a = phi_a + (K[0][i] * K[1][i] * K[2][i] * pow(delta, (K[1][i] - 1)) * pow(tau, (K[2][i] - 1)));

	for (i = limit[0]; i < limit[1]; i++)
		phi_b = phi_b + (K[0][i] * K[2][i] * pow(delta, (K[1][i] - 1)) * pow(tau, (K[2][i] - 1))
		                 * (K[1][i] - K[3][i] * pow(delta, K[3][i]))
		                 * exp(-pow(delta, K[3][i])));

	for (i = limit[1]; i < limit[2]; i++)

		phi_c = phi_c + ((K[0][i] * pow(delta, K[1][i]) * pow(tau, K[2][i])
		                  * exp(-K[10][i] * ((delta - K[13][i]) * (delta - K[13][i]))
		                        - K[11][i] * ((tau - K[12][i]) * (tau - K[12][i]))))
		                 *

		                 (K[1][i] / delta - 2 * K[10][i] * (delta - K[13][i]))
		                 * (K[2][i] / tau - 2 * K[11][i] * (tau - K[12][i])));

	for (i = limit[2]; i < limit[3]; i++)
	{
		phifn = phi_fn(K[8][i], delta, K[9][i], tau);
		thetafn = theta_fn(tau, K[6][i], delta, K[11][i]);
		deltafn = delta_fn(thetafn, K[7][i], delta, K[4][i]);
		d2phideriv = d2phi_ddeltadtau(K[8][i], K[9][i], delta, tau, phifn);

		dDELTAderiv = dDELTA_ddelta(delta, K[6][i], thetafn, K[11][i], K[7][i], K[4][i]);

		dDELTApowbddelta = dDELTApowb_ddelta(K[5][i], deltafn, dDELTAderiv);
		dphidtau = dphi_dtau(K[9][i], tau, phifn);
		dDELTApowbdtau = dDELTA_dtau(thetafn, K[5][i], deltafn);
		dphiddelta = dphi_ddelta(K[8][i], delta, phifn);
		d2DELTApowbddeltadtau
		    = d2DELTApowb_ddeltadtau(K[6][i], K[5][i], K[11][i], deltafn, delta, thetafn, dDELTAderiv);

		phi_d = phi_d + (K[0][i] * (pow(deltafn, K[5][i]) * (dphidtau + delta * d2phideriv)
		                            + delta * dDELTApowbddelta * dphidtau
		                            + dDELTApowbdtau * (phifn + delta * dphiddelta)
		                            + d2DELTApowbddeltadtau * delta * phifn));
	}

	return phi_a + phi_b + phi_c + phi_d;
}

/**********************************************************************
   A derivation of the free energy function phi
   last change: NB 4.9.05
 ***********************************************************************/
double CFluidProperties::phi_r_dd(double rho, double T) const
{
	double phi_a = 0, phi_b = 0, phi_c = 0, phi_d = 0;
	int i;
	double delta, tau;
	double thetafn, deltafn, phifn, dphiddelta, d2phiddelta2, dDELTA_deriv, dDELTApowbddelta, d2DELTA_deriv,
	    d2DELTApowbddelta2;

	tau = Tc / T;
	delta = rho / rhoc;

	for (i = 0; i < limit[0]; i++)
		phi_a = phi_a + (K[0][i] * K[1][i] * (K[1][i] - 1) * pow(delta, (K[1][i] - 2)) * pow(tau, K[2][i]));
	for (i = limit[0]; i < limit[1]; i++)
		phi_b
		    = phi_b + ((K[0][i] * exp(-(pow(delta, K[3][i]))))
		               * (((pow(delta, (K[1][i] - 2)) * pow(tau, K[2][i])))
		                  * (((K[1][i] - K[3][i] * pow(delta, K[3][i]))) * (K[1][i] - 1 - K[3][i] * pow(delta, K[3][i]))
		                     - (pow(K[3][i], 2) * pow(delta, K[3][i])))));

	for (i = limit[1]; i < limit[2]; i++)
		phi_c = phi_c
		        + ((K[0][i] * pow(tau, K[2][i]))
		           * exp(-K[10][i] * pow((delta - K[13][i]), 2) - K[11][i] * ((tau - K[12][i]) * (tau - K[12][i])))
		           * ((-2 * K[10][i] * pow(delta, K[1][i])
		               + 4 * (K[10][i] * K[10][i]) * pow(delta, K[1][i]) * ((delta - K[13][i]) * (delta - K[13][i])))
		              + (-4 * K[1][i] * K[10][i] * pow(delta, (K[1][i] - 1)) * (delta - K[13][i])
		                 + K[1][i] * (K[1][i] - 1) * pow(delta, (K[1][i] - 2)))));

	for (i = limit[2]; i < limit[3]; i++)
	{
		thetafn = theta_fn(tau, K[6][i], delta, K[11][i]);
		deltafn = delta_fn(thetafn, K[7][i], delta, K[4][i]);
		phifn = phi_fn(K[8][i], delta, K[9][i], tau);
		dphiddelta = dphi_ddelta(K[8][i], delta, phifn);
		d2phiddelta2 = d2phi_ddelta2(K[8][i], delta, phifn);
		dDELTA_deriv = dDELTA_ddelta(delta, K[6][i], thetafn, K[11][i], K[7][i], K[4][i]);
		dDELTApowbddelta = dDELTApowb_ddelta(K[5][i], deltafn, dDELTA_deriv);
		dphiddelta = dphi_ddelta(K[8][i], delta, phifn);
		d2DELTA_deriv = d2DELTA_ddelta2(delta, dDELTA_deriv, K[6][i], thetafn, K[11][i], K[7][i], K[4][i]);
		d2DELTApowbddelta2 = d2DELTApowb_ddelta2(K[5][i], deltafn, d2DELTA_deriv, dDELTA_deriv);

		phi_d = phi_d + (K[0][i] * (pow(deltafn, K[5][i]) * (2 * dphiddelta + delta * d2phiddelta2)
		                            + 2 * dDELTApowbddelta * (phifn + delta * dphiddelta)
		                            + d2DELTApowbddelta2 * delta * phifn));
	}

	return phi_a + phi_b + phi_c + phi_d;
}

/**********************************************************************
   Function for calculating the Pressure of a gas/liquid density on density
   and temperature.

   Parameters:
         rho  - density
         rhoc - density at the critical point
         T    - temperature
         Tc   - critical temperature
         R    - specific gas constant

   Programming: NB
   Aug 2008
 ***********************************************************************/
double pressure(double rho, double T, int c)
{
	CFluidProperties const* const fluid_prop(MFPGet(c));
	const double rhoc(fluid_prop->getCriticalDensity());
	const double R(fluid_prop->getSpecificGasConstant());

	return (1 + (rho / rhoc) * fluid_prop->phi_r_d(rho, T)) * rho * R * T;
}

/**********************************************************************
   This function calculates the density depending on pressure and temperature
   by iteration. The iteration process may take a long time, so it's not re-
   comended to use this function for simulations. Better use the GetMatrixValue
   function and a table with rho-p-T data.
   This function does not work for every pressure/temperature Range!

   Parameters:
         P    - pressure
         rho0 - initial density for iteration
         rhoc - density at the critical point
   T    - temperature
   Tc   - critical temperature
   R    - specific gas constant
   prec - precision for iteration criteria

   Programming: NB
   Aug 2008
   last change: NB 4.9.05
 ***********************************************************************/
double density(double P, double rho0, double T, double prec, int c)
{
	CFluidProperties const* const fluid_prop(MFPGet(c));
	int iterations = 0;
	double rho = 0.0, p0; // OK411
	const double rhoc(fluid_prop->getCriticalDensity());
	const double R(fluid_prop->getSpecificGasConstant());

	p0 = 0;
	while (fabs(P - p0) > prec) // change to fabs. 15.09.2008 WW
	{
		rho = P / ((1 + (rho0 / rhoc) * fluid_prop->phi_r_d(rho0, T)) * R * T);
		p0 = pressure(rho0, T, c);
		rho0 = rho;
		iterations++;
		if (iterations > 50)
			return 0;
	}

	return rho;
}

/**********************************************************************
   Function for calculating enthalpy depending on density and
   Temperature.

   Parameters:
         rho  - density
         rhoc - density at the critical point
         T    - temperature
         Tc   - critical temperature
         R    - specific gas constant

   Programming: NB
   Aug 2008
 ***********************************************************************/
double enthalpy(double rho, double T, int c)
{
	CFluidProperties const* const fluid_prop(MFPGet(c));
	const double tau(fluid_prop->getCriticalTemperature() / T);
	const double delta(rho / fluid_prop->getCriticalDensity());
	const double R(fluid_prop->getSpecificGasConstant());

	double h = (1 + tau * (fluid_prop->phi_0_t(T) + fluid_prop->phi_r_t(rho, T)) + delta * fluid_prop->phi_r_d(rho, T))
	           * R * T;
	return h;
}

/**********************************************************************
   Function for calculating isobaric heat capacity depending on Temperature.

   linearised in a specified temperature interval for typical technological applications
 ***********************************************************************/

double linear_heat_capacity(double T, int c)
{
	double Tl[2], cpl[2]; // temperature limits in K, capacity limits in J/kg/K

	switch (c)
	{
		case 1: // H2O
			Tl[0] = 275.0;
			Tl[1] = 1000.0;
			cpl[0] = 1859.0;
			cpl[1] = 2288.0;
			break;
		case 3: // N2  1056.8+(1146.4-1056.8)/(900.0-500.0)*(variables[1]-500.0);
			Tl[0] = 275.0;
			Tl[1] = 1200.0;
			cpl[0] = 1039.0;
			cpl[1] = 1204.0;
			break;
		case 5: // O2
			Tl[0] = 275.0;
			Tl[1] = 1200.0;
			cpl[0] = 915.0;
			cpl[1] = 1115.0;
			break;
		default:
			std::cout << "WARNING: Fluid not specified in linear_heat_capacity. Setting cp to 1000.\n";
			return 1000.0;
	}

	return (cpl[0] + (cpl[1] - cpl[0]) / (Tl[1] - Tl[0]) * (T - Tl[0]));
}

/**********************************************************************
   Function for calculating isobaric heat capacity depending on Temperature.

   polynomials from Ostermayer for model comparison
 ***********************************************************************/

double polynomial_heat_capacity(double T, int c)
{
    //TODO: input file control: polynomial degree and coefficient array
    //followed by a general loop here
    double res(0.);
    switch (c) {
        case 1: //H2O
            res = 1.5630770E+03 + T*(1.6037550E+00 + T *(-2.9327840E-03 + T * (3.2161010E-06 - 1.1568270E-09*T)));
            break;
        case 3: //N2
            res = 9.7904300E+02 + T*(4.1796390E-01 + T *(-1.1762790E-03 + T * (1.6743940E-06 -7.2562970E-10*T)));
            break;
        default:
            std::cout << "WARNING: Fluid not specified in polynomial_heat_capacity. Setting cp to 1000.\n";
            return 1000.0;
    }

    return res;
}
/**********************************************************************
   Function for calculating isochoric heat capacity depending on density and
   Temperature.

   Parameters:
         rho  - density
         rhoc - density at the critical point
         T    - temperature
         Tc   - critical temperature
         R    - specific gas constant

   Programming: NB
   Aug 2008
   last change: NB 4.9.05
 ***********************************************************************/
double isochoric_heat_capacity(double rho, double T, int c)
{
	CFluidProperties const* const fluid_prop(MFPGet(c));
	//	thermal_properties (fluid, rhoc, Tc, R);
	const double h((fluid_prop->getCriticalTemperature() / T));
	double cv
	    = -(h * h * (fluid_prop->phi_0_tt(T) + fluid_prop->phi_r_tt(rho, T))) * fluid_prop->getSpecificGasConstant();

	return cv;
}

/**********************************************************************
   Function for calculating isobaric heat capacity depending on density and
   Temperature.

   Parameters:
         rho  - density
         rhoc - density at the critical point
         T    - temperature
         Tc   - critical temperature
         R    - specific gas constant

   Programming: NB
   Aug 2008
 ***********************************************************************/

double isobaric_heat_capacity(double rho, double T, int c)
{
	CFluidProperties const* const fluid_prop(MFPGet(c));
	const double tau(fluid_prop->getCriticalTemperature() / T);
	const double delta(rho / fluid_prop->getCriticalDensity());

	double cp = (-tau * tau * (fluid_prop->phi_0_tt(T) + fluid_prop->phi_r_tt(rho, T))
	             + (pow((1 + delta * fluid_prop->phi_r_d(rho, T) - delta * tau * fluid_prop->phi_r_dt(rho, T)), 2))
	                   / ((1 + 2 * delta * fluid_prop->phi_r_d(rho, T) + delta * delta * fluid_prop->phi_r_dd(rho, T))))
	            * fluid_prop->getSpecificGasConstant();

	return cp;
}

/**********************************************************************
   Function for calculating viscosity of CO2 depending on density and Temperature.
   Programming: NB
          Aug 2008
 ***********************************************************************/
double co2_viscosity(double rho, double T)
{
	double eta, eta_0, d_eta;
	double psi, t_r;
	double a[5], b[8], d[8][5];
	int i, j;

	// coefficients of the representation of the zero-density viscosity of co2
	a[0] = 0.235156;
	a[1] = -0.491266;
	a[2] = 0.05211155;
	a[3] = 0.05347906;
	a[4] = -0.01537102;

	psi = 0; // variable for reduced effective cross-section

	// coefficients of the representation of the excess viscosity of co2

	for (i = 0; i < 8; i++)
		for (j = 0; j < 5; j++)
		{
			b[i] = 0;
			d[i][j] = 0;
		}

	d[0][0] = 0.4071119e-02;
	d[1][0] = 0.7198037e-04;
	d[5][3] = 0.2411697e-16;
	d[7][0] = 0.2971072e-22;
	d[7][1] = -0.1627888e-22;

	t_r = T / 251.196; // reduced temperature

	// deriving the zero-density viscosity eta_0(T)

	for (i = 0; i < 5; i++)
		psi = psi + a[i] * MathLib::fastpow(log(t_r), (i));
	psi = exp(psi);

	// TF   eta_0 = 1.00697 * pow (T,0.5) / psi;
	eta_0 = 1.00697 * sqrt(T) / psi;
	d_eta = 0;

	// deriving the excess viscosity d_eta(rho,T)

	for (i = 0; i < 8; i++)
	{
		for (j = 0; j < 5; j++)
			b[i] = b[i] + d[i][j] / MathLib::fastpow(t_r, j);
		d_eta = d_eta + b[i] * MathLib::fastpow(rho, i + 1);
	}

	// deriving dynamic viscosity as sum of eta_o(T) and d_eta(rho,T)

	eta = (eta_0 + d_eta) * 1e-06; // eta in [Pa s]

	return eta;
}

/**********************************************************************
   Function for calculating heat conductivity of CO2 depending on density and Temperature.
   (Vesovic&Wakeham)
   Programming: NB 4.8.01
          Nov 2008
 ***********************************************************************/
double co2_heat_conductivity(double rho, double T)
{
	double Tc = 304.1282;
	// double    pc = 7377300;
	double rho_c = 467.6;
	// double     M = 0.044098;
	// double     R = 8.314472;
	// double    NA = 6.0221353e23;
	double rho_r = rho / rho_c;
	double Tr = T / Tc;
	double lamda_r;
	double lamda;

	double a[13] = {0, 3, 6.70697, 0.94604, 0.3, 0.3, 0.39751, 0.33791, 0.77963, 0.79857, 0.9, 0.02, 0.2};
	double g[11] = {0, 0, 0, 1.5, 0, 1, 1.5, 1.5, 1.5, 3.5, 5.5};
	double h[11] = {0, 1, 5, 1, 1, 2, 0, 5, 9, 0, 0};
	double n[11] = {0,
	                7.69857587E+00,
	                1.59885811E-01,
	                1.56918621E+00,
	                -6.73400790E+00,
	                1.63890156E+01,
	                3.69415242E+00,
	                2.23205514E+01,
	                6.61420950E+01,
	                -1.71779133E-01,
	                4.33043347E-03};
	double nc = 0.775547504;
	double lamda_r_ce;
	double temp1, temp2, temp3;

	double var = 1.0 + a[11] * pow(pow(1 - Tr, 2), a[12]);
	double alpha = 1 - a[10] * log(var + sqrt(var * var - 1.0));
	temp1 = rho_r * exp(-pow(rho_r, a[1]) / a[1] - pow(a[2] * (Tr - 1), 2) - pow(a[3] * (rho_r - 1), 2));
	temp2 = pow(pow((1 - 1 / Tr) + a[4] * pow(pow(rho_r - 1, 2), 1. / (2. * a[5])), 2), a[6]);
	temp3 = pow(pow(a[7] * (rho_r - alpha), 2), a[8]);
	lamda_r_ce = temp1 / pow(temp2 + temp3, a[9]);
	temp1 = 0;
	temp2 = 0;

	for (int i = 1; i < 4; i++)
	{
		temp1 += n[i] * pow(Tr, g[i]) * pow(rho_r, h[i]);
	}
	for (int i = 4; i < 11; i++)
	{
		temp2 += n[i] * pow(Tr, g[i]) * pow(rho_r, h[i]);
	}

	temp2 *= exp(-5. * rho_r * rho_r);

	lamda_r = temp1 + temp2 + nc * lamda_r_ce;

	//    double Lamda_C=pow(R,5./6.)*pow(pc,2./3.)/pow(Tc,1./6.)/pow(M,0.5)/pow(NA,1./3.);
	double Lamda_C = 0.00481384; // W/m/K

	lamda = lamda_r * Lamda_C;

	return lamda;
}

/**********************************************************************
   Function for calculating viscosity of CH4 at 295K depending on pressure.
   (Gulik,Mostert,Van den Berg)
   Programming: NB 4.8.01
          Nov 2008
 ***********************************************************************/
double ch4_viscosity_295K(double p)
{
	double h;

	p = p / 100000;
	h = (-3.7091411E-14 * MathLib::fastpow(p, 4) + 9.1937114E-10 * MathLib::fastpow(p, 3) - 6.6099446E-06 * p * p
	     + 4.8400147E-02 * p
	     + 1.0934694E+01)
	    / 1.0e6;

	return h;
}

/**********************************************************************
   Function for calculating viscosity of pure water at a density rho and
   a temperature T.
   Programming: NB
          Apr 2009
 ***********************************************************************/
double h2o_viscosity_IAPWS(double rho, double T)
{
	double my, my_0, my_1;
	double H[4], h[6][7];
	double sum1 = 0, sum2 = 0, sum3 = 0;
	int i, j;

	T = T / 647.096;
	rho = rho / 322.0;

	H[0] = 1.67752;
	H[1] = 2.20462;
	H[2] = 0.6366564;
	H[3] = -0.241605;

	for (i = 0; i < 6; i++)
		for (j = 0; j < 7; j++)
			h[i][j] = 0;
	h[0][0] = 0.520094000;
	h[1][0] = 0.085089500;
	h[2][0] = -1.083740000;
	h[3][0] = -0.289555000;
	h[0][1] = 0.222531000;
	h[1][1] = 0.999115000;
	h[2][1] = 1.887970000;
	h[3][1] = 1.266130000;
	h[5][1] = 0.120573000;
	h[0][2] = -0.281378000;
	h[1][2] = -0.906851000;
	h[2][2] = -0.772479000;
	h[3][2] = -0.489837000;
	h[4][2] = -0.257040000;
	h[0][3] = 0.161913000;
	h[1][3] = 0.257399000;
	h[0][4] = -0.032537200;
	h[3][4] = 0.069845200;
	h[4][5] = 0.008721020;
	h[3][6] = -0.004356730;
	h[5][6] = -0.000593264;

	for (i = 0; i < 4; i++)
		sum1 = sum1 + (H[i] / (MathLib::fastpow(T, i)));

	my_0 = 100 * sqrt(T) / sum1;

	for (i = 0; i < 6; i++)
	{
		for (j = 0; j < 7; j++)
			sum3 = sum3 + h[i][j] * MathLib::fastpow(rho - 1, j);
		sum2 = sum2 + (MathLib::fastpow(1 / T - 1, i) * sum3);
		sum3 = 0;
	}

	my_1 = exp(rho * sum2);

	my = (my_0 * my_1) / 1e6;
	return my;
}

/**********************************************************************
   Function for calculating viscosity of oxygen at a density rho and
   a temperature T.
   Programming: NB/TN
   see Lemmon & Jacobson, Int J of Thermophys, 25(1), 2004
 ***********************************************************************/
double o2_viscosity(double rho, double T)
{
	const double MM(32.0); //
	const double sigma_squared(0.11751184); // nm^2
	const double ek(118.5); // K
	const double T_c(154.581); // K
	const double rho_c(13.63);

	rho /= MM; // convert density in mol/dm^3

	double b[5];

	b[0] = 0.431;
	b[1] = -0.4623;
	b[2] = 0.08406;
	b[3] = 0.005341;
	b[4] = -0.00331;

	double exponent(0);

	for (unsigned i = 0; i < 5; i++)
		exponent += b[i] * MathLib::fastpow(log(T / ek), i);

	const double Omega = exp(exponent);
	const double eta_0 = 0.0266958 * pow(MM * T, 0.5) / sigma_squared / Omega;
	double rho_min(0.2); // densities lower than that do not influence eta_r (less than 1%)
	if (rho < rho_min)
		return eta_0 * 1e-6; // for better efficiency near dilute gas (eta in Pa*s)

	double N[5];
	double t[5];
	int d[5];
	int l[5];
	int gamma[5];

	N[0] = 17.67;
	N[1] = 0.4042;
	N[2] = 0.0001077;
	N[3] = 0.3510;
	N[4] = -13.67;

	t[0] = 0.05;
	t[1] = 0.0;
	t[2] = 2.1;
	t[3] = 0.0;
	t[4] = 0.5;

	d[0] = 1;
	d[1] = 5;
	d[2] = 12;
	d[3] = 8;
	d[4] = 1;

	l[0] = 0;
	l[1] = 0;
	l[2] = 0;
	l[3] = 1;
	l[4] = 2;

	gamma[0] = 0;
	gamma[1] = 0;
	gamma[2] = 0;
	gamma[3] = 1;
	gamma[4] = 1;

	const double tau = T_c / T;
	const double delta = rho / rho_c;
	double eta_r(0);

	for (unsigned i = 0; i < 5; i++)
	{
		eta_r += N[i] * pow(tau, t[i]) * MathLib::fastpow(delta, d[i]) * exp(-gamma[i] * MathLib::fastpow(delta, l[i]));
	}

	return (eta_0 + eta_r) * 1e-6; // returns viscosity in Pa*s
}

/**********************************************************************
   Viscosity for different Fluids

   Programming: NB 4.8.01
          Nov 2008
 ***********************************************************************/
double Fluid_Viscosity(double rho, double T, double p, int fluid)
{
	p = p; // OK411
	// TODO: make a global function for all properties: Fluid_Property(int c, int property, double T, double P) (NB)

	double h;

	switch (fluid)
	{
		case 0: // CARBON_DIOXIDE
			h = co2_viscosity(rho, T);
			break;
		case 1: // WATER
			//h = h2o_viscosity_IAPWS(rho, T);
			{
				return MaterialLib::Fluid::WaterViscosityIAPWS::getValue(T, rho);
			}
			break;
		case 2: // METHANE
			// h = ch4_viscosity_295K(p);
			h = ch4_viscosity(rho, T);
			break;

		case 3: // Nitrogen
			h = n2_viscosity(rho, T);
			break;
		case 5:
			h = o2_viscosity(rho, T); // Oxygen
			break;
		default:
			h = 1E-3;
	}

	return h;
}

/**********************************************************************
   Heat conductivity for different Fluids

   Programming: NB 4.8.01
          Nov 2008
 ***********************************************************************/
double Fluid_Heat_Conductivity(double rho, double T, int fluid)
{
	double h;

	switch (fluid)
	{
		case 0: // CARBON_DIOXIDE
			h = co2_heat_conductivity(rho, T);
			break;
		case 1: // WATER
			// h = 0.598; // [W/m/K] at 293K and 1 bar
			h = h2o_heat_conductivity_IAPWS_ind(rho, T);
			break;
		case 2: // METHANE
			h = ch4_heat_conductivity(rho, T);
			break;
		case 3: // NITROGEN
			h = n2_heat_conductivity(rho, T);
			break;
		case 4: // NITROGEN
			h = n2_heat_conductivity(rho, T);
			break;
		case 5:
			h = o2_heat_conductivity(rho, T); // OXYGEN
			break;

		default:
			h = 0.5;
	}

	// if(caption.find("CARBON_DIOXIDE")!=string::npos)
	//
	// if(caption.find("METHANE")!=string::npos)
	//      h = 0.0338; // [W/m/K] at 298K and 1 bar
	// if(caption.find("WATER")!=string::npos)
	//      h = 0.598; // [W/m/K] at 293K and 1 bar

	return h;
}

/*************************************************************
* vapour pressure for co2
* NB, Dec 08
   last change: NB 4.9.05
*************************************************************/
double vapour_pressure_co2(double T)
{
	double Tc = 304.128;
	double pc = 7377300;
	double p = 0;
	double a[4], t[4];
	int i;
	a[0] = -7.0602087;
	t[0] = 1;
	a[1] = 1.9391218;
	t[1] = 1.5;
	a[2] = -1.6463597;
	t[2] = 2;
	a[3] = -3.2995634;
	t[3] = 4;

	for (i = 0; i < 4; i++)
		p = p + a[i] * pow(1 - (T / Tc), t[i]);

	p = exp(Tc / T * p) * pc;

	return p;
}

/*************************************************************
* sublime pressure for co2
* NB, Dec 08
   last change: NB 4.9.05
*************************************************************/
double sublime_pressure_co2(double T, double Tt, double pt)
{
	double p;
	double a[3];

	a[0] = -14.740846;
	a[1] = 2.4327015;
	a[2] = -5.3061778;

	p = exp(Tt / T * (a[0] * (1 - T / Tt) + a[1] * pow((1 - T / Tt), 1.9) + a[2] * pow((1 - T / Tt), 2.9))) * pt;

	return p;
}

/*************************************************************
* melting pressure for co2
* NB, Dec 08
   last change: NB 4.9.05
*************************************************************/
// just for CO2
double melting_pressure_co2(double T, double Tt, double pt)
{
	double p;
	double a[2];
	a[0] = 1955.539; // CO2
	a[1] = 2055.4593;

	p = (1 + a[0] * (T / Tt - 1) + a[1] * (T / Tt - 1) * (T / Tt - 1)) * pt;

	return p;
}

/*************************************************************
 * Peng&Robinson Equation of State
 * Analytical solving of third grade polynomial
 *
 * Parameters: temperature and pressure
 *             caption defines fluid
 * Programming: NB, Dec 08
 **************************************************************/
double preos(const CFluidProperties* mfp, double T, double P)
{
	double z1, z2, z3, h;
	vector<double> roots;

	// WW CFluidProperties* mfp_prop;

	double a, b, Tc, pc, MM, Ru;
	double omega, alpha;

	// int i;
	// mfp_prop = MFPGet (c); //14.11.2012. WW

	Ru = mfp->getUniversalGasConstant();
	MM = mfp->getMolarMass();
	Tc = mfp->getCriticalTemperature();
	pc = mfp->getCriticalPressure() / 1000;
	omega = mfp->getAzentricFactor(); // azentric factor

	// Peng Robinson EOS:
	// P= R*T / (V-b) - a*alpha / (V^2+2bV-b^2)   where V = MM/rho

	a = 0.457235 * Ru * Ru * Tc * Tc / pc;
	b = 0.077796 * Ru * Tc / pc;
	P = P / 1000; // P in kPa
	alpha = ((1 + (0.37464 + 1.5422 * omega - 0.26992 * omega * omega) * (1 - sqrt(T / Tc))))
	        * ((1 + (0.37464 + 1.5422 * omega - 0.26992 * omega * omega) * (1 - sqrt(T / Tc))));

	// EOS in the form: 0 = rho^3 + z1*rho^2 + z2*rho + z3

	z1 = (MM * a * alpha - 3 * MM * pow(b, 2) * P - 2 * MM * Ru * T * b)
	     / (b * (P * pow(b, 2) + b * Ru * T - a * alpha));
	z2 = (pow(MM, 2) * (b * P - Ru * T)) / (b * (P * pow(b, 2) + b * Ru * T - a * alpha));
	z3 = (MathLib::fastpow(MM, 3) * P) / (b * (P * b * b + b * Ru * T - a * alpha));

	NsPol3(z1, z2, z3, &roots); // derives the roots of the polynomial

	h = FindMin(roots);
	return h; // returns the lowest positive root
}

/////****************************************************************************
////* Finds and returns the positive minimum of a vector.
////* Programming: NB Dec 08
////*****************************************************************************/
////double FindMin (vector<double>Vec)
////{
////double x=DBL_MAX;
////int unsigned i;
////
////for(i=0;i<Vec.size();i++) {if ((Vec[i]>=0)&&(Vec[i]<x)) x=Vec[i];}
////
////return x;
////
////}

/*************************************************************
 * Redlich&Kwong Equation of State
 * Analytical solving of third grade polynomial
 *
 * Parameters: temperature and pressure
 *             caption defines fluid
 * Programming: NB, Jan 09
 **************************************************************/
double rkeos(double T, double P, int c)
{
	double z1, z2, z3, h;
	vector<double> roots;

	CFluidProperties const* const mfp_prop(MFPGet(c));

	double a, b, Tc, pc, MM;
	double Ru;

	if (P < 0)
		P = 100000; // set pressure to 1atm if unstable NB

	Ru = mfp_prop->getUniversalGasConstant() * 10; // universal gas constant [bar cm3/mol/K]
	MM = mfp_prop->getMolarMass();
	Tc = mfp_prop->getCriticalTemperature();
	pc = mfp_prop->getCriticalPressure() / 100000; // critical pressure

	// Redlich-Kwong EOS:
	// P= R*T (1+y+y^2-y^3)/ v(1-y^3) - a / (T^0.5*v(cv+b)   where V = MM/rho and y = b / (4v)

	a = 27 * Ru * Ru * pow(Tc, 2.5) / (64 * pc);
	b = 0.0866 * Ru * Tc / pc;
	P = P / 100000; // P in bar

	// EOS in the form: 0 = vm^3 + z1*vm^2 + z2*vm + z3

	z1 = -Ru * T / P;
	z2 = -(Ru * T * b / P - a / (sqrt(T) * P) + b * b);
	z3 = -a * b / (sqrt(T) * P);

	NsPol3(z1, z2, z3, &roots); // derives the roots of the polynomial

	h = FindMax(roots); // returns the lowest positive root (molar volume)
	h = MM / h * 1000; // density in kg/m3
	return h;
}

/*************************************************************
 * Redlich&Kwong Equation of State
 * Analytical solving of third grade polynomial
 *
 * This form of RKEOS requires adjusting parameters a and b as
   input values, also the molecular mass

 * Parameters: temperature and pressure
 *             adj. parameters a and b
 *             molecular mass

 * Programming: NB, May 09
 **************************************************************/
double rkeos(double T, double P, double MM, double a, double b)
{
	double z1, z2, z3, h;
	vector<double> roots;

	double Ru;
	P = P / 100000; // P in bar
	Ru = 83.14472; // universal gas constant [bar cm3/mol/K]
	// Redlich-Kwong EOS:
	// 0 = vm^3 + z1*vm^2 + z2*vm + z3

	z1 = -Ru * T / P;
	z2 = -(Ru * T * b / P - a / (sqrt(T) * P) + b * b);
	z3 = -a * b / (sqrt(T) * P);

	NsPol3(z1, z2, z3, &roots); // derives the roots of the polynomial

	h = FindMax(roots); // returns the lowest positive root (molar volume)
	h = MM / h * 1000; // density in kg/m3
	return h;
}

/**********************************************************************
   Function for calculating thermal conductivity of pure water at a density
   rho and a temperature T.
   IAPWS Formulation 1997 for industrial use
   Programming: NB
          Apr 2009
 ***********************************************************************/
double h2o_heat_conductivity_IAPWS_ind(double rho, double T)
{
	double lamda, lamda_0, lamda_1, lamda_2;
	double sum1 = 0;
	double S, Q, dT;
	double a[4], b[3], B[2], d[4], C[6];
	int i;

	T = T / 647.096;
	rho = rho / 317.11;

	a[0] = 0.0102811;
	a[1] = 0.0299621;
	a[2] = 0.0156146;
	a[3] = -0.00422464;

	b[0] = -0.397070;
	b[1] = 0.400302;
	b[2] = 1.060000;

	B[0] = -0.171587;
	B[1] = 2.392190;

	d[0] = 0.0701309;
	d[1] = 0.0118520;
	d[2] = 0.00169937;
	d[3] = -1.0200;

	C[0] = 0.642857;
	C[1] = -4.11717;
	C[2] = -6.17937;
	C[3] = 0.00308976;
	C[4] = 0.0822994;
	C[5] = 10.0932;

	for (i = 0; i < 4; i++)
		sum1 = sum1 + a[i] * MathLib::fastpow(T, i);

	lamda_0 = sqrt(T) * sum1;
	lamda_1 = b[0] + b[1] * rho + b[2] * exp(B[0] * (rho + B[1]) * (rho + B[1]));

	dT = fabs(T - 1) + C[3];
	Q = 2 + (C[4] / pow(dT, 3. / 5.));

	if (T >= 1)
		S = 1 / dT;
	else
		S = C[5] / pow(dT, 3. / 5.);

	lamda_2 = (d[0] / MathLib::fastpow(T, 10) + d[1]) * pow(rho, 9. / 5.) * exp(C[0] * (1 - pow(rho, 14. / 5.)))
	          + d[2] * S * pow(rho, Q) * exp((Q / (1. + Q)) * (1 - pow(rho, (1. + Q))))
	          + d[3] * exp(C[1] * pow(T, 3. / 2.) + C[2] / MathLib::fastpow(rho, 5));

	lamda = (lamda_0 + lamda_1 + lamda_2); // lamda in [W/m/K]

	return lamda;
}

/**********************************************************************
   Function for calculating viscosity of methane (CH4) depending on density
   and Temperature.

   see
   Friend, Ely, Ingham: The Transport Properties of Methane,
   J. Chem. Phys. Ref. Data,1989.

   Programming: NB
          Apr 2009
 ***********************************************************************/
double ch4_viscosity(double rho, double T)
{
	double eta, eta_0, eta_ex;
	double t, delta, tau, Omega = 0;
	double C[9], s[11], g[11];
	size_t r[11]; // TF
	double sum1 = 0, sum2 = 0;
	int unsigned i;

	rho = rho / 16.043; // rho in [mol/dm3]
	t = T / 174.;

	C[0] = -3.0328138281;
	C[1] = 16.918880086;
	C[2] = -37.189364917;
	C[3] = 41.288861858;
	C[4] = -24.615921140;
	C[5] = 8.9488430959;
	C[6] = -1.8739240542;
	C[7] = 0.20966101390;
	C[8] = -9.6570437074e-03;

	r[0] = 1;
	r[1] = 1;
	r[2] = 2;
	r[3] = 2;
	r[4] = 2;
	r[5] = 3;
	r[6] = 3;
	r[7] = 4;
	r[8] = 4;
	r[9] = 1;
	r[10] = 1;
	s[0] = 0;
	s[1] = 1;
	s[2] = 0;
	s[3] = 1;
	s[4] = 1.5;
	s[5] = 0;
	s[6] = 2;
	s[7] = 0;
	s[8] = 1;
	s[9] = 0;
	s[10] = 1;

	g[0] = 0.41250137;
	g[1] = -0.14390912;
	g[2] = 0.10366993;
	g[3] = 0.40287464;
	g[4] = -0.24903524;
	g[5] = -0.12953131;
	g[6] = 0.06575776;
	g[7] = 0.02566628;
	g[8] = -0.03716526;
	g[9] = -0.38798341;
	g[10] = 0.03533815;

	for (i = 0; i < 9; i++)
		Omega = Omega + C[i] * pow(t, (i / 3. - 1));
	Omega = 1 / Omega;

	eta_0 = 10.5 * sqrt(t) / Omega;

	delta = rho / 10.139;
	tau = 190.551 / T;

	for (i = 0; i < 9; i++)
		sum1 = sum1 + g[i] * MathLib::fastpow(delta, r[i]) * pow(tau, s[i]);

	for (i = 9; i < 11; i++)
		sum2 = sum2 + g[i] * MathLib::fastpow(delta, r[i]) * pow(tau, s[i]);

	eta_ex = 12.149 * sum1 * pow(1 + sum2, -1.0);

	eta = (eta_0 + eta_ex) * 1e-6;
	return eta;
}

/**********************************************************************
   Function for calculating viscosity of methane (CH4) depending on density
   and Temperature.

   see
   Friend, Ely, Ingham: The Transport Properties of Methane,
   J. Chem. Phys. Ref. Data,1989.

   Programming: NB
          Apr 2009
 ***********************************************************************/
double ch4_heat_conductivity(double rho, double T)
{
	double lamda, eta_0, lamda_0, lamda_ex;
	double t, delta, tau, Omega = 0, epsilon, beta;
	double C[9], Q[6], j[7], f[2], H[6], J[5];
	size_t r[7], s[7];
	double sum = 0, phi_id_tt, f_int;
	double P_sigma, Z_c, delta_star, T_star, rho_sigma_v;

	int unsigned i;
	const double T_c = 190.551;
	const double rho_c = 10.139;
	const double P_c = 4.599200;
	const double R = 8.314510; // [J/mol/K]

	t = T / 174.;
	tau = T_c / T;
	T_star = (T_c - T) / T_c;
	if (T_star < 0)
		T_star = 0;
	rho = rho / 16.043; // rho in [mol/dm3]
	delta = rho / rho_c;

	C[0] = -3.0328138281;
	C[1] = 16.918880086;
	C[2] = -37.189364917;
	C[3] = 41.288861858;
	C[4] = -24.615921140;
	C[5] = 8.9488430959;
	C[6] = -1.8739240542;
	C[7] = 0.20966101390;
	C[8] = -9.6570437074e-03;
	Q[0] = 2.5998324;
	Q[1] = -3.3854083;
	Q[2] = 1.6900979;
	Q[3] = -0.3911541;
	Q[4] = 4.7206715;
	Q[5] = -10.543907;

	f[0] = 1.458850;
	f[1] = -0.4377162;
	r[0] = 1;
	r[1] = 3;
	r[2] = 4;
	r[3] = 4;
	r[4] = 5;
	r[5] = 5;
	r[6] = 2;
	s[0] = 0;
	s[1] = 0;
	s[2] = 0;
	s[3] = 1;
	s[4] = 0;
	s[5] = 1;
	s[6] = 0;

	j[0] = 2.414927;
	j[1] = 0.55166331;
	j[2] = -0.52837734;
	j[3] = 0.073809553;
	j[4] = 0.24465507;
	j[5] = -0.047613626;
	j[6] = 1.5554612;
	H[0] = -6.589879;
	H[1] = 0.6355175;
	H[2] = 11.31028;
	H[3] = -10.38720;
	H[4] = 3.393075;
	J[0] = -0.7377483;
	J[1] = -1.241532;
	J[2] = -1.649972;
	J[3] = 2.281949;
	J[4] = 1.439570;

	beta = 0.355;
	epsilon = 1.9;

	for (i = 0; i < 9; i++)
		Omega = Omega + C[i] * pow(t, (i / 3. - 1));
	Omega = 1 / Omega;

	eta_0 = 10.5 * sqrt(t) / Omega; // (Eq. 10a)

	phi_id_tt = -Q[0] + 4 * Q[1] / 9 * pow(tau, -1. / 3.) + 10 * Q[2] / 9 * pow(tau, -2. / 3.) + 2 * Q[3] * pow(tau, -1)
	            - Q[4] * Q[5] * Q[5] * tau * tau * exp(Q[5] * tau) * pow((exp(Q[5] * tau) - 1), -2);

	f_int = f[0] + (f[1] / t);

	//(Eq. 14)
	lamda_0 = 0.51826 * eta_0 * (3.75 - f_int * (phi_id_tt + 1.5));

	//(Eq. 3)
	P_sigma
	    = P_c * exp(H[0] * T_star / (1 - T_star) + H[1] * T_star + H[2] * pow(T_star, epsilon) + H[3] * pow(T_star, 2)
	                + H[4] * pow(T_star, 3));
	Z_c = P_c / (R * T_c * rho_c);

	rho_sigma_v = P_sigma / (R * T) * (1
	                                   + P_sigma * MathLib::fastpow(tau, 8) * (Z_c - 1) / P_c
	                                         * (1
	                                            + (J[0] * pow(T_star, beta) + J[1] * pow(T_star, 2 * beta)
	                                               + J[2] * (T_star + MathLib::fastpow(T_star, 4))
	                                               + J[3] * MathLib::fastpow(T_star, 2))
	                                                  / (1 + J[4] * T_star))); // (Eq. 5)

	if ((T < T_c) && (rho < rho_c)) // (Eq. 16)
		delta_star = rho_sigma_v / rho_c;
	else
		delta_star = 11;

	for (i = 0; i < 7; i++) // (Eq. 17)

		sum = sum + j[i] * MathLib::fastpow(delta, r[i]) * MathLib::fastpow(tau, s[i]);

	// (Eq. 17)
	lamda_ex = 6.29638 * (sum + j[6] * delta * delta / delta_star);

	lamda = (lamda_0 + lamda_ex) / 1000; // lamda in [W/m/K]

	return lamda;
}

/**********************************************************************
   Function for calculating viscosity of nitrogen (N2) depending on density
   and Temperature.

   see
   Stephan, Krauss and Laeseke: Viscosity and themal conductivity of fluid
   Nitrogen, J. Chem. Phys. Ref. Data,Vol. 16, No. 4, 1987.

   Programming: NB
          Apr 2009
 ***********************************************************************/
double n2_viscosity(double rho, double T)
{
	// OK411 const double M = 28.013;
	// OK411 const double T_c = 126.2; //[K]
	// OK411 const double P_c = 3.4; // [MPa]
	const double rho_c = 314; // [kg/m3]
	// OK411 const double T_t = 63.1; //[K]
	const double CVF = 14.058; // [1e-3 Pa-s]

	const double sigma = 0.36502496e-09;
	const double k = 1.38062e-23;
	const double eps = 138.08483e-23;
	// OK411 const double N_A = 6.02213E26;
	// OK411 const double pi = 3.14159;
	const double c1 = 0.3125;
	const double c2 = 2.0442e-49;

	double A[5], C[5];

	double sum = 0, eta, eta_0, eta_r, T_star, Omega = 0;
	int i;

	A[0] = 0.46649;
	A[1] = -0.57015;
	A[2] = 0.19164;
	A[3] = -0.03708;
	A[4] = 0.00241;

	C[0] = -20.09997;
	C[1] = 3.4376416;
	C[2] = -1.4470051;
	C[3] = -0.027766561;
	C[4] = -0.21662362;

	T_star = T * k / eps;
	rho = rho / rho_c;

	for (i = 0; i < 5; i++)
		Omega = Omega + A[i] * MathLib::fastpow(log(T_star), i);

	Omega = exp(Omega);

	// eta in [Pa*s]
	eta_0 = c1 * sqrt(c2 * T) / (sigma * sigma * Omega);

	for (i = 2; i < 5; i++)
		sum = sum + C[i] * MathLib::fastpow(rho, i - 1);

	//
	eta_r = CVF * 1e-6 * (C[0] / (rho - C[1]) + C[0] / C[1] + sum);

	eta = (eta_0 + eta_r); // [Pa*s]

	return eta;
}

/**********************************************************************
   Function for calculating thermal conductivity of oxygen (O2) depending
   on density and Temperature.

   see
   Stephan, Krauss and Laeseke: Viscosity and themal conductivity of fluid
   Nitrogen, J. Chem. Phys. Ref. Data,Vol. 16, No. 4, 1987.

   Programming: TN
   see Lemmon & Jacobson, Int J of Thermophys, 25(1), 2004
 ***********************************************************************/
double o2_heat_conductivity(double rho, double T)
{
	const double MM(32.0); //
	const double sigma_squared(0.11751184); // nm^2
	const double ek(118.5); // K
	const double T_c(154.581); // K
	const double rho_c(13.63);

	rho /= MM; // convert density in mol/dm^3

	const double tau = T_c / T;
	const double delta = rho / rho_c;

	double b[5];

	b[0] = 0.431;
	b[1] = -0.4623;
	b[2] = 0.08406;
	b[3] = 0.005341;
	b[4] = -0.00331;

	double exponent(0);

	for (unsigned i = 0; i < 5; i++)
		exponent += b[i] * MathLib::fastpow(log(T / ek), i);

	const double Omega = exp(exponent);
	const double eta_0 = 0.0266958 * pow(MM * T, 0.5) / sigma_squared / Omega;

	double N[9], t[8];
	N[0] = 1.036;
	N[1] = 6.283;
	N[2] = -4.262;
	N[3] = 15.31;
	N[4] = 8.898;
	N[5] = -0.7336;
	N[6] = 6.728;
	N[7] = -4.374;
	N[8] = -0.4747;

	t[0] = -0.9;
	t[1] = -0.6;
	t[2] = 0.0;
	t[3] = 0.0;
	t[4] = 0.3;
	t[5] = 4.3;
	t[6] = 0.5;
	t[7] = 1.8;

	const double lam_0 = N[0] * eta_0 + N[1] * pow(tau, t[0])
	                     + N[2] * pow(tau, t[1]); // note counting is different due to different array sizes

	int d[6], l[6], gamma[6];

	d[0] = 1;
	d[1] = 3;
	d[2] = 4;
	d[3] = 5;
	d[4] = 7;
	d[5] = 10;

	l[0] = 0;
	l[1] = 0;
	l[2] = 0;
	l[3] = 2;
	l[4] = 2;
	l[5] = 2;

	gamma[0] = 0;
	gamma[1] = 0;
	gamma[2] = 0;
	gamma[3] = 1;
	gamma[4] = 1;
	gamma[5] = 1;

	double lam_r(0);

	for (unsigned i = 0; i < 6; i++)
	{
		lam_r += N[i + 3] * pow(tau, t[i + 2]) * MathLib::fastpow(delta, d[i])
		         * exp(-gamma[i] * MathLib::fastpow(delta, l[i]));
	}

	return (lam_0 + lam_r) * 1e-3; // convert mW/m/K into W/m/K
}

/**********************************************************************
   Function for calculating thermal conductivity of nitrogen (N2) depending
   on density and Temperature.

   see
   Stephan, Krauss and Laeseke: Viscosity and thermal conductivity of fluid
   Nitrogen, J. Chem. Phys. Ref. Data,Vol. 16, No. 4, 1987.

   Programming: NB
          Apr 2009
 ***********************************************************************/
double n2_heat_conductivity(double rho, double T)
{
	const double X1 = 0.95185202;
	const double X2 = 1.0205422;

	const double rho_c = 314; // [kg/m3]
	const double M = 28.013;
	const double k = 1.38062e-23;
	const double eps = 138.08483e-23;
	const double N_A = 6.02213E26;
	const double R = 8.31434;
	const double CCF = 4.173; // mW/m/K

	const double c1 = 0.3125;
	const double c2 = 2.0442e-49;
	const double sigma = 0.36502496e-09;

	double F;
	double A[5], f[9], C[4];
	double sum = 0, eta_0, c_v0, T_star, Omega = 0;
	double lamda_tr, lamda_in, lamda_r, lamda_0, lamda;

	int i;

	T_star = T * k / eps;
	rho = rho / rho_c;

	A[0] = 0.46649;
	A[1] = -0.57015;
	A[2] = 0.19164;
	A[3] = -0.03708;
	A[4] = 0.00241;

	f[0] = -0.837079888737e3;
	f[1] = 0.37914711487e2;
	f[2] = -0.601737844275;
	f[3] = 0.350418363823e1;
	f[4] = -0.874955653028e-5;
	f[5] = 0.148968607239e-7;
	f[6] = -0.256370354277e-11;
	f[7] = 0.100773735767e1;
	f[8] = 0.335340610e4;

	C[0] = 3.3373542;
	C[1] = 0.37098251;
	C[2] = 0.89913456;
	C[3] = 0.16972505;

	// dilute heat conductivity
	for (i = 0; i < 7; i++)
		sum = sum + f[i] * pow(T, (i - 3));
	const double temp(exp((f[8] / T)) - 1);
	c_v0 = R * (sum + ((f[7] * (f[8] / T) * (f[8] / T) * (exp((f[8] / T)))) / (temp * temp) - 1));
	sum = 0;

	double cvint;
	cvint = c_v0 * 1000 / N_A;

	// dilute gas viscosity
	for (i = 0; i < 5; i++)
		Omega = Omega + A[i] * MathLib::fastpow(log(T_star), i);
	Omega = exp(Omega);

	// eta in [Pa*s]
	eta_0 = 1e6 * (c1 * sqrt(c2 * T) / (sigma * sigma * Omega));

	F = eta_0 * k * N_A / (M * 1000);

	lamda_tr = 2.5 * (1.5 - X1);
	lamda_in = X2 * (cvint / k + X1);

	lamda_0 = F * (lamda_tr + lamda_in);
	sum = 0;
	for (i = 0; i < 4; i++)
		sum = sum + C[i] * MathLib::fastpow(rho, (i + 1));

	lamda_r = sum * CCF;

	lamda = (lamda_0 + lamda_r) / 1000; // lamda in [W/m/K]

	return lamda;
}

/*************************************************************
   saturation vapour pressure at a temperature.

   Wagner, W. and Kretzschmar, H.-J., 1998: International Steam
   Tables. Properties ofWater and Steam Based on the Industrial
   Formulation IAPWS-IF97. Second edition. Springer-Verlag.
   P. 25

   Input: T in K
   Output: p in bar

   Programming: Dedong Li

*************************************************************/
double vapour_pressure_IF97(double T)
{
	double A, B, C, theta, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10;
	n1 = 0.11670521452767e+04;
	n2 = -0.72421316703206e+06;
	n3 = -0.17073846940092e+02;
	n4 = 0.12020824702470e+05;
	n5 = -0.32325550322333e+07;
	n6 = 0.14915108613530e+02;
	n7 = -0.48232657361591e+04;
	n8 = 0.40511340542057e+06;
	n9 = -0.23855557567849e+00;
	n10 = 0.65017534844798e+03;

	theta = T + n9 / (T - n10);

	double theta_sqr(theta * theta);
	A = theta_sqr + n1 * theta + n2;
	B = n3 * theta_sqr + n4 * theta + n5;
	C = n6 * theta_sqr + n7 * theta + n8;
	double h;
	h = MathLib::fastpow((2.0 * C / (-B + sqrt(B * B - 4.0 * A * C))), 4) * 10.0;
	return h;
}

/*************************************************************
   Methane saturation vapour pressure at a temperature, Setzmann,1991.

   Input: T in K
   Output: p in Pa

   Programming: NB
*************************************************************/
double vapour_pressure_ch4(double T)
{
	const double Tc = 190.564;
	const double pc = 4592200;

	double theta;
	double n1, n2, n3, n4, h;

	n1 = -6.036219;
	n2 = 1.409359;
	n3 = -0.4945199;
	n4 = -1.443048;

	theta = (1 - T / Tc);
	h = Tc / T * (n1 * theta + n2 * pow(theta, 1.5) + n3 * theta * theta + n4 * pow(theta, 4.5));

	return exp(h) * pc;
}

/*************************************************************
   Water saturation vapour pressure at a temperature, Wagner,2002.

   Input: T in K
   Output: p in Pa

   Programming: NB
*************************************************************/
double vapour_pressure_h2o(double T)
{
	double theta;
	const double Tc = 647.096;
	const double pc = 22064000;

	double a1, a2, a3, a4, a5, a6, h;
	a1 = -7.85951783;
	a2 = 1.84408259;
	a3 = -11.7866497;
	a4 = 22.6807411;
	a5 = -15.9618719;
	a6 = 1.80122502;

	theta = (1 - T / Tc);

	h = Tc / T * (a1 * theta + a2 * pow(theta, 1.5) + a3 * MathLib::fastpow(theta, 3) + a4 * pow(theta, 3.5)
	              + a5 * MathLib::fastpow(theta, 4)
	              + a6 * MathLib::fastpow(theta, 7));

	return exp(h) * pc;
}

/*************************************************************
   returns min and max interpolation boundaries for zbrent

   Programming: Dedong Li
             May 2009
   Last modification: NB Jun 2009
*************************************************************/
double dpressure(double TT, double PP, int fluid, double ds)
{
	// return PP-pressure(ds, TT, fluid)*1.0E-5;
	return PP - pressure(ds, TT, fluid);
}

/*************************************************************
   Programming: Dedong Li
             May 2009
   Last modification:
*************************************************************/
inline float SIGN(const double& a, const float& b)
{
	return (float)(b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a));
}

/*************************************************************
   returns min and max interpolation boundaries for zbrent

   Programming: Dedong Li
             May 2009
   Last modification: NB Jun 2009
*************************************************************/
void auswahl(double T, double P, int fluid, double* x1, double* x2, double eps)
{
	eps = eps; // OK411
	switch (fluid)
	{
		case 0: // CO2
			*x1 = 0.0005;
			*x2 = 2000.0;

			if (T < 304.128) // gas or liquid
			{
				if (P < (vapour_pressure_co2(T))) // gas
				{
					*x1 = 5e-4; // min_density;
					*x2 = vapour_saturation_density_co2(T);
				}
				else // liquid
				{
					*x1 = liquid_saturation_density_co2(T);
					*x2 = 2000; // max_density;
				}
			}

			break;
		case 1: // Water
			if (T < 647.096 && P > (vapour_pressure_h2o(T)))
			{
				*x2 = 1500.0;
				if (T <= 440.0)
					*x1 = 9.0e+2;
				if (T > 440.0 && T < 573.15)
					*x1 = 7.1e+2;
				if (T >= 573.15)
					*x1 = 3.0e+2;
			}
			else
			{
				*x1 = 1.0e-8;
				*x2 = 2000;
				if (T < 423.0)
					*x2 = 2.6;
				if (T < 400.0)
					*x2 = 1.4;
				if (T < 335.0)
					*x2 = 0.15;
			}

			break;

		case 2: // Methane
			*x1 = 0.0005; // min_density
			*x2 = 600; // max_density

			if (T < 190.564) // gas or liquid
			{
				if (P < (vapour_pressure_ch4(T))) // gas
				{
					*x1 = 5e-4; // min_density;
					*x2 = vapour_saturation_density_ch4(T);
				}
				else // liquid
				{
					*x1 = liquid_saturation_density_ch4(T);
					*x2 = 600; // max_density;
				}
			}
			break;
		case 3: // Nitrogen

		{
			*x1 = 0.005; // min_density
			*x2 = 2000; // max_density

			if (T < 126.192) // gas or liquid
			{
				if (P < (vapour_pressure_n2(T))) // gas
				{
					*x1 = 5e-4; // min_density;
					*x2 = vapour_saturation_density_n2(T);
				}
				else // liquid
				{
					*x1 = liquid_saturation_density_n2(T);
					*x2 = 2000; // max_density;
				}
			}
		}

		break;
	}
}

//****************************************************************************80

double zero(double T, double P, int fluid, double t)

//****************************************************************************80
//
//  Purpose:
//
//    (Original):
//    ZERO seeks the root of a function F(X) in an interval [A,B].
//
//  Discussion:
//
//    The interval [A,B] must be a change of sign interval for F.
//    That is, F(A) and F(B) must be of opposite signs.  Then
//    assuming that F is continuous implies the existence of at least
//    one value C between A and B for which F(C) = 0.
//
//    The location of the zero is determined to within an accuracy
//    of 6 * MACHEPS * r8_abs ( C ) + 2 * T.
//
//    Thanks to Thomas Secretin for pointing out a transcription error in the
//    setting of the value of P, 11 February 2013.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modification:
//    This implementation of zero seeks the density for a specific pressure
//    and temperature condition according to a specified fluid
//
//  Modified:
//
//    11 February 2013
//    31 July 2015
//
//  Author:
//
//    Original FORTRAN77 version by Richard Brent.
//    C++ version by John Burkardt.
//    Modified for OGS purpose by N Boettcher
//
//  Reference:
//
//    Richard Brent,
//    Algorithms for Minimization Without Derivatives,
//    Dover, 2002,
//    ISBN: 0-486-41998-3,
//    LC: QA402.5.B74.
//
//  Parameters:
//
//    Input: double T, P, --> temperature and pressure
//           int fluid --> an index specifying a fluid
//
//    Output: estimated density value according to specified fluid, pressure,
//    and temperature.
//
{
	double c;
	double d;
	double e;
	double fa;
	double fb;
	double fc;
	double m;
	double macheps;
	double p;
	double q;
	double r;
	double s;
	double sa;
	double sb;
	double tol;
	//
	//  Make local copies of A and B.
	//
	double a, b;
	auswahl(T, P, fluid, &a, &b, t);

	sa = a;
	sb = b;
	fa = dpressure(T, P, fluid, sa);
	fb = dpressure(T, P, fluid, sb);

	c = sa;
	fc = fa;
	e = sb - sa;
	d = e;

	macheps = std::numeric_limits<double>::epsilon();

	for (;;)
	{
		if (std::abs(fc) < std::abs(fb))
		{
			sa = sb;
			sb = c;
			c = sa;
			fa = fb;
			fb = fc;
			fc = fa;
		}

		tol = 2.0 * macheps * std::abs(sb) + t;
		m = 0.5 * (c - sb);

		if (std::abs(m) <= tol || fb == 0.0)
		{
			break;
		}

		if (std::abs(e) < tol || std::abs(fa) <= std::abs(fb))
		{
			e = m;
			d = e;
		}
		else
		{
			s = fb / fa;

			if (sa == c)
			{
				p = 2.0 * m * s;
				q = 1.0 - s;
			}
			else
			{
				q = fa / fc;
				r = fb / fc;
				p = s * (2.0 * m * q * (q - r) - (sb - sa) * (r - 1.0));
				q = (q - 1.0) * (r - 1.0) * (s - 1.0);
			}

			if (0.0 < p)
			{
				q = -q;
			}
			else
			{
				p = -p;
			}

			s = e;
			e = d;

			if (2.0 * p < 3.0 * m * q - std::abs(tol * q) && p < std::abs(0.5 * s * q))
			{
				d = p / q;
			}
			else
			{
				e = m;
				d = e;
			}
		}
		sa = sb;
		fa = fb;

		if (tol < std::abs(d))
		{
			sb = sb + d;
		}
		else if (0.0 < m)
		{
			sb = sb + tol;
		}
		else
		{
			sb = sb - tol;
		}

		fb = dpressure(T, P, fluid, sb);

		if ((0.0 < fb && 0.0 < fc) || (fb <= 0.0 && fc <= 0.0))
		{
			c = sa;
			fc = fa;
			e = sb - sa;
			d = e;
		}
	}
	return sb;
}

/**********************************************************************
   Function for the mixing of adjusting parameters for cubic equations of state.
   For ternary mixtures.

   Fluid List:  0 - carbon dioxide
             1 - water
             2 - methane
             3 - nitrogen
             ...

   Arguments: x[i] ... mole fraction of component i
   a[i],b[i] ... parameters for pure substance
   MM[i]... molar mass of pure substance
   ra,rb,MM   ... mixed parameters

   Programming: NB
   May 2009
 ***********************************************************************/
double mixing_ternary(double* x, double* a, double* b, double* MM, double* ra, double* rb, double* rMM)
{
	int unsigned i, j, k;
	double a_mix = 0, b_mix = 0, MM_mix = 0;

	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			for (k = 0; k < 3; k++)
				a_mix += x[i] * x[j] * x[k] * pow(a[i] * a[j] * a[k], (1. / 3.));

	for (i = 0; i < 3; i++)
		b_mix += x[i] * b[i];

	for (i = 0; i < 3; i++)
		MM_mix += x[i] * MM[i];

	*ra = a_mix;
	*rb = b_mix;
	*rMM = MM_mix;

	return 0;
}

/*************************************************************
   Methane liquid saturation density at a temperature, Setzmann,1991.

   Input: T in K
   Output: p in Pa

   Programming: NB
*************************************************************/
double liquid_saturation_density_ch4(double T)
{
	const double Tc = 190.564;
	const double rhoc = 162.66;

	double theta;
	double n1, n2, n3, h;

	n1 = 1.9906389;
	n2 = -0.78756197;
	n3 = 0.036976723;

	theta = (1 - T / Tc);
	h = (n1 * pow(theta, 0.354) + n2 * sqrt(theta) + n3 * pow(theta, (5. / 2.)));

	return exp(h) * rhoc;
}

/*************************************************************
   Methane vapour saturation density at a temperature, Setzmann,1991.

   Input: T in K
   Output: p in Pa

   Programming: NB
*************************************************************/
double vapour_saturation_density_ch4(double T)
{
	const double Tc = 190.564;
	const double rhoc = 162.66;

	double theta;
	double n1, n2, n3, n4, n5, n6, h;

	n1 = -1.880284;
	n2 = -2.8526531;
	n3 = -3.000648;
	n4 = -5.251169;
	n5 = -13.191859;
	n6 = -37.553961;

	theta = (1 - T / Tc);
	h = (n1 * pow(theta, 0.354) + n2 * pow(theta, (5. / 6.)) + n3 * pow(theta, (3. / 2.)) + n4 * pow(theta, (5. / 2.))
	     + n5 * pow(theta, (25. / 6.))
	     + n6 * pow(theta, (47. / 6.)));

	return exp(h) * rhoc;
}

/*************************************************************
   Nitrogen vapour pressure at a temperature, Setzmann,1991.

   Input: T in K
   Output: p in Pa

   Programming: NB
*************************************************************/
double vapour_pressure_n2(double T)
{
	const double Tc = 126.192;
	const double pc = 3395800.;

	double theta;
	double n1, n2, n3, n4, h;

	n1 = -6.12445284;
	n2 = 1.26327220;
	n3 = -0.765910082;
	n4 = -1.77570564;

	theta = (1 - T / Tc);
	h = (Tc / T) * (n1 * theta + n2 * pow(theta, (1.5)) + n3 * pow(theta, (2.5)) + n4 * MathLib::fastpow(theta, (5)));

	return exp(h) * pc;
}

/*************************************************************
   Nitrogen liquid saturation density at a temperature, Setzmann,1991.

   Input: T in K
   Output: density in kg/m3

   Programming: NB
*************************************************************/
double liquid_saturation_density_n2(double T)
{
	const double Tc = 126.192;
	const double rhoc = 313.3;

	double theta;
	double n1, n2, n3, n4, h; // OK411 n5,n6,

	n1 = 1.48654237;
	n2 = -0.280476066;
	n3 = 0.0894143085;
	n4 = -0.119879866;

	theta = (1 - T / Tc);
	h = (n1 * pow(theta, 0.3294) + n2 * pow(theta, (2. / 3.)) + n3 * pow(theta, (8. / 3.))
	     + n4 * pow(theta, (35. / 6.)));

	return exp(h) * rhoc;
}

/*************************************************************
   Nitrogen vapour saturation density at a temperature, Setzmann,1991.

   Input: T in K
   Output: density in kg/m3

   Programming: NB
*************************************************************/
double vapour_saturation_density_n2(double T)
{
	const double Tc = 126.192;
	const double rhoc = 313.3;

	double theta;
	double n1, n2, n3, n4, n5, h; // OK411 n6,

	n1 = -1.70127164;
	n2 = -3.70402649;
	n3 = 1.29859383;
	n4 = -0.561424977;
	n5 = -2.68505381;

	theta = (1 - T / Tc);
	h = (Tc / T)
	    * (n1 * pow(theta, 0.34) + n2 * pow(theta, (5. / 6.)) + n3 * pow(theta, (7. / 6.)) + n4 * pow(theta, (13. / 6.))
	       + n5 * pow(theta, (14. / 3.)));

	return exp(h) * rhoc;
}

/*************************************************************
   Carbon dioxide liquid saturation density at a temperature, Span,1996.

   Input: T in K
   Output: density in kg/m3

   Programming: NB
*************************************************************/
double liquid_saturation_density_co2(double T)
{
	const double Tc = 304.128;
	const double rhoc = 467.6;

	double n[4], t[4], h = 0;
	int i;

	n[0] = 1.9245108;
	t[0] = 0.34;
	n[1] = -0.62385555;
	t[1] = 0.5;
	n[2] = -0.32731127;
	t[2] = (10. / 6.);
	n[3] = 0.39245142;
	t[3] = (11. / 6.);

	for (i = 0; i < 4; i++)
		h += n[i] * pow((1 - T / Tc), t[i]);

	return exp(h) * rhoc;
}

/*************************************************************
   Carbon dioxide vapour saturation density at a temperature, Span,1996.

   Input: T in K
   Output: density in kg/m

   Programming: NB
*************************************************************/
double vapour_saturation_density_co2(double T)
{
	const double Tc = 304.128;
	const double rhoc = 467.6;

	double n[5], t[5], h = 0;
	int i;

	n[0] = -1.7074879;
	t[0] = 0.34;
	n[1] = -0.8227467;
	t[1] = 0.5;
	n[2] = -4.6008549;
	t[2] = 1;
	n[3] = -10.111178;
	t[3] = (7. / 3.);
	n[4] = -29.742252;
	t[4] = (14. / 3.);

	for (i = 0; i < 5; i++)
		h += n[i] * pow((1 - T / Tc), t[i]);

	return exp(h) * rhoc;
}

/**********************************************************************
   Function thermal_properties (fluid, critical_density, critical_temperature, specific_gas_constant)
   returns the thermal properties of a given fluid
   Programming: NB Mar09
**********************************************************************/
void CFluidProperties::therm_prop(string caption)
{
	// CFluidProperties fpc;
	char letter;
	int i, j;
	// TODO: Change first letter approach to an integer (NB) - done (fluid_number) NB, 2009-06-02
	letter = caption[0];
	Ru = 8.314472;
	switch (letter)
	{
		case 'C': // CARBON DIOXIDE
		{
			fluid_id = 0;
			rhoc = 467.6; // critical density [kg/m3]
			Tc = 304.1282; // critical temperature [K]
			pc = 7377300; // critical pressure [Pa]
			Tt = 216.592; // triple point temperature [K]
			pt = 517950; //  triple point pressure [Pa]
			Rs = 188.9241; // specific gas constant [J/kg/K]
			molar_mass = 44.0099; // [g/mol]
			omega = 0.22491; // azentric factor, see PREOS
			Vd = 26.9;
			Zc = 0.27468; // critical super-compressibility, see PREOS
			n0 = 0.11333;
			k3 = 0.28996;
			m0 = 0.38493912223725674;
			a = 383766.38336903340;
			b = 0.026666761740035457;
			k1 = 0.012304583901235679;
			k2 = -0.14826846628826726;
			// Limits sums in FHE-derivations

			limit[0] = 7;
			limit[1] = 34;
			limit[2] = 39;
			limit[3] = 42;

			// Coefficients for FHE-derivations

			for (i = 0; i < 2; i++)
				for (j = 0; j < 8; j++)
					k[i][j] = 0;

			for (i = 0; i < 14; i++)
				for (j = 0; j < 56; j++)
					K[i][j] = 0;

			// ideal gas part
			// a
			k[0][0] = 8.37304456;
			k[0][1] = -3.70454304;
			k[0][2] = 2.5;
			k[0][3] = 1.99427042;
			k[0][4] = 0.62105248;
			k[0][5] = 0.41195293;
			k[0][6] = 1.04028922;
			k[0][7] = 0.08327678;
			// theta
			k[1][0] = 0;
			k[1][1] = 0;
			k[1][2] = 0;
			k[1][3] = 3.15163;
			k[1][4] = 6.1119;
			k[1][5] = 6.77708;
			k[1][6] = 11.32384;
			k[1][7] = 27.08792;

			// real gas part

			// n
			K[0][0] = 3.8856823203161E-01;
			K[0][1] = 2.9385475942740E+00;
			K[0][2] = -5.5867188534934E+00;
			K[0][3] = -7.6753199592477E-01;
			K[0][4] = 3.1729005580416E-01;
			K[0][5] = 5.4803315897767E-01;
			K[0][6] = 1.2279411220335E-01;
			K[0][7] = 2.1658961543220E+00;
			K[0][8] = 1.5841735109724E+00;
			K[0][9] = -2.3132705405503E-01;
			K[0][10] = 5.8116916431436E-02;
			K[0][11] = -5.5369137205382E-01;
			K[0][12] = 4.8946615909422E-01;
			K[0][13] = -2.4275739843501E-02;
			K[0][14] = 6.2494790501678E-02;
			K[0][15] = -1.2175860225246E-01;
			K[0][16] = -3.7055685270086E-01;
			K[0][17] = -1.6775879700426E-02;
			K[0][18] = -1.1960736637987E-01;
			K[0][19] = -4.5619362508778E-02;
			K[0][20] = 3.5612789270346E-02;
			K[0][21] = -7.4427727132052E-03;
			K[0][22] = -1.7395704902432E-03;
			K[0][23] = -2.1810121289527E-02;
			K[0][24] = 2.4332166559236E-02;
			K[0][25] = -3.7440133423463E-02;
			K[0][26] = 1.4338715756878E-01;
			K[0][27] = -1.3491969083286E-01;
			K[0][28] = -2.3151225053480E-02;
			K[0][29] = 1.2363125492901E-02;
			K[0][30] = 2.1058321972940E-03;
			K[0][31] = -3.3958519026368E-04;
			K[0][32] = 5.5993651771592E-03;
			K[0][33] = -3.0335118055646E-04;
			K[0][34] = -2.1365488688320E+02;
			K[0][35] = 2.6641569149272E+04;
			K[0][36] = -2.4027212204557E+04;
			K[0][37] = -2.8341603423999E+02;
			K[0][38] = 2.1247284400179E+02;
			K[0][39] = -6.6642276540751E-01;
			K[0][40] = 7.2608632349897E-01;
			K[0][41] = 5.5068668612842E-02;

			// d
			K[1][0] = 1;
			K[1][1] = 1;
			K[1][2] = 1;
			K[1][3] = 1;
			K[1][4] = 2;
			K[1][5] = 2;
			K[1][6] = 3;
			K[1][7] = 1;
			K[1][8] = 2;
			K[1][9] = 4;
			K[1][10] = 5;
			K[1][11] = 5;
			K[1][12] = 5;
			K[1][13] = 6;
			K[1][14] = 6;
			K[1][15] = 6;
			K[1][16] = 1;
			K[1][17] = 1;
			K[1][18] = 4;
			K[1][19] = 4;
			K[1][20] = 4;
			K[1][21] = 7;
			K[1][22] = 8;
			K[1][23] = 2;
			K[1][24] = 3;
			K[1][25] = 3;
			K[1][26] = 5;
			K[1][27] = 5;
			K[1][28] = 6;
			K[1][29] = 7;
			K[1][30] = 8;
			K[1][31] = 10;
			K[1][32] = 4;
			K[1][33] = 8;
			K[1][34] = 2;
			K[1][35] = 2;
			K[1][36] = 2;
			K[1][37] = 3;
			K[1][38] = 3;

			// t
			K[2][0] = 0;
			K[2][1] = 0.75;
			K[2][2] = 1;
			K[2][3] = 2;
			K[2][4] = 0.75;
			K[2][5] = 2;
			K[2][6] = 0.75;
			K[2][7] = 1.5;
			K[2][8] = 1.5;
			K[2][9] = 2.5;
			K[2][10] = 0;
			K[2][11] = 1.5;
			K[2][12] = 2;
			K[2][13] = 0;
			K[2][14] = 1;
			K[2][15] = 2;
			K[2][16] = 3;
			K[2][17] = 6;
			K[2][18] = 3;
			K[2][19] = 6;
			K[2][20] = 8;
			K[2][21] = 6;
			K[2][22] = 0;
			K[2][23] = 7;
			K[2][24] = 12;
			K[2][25] = 16;
			K[2][26] = 22;
			K[2][27] = 24;
			K[2][28] = 16;
			K[2][29] = 24;
			K[2][30] = 8;
			K[2][31] = 2;
			K[2][32] = 28;
			K[2][33] = 14;
			K[2][34] = 1;
			K[2][35] = 0;
			K[2][36] = 1;
			K[2][37] = 3;
			K[2][38] = 3;

			// c
			K[3][7] = 1;
			K[3][8] = 1;
			K[3][9] = 1;
			K[3][10] = 1;
			K[3][11] = 1;
			K[3][12] = 1;
			K[3][13] = 1;
			K[3][14] = 1;
			K[3][15] = 1;
			K[3][16] = 2;
			K[3][17] = 2;
			K[3][18] = 2;
			K[3][19] = 2;
			K[3][20] = 2;
			K[3][21] = 2;
			K[3][22] = 2;
			K[3][23] = 3;
			K[3][24] = 3;
			K[3][25] = 3;
			K[3][26] = 4;
			K[3][27] = 4;
			K[3][28] = 4;
			K[3][29] = 4;
			K[3][30] = 4;
			K[3][31] = 4;
			K[3][32] = 5;
			K[3][33] = 6;

			// a
			K[4][39] = 3.5;
			K[4][40] = 3.5;
			K[4][41] = 3;
			// b
			K[5][39] = 0.875;
			K[5][40] = 0.925;
			K[5][41] = 0.875;
			// A
			K[6][39] = 0.7;
			K[6][40] = 0.7;
			K[6][41] = 0.7;
			// B
			K[7][39] = 0.3;
			K[7][40] = 0.3;
			K[7][41] = 1;
			// C
			K[8][39] = 10;
			K[8][40] = 10;
			K[8][41] = 12.5;
			// D
			K[9][39] = 275;
			K[9][40] = 275;
			K[9][41] = 275;

			// alpha
			K[10][34] = 25;
			K[10][35] = 25;
			K[10][36] = 25;
			K[10][37] = 15;
			K[10][38] = 20;
			// beta
			K[11][34] = 325;
			K[11][35] = 300;
			K[11][36] = 300;
			K[11][37] = 275;
			K[11][38] = 275;
			K[11][39] = 0.3;
			K[11][40] = 0.3;
			K[11][41] = 0.3;
			// gamma
			K[12][34] = 1.16;
			K[12][35] = 1.19;
			K[12][36] = 1.19;
			K[12][37] = 1.25;
			K[12][38] = 1.22;
			// epsilon
			K[13][34] = 1;
			K[13][35] = 1;
			K[13][36] = 1;
			K[13][37] = 1;
			K[13][38] = 1;

			break;
		}
		case 'W': // WATER
		{
			fluid_id = 1;
			rhoc = 322; //[kg/m3]
			Tc = 647.096; //[K]
			pc = 22064000; // [Pa]
			Tt = PhysicalConstant::CelsiusZeroInKelvin; //  [K]
			pt = 611.657; //  [Pa]
			Rs = 461.51805; //  [J/kg/K]
			molar_mass = 18.01528; //  [g/mol]
			omega = 0.344; // azentric factor, see PREOS
			Vd = 25.14;
			Zc = 0.22944; // critical super-compressibility, see PREOS
			n0 = 0.1156;
			k3 = 0.0471;
			m0 = 0.47568277359898614;
			a = 943391.02482869523;
			b = 0.018971230469153735;
			k1 = 0.017189421358489602;
			k2 = -0.029385598856191408;

			// Limits for Sums in FHE-derivations

			limit[0] = 7;
			limit[1] = 51;
			limit[2] = 54;
			limit[3] = 56;

			// Coefficients for FHE-derivations

			for (i = 0; i < 2; i++)
				for (j = 0; j < 8; j++)
					k[i][j] = 0;

			for (i = 0; i < 14; i++)
				for (j = 0; j < 56; j++)
					K[i][j] = 0;

			// ideal gas part

			k[0][0] = -8.32044648201;
			k[0][1] = 6.6832105268;
			k[0][2] = 3.00632;
			k[0][3] = 0.012436;
			k[0][4] = 0.97315;
			k[0][5] = 1.27950;
			k[0][6] = 0.96956;
			k[0][7] = 0.24873;

			k[1][0] = 0;
			k[1][1] = 0;
			k[1][2] = 0;
			k[1][3] = 1.28728967;
			k[1][4] = 3.53734222;
			k[1][5] = 7.74073708;
			k[1][6] = 9.24437796;
			k[1][7] = 27.5075105;

			// real gas part

			K[0][0] = 1.2533547935523E-02;
			K[0][1] = 7.8957634722828E+00;
			K[0][2] = -8.7803203303561E+00;
			K[0][3] = 3.1802509345418E-01;
			K[0][4] = -2.6145533859358E-01;
			K[0][5] = -7.8199751687981E-03;
			K[0][6] = 8.8089493102134E-03;
			K[0][7] = -6.6856572307965E-01;
			K[0][8] = 2.0433810950965E-01;
			K[0][9] = -6.6212605039687E-05;
			K[0][10] = -1.9232721156002E-01;
			K[0][11] = -2.5709043003438E-01;
			K[0][12] = 1.6074868486251E-01;
			K[0][13] = -4.0092828925807E-02;
			K[0][14] = 3.9343422603254E-07;

			K[0][15] = -7.5941377088144E-06;
			K[0][16] = 5.6250979351888E-04;
			K[0][17] = -1.5608652257135E-05;
			K[0][18] = 1.1537996422951E-09;
			K[0][19] = 3.6582165144204E-07;
			K[0][20] = -1.3251180074668E-12;
			K[0][21] = -6.2639586912454E-10;
			K[0][22] = -1.0793600908932E-01;
			K[0][23] = 1.7611491008752E-02;
			K[0][24] = 2.2132295167546E-01;
			K[0][25] = -4.0247669763528E-01;
			K[0][26] = 5.8083399985759E-01;
			K[0][27] = 4.9969146990806E-03;
			K[0][28] = -3.1358700712549E-02;
			K[0][29] = -7.4315929710341E-01;

			K[0][30] = 4.7807329915480E-01;
			K[0][31] = 2.0527940895948E-02;
			K[0][32] = -1.3636435110343E-01;
			K[0][33] = 1.4180634400617E-02;
			K[0][34] = 8.3326504880713E-03;
			K[0][35] = -2.9052336009585E-02;
			K[0][36] = 3.8615085574206E-02;
			K[0][37] = -2.0393486513704E-02;
			K[0][38] = -1.6554050063734E-03;
			K[0][39] = 1.9955571979541E-03;
			K[0][40] = 1.5870308324157E-04;
			K[0][41] = -1.6388568342530E-05;
			K[0][42] = 4.3613615723811E-02;
			K[0][43] = 3.4994005463765E-02;
			K[0][44] = -7.6788197844621E-02;

			K[0][45] = 2.2446277332006E-02;
			K[0][46] = -6.2689710414685E-05;
			K[0][47] = -5.5711118565645E-10;
			K[0][48] = -1.9905718354408E-01;
			K[0][49] = 3.1777497330738E-01;
			K[0][50] = -1.1841182425981E-01;
			K[0][51] = -3.1306260323435E+01;
			K[0][52] = 3.1546140237781E+01;
			K[0][53] = -2.5213154341695E+03;
			K[0][54] = -1.4874640856724E-01;
			K[0][55] = 3.1806110878444E-01;

			K[1][0] = 1;
			K[1][1] = 1;
			K[1][2] = 1;
			K[1][3] = 2;
			K[1][4] = 2;
			K[1][5] = 3;
			K[1][6] = 4;
			K[1][7] = 1;
			K[1][8] = 1;
			K[1][9] = 1;
			K[1][10] = 2;
			K[1][11] = 2;
			K[1][12] = 3;
			K[1][13] = 4;
			K[1][14] = 4;
			K[1][15] = 5;
			K[1][16] = 7;
			K[1][17] = 9;
			K[1][18] = 10;
			K[1][19] = 11;
			K[1][20] = 13;
			K[1][21] = 15;
			K[1][22] = 1;
			K[1][23] = 2;
			K[1][24] = 2;
			K[1][25] = 2;
			K[1][26] = 3;
			K[1][27] = 4;
			K[1][28] = 4;
			K[1][29] = 4;
			K[1][30] = 5;
			K[1][31] = 6;
			K[1][32] = 6;
			K[1][33] = 7;
			K[1][34] = 9;
			K[1][35] = 9;
			K[1][36] = 9;
			K[1][37] = 9;
			K[1][38] = 9;
			K[1][39] = 10;
			K[1][40] = 10;
			K[1][41] = 12;
			K[1][42] = 3;
			K[1][43] = 4;
			K[1][44] = 4;
			K[1][45] = 5;
			K[1][46] = 14;
			K[1][47] = 3;
			K[1][48] = 6;
			K[1][49] = 6;
			K[1][50] = 6;
			K[1][51] = 3;
			K[1][52] = 3;
			K[1][53] = 3;

			K[2][0] = -0.5;
			K[2][1] = 0.875;
			K[2][2] = 1;
			K[2][3] = 0.5;
			K[2][4] = 0.75;
			K[2][5] = 0.375;
			K[2][6] = 1;
			K[2][7] = 4;
			K[2][8] = 6;
			K[2][9] = 12;
			K[2][10] = 1;
			K[2][11] = 5;
			K[2][12] = 4;
			K[2][13] = 2;
			K[2][14] = 13;
			K[2][15] = 9;
			K[2][16] = 3;
			K[2][17] = 4;
			K[2][18] = 11;
			K[2][19] = 4;
			K[2][20] = 13;
			K[2][21] = 1;
			K[2][22] = 7;
			K[2][23] = 1;
			K[2][24] = 9;
			K[2][25] = 10;
			K[2][26] = 10;
			K[2][27] = 3;
			K[2][28] = 7;
			K[2][29] = 10;
			K[2][30] = 10;
			K[2][31] = 6;
			K[2][32] = 10;
			K[2][33] = 10;
			K[2][34] = 1;
			K[2][35] = 2;
			K[2][36] = 3;
			K[2][37] = 4;
			K[2][38] = 8;
			K[2][39] = 6;
			K[2][40] = 9;
			K[2][41] = 8;
			K[2][42] = 16;
			K[2][43] = 22;
			K[2][44] = 23;
			K[2][45] = 23;
			K[2][46] = 10;
			K[2][47] = 50;
			K[2][48] = 44;
			K[2][49] = 46;
			K[2][50] = 50;
			K[2][51] = 0;
			K[2][52] = 1;
			K[2][53] = 4;

			K[3][7] = 1;
			K[3][8] = 1;
			K[3][9] = 1;
			K[3][10] = 1;
			K[3][11] = 1;
			K[3][12] = 1;
			K[3][13] = 1;
			K[3][14] = 1;
			K[3][15] = 1;
			K[3][16] = 1;
			K[3][17] = 1;
			K[3][18] = 1;
			K[3][19] = 1;
			K[3][20] = 1;
			K[3][21] = 1;
			K[3][22] = 2;
			K[3][23] = 2;
			K[3][24] = 2;
			K[3][25] = 2;
			K[3][26] = 2;
			K[3][27] = 2;
			K[3][28] = 2;
			K[3][29] = 2;
			K[3][30] = 2;
			K[3][31] = 2;
			K[3][32] = 2;
			K[3][33] = 2;
			K[3][34] = 2;
			K[3][35] = 2;
			K[3][36] = 2;
			K[3][37] = 2;
			K[3][38] = 2;
			K[3][39] = 2;
			K[3][40] = 2;
			K[3][41] = 2;
			K[3][42] = 3;
			K[3][43] = 3;
			K[3][44] = 3;
			K[3][45] = 3;
			K[3][46] = 4;
			K[3][47] = 6;
			K[3][48] = 6;
			K[3][49] = 6;
			K[3][50] = 6;

			K[4][54] = 3.5;
			K[4][55] = 3.5;
			K[5][54] = 0.85;
			K[5][55] = 0.95;
			K[6][54] = 0.32;
			K[6][55] = 0.32;
			K[7][54] = 0.2;
			K[7][55] = 0.2;
			K[8][54] = 28;
			K[8][55] = 32;
			K[9][54] = 700;
			K[9][55] = 800;

			K[10][51] = 20;
			K[10][52] = 20;
			K[10][53] = 20;
			K[11][51] = 150;
			K[11][52] = 150;
			K[11][53] = 250;
			K[11][54] = 0.3;
			K[11][55] = 0.3;
			K[12][51] = 1.21;
			K[12][52] = 1.21;
			K[12][53] = 1.25;
			K[13][51] = 1;
			K[13][52] = 1;
			K[13][53] = 1;

			break;
		}
		case 'M': // METHANE
		{
			fluid_id = 2;
			rhoc = 162.66; //[kg/m3]
			Tc = 190.551; //[K]
			pc = 4599200; // [Pa]
			Tt = 90.685; //  [K]
			pt = 11696; //  [Pa]
			Rs = 518.3; //  [J/kg/K]
			molar_mass = 16.04; //  [g/mol]
			omega = 0.011; // azentric factor, see PREOS
			Vd = 25.14;
			Zc = 0.286060; // critical super-compressibility, see PREOS
			n0 = 0.08248;
			k3 = 0.20978;
			m0 = 0.21389815179757277;
			a = 206793.24123880462;
			b = 0.026800319421855536;
			k1 = 0.0019409288415128503;
			k2 = -0.11003457942904071;
			// Limits sums in FHE-derivations

			limit[0] = 13;
			limit[1] = 36;
			limit[2] = 40;
			limit[3] = 40;

			// Coefficients for FHE-derivations

			for (i = 0; i < 2; i++)
				for (j = 0; j < 8; j++)
					k[i][j] = 0;

			for (i = 0; i < 14; i++)
				for (j = 0; j < 56; j++)
					K[i][j] = 0;

			// ideal gas part

			k[0][0] = 9.91243972;
			k[0][1] = -6.33270087;
			k[0][2] = 3.0016;
			k[0][3] = 0.008449;
			k[0][4] = 4.6942;
			k[0][5] = 3.4865;
			k[0][6] = 1.6572;
			k[0][7] = 1.4115;

			k[1][0] = 0;
			k[1][1] = 0;
			k[1][2] = 0;
			k[1][3] = 3.4004324;
			k[1][4] = 10.26951575;
			k[1][5] = 20.43932747;
			k[1][6] = 29.93744884;
			k[1][7] = 79.13351945;

			// real gas part

			K[0][0] = 4.368E-02;
			K[0][1] = 6.709E-01;
			K[0][2] = -1.766E+00;
			K[0][3] = 8.582E-01;
			K[0][4] = -1.207E+00;
			K[0][5] = 5.120E-01;
			K[0][6] = -4.000E-04;
			K[0][7] = -1.248E-02;
			K[0][8] = 3.100E-02;
			K[0][9] = 1.755E-03;
			K[0][10] = -3.172E-06;
			K[0][11] = -2.240E-06;
			K[0][12] = 2.947E-07;
			K[0][13] = 1.830E-01;
			K[0][14] = 1.512E-01;
			K[0][15] = -4.289E-01;
			K[0][16] = 6.894E-02;
			K[0][17] = -1.408E-02;
			K[0][18] = -3.063E-02;
			K[0][19] = -2.970E-02;
			K[0][20] = -1.932E-02;
			K[0][21] = -1.106E-01;
			K[0][22] = 9.953E-02;
			K[0][23] = 8.548E-03;
			K[0][24] = -6.151E-02;
			K[0][25] = -4.292E-02;
			K[0][26] = -1.813E-02;
			K[0][27] = 3.446E-02;
			K[0][28] = -2.386E-03;
			K[0][29] = -1.159E-02;
			K[0][30] = 6.642E-02;
			K[0][31] = -2.372E-02;
			K[0][32] = -3.962E-02;
			K[0][33] = -1.387E-02;
			K[0][34] = 3.389E-02;
			K[0][35] = -2.927E-03;
			K[0][36] = 9.325E-05;
			K[0][37] = -6.287E+00;
			K[0][38] = 1.271E+01;
			K[0][39] = -6.424E+00;

			K[3][13] = 1;
			K[3][14] = 1;
			K[3][15] = 1;
			K[3][16] = 1;
			K[3][17] = 1;
			K[3][18] = 1;
			K[3][19] = 1;
			K[3][20] = 2;
			K[3][21] = 2;
			K[3][22] = 2;
			K[3][23] = 2;
			K[3][24] = 2;
			K[3][25] = 3;
			K[3][26] = 3;
			K[3][27] = 3;
			K[3][28] = 3;
			K[3][29] = 4;
			K[3][30] = 4;
			K[3][31] = 4;
			K[3][32] = 4;
			K[3][33] = 4;
			K[3][34] = 4;
			K[3][35] = 4;

			K[1][0] = 1;
			K[1][1] = 1;
			K[1][2] = 1;
			K[1][3] = 2;
			K[1][4] = 2;
			K[1][5] = 2;
			K[1][6] = 2;
			K[1][7] = 3;
			K[1][8] = 4;
			K[1][9] = 4;
			K[1][10] = 8;
			K[1][11] = 9;
			K[1][12] = 10;
			K[1][13] = 1;
			K[1][14] = 1;
			K[1][15] = 1;
			K[1][16] = 2;
			K[1][17] = 4;
			K[1][18] = 5;
			K[1][19] = 6;
			K[1][20] = 1;
			K[1][21] = 2;
			K[1][22] = 3;
			K[1][23] = 4;
			K[1][24] = 4;
			K[1][25] = 3;
			K[1][26] = 5;
			K[1][27] = 5;
			K[1][28] = 8;
			K[1][29] = 2;
			K[1][30] = 3;
			K[1][31] = 4;
			K[1][32] = 4;
			K[1][33] = 4;
			K[1][34] = 5;
			K[1][35] = 6;
			K[1][36] = 2;
			K[1][37] = 0;
			K[1][38] = 0;
			K[1][39] = 0;

			K[2][0] = -0.5;
			K[2][1] = 0.5;
			K[2][2] = 1;
			K[2][3] = 0.5;
			K[2][4] = 1;
			K[2][5] = 1.5;
			K[2][6] = 4.5;
			K[2][7] = 0;
			K[2][8] = 1;
			K[2][9] = 3;
			K[2][10] = 1;
			K[2][11] = 3;
			K[2][12] = 3;
			K[2][13] = 0;
			K[2][14] = 1;
			K[2][15] = 2;
			K[2][16] = 0;
			K[2][17] = 0;
			K[2][18] = 2;
			K[2][19] = 2;
			K[2][20] = 5;
			K[2][21] = 5;
			K[2][22] = 5;
			K[2][23] = 2;
			K[2][24] = 4;
			K[2][25] = 12;
			K[2][26] = 8;
			K[2][27] = 10;
			K[2][28] = 10;
			K[2][29] = 10;
			K[2][30] = 14;
			K[2][31] = 12;
			K[2][32] = 18;
			K[2][33] = 22;
			K[2][34] = 18;
			K[2][35] = 14;
			K[2][36] = 2;
			K[2][37] = 0;
			K[2][38] = 1;
			K[2][39] = 2;
			K[10][36] = 20;
			K[10][37] = 40;
			K[10][38] = 40;
			K[10][39] = 40;
			K[11][36] = 200;
			K[11][37] = 250;
			K[11][38] = 250;
			K[11][39] = 250;
			K[12][36] = 1.07;
			K[12][37] = 1.11;
			K[12][38] = 1.11;
			K[12][39] = 1.11;
			K[13][36] = 1;
			K[13][37] = 1;
			K[13][38] = 1;
			K[13][39] = 1;
			break;
		}
		case 'N': // Nitrogen
		{
			fluid_id = 3;
			rhoc = 314.0; //[kg/m3]
			Tc = 126.20; //[K]
			pc = 3383000; // [Pa]
			Tt = 63.148; //  [K]
			pt = 12500; //  [Pa]
			Rs = 296.8; //  [J/kg/K]
			molar_mass = 28.013; //  [g/mol]
			omega = 0.039; // azentric factor, see PREOS
			Vd = 18.5;
			Zc = 0.287634; // critical super-compressibility, see PREOS
			n0 = 0.09967;
			k3 = 0.24086;
			m0 = 0.23704245415214481;
			a = 0.91551836047871149;
			b = 0.024130615006680459;
			k1 = 0.0025206904456128499;
			k2 = -0.12498696464510012;
			// Limits sums in FHE-derivations
			limit[0] = 6;
			limit[1] = 32;
			limit[2] = 36;
			limit[3] = 36;

			// Coefficients for FHE-derivations
			// ideal gas part
			k[0][0] = 9.912644;
			k[0][1] = -6.333133;
			k[0][2] = 3.0016;
			k[0][3] = 0.008449;
			k[0][4] = 4.6942;
			k[0][5] = 3.4865;
			k[0][6] = 1.6572;
			k[0][7] = 1.4115;

			k[1][0] = 0;
			k[1][1] = 0;
			k[1][2] = 0;
			k[1][3] = 3.400664;
			k[1][4] = 10.27022;
			k[1][5] = 20.44072;
			k[1][6] = 29.93949;
			k[1][7] = 79.13892;

			// real gas part
			K[0][0] = 0.924803575275;
			K[0][1] = -0.492448489428;
			K[0][2] = 0.661883336938;
			K[0][3] = -0.192902649201e1;
			K[0][4] = -0.622469309629e-1;
			K[0][5] = 0.349943957581;
			K[0][6] = 0.564857472498;
			K[0][7] = -0.161720005987e1;
			K[0][8] = -0.481395031883;
			K[0][9] = 0.421150636384;
			K[0][10] = -0.161962230825e-1;
			K[0][11] = 0.172100994165;
			K[0][12] = 0.735448924933e-2;
			K[0][13] = 0.168077305479e-1;
			K[0][14] = -0.107626664179e-2;
			K[0][15] = -0.137318088513e-1;
			K[0][16] = 0.635466899859e-3;
			K[0][17] = 0.304432279419e-2;
			K[0][18] = -0.435762336045e-1;
			K[0][19] = -0.723174889316e-1;
			K[0][20] = 0.389644315272e-1;
			K[0][21] = -0.212201363910e-1;
			K[0][22] = 0.408822981509e-2;
			K[0][23] = -0.551990017984e-4;
			K[0][24] = -0.462016716479e-1;
			K[0][25] = -0.300311716011e-2;
			K[0][26] = 0.368825891208e-1;
			K[0][27] = -0.255856846220e-2;
			K[0][28] = 0.896915264558e-2;
			K[0][29] = -0.441513370350e-2;
			K[0][30] = 0.133722924858e-2;
			K[0][31] = 0.264832491957e-3;
			K[0][32] = 0.196688194015e2;
			K[0][33] = -0.209115600730e2;
			K[0][34] = 0.167788306989e-1;
			K[0][35] = 0.262767566274e4;

			K[1][0] = 1;
			K[1][1] = 1;
			K[1][2] = 2;
			K[1][3] = 2;
			K[1][4] = 3;
			K[1][5] = 3;
			K[1][6] = 1;
			K[1][7] = 1;
			K[1][8] = 1;
			K[1][9] = 3;
			K[1][10] = 3;
			K[1][11] = 4;
			K[1][12] = 6;
			K[1][13] = 6;
			K[1][14] = 7;
			K[1][15] = 7;
			K[1][16] = 8;
			K[1][17] = 8;
			K[1][18] = 1;
			K[1][19] = 2;
			K[1][20] = 3;
			K[1][21] = 4;
			K[1][22] = 5;
			K[1][23] = 8;
			K[1][24] = 4;
			K[1][25] = 5;
			K[1][26] = 5;
			K[1][27] = 8;
			K[1][28] = 3;
			K[1][29] = 5;
			K[1][30] = 6;
			K[1][31] = 9;
			K[1][32] = 1;
			K[1][33] = 1;
			K[1][34] = 3;
			K[1][35] = 2;

			K[2][00] = 0.25;
			K[2][01] = 0.875;
			K[2][02] = 0.5;
			K[2][03] = 0.875;
			K[2][04] = 0.375;
			K[2][05] = 0.75;
			K[2][6] = 0.5;
			K[2][7] = 0.75;
			K[2][8] = 2;
			K[2][9] = 1.25;
			K[2][10] = 3.5;
			K[2][11] = 1;

			K[2][12] = 0.5;
			K[2][13] = 3;
			K[2][14] = 0;
			K[2][15] = 2.75;
			K[2][16] = 0.75;
			K[2][17] = 2.5;
			K[2][18] = 4;
			K[2][19] = 6;
			K[2][20] = 6;
			K[2][21] = 3;
			K[2][22] = 3;
			K[2][23] = 6;
			K[2][24] = 16;
			K[2][25] = 11;
			K[2][26] = 15;
			K[2][27] = 12;
			K[2][28] = 12;
			K[2][29] = 7;
			K[2][30] = 4;
			K[2][31] = 16;
			K[2][32] = 0;
			K[2][33] = 1;
			K[2][34] = 2;
			K[2][35] = 3;

			K[3][00] = 0;
			K[3][01] = 0;
			K[3][02] = 0;
			K[3][03] = 0;
			K[3][04] = 0;
			K[3][05] = 0;
			K[3][06] = 1;
			K[3][07] = 1;
			K[3][8] = 1;
			K[3][9] = 1;
			K[3][10] = 1;
			K[3][11] = 1;
			K[3][12] = 1;
			K[3][13] = 1;
			K[3][14] = 1;
			K[3][15] = 1;
			K[3][16] = 1;
			K[3][17] = 1;
			K[3][18] = 2;
			K[3][19] = 2;
			K[3][20] = 2;
			K[3][21] = 2;
			K[3][22] = 2;
			K[3][23] = 2;
			K[3][24] = 3;
			K[3][25] = 3;
			K[3][26] = 3;
			K[3][27] = 3;
			K[3][28] = 4;
			K[3][29] = 4;
			K[3][30] = 4;
			K[3][31] = 4;
			K[3][32] = 2;
			K[3][33] = 2;
			K[3][34] = 2;
			K[3][35] = 2;

			K[10][32] = 20;
			K[10][33] = 20;
			K[10][34] = 15;
			K[10][35] = 25;
			K[11][32] = 325;
			K[11][33] = 325;
			K[11][34] = 300;
			K[11][35] = 275;
			K[12][32] = 1.16;
			K[12][33] = 1.16;
			K[12][34] = 1.13;
			K[12][35] = 1.25;
			K[13][32] = 1;
			K[13][33] = 1;
			K[13][34] = 1;
			K[13][35] = 1;
			break;
		}

		case 'H': // Hydrogen; BG, 03/2012
		{
			fluid_id = 4;
			rhoc = 313.6; // [kg/m3]
			Tc = 33.30; // [K]
			pc = 1297000; // [Pa]
			omega = -0.215; // azentric factor, see PREOS
			molar_mass = 2.015894;
			// Limits sums in FHE-derivations
			break;
		}
		case 'O': // Oxygen
		{
			fluid_id = 5;
			// rhoc = 313.6;             // [kg/m3]
			// Tc = 33.30;               // [K]
			// pc = 1297000;             // [Pa]
			// omega = -0.215;           // azentric factor, see PREOS
			// molar_mass = 2.015894;
			// Limits sums in FHE-derivations
			break;
		}

		default:
			cout << "Error in eos.cpp: no fluid name specified!"
			     << "\n";
			break;
	}
}

//****************************************************************************
//* task: find compressibility factor of a mixture component
//*  note: P in bar! (1e5 Pa)
//* Programming: NB, Sep10
//*****************************************************************************/
// double DuanMixCompressibility(double T,double P,double V,CVirialCoefficients w)
double DuanMixCompressibility(double T, double P, double V, VirialCoefficients w)
{
	double R = 83.14472;
	return P * V / R / T - 1 - w.B / V - w.C / (V * V) - w.D / (V * V * V * V) - w.E / (V * V * V * V * V)
	       - w.F / (V * V) * (w.b + w.G / (V * V)) * exp(-w.G / (V * V));
}

//****************************************************************************
//* returns the third root of a number x, -inf < x < inf
//* Programming: NB, Sep10
//*****************************************************************************/
inline double W3(double x)
{
	if (x < 0)
		return -pow(fabs(x), 1. / 3.);
	else
		return pow(x, 1. / 3.);
}

//****************************************************************************
//* binary mixing coefficient for water and co2 , Duan 1992
//* Programming: NB, Sep10
//*****************************************************************************/
double k_co2_h20(int number, double T)
{
	if (T < 373.15)
	{
		switch (number)
		{
			case 1:
				return 0.20611 + 0.0006 * T;
			case 2:
				return 0.8023278 - 0.0022206 * T + 184.76824 / T;
			case 3:
				return 1.80544 - 0.0032605 * T;
			default:
				return 1;
		}
	}
	else if (T > 673.15)
	{
		switch (number)
		{
			case 1:
				return 3.131 - 5.0624e-03 * T + 1.8641e-06 * T * T - 31.409 / T;
			case 2:
				return -46.646 + 4.2877e-02 * T - 1.0892e-05 * T * T + 1.5782e+04 / T;
			case 3:
				return 0.9;
			default:
				return 1;
		}
	}
	else if (T < 495.15)
	{
		switch (number)
		{
			case 1:
				return -10084.5042 - 4.27134485 * T + 256477.783 / T + 0.00166997474 * T * T + 1816.78 * log(T);
			case 2:
				return 9.000263 - 0.00623494 * T - 2307.7125 / T;
			case 3:
				return -74.1163 + 0.1800496 * T - 1.40904946e-4 * T * T + 101305246 / T;
			default:
				return 1;
		}
	}
	else
		switch (number)
		{
			case 1:
				return -0.3568 + 7.8888e-4 * T + 333.399 / T;
			case 2:
				return -19.97444 + 0.0192515 * T + 5707.4229 / T;
			case 3:
				return 12.1308 - 0.0099489 * T - 3042.09583 / T;
			default:
				return 1;
		}
	cout << " This text should not be printed, something is wrong!"
	     << "\n";
	return 0;
}

//****************************************************************************
//* binary mixing coefficient for co2 and methane , Duan 1992
//* Programming: NB, Sep10
//*****************************************************************************/
double k_co2_ch4(int number, double T)
{
	if (T < 304)
	{
		switch (number)
		{
			case 1:
				return 0.38;
			case 2:
				return 1.74094 - 0.0058903 * T;
			case 3:
				return 1.59;
			default:
				return 1;
		}
	}
	else if (T > 498)
		return 1;
	else
	{
		switch (number)
		{
			case 1:
				return 1.1;
			case 2:
				return 3.211 - 0.00158 * T - 537.814 / T;
			// case 3: return -0.7;
			case 3:
				return 1.80544 - 0.0032605 * T; // aus F90
				// case 3: return 12.1308 -0.0099489*T - 3042.09583/T; aus F90
		}
	}
	cout << " This text should not be printed, something is wrong!"
	     << "\n";
	return 0;
}

//****************************************************************************
//* binary mixing coefficient for water and methane , Duan 1992
//* Programming: NB, Sep10
//*****************************************************************************/
double k_ch4_h2o(int number, double T)
{
	(void)number;
	(void)T;
	// Not implemented yet! If you feel constrained to change this, you'll find the correlation in Duan, Moller and
	// Weare ,1992.
	return 1;
}

//****************************************************************************
//* Parameters for Water, CO2 and Methane for Duan EOS
//*
//* Programming: NB, Sep10
//*****************************************************************************/
void DuansParameter(int fluid, double a[15], double* Tc, double* Pc, double* M)
{
	switch (fluid)
	{
		case 0: // CO2, Duan 1992
			*Tc = 304.1282;
			*Pc = 73.77300; // bar
			*M = 44.01; // g/mol
			a[0] = 8.99288497e-2;
			a[1] = -4.94783127e-1;
			a[2] = 4.77922245e-2;
			a[3] = 1.03808883e-2;
			a[4] = -2.82516861e-2;
			a[5] = 9.49887563e-2;
			a[6] = 5.20600880e-4;
			a[7] = -2.93540971e-4;
			a[8] = -1.77265112e-3;
			a[9] = -2.51101973e-5;
			a[10] = 8.93353441e-5;
			a[11] = 7.88998563e-5;
			a[12] = -1.66727022e-2;
			a[13] = 1.398;
			a[14] = 2.96e-2;
			break;
		case 1: // H2O Duan 1992
			*Tc = 647.25;
			*Pc = 221.19000; // bar
			*M = 18.01; // g/mol
			a[0] = 8.64449220E-02;
			a[1] = -3.96918955E-01;
			a[2] = -5.73334886E-02;
			a[3] = -2.93893000E-04;
			a[4] = -4.15775512E-03;
			a[5] = 1.99496791E-02;
			a[6] = 1.18901426E-04;
			a[7] = 1.55212063E-04;
			a[8] = -1.06855859E-04;
			a[9] = -4.93197687E-06;
			a[10] = -2.73739155E-06;
			a[11] = 2.65571238E-06;
			a[12] = 8.96079018E-03;
			a[13] = 4.02000000E+00;
			a[14] = 2.57000000E-02;
			break;
		case 2: // CH4 Duan 1992
			*Tc = 190.6;
			*Pc = 46.41000; // bar
			*M = 16.04; // g/mol
			a[0] = 8.72553928E-02;
			a[1] = -7.52599476E-01;
			a[2] = 3.75419887E-01;
			a[3] = 1.07291342E-02;
			a[4] = 5.49626360E-03;
			a[5] = -1.84772802E-02;
			a[6] = 3.18993183E-04;
			a[7] = 2.11079375E-04;
			a[8] = 2.01682801E-05;
			a[9] = -1.65606189E-05;
			a[10] = 1.19614546E-04;
			a[11] = -1.08087289E-04;
			a[12] = 4.48262295E-02;
			a[13] = 7.5397E-01;
			a[14] = 7.7167E-02;
			break;
		case 5: // h20 Duan 2006
			*Tc = 647.25;
			*Pc = 22119000; // bar
			*M = 18.01; // g/mol
			a[0] = 4.38269941e-02;
			a[1] = -1.68244362e-01;
			a[2] = -2.36923373e-01;
			a[3] = 1.13027462e-02;
			a[4] = -7.67764181e-02;
			a[5] = 9.71820593e-02;
			a[6] = 6.62674916e-05;
			a[7] = 1.06637349e-03;
			a[8] = -1.23265258e-03;
			a[9] = -8.93953948e-06;
			a[10] = -3.88124606e-05;
			a[11] = 5.61510206e-05;
			a[12] = 7.51274488e-03;
			a[13] = 2.51598931e+00;
			a[14] = 3.94000000e-02;
			break;
		case 6: // co2 Duan 2006
			*Tc = 304.1282;
			*Pc = 7377300; // bar
			*M = 18.01; // g/mol
			a[0] = 1.14400435e-01;
			a[1] = -9.38526684e-01;
			a[2] = 7.21857006e-01;
			a[3] = 8.81072902e-03;
			a[4] = 6.36473911e-02;
			a[5] = -7.70822213e-02;
			a[6] = 9.01506064e-04;
			a[7] = -6.81834166e-03;
			a[8] = 7.32364258e-03;
			a[9] = -1.10288237e-04;
			a[10] = 1.26524193e-03;
			a[11] = -1.49730823e-03;
			a[12] = 7.81940730e-03;
			a[13] = -4.22918013e+00;
			a[14] = 1.58500000e-01;
			break;

		default:
			break;
	}
}

//****************************************************************************
//* this function calls the binary mixing parameter for a certain condition
//*
//* Programming: NB, Sep10
//*****************************************************************************/
double bip(int number, int fluid_a, int fluid_b, int i, int j, int k, double T)
{
	// number : k1, k2 or k3
	// fluid_a, fluid_b: a-b, b-c, a-c have different binary interaction parameters
	// i,j,k  : needed for mixing rule
	// T : bip depends on Temperature

	if (((number == 1) && (i == j)) || (((number == 2) || (number == 3)) && ((i == j) && (i == k))))
		return 1;
	else
		switch (fluid_a + fluid_b)
		{
			case (1): // 2+3=5 --> CO2+H2O (duan 1992) or 0+1=1 --> CO2+H20 fo high pressure range, duan 2006
				return k_co2_h20(number, T);
			case (5): // 2+3=5 --> CO2+H2O (duan 1992) or 0+1=1 --> CO2+H20 fo high pressure range, duan 2006
				return k_co2_h20(number, T);

			case 6: // 2+4=6 --> CO2+CH4
				return k_co2_ch4(number, T);
			case 7: // 3+4=7 --> CH4+H2O
				return k_co2_ch4(number, T);
			default:
				return 1;
		}
	cout << " This text should not be printed, something is wrong!"
	     << "\n";
	return 0;
}

//****************************************************************************
//* task: calculate the virial coefficients of Duans EOS for a certain fluid at
//* a temperature T
//* Programming: NB, Sep10
//*****************************************************************************/
// CVirialCoefficients DuansVirialCoefficients(int Fluid, double T)
VirialCoefficients DuansVirialCoefficients(int Fluid, double T)
{
	double a[15] = {};
	double Tc = 0.0, Pc = 0.0, M = 0.0;
	double R = 83.14467; // cm��?bar/(K* mol)
	// CVirialCoefficients x;	//BG
	VirialCoefficients x; // BG
	DuansParameter(Fluid, a, &Tc, &Pc, &M);
	double Tr = T / Tc;

	x.B = a[0] + a[1] / (Tr * Tr) + a[2] / (Tr * Tr * Tr);
	x.C = a[3] + a[4] / (Tr * Tr) + a[5] / (Tr * Tr * Tr);
	x.D = a[6] + a[7] / (Tr * Tr) + a[8] / (Tr * Tr * Tr);
	x.E = a[9] + a[10] / (Tr * Tr) + a[11] / (Tr * Tr * Tr);
	x.F = a[12] / (Tr * Tr * Tr);
	x.b = a[13];
	x.G = a[14];
	x.Tc = Tc;
	x.Pc = Pc;
	x.Vc = R * Tc / Pc;
	x.M = M;
	x.id = Fluid;

	return x;
}

//****************************************************************************
//* this function mixes all virial coefficients of two fluids for DUAN EOS with
//* respect to temperature and mole fraction
//*
//* Input: fluid number A and B (see DuanParameter() for respective numbers)
//*              T:   Temparature
//*        x_0: mole fraction of fluid A
//*
//* Programming: NB, Sep10
//*****************************************************************************/
// void MixDuansVirialCoefficients(CVirialCoefficients fluid_a, CVirialCoefficients fluid_b, double T, double x_0,
// CVirialCoefficients*mix)
void MixDuansVirialCoefficients(
    VirialCoefficients fluid_a, VirialCoefficients fluid_b, double T, double x_0, VirialCoefficients* mix)
{
	// this code may look weird, I tried to avoid functions like pow() to increase speed

	double Bij, Cij, Dij, Eij, Fij, Gij;
	double Vcij;
	double B[2], C[2], D[2], E[2], F[2], G[2], bi[2], MM[2];
	double Vc[2];
	double x[2];

	double BVc = 0;
	double CVc = 0;
	double DVc = 0;
	double EVc = 0;
	double FVc = 0;
	double b = 0;
	double GVc = 0;

	double Vcr = 0;
	double M = 0;

	double V2, V3, V5, V6;

	Vc[0] = fluid_a.Vc;
	Vc[1] = fluid_b.Vc;
	MM[0] = fluid_a.M;
	MM[1] = fluid_b.M;

	B[0] = fluid_a.B;
	B[1] = fluid_b.B;
	C[0] = fluid_a.C;
	C[1] = fluid_b.C;
	D[0] = fluid_a.D;
	D[1] = fluid_b.D;
	E[0] = fluid_a.E;
	E[1] = fluid_b.E;
	F[0] = fluid_a.F;
	F[1] = fluid_b.F;
	bi[0] = fluid_a.b;
	bi[1] = fluid_b.b;
	G[0] = fluid_a.G;
	G[1] = fluid_b.G;

	x[0] = x_0;
	x[1] = 1 - x[0];

	for (int i = 0; i < 2; i++)
	{
		b += x[i] * bi[i];
		Vcr += x[i] * Vc[i];
		M += x[i] * MM[i];
		for (int j = 0; j < 2; j++)
		{
			double B1 = ((W3(B[i]) + W3(B[j])) / 2);
			double F1 = ((W3(F[i]) + W3(F[j])) / 2);
			V2 = ((W3(Vc[i]) + W3(Vc[j])) / 2);

			Bij = B1 * B1 * B1 * bip(1, fluid_a.id, fluid_b.id, i, j, 0, T);
			Fij = F1 * F1 * F1;
			Vcij = V2 * V2 * V2;
			BVc += x[i] * x[j] * Bij * Vcij;
			FVc += x[i] * x[j] * Fij * Vcij * Vcij;
			for (int k = 0; k < 2; k++)
			{
				double C1 = ((W3(C[i]) + W3(C[j]) + W3(C[k])) / 3);
				double G1 = ((W3(G[i]) + W3(G[j]) + W3(G[k])) / 3);
				V3 = ((W3(Vc[i]) + W3(Vc[j]) + W3(Vc[k])) / 3);
				Cij = C1 * C1 * C1 * bip(2, fluid_a.id, fluid_b.id, i, j, k, T);
				Gij = G1 * G1 * G1 * bip(3, fluid_a.id, fluid_b.id, i, j, k, T);
				Vcij = V3 * V3 * V3;
				CVc += x[i] * x[j] * x[k] * Cij * Vcij * Vcij;
				GVc += x[i] * x[j] * x[k] * Gij * Vcij * Vcij;
				for (int l = 0; l < 2; l++)
					for (int m = 0; m < 2; m++)
					{
						double D1 = ((W3(D[i]) + W3(D[j]) + W3(D[k]) + W3(D[l]) + W3(D[m])) / 5);
						V5 = ((W3(Vc[i]) + W3(Vc[j]) + W3(Vc[k]) + W3(Vc[l]) + W3(Vc[m])) / 5);
						Dij = D1 * D1 * D1;
						Vcij = V5 * V5 * V5;
						DVc += x[i] * x[j] * x[k] * x[l] * x[m] * Dij * Vcij * Vcij * Vcij * Vcij;
						for (int n = 0; n < 2; n++)
						{
							double E1 = ((W3(E[i]) + W3(E[j]) + W3(E[k]) + W3(E[l]) + W3(E[m]) + W3(E[n])) / 6);
							V6 = ((W3(Vc[i]) + W3(Vc[j]) + W3(Vc[k]) + W3(Vc[l]) + W3(Vc[m]) + W3(Vc[n])) / 6);
							Eij = E1 * E1 * E1;
							Vcij = V6 * V6 * V6;
							EVc += x[i] * x[j] * x[k] * x[l] * x[m] * x[n] * Eij * Vcij * Vcij * Vcij * Vcij * Vcij;
						} // ijklmn
					} // ijklm
			} // ijk
		} // ij
	} // i

	mix->B = BVc;
	mix->C = CVc;
	mix->D = DVc;
	mix->E = EVc;
	mix->F = FVc;
	mix->b = b;
	mix->G = GVc;
	mix->M = M;
	mix->Vc = Vcr;
}

//****************************************************************************
//* task: returns the density of a mixture of two gases (fluid1 and fluid2),
//*  where:
//*     T : Temperature [K]
//*   P : Pressure [Pa]
//*   x : mole fraction of fluid1
//*
//*  Programming: NB, Sep10
//*****************************************************************************/
double DuansMixingRule(double T, double P, double x, int fluid1, int fluid2, bool neu)
{
	(void)neu; // unused
	P /= 1e5;
	// CVirialCoefficients u,v,w;
	VirialCoefficients u, v, w;

	u = DuansVirialCoefficients(fluid1, T);
	v = DuansVirialCoefficients(fluid2, T);
	MixDuansVirialCoefficients(u, v, T, x, &w);

	double V1, V2;
	int n = 0;
	double V;

	double R = 83.14467;

	double dev_1, dev_2, dev;

	V1 = 6.0;
	dev_1 = DuanMixCompressibility(T, P, V1, w);
	V = V1;

	// while (true) {
	//	V = V + 2.0;
	//	n++;
	//	dev = DuansCompressibility(T,P,V,w);

	//	if (V<=50){
	//		if (dev <= dev_1)
	//		{
	//			V1=V;
	//			dev_1=dev;
	//		} }else
	//		{
	//			//cout << " loop 1 left after " << n << " circles "<< "\n";
	//			break;
	//		}
	//}*/

	for (int i = (int)V1; i < 51; i += 2)
	{
		V += 2;
		dev = DuanMixCompressibility(T, P, V, w);
		if (dev <= dev_1)
		{
			V1 = V;
			dev_1 = dev;
		}
	}

	// Find the first maximum point as V_2
	V2 = V1;
	dev_2 = DuanMixCompressibility(T, P, V2, w);
	V = V2;

	n = 0;
	while (true)
	{
		n++;
		V = V + 2.0;
		dev = DuanMixCompressibility(T, P, V, w);

		if (V <= max(R * T / P, 1000.0))
		{
			if (dev > dev_2)
			{
				V2 = V;
				dev_2 = dev;
			}
			else
				// cout << " loop 2 left after " << n << " circles "<< "\n";
				break;
			// goto l400;
		}
		else
			// cout << " loop 2 left after " << n << " circles "<< "\n";
			break;
	}

	// Get the solution of volume with divition method
	n = 0;
	while (true)
	{
		dev_1 = DuanMixCompressibility(T, P, V1, w);
		dev_2 = DuanMixCompressibility(T, P, V2, w);

		if ((dev_1 * dev_2) > 0.0)
			V2 = V2 + 1;

		else
		{
			V = (V1 + V2) / 2.0;

			dev = DuanMixCompressibility(T, P, V, w);

			if (fabs(dev) < 1.0e-5)
				// cout << " returning after " << n << " iterations " << "\n";
				return w.M / V * 1000;

			if ((dev_1 * dev) > 0.0)
				V1 = V;
			else if (dev_2 * dev > 0.0)
				V2 = V;
			n++;
			if (n > 1000)
			{
				cout << " max it reached "
				     << "\n";
				return -1; // Max it
			}
		}
	}

	cout << "This should not happen. Call NB at TUD!"
	     << "\n";
	return -1;
}
