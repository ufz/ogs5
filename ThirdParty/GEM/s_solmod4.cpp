//-------------------------------------------------------------------
// $Id: s_fgl1.cpp 724 2012-10-02 14:25:25Z kulik $
//
/// \file s_solmod4.cpp
/// Implementation of TSolMod derived classes for ion-association aqueous
/// activity models (THelgeson, TDavies, TLimitingLaw, TDebyeHueckel,
/// TKarpov, TShvarov)
//
// Copyright (c) 2008-2012  T.Wagner, S.Dmitrieva, D.Kulik
// <GEMS Development Team, mailto:gems2.support@psi.ch>
//
// This file is part of the GEMS3K code for thermodynamic modelling
// by Gibbs energy minimization <http://gems.web.psi.ch/GEMS3K/>
//
// GEMS3K is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.

// GEMS3K is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with GEMS3K code. If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------
//

#include <cmath>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;
#include "verror.h"
#include "s_solmod.h"





//=============================================================================================
// Extended Debye-Hueckel (EDH) model for aqueous solutions, Helgeson version
// References: Helgeson et al. (1981); Oelkers and Helgeson (1990);
// Pokrovskii and Helgeson (1995; 1997a; 1997b)
// (c) TW May 2009
//=============================================================================================


// Generic constructor for the THelgeson class
THelgeson::THelgeson( SolutionData *sd, double *arM, double *arZ, double *dW, double *eW ):
                TSolMod( sd )
{
    alloc_internal();
    m = arM;
    z = arZ;
    RhoW = dW;
    EpsW = eW;
    ac = aIPc[1];   // common ion size parameter
    bc = aIPc[0];   // common b_gamma
    flagNeut = (long int)aIPc[2];   // 0: unity, 1: calculated
    flagH2O = (long int)aIPc[3];   // 0: unity, 1: calculated
    flagElect = (long int)aIPc[4];  // 0: constant, 1: NaCl, 2: KCl, 3: NaOH, 4: KOH
}


THelgeson::~THelgeson()
{
	free_internal();
}


void THelgeson::alloc_internal()
{
	LnG = new double [NComp];
	dLnGdT = new double [NComp];
	d2LnGdT2 = new double [NComp];
	dLnGdP = new double [NComp];
}


void THelgeson::free_internal()
{
	delete[]LnG;
	delete[]dLnGdT;
	delete[]d2LnGdT2;
	delete[]dLnGdP;
}


/// Calculates T,P corrected parameters
long int THelgeson::PTparam()
{
	double alp, bet, dal, rho, eps, dedt, d2edt2, dedp;

	// read and convert rho and eps
	rho = RhoW[0];
	alp = - 1./rho*RhoW[1];
	dal = pow(alp,2.) - 1./rho*RhoW[2];
	bet = 1./rho*RhoW[3];
	eps = EpsW[0];
	dedt = 1./eps*EpsW[1];
	d2edt2 = - 1./pow(eps,2.)*pow(EpsW[1],2.) + 1./eps*EpsW[2];
	dedp = 1./eps*EpsW[3];

	// calculate A and B terms of Debye-Huckel equation (and derivatives)
	A = (1.82483e6)*sqrt(rho) / pow(Tk*eps,1.5);
	dAdT = - 3./2.*A*( dedt + 1./Tk + alp/3. );
	d2AdT2 = 1./A*pow(dAdT,2.) - 3./2.*A*( d2edt2 - 1/pow(Tk,2.) + 1/3.*dal );
	dAdP = 1./2.*A*( bet - 3.*dedp);
	B = (50.2916)*sqrt(rho) / sqrt(Tk*eps);
	dBdT = - 1./2.*B*( dedt +1./Tk + alp );
	d2BdT2 = 1./B*pow(dBdT,2.) - 1./2.*B*( d2edt2 - 1./pow(Tk,2.) + dal );
	dBdP = 1./2.*B*( bet - dedp );

 	ErrorIf( fabs(A) < 1e-9 || fabs(B) < 1e-9, "Helgeson EDH model",
 			"Error: DH parameter A or B was not calculated - no values of RhoW and EpsW !" );

	// b_gamma constant
	if ( flagElect == 0)
	{
		bgam = bc;
		dbgdT = 0.0;
		d2bgdT2 = 0.0;
		dbgdP = 0.0;
		ao = ac;
		daodT = 0.0;
		d2aodT2 = 0.0;
		daodP = 0.0;
	}

	// b_gamma TP-dependent
	else
	{
		Gfunction();
		BgammaTP();
		IonsizeTP();
		// ao = ac;
		// daodT = 0.0;
		// d2aodT2 = 0.0;
		// daodP = 0.0;
	}

	return 0;
}


/// Calculates activity coefficients
long int THelgeson::MixMod()
{
	long int j, k, w;
	double sqI, Z2, lgGam, lnGam, Nw, Lgam, lnxw, xw, Lam, sig,
			Phi, Phit, zc, za, psi, lnaw, lg_to_ln;
	zc = 1.; za = 1.; psi = 1.; lg_to_ln = 2.302585093;

	// get index of water (assumes water is last species in phase)
	w = NComp - 1;

	// calculate ionic strength and total molaities
	IonicStrength();

	xw = x[w];
	Nw = 1000./18.01528;
	// Lgam = -log10(1.+0.0180153*molT);
	Lgam = log10(xw);  // Helgeson large gamma simplified
	if( Lgam < -0.7 )
		Lgam = -0.7;  // experimental truncation of Lgam to min ln(0.5)
	lnxw = log(xw);
	sqI = sqrt(IS);

	// loop over species
	for( j=0; j<NComp; j++ )
	{
		lgGam = 0.0;
		lnGam = 0.0;

		// charged species
		if ( z[j] )
		{
			lgGam = 0.0;
			Z2 = z[j]*z[j];
			lgGam = ( - A * sqI * Z2 ) / ( 1. + B * ao * sqI ) + bgam * IS ;
			lnGamma[j] = (lgGam + Lgam) * lg_to_ln;
		}

		// neutral species and water solvent
		else
		{
			// neutral species
			if ( j != (NComp-1) )
			{
				lgGam = 0.0;
				if ( flagNeut == 1 )
					lgGam = bgam * IS;
				else
					lgGam = 0.0;
				lnGamma[j] = (lgGam + Lgam) * lg_to_ln;
				continue;
			}

			// water solvent
			else
			{
				lgGam = 0.0;
				lnGam = 0.0;
				if ( flagH2O == 1 )
				{
					Phit = 0.0;
					// Phi corrected using eq. (190) from Helgeson et al. (1981)
					Lam = 1. + (ao*B) * sqI;
					sig = 3./(pow(ao,3.)*pow(B,3.)*pow(IS,(3./2.))) * (Lam-1./Lam-2*log(Lam));
					// Phi = -2.3025851*(A*sqI*SigTerm/3. + Lgam/(0.0180153*2.*IS) - bgam*IS/2.);
					// Phi = - log(10.) * (molZ/molT) * ( (zc*za*A*sqI*SigTerm)/3. + (Lgam*psi)/(0.0180153*2.*IS) - bgam*IS/2. );

					for (k=0; k<(NComp-1); k++)
					{
						if ( (z[k] == 0) && (flagNeut == 0) )
							Phit += - log(10.) * m[k] * ( (pow(z[k],2.)*A*sqI*sig)/3. + Lgam/(0.0180153*molT) );
						else
							Phit += - log(10.) * m[k] * ( (pow(z[k],2.)*A*sqI*sig)/3. + Lgam/(0.0180153*molT) - bgam*IS/2. );
					}

					Phi = Phit/molT;
					lnaw = - Phi*molT/Nw;
					lnGam = lnaw - lnxw;
				}
				else
					lnGam = 0.0;
				lnGamma[j] = lnGam;
			}
		}
	}  // j

	return 0;
}


/// calculates excess properties
long int THelgeson::ExcessProp( double *Zex )
{
	long int j, k, w;
	double sqI, Z2, Nw, Lgam, lnxw, xw, Lam, sig, Phi, dPhidT, d2PhidT2, dPhidP,
			Phit, dPhitdT, d2PhitdT2, dPhitdP, lnaw, lg_to_ln, g, dgt, d2gt, dgp;
	double U, V, dUdT, dVdT, d2UdT2, d2VdT2, dUdP, dVdP, U1, U2, U3, V1, V2, V3,
			dU1dT, dU2dT, dU3dT, dV1dT, dV2dT, dV3dT, d2U1dT2, d2U2dT2, d2U3dT2,
			d2V1dT2, d2V2dT2, d2V3dT2, dU1dP, dU2dP, dU3dP, dV1dP, dV2dP, dV3dP,
			L, dLdT, d2LdT2, dLdP, Z, dZdT, d2ZdT2, dZdP;
	lg_to_ln = 2.302585093;
	g = 0.; dgt = 0.; d2gt = 0.; dgp = 0.;

	// get index of water (assumes water is last species in phase)
	w = NComp - 1;

	// calculate ionic strength and total molalities
	IonicStrength();

	xw = x[w];
	Nw = 1000./18.01528;
	// Lgam = -log10(1.+0.0180153*molT);
	Lgam = log10(xw);  // Helgeson large gamma simplified
	if( Lgam < -0.7 )
		Lgam = -0.7;  // experimental truncation of Lgam to min ln(0.5)
	lnxw = log(xw);
	sqI = sqrt(IS);

	// loop over species
	for( j=0; j<NComp; j++ )
	{
		// charged species
		if ( z[j] )
		{
			Z2 = z[j]*z[j];
			U = - (Z2*A) * sqI;
			dUdT = - (Z2*dAdT) * sqI;
			d2UdT2 = - (Z2*d2AdT2) * sqI;
			dUdP = - (Z2*dAdP) * sqI;
			V = 1. + (ao*B) * sqI;
			dVdT = ( daodT*B + ao*dBdT ) * sqI;
			d2VdT2 = ( d2aodT2*B + 2.*daodT*dBdT + ao*d2BdT2 ) * sqI;
			dVdP = ( daodP*B + ao*dBdP ) * sqI;
			LnG[j] = ( U/V + bgam*IS ) * lg_to_ln;
			dLnGdT[j] = ( (dUdT*V - U*dVdT)/pow(V,2.) + dbgdT*IS ) * lg_to_ln;
			d2LnGdT2[j] = ( (d2UdT2*V + dUdT*dVdT)/pow(V,2.) - (dUdT*V)*(2.*dVdT)/pow(V,3.)
				- (dUdT*dVdT + U*d2VdT2)/pow(V,2.) + (U*dVdT)*(2.*dVdT)/pow(V,3.)
				+ d2bgdT2*IS ) * lg_to_ln;
			dLnGdP[j] = ( (dUdP*V - U*dVdP)/pow(V,2.) + dbgdP*IS ) * lg_to_ln;
		}

		// neutral species and water solvent
		else
		{
			// neutral species
			if ( j != (NComp-1) )
			{
				if ( flagNeut == 1 )
				{
					LnG[j] = ( bgam*IS ) * lg_to_ln;
					dLnGdT[j] = ( dbgdT*IS ) * lg_to_ln;
					d2LnGdT2[j] = ( d2bgdT2*IS ) * lg_to_ln;
					dLnGdP[j] = ( dbgdP*IS ) * lg_to_ln;
				}

				else
				{
					LnG[j] = 0.;
					dLnGdT[j] = 0.;
					d2LnGdT2[j] = 0.;
					dLnGdP[j] = 0.;
				}
				continue;
			}

			// water solvent
			else
			{
				// H2O activity coefficient calculated
				if ( flagH2O == 1 )
				{
					Phit = 0.; dPhitdT = 0.; d2PhitdT2 = 0.; dPhitdP = 0.;

					// derivatives of lambda and sigma terms
					L = 1. + (ao*B) * sqI;
					dLdT = ( daodT*B + ao*dBdT ) * sqI;
					d2LdT2 = ( d2aodT2*B + 2.*daodT*dBdT + ao*d2BdT2 ) * sqI;
					dLdP = ( daodP*B + ao*dBdP ) * sqI;

					U1 = (A*L);
					dU1dT = (dAdT*L + A*dLdT);
					d2U1dT2 = ( d2AdT2*L + 2.*dAdT*dLdT + A*d2LdT2 );
					dU1dP = ( dAdP*L + A*dLdP );
					V1 = pow(ao,3.)*pow(B,3.) * IS;
					dV1dT = ( 3.*pow(ao,2.)*daodT*pow(B,3.) + 3.*pow(ao,3.)*pow(B,2.)*dBdT ) * IS;
					d2V1dT2 = ( 6.*ao*pow(daodT,2.)*pow(B,3.) + 3.*pow(ao,2.)*d2aodT2*pow(B,3.)
								+ 18.*pow(ao,2.)*daodT*pow(B,2.)*dBdT + 6.*pow(ao,3.)*B*pow(dBdT,2.)
								+ 3.*pow(ao,3.)*pow(B,2.)*d2BdT2 ) * IS;
					dV1dP = ( 3.*pow(ao,2.)*daodP*pow(B,3.) + 3.*pow(ao,3.)*pow(B,2.)*dBdP ) * IS;

					U2 = A;
					dU2dT = dAdT;
					d2U2dT2 = d2AdT2;
					dU2dP = dAdP;
					V2 = pow(ao,3.)*pow(B,3.)*L * IS;
					dV2dT = ( 3.*pow(ao,2.)*daodT*pow(B,3.)*L + 3.*pow(ao,3.)*pow(B,2.)*dBdT*L
								+ pow(ao,3.)*pow(B,3.)*dLdT ) * IS;
					d2V2dT2 = ( 6.*ao*pow(daodT,2.)*pow(B,3.)*L + 3.*pow(ao,2.)*d2aodT2*pow(B,3.)*L
								+ 18.*pow(ao,2.)*daodT*pow(B,2.)*dBdT*L + 6.*pow(ao,2.)*daodT*pow(B,3.)*dLdT
								+ 6.*pow(ao,3.)*B*pow(dBdT,2.)*L + 3.*pow(ao,3.)*pow(B,2.)*d2BdT2*L
								+ 6.*pow(ao,3.)*pow(B,2.)*dBdT*dLdT + pow(ao,3.)*pow(B,3.)*d2LdT2 ) * IS;
					dV2dP = ( 3.*pow(ao,2.)*daodP*pow(B,3.)*L + 3.*pow(ao,3.)*pow(B,2.)*dBdP*L
								+ pow(ao,3.)*pow(B,3.)*dLdP ) * IS;

					U3 = 2.*( A*log(L) );
					dU3dT = 2.*( dAdT*log(L) + A*(1./L)*dLdT );
					d2U3dT2 = 2.*( d2AdT2*log(L) + 2.*dAdT*(1./L)*dLdT
								- A*(1./pow(L,2.))*pow(dLdT,2.) + A*(1./L)*d2LdT2 );
					dU3dP = 2.*( dAdP*log(L) + A*(1./L)*dLdP );
					V3 = pow(ao,3.)*pow(B,3.) * IS;
					dV3dT = ( 3.*pow(ao,2.)*daodT*pow(B,3.) + 3.*pow(ao,3.)*pow(B,2.)*dBdT ) * IS;
					d2V3dT2 = ( 6.*ao*pow(daodT,2.)*pow(B,3.) + 3.*pow(ao,2.)*d2aodT2*pow(B,3.)
								+ 18.*pow(ao,2.)*daodT*pow(B,2.)*dBdT + 6.*pow(ao,3.)*B*pow(dBdT,2.)
								+ 3.*pow(ao,3.)*pow(B,2.)*d2BdT2 ) * IS;
					dV3dP = ( 3.*pow(ao,2.)*daodP*pow(B,3.) + 3.*pow(ao,3.)*pow(B,2.)*dBdP ) * IS;

					Z = U1/V1 - U2/V2 - U3/V3;
					dZdT = (dU1dT*V1 - U1*dV1dT)/pow(V1,2.) - (dU2dT*V2 - U2*dV2dT)/pow(V2,2.)
								- (dU3dT*V3 - U3*dV3dT)/pow(V3,2.);
					d2ZdT2 = (d2U1dT2*V1 + dU1dT*dV1dT)/pow(V1,2.) - (dU1dT*V1)*(2.*dV1dT)/pow(V1,3.)
								- (dU1dT*dV1dT + U1*d2V1dT2)/pow(V1,2.) + (U1*dV1dT)*(2.*dV1dT)/pow(V1,3.)
								- (d2U2dT2*V2 + dU2dT*dV2dT)/pow(V2,2.) + (dU2dT*V2)*(2.*dV2dT)/pow(V2,3.)
								+ (dU2dT*dV2dT + U2*d2V2dT2)/pow(V2,2.) - (U2*dV2dT)*(2.*dV2dT)/pow(V2,3.)
								- (d2U3dT2*V3 + dU3dT*dV3dT)/pow(V3,2.) + (dU3dT*V3)*(2.*dV3dT)/pow(V3,3.)
								+ (dU3dT*dV3dT + U3*d2V3dT2)/pow(V3,2.) - (U3*dV3dT)*(2.*dV3dT)/pow(V3,3.);
					dZdP = (dU1dP*V1 - U1*dV1dP)/pow(V1,2.) - (dU2dP*V2 - U2*dV2dP)/pow(V2,2.)
								- (dU3dP*V3 - U3*dV3dP)/pow(V3,2.);

					// increments to osmotic coefficient (and derivatives)
					Lam = 1. + (ao*B) * sqI;
					sig = 3./(pow(ao,3.)*pow(B,3.)*pow(IS,(3./2.))) * (Lam-1./Lam-2*log(Lam));

					for (k=0; k<(NComp-1); k++)
					{
						if ( (z[k] == 0) && (flagNeut == 0) )
						{
							Phit += - log(10.) * m[k] * ( (pow(z[k],2.)*A*sqI*sig)/3. + Lgam/(0.0180153*molT) );
							dPhitdT  += - log(10.) * m[k] * ( pow(z[k],2.)*dZdT );
							d2PhitdT2 += - log(10.) * m[k] * ( pow(z[k],2.)*d2ZdT2 );
							dPhitdP += - log(10.) * m[k] * ( pow(z[k],2.)*dZdP );
						}

						else
						{
							Phit += - log(10.) * m[k] * ( (pow(z[k],2.)*A*sqI*sig)/3. + Lgam/(0.0180153*molT) - bgam*IS/2. );
							dPhitdT  += - log(10.) * m[k] * ( pow(z[k],2.)*dZdT - dbgdT*IS/2. );
							d2PhitdT2 += - log(10.) * m[k] * ( pow(z[k],2.)*d2ZdT2 - d2bgdT2*IS/2. );
							dPhitdP += - log(10.) * m[k] * ( pow(z[k],2.)*dZdP - dbgdP*IS/2. );
						}
					}

					Phi = Phit/molT;
					dPhidT = dPhitdT/molT;
					d2PhidT2 = d2PhitdT2/molT;
					dPhidP = dPhitdP/molT;

					// activity coefficient (and derivatives)
					lnaw = - Phi*molT/Nw;
					LnG[j] = lnaw - lnxw;
					dLnGdT[j] = - (molT/Nw) * dPhidT;
					d2LnGdT2[j] = - (molT/Nw) * d2PhidT2;
					dLnGdP[j] = - (molT/Nw) * dPhidP;

				}

				// H2O activity coefficient 0
				else
				{
					LnG[j] = 0.;
					dLnGdT[j] = 0.;
					d2LnGdT2[j] = 0.;
					dLnGdP[j] = 0.;
				}
			}
		}

		g += x[j]*LnG[j];
		dgt += x[j]*dLnGdT[j];
		d2gt += x[j]*d2LnGdT2[j];
		dgp += x[j]*dLnGdP[j];

	} // j

	// increment thermodynamic properties
	Gex = (R_CONST*Tk) * g;
	Hex = - R_CONST*pow(Tk,2.) * dgt;
	// Sex = - R_CONST * ( g + Tk*dgt );
	Sex = (Hex-Gex)/Tk;
	CPex = - R_CONST * ( 2.*Tk*dgt + pow(Tk,2.)*d2gt );
	Vex = (R_CONST*Tk) * dgp;
	Aex = Gex - Vex*Pbar;
	Uex = Hex - Vex*Pbar;

	// assigments (excess properties)
	Zex[0] = Gex;
	Zex[1] = Hex;
	Zex[2] = Sex;
	Zex[3] = CPex;
	Zex[4] = Vex;
	Zex[5] = Aex;
	Zex[6] = Uex;

	return 0;
}


/// calculates ideal mixing properties
long int THelgeson::IdealProp( double *Zid )
{
	long int j;
	double si;
	si = 0.0;
	for (j=0; j<NComp; j++)
	{
		if ( x[j] > 1.0e-32 )
			si += x[j]*log(x[j]);
	}
	Hid = 0.0;
	CPid = 0.0;
	Vid = 0.0;
	Sid = (-1.)*R_CONST*si;
	Gid = Hid - Sid*Tk;
	Aid = Gid - Vid*Pbar;
	Uid = Hid - Vid*Pbar;

	// assignments (ideal mixing properties)
	Zid[0] = Gid;
	Zid[1] = Hid;
	Zid[2] = Sid;
	Zid[3] = CPid;
	Zid[4] = Vid;
	Zid[5] = Aid;
	Zid[6] = Uid;
	return 0;
}


/// calculates true ionic strength
long int THelgeson::IonicStrength()
{
	long int j;
	double is, mt, mz;
	is = 0.0; mt = 0.0; mz = 0.0;

	// calculate ionic strength and total molalities
	// needs to be modified when nonaqueous solvents are present
	for (j=0; j<(NComp-1); j++)
	{
		is += 0.5*m[j]*z[j]*z[j];
		mt += m[j];
		if ( z[j] )
			mz += m[j];
	}

	// assignments
	IS = is;
	molT = mt;
	molZ = mz;

	return 0;
}


/// calculates TP dependence of b_gamma (and derivatives)
long int THelgeson::BgammaTP()
{
	// ni: stoichiometric number of moles of ions in one mole of electrolyte
	// rc, ra: radius of cation and anion, respectively at 298 K/1 bar
	// units are cal, kg, K, mol, bar
	double ni, nc, na, zc, za, rc, ra, a1, a2, a3, a4, a5, c1, c2, omg, bg, bs, bh, rec, rea,
			omgpt, domdt, d2omdt2, domdp, nbg, nbv, nbj, nbh;
	double eps, eta, xborn, yborn, qborn, X1, X2;

	// set parameters
	eps = EpsW[0];
	yborn = 1./pow(EpsW[0],2.)*EpsW[1];
	xborn = 1./pow(EpsW[0],2.) * ( EpsW[2] - 2./EpsW[0]*pow(EpsW[1],2.) );
	qborn = 1./pow(EpsW[0],2.)*EpsW[3];
	eta = (1.66027e5);

	switch ( flagElect )
	{
		case 1:  // NaCl
			ni = 2.; nc = 1.; na = 1.; zc = 1.; za = -1.;
			a1 = 0.030056; a2 = -202.55; a3 = -2.9092; a4 = 20302;
			a5 = -0.206; c1 = -1.50; c2 = 53300.; omg = 178650.;
			bg = -174.623; bs = 2.164; rc = 0.97; ra = 1.81;
			break;
		case 2:  // KCl
			ni = 2.; nc = 1.; na = 1.; zc = 1.; za = -1.;
			a1 = 0.0172; a2 = -115.36; a3 = -1.1857; a4 = 13854.2;
			a5 = -0.262; c1 = -2.53; c2 = 38628.4; omg = 164870.;
			bg = -70.0; bs = 1.727; rc = 1.33; ra = 1.81;
			break;
		case 3:  // NaOH
			ni = 2.; nc = 1.; na = 1.; zc = 1.; za = -1.;
			a1 = 0.030056; a2 = -202.55; a3 = -2.9092; a4 = 20302;
			a5 = -0.206; c1 = -1.50; c2 = 53300.; omg = 205520.;
			bg = -267.4; bs = 1.836; rc = 0.97; ra = 1.40;
			break;
		case 4:  // KOH
			ni = 2.; nc = 1.; na = 1.; zc = 1.; za = -1.;
			a1 = 0.0172; a2 = -115.36; a3 = -1.1857; a4 = 13854.2;
			a5 = -0.262; c1 = -2.53; c2 = 38628.4; omg = 191730.;
			bg = -335.7; bs = 1.26; rc = 1.33; ra = 1.40;
			break;
		default:  // wrong mode
			return -1;
	}

	// calculation part, extended 06.06.2009 (TW)
	bh = bg + (298.15)*bs;
	rec = rc + fabs(zc)*(0.94+Gf);
	rea = ra + fabs(za)*Gf;
	X1 = - eta*nc*( fabs(pow(zc,3.))/pow(rec,2.) - zc/pow((3.082+Gf),2.) )
			- eta*na*( fabs(pow(za,3.))/pow(rea,2.) - za/pow((3.082+Gf),2.) );
	X2 = 2.*eta*nc*( fabs(pow(zc,4.))/pow(rec,3.) - zc/pow((3.082+Gf),3.) )
			+ 2.*eta*nc * ( fabs(pow(zc,4.))/pow(rec,3.) - zc/pow((3.082+Gf),3.) );
	omgpt = eta*( nc*pow(zc,2.)/rec + na*pow(za,2.)/rea );
	// omgpt = (1.66027e5)*(1./(0.94+rc+Gf)+1./(ra+Gf));
	domdt = X1*dGfdT;
	d2omdt2 = X2*pow(dGfdT,2.) + X1*d2GfdT2;
	domdp = X1*dGfdP;
	nbg = - ni*bg/2. + ni*bs*(Tk-298.15)/2. - c1*(Tk*log(Tk/298.15)-Tk+298.15)
				+ a1*(Pbar-1.) + a2*log((2600.+Pbar)/(2600.+1.))
				- c2*((1./(Tk-228.)-1./(298.15-228.))*(228.-Tk)/228.-Tk/(228.*228.)
				* log((298.15*(Tk-228.))/(Tk*(298.15-228.))))
				+ 1./(Tk-228.)*(a3*(Pbar-1.) + a4*log((2600.+Pbar)/(2600.+1.)))
				+ a5*(omgpt*(1./eps-1.)-omg*(1./(78.24513795)-1.)+(-5.798650444e-5)*omg*(Tk-298.15));
	nbh = - ni/2.*bh + c1*(Tk-298.15) - c2*(1./(Tk-228.)-1./(298.15-228.))
				+ a1*(Pbar-1.) + a2*log((2600.+Pbar)/(2600.+1.))
				+ ((2.*Tk-228.)/pow((Tk-228.),2.))*(a3*(Pbar-1.)+a4*log((2600.+Pbar)/(2600.+1.)))
				+ a5*( omgpt*(1./eps-1.) + omgpt*Tk*yborn - Tk*(1./eps-1.)*domdt
				- omg*(1./(78.24513795)-1.) - omg*(298.15)*(-5.798650444e-5) );
	nbj = c1 + c2/pow((Tk-228.),2.)  - ( (2.*Tk)/pow((Tk-228.),3.) )
				* ( a3*(Pbar-1.) + a4*log((2600.+Pbar)/(2600.+1.)) )
				+ a5*( omgpt*Tk*xborn + 2.*Tk*yborn*domdt - Tk*(1./eps-1.)*d2omdt2 );
	nbv = a1 + a2/(2600.+Pbar) + a3/(Tk-228.) + a4/((2600.+Pbar)*(Tk-228.))
				+ a5*(-omgpt*qborn + (1./eps-1.)*domdp);
	// b_gamma = nbg/(2.*log(10.)*1.98721*Tk);
	// bgam = b_gamma;
	bgam = nbg/(2.*log(10.)*(1.98721)*Tk)*2./ni;
	dbgdT = - nbh/(2.*log(10.)*(1.98721)*pow(Tk,2.))*2./ni;
	d2bgdT2 = - nbj/( 2.*log(10.)*(1.98721)*pow(Tk,2.))*2./ni - 2./Tk*dbgdT;
	dbgdP = nbv/(2.*log(10.)*(1.98721)*Tk)*2./ni;

	return 0;
}


/// calculates TP dependence of a_not (and derivatives)
long int THelgeson::IonsizeTP()
{
	double nc, na, ni, zc, za, c;

	switch ( flagElect )
	{
		case 1:  // NaCl
			nc = 1.; na = 1.; ni = 2.;
			zc = 1.; za = -1.;
			break;
		case 2:  // KCl
			nc = 1.; na = 1.; ni = 2.;
			zc = 1.; za = -1.;
			break;
		case 3:  // NaOH
			nc = 1.; na = 1.; ni = 2.;
			zc = 1.; za = -1.;
			break;
		case 4:  // KOH
			nc = 1.; na = 1.; ni = 2.;
			zc = 1.; za = -1.;
			break;
		default:  // wrong mode
			return -1;
	}

	c = 2./ni * ( nc*fabs(zc) + na*fabs(za) );
	ao = ac + c*Gf;
	daodT = c*dGfdT;
	d2aodT2 = c*d2GfdT2;
	daodP = c*dGfdP;

	return 0;
}


/// wrapper for g-function
long int THelgeson::Gfunction()
{
	double T, P, D, beta, alpha, daldT, g, dgdP, dgdT, d2gdT2;
	double TMAX = 1000., PMAX = 5000., TOL = 1.0e-4;

	// convert parameters
	T = Tk - 273.15;
	P = Pbar;
	D = RhoW[0];
	alpha = - RhoW[1]/RhoW[0];
	daldT = pow(alpha, 2.) - RhoW[2]/RhoW[0];
	beta = RhoW[3]/RhoW[0];

	// initialize g and derivatives to zero
	g = 0.0;
	dgdP = 0.0;
	dgdT = 0.0;
	d2gdT2 = 0.0;

	if ((T > TMAX+TOL) || (P > PMAX+TOL))
		return -1;
	else
		GShok2( T, P, D, beta, alpha, daldT, g, dgdP, dgdT, d2gdT2 );

	// assignments
	Gf = g;
	dGfdT = dgdT;
	d2GfdT2 = d2gdT2;
	dGfdP = dgdP;

	return 0;
}


/// calculates g-function and derivatives
long int THelgeson::GShok2( double T, double P, double D, double beta,
		double alpha, double daldT, double &g, double &dgdP, double &dgdT, double &d2gdT2 )
{
	double a, b, dgdD, /*dgdD2,*/ dadT, dadTT, dbdT, dbdTT, dDdT, dDdP,
			dDdTT, Db, dDbdT, dDbdTT, ft, dftdT, dftdTT, fp, dfpdP,
			f, dfdP, dfdT, d2fdT2, tempy;
	double C[6]  = {-0.2037662e+01,  0.5747000e-02, -0.6557892e-05,
			0.6107361e+01, -0.1074377e-01,  0.1268348e-04 };
	double cC[3] = { 0.3666666e+02, -0.1504956e-09,  0.5017997e-13 };

	if ( D >= 1.3 )
		return -1;
	// Sveta 19/02/2000 1-D < 0 pri D>1 => pow(-number, b) = -NaN0000 error
	double pw = fabs(1.0 - D); // insert Sveta 19/02/2000

	// calculation part
	a = C[0] + C[1]*T + C[2]*pow(T,2.);
	b = C[3] + C[4]*T + C[5]*pow(T,2.);
	g = a * pow(pw, b);

	dgdD = - a*b* pow(pw, (b - 1.0));
	// dgdD2 = a * b * (b - 1.0) * pow((1.0 - D),(b - 2.0));

	dadT = C[1] + 2.0*C[2]*T;
	dadTT = 2.0*C[2];
	dbdT = C[4] + 2.0*C[5]*T;
	dbdTT = 2.0*C[5];
	dDdT = - D * alpha;
	dDdP = D * beta;
	dDdTT = - D * (daldT - pow(alpha,2.));

	// Db = pow((1.0 - D), b);  Fixed by DAK 01.11.00
	Db = pow(pw , b);
	dDbdT = - b*pow(pw,(b - 1.0))*dDdT + log(pw)*Db*dbdT;

	dDbdTT = - (b*pow(pw,(b-1.0)) *dDdTT + pow(pw,(b - 1.0)) * dDdT * dbdT
				+ b * dDdT * (-(b - 1.0) * pow(pw,(b - 2.0)) * dDdT
				+ log(pw) * pow(pw,(b - 1.0)) * dbdT))
				+ log(pw) * pow(pw,b) * dbdTT
				- pow(pw,b) * dbdT * dDdT / (1.0 - D)
				+ log(pw) * dbdT * dDbdT;

	dgdP = dgdD * dDdP;
	dgdT = a * dDbdT + Db * dadT;
	d2gdT2 = a * dDbdTT + 2.0 * dDbdT * dadT + Db * dadTT;

	if((T < 155.0) || (P > 1000.0) || (T > 355.0))
		return 0;

	tempy = ((T - 155.0) / 300.0);
	ft = pow(tempy,4.8) + cC[0] * pow(tempy,16.);

	dftdT = 4.8 / 300.0 * pow(tempy, 3.8) + 16.0 / 300.0 * cC[0] * pow(tempy, 15.);

	dftdTT = 3.8 * 4.8 / (300.0 * 300.0) * pow(tempy, 2.8)
		+ 15.0 * 16.0 / (300.0*300.0) * cC[0] * pow(tempy,14.);

	fp = cC[1] * pow((1000.0 - P),3.) + cC[2] * pow((1000.0 - P),4.);

	dfpdP = -3.0 * cC[1] * pow((1000.0 - P),2.) - 4.0 * cC[2] * pow((1000.0 - P),3.);

	f = ft * fp;
	dfdP = ft * dfpdP;
	dfdT = fp * dftdT;
	d2fdT2 = fp * dftdTT;

	g -= f;
	dgdP -= dfdP;
	dgdT -= dfdT;
	d2gdT2 -= d2fdT2;

	return 0;
}





//=============================================================================================
// Extended Debye-Hueckel (EDH) model for aqueous solutions, Davies version
// References: Langmuir (1997)
// (c) TW May 2009
//=============================================================================================


// Generic constructor for the TDavies class
TDavies::TDavies( SolutionData *sd, double *arM, double *arZ, double *dW, double *eW ):
                TSolMod( sd )
{
    alloc_internal();
    m = arM;
    z = arZ;
    RhoW = dW;
    EpsW = eW;
    flagNeut = 0;
    flagH2O = (long int)aIPc[3];  // 0: unity, 1: calculated
    flagMol = (long int)aIPc[5];  // 0: no scale correction, 1: scale correction
}


TDavies::~TDavies()
{
	free_internal();
}


void TDavies::alloc_internal()
{
	LnG = new double [NComp];
	dLnGdT = new double [NComp];
	d2LnGdT2 = new double [NComp];
	dLnGdP = new double [NComp];
}


void TDavies::free_internal()
{
	delete[]LnG;
	delete[]dLnGdT;
	delete[]d2LnGdT2;
	delete[]dLnGdP;
}


/// Calculates T,P corrected parameters
long int TDavies::PTparam()
{
	double alp, bet, dal, rho, eps, dedt, d2edt2, dedp;

	// read and convert rho and eps
	rho = RhoW[0];
	alp = - 1./rho*RhoW[1];
	dal = pow(alp,2.) - 1./rho*RhoW[2];
	bet = 1./rho*RhoW[3];
	eps = EpsW[0];
	dedt = 1./eps*EpsW[1];
	d2edt2 = - 1./pow(eps,2.)*pow(EpsW[1],2.) + 1./eps*EpsW[2];
	dedp = 1./eps*EpsW[3];

	// calculate A term of Debye-Huckel equation (and derivatives)
	A = (1.82483e6)*sqrt(rho) / pow(Tk*eps,1.5);
	dAdT = - 3./2.*A*( dedt + 1./Tk + alp/3. );
	d2AdT2 = 1./A*pow(dAdT,2.) - 3./2.*A*( d2edt2 - 1/pow(Tk,2.) + 1/3.*dal );
	dAdP = 1./2.*A*( bet - 3.*dedp);

 	ErrorIf( fabs(A) < 1e-9, "Davies model",
 			"Error: DH parameter A was not calculated - no values of RhoW and EpsW !" );

	return 0;
}


/// Calculates activity coefficients
long int TDavies::MixMod()
{
	long int j, k, w;
	double sqI, Z2, lgGam, lnGam, Nw, Lgam, lnxw, xw, lgaw, lnaw,
				sig, zt, zz, lg_to_ln;
	lg_to_ln = 2.302585093;

	// get index of water (assumes water is last species in phase)
	w = NComp - 1;

	// calculate ionic strength and total molalities
	IonicStrength();

	xw = x[w];
	Nw = 1000./18.01528;
	// Lgam = -log10(1.+0.0180153*molT);
	Lgam = log10(xw);  // Helgeson large gamma simplified
	if( Lgam < -0.7 )
		Lgam = -0.7;  // experimental truncation of Lgam to min ln(0.5)
	lnxw = log(xw);
	sqI = sqrt(IS);

	// loop over species
	for( j=0; j<NComp; j++ )
	{
		lgGam = 0.0;
		lnGam = 0.0;

		// charged species
		if ( z[j] )
		{
			lgGam = 0.0;
			Z2 = z[j]*z[j];
			lgGam = ( -A * Z2 ) * ( sqI/( 1. + sqI ) - 0.3 * IS );
			if ( flagMol == 0 )
				lnGamma[j] = (lgGam + Lgam) * lg_to_ln;
			else
				lnGamma[j] = lgGam * lg_to_ln;
		}

		// neutral species and water solvent
		else
		{
			// neutral species
			if ( j != (NComp-1) )
			{
				lgGam = 0.0;
				if ( flagMol == 0 )
					lnGamma[j] = (lgGam + Lgam) * lg_to_ln;
				else
					lnGamma[j] = lgGam * lg_to_ln;
				continue;
			}

			// water solvent
			else
			{
				lgGam = 0.0;
				lnGam = 0.0;
				if ( flagH2O == 1 )
				{
					// equation from Wolery (1990), corrected for species charges terms
					zt = 0.;
					sig = 3./sqI * ( 1. + sqI - 1./(1.+sqI) - 2.*log(1.+sqI) );

					for (k=0; k<(NComp-1); k++)
					{
						zt += m[k]*pow(z[k],2.);
					}

					zz = zt/molT;
					lgaw = (1./Nw) * ( - molT/log(10.) + 2./3.*A*sqI*sig - 2.*(-0.3/zz)*A*pow(IS,2.) );
					// lgaw = (1./Nw) * ( - molT/log(10.) + 2./3.*A*sqI*sig - 2.*(-0.3/zz)*A*pow(IS,2.) );
					if ( flagMol == 0 )
						lgaw += molT/(Nw)/log(10.) + lnxw/log(10.);
					lnaw = lgaw * lg_to_ln;
					lnGam = lnaw - lnxw;
				}
				else
					lnGam = 0.0;
				lnGamma[j] = lnGam;
			}
		}
	} // j

	return 0;
}


/// calculates excess properties
long int TDavies::ExcessProp( double *Zex )
{
	long int j, k, w;
	double sqI, Z2, Nw, Lgam, lnxw, xw, lnaw, lgGam, lgaw, sig, zt, zz,
				dawdt, d2awdt2, dawdp, lg_to_ln, g, dgt, d2gt, dgp;
	lg_to_ln = 2.302585093;
	g = 0.; dgt = 0.; d2gt = 0.; dgp = 0.;

	// get index of water (assumes water is last species in phase)
	w = NComp - 1;

	// calculate ionic strength and total molalities
	IonicStrength();

	xw = x[w];
	Nw = 1000./18.01528;
	// Lgam = -log10(1.+0.0180153*molT);
	Lgam = log10(xw);  // Helgeson large gamma simplified
	if( Lgam < -0.7 )
		Lgam = -0.7;  // experimental truncation of Lgam to min ln(0.5)
	lnxw = log(xw);
	sqI = sqrt(IS);

	// loop over species
	for( j=0; j<NComp; j++ )
	{
		// charged species
		if ( z[j] )
		{
			lgGam = 0.0;
			Z2 = z[j]*z[j];
			lgGam = - ( A * Z2 ) * ( sqI/( 1. + sqI ) - 0.3 * IS );
			if ( flagMol == 0 )
				LnG[j] = lgGam * lg_to_ln;
			else
				LnG[j] = (lgGam - Lgam) * lg_to_ln;
			dLnGdT[j] = - ( dAdT * Z2 ) * ( sqI/( 1. + sqI ) - 0.3 * IS ) * lg_to_ln;
			d2LnGdT2[j] = - ( d2AdT2 * Z2 ) * ( sqI/( 1. + sqI ) - 0.3 * IS ) * lg_to_ln;
			dLnGdP[j] = - ( dAdP * Z2 ) * ( sqI/( 1. + sqI ) - 0.3 * IS ) * lg_to_ln;
		}

		// neutral species and water solvent
		else
		{
			// neutral species
			if ( j != (NComp-1) )
			{
				lgGam = 0.0;
				if ( flagMol == 0 )
					LnG[j] = 0.0;
				else
					LnG[j] = (lgGam - Lgam) * lg_to_ln;
				dLnGdT[j] = 0.;
				d2LnGdT2[j] = 0.;
				dLnGdP[j] = 0.;
			}

			// water solvent
			else
			{
				// H2O activity coefficient calculated
				if ( flagH2O == 1 )
				{
					// equation from Wolery (1990), corrected for species charges terms
					zt = 0.;
					sig = 3./sqI * ( 1. + sqI - 1./(1.+sqI) - 2.*log(1.+sqI) );

					for (k=0; k<(NComp-1); k++)
					{
						zt += m[k]*pow(z[k],2.);
					}

					zz = zt/molT;
					// lgaw = (1./Nw) * ( - molT/log(10.) + 2./3.*A*sqI*sig - 2.*(-0.3)*A*pow(IS,2.) );
					lgaw = (1./Nw) * ( - molT/log(10.) + 2./3.*A*sqI*sig - 2.*(-0.3/zz)*A*pow(IS,2.) );
					if ( flagMol == 0 )
						lgaw += molT/(Nw)/log(10.) + lnxw/log(10.);
					lnaw = lgaw * lg_to_ln;
					dawdt = (1./Nw) * ( 2./3.*dAdT*sqI*sig - 2.*(-0.3/zz)*dAdT*pow(IS,2.) ) * lg_to_ln;
					d2awdt2 = (1./Nw) * ( 2./3.*d2AdT2*sqI*sig - 2.*(-0.3/zz)*d2AdT2*pow(IS,2.) ) * lg_to_ln;
					dawdp = (1./Nw) * ( 2./3.*dAdP*sqI*sig - 2.*(-0.3/zz)*dAdP*pow(IS,2.) ) * lg_to_ln;
					LnG[j] = lnaw - lnxw;
					dLnGdT[j] = dawdt;
					d2LnGdT2[j] = d2awdt2;
					dLnGdP[j] = dawdp;
				}

				// H2O activity coefficient 0
				else
				{
					LnG[j] = 0.;
					dLnGdT[j] = 0.;
					d2LnGdT2[j] = 0.;
					dLnGdP[j] = 0.;
				}
			}
		}

		g += x[j]*LnG[j];
		dgt += x[j]*dLnGdT[j];
		d2gt += x[j]*d2LnGdT2[j];
		dgp += x[j]*dLnGdP[j];

	} // j

	// increment thermodynamic properties
	Gex = (R_CONST*Tk) * g;
	Hex = - R_CONST*pow(Tk,2.) * dgt;
	// Sex = - R_CONST * ( g + Tk*dgt );
	Sex = (Hex-Gex)/Tk;
	CPex = - R_CONST * ( 2.*Tk*dgt + pow(Tk,2.)*d2gt );
	Vex = (R_CONST*Tk) * dgp;
	Aex = Gex - Vex*Pbar;
	Uex = Hex - Vex*Pbar;

	// assigments (excess properties)
	Zex[0] = Gex;
	Zex[1] = Hex;
	Zex[2] = Sex;
	Zex[3] = CPex;
	Zex[4] = Vex;
	Zex[5] = Aex;
	Zex[6] = Uex;

	return 0;
}


/// calculates ideal mixing properties
long int TDavies::IdealProp( double *Zid )
{
	long int j;
	double si;
	si = 0.0;
	for (j=0; j<NComp; j++)
	{
		if ( x[j] > 1.0e-32 )
			si += x[j]*log(x[j]);
	}
	Hid = 0.0;
	CPid = 0.0;
	Vid = 0.0;
	Sid = (-1.)*R_CONST*si;
	Gid = Hid - Sid*Tk;
	Aid = Gid - Vid*Pbar;
	Uid = Hid - Vid*Pbar;

	// assignments (ideal mixing properties)
	Zid[0] = Gid;
	Zid[1] = Hid;
	Zid[2] = Sid;
	Zid[3] = CPid;
	Zid[4] = Vid;
	Zid[5] = Aid;
	Zid[6] = Uid;
	return 0;
}


/// calculates true ionic strength
long int TDavies::IonicStrength()
{
	long int j;
	double is, mt;
	is = 0.0; mt = 0.0;

	// calculate ionic strength and total molalities
	// needs to be modified if nonaqueous solvents are present
	for (j=0; j<(NComp-1); j++)
	{
		is += 0.5*m[j]*z[j]*z[j];
		mt += m[j];
	}

	// assignments
	IS = is;
	molT = mt;

	return 0;
}





//=============================================================================================
// Debye-Hueckel (DH) limiting law for aqueous solutions
// References: Langmuir (1997)
// ((c) TW May 2009
//=============================================================================================


// Generic constructor for the TDLimitingLaw class
TLimitingLaw::TLimitingLaw( SolutionData *sd, double *arM, double *arZ, double *dW, double *eW ):
                TSolMod( sd )
{
    alloc_internal();
    m = arM;
    z = arZ;
    RhoW = dW;
    EpsW = eW;
    flagNeut = (long int)aIPc[2];  // 0: unity
    flagH2O = (long int)aIPc[3];  // 0: unity, 1: calculated
}


TLimitingLaw::~TLimitingLaw()
{
	free_internal();
}


void TLimitingLaw::alloc_internal()
{
	LnG = new double [NComp];
	dLnGdT = new double [NComp];
	d2LnGdT2 = new double [NComp];
	dLnGdP = new double [NComp];
}


void TLimitingLaw::free_internal()
{
	delete[]LnG;
	delete[]dLnGdT;
	delete[]d2LnGdT2;
	delete[]dLnGdP;
}


/// Calculates T,P corrected parameters
long int TLimitingLaw::PTparam()
{
	double alp, bet, dal, rho, eps, dedt, d2edt2, dedp;

	// read and convert rho and eps
	rho = RhoW[0];
	alp = - 1./rho*RhoW[1];
	dal = pow(alp,2.) - 1./rho*RhoW[2];
	bet = 1./rho*RhoW[3];
	eps = EpsW[0];
	dedt = 1./eps*EpsW[1];
	d2edt2 = - 1./pow(eps,2.)*pow(EpsW[1],2.) + 1./eps*EpsW[2];
	dedp = 1./eps*EpsW[3];

	// calculate A term of Debye-Huckel equation (and derivatives)
	A = (1.82483e6)*sqrt(rho) / pow(Tk*eps,1.5);
	dAdT = - 3./2.*A*( dedt + 1./Tk + alp/3. );
	d2AdT2 = 1./A*pow(dAdT,2.) - 3./2.*A*( d2edt2 - 1/pow(Tk,2.) + 1/3.*dal );
	dAdP = 1./2.*A*( bet - 3.*dedp);

 	ErrorIf( fabs(A) < 1e-9, "DH limiting law model",
 			"Error: DH parameter A was not calculated - no values of RhoW and EpsW !" );

	return 0;
}


/// Calculates activity coefficients
long int TLimitingLaw::MixMod()
{
	long int j, k, w;
	double sqI, Z2, lgGam, lnGam, Nw, Lgam, lnxw, xw, lnaw, Phi, Phit, lg_to_ln;
	lg_to_ln = 2.302585093;

	// get index of water (assumes water is last species in phase)
	w = NComp - 1;

	// calculate ionic strength and total molalities
	IonicStrength();

	xw = x[w];
	Nw = 1000./18.01528;
	// Lgam = -log10(1.+0.0180153*molT);
	Lgam = log10(xw);  // Helgeson large gamma simplified
	if( Lgam < -0.7 )
		Lgam = -0.7;  // experimental truncation of Lgam to min ln(0.5)
	lnxw = log(xw);
	sqI = sqrt(IS);

	// loop over species
	for( j=0; j<NComp; j++ )
	{
		lgGam = 0.0;
		lnGam = 0.0;

		// charged species
		if ( z[j] )
		{
			lgGam = 0.0;
			Z2 = z[j]*z[j];
			lgGam = ( -A * Z2 * sqI );
			lnGamma[j] = (lgGam + Lgam) * lg_to_ln;
		}

		// neutral species and water solvent
		else
		{
			// neutral species
			if ( j != (NComp-1) )
			{
				lgGam = 0.0;
				lnGamma[j] = (lgGam + Lgam) * lg_to_ln;
				continue;
			}

			// water solvent
			else
			{
				lgGam = 0.0;
				lnGam = 0.0;
				if ( flagH2O == 1 )
				{
					Phit = 0.0;
					// Phi corrected using eq. (190) from Helgeson et al. (1981)

					for (k=0; k<(NComp-1); k++)
					{
							Phit += - log(10.) * m[k] * ( (pow(z[k],2.)*A*sqI)/3. + Lgam/(0.0180153*molT) );
					}

					Phi = Phit/molT;
					lnaw = - Phi*molT/Nw;
					lnGam = lnaw - lnxw;
				}
				else
					lnGam = 0.0;
				lnGamma[j] = lnGam;
			}
		}
	} // j

	return 0;
}


/// calculates excess properties
long int TLimitingLaw::ExcessProp( double *Zex )
{
	long int j, k, w;
	double sqI, Z2, Nw, Lgam, lnxw, xw, Phi, dPhidT, d2PhidT2, dPhidP,
			Phit, dPhitdT, d2PhitdT2, dPhitdP, lnaw, lg_to_ln, g, dgt, d2gt, dgp;
	lg_to_ln = 2.302585093;
	g = 0.; dgt = 0.; d2gt = 0.; dgp = 0.;

	// get index of water (assumes water is last species in phase)
	w = NComp - 1;

	// calculate ionic strength and total molalities
	IonicStrength();

	xw = x[w];
	Nw = 1000./18.01528;
	// Lgam = -log10(1.+0.0180153*molT);
	Lgam = log10(xw);  // Helgeson large gamma simplified
	if( Lgam < -0.7 )
		Lgam = -0.7;  // experimental truncation of Lgam to min ln(0.5)
	lnxw = log(xw);
	sqI = sqrt(IS);

	// loop over species
	for( j=0; j<NComp; j++ )
	{
		// charged species
		if ( z[j] )
		{
			Z2 = z[j]*z[j];
			LnG[j] = - ( A * Z2 * sqI ) * lg_to_ln;
			dLnGdT[j] = - ( dAdT * Z2 * sqI ) * lg_to_ln;
			d2LnGdT2[j] = - ( d2AdT2 * Z2 * sqI ) * lg_to_ln;
			dLnGdP[j] = - ( dAdP * Z2 * sqI ) * lg_to_ln;
		}

		// neutral species and water solvent
		else
		{
			// neutral species
			if ( j != (NComp-1) )
			{
				LnG[j] = 0.;
				dLnGdT[j] = 0.;
				d2LnGdT2[j] = 0.;
				dLnGdP[j] = 0.;
			}

			// water solvent
			else
			{
				if ( flagH2O == 1 )
				{
					Phit = 0.; dPhitdT = 0.; d2PhitdT2 = 0.; dPhitdP = 0.;
					// Phi corrected using eq. (190) from Helgeson et al. (1981)

					for (k=0; k<(NComp-1); k++)
					{
						Phit += - log(10.) * m[k] * ( (pow(z[k],2.)*A*sqI)/3. + Lgam/(0.0180153*molT) );
						dPhitdT += - log(10.) * m[k] * ( (pow(z[k],2.)*dAdT*sqI)/3. );
						d2PhitdT2 += - log(10.) * m[k] * ( (pow(z[k],2.)*d2AdT2*sqI)/3. );
						dPhitdP += - log(10.) * m[k] * ( (pow(z[k],2.)*dAdP*sqI)/3. );
					}

					Phi = Phit/molT;
					dPhidT = dPhitdT/molT;
					d2PhidT2 = d2PhitdT2/molT;
					dPhidP = dPhitdP/molT;

					// activity coefficient (and derivatives)
					lnaw = - Phi*molT/Nw;
					LnG[j] = lnaw - lnxw;
					dLnGdT[j] = - (molT/Nw) * dPhidT;
					d2LnGdT2[j] = - (molT/Nw) * d2PhidT2;
					dLnGdP[j] = - (molT/Nw) * dPhidP;
				}

				else
				{
					LnG[j] = 0.;
					dLnGdT[j] = 0.;
					d2LnGdT2[j] = 0.;
					dLnGdP[j] = 0.;
				}
			}
		}

		g += x[j]*LnG[j];
		dgt += x[j]*dLnGdT[j];
		d2gt += x[j]*d2LnGdT2[j];
		dgp += x[j]*dLnGdP[j];

	} // j

	// increment thermodynamic properties
	Gex = (R_CONST*Tk) * g;
	Hex = - R_CONST*pow(Tk,2.) * dgt;
	// Sex = - R_CONST * ( g + Tk*dgt );
	Sex = (Hex-Gex)/Tk;
	CPex = - R_CONST * ( 2.*Tk*dgt + pow(Tk,2.)*d2gt );
	Vex = (R_CONST*Tk) * dgp;
	Aex = Gex - Vex*Pbar;
	Uex = Hex - Vex*Pbar;

	// assigments (excess properties)
	Zex[0] = Gex;
	Zex[1] = Hex;
	Zex[2] = Sex;
	Zex[3] = CPex;
	Zex[4] = Vex;
	Zex[5] = Aex;
	Zex[6] = Uex;

	return 0;
}


/// calculates ideal mixing properties
long int TLimitingLaw::IdealProp( double *Zid )
{
	long int j;
	double si;
	si = 0.0;
	for (j=0; j<NComp; j++)
	{
		if ( x[j] > 1.0e-32 )
			si += x[j]*log(x[j]);
	}
	Hid = 0.0;
	CPid = 0.0;
	Vid = 0.0;
	Sid = (-1.)*R_CONST*si;
	Gid = Hid - Sid*Tk;
	Aid = Gid - Vid*Pbar;
	Uid = Hid - Vid*Pbar;

	// assignments (ideal mixing properties)
	Zid[0] = Gid;
	Zid[1] = Hid;
	Zid[2] = Sid;
	Zid[3] = CPid;
	Zid[4] = Vid;
	Zid[5] = Aid;
	Zid[6] = Uid;
	return 0;
}


/// calculates true ionic strength
long int TLimitingLaw::IonicStrength()
{
	long int j;
	double is, mt;
	is = 0.0; mt = 0.0;

	// calculate ionic strength and total molalities
	// needs to be modified when nonaqueous solvents are present
	for (j=0; j<(NComp-1); j++)
	{
		is += 0.5*m[j]*z[j]*z[j];
		mt += m[j];
	}

	// assignments
	IS = is;
	molT = mt;

	return 0;
}





//=============================================================================================
// Debye-Hueckel (DH) two term model for aqueous solutions
// References: Helgeson et al. (1981)
// uses individual ion-size parameters, optionally individual salting-out coefficients
// (c) TW May 2009
//=============================================================================================


// Generic constructor for the TDebyeHueckel class
TDebyeHueckel::TDebyeHueckel( SolutionData *sd, double *arM, double *arZ, double *dW, double *eW ):
                TSolMod( sd )
{
    alloc_internal();
    m = arM;
    z = arZ;
    RhoW = dW;
    EpsW = eW;
    ac = 3.72;   // common ion size parameter
    bc = 0.064;   // common b_setch
    flagNeut = (long int)aIPc[2];   // 0: unity, 1: calculated from bg
    flagH2O = (long int)aIPc[3];   // 0: unity, 1: calculated from bg
}


TDebyeHueckel::~TDebyeHueckel()
{
	free_internal();
}


void TDebyeHueckel::alloc_internal()
{
	LnG = new double [NComp];
	dLnGdT = new double [NComp];
	d2LnGdT2 = new double [NComp];
	dLnGdP = new double [NComp];
	an = new double [NComp];
	bg = new double [NComp];
}


void TDebyeHueckel::free_internal()
{
	delete[]LnG;
	delete[]dLnGdT;
	delete[]d2LnGdT2;
	delete[]dLnGdP;
	delete[]an;
	delete[]bg;
}


/// Calculates T,P corrected parameters
long int TDebyeHueckel::PTparam()
{
	long int j;
	double alp, bet, dal, rho, eps, dedt, d2edt2, dedp;

	// read and copy individual an and bg
	for (j=0; j<NComp; j++)
	{
		an[j] = aDCc[NP_DC*j];   // individual an
		bg[j] = aDCc[NP_DC*j+1];   // individual bg
	}

	// read and convert rho and eps
	rho = RhoW[0];
	alp = - 1./rho*RhoW[1];
	dal = pow(alp,2.) - 1./rho*RhoW[2];
	bet = 1./rho*RhoW[3];
	eps = EpsW[0];
	dedt = 1./eps*EpsW[1];
	d2edt2 = - 1./pow(eps,2.)*pow(EpsW[1],2.) + 1./eps*EpsW[2];
	dedp = 1./eps*EpsW[3];

	// calculate A and B terms of Debye-Huckel equation (and derivatives)
	A = (1.82483e6)*sqrt(rho) / pow(Tk*eps,1.5);
	dAdT = - 3./2.*A*( dedt + 1./Tk + alp/3. );
	d2AdT2 = 1./A*pow(dAdT,2.) - 3./2.*A*( d2edt2 - 1/pow(Tk,2.) + 1/3.*dal );
	dAdP = 1./2.*A*( bet - 3.*dedp);
	B = (50.2916)*sqrt(rho) / sqrt(Tk*eps);
	dBdT = - 1./2.*B*( dedt +1./Tk + alp );
	d2BdT2 = 1./B*pow(dBdT,2.) - 1./2.*B*( d2edt2 - 1./pow(Tk,2.) + dal );
	dBdP = 1./2.*B*( bet - dedp );

 	ErrorIf( fabs(A) < 1e-9 || fabs(B) < 1e-9, "DH two-term model",
 			"Error: DH parameter A or B was not calculated - no values of RhoW and EpsW !" );

	return 0;
}


/// Calculates activity coefficients
long int TDebyeHueckel::MixMod()
{
	long int j, k, w;
	double sqI, Z2, lgGam, lnGam, Nw, Lgam, Lam, lnaw, lnxw, xw, lg_to_ln,
				sig, Phi, Phit;
	lg_to_ln = 2.302585093;

	// get index of water (assumes water is last species in phase)
	w = NComp - 1;

	// calculate ionic strength and total molalities
	IonicStrength();

	xw = x[w];
	Nw = 1000./18.01528;
	// Lgam = -log10(1.+0.0180153*molT);
	Lgam = log10(xw);  // Helgeson large gamma simplified
	if( Lgam < -0.7 )
		Lgam = -0.7;  // experimental truncation of Lgam to min ln(0.5)
	lnxw = log(xw);
	sqI = sqrt(IS);

	// loop over species
	for( j=0; j<NComp; j++ )
	{
		lgGam = 0.0;
		lnGam = 0.0;

		// charged species
		if ( z[j] )
		{
			lgGam = 0.0;
			Z2 = z[j]*z[j];
			lgGam = ( - A * sqI * Z2 ) / ( 1. + B * an[j] * sqI ) ;
			lnGamma[j] = (lgGam + Lgam) * lg_to_ln;
		}

		// neutral species and water solvent
		else
		{
			// neutral species
			if ( j != (NComp-1) )
			{
				lgGam = 0.0;

				// rational Setchenow coefficient
				if ( flagNeut == 1 )
					lgGam = bg[j] * IS;

				else
					lgGam = 0.0;

				lnGamma[j] = (lgGam + Lgam) * lg_to_ln;
				continue;
			}

			// water solvent
			else
			{
				lgGam = 0.0;
				lnGam = 0.0;

				// water activity coefficient calculated
				if ( flagH2O == 1 )
				{
					Phit = 0.0;
					// Phi corrected using eq. (190) from Helgeson et al. (1981)
					Lam = 1. + (ao*B) * sqI;
					sig = 3./(pow(ao,3.)*pow(B,3.)*pow(IS,(3./2.))) * (Lam-1./Lam-2*log(Lam));
					// Phi = -2.3025851*(A*sqI*SigTerm/3. + Lgam/(0.0180153*2.*IS) - bgam*IS/2.);
					// Phi = - log(10.) * (molZ/molT) * ( (zc*za*A*sqI*SigTerm)/3. + (Lgam*psi)/(0.0180153*2.*IS) - bgam*IS/2. );

					for (k=0; k<(NComp-1); k++)
					{
						if ( (z[k] == 0) && (flagNeut == 1) )
							Phit += - log(10.) * m[k] * ( (pow(z[k],2.)*A*sqI*sig)/3. + Lgam/(0.0180153*molT) - bg[k]*IS/2. );
						else
							Phit += - log(10.) * m[k] * ( (pow(z[k],2.)*A*sqI*sig)/3. + Lgam/(0.0180153*molT) );
					}

					Phi = Phit/molT;
					lnaw = - Phi*molT/Nw;
					lnGam = lnaw - lnxw;
				}
				else
					lnGam = 0.0;
				lnGamma[j] = lnGam;
			}
		}
	} // j

	return 0;
}


/// calculates excess properties
long int TDebyeHueckel::ExcessProp( double *Zex )
{
	// (under construction)
	long int j, k, w;
	double sqI, Z2, Nw, Lgam, lnxw, xw, Lam, sig, Phi, dPhidT, d2PhidT2, dPhidP,
			Phit, dPhitdT, d2PhitdT2, dPhitdP, lnaw, lg_to_ln, g, dgt, d2gt, dgp;
	double U, V, dUdT, dVdT, d2UdT2, d2VdT2, dUdP, dVdP, U1, U2, U3, V1, V2, V3,
			dU1dT, dU2dT, dU3dT, dV1dT, dV2dT, dV3dT, d2U1dT2, d2U2dT2, d2U3dT2,
			d2V1dT2, d2V2dT2, d2V3dT2, dU1dP, dU2dP, dU3dP, dV1dP, dV2dP, dV3dP,
			L, dLdT, d2LdT2, dLdP, Z, dZdT, d2ZdT2, dZdP;
	lg_to_ln = 2.302585093;
	g = 0.; dgt = 0.; d2gt = 0.; dgp = 0.;

	// get index of water (assumes water is last species in phase)
	w = NComp - 1;

	// calculate ionic strength and total molalities
	IonicStrength();

	xw = x[w];
	Nw = 1000./18.01528;
	// Lgam = -log10(1.+0.0180153*molT);
	Lgam = log10(xw);  // Helgeson large gamma simplified
	if( Lgam < -0.7 )
		Lgam = -0.7;  // experimental truncation of Lgam to min ln(0.5)
	lnxw = log(xw);
	sqI = sqrt(IS);

	// loop over species
	for( j=0; j<NComp; j++ )
	{
		// charged species
		if ( z[j] )
		{
			Z2 = z[j]*z[j];
			U = - (Z2*A) * sqI;
			dUdT = - (Z2*dAdT) * sqI;
			d2UdT2 = - (Z2*d2AdT2) * sqI;
			dUdP = - (Z2*dAdP) * sqI;
			V = 1. + (an[j]*B) * sqI;
			dVdT = ( an[j]*dBdT ) * sqI;
			d2VdT2 = ( an[j]*d2BdT2 ) * sqI;
			dVdP = ( an[j]*dBdP ) * sqI;
			LnG[j] = ( ( - A * sqI * Z2 ) / ( 1. + B * an[j] * sqI ) ) * lg_to_ln;
			dLnGdT[j] = ( (dUdT*V - U*dVdT)/pow(V,2.) ) * lg_to_ln;
			d2LnGdT2[j] = ( (d2UdT2*V + dUdT*dVdT)/pow(V,2.) - (dUdT*V)*(2.*dVdT)/pow(V,3.)
				- (dUdT*dVdT + U*d2VdT2)/pow(V,2.) + (U*dVdT)*(2.*dVdT)/pow(V,3.) ) * lg_to_ln;
			dLnGdP[j] = ( (dUdP*V - U*dVdP)/pow(V,2.) ) * lg_to_ln;
		}

		// neutral species and water solvent
		else
		{
			// neutral species
			if ( j != (NComp-1) )
			{
				// rational Setchenow coefficient
				if ( flagNeut == 1 )
				{
					LnG[j] = ( bg[j] * IS ) * lg_to_ln;
					dLnGdT[j] = 0.;
					d2LnGdT2[j] = 0.;
					dLnGdP[j] = 0.;
				}

				else
				{
					LnG[j] = 0.;
					dLnGdT[j] = 0.;
					d2LnGdT2[j] = 0.;
					dLnGdP[j] = 0.;
				}
				continue;
			}

			// water solvent
			else
			{
				// H2O activity coefficient calculated
				if ( flagH2O == 1 )
				{
					Phit = 0.; dPhitdT = 0.; d2PhitdT2 = 0.; dPhitdP = 0.;

					// derivatives of lambda and sigma terms
					L = 1. + (ao*B) * sqI;
					dLdT = ( ao*dBdT ) * sqI;
					d2LdT2 = ( ao*d2BdT2 ) * sqI;
					dLdP = ( ao*dBdP ) * sqI;

					U1 = (A*L);
					dU1dT = (dAdT*L + A*dLdT);
					d2U1dT2 = ( d2AdT2*L + 2.*dAdT*dLdT + A*d2LdT2 );
					dU1dP = ( dAdP*L + A*dLdP );
					V1 = pow(ao,3.)*pow(B,3.) * IS;
					dV1dT = ( 3.*pow(ao,3.)*pow(B,2.)*dBdT ) * IS;
					d2V1dT2 = ( 6.*pow(ao,3.)*B*pow(dBdT,2.) + 3.*pow(ao,3.)*pow(B,2.)*d2BdT2 ) * IS;
					dV1dP = ( 3.*pow(ao,3.)*pow(B,2.)*dBdP ) * IS;

					U2 = A;
					dU2dT = dAdT;
					d2U2dT2 = d2AdT2;
					dU2dP = dAdP;
					V2 = pow(ao,3.)*pow(B,3.)*L * IS;
					dV2dT = ( 3.*pow(ao,3.)*pow(B,2.)*dBdT*L + pow(ao,3.)*pow(B,3.)*dLdT ) * IS;
					d2V2dT2 = ( 6.*pow(ao,3.)*B*pow(dBdT,2.)*L + 3.*pow(ao,3.)*pow(B,2.)*d2BdT2*L
								+ 6.*pow(ao,3.)*pow(B,2.)*dBdT*dLdT + pow(ao,3.)*pow(B,3.)*d2LdT2 ) * IS;
					dV2dP = ( 3.*pow(ao,3.)*pow(B,2.)*dBdP*L + pow(ao,3.)*pow(B,3.)*dLdP ) * IS;

					U3 = 2.*( A*log(L) );
					dU3dT = 2.*( dAdT*log(L) + A*(1./L)*dLdT );
					d2U3dT2 = 2.*( d2AdT2*log(L) + 2.*dAdT*(1./L)*dLdT
								- A*(1./pow(L,2.))*pow(dLdT,2.) + A*(1./L)*d2LdT2 );
					dU3dP = 2.*( dAdP*log(L) + A*(1./L)*dLdP );
					V3 = pow(ao,3.)*pow(B,3.) * IS;
					dV3dT = ( 3.*pow(ao,3.)*pow(B,2.)*dBdT ) * IS;
					d2V3dT2 = ( 6.*pow(ao,3.)*B*pow(dBdT,2.) + 3.*pow(ao,3.)*pow(B,2.)*d2BdT2 ) * IS;
					dV3dP = ( 3.*pow(ao,3.)*pow(B,2.)*dBdP ) * IS;

					Z = U1/V1 - U2/V2 - U3/V3;
					dZdT = (dU1dT*V1 - U1*dV1dT)/pow(V1,2.) - (dU2dT*V2 - U2*dV2dT)/pow(V2,2.)
								- (dU3dT*V3 - U3*dV3dT)/pow(V3,2.);
					d2ZdT2 = (d2U1dT2*V1 + dU1dT*dV1dT)/pow(V1,2.) - (dU1dT*V1)*(2.*dV1dT)/pow(V1,3.)
								- (dU1dT*dV1dT + U1*d2V1dT2)/pow(V1,2.) + (U1*dV1dT)*(2.*dV1dT)/pow(V1,3.)
								- (d2U2dT2*V2 + dU2dT*dV2dT)/pow(V2,2.) + (dU2dT*V2)*(2.*dV2dT)/pow(V2,3.)
								+ (dU2dT*dV2dT + U2*d2V2dT2)/pow(V2,2.) - (U2*dV2dT)*(2.*dV2dT)/pow(V2,3.)
								- (d2U3dT2*V3 + dU3dT*dV3dT)/pow(V3,2.) + (dU3dT*V3)*(2.*dV3dT)/pow(V3,3.)
								+ (dU3dT*dV3dT + U3*d2V3dT2)/pow(V3,2.) - (U3*dV3dT)*(2.*dV3dT)/pow(V3,3.);
					dZdP = (dU1dP*V1 - U1*dV1dP)/pow(V1,2.) - (dU2dP*V2 - U2*dV2dP)/pow(V2,2.)
								- (dU3dP*V3 - U3*dV3dP)/pow(V3,2.);

					// increments to osmotic coefficient (and derivatives)
					Lam = 1. + (ao*B) * sqI;
					sig = 3./(pow(ao,3.)*pow(B,3.)*pow(IS,(3./2.))) * (Lam-1./Lam-2*log(Lam));

					for (k=0; k<(NComp-1); k++)
					{
						if ( (z[k] == 0) && (flagNeut == 1) )
						{
							Phit += - log(10.) * m[k] * ( (pow(z[k],2.)*A*sqI*sig)/3. + Lgam/(0.0180153*molT) - bg[k]*IS/2. );
							dPhitdT  += - log(10.) * m[k] * ( pow(z[k],2.)*dZdT );
							d2PhitdT2 += - log(10.) * m[k] * ( pow(z[k],2.)*d2ZdT2 );
							dPhitdP += - log(10.) * m[k] * ( pow(z[k],2.)*dZdP );

						}

						else
						{
							Phit += - log(10.) * m[k] * ( (pow(z[k],2.)*A*sqI*sig)/3. + Lgam/(0.0180153*molT) );
							dPhitdT  += - log(10.) * m[k] * ( pow(z[k],2.)*dZdT );
							d2PhitdT2 += - log(10.) * m[k] * ( pow(z[k],2.)*d2ZdT2 );
							dPhitdP += - log(10.) * m[k] * ( pow(z[k],2.)*dZdP );
						}
					}

					Phi = Phit/molT;
					dPhidT = dPhitdT/molT;
					d2PhidT2 = d2PhitdT2/molT;
					dPhidP = dPhitdP/molT;

					// activity coefficient (and derivatives)
					lnaw = - Phi*molT/Nw;
					LnG[j] = lnaw - lnxw;
					dLnGdT[j] = - (molT/Nw) * dPhidT;
					d2LnGdT2[j] = - (molT/Nw) * d2PhidT2;
					dLnGdP[j] = - (molT/Nw) * dPhidP;

				}

				// H2O activity coefficient 0
				else
				{
					LnG[j] = 0.;
					dLnGdT[j] = 0.;
					d2LnGdT2[j] = 0.;
					dLnGdP[j] = 0.;
				}
			}
		}

		g += x[j]*LnG[j];
		dgt += x[j]*dLnGdT[j];
		d2gt += x[j]*d2LnGdT2[j];
		dgp += x[j]*dLnGdP[j];

	} // j

	// increment thermodynamic properties
	Gex = (R_CONST*Tk) * g;
	Hex = - R_CONST*pow(Tk,2.) * dgt;
	// Sex = - R_CONST * ( g + Tk*dgt );
	Sex = (Hex-Gex)/Tk;
	CPex = - R_CONST * ( 2.*Tk*dgt + pow(Tk,2.)*d2gt );
	Vex = (R_CONST*Tk) * dgp;
	Aex = Gex - Vex*Pbar;
	Uex = Hex - Vex*Pbar;

	// assigments (excess properties)
	Zex[0] = Gex;
	Zex[1] = Hex;
	Zex[2] = Sex;
	Zex[3] = CPex;
	Zex[4] = Vex;
	Zex[5] = Aex;
	Zex[6] = Uex;

	return 0;
}


/// calculates ideal mixing properties
long int TDebyeHueckel::IdealProp( double *Zid )
{
	long int j;
	double si;
	si = 0.0;
	for (j=0; j<NComp; j++)
	{
		if ( x[j] > 1.0e-32 )
			si += x[j]*log(x[j]);
	}
	Hid = 0.0;
	CPid = 0.0;
	Vid = 0.0;
	Sid = (-1.)*R_CONST*si;
	Gid = Hid - Sid*Tk;
	Aid = Gid - Vid*Pbar;
	Uid = Hid - Vid*Pbar;

	// assignments (ideal mixing properties)
	Zid[0] = Gid;
	Zid[1] = Hid;
	Zid[2] = Sid;
	Zid[3] = CPid;
	Zid[4] = Vid;
	Zid[5] = Aid;
	Zid[6] = Uid;
	return 0;
}


/// calculates true ionic strength
long int TDebyeHueckel::IonicStrength()
{
	long int j;
	double is, mt, mz, as;
	is = 0.0; mt = 0.0; mz = 0.0; as = 0.0;

	// calculate ionic strength and total molalities
	// needs to be modified when nonaqueous solvents are present
	for (j=0; j<(NComp-1); j++)
	{
		is += 0.5*m[j]*z[j]*z[j];
		mt += m[j];
		if ( z[j] )
		{
			mz += m[j];
			as += m[j]*an[j];
		}
	}

	// assignments
	IS = is;
	molT = mt;
	ao = as/mz;

	return 0;
}





//=============================================================================================
// Extended Debye-Hueckel (EDH) model for aqueous solutions, Karpov version
// References: Karpov et al. (1997); Helgeson et al. (1981); Oelkers and Helgeson (1990);
// Pokrovskii and Helgeson (1995; 1997a; 1997b)
// (c) TW May 2009
//=============================================================================================


// Generic constructor for the TKarpov class
TKarpov::TKarpov( SolutionData *sd, double *arM, double *arZ, double *dW, double *eW ):
                TSolMod( sd )
{
    alloc_internal();
    m = arM;
    z = arZ;
    RhoW = dW;
    EpsW = eW;
    ac = aIPc[1];   // common ion size parameter
    bc = aIPc[0];   // common b_gamma
    flagNeut = (long int)aIPc[2];   // 0: unity, 1: calculated
    flagH2O = (long int)aIPc[3];   // 0: unity, 1: calculated
    flagElect = (long int)aIPc[4];  // 0: constant, 1: NaCl, 2: KCl, 3: NaOH, 4: KOH
}


TKarpov::~TKarpov()
{
	free_internal();
}


void TKarpov::alloc_internal()
{
	LnG = new double [NComp];
	dLnGdT = new double [NComp];
	d2LnGdT2 = new double [NComp];
	dLnGdP = new double [NComp];
	an = new double [NComp];
	bg = new double [NComp];
}


void TKarpov::free_internal()
{
	delete[]LnG;
	delete[]dLnGdT;
	delete[]d2LnGdT2;
	delete[]dLnGdP;
	delete[]an;
	delete[]bg;
}


/// Calculates T,P corrected parameters
long int TKarpov::PTparam()
{
	long int j;
	double alp, bet, dal, rho, eps, dedt, d2edt2, dedp;

	// read and copy individual an and bg
	for (j=0; j<NComp; j++)
	{
		an[j] = aDCc[NP_DC*j];
		bg[j] = aDCc[NP_DC*j+1];  // individual bg (not used)
	}

	// read and convert rho and eps
	rho = RhoW[0];
	alp = - 1./rho*RhoW[1];
	dal = pow(alp,2.) - 1./rho*RhoW[2];
	bet = 1./rho*RhoW[3];
	eps = EpsW[0];
	dedt = 1./eps*EpsW[1];
	d2edt2 = - 1./pow(eps,2.)*pow(EpsW[1],2.) + 1./eps*EpsW[2];
	dedp = 1./eps*EpsW[3];

	// calculate A and B terms of Debye-Huckel equation (and derivatives)
	A = (1.82483e6)*sqrt(rho) / pow(Tk*eps,1.5);
	dAdT = - 3./2.*A*( dedt + 1./Tk + alp/3. );
	d2AdT2 = 1./A*pow(dAdT,2.) - 3./2.*A*( d2edt2 - 1/pow(Tk,2.) + 1/3.*dal );
	dAdP = 1./2.*A*( bet - 3.*dedp);
	B = (50.2916)*sqrt(rho) / sqrt(Tk*eps);
	dBdT = - 1./2.*B*( dedt +1./Tk + alp );
	d2BdT2 = 1./B*pow(dBdT,2.) - 1./2.*B*( d2edt2 - 1./pow(Tk,2.) + dal );
	dBdP = 1./2.*B*( bet - dedp );

 	ErrorIf( fabs(A) < 1e-9 || fabs(B) < 1e-9, "Karpov EDH model",
 			"Error: DH parameter A or B was not calculated - no values of RhoW and EpsW !" );

	// b_gamma constant
	if ( flagElect == 0)
	{
		bgam = bc;
		dbgdT = 0.0;
		d2bgdT2 = 0.0;
		dbgdP = 0.0;
		ao = ac;  // constant
	}

	// b_gamma TP-dependent
	else
	{
		Gfunction();
		BgammaTP();
		ao = ac;  // constant
	}

	return 0;
}


/// Calculates activity coefficients
long int TKarpov::MixMod()
{
	long int j, k, w;
	double sqI, Z2, lgGam, lnGam, Nw, Lgam, lnxw, xw, Lam, sig,
			Phi, Phit, psi, zc, za, lnaw, lg_to_ln;
	zc = 1.; za = 1.; psi = 1.; lg_to_ln = 2.302585093;

	// get index of water (assumes water is last species in phase)
	w = NComp - 1;

	// calculate ionic strength and total molalities
	IonicStrength();

	xw = x[w];
	Nw = 1000./18.01528;
	// Lgam = -log10(1.+0.0180153*molT);
	Lgam = log10(xw);  // Helgeson large gamma simplified
	if( Lgam < -0.7 )
		Lgam = -0.7;  // experimental truncation of Lgam to min ln(0.5)
	lnxw = log(xw);
	sqI = sqrt(IS);

	// loop over species
	for( j=0; j<NComp; j++ )
	{
		lgGam = 0.0;
		lnGam = 0.0;

		// charged species (individual ion size parameters)
		if ( z[j] )
		{
			lgGam = 0.0;
			Z2 = z[j]*z[j];
			lgGam = ( - A * sqI * Z2 ) / ( 1. + B * an[j] * sqI ) + bgam * IS ;
			lnGamma[j] = (lgGam + Lgam) * lg_to_ln;
		}

		// neutral species and water solvent
		else
		{
			// neutral species
			if ( j != (NComp-1) )
			{
				lgGam = 0.0;
				if ( flagNeut == 1 )
					lgGam = bgam * IS;
				else
					lgGam = 0.0;
				lnGamma[j] = (lgGam + Lgam) * lg_to_ln;
				continue;
			}

			// water solvent
			else
			{
				lgGam = 0.0;
				lnGam = 0.0;
				if ( flagH2O == 1 )
				{
					Phit = 0.0;
					// Phi corrected using eq. (190) from Helgeson et al. (1981)
					Lam = 1. + ao*B*sqI;
					sig = 3./(pow(ao,3.)*pow(B,3.)*pow(IS,(3./2.))) * (Lam-1./Lam-2*log(Lam));
					// Phi = -2.3025851*(A*sqI*SigTerm/3. + Lgam/(0.0180153*2.*IS) - bgam*IS/2.);
					// Phi = - log(10.) * (molZ/molT) * ( (zc*za*A*sqI*SigTerm)/3. + (psi*Lgam)/(0.0180153*2.*IS) - bgam*IS/2. );

					for (k=0; k<(NComp-1); k++)
					{
						if ( (z[k] == 0) && (flagNeut == 0) )
							Phit += - log(10.) * m[k] * ( (pow(z[k],2.)*A*sqI*sig)/3. + Lgam/(0.0180153*molT) );
						else
							Phit += - log(10.) * m[k] * ( (pow(z[k],2.)*A*sqI*sig)/3. + Lgam/(0.0180153*molT) - bgam*IS/2. );
					}

					Phi = Phit/molT;
					lnaw = - Phi*molT/Nw;
					lnGam = lnaw - lnxw;
				}
				else
					lnGam = 0.0;
				lnGamma[j] = lnGam;
			}
		}
	} // j

	return 0;
}


/// calculates excess properties
long int TKarpov::ExcessProp( double *Zex )
{
	long int j, k, w;
	double sqI, Z2, Nw, Lgam, lnxw, xw, Lam, sig, Phi, dPhidT, d2PhidT2, dPhidP,
			Phit, dPhitdT, d2PhitdT2, dPhitdP, lnaw, lg_to_ln, g, dgt, d2gt, dgp;
	double U, V, dUdT, dVdT, d2UdT2, d2VdT2, dUdP, dVdP, U1, U2, U3, V1, V2, V3,
			dU1dT, dU2dT, dU3dT, dV1dT, dV2dT, dV3dT, d2U1dT2, d2U2dT2, d2U3dT2,
			d2V1dT2, d2V2dT2, d2V3dT2, dU1dP, dU2dP, dU3dP, dV1dP, dV2dP, dV3dP,
			L, dLdT, d2LdT2, dLdP, Z, dZdT, d2ZdT2, dZdP;
	lg_to_ln = 2.302585093;
	g = 0.; dgt = 0.; d2gt = 0.; dgp = 0.;

	// get index of water (assumes water is last species in phase)
	w = NComp - 1;

	// calculate ionic strength and total molalities
	IonicStrength();

	xw = x[w];
	Nw = 1000./18.01528;
	// Lgam = -log10(1.+0.0180153*molT);
	Lgam = log10(xw);  // Helgeson large gamma simplified
	if( Lgam < -0.7 )
		Lgam = -0.7;  // experimental truncation of Lgam to min ln(0.5)
	lnxw = log(xw);
	sqI = sqrt(IS);

	// loop over species
	for( j=0; j<NComp; j++ )
	{
		// charged species
		if ( z[j] )
		{
			Z2 = z[j]*z[j];
			U = - (Z2*A) * sqI;
			dUdT = - (Z2*dAdT) * sqI;
			d2UdT2 = - (Z2*d2AdT2) * sqI;
			dUdP = - (Z2*dAdP) * sqI;
			V = 1. + (an[j]*B) * sqI;
			dVdT = ( an[j]*dBdT ) * sqI;
			d2VdT2 = ( an[j]*d2BdT2 ) * sqI;
			dVdP = ( an[j]*dBdP ) * sqI;
			LnG[j] = ( U/V + bgam*IS ) * lg_to_ln;
			dLnGdT[j] = ( (dUdT*V - U*dVdT)/pow(V,2.) + dbgdT*IS ) * lg_to_ln;
			d2LnGdT2[j] = ( (d2UdT2*V + dUdT*dVdT)/pow(V,2.) - (dUdT*V)*(2.*dVdT)/pow(V,3.)
				- (dUdT*dVdT + U*d2VdT2)/pow(V,2.) + (U*dVdT)*(2.*dVdT)/pow(V,3.)
				+ d2bgdT2*IS ) * lg_to_ln;
			dLnGdP[j] = ( (dUdP*V - U*dVdP)/pow(V,2.) + dbgdP*IS ) * lg_to_ln;
		}

		// neutral species and water solvent
		else
		{
			// neutral species
			if ( j != (NComp-1) )
			{
				if ( flagNeut == 1 )
				{
					LnG[j] = ( bgam*IS ) * lg_to_ln;
					dLnGdT[j] = ( dbgdT*IS ) * lg_to_ln;
					d2LnGdT2[j] = ( d2bgdT2*IS ) * lg_to_ln;
					dLnGdP[j] = ( dbgdP*IS ) * lg_to_ln;
				}

				else
				{
					LnG[j] = 0.;
					dLnGdT[j] = 0.;
					d2LnGdT2[j] = 0.;
					dLnGdP[j] = 0.;
				}
				continue;
			}

			// water solvent
			else
			{
				// H2O activity coefficient calculated
				if ( flagH2O == 1 )
				{
					Phit = 0.; dPhitdT = 0.; d2PhitdT2 = 0.; dPhitdP = 0.;

					// derivatives of lambda and sigma terms
					L = 1. + (ao*B) * sqI;
					dLdT = ( ao*dBdT ) * sqI;
					d2LdT2 = ( ao*d2BdT2 ) * sqI;
					dLdP = ( ao*dBdP ) * sqI;

					U1 = (A*L);
					dU1dT = (dAdT*L + A*dLdT);
					d2U1dT2 = ( d2AdT2*L + 2.*dAdT*dLdT + A*d2LdT2 );
					dU1dP = ( dAdP*L + A*dLdP );
					V1 = pow(ao,3.)*pow(B,3.) * IS;
					dV1dT = ( 3.*pow(ao,3.)*pow(B,2.)*dBdT ) * IS;
					d2V1dT2 = ( 6.*pow(ao,3.)*B*pow(dBdT,2.) + 3.*pow(ao,3.)*pow(B,2.)*d2BdT2 ) * IS;
					dV1dP = ( 3.*pow(ao,3.)*pow(B,2.)*dBdP ) * IS;

					U2 = A;
					dU2dT = dAdT;
					d2U2dT2 = d2AdT2;
					dU2dP = dAdP;
					V2 = pow(ao,3.)*pow(B,3.)*L * IS;
					dV2dT = ( 3.*pow(ao,3.)*pow(B,2.)*dBdT*L + pow(ao,3.)*pow(B,3.)*dLdT ) * IS;
					d2V2dT2 = ( 6.*pow(ao,3.)*B*pow(dBdT,2.)*L + 3.*pow(ao,3.)*pow(B,2.)*d2BdT2*L
								+ 6.*pow(ao,3.)*pow(B,2.)*dBdT*dLdT + pow(ao,3.)*pow(B,3.)*d2LdT2 ) * IS;
					dV2dP = ( 3.*pow(ao,3.)*pow(B,2.)*dBdP*L + pow(ao,3.)*pow(B,3.)*dLdP ) * IS;

					U3 = 2.*( A*log(L) );
					dU3dT = 2.*( dAdT*log(L) + A*(1./L)*dLdT );
					d2U3dT2 = 2.*( d2AdT2*log(L) + 2.*dAdT*(1./L)*dLdT
								- A*(1./pow(L,2.))*pow(dLdT,2.) + A*(1./L)*d2LdT2 );
					dU3dP = 2.*( dAdP*log(L) + A*(1./L)*dLdP );
					V3 = pow(ao,3.)*pow(B,3.) * IS;
					dV3dT = ( 3.*pow(ao,3.)*pow(B,2.)*dBdT ) * IS;
					d2V3dT2 = ( 6.*pow(ao,3.)*B*pow(dBdT,2.) + 3.*pow(ao,3.)*pow(B,2.)*d2BdT2 ) * IS;
					dV3dP = ( 3.*pow(ao,3.)*pow(B,2.)*dBdP ) * IS;

					Z = U1/V1 - U2/V2 - U3/V3;
					dZdT = (dU1dT*V1 - U1*dV1dT)/pow(V1,2.) - (dU2dT*V2 - U2*dV2dT)/pow(V2,2.)
								- (dU3dT*V3 - U3*dV3dT)/pow(V3,2.);
					d2ZdT2 = (d2U1dT2*V1 + dU1dT*dV1dT)/pow(V1,2.) - (dU1dT*V1)*(2.*dV1dT)/pow(V1,3.)
								- (dU1dT*dV1dT + U1*d2V1dT2)/pow(V1,2.) + (U1*dV1dT)*(2.*dV1dT)/pow(V1,3.)
								- (d2U2dT2*V2 + dU2dT*dV2dT)/pow(V2,2.) + (dU2dT*V2)*(2.*dV2dT)/pow(V2,3.)
								+ (dU2dT*dV2dT + U2*d2V2dT2)/pow(V2,2.) - (U2*dV2dT)*(2.*dV2dT)/pow(V2,3.)
								- (d2U3dT2*V3 + dU3dT*dV3dT)/pow(V3,2.) + (dU3dT*V3)*(2.*dV3dT)/pow(V3,3.)
								+ (dU3dT*dV3dT + U3*d2V3dT2)/pow(V3,2.) - (U3*dV3dT)*(2.*dV3dT)/pow(V3,3.);
					dZdP = (dU1dP*V1 - U1*dV1dP)/pow(V1,2.) - (dU2dP*V2 - U2*dV2dP)/pow(V2,2.)
								- (dU3dP*V3 - U3*dV3dP)/pow(V3,2.);

					// increments to osmotic coefficient (and derivatives)
					Lam = 1. + (ao*B) * sqI;
					sig = 3./(pow(ao,3.)*pow(B,3.)*pow(IS,(3./2.))) * (Lam-1./Lam-2*log(Lam));

					for (k=0; k<(NComp-1); k++)
					{
						if ( (z[k] == 0) && (flagNeut == 0) )
						{
							Phit += - log(10.) * m[k] * ( (pow(z[k],2.)*A*sqI*sig)/3. + Lgam/(0.0180153*molT) );
							dPhitdT  += - log(10.) * m[k] * ( pow(z[k],2.)*dZdT );
							d2PhitdT2 += - log(10.) * m[k] * ( pow(z[k],2.)*d2ZdT2 );
							dPhitdP += - log(10.) * m[k] * ( pow(z[k],2.)*dZdP );

						}

						else
						{
							Phit += - log(10.) * m[k] * ( (pow(z[k],2.)*A*sqI*sig)/3. + Lgam/(0.0180153*molT) - bgam*IS/2. );
							dPhitdT  += - log(10.) * m[k] * ( pow(z[k],2.)*dZdT - dbgdT*IS/2. );
							d2PhitdT2 += - log(10.) * m[k] * ( pow(z[k],2.)*d2ZdT2 - d2bgdT2*IS/2. );
							dPhitdP += - log(10.) * m[k] * ( pow(z[k],2.)*dZdP - dbgdP*IS/2. );
						}
					}

					Phi = Phit/molT;
					dPhidT = dPhitdT/molT;
					d2PhidT2 = d2PhitdT2/molT;
					dPhidP = dPhitdP/molT;

					// activity coefficient (and derivatives)
					lnaw = - Phi*molT/Nw;
					LnG[j] = lnaw - lnxw;
					dLnGdT[j] = - (molT/Nw) * dPhidT;
					d2LnGdT2[j] = - (molT/Nw) * d2PhidT2;
					dLnGdP[j] = - (molT/Nw) * dPhidP;
				}

				// H2O activity coefficient 0
				else
				{
					LnG[j] = 0.;
					dLnGdT[j] = 0.;
					d2LnGdT2[j] = 0.;
					dLnGdP[j] = 0.;
				}
			}
		}

		g += x[j]*LnG[j];
		dgt += x[j]*dLnGdT[j];
		d2gt += x[j]*d2LnGdT2[j];
		dgp += x[j]*dLnGdP[j];

	} // j

	// increment thermodynamic properties
	Gex = (R_CONST*Tk) * g;
	Hex = - R_CONST*pow(Tk,2.) * dgt;
	// Sex = - R_CONST * ( g + Tk*dgt );
	Sex = (Hex-Gex)/Tk;
	CPex = - R_CONST * ( 2.*Tk*dgt + pow(Tk,2.)*d2gt );
	Vex = (R_CONST*Tk) * dgp;
	Aex = Gex - Vex*Pbar;
	Uex = Hex - Vex*Pbar;

	// assigments (excess properties)
	Zex[0] = Gex;
	Zex[1] = Hex;
	Zex[2] = Sex;
	Zex[3] = CPex;
	Zex[4] = Vex;
	Zex[5] = Aex;
	Zex[6] = Uex;

	return 0;
}


/// calculates ideal mixing properties
long int TKarpov::IdealProp( double *Zid )
{
	long int j;
	double si;
	si = 0.0;
	for (j=0; j<NComp; j++)
	{
		if ( x[j] > 1.0e-32 )
			si += x[j]*log(x[j]);
	}
	Hid = 0.0;
	CPid = 0.0;
	Vid = 0.0;
	Sid = (-1.)*R_CONST*si;
	Gid = Hid - Sid*Tk;
	Aid = Gid - Vid*Pbar;
	Uid = Hid - Vid*Pbar;

	// assignments (ideal mixing properties)
	Zid[0] = Gid;
	Zid[1] = Hid;
	Zid[2] = Sid;
	Zid[3] = CPid;
	Zid[4] = Vid;
	Zid[5] = Aid;
	Zid[6] = Uid;
	return 0;
}


/// calculates true ionic strength
long int TKarpov::IonicStrength()
{
	long int j;
	double is, mt, mz, as;
	is = 0.0; mt = 0.0; mz = 0.0; as = 0.0;

	// calculate ionic strength and total molalities
	// needs to be modified when nonaqueous solvents are present
	for (j=0; j<(NComp-1); j++)
	{
		is += 0.5*m[j]*z[j]*z[j];
		mt += m[j];
		if ( z[j] )
		{
			mz += m[j];
			as += m[j]*an[j];
		}
	}

	// conversions and assignments
	ao = as/mz;
	IS = is;
	molT = mt;
	molZ = mz;

	return 0;
}


/// calculates TP dependence of b_gamma (and derivatives)
long int TKarpov::BgammaTP()
{
	// ni: stoichiometric number of moles of ions in one mole of electrolyte
	// rc, ra: radius of cation and anion, respectively at 298 K/1 bar
	// units are cal, kg, K, mol, bar
	double ni, nc, na, zc, za, rc, ra, a1, a2, a3, a4, a5, c1, c2, omg, bg, bs, bh, rec, rea,
			omgpt, domdt, d2omdt2, domdp, nbg, nbv, nbj, nbh;
	double eps, eta, xborn, yborn, qborn, X1, X2;

	// set parameters
	eps = EpsW[0];
	yborn = 1./pow(EpsW[0],2.)*EpsW[1];
	xborn = 1./pow(EpsW[0],2.) * ( EpsW[2] - 2./EpsW[0]*pow(EpsW[1],2.) );
	qborn = 1./pow(EpsW[0],2.)*EpsW[3];
	eta = (1.66027e5);

	switch ( flagElect )
	{
		case 1:  // NaCl
			ni = 2.; nc = 1.; na = 1.; zc = 1.; za = -1.;
			a1 = 0.030056; a2 = -202.55; a3 = -2.9092; a4 = 20302;
			a5 = -0.206; c1 = -1.50; c2 = 53300.; omg = 178650.;
			bg = -174.623; bs = 2.164; rc = 0.97; ra = 1.81;
			break;
		case 2:  // KCl
			ni = 2.; nc = 1.; na = 1.; zc = 1.; za = -1.;
			a1 = 0.0172; a2 = -115.36; a3 = -1.1857; a4 = 13854.2;
			a5 = -0.262; c1 = -2.53; c2 = 38628.4; omg = 164870.;
			bg = -70.0; bs = 1.727; rc = 1.33; ra = 1.81;
			break;
		case 3:  // NaOH
			ni = 2.; nc = 1.; na = 1.; zc = 1.; za = -1.;
			a1 = 0.030056; a2 = -202.55; a3 = -2.9092; a4 = 20302;
			a5 = -0.206; c1 = -1.50; c2 = 53300.; omg = 205520.;
			bg = -267.4; bs = 1.836; rc = 0.97; ra = 1.40;
			break;
		case 4:  // KOH
			ni = 2.; nc = 1.; na = 1.; zc = 1.; za = -1.;
			a1 = 0.0172; a2 = -115.36; a3 = -1.1857; a4 = 13854.2;
			a5 = -0.262; c1 = -2.53; c2 = 38628.4; omg = 191730.;
			bg = -335.7; bs = 1.26; rc = 1.33; ra = 1.40;
			break;
		default:  // wrong mode
			return -1;
	}

	// calculation part, extended 06.06.2009 (TW)
	bh = bg + (298.15)*bs;
	rec = rc + fabs(zc)*(0.94+Gf);
	rea = ra + fabs(za)*Gf;
	X1 = - eta*nc*( fabs(pow(zc,3.))/pow(rec,2.) - zc/pow((3.082+Gf),2.) )
			- eta*na*( fabs(pow(za,3.))/pow(rea,2.) - za/pow((3.082+Gf),2.) );
	X2 = 2.*eta*nc*( fabs(pow(zc,4.))/pow(rec,3.) - zc/pow((3.082+Gf),3.) )
			+ 2.*eta*nc * ( fabs(pow(zc,4.))/pow(rec,3.) - zc/pow((3.082+Gf),3.) );
	omgpt = eta*( nc*pow(zc,2.)/rec + na*pow(za,2.)/rea );
	// omgpt = (1.66027e5)*(1./(0.94+rc+Gf)+1./(ra+Gf));
	domdt = X1*dGfdT;
	d2omdt2 = X2*pow(dGfdT,2.) + X1*d2GfdT2;
	domdp = X1*dGfdP;
	nbg = - ni*bg/2. + ni*bs*(Tk-298.15)/2. - c1*(Tk*log(Tk/298.15)-Tk+298.15)
				+ a1*(Pbar-1.) + a2*log((2600.+Pbar)/(2600.+1.))
				- c2*((1./(Tk-228.)-1./(298.15-228.))*(228.-Tk)/228.-Tk/(228.*228.)
				* log((298.15*(Tk-228.))/(Tk*(298.15-228.))))
				+ 1./(Tk-228.)*(a3*(Pbar-1.) + a4*log((2600.+Pbar)/(2600.+1.)))
				+ a5*(omgpt*(1./eps-1.)-omg*(1./(78.24513795)-1.)+(-5.798650444e-5)*omg*(Tk-298.15));
	nbh = - ni/2.*bh + c1*(Tk-298.15) - c2*(1./(Tk-228.)-1./(298.15-228.))
				+ a1*(Pbar-1.) + a2*log((2600.+Pbar)/(2600.+1.))
				+ ((2.*Tk-228.)/pow((Tk-228.),2.))*(a3*(Pbar-1.)+a4*log((2600.+Pbar)/(2600.+1.)))
				+ a5*( omgpt*(1./eps-1.) + omgpt*Tk*yborn - Tk*(1./eps-1.)*domdt
				- omg*(1./(78.24513795)-1.) - omg*(298.15)*(-5.798650444e-5) );
	nbj = c1 + c2/pow((Tk-228.),2.)  - ( (2.*Tk)/pow((Tk-228.),3.) )
				* ( a3*(Pbar-1.) + a4*log((2600.+Pbar)/(2600.+1.)) )
				+ a5*( omgpt*Tk*xborn + 2.*Tk*yborn*domdt - Tk*(1./eps-1.)*d2omdt2 );
	nbv = a1 + a2/(2600.+Pbar) + a3/(Tk-228.) + a4/((2600.+Pbar)*(Tk-228.))
				+ a5*(-omgpt*qborn + (1./eps-1.)*domdp);
	// b_gamma = nbg/(2.*log(10.)*1.98721*Tk);
	// bgam = b_gamma;
	bgam = nbg/(2.*log(10.)*(1.98721)*Tk)*2./ni;
	dbgdT = - nbh/(2.*log(10.)*(1.98721)*pow(Tk,2.))*2./ni;
	d2bgdT2 = - nbj/( 2.*log(10.)*(1.98721)*pow(Tk,2.))*2./ni - 2./Tk*dbgdT;
	dbgdP = nbv/(2.*log(10.)*(1.98721)*Tk)*2./ni;

	return 0;
}


/// wrapper for g-function
long int TKarpov::Gfunction()
{
	double T, P, D, beta, alpha, daldT, g, dgdP, dgdT, d2gdT2;
	double TMAX = 1000., PMAX = 5000., TOL = 1.0e-4;

	// convert parameters
	T = Tk - 273.15;
	P = Pbar;
	D = RhoW[0];
	alpha = - RhoW[1]/RhoW[0];
	daldT = pow(alpha, 2.) - RhoW[2]/RhoW[0];
	beta = RhoW[3]/RhoW[0];

	// initialize g and derivatives to zero
	g = 0.0;
	dgdP = 0.0;
	dgdT = 0.0;
	d2gdT2 = 0.0;

	if ((T > TMAX+TOL) || (P > PMAX+TOL))
		return -1;
	else
		GShok2( T, P, D, beta, alpha, daldT, g, dgdP, dgdT, d2gdT2 );

	// assignments
	Gf = g;
	dGfdT = dgdT;
	d2GfdT2 = d2gdT2;
	dGfdP = dgdP;

	return 0;
}


/// calculates g-function and derivatives
long int TKarpov::GShok2( double T, double P, double D, double beta,
		double alpha, double daldT, double &g, double &dgdP, double &dgdT, double &d2gdT2 )
{
	double a, b, dgdD, /*dgdD2,*/ dadT, dadTT, dbdT, dbdTT, dDdT, dDdP,
		dDdTT, Db, dDbdT, dDbdTT, ft, dftdT, dftdTT, fp, dfpdP,
		f, dfdP, dfdT, d2fdT2, tempy;
	double C[6]  = {-0.2037662e+01,  0.5747000e-02, -0.6557892e-05,
			0.6107361e+01, -0.1074377e-01,  0.1268348e-04 };
	double cC[3] = { 0.3666666e+02, -0.1504956e-09,  0.5017997e-13 };

	if ( D >= 1.3 )
		return -1;
	// Sveta 19/02/2000 1-D < 0 pri D>1 => pow(-number, b) = -NaN0000 error
	double pw = fabs(1.0 - D); // insert Sveta 19/02/2000

	// calculation part
	a = C[0] + C[1]*T + C[2]*pow(T,2.);
	b = C[3] + C[4]*T + C[5]*pow(T,2.);
	g = a * pow(pw, b);

	dgdD = - a*b* pow(pw, (b - 1.0));
	// dgdD2 = a * b * (b - 1.0) * pow((1.0 - D),(b - 2.0));

	dadT = C[1] + 2.0*C[2]*T;
	dadTT = 2.0*C[2];
	dbdT = C[4] + 2.0*C[5]*T;
	dbdTT = 2.0*C[5];
	dDdT = - D * alpha;
	dDdP = D * beta;
	dDdTT = - D * (daldT - pow(alpha,2.));

	// Db = pow((1.0 - D), b);  Fixed by DAK 01.11.00
	Db = pow(pw , b);
	dDbdT = - b*pow(pw,(b - 1.0))*dDdT + log(pw)*Db*dbdT;

	dDbdTT = - (b*pow(pw,(b-1.0)) *dDdTT + pow(pw,(b - 1.0)) * dDdT * dbdT
				+ b * dDdT * (-(b - 1.0) * pow(pw,(b - 2.0)) * dDdT
				+ log(pw) * pow(pw,(b - 1.0)) * dbdT))
				+ log(pw) * pow(pw,b) * dbdTT
				- pow(pw,b) * dbdT * dDdT / (1.0 - D)
				+ log(pw) * dbdT * dDbdT;

	dgdP = dgdD * dDdP;
	dgdT = a * dDbdT + Db * dadT;
	d2gdT2 = a * dDbdTT + 2.0 * dDbdT * dadT + Db * dadTT;

	if((T < 155.0) || (P > 1000.0) || (T > 355.0))
		return 0;

	tempy = ((T - 155.0) / 300.0);
	ft = pow(tempy,4.8) + cC[0] * pow(tempy,16.);

	dftdT = 4.8 / 300.0 * pow(tempy, 3.8) + 16.0 / 300.0 * cC[0] * pow(tempy, 15.);

	dftdTT = 3.8 * 4.8 / (300.0 * 300.0) * pow(tempy, 2.8)
		+ 15.0 * 16.0 / (300.0*300.0) * cC[0] * pow(tempy,14.);

	fp = cC[1] * pow((1000.0 - P),3.) + cC[2] * pow((1000.0 - P),4.);

	dfpdP = -3.0 * cC[1] * pow((1000.0 - P),2.) - 4.0 * cC[2] * pow((1000.0 - P),3.);

	f = ft * fp;
	dfdP = ft * dfpdP;
	dfdT = fp * dftdT;
	d2fdT2 = fp * dftdTT;

	g -= f;
	dgdP -= dfdP;
	dgdT -= dfdT;
	d2gdT2 -= d2fdT2;

	return 0;
}





//=============================================================================================
// Extended Debye-Hueckel (EDH) model for aqueous solutions, Shvarov version
// References: Shvarov (2007); Oelkers and Helgeson (1990);
// Pokrovskii and Helgeson (1995; 1997a; 1997b)
// (c) TW June 2009
//=============================================================================================


// Generic constructor for the TShvarov class
TShvarov::TShvarov( SolutionData *sd, double *arM, double *arZ, double *dW, double *eW ):
                TSolMod( sd )
{
    alloc_internal();
    m = arM;
    z = arZ;
    RhoW = dW;
    EpsW = eW;
    ac = aIPc[1];   // common ion size parameter
    bc = aIPc[0];   // common b_gamma
    flagNeut = (long int)aIPc[2];   // 0: unity, 1: calculated
    flagH2O = (long int)aIPc[3];   // 0: unity, 1: calculated
    flagElect = (long int)aIPc[4];  // 0: constant, 1: NaCl, 2: KCl, 3: NaOH, 4: KOH
}


TShvarov::~TShvarov()
{
	free_internal();
}


void TShvarov::alloc_internal()
{
	LnG = new double [NComp];
	dLnGdT = new double [NComp];
	d2LnGdT2 = new double [NComp];
	dLnGdP = new double [NComp];
	bj = new double [NComp];
}


void TShvarov::free_internal()
{
	delete[]LnG;
	delete[]dLnGdT;
	delete[]d2LnGdT2;
	delete[]dLnGdP;
	delete[]bj;
}


long int TShvarov::PTparam()
{
	long int j;
	double alp, bet, dal, rho, eps, dedt, d2edt2, dedp;

	// read and copy individual bj parameters
	for (j=0; j<NComp; j++)
	{
		if ( (aDCc[NP_DC*j+1]) > 0 )
			bj[j] = aDCc[NP_DC*j+1];  // individual bj
		else
			bj[j] = 1.;
	}

	// read and convert rho and eps
	rho = RhoW[0];
	alp = - 1./rho*RhoW[1];
	dal = pow(alp,2.) - 1./rho*RhoW[2];
	bet = 1./rho*RhoW[3];
	eps = EpsW[0];
	dedt = 1./eps*EpsW[1];
	d2edt2 = - 1./pow(eps,2.)*pow(EpsW[1],2.) + 1./eps*EpsW[2];  // corrected 23.05.2009 (TW)
	dedp = 1./eps*EpsW[3];

	// calculate A and B terms of Debye-Huckel equation (and derivatives)
	A = (1.82483e6)*sqrt(rho) / pow(Tk*eps,1.5);
	dAdT = - 3./2.*A*( dedt + 1./Tk + alp/3. );
	d2AdT2 = 1./A*pow(dAdT,2.) - 3./2.*A*( d2edt2 - 1/pow(Tk,2.) + 1/3.*dal );
	dAdP = 1./2.*A*( bet - 3.*dedp);
	B = (50.2916)*sqrt(rho) / sqrt(Tk*eps);
	dBdT = - 1./2.*B*( dedt +1./Tk + alp );
	d2BdT2 = 1./B*pow(dBdT,2.) - 1./2.*B*( d2edt2 - 1./pow(Tk,2.) + dal );
	dBdP = 1./2.*B*( bet - dedp );

 	ErrorIf( fabs(A) < 1e-9 || fabs(B) < 1e-9, "Shvarov EDH model",
 			"Error: DH parameter A or B was not calculated - no values of RhoW and EpsW !" );

	// b_gamma constant
	if ( flagElect == 0)
	{
		bgam = bc;
		dbgdT = 0.0;
		d2bgdT2 = 0.0;
		dbgdP = 0.0;
		ao = ac;
		daodT = 0.0;
		d2aodT2 = 0.0;
		daodP = 0.0;
	}

	// b_gamma TP-dependent
	else
	{
		Gfunction();
		BgammaTP();
		IonsizeTP();
		// ao = ac;
		// daodT = 0.0;
		// d2aodT2 = 0.0;
		// daodP = 0.0;
	}

	return 0;
}


long int TShvarov::MixMod()
{
	long int j, k, w;
	double sqI, Z2, lgGam, lnGam, Nw, xw, lnxw, Lgam,
				msum, C, lg_to_ln;
	lg_to_ln = 2.302585093;

	// get index of water (assumes water is last species in phase)
	w = NComp - 1;

	// calculate ionic strength and total molaities (molT and molZ)
	IonicStrength();

	xw = x[w];
	Nw = 1000./18.01528;
	// Lgam = -log10(1.+0.0180153*molT);
	Lgam = log10(xw);  // Helgeson large gamma simplified
	if( Lgam < -0.7 )
		Lgam = -0.7;  // experimental truncation of Lgam to min ln(0.5)
	lnxw = log(xw);
	sqI = sqrt(IS);
	C = (0.5*bgam);

	// loop over species
	for( j=0; j<NComp; j++ )
	{
		lgGam = 0.0;
		lnGam = 0.0;
		msum = 0.0;

		// calculate (bj*mj) sum
		for (k=0; k<(NComp-1); k++)
		{
			msum += bj[k]*m[k];
		}

		// charged species
		if ( z[j] )
		{
			lgGam = 0.0;
			Z2 = z[j]*z[j];
			lgGam = - (A*Z2*sqI) / (1.+ao*B*sqI) + C*bj[j]*msum;
			lnGamma[j] = (lgGam + Lgam) * lg_to_ln;
		}

		// neutral species and water solvent
		else
		{
			// neutral species
			if ( j != (NComp-1) )
			{
				lgGam = 0.0;
				if ( flagNeut == 1 )
					lgGam = C*bj[j]*msum;
				else
					lgGam = 0.0;
				lnGamma[j] = (lgGam + Lgam) * lg_to_ln;
				continue;
			}

			// water solvent
			else
			{
				lgGam = 0.0;
				lnGam = 0.0;
				if ( flagH2O == 1 )
				{
					lgGam = - A/(ao*B) * (2./Nw) * ( IS/(1.+ao*B*sqI) - 2.*sqI/(ao*B) + 2./pow((ao*B),2.) * log(1.+ao*B*sqI) )
								 - C*pow(msum,2.)/(2.*Nw);
				}
				else
					lgGam = 0.0;
				lnGamma[j] = (lgGam) * lg_to_ln;
			}
		}
	}  // j

	return 0;
}


long int TShvarov::ExcessProp( double *Zex )
{
	long int j, k, w;
	double sqI, Z2, Nw, xw, msum, C, dCdT, d2CdT2, dCdP, lg_to_ln,
				g, dgt, d2gt, dgp;
	double U, V, dUdT, dVdT, d2UdT2, d2VdT2, dUdP, dVdP, U1, U2, U3, V1, V2, V3,
				dU1dT, dU2dT, dU3dT, dV1dT, dV2dT, dV3dT, d2U1dT2, d2U2dT2, d2U3dT2,
				d2V1dT2, d2V2dT2, d2V3dT2, dU1dP, dU2dP, dU3dP, dV1dP, dV2dP, dV3dP,
				X, dXdT, d2XdT2, dXdP, Z, dZdT, d2ZdT2, dZdP;
	lg_to_ln = 2.302585093;
	g = 0.; dgt = 0.; d2gt = 0.; dgp = 0.;

	// get index of water (assumes water is last species in phase)
	w = NComp - 1;

	// calculate ionic strength and total molalities
	IonicStrength();

	xw = x[w];
	Nw = 1000./18.01528;
	sqI = sqrt(IS);
	C = (0.5*bgam);
	dCdT = (0.5*dbgdT);
	d2CdT2 = (0.5*d2bgdT2);
	dCdP = (0.5*dbgdP);

	// loop over species
	for( j=0; j<NComp; j++ )
	{
		msum = 0.0;

		// calculate bj*mj sum
		for (k=0; k<(NComp-1); k++)
		{
			msum += bj[k]*m[k];
		}

		// charged species
		if ( z[j] )
		{
			Z2 = z[j]*z[j];
			U = - (Z2*A) * sqI;
			dUdT = - (Z2*dAdT) * sqI;
			d2UdT2 = - (Z2*d2AdT2) * sqI;
			dUdP = - (Z2*dAdP) * sqI;
			V = 1. + (ao*B) * sqI;
			dVdT = ( daodT*B + ao*dBdT ) * sqI;
			d2VdT2 = ( d2aodT2*B + 2.*daodT*dBdT + ao*d2BdT2 ) * sqI;
			dVdP = ( daodP*B + ao*dBdP ) * sqI;
			LnG[j] = ( U/V + C*bj[j]*msum ) * lg_to_ln;
			dLnGdT[j] = ( (dUdT*V - U*dVdT)/pow(V,2.) + dCdT*bj[j]*msum ) * lg_to_ln;
			d2LnGdT2[j] = ( (d2UdT2*V + dUdT*dVdT)/pow(V,2.) - (dUdT*V)*(2.*dVdT)/pow(V,3.)
				- (dUdT*dVdT + U*d2VdT2)/pow(V,2.) + (U*dVdT)*(2.*dVdT)/pow(V,3.)
				+ d2CdT2*bj[j]*msum ) * lg_to_ln;
			dLnGdP[j] = ( (dUdP*V - U*dVdP)/pow(V,2.) + dCdP*bj[j]*msum ) * lg_to_ln;
		}

		// neutral species and water solvent
		else
		{
			// neutral species
			if ( j != (NComp-1) )
			{
				if ( flagNeut == 1 )
				{
					LnG[j] = ( dCdT*bj[j]*msum ) * lg_to_ln;
					dLnGdT[j] = ( dCdT*bj[j]*msum ) * lg_to_ln;
					d2LnGdT2[j] = ( d2CdT2*bj[j]*msum ) * lg_to_ln;
					dLnGdP[j] = ( dCdP*bj[j]*msum ) * lg_to_ln;
				}

				else
				{
					LnG[j] = 0.;
					dLnGdT[j] = 0.;
					d2LnGdT2[j] = 0.;
					dLnGdP[j] = 0.;
				}
				continue;
			}

			// water solvent
			else
			{
				// H2O activity coefficient calculated
				if ( flagH2O == 1 )
				{
					// derivatives of lambda and sigma terms
					U1 = A * IS;
					dU1dT = dAdT * IS;
					d2U1dT2 = d2AdT2 * IS;
					dU1dP = dAdP * IS;
					V1 = (ao*B) + (pow(ao,2.)*pow(B,2.)) * sqI;
					dV1dT = ( daodT*B + ao*dBdT ) + 2.*( ao*daodT*pow(B,2.) + pow(ao,2.)*B*dBdT ) * sqI;
					d2V1dT2 = ( d2aodT2*B + 2.*daodT*dBdT + ao*d2BdT2 )
								+ 2. * ( pow(daodT,2.)*pow(B,2.) + ao*d2aodT2*pow(B,2.) + 4.*ao*daodT*B*dBdT
								+ pow(ao,2.)*pow(dBdT,2.) + pow(ao,2.)*B*d2BdT2 ) * sqI;
					dV1dP = ( daodP*B + ao*dBdP ) + 2.*( ao*daodP*pow(B,2.) + pow(ao,2.)*B*dBdP ) * sqI;

					U2 = (2.*A) * sqI;
					dU2dT = (2*dAdT) * sqI;
					d2U2dT2 = (2.*d2AdT2) * sqI;
					dU2dP = (2.*dAdP) * sqI;
					V2 = pow(ao,2.)*pow(B,2.);
					dV2dT = 2.*( ao*daodT*pow(B,2.) + pow(ao,2.)*B*dBdT );
					d2V2dT2 = 2.*( pow(daodT,2.)*pow(B,2.) + ao*d2aodT2*pow(B,2.) + 4.*ao*daodT*B*dBdT
								+ pow(ao,2.)*pow(dBdT,2.) + pow(ao,2.)*B*d2BdT2 );
					dV2dP = 2.*( ao*daodP*pow(B,2.) + pow(ao,2.)*B*dBdP );

					U3 = 2.*A*log(1.+ao*B*sqI);
					X = log(1.+ao*B*sqI);
					dXdT = pow((1.+ao*B*sqI),-1.) * ( daodT*B + ao*dBdT ) * sqI;
					d2XdT2 = - pow((1.+ao*B*sqI),-2.) * pow(( daodT*B + ao*dBdT ),2.) * IS
								+ pow((1.+ao*B*sqI),-1.) * ( d2aodT2*B + 2.*daodT*dBdT + ao*d2BdT2 ) * sqI;
					dXdP = pow((1.+ao*B*sqI),-1.) * ( daodP*B + ao*dBdP ) * sqI;
					dU3dT = 2.*( dAdT*X + A*dXdT );
					d2U3dT2 = 2.*( d2AdT2*X + 2.*dAdT*dXdT + A*d2XdT2 );
					dU3dP = 2.*( dAdP*X + A*dXdP );
					V3 = pow(ao,3.)*pow(B,3.);
					dV3dT = 3.*( pow(ao,2.)*daodT*pow(B,3.) + pow(ao,3.)*pow(B,2.)*dBdT );
					d2V3dT2 = 3.*( 2.*ao*pow(daodT,2.)*pow(B,3.) + pow(ao,2.)*d2aodT2*pow(B,3.)
								+ 6.*pow(ao,2.)*daodT*pow(B,2.)*dBdT + 2.*pow(ao,3.)*B*pow(dBdT,2.)
								+ pow(ao,3.)*pow(B,2.)*d2BdT2 );
					dV3dP = 3.*( pow(ao,2.)*daodP*pow(B,3.) + pow(ao,3.)*pow(B,2.)*dBdP );

					Z = U1/V1 - U2/V2 + U3/V3;
					dZdT = (dU1dT*V1 - U1*dV1dT)/pow(V1,2.) - (dU2dT*V2 - U2*dV2dT)/pow(V2,2.)
								+ (dU3dT*V3 - U3*dV3dT)/pow(V3,2.);
					d2ZdT2 = (d2U1dT2*V1 + dU1dT*dV1dT)/pow(V1,2.) - (dU1dT*V1)*(2.*dV1dT)/pow(V1,3.)
								- (dU1dT*dV1dT + U1*d2V1dT2)/pow(V1,2.) + (U1*dV1dT)*(2.*dV1dT)/pow(V1,3.)
								- (d2U2dT2*V2 + dU2dT*dV2dT)/pow(V2,2.) + (dU2dT*V2)*(2.*dV2dT)/pow(V2,3.)
								+ (dU2dT*dV2dT + U2*d2V2dT2)/pow(V2,2.) - (U2*dV2dT)*(2.*dV2dT)/pow(V2,3.)
								+ (d2U3dT2*V3 + dU3dT*dV3dT)/pow(V3,2.) - (dU3dT*V3)*(2.*dV3dT)/pow(V3,3.)
								- (dU3dT*dV3dT + U3*d2V3dT2)/pow(V3,2.) + (U3*dV3dT)*(2.*dV3dT)/pow(V3,3.);
					dZdP = (dU1dP*V1 - U1*dV1dP)/pow(V1,2.) - (dU2dP*V2 - U2*dV2dP)/pow(V2,2.)
								+ (dU3dP*V3 - U3*dV3dP)/pow(V3,2.);

					// activity coefficient (and derivatives)
					LnG[j] = ( - (2./Nw)*Z - C*pow(msum,2.)/(2.*Nw) ) * lg_to_ln;
					dLnGdT[j] = ( - (2./Nw)*dZdT - dCdT*pow(msum,2.)/(2.*Nw) ) * lg_to_ln;
					d2LnGdT2[j] = ( - (2./Nw)*d2ZdT2 - d2CdT2*pow(msum,2.)/(2.*Nw) ) * lg_to_ln;
					dLnGdP[j] = ( - (2./Nw)*dZdP - dCdP*pow(msum,2.)/(2.*Nw) ) * lg_to_ln;
				}

				// H2O activity coefficient 0
				else
				{
					LnG[j] = 0.;
					dLnGdT[j] = 0.;
					d2LnGdT2[j] = 0.;
					dLnGdP[j] = 0.;
				}
			}

		}

		g += x[j]*LnG[j];
		dgt += x[j]*dLnGdT[j];
		d2gt += x[j]*d2LnGdT2[j];
		dgp += x[j]*dLnGdP[j];

	}  // j

	// increment thermodynamic properties
	Gex = (R_CONST*Tk) * g;
	Hex = - R_CONST*pow(Tk,2.) * dgt;
	// Sex = - R_CONST * ( g + Tk*dgt );
	Sex = (Hex-Gex)/Tk;
	CPex = - R_CONST * ( 2.*Tk*dgt + pow(Tk,2.)*d2gt );
	Vex = (R_CONST*Tk) * dgp;
	Aex = Gex - Vex*Pbar;
	Uex = Hex - Vex*Pbar;

	// assigments (excess properties)
	Zex[0] = Gex;
	Zex[1] = Hex;
	Zex[2] = Sex;
	Zex[3] = CPex;
	Zex[4] = Vex;
	Zex[5] = Aex;
	Zex[6] = Uex;

	return 0;
}


long int TShvarov::IdealProp( double *Zid )
{
	long int j;
	double si;
	si = 0.0;
	for (j=0; j<NComp; j++)
	{
		if ( x[j] > 1.0e-32 )
			si += x[j]*log(x[j]);
	}
	Hid = 0.0;
	CPid = 0.0;
	Vid = 0.0;
	Sid = (-1.)*R_CONST*si;
	Gid = Hid - Sid*Tk;
	Aid = Gid - Vid*Pbar;
	Uid = Hid - Vid*Pbar;

	// assignments (ideal mixing properties)
	Zid[0] = Gid;
	Zid[1] = Hid;
	Zid[2] = Sid;
	Zid[3] = CPid;
	Zid[4] = Vid;
	Zid[5] = Aid;
	Zid[6] = Uid;
	return 0;
}


long int TShvarov::IonicStrength()
{
	long int j;
	double is, mt;
	is = 0.0; mt = 0.0;

	// calculate ionic strength and total molalities
	// needs to be modified when nonaqueous solvents are present
	for (j=0; j<(NComp-1); j++)
	{
		is += 0.5*m[j]*z[j]*z[j];
		mt += m[j];
	}

	// assignments
	IS = is;
	molT = mt;

	return 0;
}


long int TShvarov::BgammaTP()
{
	// ni: stoichiometric number of moles of ions in one mole of electrolyte
	// rc, ra: radius of cation and anion, respectively at 298 K/1 bar
	// units are cal, kg, K, mol, bar
	double ni, nc, na, zc, za, rc, ra, a1, a2, a3, a4, a5, c1, c2, omg, bg, bs, bh;
	double rec, rea, omgpt, domdt, d2omdt2, domdp, nbg, nbv, nbj, nbh;
	double eps, eta, xborn, yborn, qborn, X1, X2;

	// set parameters
	eps = EpsW[0];
	yborn = 1./pow(EpsW[0],2.)*EpsW[1];
	xborn = 1./pow(EpsW[0],2.) * ( EpsW[2] - 2./EpsW[0]*pow(EpsW[1],2.) );
	qborn = 1./pow(EpsW[0],2.)*EpsW[3];
	eta = (1.66027e5);

	switch ( flagElect )
	{
		case 1:  // NaCl
			ni = 2.; nc = 1.; na = 1.; zc = 1.; za = -1.;
			a1 = 0.030056; a2 = -202.55; a3 = -2.9092; a4 = 20302;
			a5 = -0.206; c1 = -1.50; c2 = 53300.; omg = 178650.;
			bg = -174.623; bs = 2.164; rc = 0.97; ra = 1.81;
			break;
		case 2:  // KCl
			ni = 2.; nc = 1.; na = 1.; zc = 1.; za = -1.;
			a1 = 0.0172; a2 = -115.36; a3 = -1.1857; a4 = 13854.2;
			a5 = -0.262; c1 = -2.53; c2 = 38628.4; omg = 164870.;
			bg = -70.0; bs = 1.727; rc = 1.33; ra = 1.81;
			break;
		case 3:  // NaOH
			ni = 2.; nc = 1.; na = 1.; zc = 1.; za = -1.;
			a1 = 0.030056; a2 = -202.55; a3 = -2.9092; a4 = 20302;
			a5 = -0.206; c1 = -1.50; c2 = 53300.; omg = 205520.;
			bg = -267.4; bs = 1.836; rc = 0.97; ra = 1.40;
			break;
		case 4:  // KOH
			ni = 2.; nc = 1.; na = 1.; zc = 1.; za = -1.;
			a1 = 0.0172; a2 = -115.36; a3 = -1.1857; a4 = 13854.2;
			a5 = -0.262; c1 = -2.53; c2 = 38628.4; omg = 191730.;
			bg = -335.7; bs = 1.26; rc = 1.33; ra = 1.40;
			break;
		default:  // wrong mode
			return -1;
	}

	// calculation part, extended 06.06.2009 (TW)
	bh = bg + (298.15)*bs;
	rec = rc + fabs(zc)*(0.94+Gf);
	rea = ra + fabs(za)*Gf;
	X1 = - eta*nc*( fabs(pow(zc,3.))/pow(rec,2.) - zc/pow((3.082+Gf),2.) )
			- eta*na*( fabs(pow(za,3.))/pow(rea,2.) - za/pow((3.082+Gf),2.) );
	X2 = 2.*eta*nc*( fabs(pow(zc,4.))/pow(rec,3.) - zc/pow((3.082+Gf),3.) )
			+ 2.*eta*nc * ( fabs(pow(zc,4.))/pow(rec,3.) - zc/pow((3.082+Gf),3.) );
	omgpt = eta*( nc*pow(zc,2.)/rec + na*pow(za,2.)/rea );
	// omgpt = (1.66027e5)*(1./(0.94+rc+Gf)+1./(ra+Gf));
	domdt = X1*dGfdT;
	d2omdt2 = X2*pow(dGfdT,2.) + X1*d2GfdT2;
	domdp = X1*dGfdP;
	nbg = - ni*bg/2. + ni*bs*(Tk-298.15)/2. - c1*(Tk*log(Tk/298.15)-Tk+298.15)
				+ a1*(Pbar-1.) + a2*log((2600.+Pbar)/(2600.+1.))
				- c2*((1./(Tk-228.)-1./(298.15-228.))*(228.-Tk)/228.-Tk/(228.*228.)
				* log((298.15*(Tk-228.))/(Tk*(298.15-228.))))
				+ 1./(Tk-228.)*(a3*(Pbar-1.) + a4*log((2600.+Pbar)/(2600.+1.)))
				+ a5*(omgpt*(1./eps-1.)-omg*(1./(78.24513795)-1.)+(-5.798650444e-5)*omg*(Tk-298.15));
	nbh = - ni/2.*bh + c1*(Tk-298.15) - c2*(1./(Tk-228.)-1./(298.15-228.))
				+ a1*(Pbar-1.) + a2*log((2600.+Pbar)/(2600.+1.))
				+ ((2.*Tk-228.)/pow((Tk-228.),2.))*(a3*(Pbar-1.)+a4*log((2600.+Pbar)/(2600.+1.)))
				+ a5*( omgpt*(1./eps-1.) + omgpt*Tk*yborn - Tk*(1./eps-1.)*domdt
				- omg*(1./(78.24513795)-1.) - omg*(298.15)*(-5.798650444e-5) );
	nbj = c1 + c2/pow((Tk-228.),2.)  - ( (2.*Tk)/pow((Tk-228.),3.) )
				* ( a3*(Pbar-1.) + a4*log((2600.+Pbar)/(2600.+1.)) )
				+ a5*( omgpt*Tk*xborn + 2.*Tk*yborn*domdt - Tk*(1./eps-1.)*d2omdt2 );
	nbv = a1 + a2/(2600.+Pbar) + a3/(Tk-228.) + a4/((2600.+Pbar)*(Tk-228.))
				+ a5*(-omgpt*qborn + (1./eps-1.)*domdp);
	// b_gamma = nbg/(2.*log(10.)*1.98721*Tk);
	// bgam = b_gamma;
	bgam = nbg/(2.*log(10.)*(1.98721)*Tk)*2./ni;
	dbgdT = - nbh/(2.*log(10.)*(1.98721)*pow(Tk,2.))*2./ni;
	d2bgdT2 = - nbj/( 2.*log(10.)*(1.98721)*pow(Tk,2.))*2./ni - 2./Tk*dbgdT;
	dbgdP = nbv/(2.*log(10.)*(1.98721)*Tk)*2./ni;

	return 0;
}


long int TShvarov::IonsizeTP()
{
	double nc, na, ni, zc, za, c;

	switch ( flagElect )
	{
		case 1:  // NaCl
			nc = 1.; na = 1.; ni = 2.;
			zc = 1.; za = -1.;
			break;
		case 2:  // KCl
			nc = 1.; na = 1.; ni = 2.;
			zc = 1.; za = -1.;
			break;
		case 3:  // NaOH
			nc = 1.; na = 1.; ni = 2.;
			zc = 1.; za = -1.;
			break;
		case 4:  // KOH
			nc = 1.; na = 1.; ni = 2.;
			zc = 1.; za = -1.;
			break;
		default:  // wrong mode
			return -1;
	}

	c = 2./ni * ( nc*fabs(zc) + na*fabs(za) );
	ao = ac + c*Gf;
	daodT = c*dGfdT;
	d2aodT2 = c*d2GfdT2;
	daodP = c*dGfdP;

	return 0;
}


long int TShvarov::Gfunction()
{
	double T, P, D, beta, alpha, daldT;
	double g, dgdP, dgdT, d2gdT2;
	double TMAX = 1000., PMAX = 5000., TOL = 1.0e-4;

	// convert parameters
	T = Tk - 273.15;
	P = Pbar;
	D = RhoW[0];
	alpha = - RhoW[1]/RhoW[0];
	daldT = pow(alpha, 2.) - RhoW[2]/RhoW[0];
	beta = RhoW[3]/RhoW[0];

	// initialize g and derivatives to zero
	g = 0.0;
	dgdP = 0.0;
	dgdT = 0.0;
	d2gdT2 = 0.0;

	if ((T > TMAX+TOL) || (P > PMAX+TOL))
		return -1;
	else
		GShok2( T, P, D, beta, alpha, daldT, g, dgdP, dgdT, d2gdT2 );

	// assignments
	Gf = g;
	dGfdT = dgdT;
	d2GfdT2 = d2gdT2;
	dGfdP = dgdP;

	return 0;
}


long int TShvarov::GShok2( double T, double P, double D, double beta,
		double alpha, double daldT, double &g, double &dgdP, double &dgdT, double &d2gdT2 )
{
	double a, b, dgdD, /*dgdD2,*/ dadT, dadTT, dbdT, dbdTT, dDdT, dDdP,
		dDdTT, Db, dDbdT, dDbdTT, ft, dftdT, dftdTT, fp, dfpdP,
		f, dfdP, dfdT, d2fdT2, tempy;
	double C[6]  = {-0.2037662e+01,  0.5747000e-02, -0.6557892e-05,
			0.6107361e+01, -0.1074377e-01,  0.1268348e-04 };
	double cC[3] = { 0.3666666e+02, -0.1504956e-09,  0.5017997e-13 };

	if ( D >= 1.3 )
		return -1;
	// Sveta 19/02/2000 1-D < 0 pri D>1 => pow(-number, b) = -NaN0000 error
	double pw = fabs(1.0 - D); // insert Sveta 19/02/2000

	// calculation part
	a = C[0] + C[1]*T + C[2]*pow(T,2.);
	b = C[3] + C[4]*T + C[5]*pow(T,2.);
	g = a * pow(pw, b);

	dgdD = - a*b* pow(pw, (b - 1.0));
		// dgdD2 = a * b * (b - 1.0) * pow((1.0 - D),(b - 2.0));

	dadT = C[1] + 2.0*C[2]*T;
	dadTT = 2.0*C[2];
	dbdT = C[4] + 2.0*C[5]*T;
	dbdTT = 2.0*C[5];
	dDdT = - D * alpha;
	dDdP = D * beta;
	dDdTT = - D * (daldT - pow(alpha,2.));

		// Db = pow((1.0 - D), b);  Fixed by DAK 01.11.00
	Db = pow(pw , b);
	dDbdT = - b*pow(pw,(b - 1.0))*dDdT + log(pw)*Db*dbdT;

	dDbdTT = - (b*pow(pw,(b-1.0)) *dDdTT + pow(pw,(b - 1.0)) * dDdT * dbdT
				+ b * dDdT * (-(b - 1.0) * pow(pw,(b - 2.0)) * dDdT
				+ log(pw) * pow(pw,(b - 1.0)) * dbdT))
				+ log(pw) * pow(pw,b) * dbdTT
				- pow(pw,b) * dbdT * dDdT / (1.0 - D)
				+ log(pw) * dbdT * dDbdT;

	dgdP = dgdD * dDdP;
	dgdT = a * dDbdT + Db * dadT;
	d2gdT2 = a * dDbdTT + 2.0 * dDbdT * dadT + Db * dadTT;

	if((T < 155.0) || (P > 1000.0) || (T > 355.0))
		return 0;

	tempy = ((T - 155.0) / 300.0);
	ft = pow(tempy,4.8) + cC[0] * pow(tempy,16.);

	dftdT = 4.8 / 300.0 * pow(tempy, 3.8) + 16.0 / 300.0 * cC[0] * pow(tempy, 15.);

	dftdTT = 3.8 * 4.8 / (300.0 * 300.0) * pow(tempy, 2.8)
		+ 15.0 * 16.0 / (300.0*300.0) * cC[0] * pow(tempy,14.);

	fp = cC[1] * pow((1000.0 - P),3.) + cC[2] * pow((1000.0 - P),4.);

	dfpdP = -3.0 * cC[1] * pow((1000.0 - P),2.) - 4.0 * cC[2] * pow((1000.0 - P),3.);

	f = ft * fp;
	dfdP = ft * dfpdP;
	dfdT = fp * dftdT;
	d2fdT2 = fp * dftdTT;

	g -= f;
	dgdP -= dfdP;
	dgdT -= dfdT;
	d2gdT2 -= d2fdT2;

	return 0;
}


//--------------------- End of s_solmod4.cpp ---------------------------


