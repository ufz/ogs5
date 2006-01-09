/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <iostream>
#include <limits>
#include <cstdlib>
#include <cmath>
#include "IAPWS-IF97.h"
#include "NR.h"
#include "Brent/brent.hpp"
using namespace std;


double IF97::R, IF97::Rm, IF97::M, IF97::Tc, IF97::Pc, IF97::Dc, IF97::Tt, IF97::Pt, IF97::Tb;
double IF97::TT, IF97::PP;
//-----------------------------
IF97::IF97(void){}
IF97::~IF97(void){}

void IF97::ReferenceConstants(void){
	R  = 0.461526  ;//kJ kg-1 K-1,   (1.1)
	Rm = 8.31451   ;//kJ kmol-1 K-1, (1.2)
	M  = 18.015257 ;//kg kmol-1,     (1.3)
	Tc = 647.096   ;//K ,            (1.4)
	Pc = 22.064    ;//MPa ,          (1.5)
	Dc = 322.0     ;//kg m-3,        (1.6)
	Tt = 273.16    ;//K ,            (1.7)
	Pt = 611.657   ;//Pa ,           (1.8)
	Tb = 373.1243  ;//K ,            (1.9)
}


double IF97::Psat(double T){
    double A,B,C,theta,n1,n2,n3,n4,n5,n6,n7,n8,n9,n10;
	n1 =  0.11670521452767e+04;
	n2 = -0.72421316703206e+06;
	n3 = -0.17073846940092e+02;
	n4 =  0.12020824702470e+05;
	n5 = -0.32325550322333e+07;
	n6 =  0.14915108613530e+02;
	n7 = -0.48232657361591e+04;
	n8 =  0.40511340542057e+06;
	n9 = -0.23855557567849e+00;
	n10=  0.65017534844798e+03;
	theta = T+n9/(T-n10);
	A =    pow(theta,2.0)+n1*theta+n2;
	B = n3*pow(theta,2.0)+n4*theta+n5;
	C = n6*pow(theta,2.0)+n7*theta+n8;
	return pow(2.0*C/(-B+pow(B*B-4.0*A*C,0.5)),4.0);
}

double IF97::Tsat(double P){
    double D,E,F,G,beta,n1,n2,n3,n4,n5,n6,n7,n8,n9,n10;
	n1 =  0.11670521452767e+04;
	n2 = -0.72421316703206e+06;
	n3 = -0.17073846940092e+02;
	n4 =  0.12020824702470e+05;
	n5 = -0.32325550322333e+07;
	n6 =  0.14915108613530e+02;
	n7 = -0.48232657361591e+04;
	n8 =  0.40511340542057e+06;
	n9 = -0.23855557567849e+00;
	n10=  0.65017534844798e+03;
	beta=pow(P,0.25);
	E =    pow(beta,2.0)+n3*beta+n6;
	F = n1*pow(beta,2.0)+n4*beta+n7;
	G = n2*pow(beta,2.0)+n5*beta+n8;
	D = 2.0*G/(-1.0*F-pow(pow(F,2.0)-4.0*E*G,0.5));
	return 0.5*(n10+D-pow(pow(n10+D,2.0)-4.0*(n9+n10*D),0.5));
}

double IF97::Pb23(double T){
	const double B23_N[6]
	= { 0, 0.34805185628969E+03, -0.11671859879975E+01, 0.10192970039326E-02, 0.57254459862746E+03, 0.13918839778870E+02 };
	return B23_N[1]+B23_N[2]*T+B23_N[3]*T*T;
}

double IF97::Tb23(double P){
	const double B23_N[6]
	= { 0, 0.34805185628969E+03, -0.11671859879975E+01, 0.10192970039326E-02, 0.57254459862746E+03, 0.13918839778870E+02 };
	return B23_N[4]+pow((P-B23_N[5])/B23_N[3],0.5);
}

int IF97::region(double T, double P){
	double Ps,P23;
	if (T<=623.15){
		Ps=Psat(T);
		if(P==Ps) return 4;
		else if(P>Ps) return 1;
		         else return 2;
	}
	else{
		P23=Pb23(T);
		if(P>P23) return 3;
		else return 2;
	}
}




double IF97::g1PT(double P, double T){

	// REGION 1 CORRELATION DATA
	#define REG1_COUNT 34
	const int REGION1_I[REG1_COUNT] = {
										  0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4,
										  4, 4, 5, 8, 8, 21, 23, 29, 30, 31, 32
									  };
	const int REGION1_J[REG1_COUNT] = {
										  -2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0, 1, 3, 17, -4,
										  0, 6, -5, -2, 10, -8, -11, -6, -29, -31, -38, -39, -40, -41
									  };
	const double REGION1_N[REG1_COUNT] = {
				0.14632971213167E+00, -0.84548187169114E+00, -0.37563603672040E+01,
				0.33855169168385E+01, -0.95791963387872E+00, 0.15772038513228E+00, -0.16616417199501E-01,
				0.81214629983568E-03, 0.28319080123804E-03, -0.60706301565874E-03, -0.18990068218419E-01,
				-0.32529748770505E-01, -0.21841717175414E-01, -0.52838357969930E-04, -0.47184321073267E-03,
				-0.30001780793026E-03, 0.47661393906987E-04, -0.44141845330846E-05, -0.72694996297594E-15,
				-0.31679644845054E-04, -0.28270797985312E-05, -0.85205128120103E-09, -0.22425281908000E-05,
				-0.65171222895601E-06, -0.14341729937924E-12, -0.40516996860117E-06, -0.12734301741641E-08,
				-0.17424871230634E-09, -0.68762131295531E-18, 0.14478307828521E-19, 0.26335781662795E-22,
				-0.11947622640071E-22, 0.18228094581404E-23, -0.93537087292458E-25
			};

	int i;
	double R=0.461526,Px=16.53,Tx=1386.0,Pi,Tau,res=0;
	Pi = P/Px;
	Tau= Tx/T;
	for(i=0;i<34;i++) res += REGION1_N[i]*pow(7.1-Pi,REGION1_I[i])*pow(Tau-1.222,REGION1_J[i]);
	return res*R*T;
}

double IF97::g2PT(double P, double T){

	// REGION 2 CORRELATION DATA
	// Ideal Gas Series Data
	#define REG2I_COUNT 9
	const int REGION2_J0[REG2I_COUNT] = { 0, 1, -5, -4, -3, -2, -1, 2, 3 };
	const double REGION2_N0[REG2I_COUNT] = {
				-0.96927686500217E+01, 0.10086655968018E+02, -0.56087911283020E-02,
				0.71452738081455E-01, -0.40710498223928E+00, 0.14240819171444E+01, -0.43839511319450E+01,
				-0.28408632460772E+00, 0.21268463753307E-01
			};
	// Residual Series Data
	#define REG2R_COUNT 43
	const int REGION2_I[REG2R_COUNT] = {
										   1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 5, 6, 6, 6, 7,
										   7, 7, 8, 8, 9, 10, 10, 10, 16, 16, 18, 20, 20, 20, 21, 22, 23, 24, 24, 24
									   };
	const int REGION2_J[REG2R_COUNT] = {
										   0, 1, 2, 3, 6, 1, 2, 4, 7, 36, 0, 1, 3, 6, 35, 1, 2, 3, 7, 3, 16, 35,
										   0, 11, 25, 8, 36, 13, 4, 10, 14, 29, 50, 57, 20, 35, 48, 21, 53, 39, 26, 40, 58
									   };
	const double REGION2_N[REG2R_COUNT] = {
				-0.17731742473213E-02, -0.17834862292358E-01, -0.45996013696365E-01,
				-0.57581259083432E-01, -0.50325278727930E-01, -0.33032641670203E-04, -0.18948987516315E-03,
				-0.39392777243355E-02, -0.43797295650573E-01, -0.26674547914087E-04, 0.20481737692309E-07,
				0.43870667284435E-06, -0.32277677238570E-04, -0.15033924542148E-02, -0.40668253562649E-01,
				-0.78847309559367E-09, 0.12790717852285E-07, 0.48225372718507E-06, 0.22922076337661E-05,
				-0.16714766451061E-10, -0.21171472321355E-02, -0.23895741934104E+02, -0.59059564324270E-17,
				-0.12621808899101E-05, -0.38946842435739E-01, 0.11256211360459E-10, -0.82311340897998E+01,
				0.19809712802088E-07, 0.10406965210174E-18, -0.10234747095929E-12, -0.10018179379511E-08,
				-0.80882908646985E-10, 0.10693031879409E+00, -0.33662250574171E+00, 0.89185845355421E-24,
				0.30629316876232E-12, -0.42002467698208E-05, -0.59056029685639E-25, 0.37826947613457E-05,
				-0.12768608934681E-14, 0.73087610595061E-28, 0.55414715350778E-16, -0.94369707241210E-06
			};

	int i;
	double R=0.461526,Px=1.0,Tx=540.0,Pi,Tau,resI=0,resR=0;
	Pi = P/Px;
	Tau= Tx/T;
	resI=log(Pi);
	for(i=0;i<9; i++) resI += REGION2_N0[i]*pow(Tau,REGION2_J0[i]);
	for(i=0;i<43;i++) resR += REGION2_N[i] *pow(Pi,REGION2_I[i])*pow(Tau-0.5,REGION2_J[i]);
	return (resI+resR)*R*T;
}

double IF97::f3DT(double D, double T){

	// CORRELATION DATA
	#define REG3_COUNT 40
	const int REGION3_I[REG3_COUNT] = {
										  0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3,
										  4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 8, 9, 9, 10, 10, 11
									  };
	const int REGION3_J[REG3_COUNT] = {
										  0, 0, 1, 2, 7, 10, 12, 23, 2, 6, 15, 17, 0, 2, 6, 7, 22, 26, 0, 2, 4,
										  16, 26, 0, 2, 4, 26, 1, 3, 26, 0, 2, 26, 2, 26, 2, 26, 0, 1, 26
									  };
	const double REGION3_N[REG3_COUNT] = {
				0.10658070028513E+01, -0.15732845290239E+02, 0.20944396974307E+02,
				-0.76867707878716E+01, 0.26185947787954E+01, -0.28080781148620E+01, 0.12053369696517E+01,
				-0.84566812812502E-02, -0.12654315477714E+01, -0.11524407806681E+01, 0.88521043984318E+00,
				-0.64207765181607E+00, 0.38493460186671E+00, -0.85214708824206E+00, 0.48972281541877E+01,
				-0.30502617256965E+01, 0.39420536879154E-01, 0.12558408424308E+00, -0.27999329698710E+00,
				0.13899799569460E+01, -0.20189915023570E+01, -0.82147637173963E-02, -0.47596035734923E+00,
				0.43984074473500E-01, -0.44476435428739E+00, 0.90572070719733E+00, 0.70522450087967E+00,
				0.10770512626332E+00, -0.32913623258954E+00, -0.50871062041158E+00, -0.22175400873096E-01,
				0.94260751665092E-01, 0.16436278447961E+00, -0.13503372241348E-01, -0.14834345352472E-01,
				0.57922953628084E-03, 0.32308904703711E-02, 0.80964802996215E-04, -0.16557679795037E-03,
				-0.44923899061815E-04
			};

	int i;
	double R=0.461526,Dx=322.0,Tx=647.096,Delta,Tau,res=0;
	Delta=D/Dx;
	Tau=Tx/T;
	res=REGION3_N[0]*log(Delta);
	for(i=1;i<40;i++) res += REGION3_N[i]*pow(Delta,REGION3_I[i])*pow(Tau,REGION3_J[i]);
	return res*R*T;
}

double IF97::dpressure(double ds){
	return IF97::PP-ds*ds*NR::dfridrX(IF97::f3DT, ds, IF97::TT)*1.0e-3;
}


//Gibbs energy (kJ/kg) (T K, P Mpa)
double IF97::G(double T, double P){
	int rg=region(T,P);
	double dens, res=0.0;
	if(rg==1 || rg==4) res=IF97::g1PT(P,T);
	if(rg==2)          res=IF97::g2PT(P,T);
	if(rg==3) {
		dens=IF97::density(T,P);
		res=IF97::f3DT(dens,T)+dens*NR::dfridrX(f3DT, dens, T);
	}
	return res;
}
//Specific enthalpy (kJ/kg) (T K, P Mpa)
double IF97::H(double T, double P){
	int rg=region(T,P);
	double dens,res=0.0;
	if(rg==1 || rg==4) res=IF97::g1PT(P,T)-T*NR::dfridrY(g1PT, P, T);
	if(rg==2)          res=IF97::g2PT(P,T)-T*NR::dfridrY(g2PT, P, T);
	if(rg==3) {
		dens=IF97::density(T,P);
		res=IF97::f3DT(dens,T)-T*NR::dfridrY(f3DT, dens, T)+dens*NR::dfridrX(f3DT, dens, T);
	}
	return res;
}
//Specific entropy (kJ/kg K) (T K, P Mpa)
double IF97::S(double T, double P){
	int rg=region(T,P);
	double dens,res=0.0;
	if(rg==1 || rg==4) res= -1.0*NR::dfridrY(g1PT, P, T);
	if(rg==2)          res= -1.0*NR::dfridrY(g2PT, P, T);
	if(rg==3) {
		dens=IF97::density(T,P);
		res= -1.0*NR::dfridrY(f3DT, dens, T);
	}
	return res;
}


//density (kg m^-3) (T K, P Mpa)
double IF97::density(double T, double P){
	int rg=region(T,P);
	double res=0.0;
	if(rg==1 || rg==4) res=1.0/NR::dfridrX(g1PT, P, T)*1.0e3;
	if(rg==2)          res=1.0/NR::dfridrX(g2PT, P, T)*1.0e3;
	if(rg==3) {
		IF97::TT=T;
		IF97::PP=P;
		brent::local_min(85.0, 750.0, 1e-6, dpressure, res);
	}
	return res;
}

//viscosity (Pa s)  (T K, P Mpa)
double IF97::viscosity(double T, double P){

	const int VISC_I[21] = {  0, 0, 0, 0, 1 ,1 ,1, 1, 1, 2, 2, 2 ,2, 2, 3, 3, 4, 4, 5, 6, 6  };
	const int VISC_J[21] = {  0, 1, 2, 3, 0, 1, 2, 3, 5, 0, 1, 2, 3, 4, 0, 1, 0, 3, 4, 3, 5  };
	const double VISC_N[21] = {
		0.520094, 0.850895e-1, -0.108374e1, -0.289555, 0.222531, 0.999115, 0.188797e1, 0.126613e1,
		0.120573, -0.281378, -0.906851, -0.772479, -0.489837, -0.257040, 0.161913, 0.257399, -0.325372e-1,
		0.698452e-1, 0.872102e-2, -0.435673e-2, -0.593264e-3
	};
	const double VISC_N0[5] = { 0, 0.167752e-1, 0.220462e-1, 0.6366564e-2, -0.241605e-2 };

	ReferenceConstants();
	int i;
	double Delta, Theta, res0=0, res1=0;
	Delta=density(T,P)/Dc;
	Theta=T/Tc;
	for(i=1;i<=4;i++) res0 += VISC_N0[i]*pow(Theta, 1-i);
	res0 = pow(Theta,0.5)*pow(res0,-1.0);
	Theta=Tc/T;
	for(i=0;i<21;i++) res1 += VISC_N[i]*pow(Delta-1.0, VISC_I[i])*pow(Theta-1.0,VISC_J[i]);
	res1 = exp(Delta*res1);
	return res0*res1*1.0e-6;
}

//dielectric constant ()
double IF97::dielectric(double T, double P){
	const int DIEL_I[11] = {  1, 1, 1, 2, 3, 3, 4, 5, 6, 7, 10  };
	const double DIEL_J[11] = {  0.25, 1, 2.5, 1.5, 1.5, 2.5, 2, 2, 5, 0.5, 10  };
	const double DIEL_N[12] = {
		0.978224486826, -0.957771379375, 0.237511794148, 0.714692244396, -0.298217036956, -0.108863472196,
		0.949327488264e-1, -0.980469816509e-2, 0.165167634970e-4, 0.937359795772e-4, -0.123179218720e-9, 0.196096504426e-2
	};
	int i;
	double k, Na, a, eps0, u, Mw;
	double A, B, g, rho, Delta, Tau, eps;
	k    = 1.380658e-23; //J K-1  Boltzmann s constant k
	Na   = 6.0221367e23; //mol-1 Avogadro s number NA
	a    = 1.636e-40; //C2 J-1 m2 Mean molecular polarizability a
	eps0 = 8.854187817e-12; //C2 J-1 m-1 Permittivity of vacuum e 0
	u    = 6.138e-30; //C m Molecular dipole moment µ
	Mw   = 0.018015268; //kg mol-1 Molar mass M
	rho = density(T,P);
	ReferenceConstants();
	Delta= rho/Dc;
	Tau  = Tc/T;
	g=1.0+DIEL_N[11]*Delta*pow((Tc/228.0/Tau-1),-1.2);
	for(i=0;i<11;i++) g += DIEL_N[i]*pow(Delta,DIEL_I[i])*pow(Tau,DIEL_J[i]);
	A = Na*pow(u,2.0)*rho*g/(Mw*eps0*k*T);
	B = Na*a*rho/(3.0*Mw*eps0);
	eps=(1.0+A+5.0*B+pow((9.0+2.0*A+18.0*B+pow(A,2.0)+10.0*A*B+9.0*pow(B,2.0)),0.5))/(4.0*(1.0-B));
	return eps;
}
