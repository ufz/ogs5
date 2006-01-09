/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <iostream>
#include <math.h>
#include <vector>
#include <cstdlib>
#include <cmath>
#include "IAPWS-IF97.h"
#include "Density.h"
using namespace std;


double density::Ds, density::Vs;

density::density(void){}
density::~density(void){}




void density::MaoModel(double T, double P, double m, int f){

	double R=83.14472, N0=6.0221367E+23, pai=3.14159265, qe=1.60217733E-19, kB=1.380658E-23, D0=8.85418781762E-12;
	double dP=1.0e-2, dH2O, dH2Oa, dH2Ob, Av, D1000, C, B, DD;

	dH2O =IF97::density(T,P/10.0);
	dH2Oa=IF97::density(T,(P-dP)/10.0);
	dH2Ob=IF97::density(T,(P+dP)/10.0);

	const double U[10]=
		{ 0.0, 3.4279e+2, -5.0866e-3, 9.4690e-7, -2.0525, 3.1159e+3, -1.8289e+2, -8.0325e+3, 4.2142e+6, 2.1417 };
	D1000 = U[1]*exp(U[2]*T + U[3]*T*T);
	C = U[4]+U[5]/(U[6]+T);
	B = U[7]+U[8]/T+U[9]*T;
	DD= D1000+C*log((B+P)/(B+1000));

	//Aphia = 1.0/3.0*pow(2.0*pai*N0*dH2Oa,0.5) * pow(qe*qe/(4.0*pai*D0*Da*kB*T),1.5);
	//Aphib = 1.0/3.0*pow(2.0*pai*N0*dH2Ob,0.5) * pow(qe*qe/(4.0*pai*D0*Db*kB*T),1.5);
	//Av = -2.0*R*T*(Aphib-Aphia)/dP/10.0;
	Av=  -4.0*R*T  /3.0* pow(2.0*pai*N0,0.5) * pow(qe*qe/(4.0*pai*D0*kB*T),1.5)
	    *(pow(dH2O,0.5)*(-1.5*pow(DD,-2.5)*C/(B+P))+pow(DD,-1.5)*0.5*pow(dH2O,-0.5)*(dH2Ob-dH2Oa)/(2.0*dP));
	//cout << " Av---- " << Av << endl;

	const double c[7][24]={
	{0.0, 1.17271480E+03,  1.40527916E-01, -5.53962649E-04,  1.72402126E-06, -1.57184556E-01,  8.89959461E-04,
	     -1.52090064E-06,  0.0,             7.33553879E-04,  4.06701494E-05, -2.38873863E-07,  3.94900863E-10,
		 -3.56131664E-03, -6.18472877E-07,  3.00484214E-08, -1.02229075E-10,  0.0,             2.35592818E-06,
		 -2.68117086E-05, -2.17228726E-06,  1.19732313E-08, -1.51104808E-11,  3.83403994E-05},
	{0.0, 1.06607098E+03, -8.39622456E-03,  5.35429127E-04,  7.55373789E-07, -4.19512335E-02,  1.45082899E-04,
	     -3.47807732E-07,  0.0,             1.10913788E-03,  1.14498252E-04, -5.51181270E-07,  7.05483955E-10,
		 -5.05734723E-03, -1.32747828E-06,  4.77261581E-08, -1.76888377E-10,  0.0, 		       6.40541237E-06,
		  3.07698827E-05, -1.64042763E-05,  7.06784935E-08, -6.50338372E-11, -4.50906014E-05},
	{0.0, 2.90812061E+02,  6.54111195E+00, -1.61831978E-02,  1.46280384E-05,  1.41397987E+00, -1.07266230E-02,
	      2.64506021E-05, -2.19789708E-08,  3.02182158E-03, -2.15621394E-04,  9.24163206E-07, -1.10089434E-09,
		  2.87018859E-03, -6.73119697E-06,  1.68332473E-06, -7.99645640E-09,  1.11881560E-11, -6.59292385E-05,
		 -2.02369103E-04, -1.70609099E-05,  1.00510108E-07, -1.86624642E-10,  1.91919166E-03},
	{0.0, 1.18880927E+03, -1.43194546E+00,  3.87973220E-03, -2.20330377E-06,  6.38745038E-01, -5.51728055E-03,
	      1.50231562E-05, -1.35757912E-08,  8.43627549E-04,  5.25365072E-04, -1.87204100E-06,  4.20263897E-09,
		 -1.18062548E-01,  6.07424747E-06, -1.20268210E-06,  5.23784551E-09, -8.23940319E-12,  9.75167613E-05,
		 -4.92959181E-05, -2.73642775E-05, 5.42602386E-08, 	-1.95602825E-10,  1.00921935E-02},
	{0.0, 1.12080057E+03, -2.61669538E-01,  1.52042960E-03, -6.89131095E-07, -5.11802652E-02,  2.22234857E-04,
	     -5.66464544E-07,  2.92950266E-10,  2.43934633E-03, -1.42746873E-04,  7.35840529E-07, -9.43615480E-10,
		 -5.18606814E-03, -6.16536928E-07, -1.04523561E-07,  4.52637296E-10, -1.05076158E-12,  2.31544709E-05,
		 -1.09663211E-04,  1.90836111E-05, -9.25997994E-08,  1.54388261E-10, -1.29354832E-03},
	{0.0, 1.11894213E+03, -7.37321458E-01,  1.77908655E-03,  0.0,             0.0,             0.0,
	      0.0,             0.0,             2.21225680E-03,  6.62517291E-05, -2.37296050E-07,  0.0,
		  0.0,             0.0,             0.0,             0.0,             0.0,             0.0,
		 -4.21300430E-04,  0.0,             9.46738388E-09,  0.0,             0.0},
	{0.0, 1.10229139E+03, -7.53497776E-01,  1.92829036E-03,  0.0,             0.0,            -1.15406910E-04,
	      0.0,             0.0,    			2.57437715E-03,  1.64541676E-05, -9.30035886E-08,  0.0,
		  0.0,             0.0,             0.0,             0.0,             0.0,             0.0,
		  0.0,             0.0,             0.0,             0.0,             0.0},
	};
	const double mr[7]={ 10.0,6.0,6.0,2.0,5.0,2.0,1.5 };
	const double Ms[7]={ 42.394, 58.443, 74.551, 95.236, 110.986, 158.536, 206.286 };

	double Imx, Imrx, hIx, hIrx, CA, CB, CC;
	if(f==0 || f==1 || f==2){
		Imx = m;
		Imrx= mr[f];
		CA=2.0;
		CB=2.0;
		CC=2.0;
	}
	else {
		Imx = 3.0*m;
		Imrx= 3.0*mr[f];
		CA=6.0;
		CB=4.0;
		CC=8.0;
	}
	hIx =log(1.0+1.2*pow(Imx,0.5) )/(2.0*1.2);
	hIrx=log(1.0+1.2*pow(Imrx,0.5))/(2.0*1.2);

	if(m==0) density::Ds=dH2O/1000.0;
	else     density::Ds=
			(1000+m*Ms[f])/m/ ((1/m-1/mr[f])*1.0e6/dH2O + CA*Av*( hIx-hIrx )
			  +1.0/mr[f]*(c[f][1]+c[f][2]*T+c[f][3]*pow(T,2.0)+c[f][4]*pow(T,3.0)
					+P*(c[f][5]+c[f][6]*T+c[f][7]*pow(T,2.0)+c[f][8]*pow(T,3.0)))
			  +CB*R*T*(m-mr[f])*(c[f][9]/(T-227.0)+c[f][10]+c[f][11]*T
			  +c[f][12]*pow(T,2.0)+c[f][13]/(647.0-T)+P*(c[f][14]/(T-227.0)+c[f][15]+c[f][16]*T
			  +c[f][17]*pow(T,2.0)+c[f][18]/(647-T)))
			  +CC*R*T*(m*m-pow(mr[f],2.0))*(c[f][19]/(T-227.0)+c[f][20]+c[f][21]*T
			  +c[f][22]*pow(T,2.0) + c[f][23]/(647.0-T)));
	density::Vs=
			  1.0/mr[f]*(c[f][1]+c[f][2]*T+c[f][3]*pow(T,2.0)+c[f][4]*pow(T,3.0)
					+ P*(c[f][5]+c[f][6]*T+c[f][7]*pow(T,2.0)+c[f][8]*pow(T,3.0)))
			 - 1.0e6/mr[f]/dH2O-CA*Av*hIrx
			 -CB*R*T*mr[f]*(c[f][9] /(T-227.0)+c[f][10]+c[f][11]*T+c[f][12]*pow(T,2.0)+c[f][13]/(647.0-T)
			    		+P*(c[f][14]/(T-227.0)+c[f][15]+c[f][16]*T+c[f][17]*pow(T,2.0)+c[f][18]/(647.0-T)))
	-CC*R*T*pow(mr[f],2.0)*(c[f][19]/(T-227.0)+c[f][20]+c[f][21]*T+c[f][22]*pow(T,2.0)+c[f][23]/(647.0-T));


}

double density::MultiDensity(double T, double P, vector<double> mv, vector<int> fv){

	const double Ms[7]={ 42.394, 58.443, 74.551, 95.236, 110.986, 158.536, 206.286 };
	double mz=0,mc=0,mc0=0,Mt=1000.0,Vt=0,Vex,Vw;
	int i,fr=1,n=(int)fv.size();

	density::MaoModel(T,P,0.0,0);
	Vw=1000.0/density::Ds;
	Vt=Vw;

	for(i=0;i<n;i++){
		if(fv[i]==0 || fv[i]==1 || fv[i]==2) mc=mv[i];
		else mc=3.0*mv[i];
		if(mc>mc0) fr=fv[i];
		mc0=mc;
		mz+=mc;
		Mt+=mv[i]*Ms[fv[i]];
		density::MaoModel(T,P,mv[i],fv[i]);
		//cout << " Ds " << tdensity::Ds << endl;
		Vt+=mv[i]*density::Vs;
	}

	if(fr!=0 && fr!=1 && fr!=2) mz /=3.0;
	density::MaoModel(T,P,mz,fr);
	Vex=(1000.0+mz*Ms[fr])/density::Ds - Vw - mz*density::Vs;
	Vt+=Vex;
	return Mt/Vt;
}



double density::ExcessM(double mv[]){

	double mex=0.0;
	mex += mv[0];
	mex += mv[1];
	mex += mv[2];
	mex += 3.0*mv[3];
	mex += 3.0*mv[4];
	mex += mv[5];
	mex += 3.0*mv[6];
	mex += 3.0*mv[7];
	return 0.5*mex;
}


double density::Multi_density(double T, double P, double mv[]){

	double Vs[8];
	const double Ms[8]={ 6.941, 22.990, 39.098, 24.305, 40.080, 35.453, 96.058, 60.009  }; //moleculer weight
	double Mt=1000.0,Vt=0,Vex,Vw,mex;
	int i;

	//0-Li, 1-Na, 2-K, 3-Mg, 4-Ca, 5-Cl, 6-SO4, 7-CO3
	//calc Partial Moler Volume
	//ref solute KCl, set Vs_K = Vs_Cl
	density::MaoModel(T,P,0.0,2);
	Vs[2]=0.5*density::Vs;
	Vs[5]=0.5*density::Vs;

	//Vs_Li = Vs_LiCl - Vs_Cl
	density::MaoModel(T,P,0.0,0);
	Vs[0]=density::Vs-Vs[5];

	//Vs_Na
	density::MaoModel(T,P,0.0,1);
	Vs[1]=density::Vs-Vs[5];


	//Vs_Mg = Vs_Na
	Vs[3]=Vs[1];

	//Vs_Ca = Vs_K
	Vs[4]=Vs[2];

	Vs[6]=1.5*Vs[5];
	Vs[7]=1.5*Vs[5];


	//=========

	density::MaoModel(T,P,0.0,0);
	Vw=1000.0/density::Ds;
	Vt=Vw;

	for(i=0;i<=7;i++){
		Mt+=mv[i]*Ms[i];
		Vt+=mv[i]*Vs[i];
	}

	mex=density::ExcessM(mv);
	density::MaoModel(T,P,mex,1);
	Vex=(1000.0+mex*(Ms[1]+Ms[5]))/density::Ds - Vw - mex*density::Vs;

	//cout << endl << mex<< " " << (1000.0+mex*(Ms[1]+Ms[5]))/density::Ds << " " << Vw<< " " << Vex << endl;

	Vt+=Vex;

	//cout << " " << Mt << " " << Vt << endl;

	return Mt/Vt;

}




double density::CO2brine(double T, double P, double mNaCl, double mCO2)
{
	double Ps /*,mTot*/,V2,Dens;

	if(mNaCl<0.0 || mNaCl>6.0 || mCO2<0.0){
		Dens = -1.0;
		goto T01;
	}
	if(T<273.15 || T>573.15 || P<0.001 || P>1000.0){
		Dens = -1.0;
		goto T01;
	}

	Ps=IF97::Psat(T)*10.0;
	if(Ps>P){
		Dens = IF97::density(T,P/10.0)/1000.0;
		cout << " Ps " << endl;
		goto T01;
	}

	// mTot = 55.5084 +mNaCl +mCO2;

	density::MaoModel(T, P, mNaCl,1);

	V2 = 18.0153/IF97::density(T,P/10.0)*1000.0*
		(1.0 +   0.38384020e-3*T*T - 0.55953850   *T + 0.30429268e3 - 0.72044305e5/T + 0.63003388e7/T/T
	         + (-0.57709332e-5*T*T + 0.82764653e-2*T - 0.43813556e1 + 0.10144907e4/T - 0.86777045e5/T/T)*P/10.0);

	Dens = (1000.0+58.4428*mNaCl+44.0098*mCO2)/((1000.0+58.4428*mNaCl)/density::Ds+V2*mCO2);

T01:
	return Dens;
}


//0-Li, 1-Na, 2-K, 3-Mg, 4-Ca, 5-Cl, 6-SO4, 7-CO3
double density::CO2_MultiBrine_density(double T, double P, double mv[], double mCO2)
{
	double Ps,V2,Dens,Dens_B,Mt=1000.0;
	const double Ms[8]={ 6.941, 22.990, 39.098, 24.305, 40.080, 35.453, 96.058, 60.009  }; //moleculer weight
	int i;

	Ps=IF97::Psat(T)*10.0;
	if(Ps>P){
		Dens = IF97::density(T,P/10.0)/1000.0;
		goto T01;
	}

	Dens_B=density::Multi_density(T,P,mv);
	for(i=0;i<=7;i++)
		Mt += mv[i]*Ms[i];

	V2 = 18.0153/IF97::density(T,P/10.0)*1000.0*
		(1.0 +   0.38384020e-3*T*T - 0.55953850   *T + 0.30429268e3 - 0.72044305e5/T + 0.63003388e7/T/T
	         + (-0.57709332e-5*T*T + 0.82764653e-2*T - 0.43813556e1 + 0.10144907e4/T - 0.86777045e5/T/T)*P/10.0);

	Dens = (Mt+44.0098*mCO2)/(Mt/Dens_B+V2*mCO2);

T01:
	return Dens;
}


//0-Li, 1-Na, 2-K, 3-Mg, 4-Ca, 5-Cl, 6-SO4, 7-CO3
double density::concentration_water(double T, double P, double cv[], double cv_CO2) //unit mol/m3
{
	const double Ms[8]={ 6.941, 22.990, 39.098, 24.305, 40.080, 35.453, 96.058, 60.009  }; //moleculer weight
	double mv[8],mCO2,dens,d0,d1,mass,densx,err=1.0e-4;
	int i,ix;

	d0=0.8;
	d1=1.5;

	for(ix=0;ix<8;ix++){
		dens=0.5*(d0+d1); //unit g/cm3
		mass=dens*1.0e6;  //unit (g)

		for(i=0;i<8;i++)
			mass -= cv[i]*Ms[i];
		mass -= cv_CO2*44.01;

		for(i=0;i<8;i++){
			mv[i]= 1000.0*cv[i]/mass;
			//cout << mv[i] << " " ;
		}
		mCO2= 1000.0*cv_CO2/mass;
		//cout << mCO2 << endl;

		densx=density::CO2_MultiBrine_density(T,P,mv,mCO2);
		//cout << " dens " << dens << " densx " << densx << endl;
		if( abs(dens-densx) < err )
			break;
		else if(dens>densx)
			d1=dens;
		else
			d0=dens;

	}
	return mass/18.015;
}
