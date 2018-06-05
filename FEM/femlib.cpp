/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   ROCKFLOW - Modul: femlib.c

   Aufgabe:
   Funktionen fuer die FEM

   Programmaenderungen:
   02/1997     R.Kaiser     Uebertragen aus dem aTM
   04/1997     R.Kaiser     Funktionsnamen geaendert
   05/1997     R.Kaiser     jac_asm.c -> femlib.c
   16.10.1997  R.Kaiser     Funktionen zur Berechnung der Elementvolumen
   09.03.1998  R.Kaiser     Funktionen zur Berechnung eines
   Einheitsnormalenvektors zur Elementkante (2D)
   02/2000     RK           Funktionen fuer zweidim. 9-Knoten-Elemente
   02/2000     C.Thorenz    Wrapper fuer InterpolValuexD
   11/2000     OK           Volumenberechnung fuer Tetraeder und Prismen
   03/2003     RK           Quellcode bereinigt, Globalvariablen entfernt

**************************************************************************/
#include "femlib.h"

/***************************************************************************
   ROCKFLOW - Funktion: MXPGaussPktTri
   Aufgabe:
   X Punkt Gauss-Integration Dreieck
           bestimmtes Integral
           1/           X
 | P(x)dx = Sigma  Fi*P(xi)
 |
          -1/          i=1
   Formalparameter:
           E: int anzgp:    Anzahl der Gauspunkte im Dreieck
   E: long xx:      Nr des Gauspunktes
   E: long coord:   coordr
   return : double coord des Gauspunktes
   Ergebnis:

   Aenderungen/Korrekturen:
   07/2003     mb        Erste Version

 **************************************************************************/
double MXPGaussPktTri(int anzgp, long xx, long coordr)
{
	double a[3] = {0.66666666666666, 0.16666666666666, 0.16666666666666};
	double b[3] = {0.16666666666666, 0.66666666666666, 0.16666666666666};
	double c[3] = {0.16666666666666, 0.16666666666666, 0.66666666666666};
	if (anzgp == 1)
		return 0.3333333333333333;
	else if (anzgp == 3)
	{
		switch (xx)
		{
			case 0:
				return a[coordr];
			case 1:
				return b[coordr];
			case 2:
				return c[coordr];
		}
	}
	else if (anzgp == 4)
	{
		switch (xx)
		{
			case 0:
				switch (coordr)
				{
					case 0:
						return 0.333333333333333;
					case 1:
						return 0.333333333333333;
				}
				break;

			case 1:
				switch (coordr)
				{
					case 0:
						return 0.600000000000000;
					case 1:
						return 0.200000000000000;
				}
				break;

			case 2:
				switch (coordr)
				{
					case 0:
						return 0.200000000000000;
					case 1:
						return 0.600000000000000;
				}
				break;

			case 3:
				switch (coordr)
				{
					case 0:
						return 0.200000000000000;
					case 1:
						return 0.200000000000000;
				}
				break;
		}
	} /*else if*/ /* switch grd */
	return 0.0;
}

/***************************************************************************
   ROCKFLOW - Funktion: MXPGaussFktTri
   Aufgabe:
   X Punkt Gauss-Integration Dreieck
           bestimmtes Integral
           1/           X
 | P(x)dx = Sigma  Fi*P(xi)
 |
          -1/          i=1
   Formalparameter:
           E: int anzgp:     Anzahl der Gauspunkte im Dreieck
   E: long pkt:      Nr des Gauspunktes
   return : Wichtung des Gauspunktes
   Ergebnis:  double: Wichtung des Gauspunktes

   Aenderungen/Korrekturen:
   07/2003     mb        Erste Version

 **************************************************************************/
double MXPGaussFktTri(int anzgp, long pkt)
{
	if (anzgp == 1)
		return 0.5;
	else if (anzgp == 3)
		return 0.166666666666666;
	else if (anzgp == 4)
	{
		switch (pkt)
		{
			case 0:
				return -0.281250000000000;
				break;
			case 1:
				return 0.260416666666667;
				break;
			case 2:
				return 0.260416666666667;
				break;
			case 3:
				return 0.260416666666667;
				break;
		}
	}
	return 0.0;
}

/**************************************************************************
   ROCKFLOW - Function: Get_Nt_x_Nt

   Task:
   Gets the linear component of the capacitance matrix for triangular prisms:
   integral Nt x Nt

   Parameter: (I: Input; R: Return; X: Both)
           I: void

   Return:
  *Ct

   Programming:
   07/2003   MB   First Version

**************************************************************************/
int Get_Nt_x_Nt(double* Ct)
{
	double Ct1;
	double Ct2;

	Ct1 = 2.0 / 3.0;
	Ct2 = 1.0 / 3.0;

	Ct[0] = Ct[21] = Ct1;
	Ct[1] = Ct[22] = Ct1;
	Ct[2] = Ct[23] = Ct1;
	Ct[6] = Ct[27] = Ct1;
	Ct[7] = Ct[28] = Ct1;
	Ct[8] = Ct[29] = Ct1;
	Ct[12] = Ct[33] = Ct1;
	Ct[13] = Ct[34] = Ct1;
	Ct[14] = Ct[35] = Ct1;

	Ct[3] = Ct[18] = Ct2;
	Ct[4] = Ct[19] = Ct2;
	Ct[5] = Ct[20] = Ct2;
	Ct[9] = Ct[24] = Ct2;
	Ct[10] = Ct[25] = Ct2;
	Ct[11] = Ct[26] = Ct2;
	Ct[15] = Ct[30] = Ct2;
	Ct[16] = Ct[31] = Ct2;
	Ct[17] = Ct[32] = Ct2;

	return 1;
}

/**************************************************************************
   ROCKFLOW - Function: Get_Nt_x_gradNt

   Task:
   Gets the linear component for triangular prisms:
   integral Nt x gradNt

   Parameter: (I: Input; R: Return; X: Both)
           I: void

   Return:
  *Ct

   Programming:
   07/2003   MB   First Version

**************************************************************************/
int Get_Nt_x_gradNt(double* Nt_x_gradNt)
{
	double Ct1;
	double Ct2;

	Ct1 = 1.0 / 2.0;
	Ct2 = -1.0 / 2.0;

	Nt_x_gradNt[0] = Nt_x_gradNt[18] = Ct1;
	Nt_x_gradNt[1] = Nt_x_gradNt[19] = Ct1;
	Nt_x_gradNt[2] = Nt_x_gradNt[20] = Ct1;
	Nt_x_gradNt[6] = Nt_x_gradNt[24] = Ct1;
	Nt_x_gradNt[7] = Nt_x_gradNt[25] = Ct1;
	Nt_x_gradNt[8] = Nt_x_gradNt[26] = Ct1;
	Nt_x_gradNt[12] = Nt_x_gradNt[30] = Ct1;
	Nt_x_gradNt[13] = Nt_x_gradNt[31] = Ct1;
	Nt_x_gradNt[14] = Nt_x_gradNt[32] = Ct1;

	Nt_x_gradNt[3] = Nt_x_gradNt[21] = Ct2;
	Nt_x_gradNt[4] = Nt_x_gradNt[22] = Ct2;
	Nt_x_gradNt[5] = Nt_x_gradNt[23] = Ct2;
	Nt_x_gradNt[9] = Nt_x_gradNt[27] = Ct2;
	Nt_x_gradNt[10] = Nt_x_gradNt[28] = Ct2;
	Nt_x_gradNt[11] = Nt_x_gradNt[29] = Ct2;
	Nt_x_gradNt[15] = Nt_x_gradNt[33] = Ct2;
	Nt_x_gradNt[16] = Nt_x_gradNt[34] = Ct2;
	Nt_x_gradNt[17] = Nt_x_gradNt[35] = Ct2;

	return 1;
}

/**************************************************************************
   ROCKFLOW - Function: Get_gradNt_x_Nt

   Task:
   Gets the linear component for triangular prisms:
   integral Nt x gradNt

   Parameter: (I: Input; R: Return; X: Both)
           I: void

   Return:
  *Ct

   Programming:
   07/2003   MB   First Version

**************************************************************************/
int Get_gradNt_x_Nt(double* gradNt_x_Nt)
{
	double Ct1;
	double Ct2;

	Ct1 = 1.0 / 2.0;
	Ct2 = -1.0 / 2.0;

	gradNt_x_Nt[0] = gradNt_x_Nt[3] = Ct1;
	gradNt_x_Nt[1] = gradNt_x_Nt[4] = Ct1;
	gradNt_x_Nt[2] = gradNt_x_Nt[5] = Ct1;
	gradNt_x_Nt[6] = gradNt_x_Nt[9] = Ct1;
	gradNt_x_Nt[7] = gradNt_x_Nt[10] = Ct1;
	gradNt_x_Nt[8] = gradNt_x_Nt[11] = Ct1;
	gradNt_x_Nt[12] = gradNt_x_Nt[15] = Ct1;
	gradNt_x_Nt[13] = gradNt_x_Nt[16] = Ct1;
	gradNt_x_Nt[14] = gradNt_x_Nt[17] = Ct1;

	gradNt_x_Nt[18] = gradNt_x_Nt[21] = Ct2;
	gradNt_x_Nt[19] = gradNt_x_Nt[22] = Ct2;
	gradNt_x_Nt[20] = gradNt_x_Nt[23] = Ct2;
	gradNt_x_Nt[24] = gradNt_x_Nt[27] = Ct2;
	gradNt_x_Nt[25] = gradNt_x_Nt[28] = Ct2;
	gradNt_x_Nt[26] = gradNt_x_Nt[29] = Ct2;
	gradNt_x_Nt[30] = gradNt_x_Nt[33] = Ct2;
	gradNt_x_Nt[31] = gradNt_x_Nt[34] = Ct2;
	gradNt_x_Nt[32] = gradNt_x_Nt[35] = Ct2;

	return 1;
}

/**************************************************************************
   ROCKFLOW - Function: Get_gradNt_x_gradNt

   Task:
   Gets the linear component for triangular prisms:
   integral gradNt x gradNt

   Parameter: (I: Input; R: Return; X: Both)
           I: void

   Return:
  *GradNGradN

   Programming:
   07/2003   MB   First Version

**************************************************************************/
int Get_gradNt_x_gradNt(double* GradNGradN)
{
	double Ct1;
	double Ct2;

	Ct1 = 1.0 / 2.0;
	Ct2 = -1.0 / 2.0;

	GradNGradN[0] = GradNGradN[21] = Ct1;
	GradNGradN[1] = GradNGradN[22] = Ct1;
	GradNGradN[2] = GradNGradN[23] = Ct1;
	GradNGradN[6] = GradNGradN[27] = Ct1;
	GradNGradN[7] = GradNGradN[28] = Ct1;
	GradNGradN[8] = GradNGradN[29] = Ct1;
	GradNGradN[12] = GradNGradN[33] = Ct1;
	GradNGradN[13] = GradNGradN[34] = Ct1;
	GradNGradN[14] = GradNGradN[35] = Ct1;

	GradNGradN[3] = GradNGradN[18] = Ct2;
	GradNGradN[4] = GradNGradN[19] = Ct2;
	GradNGradN[5] = GradNGradN[20] = Ct2;
	GradNGradN[9] = GradNGradN[24] = Ct2;
	GradNGradN[10] = GradNGradN[25] = Ct2;
	GradNGradN[11] = GradNGradN[26] = Ct2;
	GradNGradN[15] = GradNGradN[30] = Ct2;
	GradNGradN[16] = GradNGradN[31] = Ct2;
	GradNGradN[17] = GradNGradN[32] = Ct2;

	return 1;
}

/***************************************************************************
   ROCKFLOW - Funktion: MOmegaPrism
   Aufgabe:
           Berechnet Ansatzfunktion (r,s,t) eines Dreieckprismas.

   Formalparameter:
           Z: *vf - 1x6 Feld
           E: r,s,t (s.o.)
   Ergebnis:
           double *vf
   Aenderungen/Korrekturen:
   07/2003   MB   First Version
 **************************************************************************/

int MOmegaPrism(double* vf, double r, double s, double t)
{
	long i;
	int ok = 0;
#ifdef ERROR_CONTROL
	if (vf == NULL)
	{
		ok = 0;
		return ok;
	}
#endif
	vf[0] = (1.0 - r - s) * (1.0 + t);
	vf[1] = r * (1.0 + t);
	vf[2] = s * (1.0 + t);
	vf[3] = (1.0 - r - s) * (1.0 - t);
	vf[4] = r * (1.0 - t);
	vf[5] = s * (1.0 - t);
	for (i = 0; i < 6; i++)
		vf[i] *= 0.5;
	return ok = 1;
}

/***************************************************************************
   ROCKFLOW - Funktion: MGradOmegaPrism
   Aufgabe:
           Berechnet Gradient der Ansatzfunktionen (r,s,t).

 |  0   1   2   3   4   5   |
 |                          |
      grad (omega)= |  6   7   8   9  10  11   |
 |                          |
 | 12  13  14  15  16  17   |

   Zahlen in der Matrix stehen fuer Positionen im Feld vf.
   Gleichungen s.u..

   Formalparameter:
   Z: *vf - 3x6 Feld
   E: r,s (s.o.)
   Ergebnis:
   3x8 Matrix
   Aenderungen/Korrekturen:
   06/2003     mb        Erste Version

 **************************************************************************/

int MGradOmegaPrism(double r, double s, double t, double* vf)
{
	long i;
	int ok = 0;
#ifdef ERROR_CONTROL
	if (vf == NULL)
	{
		ok = 0;
		return ok;
	}
#endif
	/*grad r N */
	vf[0] = -(1.0 + t);
	vf[1] = (1.0 + t);
	vf[2] = 0.0;
	vf[3] = -(1.0 - t);
	vf[4] = (1.0 - t);
	vf[5] = 0.0;
	/*grad s N */
	vf[6] = -(1.0 + t);
	vf[7] = 0.0;
	vf[8] = (1.0 + t);
	vf[9] = -(1.0 - t);
	vf[10] = 0.0;
	vf[11] = (1.0 - t);
	/*grad t N */
	vf[12] = (1.0 - r - s);
	vf[13] = r;
	vf[14] = s;
	vf[15] = -(1.0 - r - s);
	vf[16] = -r;
	vf[17] = -s;

	for (i = 0; i < 18; i++)
		vf[i] *= 0.5;
	return ok = 1;
}
