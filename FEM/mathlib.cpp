/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************

   ROCKFLOW - Modul: mathlib.c

   Aufgabe:
   ROCKFLOW-Schnittstelle fuer alle mathematischen Operationen, die nicht
   Standard (nicht in math.h enthalten) sind.
   mathlib.h benutzt externe Bibliotheken, die bei Bedarf ausgetauscht
   werden koennen, ohne die ROCKFLOW-Schnittstelle aendern zu muessen.

   Ueberblick ueber die Funktionen:
   Konstanten
   Mdrittel              - 1/3
   MKleinsteZahl         - kleinste Zahl, die groesser 0.0 fuer den
   Rechner ist
   MCalcKleinsteZahl     - "Berechnet" die kleinste Zahl
   M2GPIKoo1             - 2 Punkt Gauss-Integration Koordinate 1
   M2GPIKoo2             -  Koordinate 2
   M2GPIFk1              -  Faktor 1
   M2GPIFk2              -  Faktor 2
   M3GPIKoo1             - 3 Punkt Gauss-Integration Koordinate 1
   M3GPIKoo2             -  Koordinate 2
   M3GPIKoo3             -  Koordinate 3
   M3GPIFk1              -  Faktor 1
   M3GPIFk2              -  Faktor 2
   M3GPIFk3              -  Faktor 3
   Mathematische Funktionen
   MGleichDouble         - Vergleicht zwei double-Zahlen unter
   Beruecksichtigung einer Fehlertoleranz
   MMin                  - ermittelt die kleinere zweier Double-Werte
   MMax                  - ermittelt die groessere zweier Double-Werte
   MRange                - begrenzt eine Zahl auf ein Intervall
   MOmega1D              - Ansatzfunktionen
   MOmega2D              - Ansatzfunktionen
   MOmega2DTriangle      - Ansatzfunktionen für linear Dreiecke
   MOmega3D                (s.o.)
   MPhi2D                - Testfunktionen
   MPhi3D
   MPhi2D_SUPG           - Testfunktionen
   MPhi3D_SUPG             (s.o.)
   MGradOmega2D          - Gradient der Ansatzfunktionen
   MGradOmega3D            (s.o.)
   MGradPhi2D            - Gradient der Testfunktionen
   MGradPhi3D              (s.o.)
   MOmega2D_9N           - Ansatzfunktionen (zweidim. 9-Knoten-Elemente)
   MPhi2D_9N             - Testfunktionen (zweidim. 9-Knoten-Elemente)
   MGradOmega2D_9N       - Gradient der Ansatzfunktionen (zweidim. 9-Knoten-Elemente)
   MGradPhi2D_9N         - Gradient der Testfunktionen (zweidim. 9-Knoten-Elemente)
   M2PGaussInt           - 2 Punkt Gauss-Integration
   M3PGaussInt           - 3 Punkt Gauss-Integration
   MXPGaussPkt           - Punkte fuer die X Punkt Gauss-Integration
   MXPGaussFkt           - Faktoren fuer die X Punkt Gauss-Integration
   Funktionen fuer/mit Vektoren und Matrizen
   MNulleVec             - Setze angegebenen Vektor = 0.0
   MNulleMat             - Setze angegebene Matrix = 0.0
   MBtrgVec              - Betrag von Vektor
   MB
   MNormiere             - Normiere Vektor
   M2Determinante        - Determinante einer 2x2Matrix
   M3Determinante        - Determinante einer 3x3Matrix
   M4Determinante        - Determinante einer 4x4Matrix
   Mxg2Determinante      - Determinante einer beliebigen Matrix
   mit goesseren Ausmassen als 2x2
   (Warnung: funktioniert leider nicht richtig.)
   MTranspoVec           - Transponieren beliebiger Vektoren
   MTranspoMat           - Transponieren beliebiger Matrizen
   M2Invertiere          - Invertiert 2x2 Matrizen
   M3Invertiere          - Invertiert 3x3 Matrizen
   MAddVektoren          - Addition zweier beliebier Vektoren
   MAddSkalVektoren      - Vektoren mit Skalar multiplizieren und dann addieren
   MAddMatrizen          - Addition zweier beliebier Matrizen
   MMultVecSkalar        - Multiplikation Vektor mit Skalarwert
   MMultMatSkalar        - Multiplikation Matrix mit Skalarwert
   MSkalarprodukt        - Skalarprodukt zweier beliebiger Vektoren
   M3KreuzProdukt        - Kreuzprodukt von 3D-Vektoren
   MMultVecVec           - Multiplikation Vektor mit Vektor (-> Matrix)
   MMultVecMat           - Multiplikation Vektor mit Matrix (-> Vektor)
   MMultMatVec           - Multiplikation Matrix mit Vektor (-> Vektor)
   MMultMatMat           - Multiplikation Matrix mit Matrix (-> Matrix)
   MVekDist              - Abstand zwischen zwei Vektoren
   Auf Matrizen und Vektoren aufbauende Funktionen
   MKTF2Dr2D             - Koordinaten-Transformation 2D -> 2D (R=rotiert)
   MKTFMat2Dr2D          - Koordinaten-Transf. 2D -> 2D ; berechnet T-Matrix
   MKTF3Dr2D             - Koordinaten-Transformation 3D -> 2D (R=rotiert)
   MKTFMat3Dr2D          - Koordinaten-Transf. 3D -> 2D ; berechnet T-Matrix
   MKTF2Dt2D             - Koordinaten-Transformation 2D -> 2D (T=transvers.)
   Prueffunktion
   MBistDuDiagMat        - Prueft auf Diagonalitaet einer Matrix
   (nicht gerade schnell :-( )
   Bearbeitungsfunktionen
   MMachVec              - Erzeugt einen Vektor
   MNullVec              - Nullt einen Vektor
   MLoeschVec            - Zerstoert einen Vektor
   MKopierVec            - Kopiert einen Vektor auf einen Anderen
   MMachMat              - Erzeugt eine Matrix
   MNullMat              - Nullt eine Matrix
   MLoeschMat            - Zerstoert eine Matrix
   MKopierMat            - Kopiert eine Matrix auf eine Andere

   Geometrie-Funktionen
   MCalcDistancePointToPoint  -  Abstand Punkt - Punkt im R3
   MCalcDistancePointToLine   -  Abstand Punkt - Linie im R3
   MCalcDistancePointToPlane  -  Abstand Punkt - Ebene im R3
   MCalcProjectionOfPointOnLine  - Projektionspunkt eines Punktes auf eine Linie (Lotpunkt)
   MCalcProjectionOfPointOnPlane - Projektionspunkt eines Punktes auf einer Flaeche (Lotpunkt)

   Sortierfunktionen
   MQSort_LongDouble     - Sortiert Datensaetze aus long und double nach dem groessten double

   Aenderungen/Korrekturen:
   08/1994     Hans Herrmann        Erste Version
   09/1994     hh                   Aenderungen
   12/1994     MSR                  Display-Routinen von mathlib
   nach display portiert
   05/1995     hh                   Erweiterung um grad()
   07/1995     hh                   Erweiterung um MKTFxxxxx(),MInvertiere()
   MBtrgVec(),MGPIxxx
   08/1995     cb                   phi eingebaut
   11/1995     msr                  Vektornormen und -Gleichungen eingebaut
   11/1995     cb                   Erweiterung um phi(), grad()
   08/1996     cb                   4 Punkt Gauss Integration
   01.07.1997  R.Kaiser             MPhi2D und MPhi3D  (ohne SUPG)
   MPhi2D_SUPG und MPhi3D_SUPG (mit SUPG)
   09/97       dh                   MVekDist
   02/2000     RK                   MPhi2D_9N, MOmega2D_N
   MGradPhi2D_9N, MGradOmega2D_9N
   02/2000   C.Thorenz              MMin, MMax,    Geometrie-Funktionen
   MCalcDistancePointToPoint  -  Abstand Punkt - Punkt im R3
   MCalcDistancePointToLine   -  Abstand Punkt - Linie im R3
   MCalcDistancePointToPlane  -  Abstand Punkt - Ebene im R3
   MCalcProjectionOfPointOnLine  - Projektionspunkt eines Punktes auf eine Linie (Lotpunkt)
   MCalcProjectionOfPointOnPlane - Projektionspunkt eines Punktes auf einer Flaeche (Lotpunkt)
   06/2001   M. Kohlmeier           Sortierfunktion
   08/2001     MK                   Alle M2Invertiere in M2InvertiereUndTransponiere umbenannt
   M2Invertiere korrigiert
   04/2002     OK                   M4Determinante
   07/2003     WW                   Triangle shape funtions
   06/2004     WW                   Generalized shape functions
 ***************************************************************************/
// There is a name conflict between stdio.h and the MPI C++ binding
// with respect to the names SEEK_SET, SEEK_CUR, and SEEK_END.  MPI
// wants these in the MPI namespace, but stdio.h will #define these
// to integer values.  #undef'ing these can cause obscure problems
// with other include files (such as iostream), so we instead use
// #error to indicate a fatal error.  Users can either #undef
// the names before including mpi.h or include mpi.h *before* stdio.h
// or iostream.
#if defined(USE_MPI) || defined(USE_MPI_PARPROC) || defined(USE_MPI_REGSOIL)
#include <mpi.h>
#endif

/* Preprozessor-Definitionen */

#include <cfloat>
#include <cstdlib>

//#include <stdio.h>
#include "femlib.h" //CMCD 03 2004
#include "makros.h"
#include "memory.h"
#include "display.h"
#include "mathlib.h"
#include "geo_mathlib.h" // for M3KreuzProdukt

// WW----------------------
#include "par_ddc.h"
// WW----------------------

#include "prototyp.h"

double pai = 4.0 * atan(1.0);
VoidFuncDXCDX ShapeFunction;
VoidFuncDXCDX ShapeFunctionHQ;
VoidFuncDXCDX GradShapeFunction = NULL;

/*##########################################################################
   Mathematische Funktionen
 ######################################################################## */

/***************************************************************************
   GEO MathLib - Funktion: MBtrgVec
   Aufgabe:
           Berechnet Betrag von Vektor
   Formalparameter:
           E: *vec
           E: n
   Ergebnis:
           Betrag
   Aenderungen/Korrekturen:
   07/1995     hh        Erste Version
   11/1999     C.Thorenz Register-Variablen

 **************************************************************************/

double MBtrgVec(double* vec, long n)
{
	register long i;
	register double zwo = 0.0;
	for (i = 0; i < n; i++)
		zwo += vec[i] * vec[i];
	return sqrt(zwo);
}

#ifdef obsolete // 05.03.2010 WW
/***************************************************************************
   ROCKFLOW - Funktion: MGleichDouble
   Aufgabe:
           Vergleicht zwei double-Zahlen unter Beruecksichtigung
           einer Fehlertoleranz
   Formalparameter:
           E: zahl1, zahl2
           E: tol - Fehlertoleranz (positiv)
   Ergebnis:
            0 - ungleich
            1 - gleich
   Programmaenderungen:
   07/1994     hh        Erste Version

 **************************************************************************/

int MGleichDouble(double zahl1, double zahl2, double tol)
{
	int retval;
	if (fabs(zahl1 - zahl2) <= tol)
		retval = 1;
	else
		retval = 0;
	return retval;
}
#endif //#ifdef obsolete  //05.03.2010 WW

////////////////////////////////////////////////////////////
#ifdef obsolete // 05.03.2010 WW
/***************************************************************************
   ROCKFLOW - Funktion: MOmega1D
   Aufgabe:
           Berechnet Ansatzfunktion (r).
                        /       \
                     1  | (1-r) |
            Omega = --- |       |
                     2  | (1+r) |
 \       /
   Formalparameter:
           Z: *vf - 1x2 Feld
   E: r  (s.o.)
   Ergebnis:
   Vektor
   Aenderungen/Korrekturen:
   08/1995     cb        Erste Version

 **************************************************************************/

int MOmega1D(double* vf, double r)
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
	vf[0] = (1.0 + r);
	vf[1] = (1.0 - r);
	for (i = 0; i < 2; i++)
		vf[i] *= 0.5;
	return ok = 1;
}

/***************************************************************************
   ROCKFLOW - Funktion: MOmega2D
   Aufgabe:
           Berechnet Ansatzfunktion (r,s).
                        / (1+r)(1+s) \
 |            |
                     1  | (1-r)(1+s) |
            Omega = --- |            |
                     4  | (1-r)(1-s) |
 |            |
 \ (1+r)(1-s) /
   Formalparameter:
   Z: *vf - 1x4 Feld
   E: r,s (s.o.)
   Ergebnis:
   Vektor
   Aenderungen/Korrekturen:
   08/1995     cb        Erste Version

 **************************************************************************/

int MOmega2D(double* vf, double r, double s)
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
	vf[0] = (1.0 + r) * (1.0 + s);
	vf[1] = (1.0 - r) * (1.0 + s);
	vf[2] = (1.0 - r) * (1.0 - s);
	vf[3] = (1.0 + r) * (1.0 - s);
	for (i = 0; i < 4; i++)
		vf[i] *= 0.25;
	return ok = 1;
}

/***************************************************************************
   ROCKFLOW - Funktion: MOmega3D
   Aufgabe:
           Berechnet Ansatzfunktion (r,s).
                        / (1+r)(1+s)(1+t) \
 | (1-r)(1+s)(1+t) |
                     1  | (1-r)(1-s)(1+t) |
            Omega = --- | (1+r)(1-s)(1+t) |
                     8  | (1+r)(1+s)(1-t) |
 | (1-r)(1+s)(1-t) |
 | (1-r)(1-s)(1-t) |
 \ (1+r)(1-s)(1-t) /
   Formalparameter:
   Z: *vf - 1x8 Feld
   E: r,s,t (s.o.)
   Ergebnis:
   Vektor
   Aenderungen/Korrekturen:
   08/1995     cb        Erste Version

 **************************************************************************/

int MOmega3D(double* vf, double r, double s, double t)
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
	vf[0] = (1.0 + r) * (1.0 + s) * (1.0 + t);
	vf[1] = (1.0 - r) * (1.0 + s) * (1.0 + t);
	vf[2] = (1.0 - r) * (1.0 - s) * (1.0 + t);
	vf[3] = (1.0 + r) * (1.0 - s) * (1.0 + t);
	vf[4] = (1.0 + r) * (1.0 + s) * (1.0 - t);
	vf[5] = (1.0 - r) * (1.0 + s) * (1.0 - t);
	vf[6] = (1.0 - r) * (1.0 - s) * (1.0 - t);
	vf[7] = (1.0 + r) * (1.0 - s) * (1.0 - t);
	for (i = 0; i < 8; i++)
		vf[i] *= 0.125;
	return ok = 1;
}

/***************************************************************************
   ROCKFLOW - Funktion: MOmega2DTriangle
   Aufgabe:
           Berechnet Ansatzfunktion (r,s).
                        / (1+r)(1+s) \
 |            |
                     1  | (1-r)(1+s) |
            Omega = --- |            |
                     4  | (1-r)(1-s) |
 |            |
 \ (1+r)(1-s) /
   Formalparameter:
   Z: *vf - 1x4 Feld
   E: r,s (s.o.)
   Ergebnis:
   Vektor
   Aenderungen/Korrekturen:
   08/1995     cb        Erste Version

 **************************************************************************/
//#ifdef obsolete //WW. 06.11.2008
int MOmega2DTriangle(double* vf, double xx, double yy, long number)
{
	int i, nn = 3, ok = 0;
	double x[3], y[3];
	double area;
	long* element_nodes;

#ifdef ERROR_CONTROL
	if (vf == NULL)
	{
		ok = 0;
		return ok;
	}
#endif
	area = fabs(ElGetElementVolume(number)); // should be triangle area !
	/*  Calc2DElementCoordinatesTriangle(number,x,y); */
	element_nodes = ElGetElementNodes(number);
	for (i = 0; i < nn; i++)
	{
		x[i] = GetNodeX(element_nodes[i]);
		y[i] = GetNodeY(element_nodes[i]);
	}
	/* CMCD March 2004 The area of the triangle needs to be calculated for the 2D projected triangle
	   not the three D as above
	    area = (x[1]*y[2]-x[2]*y[1]+x[0]*y[1]-x[1]*y[0]+x[2]*y[0]-x[0]*y[2])/2; // CMCD March 2004
	   xx=0.;
	   yy=0.;
	   what is this ?
	 */
	vf[0] = (x[1] * y[2] - x[2] * y[1]) + (y[1] - y[2]) * xx + (x[2] - x[1]) * yy;
	vf[1] = (x[2] * y[0] - x[0] * y[2]) + (y[2] - y[0]) * xx + (x[0] - x[2]) * yy;
	vf[2] = (x[0] * y[1] - x[1] * y[0]) + (y[0] - y[1]) * xx + (x[1] - x[0]) * yy;
	vf[0] /= (2. * area);
	vf[1] /= (2. * area);
	vf[2] /= (2. * area);

	return ok = 1;
}

/**************************************************************************
   GeoLib-Method: MOmega3DTetrahedron
   Task: Interpolation function for tetrahedra
   Programing:
   08/2003 OK Implementation
**************************************************************************/
int MOmega3DTetrahedron(double* vf, double r, double s, double t, long number)
{
	int ok = 0;
	long* element_nodes;
	int i;
	double a1, a2, a3, a4;
	double b1, b2, b3, b4;
	double c1, c2, c3, c4;
	double d1, d2, d3, d4;
	double x[4], y[4], z[4];
	double N1, N2, N3, N4;
	double volume;
	double mat3x3[9];

	volume = fabs(ElGetElementVolume(number));
	element_nodes = ElGetElementNodes(number);
	for (i = 0; i < 4; i++)
	{
		x[i] = GetNodeX(element_nodes[i]);
		y[i] = GetNodeY(element_nodes[i]);
		z[i] = GetNodeZ(element_nodes[i]);
	}

	mat3x3[0] = x[1];
	mat3x3[1] = y[1];
	mat3x3[2] = z[1];
	mat3x3[3] = x[2];
	mat3x3[4] = y[2];
	mat3x3[5] = z[2];
	mat3x3[6] = x[3];
	mat3x3[7] = y[3];
	mat3x3[8] = z[3];
	a1 = -1.0 * M3Determinante(mat3x3);
	mat3x3[0] = x[2];
	mat3x3[1] = y[2];
	mat3x3[2] = z[2];
	mat3x3[3] = x[3];
	mat3x3[4] = y[3];
	mat3x3[5] = z[3];
	mat3x3[6] = x[0];
	mat3x3[7] = y[0];
	mat3x3[8] = z[0];
	a2 = 1.0 * M3Determinante(mat3x3);
	mat3x3[0] = x[3];
	mat3x3[1] = y[3];
	mat3x3[2] = z[3];
	mat3x3[3] = x[0];
	mat3x3[4] = y[0];
	mat3x3[5] = z[0];
	mat3x3[6] = x[1];
	mat3x3[7] = y[1];
	mat3x3[8] = z[1];
	a3 = -1.0 * M3Determinante(mat3x3);
	mat3x3[0] = x[0];
	mat3x3[1] = y[0];
	mat3x3[2] = z[0];
	mat3x3[3] = x[1];
	mat3x3[4] = y[1];
	mat3x3[5] = z[1];
	mat3x3[6] = x[2];
	mat3x3[7] = y[2];
	mat3x3[8] = z[2];
	a4 = 1.0 * M3Determinante(mat3x3);

	mat3x3[0] = 1.0;
	mat3x3[1] = y[1];
	mat3x3[2] = z[1];
	mat3x3[3] = 1.0;
	mat3x3[4] = y[2];
	mat3x3[5] = z[2];
	mat3x3[6] = 1.0;
	mat3x3[7] = y[3];
	mat3x3[8] = z[3];
	b1 = -1.0 * M3Determinante(mat3x3);
	mat3x3[0] = 1.0;
	mat3x3[1] = y[2];
	mat3x3[2] = z[2];
	mat3x3[3] = 1.0;
	mat3x3[4] = y[3];
	mat3x3[5] = z[3];
	mat3x3[6] = 1.0;
	mat3x3[7] = y[0];
	mat3x3[8] = z[0];
	b2 = 1.0 * M3Determinante(mat3x3);
	mat3x3[0] = 1.0;
	mat3x3[1] = y[3];
	mat3x3[2] = z[3];
	mat3x3[3] = 1.0;
	mat3x3[4] = y[0];
	mat3x3[5] = z[0];
	mat3x3[6] = 1.0;
	mat3x3[7] = y[1];
	mat3x3[8] = z[1];
	b3 = -1.0 * M3Determinante(mat3x3);
	mat3x3[0] = 1.0;
	mat3x3[1] = y[0];
	mat3x3[2] = z[0];
	mat3x3[3] = 1.0;
	mat3x3[4] = y[1];
	mat3x3[5] = z[1];
	mat3x3[6] = 1.0;
	mat3x3[7] = y[2];
	mat3x3[8] = z[2];
	b4 = 1.0 * M3Determinante(mat3x3);

	mat3x3[0] = x[1];
	mat3x3[1] = 1.0;
	mat3x3[2] = z[1];
	mat3x3[3] = x[2];
	mat3x3[4] = 1.0;
	mat3x3[5] = z[2];
	mat3x3[6] = x[3];
	mat3x3[7] = 1.0;
	mat3x3[8] = z[3];
	c1 = -1.0 * M3Determinante(mat3x3);
	mat3x3[0] = x[2];
	mat3x3[1] = 1.0;
	mat3x3[2] = z[2];
	mat3x3[3] = x[3];
	mat3x3[4] = 1.0;
	mat3x3[5] = z[3];
	mat3x3[6] = x[0];
	mat3x3[7] = 1.0;
	mat3x3[8] = z[0];
	c2 = 1.0 * M3Determinante(mat3x3);
	mat3x3[0] = x[3];
	mat3x3[1] = 1.0;
	mat3x3[2] = z[3];
	mat3x3[3] = x[0];
	mat3x3[4] = 1.0;
	mat3x3[5] = z[0];
	mat3x3[6] = x[1];
	mat3x3[7] = 1.0;
	mat3x3[8] = z[1];
	c3 = -1.0 * M3Determinante(mat3x3);
	mat3x3[0] = x[0];
	mat3x3[1] = 1.0;
	mat3x3[2] = z[0];
	mat3x3[3] = x[1];
	mat3x3[4] = 1.0;
	mat3x3[5] = z[1];
	mat3x3[6] = x[2];
	mat3x3[7] = 1.0;
	mat3x3[8] = z[2];
	c4 = 1.0 * M3Determinante(mat3x3);

	mat3x3[0] = x[1];
	mat3x3[1] = y[1];
	mat3x3[2] = 1.0;
	mat3x3[3] = x[2];
	mat3x3[4] = y[2];
	mat3x3[5] = 1.0;
	mat3x3[6] = x[3];
	mat3x3[7] = y[3];
	mat3x3[8] = 1.0;
	d1 = -1.0 * M3Determinante(mat3x3);
	mat3x3[0] = x[2];
	mat3x3[1] = y[2];
	mat3x3[2] = 1.0;
	mat3x3[3] = x[3];
	mat3x3[4] = y[3];
	mat3x3[5] = 1.0;
	mat3x3[6] = x[0];
	mat3x3[7] = y[0];
	mat3x3[8] = 1.0;
	d2 = 1.0 * M3Determinante(mat3x3);
	mat3x3[0] = x[3];
	mat3x3[1] = y[3];
	mat3x3[2] = 1.0;
	mat3x3[3] = x[0];
	mat3x3[4] = y[0];
	mat3x3[5] = 1.0;
	mat3x3[6] = x[1];
	mat3x3[7] = y[1];
	mat3x3[8] = 1.0;
	d3 = -1.0 * M3Determinante(mat3x3);
	mat3x3[0] = x[0];
	mat3x3[1] = y[0];
	mat3x3[2] = 1.0;
	mat3x3[3] = x[1];
	mat3x3[4] = y[1];
	mat3x3[5] = 1.0;
	mat3x3[6] = x[2];
	mat3x3[7] = y[2];
	mat3x3[8] = 1.0;
	d4 = 1.0 * M3Determinante(mat3x3);

	// Element Shape Functions
	N1 = ((a1 * 1) + (b1 * r) + (c1 * s) + (d1 * t)) / (6 * volume);
	N2 = ((a2 * 1) + (b2 * r) + (c2 * s) + (d2 * t)) / (6 * volume);
	N3 = ((a3 * 1) + (b3 * r) + (c3 * s) + (d3 * t)) / (6 * volume);
	N4 = ((a4 * 1) + (b4 * r) + (c4 * s) + (d4 * t)) / (6 * volume);

	vf[0] = N1;
	vf[1] = N2;
	vf[2] = N3;
	vf[3] = N4;

	return ok = 1;
}

//#endif //#ifndef obsolete //WW. 06.11.2008

/***************************************************************************
   ROCKFLOW - Funktion: MPhi3D
   Aufgabe:
           Berechnet Testfunktion (r,s,t) ohne Upwind-Parameter alpha.
           Phi = Omega + alpha * grad Omega
   Formalparameter:
           Z: *vf - 1x8 Feld
           E: r,s,t (s.o.)
   Ergebnis:
           Vektor
   Aenderungen/Korrekturen:
   01.07.1997  R.Kaiser  erste Version

 **************************************************************************/

int MPhi3D(double* vf, double r, double s, double t)
{
	int i;
	int ok = 0;
#ifdef ERROR_CONTROL
	if (vf == NULL)
	{
		ok = 0;
		return ok;
	}
#endif
	vf[0] = (1.0 + r) * (1.0 + s) * (1.0 + t);
	vf[1] = (1.0 - r) * (1.0 + s) * (1.0 + t);
	vf[2] = (1.0 - r) * (1.0 - s) * (1.0 + t);
	vf[3] = (1.0 + r) * (1.0 - s) * (1.0 + t);
	vf[4] = (1.0 + r) * (1.0 + s) * (1.0 - t);
	vf[5] = (1.0 - r) * (1.0 + s) * (1.0 - t);
	vf[6] = (1.0 - r) * (1.0 - s) * (1.0 - t);
	vf[7] = (1.0 + r) * (1.0 - s) * (1.0 - t);
	for (i = 0; i < 8; i++)
		vf[i] *= 0.125;
	return ok = 1;
}

/***************************************************************************
   ROCKFLOW - Funktion: MPhi2D_SUPG
   Aufgabe:
           Berechnet Testfunktion (r,s) incl. Upwind-Parameter alpha.
           Phi = Omega + alpha * grad Omega
   Formalparameter:
           Z: *vf - 1x4 Feld
           E: r,s (s.o.)
   Ergebnis:
           Vektor
   Aenderungen/Korrekturen:
   11/1995     cb        Erste Version
   01.07.1997  R.Kaiser  Umbenannt zu MPhi2D_SUPG

 **************************************************************************/

int MPhi2D_SUPG(double* vf, double r, double s, double* alpha)
{
	int i;
	int ok = 0;
#ifdef ERROR_CONTROL
	if (vf == NULL)
	{
		ok = 0;
		return ok;
	}
#endif
	vf[0] = (1.0 + r) * (1.0 + s) + alpha[0] * (1.0 + s) + alpha[1] * (1.0 + r);
	vf[1] = (1.0 - r) * (1.0 + s) - alpha[0] * (1.0 + s) + alpha[1] * (1.0 - r);
	vf[2] = (1.0 - r) * (1.0 - s) - alpha[0] * (1.0 - s) - alpha[1] * (1.0 - r);
	vf[3] = (1.0 + r) * (1.0 - s) + alpha[0] * (1.0 - s) - alpha[1] * (1.0 + r);
	for (i = 0; i < 4; i++)
		vf[i] *= 0.25;
	return ok = 1;
}

/***************************************************************************
   ROCKFLOW - Funktion: MPhi3D_SUPG
   Aufgabe:
           Berechnet Testfunktion (r,s,t) incl. Upwind-Parameter alpha.
           Phi = Omega + alpha * grad Omega
   Formalparameter:
           Z: *vf - 1x8 Feld
           E: r,s,t (s.o.)
   Ergebnis:
           Vektor
   Aenderungen/Korrekturen:
   11/1995     cb        Erste Version
   01.07.1997  R.Kaiser  Umbenannt zu MPhi3D_SUPG

 **************************************************************************/

int MPhi3D_SUPG(double* vf, double r, double s, double t, double* alpha)
{
	int i;
	int ok = 0;
#ifdef ERROR_CONTROL
	if (vf == NULL)
	{
		ok = 0;
		return ok;
	}
#endif
	vf[0] = (1.0 + r) * (1.0 + s) * (1.0 + t) + alpha[0] * (1.0 + s) * (1.0 + t) + alpha[1] * (1.0 + r) * (1.0 + t)
	        + alpha[2] * (1.0 + r) * (1.0 + s);
	vf[1] = (1.0 - r) * (1.0 + s) * (1.0 + t) - alpha[0] * (1.0 + s) * (1.0 + t) + alpha[1] * (1.0 - r) * (1.0 + t)
	        + alpha[2] * (1.0 - r) * (1.0 + s);
	vf[2] = (1.0 - r) * (1.0 - s) * (1.0 + t) - alpha[0] * (1.0 - s) * (1.0 + t) - alpha[1] * (1.0 - r) * (1.0 + t)
	        + alpha[2] * (1.0 - r) * (1.0 - s);
	vf[3] = (1.0 + r) * (1.0 - s) * (1.0 + t) + alpha[0] * (1.0 - s) * (1.0 + t) - alpha[1] * (1.0 + r) * (1.0 + t)
	        + alpha[2] * (1.0 + r) * (1.0 - s);
	vf[4] = (1.0 + r) * (1.0 + s) * (1.0 - t) + alpha[0] * (1.0 + s) * (1.0 - t) + alpha[1] * (1.0 + r) * (1.0 - t)
	        - alpha[2] * (1.0 + r) * (1.0 + s);
	vf[5] = (1.0 - r) * (1.0 + s) * (1.0 - t) - alpha[0] * (1.0 + s) * (1.0 - t) + alpha[1] * (1.0 - r) * (1.0 - t)
	        - alpha[2] * (1.0 - r) * (1.0 + s);
	vf[6] = (1.0 - r) * (1.0 - s) * (1.0 - t) - alpha[0] * (1.0 - s) * (1.0 - t) - alpha[1] * (1.0 - r) * (1.0 - t)
	        - alpha[2] * (1.0 - r) * (1.0 - s);
	vf[7] = (1.0 + r) * (1.0 - s) * (1.0 - t) + alpha[0] * (1.0 - s) * (1.0 - t) - alpha[1] * (1.0 + r) * (1.0 - t)
	        - alpha[2] * (1.0 + r) * (1.0 - s);
	for (i = 0; i < 8; i++)
		vf[i] *= 0.125;
	return ok = 1;
}

/***************************************************************************
   ROCKFLOW - Funktion: MGradOmega2D
   Aufgabe:
           Berechnet Gradient eines 2D Vektorfeldes,
           dessen Ansatzfunktion (Vektor) bekannt ist (r,s).
                        /                            \
                     1  | +(1+s) -(1+s) -(1-s) +(1-s) |
      grad (omega)= --- |                             |
                     4  | +(1+r) +(1-r) -(1-r) -(1+r) |
 \                            /
   Formalparameter:
   Z: *vf - 2x4 Feld
   E: r,s (s.o.)
   Ergebnis:
   2x4 Matrix
   Aenderungen/Korrekturen:
   05/1995     hh        Erste Version

 **************************************************************************/

int MGradOmega2D(double* vf, double r, double s)
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
	vf[0] = +(1.0 + s);
	vf[1] = -(1.0 + s);
	vf[2] = -(1.0 - s);
	vf[3] = +(1.0 - s);
	vf[4] = +(1.0 + r);
	vf[5] = +(1.0 - r);
	vf[6] = -(1.0 - r);
	vf[7] = -(1.0 + r);
	for (i = 0; i < 8; i++)
		vf[i] *= 0.25;
	return ok = 1;
}

/***************************************************************************
   ROCKFLOW - Funktion: MGradOmega3D
   Aufgabe:
           Berechnet Gradient eines 3D Vektorfeldes,
           dessen Ansatzfunktion (Vektor) bekannt ist (r,s,t).

 |  0   1   2   3   4   5   6   7  |
 |                                 |
      grad (omega)= |  8   9  10  11  12  13  14  15  |
 |                                 |
 | 16  17  18  19  20  21  22  23  |

   Zahlen in der Matrix stehen fuer Positionen im Feld vf.
   Gleichungen s.u..

   Formalparameter:
   Z: *vf - 3x8 Feld
   E: r,s (s.o.)
   Ergebnis:
   3x8 Matrix
   Aenderungen/Korrekturen:
   05/1995     hh        Erste Version

 **************************************************************************/

int MGradOmega3D(double* vf, double r, double s, double t)
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
	vf[0] = +(1.0 + s) * (1.0 + t);
	vf[1] = -(1.0 + s) * (1.0 + t);
	vf[2] = -(1.0 - s) * (1.0 + t);
	vf[3] = +(1.0 - s) * (1.0 + t);

	vf[4] = +(1.0 + s) * (1.0 - t);
	vf[5] = -(1.0 + s) * (1.0 - t);
	vf[6] = -(1.0 - s) * (1.0 - t);
	vf[7] = +(1.0 - s) * (1.0 - t);

	vf[8] = +(1.0 + r) * (1.0 + t);
	vf[9] = +(1.0 - r) * (1.0 + t);
	vf[10] = -(1.0 - r) * (1.0 + t);
	vf[11] = -(1.0 + r) * (1.0 + t);

	vf[12] = +(1.0 + r) * (1.0 - t);
	vf[13] = +(1.0 - r) * (1.0 - t);
	vf[14] = -(1.0 - r) * (1.0 - t);
	vf[15] = -(1.0 + r) * (1.0 - t);

	vf[16] = +(1.0 + r) * (1.0 + s);
	vf[17] = +(1.0 - r) * (1.0 + s);
	vf[18] = +(1.0 - r) * (1.0 - s);
	vf[19] = +(1.0 + r) * (1.0 - s);

	vf[20] = -(1.0 + r) * (1.0 + s);
	vf[21] = -(1.0 - r) * (1.0 + s);
	vf[22] = -(1.0 - r) * (1.0 - s);
	vf[23] = -(1.0 + r) * (1.0 - s);

	for (i = 0; i < 24; i++)
		vf[i] *= 0.125;
	return ok = 1;
}

/***************************************************************************
   ROCKFLOW - Funktion: MGradPhi2D
   Aufgabe:
           Berechnet Gradient der Testfunktionen (r,s)
           grad phi = grad omega
   Formalparameter:
           Z: *vf - 2x4 Feld
           E: r,s (s.o.)
   Ergebnis:
           2x4 Matrix
   Aenderungen/Korrekturen:
   11/1995     cb        Erste Version

 **************************************************************************/

int MGradPhi2D(double* vf, double r, double s)
{
	return MGradOmega2D(vf, r, s);
}

/***************************************************************************
   ROCKFLOW - Funktion: MGradPhi3D
   Aufgabe:
           Berechnet Gradient der Testfunktionen (r,s,t).
           grad phi = grad omega
   Formalparameter:
           Z: *vf - 3x8 Feld
           E: r,s (s.o.)
   Ergebnis:
           3x8 Matrix
   Aenderungen/Korrekturen:
   11/1995     hh        Erste Version

 **************************************************************************/

int MGradPhi3D(double* vf, double r, double s, double t)
{
	return MGradOmega3D(vf, r, s, t);
}

/***************************************************************************
   ROCKFLOW - Funktion: MGetCoor
   Aufgabe:
            lokale Koordinaten der Eckpunkte
   Formalparameter:
           E: typ Element-Dimension - 1
           E: nr  lokale Knotennummer
           X: r,s,t lokale Koordinaten
   Ergebnis:
           -void-
   Aenderungen/Korrekturen:
   04/1996     cb        Erste Version

 **************************************************************************/
void MGetCoor(int typ, long j, double* r, double* s, double* t)
{
	switch (typ)
	{
		case 0:
			*r = (double)(1l - j - j);
			*s = *t = 0.0;
			break;
		case 1:
			*r = 1.0 - 2.0 * (double)(((j + 1l) % 4) / 2l);
			*s = 1.0 - 2.0 * (double)(j / 2l);
			*t = 0.0;
			break;
		case 2:
			*r = 1.0 - 2.0 * (double)(((j + 1l) % 4) / 2l);
			*s = 1.0 - 2.0 * (double)((j % 4) / 2l);
			*t = 1.0 - 2.0 * (double)(j / 4l);
			break;
	} /* switch typ */
}

/**************************************************************************
   ROCKFLOW - Function: MAngleVectors

   Task:  Calculate angle between 2 vectors

   Parameter: (I: Input; R: Return; X: Both)
           I: *v1, *v2

   Return:
           Angle

   Programming:
   09/2002   MB   First Version
**************************************************************************/
double MAngleVectors(double* v1, double* v2)
{
	double Pproduct;
	double angle;

	Pproduct = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
	angle = (acos(Pproduct / (MBtrgVec(v1, 3) * MBtrgVec(v2, 3)))) * 180 / PI;

	return angle;
}

/***************************************************************************
   ROCKFLOW - Funktion: MNormiere
   Aufgabe:
           Berechnet Norm von Vektor
   Formalparameter:
           E: *vec
           E: n
   Ergebnis:
           normierter Vektor
   Aenderungen/Korrekturen:
   07/1995     hh        Erste Version
   11/1999     C.Thorenz Register-Variablen
 **************************************************************************/

void MNormiere(double* vec, long n)
{
	register long i;
	register double vl;
	vl = MBtrgVec(vec, n);
	for (i = 0; i < n; i++)
		vec[i] = vec[i] / vl;
}

/***************************************************************************
   ROCKFLOW - Funktion: M2Determinante
   Aufgabe:
           Berechnung der Determinante einer 2x2-Matrix
   Formalparameter:
           E: *matrix
   Ergebnis:
           Determinante
   Aenderungen/Korrekturen:
   07/1994     hh        Erste Version

 **************************************************************************/

double M2Determinante(double* matrix)
{
	return matrix[0] * matrix[3] - matrix[1] * matrix[2];
} /* extern double M2Determinante */

/***************************************************************************
   ROCKFLOW - Funktion: Mxg2Determinante
   Aufgabe:
           Berechnung der Determinante einer mxn-Matrix (n=m und n>2)
                    --- n ---       ---- n ---------
 |*********      | 0 1 2 3 4 ... n
 |********* oder | n+1 n+2 ...
                   m*********      m
 |*********      |              ...
 |*********      |         ...  n*m-1
   Formalparameter:
   E: *matrix
   E: m,n - Ausdehnung
   Ergebnis:
   Determinante oder 0.0 bei Fehler
   Aenderungen/Korrekturen:
   07/1994     hh        Erste Version
   11/1999     C.Thorenz Register-Variablen
 **************************************************************************/

double Mxg2Determinante(double* matrix, long m, long n)
{
	register long i, k, sprung;
	register double depp = 0.0, dussel;
#ifdef ERROR_CONTROL
	if (m != n)
	{
		DisplayErrorMsg("Determinate einer nicht-quadratischen Matrix ?");
		return 0.0;
	} /* if */
	if (m < 3)
	{
		DisplayErrorMsg("Determinate zu klein m>2 !!");
		return 0.0;
	} /* if */
#endif

	dussel = 1.0;
	for (i = 0; i < m * n; i += m + 1)
		dussel *= matrix[i];
	depp += dussel;
	for (k = 1; k < m; k++)
	{
		dussel = 1.0;
		sprung = (m - k) * n - 1;
		for (i = k; i < n * m; i += m + 1)
		{
			dussel *= matrix[i];
			if (sprung == i)
				i -= m;
		} /* for */
		depp += dussel;
	} /* for */
	dussel = 1.0;
	for (i = n - 1; i < m * n - m + 1; i += m - 1)
		dussel *= matrix[i];
	depp -= dussel;
	for (k = n - 2; k > -1; k--)
	{
		dussel = 1.0;
		sprung = k * n;
		for (i = k; i < n * m; i += m - 1)
		{
			dussel *= matrix[i];
			if (sprung == i)
				i += m;
		} /* for */
		depp -= dussel;
	} /* for */
	return depp;
} /* extern double Mxg2Determinante */

/**************************************************************************
   ROCKFLOW - Funktion: MTranspoVec
   Aufgabe:
           Vektor transponieren
   Formalparameter:
           X: *vec
           E: g
   Ergebnis:
           Transponierte
   Aenderungen/Korrekturen:
   09/1994     hh        Erste Version
   11/1999     C.Thorenz Register-Variablen
**************************************************************************/

void MTranspoVec(double* vec, long g)
{
	register long i;
	register double zwiebel;
	for (i = 0; i < g / 2; i++)
	{
		zwiebel = vec[i];
		vec[i] = vec[g - 1 - i];
		vec[g - 1 - i] = zwiebel;
	}
} /* MTranspoVec */

/**************************************************************************
   ROCKFLOW - Funktion: MTranspoMat
   Aufgabe:
           Matrizen transponieren
   Formalparameter:
           E: *mat1
           E: m=Zeilenzahl
           E: n=Spaltenzahl
           X: *mat2
   Ergebnis:
           Transponierte
   Aenderungen/Korrekturen:
   08/1994     hh        Erste Version
   08/1995     cb        (so geht das irgendwie nicht)
   11/1999     C.Thorenz Register-Variablen
**************************************************************************/

void MTranspoMat(double* mat1, long m, long n, double* mat2)
{
	register long i, k;
	for (i = 0; i < n; i++)
		for (k = 0; k < m; k++)
			mat2[i * m + k] = mat1[k * n + i];
} /* MTranspoMat */

/**************************************************************************
   ROCKFLOW - Funktion: M2InvertiereUndTransponiere
   Aufgabe:
           2x2 Matrizen invertieren |a b|
  |c d|
                        1     |  d  -b |
           inv m = ---------- |        |
                     det(m)   | -c   a |
   Formalparameter:
           X: *m
   Ergebnis:
   Invertierte und Transponierte
   Aenderungen/Korrekturen:
   08/1994     hh        Erste Version
   11/1999     C.Thorenz Register-Variablen
   08/2001     MK        M2Invertiere M2InvertiereUndTransponiere umbenannt
**************************************************************************/

void M2InvertiereUndTransponiere(double* m)
{
	register double eddet, zecke;
	register int i;
	eddet = m[0] * m[3] - m[1] * m[2];
	if (fabs(eddet) > MKleinsteZahl)
		eddet = 1.0 / eddet;
	zecke = m[0];
	m[0] = m[3];
	m[3] = zecke;
	zecke = -m[1];
	m[1] = -m[2]; /* hier wird auch transponiert */
	m[2] = zecke;
	for (i = 0; i < 4; i++)
		m[i] *= eddet;
} /* M2InvertiereUndTransponiere */

/**************************************************************************
   ROCKFLOW - Funktion: M2Invertiere
   Aufgabe:
           2x2 Matrizen invertieren |a b|
  |c d|
                        1     |  d  -b |
           inv m = ---------- |        |
                     det(m)   | -c   a |
   Formalparameter:
           X: *m
   Ergebnis:
   Invertierte
   Aenderungen/Korrekturen:
   08/1994   hh           Erste Version
   11/1999   C.Thorenz    Register-Variablen
   08/2001   M. Kohlmeier Funktion transponierte und invertierte jetzt
   invertiert sie nur. Alle Aufrufe ersetzt durch
   M2InvertiereUndTransponiere

**************************************************************************/

void M2Invertiere(double* m)
{
	register double eddet, zecke;
	register int i;
	eddet = m[0] * m[3] - m[1] * m[2];
	if (fabs(eddet) > MKleinsteZahl)
		eddet = 1.0 / eddet;
	zecke = m[0];
	m[0] = m[3];
	m[3] = zecke;
	m[1] = -m[1];
	m[2] = -m[2];
	for (i = 0; i < 4; i++)
		m[i] *= eddet;
} /* M2Invertiere */

/**************************************************************************
   ROCKFLOW - Funktion: M3Invertiere
   Aufgabe:
           3x3 Matrizen invertieren

   Formalparameter:
           X: *m
   Ergebnis:
           Invertierte
   Aenderungen/Korrekturen:
   06/1996     cb        Erste Version nach RaRa
   11/1999     C.Thorenz Register-Variablen
**************************************************************************/

void M3Invertiere(double* m)
{
	double z[9];
	register double d;
	register int i;
	d = m[0] * (m[4] * m[8] - m[7] * m[5]) + m[3] * (m[7] * m[2] - m[1] * m[8]) + m[6] * (m[1] * m[5] - m[4] * m[2]);
	if (fabs(d) > MKleinsteZahl)
		d = 1.0 / d;
	z[0] = m[4] * m[8] - m[7] * m[5];
	z[3] = m[6] * m[5] - m[3] * m[8];
	z[6] = m[3] * m[7] - m[6] * m[4];
	z[1] = m[7] * m[2] - m[1] * m[8];
	z[4] = m[0] * m[8] - m[6] * m[2];
	z[7] = m[1] * m[6] - m[0] * m[7];
	z[2] = m[1] * m[5] - m[4] * m[2];
	z[5] = m[3] * m[2] - m[0] * m[5];
	z[8] = m[0] * m[4] - m[3] * m[1];
	for (i = 0; i < 9; i++)
		m[i] = d * z[i];
} /* M3Invertiere */

/**************************************************************************
   ROCKFLOW - Funktion: MInvertiere
   Aufgabe:
           Matrizen invertieren (Hauptdiagonalelemente != Null)
           (Verfahren aus Schneider)
   Formalparameter:
           X: *matrix
           E: m,n
   Ergebnis:
           Invertierte
   Aenderungen/Korrekturen:
   04/1996     cb        Erste Version
   24.06.1999  OK        Warnung bei nichtquadratischen Matrizen
   11/1999     C.Thorenz Register-Variablen, Warnung in #ifdef

**************************************************************************/
void MInvertiere(double* mat, long m, long n)
{
	register long i, j, k;

#ifdef ERROR_CONTROL
	if (m != n)
	{
		DisplayErrorMsg("Abbruch !!! Nicht-quadratische Matrix *** MInvertiere");
		abort();
	}
	for (k = 0; k < n; k++)
		if (fabs(mat[n * k + k]) < MKleinsteZahl)
		{
			DisplayErrorMsg("Abbruch !!! Diagonalelement Null in MInvertiere");
			abort();
		}
#endif
	m = n; /* ah */
	for (k = 0; k < n; k++)
	{
		for (j = 0; j < n; j++)
			if (j != k)
			{
				mat[n * k + j] /= (-mat[n * k + k]);
				for (i = 0; i < n; i++)
					if (i != k)
						mat[n * i + j] += (mat[n * k + j] * mat[n * i + k]);
			}
		mat[n * k + k] = 1.0 / mat[n * k + k];
		for (i = 0; i < n; i++)
			if (i != k)
				mat[n * i + k] *= mat[n * k + k];
	}
} /* MInvertiere */

/**************************************************************************
   ROCKFLOW - Funktion: MAddVektoren
   Aufgabe:
           Vektorenen addieren
   Formalparameter:
           E: *v1, *v2
           A: *vout
           E: g
   Ergebnis:
           Vektorsumme
   Aenderungen/Korrekturen:
   08/1994     hh        Erste Version
   11/1999     C.Thorenz Register-Variablen
   10/2005 PA Vectorization
**************************************************************************/
#ifdef SX
int MAddVektoren(double* restrict v1, double* restrict v2, double* restrict vout, long g)
#else
int MAddVektoren(double* v1, double* v2, double* vout, long g)
#endif
{
	register long i;
#ifdef SX
#pragma cdir nodep
#endif
	for (i = 0; i < g; i++)
		vout[i] = v1[i] + v2[i];
	return 1;
} /* MAddVektoren */

/**************************************************************************
   ROCKFLOW - Funktion: MAddMatrizen
   Aufgabe:
           Matrizen addieren
   Formalparameter:
           E: *m1, *m2
           A: *mout
           E: m,n
   Ergebnis:
           Matrixsumme
   Aenderungen/Korrekturen:
   08/1994     hh        Erste Version
   11/1999     C.Thorenz Register-Variablen
**************************************************************************/
#ifdef SX
int MAddMatrizen(double* restrict m1, double* restrict m2, double* restrict mout, long m, long n)
#else
int MAddMatrizen(double* m1, double* m2, double* mout, long m, long n)
#endif
{
	register long i;
#ifdef SX
#pragma cdir nodep
#endif
	for (i = 0; i < m * n; i++)
		mout[i] = m1[i] + m2[i];
	return 1;
}

/**************************************************************************
   ROCKFLOW - Funktion: MMultVecSkalar
   Aufgabe:
           Multiplikation eines Vektors mit einem Skalarwert,
           wobei der Vektor (!) veraendert wird
   Formalparameter:
           X: *vec - Zeiger auf Vektor
           E: skal - Skalarwert
           E: g - Ausdehnung
   Aenderungen/Korrekturen:
   07/1994     hh        Erste Version
   11/1999     C.Thorenz Register-Variablen
**************************************************************************/
#ifdef SX
int MMultVecSkalar(double* restrict vec, double skal, long g)
#else
int MMultVecSkalar(double* vec, double skal, long g)
#endif
{
	register long i;
#ifdef SX
#pragma cdir nodep
#endif
	for (i = 0; i < g; i++)
		vec[i] *= skal;
	return 1;
} /* extern int MMultVecSkalar */

/**************************************************************************
   ROCKFLOW - Funktion: MMultMatSkalar
   Aufgabe:
           Multiplikation einer Matrix mit einem Skalarwert,
           wobei die Matrix (!) veraendert wird
   Formalparameter:
           X: *matrix - Zeiger auf Matrix
           E: skal - Skalarwert
           E: m,n - Ausdehnung
   Aenderungen/Korrekturen:
   07/1994     hh        Erste Version
   11/1999     C.Thorenz Register-Variablen
**************************************************************************/
#ifdef SX
int MMultMatSkalar(double* restrict matrix, double skal, long m, long n)
#else
int MMultMatSkalar(double* matrix, double skal, long m, long n)
#endif
{
	register long i;
#ifdef SX
#pragma cdir nodep
#endif
	for (i = 0; i < m * n; i++)
		matrix[i] *= skal;
	return 1;
}

/*##########################################################################
   Auf Matrizen und Vektoren aufbauende Funktionen
 ######################################################################## */

/**************************************************************************
   ROCKFLOW - Funktion: MKTF2Dr2D
   Aufgabe:
           Koordinaten-Transformation 2D -> 2D (r=rotiert)
   Formalparameter:
           E: *vec1
           E: winkel
           A: *vec2
   Ergebnis:
           gedrehter Vektor vec2
           Transformation von System 1 zu System 2
   x1
   /\    / x2
|   /          um winkel verdrehtes Koordinatensystem
|  /           es wird immer nach rechts gedreht
| /
|||/
   -+----------->
|||..          y1
   ... y2

   Aenderungen/Korrekturen:
   07/1995     hh        Erste Version
**************************************************************************/

int MKTF2Dr2D(double* vec1, double winkel, double* vec2)
{
#ifdef ERROR_CONTROL
	if (vec2 == NULL)
	{
		DisplayErrorMsg("MKTF2Dr2D:vec2 muss mit malloc/MMachMat def. sein");
		return 0;
	}
	if (vec1 == NULL)
	{
		DisplayErrorMsg("MKTF2Dr2D:vec1 muss mit malloc/MMachMat def. sein");
		return 0;
	}
#endif
	vec2[0] = vec1[0] * cos(winkel) + vec1[1] * sin(winkel);
	vec2[1] = -vec1[0] * sin(winkel) + vec1[1] * cos(winkel);
	return 1;
}

/**************************************************************************
   ROCKFLOW - Funktion: MKTFMat2Dr2D
   Aufgabe:
           Berchnung einer Transformationsmatrix fuer die
           Koordinaten-Transformation 2D -> 2D (r=rotiert)
   Formalparameter:
           E: winkel
           A: *tmat
   Ergebnis:
           Transformationmatrix von System 1 zu System 2
       x1
   /\    / x2
|   /          um winkel verdrehtes Koordinatensystem
|  /           es wird immer nach rechts gedreht
| /
|||/
   -+----------->
|||..          y1
   ... y2

   Aenderungen/Korrekturen:
   07/1995     hh        Erste Version
**************************************************************************/

int MKTFMat2Dr2D(double winkel, double* tmat)
{
#ifdef ERROR_CONTROL
	if (tmat == NULL)
	{
		DisplayErrorMsg("MKTFMat2Dr2D:tmat muss mit malloc/MMachMat def. sein");
		return 0;
	}
#endif
	tmat[0] = cos(winkel);
	/* tmat[1]= sin(winkel);
	   tmat[2]=-sin(winkel); */
	tmat[1] = -sin(winkel);
	tmat[2] = sin(winkel);
	tmat[3] = cos(winkel);
	return 1;
}

/**************************************************************************
   ROCKFLOW - Funktion: MKTF2Dt2D
   Aufgabe:
           Koordinaten-Transformation 2D -> 2D (t=transversal)
   Formalparameter:
           E: *vec1
           E: dx,dy
           A: *vec2
   Ergebnis:
           verschobener Vektor vec2
           Transformation von System 1 zu System 2
   x1
   /\    /\ x2
|     |
|||<--->|
| dx  |
|     +------------->
   -+-----|----->       y2
|     |      y1

   Aenderungen/Korrekturen:
   07/1995     hh        Erste Version
**************************************************************************/

int MKTF2Dt2D(double* vec1, double dx, double dy, double* vec2)
{
	vec2[1] = vec1[1] + dx;
	vec2[2] = vec1[2] + dy;
	return 1;
}

/**************************************************************************
   ROCKFLOW - Funktion: MKTF3Dr2D
   Aufgabe:
           Koordinaten-Transformation 3D-Ebene -> 2D (r=rotiert)
   Formalparameter:
           E: *vec1
           E: *vec2
           E: winkel
           A: *vec
   Ergebnis:
           gedrehter Vektor vec
   Transformation von System 1 (in 3D-Ebene)
   zu System 2 (nur noch 2D)
   x1
   /\    / x2
|   /          um winkel verdrehtes Koordinatensystem
|  /           es wird immer nach rechts gedreht
| /
|||/
   -+----------->
|||..          y1
   ... y2

   Augegeben wird der gedrehte Vektor vec1

   Aenderungen/Korrekturen:
   07/1995     hh        Erste Version
   11/1995     msr       Malloc raus / static rein
**************************************************************************/

int MKTF3Dr2D(double* vec1, double* vec2, double winkel, double* vec)
{
	static double f;
	static double zwick[6];
	/* static double zwack[6]; */
	static double zwock[4];
	static long i;
	MNormiere(vec1, 3);
	f = MSkalarprodukt(vec1, vec2, 3);
	for (i = 0; i < 3; i++)
		vec2[i] -= f * vec1[i];
	MNormiere(vec2, 3);
	zwick[0] = vec1[0];
	zwick[1] = vec2[0];
	zwick[2] = vec1[1];
	zwick[3] = vec2[1];
	zwick[4] = vec1[2];
	zwick[5] = vec2[2];
	MKTFMat2Dr2D(winkel, zwock);
	MMultMatMat(zwick, 3, 2, zwock, 2, 2, vec, 3, 2);
	/* MMultMatMat(zwick,3,2,zwock,2,3,zwack,3,2);
	   zwock[0]=vec[0];zwock[1]=vec[1];
	   MMultMatVec(zwack,3,2,zwock,2,vec,2); */
	return 1;
}

/**************************************************************************
   ROCKFLOW - Funktion: MKTFMat3Dr2D
   Aufgabe:
          Koordinaten-Transformation 3D-Ebene -> 2D (r=rotiert)
   Formalparameter:
          E: *vec1
          E: *vec2
          E: winkel
          A: *vec
   Ergebnis:
          gedrehter Vektor vec
   Transformation von System 1 (in 3D-Ebene)
   zu System 2 (nur noch 2D)
   x1
   /\    / x2
|   /          um winkel verdrehtes Koordinatensystem
|  /           es wird immer nach rechts gedreht
| /
|||/
   -+----------->
|||..          y1
   ... y2

   Augegeben wird der gedrehte Vektor vec1

   Aenderungen/Korrekturen:
   07/1995     hh        Erste Version
   11/1995     msr       Malloc raus / static rein
**************************************************************************/

void MKTFMat3Dr2D(double* vec1, double* vec2, double winkel, double* mat)
{ /* if(winkel!=0.0) */ /* diese Verzweigung spart viel Zeit */
	{
		static double f;
		static double zwick[6];
		static double zwack[6];
		static double zwock[4];
		static long i;
		MNormiere(vec1, 3);
		f = MSkalarprodukt(vec1, vec2, 3);
		for (i = 0; i < 3; i++)
			vec2[i] -= f * vec1[i];
		MNormiere(vec2, 3);
		zwick[0] = vec1[0];
		zwick[1] = vec2[0];
		zwick[2] = vec1[1];
		zwick[3] = vec2[1];
		zwick[4] = vec1[2];
		zwick[5] = vec2[2];
		MKTFMat2Dr2D(winkel, zwock);
		MMultMatMat(zwick, 3, 2, zwock, 2, 2, zwack, 3, 2);
		for (i = 0; i < 6; i++)
			mat[i] = zwack[i];
	} /* if */

	/* else
	   {MNulleMat(mat,3,2);
	   mat[0]=1.0;mat[3]=1.0;
	   } */
	/* else */
}

/*##########################################################################
   Pruefunktion fuer Matrix
 ######################################################################## */

/**************************************************************************
   ROCKFLOW - Funktion: MBistDuDiagMat
   Aufgabe:
           Prueft Matrix auf diagonalitaet
   Formalparameter:
           E: *matrix
           E: m,n
   Ergebnis:
           -1  - keine quadartische Matrix
            0  - keine Diagonalmatrix
            1  - Diagonalmatrix
   3  - Tridiagonalmatrix
   Aenderungen/Korrekturen:
   09/1994     hh        Erste Version
**************************************************************************/

int MBistDuDiagMat(double* matrix, long m, long n)
{
	long i, j, knall = 0, diag = 0;
	if (m != n)
		return -1;
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			/*printf(" %ld %ld \n",i,j); */
			if ((fabs(matrix[i + j * n]) > MKleinsteZahl) && ((i - j < -1) || (i - j > 1)))
				return 0;

	diag = 3;
	knall = 0;
	for (i = 0; i < m * n; i += m + 1)
		for (j = i + 1; ((j < i + m) && (j < m * n)); j++)
			if (fabs(matrix[j]) > MKleinsteZahl)
				knall = 1;
	if (!knall)
		diag = 1;
	return diag;
} /* of MBistDuDiagMat */

/*##########################################################################
   Bearbeitungfunktionen fuer Vektoren
 ######################################################################## */

/**************************************************************************
   ROCKFLOW - Funktion: MLoeschVec
   Aufgabe:
           Zerstoeren eines bestehenden Vektors
   Formalparameter:
           E: *vec
   Ergebnis:
           Rueckgabewert 1
   Aenderungen/Korrekturen:
   08/1994     hh        Erste Version
**************************************************************************/

void MLoeschVec(double* vec)
{
#ifdef ERROR_CONTROL
	if (vec == NULL)
		DisplayErrorMsg("Fehler in MLoeschVec");
#endif
	vec = (double*)Free(vec);
} /* MLoeschVec */

/*##########################################################################
   Bearbeitungfunktionen fuer Matrizen
 ######################################################################## */

/**************************************************************************
   ROCKFLOW - Funktion: MMachMat
   Aufgabe:
           Erzeugen einer neuen Matrix
   Formalparameter:
           E: m,n
   Ergebnis:
           Zeiger auf die erzeugte Matrix.
           oder Rueckgabewert NULL wenn nicht genuegend Speicher da ist.
   Aenderungen/Korrekturen:
   08/1994     hh        Erste Version
   11/1999     C.Thorenz Nullung herausgenommen
**************************************************************************/

double* MMachMat(long m, long n)
{
	double* zwerg;
	zwerg = (double*)Malloc(sizeof(double) * m * n);
#ifdef ERROR_CONTROL
	if (zwerg == NULL)
		DisplayErrorMsg("zuwenig Speicher in MMachMat");
#endif
	return zwerg;
} /* MMachMat */

/**************************************************************************
   ROCKFLOW - Funktion: MNullMat
   Aufgabe:
           Fuellt Matrix mit Nullen
   Formalparameter:
           E: *zwerg,m,n
   Ergebnis:
   Aenderungen/Korrekturen:
   11/1999     C.Thorenz Erste Version
**************************************************************************/

void MNullMat(double* zwerg, long m, long n)
{
	register long i;
#ifdef ERROR_CONTROL
	if (zwerg == NULL)
		DisplayErrorMsg("Fehler in MNullMat");
#endif
	for (i = 0; i < m * n; i++)
		zwerg[i] = 0.0;
} /* MNullMat */

/**************************************************************************
   ROCKFLOW - Funktion: MLoeschMat
   Aufgabe:
           Zerstoeren einer bestehenden Matrix
   Formalparameter:
           E: *vec
   Ergebnis:
           Rueckgabewert 1
   Aenderungen/Korrekturen:
   08/1994     hh        Erste Version
**************************************************************************/

void MLoeschMat(double* mat)
{
#ifdef ERROR_CONTROL
	if (mat == NULL)
		DisplayErrorMsg("Fehler in MLoeschMat");
#endif
	mat = (double*)Free(mat);
} /* MLoeschMat */

/**************************************************************************
   ROCKFLOW - Funktion: MKopierMat
   Aufgabe:
           Kopieren einer Matrix
   Formalparameter:
           E: *matq
           A: *mat
           E: m
           E: n
   Ergebnis:
           Rueckgabewert 1
   Aenderungen/Korrekturen:
   08/1994     hh        Erste Version
**************************************************************************/

void MKopierMat(double* matq, double* matz, long m, long n)
{
	register long i;
	for (i = 0; i < m * n; i++)
		matz[i] = matq[i];
} /* MKopierMat */

/*##########################################################################
   Ausgabefunktionen
 ########################################################################### */
/**************************************************************************
   ROCKFLOW - Funktion: MZeigVec
   Aufgabe:
           Anzeigen eines Vektors
   Formalparameter:
           E: *vec
           E: grad
           E: *text
   Programmaenderungen:
   08/1994     hh        Erste Version
**************************************************************************/
void MZeigVec(double* vec, long grad, char* text)
{
	long i;
	printf("\n%s\n", text);
	for (i = 0; i < grad; i++)
		printf(" | %e | \n", vec[i]);
}

/**************************************************************************
   ROCKFLOW - Funktion: M2FileVec
   Aufgabe:
           Anzeigen eines Vektors
   Formalparameter:
           E: *vec
           E: grad
           E: *text
   Programmaenderungen:
   08/1994     hh        Erste Version
**************************************************************************/
void M2FileVec(double* vec, long grad, char* text)
{
	long i;
	FILE* fp;
	fp = fopen("fgnuplot.asc", "a");
	fprintf(fp, "\n%s\n", text);
	for (i = 0; i < grad; i++)
		fprintf(fp, "%e\n", vec[i]);
	fclose(fp);
}

/**************************************************************************
   ROCKFLOW - Funktion: MZeigMat
   Aufgabe:
           Anzeigen einer Matrix
   Formalparameter:
           E: *mat
           E: m=Zeilenzahl
           E: n=Spaltenzahl
           E: *text
   Programmaenderungen:
   08/1994     hh        Erste Version
**************************************************************************/
void MZeigMat(double* mat, long m, long n, char* text)
{
	long i, k;
	printf("\n%s\n", text);
	for (k = 0; k < m; k++)
	{
		printf("| ");
		for (i = 0; i < n; i++)
			printf(" %e ", mat[i + k * n]);
		puts(" |");
	}
} /* MZeigMat */

/**************************************************************************
   ROCKFLOW - Funktion: MVekNormMax

   Aufgabe:
   Berechnet die Maximumnorm von x

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E double *x : Zeiger auf Vektor
   E long n : Dimension von x

   Ergebnis:
   Norm

   Programmaenderungen:
   11/1995     MSR        Erste Version

**************************************************************************/
double MVekNormMax(double* x, long n)
{
	register long i;
	register double erg = fabs(x[0]);
	for (i = 1l; i < n; i++)
		if (fabs(x[i]) > erg)
			erg = fabs(x[i]);
	return erg;
}

/*##########################################################################
   Geometrische Funktionen im R3
 ######################################################################## */

/***************************************************************************
   ROCKFLOW - Funktion: MCalcDistancePointToPoint

   Aufgabe:
           Abstand zweier Punkte im R3

   Formalparameter:
           E: *pt1  - Punktkoordinaten
           E: *pt2  - Punktkoordinaten

   Ergebnis:
   Abstand

   Aenderungen/Korrekturen:
   02/2000   C. Thorenz      Erste Version

 **************************************************************************/

/***************************************************************************
   ROCKFLOW - Funktion:  MCalcProjectionOfPointOnLine

   Aufgabe:
           Ermittelt den Projektionspunkt eines Punktes auf eine Linie (Lotpunkt)

   Formalparameter:
           E: *pt   - Punktkoordinaten
           E: *l1   - Punktkoordinaten fuer ersten Punkt auf Linie
           E: *l2   - Punktkoordinaten fuer zweiten Punkt auf Linie
           R: *proj   - Punktkoordinaten fuer Projektionspunkt auf Linie
   Ergebnis:
   Abstand Punkt-Linie

   Aenderungen/Korrekturen:
   02/2000   C. Thorenz      Erste Version

 **************************************************************************/
double MCalcProjectionOfPointOnLine(double* pt, double* l1, double* l2, double* proj)
{
	int i;
	double vec1[3], vec2[3], l, projektionslaenge, abstand;

	/* Vektor Linienpunkt zu Punkt */
	for (i = 0; i < 3; i++)
		vec1[i] = pt[i] - l1[i];

	/* Vektor Linienpunkt zu Linienpunkt */
	for (i = 0; i < 3; i++)
		vec2[i] = l2[i] - l1[i];

	/* Normieren */
	l = MBtrgVec(vec2, 3);
#ifdef ERROR_CONTROL
	if (l < MKleinsteZahl)
	{
		printf("\n Fehler in MCalcProjectionOfPointOnLine: Laenge ist Null !!!");
		exit(1);
	}
#endif
	for (i = 0; i < 3; i++)
		vec2[i] /= (l + MKleinsteZahl);

	projektionslaenge = MSkalarprodukt(vec1, vec2, 3);

	/* Punkt bestimmen */
	for (i = 0; i < 3; i++)
		proj[i] = l1[i] + projektionslaenge * vec2[i];

	abstand = MCalcDistancePointToPoint(proj, pt);

	return abstand;
}

/***************************************************************************
   ROCKFLOW - Funktion:  MCalcProjectionOfPointOnPlane

   Aufgabe:
           Ermittelt den Projektionspunkt eines Punktes auf einer Flaeche (Lotpunkt)

   Formalparameter:
           E: *pt  - Punktkoordinaten
           E: *e1  - Punktkoordinaten fuer ersten Punkt auf Ebene
           E: *e2  - Punktkoordinaten fuer zweiten Punkt auf Ebene
           E: *e3  - Punktkoordinaten fuer dritten Punkt auf Ebene
   R: *proj   - Punktkoordinaten fuer Projektionspunkt auf Linie
   Ergebnis:
   Abstand-Flaeche   (Kann negativ sein, je nach Lage des Punkts zur Ebene)

   Aenderungen/Korrekturen:
   02/2000   C. Thorenz      Erste Version

 **************************************************************************/

/***************************************************************************
   ROCKFLOW - Funktion: IsNodeInsideTriangle

   Aufgabe: Ueberprueft, ob sich ein vorgegebener Knoten innerhalb
            eines vorgegebenen Dreiecks befindet

   Formalparameter:
           E long n1, n2, n3: Knoten des Dreiecks
           E long node: Knoten

   Ergebnis:
   0 bei Fehler, sonst 1

   Aenderungen/Korrekturen:
   05/2000     RK        Erste Version

 **************************************************************************/

/***************************************************************************
   ROCKFLOW - Funktion: CalcTriangleArea

   Aufgabe: Berechnet die Flaeche eines vorgegebenen Dreiecks befindet

   Formalparameter:
           E long n1, n2, n3: Knoten des Dreiecks

   Ergebnis: Dreiecksflaeche

   Aenderungen/Korrekturen:
   05/2000     RK        Erste Version

 **************************************************************************/
double CalcTriangleArea(long n1, long n2, long n3)
{
	double area = 0.0;
	n1 = n1;
	n2 = n2;
	n3 = n3; // OK411
	/*OK411
	   double vec1[3], vec2[3], vec3[3];

	   vec1[0] = GetNodeX(n2)-GetNodeX(n1);
	   vec1[1] = GetNodeY(n2)-GetNodeY(n1);
	   vec1[2] = GetNodeZ(n2)-GetNodeZ(n1);

	   vec2[0] = GetNodeX(n1)-GetNodeX(n3);
	   vec2[1] = GetNodeY(n1)-GetNodeY(n3);
	   vec2[2] = GetNodeZ(n1)-GetNodeZ(n3);

	   M3KreuzProdukt(vec1,vec2,vec3);

	   area = 0.5 * MBtrgVec(vec3,3);
	 */
	return area;
}

/***************************************************************************
   ROCKFLOW - Funktion: IsNodeInsidePlain

   Aufgabe: Ueberprueft, ob sich ein vorgegebener Knoten innerhalb
            einer vorgegebenen Flaeche befindet

   Formalparameter:
           E long n1, n2, n3: Flaeche aufspannenden Knoten
           E long node: Knoten

   Ergebnis:
   0 bei Fehler, sonst 1

   Aenderungen/Korrekturen:
   05/2000     RK        Erste Version

 **************************************************************************/

/*##########################################################################
   Funktionen fuer zweidimensionale 9-Knoten-Elemente
 ######################################################################## */

/***************************************************************************
   ROCKFLOW - Funktion: MOmega2D_9N
   Aufgabe:
           Berechnet Ansatzfunktion (r,s).

   Formalparameter:
           Z: *vf - 1x9 Feld
           E: r,s (s.o.)
   Ergebnis:
           Vektor
   Aenderungen/Korrekturen:
   12/1999     RK        Erste Version

 **************************************************************************/

int MOmega2D_9N(double* vf, double r, double s)
{
	int ok = 0;
#ifdef ERROR_CONTROL
	if (vf == NULL)
	{
		ok = 0;
		return ok;
	}
#endif
	vf[8] = (1.0 - r * r) * (1.0 - s * s);
	vf[7] = 0.5 * (1.0 - s * s) * (1.0 + r) - 0.5 * vf[8];
	vf[6] = 0.5 * (1.0 - r * r) * (1.0 - s) - 0.5 * vf[8];
	vf[5] = 0.5 * (1.0 - s * s) * (1.0 - r) - 0.5 * vf[8];
	vf[4] = 0.5 * (1.0 - r * r) * (1.0 + s) - 0.5 * vf[8];
	vf[3] = 0.25 * (1.0 + r) * (1.0 - s) - 0.5 * vf[6] - 0.5 * vf[7] - 0.25 * vf[8];
	vf[2] = 0.25 * (1.0 - r) * (1.0 - s) - 0.5 * vf[5] - 0.5 * vf[6] - 0.25 * vf[8];
	vf[1] = 0.25 * (1.0 - r) * (1.0 + s) - 0.5 * vf[4] - 0.5 * vf[5] - 0.25 * vf[8];
	vf[0] = 0.25 * (1.0 + r) * (1.0 + s) - 0.5 * vf[4] - 0.5 * vf[7] - 0.25 * vf[8];
	return ok = 1;
}

/***************************************************************************
   ROCKFLOW - Funktion: MPhi2D_9N

   Aufgabe: Berechnet Testfunktion (r,s)

   Formalparameter: Z: *vf - 1x9 Feld
                    E: r,s (s.o.)

   Ergebnis: Vektor

   Aenderungen/Korrekturen:
   12/1999     R.Kaiser        Erste Version

 **************************************************************************/

int MPhi2D_9N(double* vf, double r, double s)
{
	int ok = 0;

	vf[8] = (1.0 - r * r) * (1.0 - s * s);
	vf[7] = 0.5 * (1.0 - s * s) * (1.0 + r) - 0.5 * vf[8];
	vf[6] = 0.5 * (1.0 - r * r) * (1.0 - s) - 0.5 * vf[8];
	vf[5] = 0.5 * (1.0 - s * s) * (1.0 - r) - 0.5 * vf[8];
	vf[4] = 0.5 * (1.0 - r * r) * (1.0 + s) - 0.5 * vf[8];
	vf[3] = 0.25 * (1.0 + r) * (1.0 - s) - 0.5 * vf[6] - 0.5 * vf[7] - 0.25 * vf[8];
	vf[2] = 0.25 * (1.0 - r) * (1.0 - s) - 0.5 * vf[5] - 0.5 * vf[6] - 0.25 * vf[8];
	vf[1] = 0.25 * (1.0 - r) * (1.0 + s) - 0.5 * vf[4] - 0.5 * vf[5] - 0.25 * vf[8];
	vf[0] = 0.25 * (1.0 + r) * (1.0 + s) - 0.5 * vf[4] - 0.5 * vf[7] - 0.25 * vf[8];
	return ok = 1;
}

/***************************************************************************
   ROCKFLOW - Funktion: MGradOmega2D_9N
   Aufgabe:
           Berechnet Gradient eines 2D Vektorfeldes,
           dessen Ansatzfunktion (Vektor) bekannt ist (r,s).

   Formalparameter:
           Z: *vf - 2x9 Feld
           E: r,s (s.o.)
   Ergebnis:
           2x9 Matrix
   Aenderungen/Korrekturen:
   12/1999     RK        Erste Version

 **************************************************************************/

int MGradOmega2D_9N(double* vf, double r, double s)
{
	int ok = 0;
#ifdef ERROR_CONTROL
	if (vf == NULL)
	{
		ok = 0;
		return ok;
	}
#endif

	vf[8] = -2.0 * r * (1.0 - s * s);
	vf[7] = +0.5 * (1.0 - s * s) - 0.5 * vf[8];
	vf[6] = -1.0 * r * (1.0 - s) - 0.5 * vf[8];
	vf[5] = -0.5 * (1.0 - s * s) - 0.5 * vf[8];
	vf[4] = -1.0 * r * (1.0 + s) - 0.5 * vf[8];
	vf[3] = +0.25 * (1 - s) - 0.5 * vf[6] - 0.5 * vf[7] - 0.25 * vf[8];
	vf[2] = -0.25 * (1 - s) - 0.5 * vf[5] - 0.5 * vf[6] - 0.25 * vf[8];
	vf[1] = -0.25 * (1 + s) - 0.5 * vf[4] - 0.5 * vf[5] - 0.25 * vf[8];
	vf[0] = +0.25 * (1 + s) - 0.5 * vf[4] - 0.5 * vf[7] - 0.25 * vf[8];

	vf[17] = -2.0 * s * (1.0 - r * r);
	vf[16] = -1.0 * s * (1.0 + r) - 0.5 * vf[17];
	vf[15] = -0.5 * (1.0 - r * r) - 0.5 * vf[17];
	vf[14] = -1.0 * s * (1.0 - r) - 0.5 * vf[17];
	vf[13] = +0.5 * (1 - r * r) - 0.5 * vf[17];
	vf[12] = -0.25 * (1 + r) - 0.5 * vf[15] - 0.5 * vf[16] - 0.25 * vf[17];
	vf[11] = -0.25 * (1 - r) - 0.5 * vf[14] - 0.5 * vf[15] - 0.25 * vf[17];
	vf[10] = +0.25 * (1 - r) - 0.5 * vf[13] - 0.5 * vf[14] - 0.25 * vf[17];
	vf[9] = +0.25 * (1 + r) - 0.5 * vf[13] - 0.5 * vf[16] - 0.25 * vf[17];

	return ok = 1;
}

/***************************************************************************
   ROCKFLOW - Funktion: MGradPhi2D_9N
   Aufgabe:
           Berechnet Gradient der Testfunktionen (r,s)
           grad phi = grad omega
   Formalparameter:
           Z: *vf - 2x9 Feld
           E: r,s (s.o.)
   Ergebnis:
           2x9 Matrix
   Aenderungen/Korrekturen:
   12/1999     RK        Erste Version

 **************************************************************************/

int MGradPhi2D_9N(double* vf, double r, double s)
{
	return MGradOmega2D_9N(vf, r, s);
}

/*##########################################################################
   Funktionen fuer dreidimensionale 20-Knoten-Elemente
 ######################################################################## */

/***************************************************************************
   ROCKFLOW - Funktion: MPhi3D_20N
   Aufgabe:
           Berechnet Testfunktion (r,s,t).
           (K.- J. Bathe, Finite-Elemente-Methoden)

   Formalparameter:
           Z: *vf - 1x20 Feld
           E: r,s,t (s.o.)
   Ergebnis:
           Vektor
   Aenderungen/Korrekturen:
   10/2001     RK        Erste Version

 **************************************************************************/
int MPhi3D_20N(double* vf, double r, double s, double t)
{
	int ok = 0;
	double g[20];

#ifdef ERROR_CONTROL
	if (vf == NULL)
	{
		ok = 0;
		return ok;
	}
#endif

	g[0] = 0.125 * (1.0 + r) * (1.0 + s) * (1.0 + t);
	g[1] = 0.125 * (1.0 - r) * (1.0 + s) * (1.0 + t);
	g[2] = 0.125 * (1.0 - r) * (1.0 - s) * (1.0 + t);
	g[3] = 0.125 * (1.0 + r) * (1.0 - s) * (1.0 + t);
	g[4] = 0.125 * (1.0 + r) * (1.0 + s) * (1.0 - t);
	g[5] = 0.125 * (1.0 - r) * (1.0 + s) * (1.0 - t);
	g[6] = 0.125 * (1.0 - r) * (1.0 - s) * (1.0 - t);
	g[7] = 0.125 * (1.0 + r) * (1.0 - s) * (1.0 - t);
	g[8] = 0.25 * (1.0 - r * r) * (1.0 + s) * (1.0 + t);
	g[9] = 0.25 * (1.0 - r) * (1.0 - s * s) * (1.0 + t);
	g[10] = 0.25 * (1.0 - r * r) * (1.0 - s) * (1.0 + t);
	g[11] = 0.25 * (1.0 + r) * (1.0 - s * s) * (1.0 + t);
	g[12] = 0.25 * (1.0 - r * r) * (1.0 + s) * (1.0 - t);
	g[13] = 0.25 * (1.0 - r) * (1.0 - s * s) * (1.0 - t);
	g[14] = 0.25 * (1.0 - r * r) * (1.0 - s) * (1.0 - t);
	g[15] = 0.25 * (1.0 + r) * (1.0 - s * s) * (1.0 - t);
	g[16] = 0.25 * (1.0 + r) * (1.0 + s) * (1.0 - t * t);
	g[17] = 0.25 * (1.0 - r) * (1.0 + s) * (1.0 - t * t);
	g[18] = 0.25 * (1.0 - r) * (1.0 - s) * (1.0 - t * t);
	g[19] = 0.25 * (1.0 + r) * (1.0 - s) * (1.0 - t * t);

	vf[0] = g[0] - (g[8] + g[11] + g[16]) * 0.5;
	vf[1] = g[1] - (g[8] + g[9] + g[17]) * 0.5;
	vf[2] = g[2] - (g[9] + g[10] + g[18]) * 0.5;
	vf[3] = g[3] - (g[10] + g[11] + g[19]) * 0.5;
	vf[4] = g[4] - (g[12] + g[15] + g[16]) * 0.5;
	vf[5] = g[5] - (g[12] + g[13] + g[17]) * 0.5;
	vf[6] = g[6] - (g[13] + g[14] + g[18]) * 0.5;
	vf[7] = g[7] - (g[14] + g[15] + g[19]) * 0.5;

	vf[8] = g[8];
	vf[9] = g[9];
	vf[10] = g[10];
	vf[11] = g[11];
	vf[12] = g[12];
	vf[13] = g[13];
	vf[14] = g[14];
	vf[15] = g[15];
	vf[16] = g[16];
	vf[17] = g[17];
	vf[18] = g[18];
	vf[19] = g[19];

	return ok = 1;
}

/***************************************************************************
   ROCKFLOW - Funktion: MGradOmega3D_20N
   Aufgabe:
           Berechnet Gradient der Ansatzfunktion (r,s,t)

   Formalparameter:
           Z: *vf - 3x20 Feld
           E: r,s,t (s.o.)
   Ergebnis:
           Vektor
   Aenderungen/Korrekturen:
   10/2001     RK        Erste Version

 **************************************************************************/
int MGradOmega3D_20N(double* vf, double r, double s, double t)
{
	int ok = 0;
	double g_r[20]; /* r-Ableitung */
	double g_s[20]; /* s-Ableitung */
	double g_t[20]; /* t-Ableitung */

#ifdef ERROR_CONTROL
	if (vf == NULL)
	{
		ok = 0;
		return ok;
	}
#endif

	g_r[0] = +0.125 * (1.0 + s) * (1.0 + t);
	g_r[1] = -0.125 * (1.0 + s) * (1.0 + t);
	g_r[2] = -0.125 * (1.0 - s) * (1.0 + t);
	g_r[3] = +0.125 * (1.0 - s) * (1.0 + t);
	g_r[4] = +0.125 * (1.0 + s) * (1.0 - t);
	g_r[5] = -0.125 * (1.0 + s) * (1.0 - t);
	g_r[6] = -0.125 * (1.0 - s) * (1.0 - t);
	g_r[7] = +0.125 * (1.0 - s) * (1.0 - t);
	g_r[8] = -0.5 * r * (1.0 + s) * (1.0 + t);
	g_r[9] = -0.25 * (1.0 - s * s) * (1.0 + t);
	g_r[10] = -0.5 * r * (1.0 - s) * (1.0 + t);
	g_r[11] = +0.25 * (1.0 - s * s) * (1.0 + t);
	g_r[12] = -0.5 * r * (1.0 + s) * (1.0 - t);
	g_r[13] = -0.25 * (1.0 - s * s) * (1.0 - t);
	g_r[14] = -0.5 * r * (1.0 - s) * (1.0 - t);
	g_r[15] = +0.25 * (1.0 - s * s) * (1.0 - t);
	g_r[16] = +0.25 * (1.0 + s) * (1.0 - t * t);
	g_r[17] = -0.25 * (1.0 + s) * (1.0 - t * t);
	g_r[18] = -0.25 * (1.0 - s) * (1.0 - t * t);
	g_r[19] = +0.25 * (1.0 - s) * (1.0 - t * t);

	g_s[0] = +0.125 * (1.0 + r) * (1.0 + t);
	g_s[1] = +0.125 * (1.0 - r) * (1.0 + t);
	g_s[2] = -0.125 * (1.0 - r) * (1.0 + t);
	g_s[3] = -0.125 * (1.0 + r) * (1.0 + t);
	g_s[4] = +0.125 * (1.0 + r) * (1.0 - t);
	g_s[5] = +0.125 * (1.0 - r) * (1.0 - t);
	g_s[6] = -0.125 * (1.0 - r) * (1.0 - t);
	g_s[7] = -0.125 * (1.0 + r) * (1.0 - t);
	g_s[8] = +0.25 * (1.0 - r * r) * (1.0 + t);
	g_s[9] = -0.5 * (1.0 - r) * s * (1.0 + t);
	g_s[10] = -0.25 * (1.0 - r * r) * (1.0 + t);
	g_s[11] = -0.5 * (1.0 + r) * s * (1.0 + t);
	g_s[12] = +0.25 * (1.0 - r * r) * (1.0 - t);
	g_s[13] = -0.5 * (1.0 - r) * s * (1.0 - t);
	g_s[14] = -0.25 * (1.0 - r * r) * (1.0 - t);
	g_s[15] = -0.5 * (1.0 + r) * s * (1.0 - t);
	g_s[16] = +0.25 * (1.0 + r) * (1.0 - t * t);
	g_s[17] = +0.25 * (1.0 - r) * (1.0 - t * t);
	g_s[18] = -0.25 * (1.0 - r) * (1.0 - t * t);
	g_s[19] = -0.25 * (1.0 + r) * (1.0 - t * t);

	g_t[0] = +0.125 * (1.0 + r) * (1.0 + s);
	g_t[1] = +0.125 * (1.0 - r) * (1.0 + s);
	g_t[2] = +0.125 * (1.0 - r) * (1.0 - s);
	g_t[3] = +0.125 * (1.0 + r) * (1.0 - s);
	g_t[4] = -0.125 * (1.0 + r) * (1.0 + s);
	g_t[5] = -0.125 * (1.0 - r) * (1.0 + s);
	g_t[6] = -0.125 * (1.0 - r) * (1.0 - s);
	g_t[7] = -0.125 * (1.0 + r) * (1.0 - s);
	g_t[8] = +0.25 * (1.0 - r * r) * (1.0 + s);
	g_t[9] = +0.25 * (1.0 - r) * (1.0 - s * s);
	g_t[10] = +0.25 * (1.0 - r * r) * (1.0 - s);
	g_t[11] = +0.25 * (1.0 + r) * (1.0 - s * s);
	g_t[12] = -0.25 * (1.0 - r * r) * (1.0 + s);
	g_t[13] = -0.25 * (1.0 - r) * (1.0 - s * s);
	g_t[14] = -0.25 * (1.0 - r * r) * (1.0 - s);
	g_t[15] = -0.25 * (1.0 + r) * (1.0 - s * s);
	g_t[16] = -0.5 * (1.0 + r) * (1.0 + s) * t;
	g_t[17] = -0.5 * (1.0 - r) * (1.0 + s) * t;
	g_t[18] = -0.5 * (1.0 - r) * (1.0 - s) * t;
	g_t[19] = -0.5 * (1.0 + r) * (1.0 - s) * t;

	vf[0] = g_r[0] - (g_r[8] + g_r[11] + g_r[16]) * 0.5;
	vf[1] = g_r[1] - (g_r[8] + g_r[9] + g_r[17]) * 0.5;
	vf[2] = g_r[2] - (g_r[9] + g_r[10] + g_r[18]) * 0.5;
	vf[3] = g_r[3] - (g_r[10] + g_r[11] + g_r[19]) * 0.5;
	vf[4] = g_r[4] - (g_r[12] + g_r[15] + g_r[16]) * 0.5;
	vf[5] = g_r[5] - (g_r[12] + g_r[13] + g_r[17]) * 0.5;
	vf[6] = g_r[6] - (g_r[13] + g_r[14] + g_r[18]) * 0.5;
	vf[7] = g_r[7] - (g_r[14] + g_r[15] + g_r[19]) * 0.5;

	vf[8] = g_r[8];
	vf[9] = g_r[9];
	vf[10] = g_r[10];
	vf[11] = g_r[11];
	vf[12] = g_r[12];
	vf[13] = g_r[13];
	vf[14] = g_r[14];
	vf[15] = g_r[15];
	vf[16] = g_r[16];
	vf[17] = g_r[17];
	vf[18] = g_r[18];
	vf[19] = g_r[19];

	vf[20] = g_s[0] - (g_s[8] + g_s[11] + g_s[16]) * 0.5;
	vf[21] = g_s[1] - (g_s[8] + g_s[9] + g_s[17]) * 0.5;
	vf[22] = g_s[2] - (g_s[9] + g_s[10] + g_s[18]) * 0.5;
	vf[23] = g_s[3] - (g_s[10] + g_s[11] + g_s[19]) * 0.5;
	vf[24] = g_s[4] - (g_s[12] + g_s[15] + g_s[16]) * 0.5;
	vf[25] = g_s[5] - (g_s[12] + g_s[13] + g_s[17]) * 0.5;
	vf[26] = g_s[6] - (g_s[13] + g_s[14] + g_s[18]) * 0.5;
	vf[27] = g_s[7] - (g_s[14] + g_s[15] + g_s[19]) * 0.5;

	vf[28] = g_s[8];
	vf[29] = g_s[9];
	vf[30] = g_s[10];
	vf[31] = g_s[11];
	vf[32] = g_s[12];
	vf[33] = g_s[13];
	vf[34] = g_s[14];
	vf[35] = g_s[15];
	vf[36] = g_s[16];
	vf[37] = g_s[17];
	vf[38] = g_s[18];
	vf[39] = g_s[19];

	vf[40] = g_t[0] - (g_t[8] + g_t[11] + g_t[16]) * 0.5;
	vf[41] = g_t[1] - (g_t[8] + g_t[9] + g_t[17]) * 0.5;
	vf[42] = g_t[2] - (g_t[9] + g_t[10] + g_t[18]) * 0.5;
	vf[43] = g_t[3] - (g_t[10] + g_t[11] + g_t[19]) * 0.5;
	vf[44] = g_t[4] - (g_t[12] + g_t[15] + g_t[16]) * 0.5;
	vf[45] = g_t[5] - (g_t[12] + g_t[13] + g_t[17]) * 0.5;
	vf[46] = g_t[6] - (g_t[13] + g_t[14] + g_t[18]) * 0.5;
	vf[47] = g_t[7] - (g_t[14] + g_t[15] + g_t[19]) * 0.5;

	vf[48] = g_t[8];
	vf[49] = g_t[9];
	vf[50] = g_t[10];
	vf[51] = g_t[11];
	vf[52] = g_t[12];
	vf[53] = g_t[13];
	vf[54] = g_t[14];
	vf[55] = g_t[15];
	vf[56] = g_t[16];
	vf[57] = g_t[17];
	vf[58] = g_t[18];
	vf[59] = g_t[19];

	return ok = 1;
}

/***************************************************************************
   ROCKFLOW - Funktion: MGradPhi3D_20N
   Aufgabe:
           Berechnet Gradient der Testfunktionen (r,s,t)
           grad phi = grad omega
   Formalparameter:
           Z: *vf - 3x20 Feld
           E: r,s (s.o.)
   Ergebnis:
           3x20 Matrix
   Aenderungen/Korrekturen:
   10/2001     RK        Erste Version

 **************************************************************************/
int MGradPhi3D_20N(double* vf, double r, double s, double t)
{
	return MGradOmega3D_20N(vf, r, s, t);
}

/*##########################################################################
   Sortierfunktionen
 ########################################################################*/
/***************************************************************************
   ROCKFLOW - Funktion: MCompare_for_MQSort_LongDouble
   Aufgabe:
           Vergleichfunktion fuer zugehoerige Sortierfunktion
   Formalparameter:
   *Arg1, *Arg2

   Ergebnis:
           1,0,-1  (groesser, gleich, kleiner)
   Aenderungen/Korrekturen:
   06/2001     MK       Erste Version

 **************************************************************************/
int MCompare_for_MQSort_LongDouble(const void* Arg1, const void* Arg2)
{
	typedef struct
	{
		long index;
		double value;
	} MType_LongDouble;
	MType_LongDouble *arg1, *arg2;
	arg1 = (MType_LongDouble*)Arg1;
	arg2 = (MType_LongDouble*)Arg2;
	if ((arg1->value) > (arg2->value))
		return 1;
	if (fabs(arg1->value - arg2->value) < MKleinsteZahl)
		return 0;
	else
		return -1;
}

/***************************************************************************
   ROCKFLOW - Funktion: MQSort_LongDouble
   Aufgabe:
           Sortierfunktion
   Formalparameter:
   *DataSets
           NumberOfDataSets
           SizeOfDataSet

   Ergebnis:
           Sortierter Datensatz
   Aenderungen/Korrekturen:
   06/2001     MK       Erste Version

 **************************************************************************/
void MQSort_LongDouble(void* DataSets, const int NumberOfDataSets, const int SizeOfDataSet)
{
	qsort(DataSets, NumberOfDataSets, SizeOfDataSet, MCompare_for_MQSort_LongDouble);
}

/*##########################################################################
   eventuell noch nuetzliche Funktionen,
 ######################################################################## */

/***************************************************************************
   ROCKFLOW - Funktion: JWDMMultMatSkalar

   Aufgabe:
           Multiplikation einer Matrix mit einem Skalarwert,
           wobei die Matrix nicht veraendert wird
   Formalparameter:
           E: *matrix - Zeiger auf Matrix
           E: skal - Skalarwert
           E: m,n - Ausdehnung
   Ergebnis:
   Zeiger auf Ergebnis
   Aenderungen/Korrekturen:
   07/1994     hh        Erste Version

 **************************************************************************/

double* JWDMMultMatSkalar(double* matrix, double skal, long m, long n)
{ /* flip ist Zwischenspeicher fuer die Matrix */
	double *flip, *flippy;
	long i;
	flip = (double*)Malloc(sizeof(double) * n * m);
	flippy = flip;
	if (flip == NULL)
	{
		DisplayErrorMsg(" kein Speicher mehr !!! ");
		return 0;
	} /* if */
	for (i = 0; i < m * n; i++)
		flip[i] = matrix[i] * skal;
	/* for */
	return flippy;
	/* flippy = Free(flippy); */
} /* extern double *JWDMMultMatSkalar */

/**************************************************************************
   ROCKFLOW - Function: GetPriMatFromTriMat

   Task:
      Gets a 9x9 Matrix from 3x3 matrix according to:

                 a b c  a b c
   a b c         d e f  d e f
   d e f   ==>   g h i  g h i
   g h i         a b c  a b c
                 d e f  d e f
   g h i  g h i

   Parameter: (I: Input; R: Return; X: Both)
   I: double *mat

   Return:
  *mat2

   Programming:
   07/2003   MB   First Version
**************************************************************************/

int GetPriMatFromTriMat(double* mat1, double* mat_2)
{
	mat_2[0] = mat_2[3] = mat1[0];
	mat_2[1] = mat_2[4] = mat1[1];
	mat_2[2] = mat_2[5] = mat1[2];
	mat_2[6] = mat_2[9] = mat1[3];
	mat_2[7] = mat_2[10] = mat1[4];
	mat_2[8] = mat_2[11] = mat1[5];
	mat_2[12] = mat_2[15] = mat1[6];
	mat_2[13] = mat_2[16] = mat1[7];
	mat_2[14] = mat_2[17] = mat1[8];

	mat_2[18] = mat_2[21] = mat1[0];
	mat_2[19] = mat_2[22] = mat1[1];
	mat_2[20] = mat_2[23] = mat1[2];
	mat_2[24] = mat_2[27] = mat1[3];
	mat_2[25] = mat_2[28] = mat1[4];
	mat_2[26] = mat_2[29] = mat1[5];
	mat_2[30] = mat_2[33] = mat1[6];
	mat_2[31] = mat_2[34] = mat1[7];
	mat_2[32] = mat_2[35] = mat1[8];

	return 1;
}

/**************************************************************************
   ROCKFLOW - Funktion: MMultMatMat2
   Aufgabe:
           Multiplikation zweier Matrizen der gleichen Grösse nach:

   a b  x  e f  ==>  axe bxf
   c d     g h       cxg dxh

   Formalparameter:
           E: *mat1, *mat2, m1, n1
   Ergebnis:
  *matrix

**************************************************************************/
int MMultMatMat2(double* mat1, long m1, long n1, double* mat2, double* ergebnis)
{
	int i;

	for (i = 0; i < m1 * n1; i++)
		ergebnis[i] = mat1[i] * mat2[i];

	return 1;
}

/**************************************************************************
   ROCKFLOW - Function: TensorDrehDich

   Task: Dreht Tensor von von stromlinienorietierten in lokale physikalische
         (Element) Koordinaten.

   Parameter: (I: Input; R: Return; X: Both)
           I: double*  d Tensor
              double*  velo Geschwindigkeitsvektor

   Return:
   double* d

   Programming:
   08/2003   MB   First Version, herausgelöst aus CalcEle3D

**************************************************************************/
double* TensorDrehDich(double* d, double* velo)

{
	double zwa[18];
	double zwi[18];
	long l, ii;
	double vg;
	double fkt;
	double zwo[9];

	vg = MBtrgVec(velo, 3);

	if (d[0] > MKleinsteZahl || d[4] > MKleinsteZahl || d[8] > MKleinsteZahl)
		/* Drehen: Stromrichtung - r,s,t */
		if (vg > MKleinsteZahl)
		{
			/* 1. Zeile */
			for (l = 0; l < 3; l++)
				zwa[l] = velo[l];
			MNormiere(zwa, 3);
			/* 2. Zeile */
			fkt = fabs(velo[0]);
			ii = 0;
			for (l = 1; l < 3; l++)
				if (fabs(velo[l]) < fkt)
				{
					fkt = fabs(velo[l]);
					ii = l;
				}
			zwo[ii] = 0.0;
			zwo[(ii + 1) % 3] = velo[(ii + 2) % 3];
			zwo[(ii + 2) % 3] = -velo[(ii + 1) % 3];
			MNormiere(zwo, 3);
			/* 3. Zeile */
			M3KreuzProdukt(zwa, zwo, zwi);
			MNormiere(zwi, 3);
			for (l = 0; l < 3; l++)
			{
				zwa[3 + l] = zwo[l];
				zwa[6 + l] = zwi[l];
			}
			/* dreh^T * D * dreh */
			MMultMatMat(d, 3, 3, zwa, 3, 3, zwi, 3, 3);
			MTranspoMat(zwa, 3, 3, zwo);
			MMultMatMat(zwo, 3, 3, zwi, 3, 3, d, 3, 3);
		} /* end if (vg > MKleinsteZahl) */
	/* end if (d[0] > MKleinsteZahl || d[4] > MKleinsteZahl || d[8] > MKleinsteZahl) */

	return d;
}
#endif //#ifdef obsolete  //05.03.2010 WW

///////////////////////////////////////////////////////////////////
//
//        Moved here by WW. 05.03.2010
//
///////////////////////////////////////////////////////////////////

/***************************************************************************
   ROCKFLOW - Funktion: MCalcDistancePointToLine

   Aufgabe:
           Abstand Punkt - Linie im R3

   Formalparameter:
           E: *pt  - Punktkoordinaten
           E: *l1  - Punktkoordinaten fuer ersten Punkt auf Linie
           E: *l2  - Punktkoordinaten fuer zweiten Punkt auf Linie

   Ergebnis:
   Abstand

   Aenderungen/Korrekturen:
   02/2000   C. Thorenz      Erste Version

 **************************************************************************/
double MCalcDistancePointToLine(double* pt, double* l1, double* l2)
{
	int i;
	double vec1[3], vec2[3], erg_vec[3], l;

	/* Vektor Linienpunkt zu Punkt */
	for (i = 0; i < 3; i++)
		vec1[i] = pt[i] - l1[i];

	/* Vektor Linienpunkt zu Linienpunkt */
	for (i = 0; i < 3; i++)
		vec2[i] = l2[i] - l1[i];

	/* Normieren */
	l = MBtrgVec(vec2, 3);
#ifdef ERROR_CONTROL
	if (l < MKleinsteZahl)
	{
		printf("\n FEHLER IN CalcDistancePointToLine: Laenge ist Null !!!");
		exit(1);
	}
#endif

	for (i = 0; i < 3; i++)
		vec2[i] /= (l + MKleinsteZahl);

	M3KreuzProdukt(vec1, vec2, erg_vec);

	return MBtrgVec(erg_vec, 3);
}

/***************************************************************************
   ROCKFLOW - Funktion:  MCalcProjectionOfPointOnLine

   Aufgabe:
           Ermittelt den Projektionspunkt eines Punktes auf eine Linie (Lotpunkt)

   Formalparameter:
           E: *pt   - Punktkoordinaten
           E: *l1   - Punktkoordinaten fuer ersten Punkt auf Linie
           E: *l2   - Punktkoordinaten fuer zweiten Punkt auf Linie
           R: *proj   - Punktkoordinaten fuer Projektionspunkt auf Linie
   Ergebnis:
   Abstand Punkt-Linie

   Aenderungen/Korrekturen:
   02/2000   C. Thorenz      Erste Version

 **************************************************************************/
double MCalcProjectionOfPointOnLine(double* pt, double* l1, double* l2, double* proj)
{
	int i;
	double vec1[3], vec2[3], l, projektionslaenge, abstand;

	/* Vektor Linienpunkt zu Punkt */
	for (i = 0; i < 3; i++)
		vec1[i] = pt[i] - l1[i];

	/* Vektor Linienpunkt zu Linienpunkt */
	for (i = 0; i < 3; i++)
		vec2[i] = l2[i] - l1[i];

	/* Normieren */
	l = MBtrgVec(vec2, 3);
#ifdef ERROR_CONTROL
	if (l < MKleinsteZahl)
	{
		printf("\n Fehler in MCalcProjectionOfPointOnLine: Laenge ist Null !!!");
		exit(1);
	}
#endif
	for (i = 0; i < 3; i++)
		vec2[i] /= (l + MKleinsteZahl);

	projektionslaenge = MSkalarprodukt(vec1, vec2, 3);

	/* Punkt bestimmen */
	for (i = 0; i < 3; i++)
		proj[i] = l1[i] + projektionslaenge * vec2[i];

	abstand = MCalcDistancePointToPoint(proj, pt);

	return abstand;
}

/***************************************************************************
   ROCKFLOW - Funktion: MCalcDistancePointToPlane

   Aufgabe:
           Abstand Punkt - Ebene im R3

   Formalparameter:
           E: *pt  - Punktkoordinaten
           E: *e1  - Punktkoordinaten fuer ersten Punkt auf Ebene
           E: *e2  - Punktkoordinaten fuer zweiten Punkt auf Ebene
           E: *e3  - Punktkoordinaten fuer dritten Punkt auf Ebene

   Ergebnis:
   Abstand-Flaeche   (Kann negativ sein, je nach Lage des Punkts zur Ebene)

   Aenderungen/Korrekturen:
   02/2000   C. Thorenz      Erste Version

 **************************************************************************/
double MCalcDistancePointToPlane(double const* const pt, double* e1, double* e2, double* e3)
{
	int i;
	double vec1[3], vec2[3], vec3[3], normal[3], volume, area;

	/* 1. Ebenen-Vektor  */
	for (i = 0; i < 3; i++)
		vec1[i] = e2[i] - e1[i];

	/* 2. Ebenen-Vektor  */
	for (i = 0; i < 3; i++)
		vec2[i] = e3[i] - e1[i];

	/* Ebene-Punkt-Vektor  */
	for (i = 0; i < 3; i++)
		vec3[i] = pt[i] - e1[i];

	/* Normalenvektor */
	M3KreuzProdukt(vec1, vec2, normal);

	/* Volumen des Spats */
	volume = MSkalarprodukt(normal, vec3, 3);

	/* Flaeche auf der Ebene */
	area = MBtrgVec(normal, 3);

#ifdef ERROR_CONTROL
	if (area < MKleinsteZahl)
	{
		printf("\n FEHLER IN CalcDistancePointToPlane: Flaeche ist Null !!!");
		exit(1);
	}
#endif

	return volume / (area + MKleinsteZahl);
}

/**************************************************************************/
/* ROCKFLOW - Funktion: MMin
 */
/* Aufgabe:
   Gibt den Kleineren der Eingabewerte zurueck
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E double a,b
 */
/* Ergebnis:
   kleinere Zahl
 */
/* Programmaenderungen:
   2/2000     C.Thorenz  Erste Version                                                                          */
/**************************************************************************/
double MMin(double a, double b)
{
	if (a < b)
		return a;
	else
		return b;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: MMax
 */
/* Aufgabe:
   Gibt den Groesseren der Eingabewerte zurueck
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E double a,b
 */
/* Ergebnis:
   groessere Zahl
 */
/* Programmaenderungen:
   2/2000     C.Thorenz  Erste Version                                                                          */
/**************************************************************************/
double MMax(double a, double b)
{
	if (a > b)
		return a;
	else
		return b;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: MRange
 */
/* Aufgabe:
   Begrenzt eine Zahl auf ein Intervall
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E double a,b,c
 */
/* Ergebnis:
   Zahl im Intervall
 */
/* Programmaenderungen:
   7/2000     C.Thorenz  Erste Version                                                                          */
/**************************************************************************/
double MRange(double a, double b, double c)
{
	if (b < a)
		return a;
	if (b > c)
		return c;

	return b;
}

/*##########################################################################
   Funktionen fuer Vektoren und Matrizen
 ######################################################################## */

/***************************************************************************
   ROCKFLOW - Funktion: MNulleVec
   Aufgabe:
           Setze angegebenen Vektor = 0.0
   Formalparameter:
           E: *vec
           E: g
   Ergebnis:
           viele Nullen
   Aenderungen/Korrekturen:
   09/1994     hh        Erste Version

 **************************************************************************/

void MNulleVec(double* vec, long g)
{
	register long i;
	for (i = 0; i < g; i++)
		vec[i] = 0.0;
}

/***************************************************************************
   ROCKFLOW - Funktion: MNulleMat
   Aufgabe:
           Setze angegebene Matrix = 0.0
   Formalparameter:[B
           E: *mat
           E: g
   Ergebnis:
           viele Nullen
   Aenderungen/Korrekturen:
   09/1994     hh        Erste Version

 **************************************************************************/

void MNulleMat(double* mat, long m, long n)
{
	register long i;
	for (i = 0; i < m * n; i++)
		mat[i] = 0.0;
}

#ifndef NEW_EQS // WW. 05.03.2010
/*##########################################################################
   Funktionen fuer Gleichungsloeser (CG)
 ######################################################################## */

/**************************************************************************
   ROCKFLOW - Funktion: MVekNorm1

   Aufgabe:
   Berechnet die Spaltensummennorm von x

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E double *x : Zeiger auf Vektor
   E long n : Dimension von x

   Ergebnis:
   Norm

   Programmaenderungen:
   11/1995     MSR        Erste Version

**************************************************************************/
double MVekNorm1(double* x, long n)
{
	register long i;
	register double erg = 0.0;
#ifdef SX
#pragma cdir nodep
#endif
	for (i = 0l; i < n; i++)
		erg += fabs(x[i]);
	return erg;
}

/**************************************************************************
   ROCKFLOW - Funktion: MVekNorm2

   Aufgabe:
   Berechnet die euklidische Norm von x

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E double *x : Zeiger auf Vektor
   E long n : Dimension von x

   Ergebnis:
   Norm

   Programmaenderungen:
   11/1995     MSR        Erste Version

**************************************************************************/
double MVekNorm2(double* x, long n)
{
	register long i;
	register double erg = 0.0;
#ifdef SX
#pragma cdir nodep
#endif
	for (i = 0l; i < n; i++)
		erg += x[i] * x[i];
	return sqrt(erg);
}

/**************************************************************************
   ROCKFLOW - Funktion: MVekSum

   Aufgabe:
   Fuehrt die Operation

   x = x + alpha * y

   durch

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   X double *x : Zeiger auf Vektor x
   E double alpha : Faktor alpha
   E double *y : Zeiger auf Vektor y
   E long n : Dimension von x und y

   Ergebnis:
   - void -

   Programmaenderungen:
   11/1995     MSR        Erste Version

**************************************************************************/
#ifdef SX
void MVekSum(double* restrict x, double alpha, double* restrict y, long n)
#else
void MVekSum(double* x, double alpha, double* y, long n)
#endif
{
	register long i;
#ifdef SX
#pragma cdir nodep
#endif
	for (i = 0l; i < n; i++)
		x[i] += alpha * y[i];
}

/**************************************************************************
   ROCKFLOW - Funktion: MVekGle

   Aufgabe:
   Fuehrt die Operation

   z = alpha * x + beta * y

   durch

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E double alpha : Faktor alpha
   E double *x : Zeiger auf Vektor x
   E double beta : Faktor beta
   E double *y : Zeiger auf Vektor y
   R double *z : Zeiger auf Ergabnisvektor z
   (der Speicher muss bereits allokiert sein)
   E long n : Dimension von x, y und z

   Ergebnis:
   - void -

   Programmaenderungen:
   11/1995     MSR        Erste Version

**************************************************************************/
void MVekGle(double alpha, double* x, double beta, double* y, double* z, long n)
{
	register long i;
	for (i = 0l; i < n; i++)
		z[i] = alpha * x[i] + beta * y[i];
}

/**************************************************************************
   ROCKFLOW - Funktion: MVekDist

   Aufgabe:
   Berechnet den Abstand zwischen x und y in einer gegebenen Norm

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E double *x : Zeiger auf Vektor
   E double *y : Zeiger auf Vektor
   E long n : Dimension von x und y

   Ergebnis:
   Abstand

   Programmaenderungen:
   09/1997     AH         Einbau der Funktion

**************************************************************************/
double MVekDist(double* x, double* y, long n)
{
	double dist;
	double* d;

	d = (double*)Malloc(n * sizeof(double));
	MVekGle(1., x, -1., y, d, n);
	dist = VEKNORM_BICG(d, n);
	d = (double*)Free(d);

	return dist;
}

/**************************************************************************
   ROCKFLOW - Funktion: MMachVec
   Aufgabe:
           Erzeugen eines neuen Vektors
   Formalparameter:
           E: g
   Ergebnis:
           Zeiger auf den erzeugten Vektor.
           oder Rueckgabewert NULL wenn nicht genuegend Speicher da ist.
   Aenderungen/Korrekturen:
   08/1994     hh        Erste Version
   11/1999     C.Thorenz Nullung herausgenommen
**************************************************************************/

double* MMachVec(long g)
{
	double* zwerg;
	zwerg = (double*)Malloc(sizeof(double) * g);
#ifdef ERROR_CONTROL
	if (zwerg == NULL)
		DisplayErrorMsg("zuwenig Speicher in MMachVec");
#endif
	return zwerg;
} /* MMachVec */

/**************************************************************************
   ROCKFLOW - Funktion: MNullVec
   Aufgabe:
           Schreibt Nullen in einen Vektors
   Formalparameter:
           E: *zwerg
           E: g
   Ergebnis:
   Aenderungen/Korrekturen:
   11/1999     C.Thorenz Erste Version
**************************************************************************/
#ifdef SX
void MNullVec(double* restrict zwerg, long g)
#else
void MNullVec(double* zwerg, long g)
#endif
{
	register long i;
	zwerg = (double*)Malloc(sizeof(double) * g);
#ifdef ERROR_CONTROL
	if (zwerg == NULL)
		DisplayErrorMsg("Fehler in MNullVec");
#endif
#ifdef SX
#pragma cdir nodep
#endif
	for (i = 0; i < g; i++)
		zwerg[i] = 0.0;
} /* MNullVec */

/**************************************************************************
   ROCKFLOW - Funktion: MKopierVec
   Aufgabe:
           Kopieren eines Vektors
   Formalparameter:
           E: *vecquelle
           A: *vecziel
           E: g
   Ergebnis:
           Rueckgabewert 1
   Aenderungen/Korrekturen:
   08/1994     hh        Erste Version
   11/1995     msr       register
   11/1999     C.Thorenz Fehlerabfrage

**************************************************************************/
#ifdef SX
void MKopierVec(double* restrict vecquelle, double* restrict vecziel, long g)
#else
void MKopierVec(double* vecquelle, double* vecziel, long g)
#endif
{
	register long i;
#ifdef ERROR_CONTROL
	if ((vecquelle == NULL) || (vecziel == NULL))
		DisplayErrorMsg("Fehler in MLoeschVec");
#endif
#ifdef SX
#pragma cdir nodep
#endif
	for (i = 0; i < g; i++)
		vecziel[i] = vecquelle[i];
} /* MKopierVec */

/* Function MVectorlength
   Aufgabe:
           Calculates the length of a vector in cartesian co-ordinates.
   Formalparameter:
           delta x
           delta y
           delta z
   Ergebnis:
           Length of vector
   Aenderungen/Korrekturen:
   05/2004     CMCD        Erste Version
   01/2011	TF changed from pow(*,0.5) to sqrt(*)
 */
// double MVectorlength(double dx, double dy, double dz)
//{
////    double length;
////    length=pow(((dz*dz)+(dy*dy)+(dx*dx)),0.5);
////    return length;
//    return sqrt(dz * dz + dy * dy + dx * dx);
//}

/**************************************************************************
   ROCKFLOW - Funktion: MAddSkalVektoren
   Aufgabe:
           Vektoren mit Skalar multiplizieren und dann addieren
   Formalparameter:
           E: *v1, *v2  : Vektoren
           e: m1,  m2   : Multiplikatoren
           A: *vout     : Ergebnisvektor
           E: g         : Vektorlaenge
   Ergebnis:
           Vektorsumme
   Aenderungen/Korrekturen:
   8/2001     C.Thorenz: Erste Version
**************************************************************************/
#ifdef SX
int MAddSkalVektoren(double* restrict v1, double m1, double* restrict v2, double m2, double* restrict vout, long g)
#else
int MAddSkalVektoren(double* v1, double m1, double* v2, double m2, double* vout, long g)
#endif
{
	register long i;
	// WW    if ((m1==1.)&&(m2==1.))
	if ((fabs(m1 - 1.) < MKleinsteZahl) && (fabs(m2 - 1.) < MKleinsteZahl))
#ifdef SX
#pragma cdir nodep
#endif
		for (i = 0; i < g; i++)
			vout[i] = v1[i] + v2[i];
	else if (fabs(m1 - 1.) < MKleinsteZahl)
#ifdef SX
#pragma cdir nodep
#endif
		for (i = 0; i < g; i++)
			vout[i] = v1[i] + m2 * v2[i];
	else if (fabs(m2 - 1.) < MKleinsteZahl)
#ifdef SX
#pragma cdir nodep
#endif
		for (i = 0; i < g; i++)
			vout[i] = m1 * v1[i] + v2[i];
	else
#ifdef SX
#pragma cdir nodep
#endif
		for (i = 0; i < g; i++)
			vout[i] = m1 * v1[i] + m2 * v2[i];

	return 1;
}
#endif //   05.03.2010. WW

////
/**************************************************************************
   ROCKFLOW - Funktion: MMultVecVec
   Aufgabe:
           Multiplikation Vektor mit Vektor (zeilenweise)

                  xxxxxx <- vec1
                x oooooo
        vec1 -> x oooooo <- mato
                x oooooo
                x oooooo

   Formalparameter:
   E: *vec1 *vec2
   E: gv1, gv2
   A: *mato, m (Zeilen) ,n (Spalten)
   Ergebnis:
   Matrix
   Aenderungen/Korrekturen:
   08/1994     hh        Erste Version
   11/1999     C.Thorenz Register-Variablen
**************************************************************************/

int MMultVecVec(double* vec1, long gv1, double* vec2, long gv2, double* mato, long mo, long no)
{
	register long i, j;

#ifdef ERROR_CONTROL
	const bool check_err = true;
#else
	const bool check_err = false;
#endif
	if (check_err)
	{
		if (gv1 != mo)
			DisplayErrorMsg("MMultVecVec: Groesse von Matrix und Vektor 1 passen nicht");
		if (gv2 != no)
			DisplayErrorMsg("MMultVecVec: Groesse von Matrix und Vektor 2 passen nicht");
	}

	mo = gv1;
	no = gv2;
	for (i = 0; i < gv1; i++)
		for (j = 0; j < gv2; j++)
			mato[i * gv2 + j] = vec1[i] * vec2[j];
	return 1;
} /* MMultVecVec */

/**************************************************************************
   ROCKFLOW - Funktion: MMultVecMat
   Aufgabe:
           Multiplikation Vektor mit Matrix

                     xxxxxx
                     xxxxxx
                     xxxxxx <- mat
                     xxxxxx
          vec-> xxxx oooooo <-veco

   xxxxxx      Es wird zeilenweise aufsummiert.
   XXXXXX
   xxxxxx
   xxxxxx
   XXXX Oooooo

   Formalparameter:
   E: *vec *mat
   E: gv, gm
   A: *veco
   Ergebnis:
   Vektor
   Aenderungen/Korrekturen:
   08/1994     hh        Erste Version
   11/1999     C.Thorenz Register-Variablen
**************************************************************************/
#ifdef SX
int MMultVecMat(double* restrict vec, long gv, double* restrict mat, long m, long n, double* restrict veco, long go)
#else
int MMultVecMat(double* vec, long gv, double* mat, long m, long n, double* veco, long go)
#endif
{
	register long i, j;
#ifdef ERROR_CONTROL
	const bool check_err = true;
#else
	const bool check_err = false;
#endif
	if (check_err)
	{
		if (gv != m)
			DisplayErrorMsg("MMultVecMat: Groesse von Matrix und Vektor passen nicht");
		if (go != n)
			DisplayErrorMsg("MMultVecMat: Groesse von Ergebnis-Vektor stimmt nicht");
	}

	MNulleVec(veco, n); /* cb: Nullen nicht vergessen ! */
	for (i = 0; i < m; i++)
#ifdef SX
#pragma cdir nodep
#endif
		for (j = 0; j < n; j++)
			veco[j] += vec[i] * mat[j + i * n];
	return 1;
} /* MMultVecMat */

/**************************************************************************
   ROCKFLOW - Funktion: MMultMatVec
   Aufgabe:
           Multiplikation Matrix Vektor
                  z
                  z
                  z
                  z
            n     z
        x x x x x o
      m x x x x x o
   x x x x x o

   Formalparameter:
   E: *mat
   E: m,n     Groesse von *mat
   E: *vec
   E: g       Groesse von *vec
   A: *veco
   E: r       Groesse von *veco
   Ergebnis:
   Vektor
   Aenderungen/Korrekturen:
   08/1994     hh        Erste Version
   24.06.1999  OK        Warnung
   11/1999     C.Thorenz Warnung in #IFDEF, Register-Variablen
**************************************************************************/

int MMultMatVec(/* Matrix */
                const double* mat, long m, long n,
                /* Veltor fuer das Produkt */
                double* vec, long g,
                /* Vektor fuer das Ergebnis */
                double* veco, long r)
{
	register long i, k;

#ifdef ERROR_CONTROL
	const bool check_err = true;
#else
	const bool check_err = false;
#endif
	if (check_err)
	{
		if (g != n)
			DisplayErrorMsg("MMultMatVec: Groesse von Matrix und Vektor passen nicht");
		if (r != m)
			DisplayErrorMsg("MMultMatVec: Groesse von Ergebnis-Vektor stimmt nicht");
	}

	r = m;
	for (k = 0; k < m; k++)
	{
		veco[k] = 0.0; /* cb: Nullen nicht vergessen ! */
		for (i = 0; i < g; i++)
			veco[k] += vec[i] * mat[i + k * n];
	}
	return 1;
} /* MMultMatVec */

/**************************************************************************
   ROCKFLOW - Funktion: MMultMatMat
   Aufgabe:
           Multiplikation Matrix mit Matrix
   Erlaeuterung:              n2
                            X X X   Die Funktion arbeitet sich zeilenweise
             mat2           x x x   durch die Matrizen durch. Es wird mit
                         m2 x x x   dem ersten Element der ersten Zeile von
                            x x x   mat1 begonnen. Damit werden alle
                            x x x   Multiplikationen in der ersten Zeile von
                    n1     -------  mat2 durchgefuehrt. Das Ergebnis wird
   X X X X X | O O O   jeweils an der entsprechenden Stelle in
   mat1    x x x x x | o o o   der Ergebnismatrix mato hinzuaddiert.
   m1 x x x x x | o o o   Danach geht es mit dem zweiten Element der
   x x x x x | o o o   ersten Zeile von mat1 so weiter.
   Die Elemente der einzelnen Matrizen sind     mato mat1 mat2
   immer zeilenweise durchnummeriert            0  =  0  *  0
   gespeichert/verfuegbar. In der aus-          1  =  0  *  1
   schnittweisen Tabelle sind in der            2  =  0  *  2
   Reihenfolge der Abarbeitung die              0  =  1  *  3
   Elementnummern.                              1  =  1  *  4
   2  =  1  *  5
   0  =  2  *  6

   Kommentar C.Thorenz: Vorsicht! Die zweite Matrix wird vorher trans-
   poniert.

   Formalparameter:
   E: *mat1 *mat2
   A: *mato
   E: m1,n1, m2, n2, mo, no
   Ergebnis:
   Matrix
   Aenderungen/Korrekturen:
   08/1994     hh        Erste Version
   3/1999     C.Thorenz Optimiert (... ca. doppelt so schnell)
   1/2001     C.Thorenz Noch etwas weiter optimiert
   8/2001     C.Thorenz CBLAS eingehaengt

**************************************************************************/
#ifdef SX
int MMultMatMat(double* restrict mat1, long m1, long n1, double* restrict mat2, long m2, long n2, double* restrict mato,
                long mo, long no)
#else
int MMultMatMat(double* mat1, long m1, long n1, double* mat2, long m2, long n2, double* mato, long mo, long no)
#endif
{
#ifdef CBLAS_MMultMatMat
	enum CBLAS_ORDER Order = CblasRowMajor;
	enum CBLAS_TRANSPOSE TransA = CblasNoTrans;
	enum CBLAS_TRANSPOSE TransB = CblasNoTrans;
#else
	register int i, j, ih1, ih2;
#endif

#ifdef ERROR_CONTROL
	if ((m1 != mo) || (n2 != no) || (n1 != m2))
	{
		DisplayErrorMsg("MMultMatMat:Die Matrizen passen nicht zueinander !");
		return 0;
	}
#endif

#ifndef CBLAS_MMultMatMat
#ifdef SX
#pragma cdir nodep
#endif
	for (i = 0; i < mo * no; i++)
		mato[i] = 0.0; /*Ergebnismatrix vorsichtshalber nullen */

	for (i = 0; i < n1 * m1; i++)
	{
		ih1 = i / n1 * no;
		ih2 = (i - i / n1 * m2) * n2;
#ifdef SX
#pragma cdir nodep
#endif
		for (j = 0; j < n2; j++)
			mato[j + ih1] += mat1[i] * mat2[ih2 + j];
	}
#else // ifndef CBLAS_MMultMatMat
	cblas_dgemm(Order, TransA, TransB, m1, n2, n1, 1., mat1, n1, mat2, n2, 0., mato, no);
#endif

	return 1;
} /* MMultMatMat */

/***************************************************************************
   ROCKFLOW - Funktion: MXPGaussPkt
   Aufgabe:
   X Punkt Gauss-Integration
           bestimmtes Integral
           1/           X
 | P(x)dx = Sigma  Fi*P(xi)
 |
          -1/          i=1
   Formalparameter:
           E: grd  (das X aus der Gleichung oben)
   E: pkt  (welcher von den vielen es sein soll)
   return :
   Ergebnis:
   bestimmtes Integral
   Aenderungen/Korrekturen:
   07/1995     hh        Erste Version

 **************************************************************************/
double MXPGaussPkt(long grd, long pkt)
{
	switch (grd)
	{
		case 1:
			return 0.0;
		case 2:
			switch (pkt)
			{
				case 0:
					return 0.577350269189626;
				case 1:
					return -0.577350269189626;
			}
			break;
		case 3:
			switch (pkt)
			{
				case 0:
					return 0.774596669241483;
				case 1:
					return 0.0;
				case 2:
					return -0.774596669241483;
			}
			break;
		case 4:
			switch (pkt)
			{
				case 0:
					return 0.861136311594053;
				case 1:
					return 0.339981043584856;
				case 2:
					return -0.339981043584856;
				case 3:
					return -0.861136311594053;
			}
			break;
	} /* switch grd */
	return 0.0;
}

/***************************************************************************
   ROCKFLOW - Funktion: MXPGaussFkt
   Aufgabe:
   X Punkt Gauss-Integration
           bestimmtes Integral
           1/           X
 | P(x)dx = Sigma  Fi*P(xi)
 |
          -1/          i=1
   Formalparameter:
           E: grd  (das X uas der Gleichung oben)
   E: pkt  (welcher von den vielen es sein soll)
   return : Faktor fuer Gauss-Integration
   Ergebnis:
   bestimmtes Integral
   Aenderungen/Korrekturen:
   07/1995     hh        Erste Version

 **************************************************************************/
double MXPGaussFkt(long grd, long pkt)
{
	switch (grd)
	{
		case 1:
			return 2.0;
		case 2:
			switch (pkt)
			{
				case 0:
					return 1.0;
				case 1:
					return 1.0;
			}
			break;
		case 3:
			switch (pkt)
			{
				case 0:
					return 0.555555555555556;
				case 1:
					return 0.888888888888889;
				case 2:
					return 0.555555555555556;
			}
			break;
		case 4:
			switch (pkt)
			{
				case 0:
					return 0.347854845137454;
				case 1:
					return 0.652145154862546;
				case 2:
					return 0.652145154862546;
				case 3:
					return 0.347854845137454;
			}
			break;
	} /* switch grd */
	return 0.0;
}

///////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
// NEW. WW

/***************************************************************************
   ROCKFLOW - Funktion: SamplePointTriHQ(const int nsample)

   Aufgabe:
        Provide the sample point for numerical intergral of quadratic
   triangle element. Totally there are three sample points.
   Formalparameter:
         const int nsample: Index of integration points
         double *SPoint   : 0-1 unit coordinates
                            2   Weight
   Programming:
   06/2003     WW        Erste Version
   0.166666666666667 0.166666666666667  0.166666666666667
   0.666666666666667 0.166666666666667  0.166666666666667
   0.166666666666667 0.666666666666667  0.166666666666667
   -------------------------------------------------------
   On Edges
   POINTS & WEIGHTS:
   0.5 0.5 0.166666666666667
   0.0 0.5 0.166666666666667
   0.5 0.0 0.166666666666667
 **************************************************************************/
void SamplePointTriHQ(const int nsample, double* SPoint)
{
	switch (nsample)
	{
		case 0:
			SPoint[0] = 0.166666666666667;
			SPoint[1] = 0.166666666666667;
			SPoint[2] = 0.1666666666667; // Weight
			break;
		case 1:
			SPoint[0] = 0.666666666666667;
			SPoint[1] = 0.166666666666667;
			SPoint[2] = 0.1666666666667; // Weight
			break;
		case 2:
			SPoint[0] = 0.166666666666667;
			SPoint[1] = 0.666666666666667;
			SPoint[2] = 0.1666666666667; // Weight
			break;
		default:
			break;
	}
}

/***************************************************************************
   ROCKFLOW - Funktion: SamplePointTet(const int nsample, double *SPoint)

   Aufgabe:
        Provide the sample point for numerical intergral of tedeahedra element.
        Totally there are three sample points. 5 sample points
   Formalparameter:
         const int nsample: Index of integration points
         double *SPoint   : 0-2 unit coordinates
                            3   Weight
   Programming:
   06/2003     WW        Erste Version

 **************************************************************************/
void SamplePointTet5(const int nsample, double* SPoint)
{
	switch (nsample)
	{
		case 0:
			SPoint[0] = 0.25;
			SPoint[1] = 0.25;
			SPoint[2] = 0.25;
			SPoint[3] = -0.133333333333333; // Weight
			break;
		case 1:
			SPoint[0] = 0.166666666666667;
			SPoint[1] = 0.166666666666667;
			SPoint[2] = 0.166666666666667;
			SPoint[3] = 0.07500000000000; // Weight
			break;
		case 2:
			SPoint[0] = 0.5;
			SPoint[1] = 0.166666666666667;
			SPoint[2] = 0.166666666666667;
			SPoint[3] = 0.07500000000000; // Weight
			break;
		case 3:
			SPoint[0] = 0.166666666666667;
			SPoint[1] = 0.5;
			SPoint[2] = 0.166666666666667;
			SPoint[3] = 0.07500000000000; // Weight
			break;
		case 4:
			SPoint[0] = 0.166666666666667;
			SPoint[1] = 0.166666666666667;
			SPoint[2] = 0.5;
			SPoint[3] = 0.07500000000000; // Weight
			break;
		default:
			break;
	}
}

/***************************************************************************
   ROCKFLOW - Funktion: SamplePointTetHQ(const int nsample, double *SPoint)

   Aufgabe:
        Provide the sample point for numerical intergral of tedeahedra element.
        Totally there are three sample points. 15 sample points
   Formalparameter:
         const int nsample: Index of integration points
         double *SPoint   : 0-2 unit coordinates
                            3   Weight
   Programming:
   06/2003     WW        Erste Version

 **************************************************************************/
void SamplePointTet15(const int nsample, double* SPoint)
{
	switch (nsample)
	{
		case 0:
			SPoint[0] = 0.25;
			SPoint[1] = 0.25;
			SPoint[2] = 0.25;
			SPoint[3] = 0.019753086419753086; // Weight
			break;
		case 1:
			SPoint[0] = 0.09197107805272303;
			SPoint[1] = 0.09197107805272303;
			SPoint[2] = 0.09197107805272303;
			SPoint[3] = 0.011989513963169772; // Weight
			break;
		case 2:
			SPoint[0] = 0.72408676584183096;
			SPoint[1] = 0.09197107805272303;
			SPoint[2] = 0.09197107805272303;
			SPoint[3] = 0.011989513963169772; // Weight
			break;
		case 3:
			SPoint[0] = 0.09197107805272303;
			SPoint[1] = 0.72408676584183096;
			SPoint[2] = 0.09197107805272303;
			SPoint[3] = 0.011989513963169772; // Weight
			break;
		case 4:
			SPoint[0] = 0.09197107805272303;
			SPoint[1] = 0.09197107805272303;
			SPoint[2] = 0.72408676584183096;
			SPoint[3] = 0.011989513963169772; // Weight
			break;
		case 5:
			SPoint[0] = 0.44364916731037080;
			SPoint[1] = 0.05635083268962915;
			SPoint[2] = 0.05635083268962915;
			SPoint[3] = 0.008818342151675485; // Weight
			break;
		case 6:
			SPoint[0] = 0.05635083268962915;
			SPoint[1] = 0.44364916731037080;
			SPoint[2] = 0.05635083268962915;
			SPoint[3] = 0.008818342151675485; // Weight
			break;
		case 7:
			SPoint[0] = 0.05635083268962915;
			SPoint[1] = 0.05635083268962915;
			SPoint[2] = 0.44364916731037080;
			SPoint[3] = 0.008818342151675485; // Weight
			break;
		case 8:
			SPoint[0] = 0.05635083268962915;
			SPoint[1] = 0.44364916731037080;
			SPoint[2] = 0.44364916731037080;
			SPoint[3] = 0.008818342151675485; // Weight
			break;
		case 9:
			SPoint[0] = 0.44364916731037080;
			SPoint[1] = 0.05635083268962915;
			SPoint[2] = 0.44364916731037080;
			SPoint[3] = 0.008818342151675485; // Weight
			break;
		case 10:
			SPoint[0] = 0.44364916731037080;
			SPoint[1] = 0.44364916731037080;
			SPoint[2] = 0.05635083268962915;
			SPoint[3] = 0.008818342151675485; // Weight
			break;
		case 11:
			SPoint[0] = 0.31979362782962989;
			SPoint[1] = 0.31979362782962989;
			SPoint[2] = 0.31979362782962989;
			SPoint[3] = 0.011511367871045397; // Weight
			break;
		case 12:
			SPoint[0] = 0.04061911651111023;
			SPoint[1] = 0.31979362782962989;
			SPoint[2] = 0.31979362782962989;
			SPoint[3] = 0.011511367871045397; // Weight
			break;
		case 13:
			SPoint[0] = 0.31979362782962989;
			SPoint[1] = 0.04061911651111023;
			SPoint[2] = 0.31979362782962989;
			SPoint[3] = 0.011511367871045397; // Weight
			break;
		case 14:
			SPoint[0] = 0.31979362782962989;
			SPoint[1] = 0.31979362782962989;
			SPoint[2] = 0.04061911651111023;
			SPoint[3] = 0.011511367871045397; // Weight
			break;
		default:
			break;
	}
}

/**************************************************************************
   MathLib-Method:
   Task:
   Programing:
   05/2011 NW Implementation
   Last modified:
**************************************************************************/
void SamplePointPyramid5(const int nsample, double* SPoint)
{
	static const double g1 = 0.584237394672177188; //=8/5*sqrt(2/15)
	static const double g2 = -2.0 / 3.0;
	static const double g3 = 2.0 / 5.0;
	static const double w1 = 81.0 / 100.0;
	static const double w2 = 125.0 / 27.0;

	switch (nsample)
	{
		case 0:
			SPoint[0] = -g1;
			SPoint[1] = -g1;
			SPoint[2] = g2;
			SPoint[3] = w1; // Weight
			break;
		case 1:
			SPoint[0] = g1;
			SPoint[1] = -g1;
			SPoint[2] = g2;
			SPoint[3] = w1; // Weight
			break;
		case 2:
			SPoint[0] = g1;
			SPoint[1] = g1;
			SPoint[2] = g2;
			SPoint[3] = w1; // Weight
			break;
		case 3:
			SPoint[0] = -g1;
			SPoint[1] = g1;
			SPoint[2] = g2;
			SPoint[3] = w1; // Weight
			break;
		case 4:
			SPoint[0] = .0;
			SPoint[1] = .0;
			SPoint[2] = g3;
			SPoint[3] = w2; // Weight
			break;
		default:
			break;
	}
}

/**************************************************************************
   MathLib-Method:
   Task:
   Programing:
   05/2011 NW Implementation
   Last modified:
**************************************************************************/
void SamplePointPyramid13(const int nsample, double* SPoint)
{
	static const double g1 = 0.673931986207731726;
	static const double g2 = 0.610639618865075532;
	static const double g3 = 0.580939660561084423;
	static const double g4 = -0.142857142857142857;
	static const double g5 = -0.321428571428571429;
	static const double g6 = 0.524394036075370072;
	static const double g7 = -0.830065359477124183;
	static const double w1 = 0.515003019323671498;
	static const double w2 = 0.257183745242064659;
	static const double w3 = 2.474004977113405936;
	static const double w4 = 0.419515737191525950;

	switch (nsample)
	{
		case 0:
			SPoint[0] = -g1;
			SPoint[1] = -g1;
			SPoint[2] = g4;
			SPoint[3] = w1; // Weight
			break;
		case 1:
			SPoint[0] = g1;
			SPoint[1] = -g1;
			SPoint[2] = g4;
			SPoint[3] = w1; // Weight
			break;
		case 2:
			SPoint[0] = g1;
			SPoint[1] = g1;
			SPoint[2] = g4;
			SPoint[3] = w1; // Weight
			break;
		case 3:
			SPoint[0] = -g1;
			SPoint[1] = g1;
			SPoint[2] = g4;
			SPoint[3] = w1; // Weight
			break;
		case 4:
			SPoint[0] = -g2;
			SPoint[1] = .0;
			SPoint[2] = g5;
			SPoint[3] = w2; // Weight
			break;
		case 5:
			SPoint[0] = g2;
			SPoint[1] = .0;
			SPoint[2] = g5;
			SPoint[3] = w2; // Weight
			break;
		case 6:
			SPoint[0] = .0;
			SPoint[1] = -g2;
			SPoint[2] = g5;
			SPoint[3] = w2; // Weight
			break;
		case 7:
			SPoint[0] = .0;
			SPoint[1] = g2;
			SPoint[2] = g5;
			SPoint[3] = w2; // Weight
			break;
		case 8:
			SPoint[0] = .0;
			SPoint[1] = .0;
			SPoint[2] = g6;
			SPoint[3] = w3; // Weight
			break;
		case 9:
			SPoint[0] = -g3;
			SPoint[1] = -g3;
			SPoint[2] = g7;
			SPoint[3] = w4; // Weight
			break;
		case 10:
			SPoint[0] = g3;
			SPoint[1] = -g3;
			SPoint[2] = g7;
			SPoint[3] = w4; // Weight
			break;
		case 11:
			SPoint[0] = g3;
			SPoint[1] = g3;
			SPoint[2] = g7;
			SPoint[3] = w4; // Weight
			break;
		case 12:
			SPoint[0] = -g3;
			SPoint[1] = g3;
			SPoint[2] = g7;
			SPoint[3] = w4; // Weight
			break;
		default:
			break;
	}
}
void SamplePointPyramid8(const int i, double* SPoint)
{
	static const double g1 = sqrt((double)1.0 / 3.0);
	static const double g2 = (2.0 * sqrt(10.0) - 5.0) / 15.0;
	static const double g3 = -2.0 / 3.0 - g2;
	static const double w1 = 5.0 * (68.0 + 5.0 * sqrt(10.0)) / 432.0;
	static const double w2 = 85.0 / 54.0 - w1;

	switch (i)
	{
		case 0:
			SPoint[0] = -g1;
			SPoint[1] = -g1;
			SPoint[2] = g2;
			SPoint[3] = w1; // Weight
			break;
		case 1:
			SPoint[0] = g1;
			SPoint[1] = -g1;
			SPoint[2] = g2;
			SPoint[3] = w1; // Weight
			break;
		case 2:
			SPoint[0] = g1;
			SPoint[1] = g1;
			SPoint[2] = g2;
			SPoint[3] = w1; // Weight
			break;
		case 3:
			SPoint[0] = -g1;
			SPoint[1] = g1;
			SPoint[2] = g2;
			SPoint[3] = w1; // Weight
			break;
		case 4:
			SPoint[0] = -g1;
			SPoint[1] = -g1;
			SPoint[2] = g3;
			SPoint[3] = w2; // Weight
			break;
		case 5:
			SPoint[0] = g1;
			SPoint[1] = -g1;
			SPoint[2] = g3;
			SPoint[3] = w2; // Weight
			break;
		case 6:
			SPoint[0] = g1;
			SPoint[1] = g1;
			SPoint[2] = g3;
			SPoint[3] = w2; // Weight
			break;
		case 7:
			SPoint[0] = -g1;
			SPoint[1] = g1;
			SPoint[2] = g3;
			SPoint[3] = w2; // Weight
			break;
	}
}
/**************************************************************************
   MathLib-Method:
   Task:
   Programing:
   02/2005 OK Implementation
   Last modified:
**************************************************************************/
void ShapeFunctionLine(double* N1, const double* u)
{
	int i;
	N1[0] = 1.0 - u[0];
	N1[1] = 1.0 + u[0];
	//  N1[0] = 1.0 + u[0];
	//  N1[1] = 1.0 - u[0];
	for (i = 0; i < 2; i++)
		N1[i] *= 0.5;
}

/**************************************************************************
   MathLib-Method:
   Task:
   Programing:
   07/2003 WW Implementation
   Last modified:
**************************************************************************/
void ShapeFunctionLineHQ(double* N1, const double* u)
{
	N1[0] = 0.5 * u[0] * (u[0] - 1.0);
	N1[1] = 0.5 * u[0] * (u[0] + 1.0);
	N1[2] = 1.0 - u[0] * u[0];
}

/**************************************************************************
   MathLib-Method:
   Task:
   Programing:
   02/2005 WW Implementation
   01/2006 PCH Correct the sign.
   Last modified:
**************************************************************************/
void GradShapeFunctionLine(double* dN1, const double* u)
{
	u = u;
	dN1[0] = -0.5;
	dN1[1] = 0.5;
}

/**************************************************************************
   MathLib-Method:
   Task:
   Programing:
   04/2009 NW Implementation
   Last modified:
**************************************************************************/
void GradShapeFunctionLineHQ(double* dN1, const double* u)
{
	dN1[0] = u[0] - 0.5;
	dN1[1] = u[0] + 0.5;
	dN1[2] = -2.0 * u[0];
}

/***************************************************************************
   ROCKFLOW - Funktion:  ShapeFunctionTri

   Aufgabe:
        Compute the shape function for numerical intergral of linear
   triangle element.
   Formalparameter:
      E:
        const  double *u  : Pointer to a 2 dimension array for
                               the unit coordinates
      R:
   double * N3       : Array of size 3, to store the value of shape
   function

   Programming:
   06/2003     WW        Erste Version

 **************************************************************************/
void ShapeFunctionTri(double* N3, const double* u)
{
	N3[0] = 1. - u[0] - u[1];
	N3[1] = u[0];
	N3[2] = u[1];
}

/***************************************************************************
   ROCKFLOW - Funktion:  *GradShapeFunctionTri

   Aufgabe:
        Compute the shape function for numerical intergral of linear
   triangle element.
   Formalparameter:
      E:
       double * N3       : Array of size 6, gradient of shape function
                          0--2 d()/dL_1
                          3--5 d()/dL_2
   const  double *u  : Pointer to a 2 dimension array for
   the unit coordinates

   Programming:
   06/2006     WW        Erste Version

 **************************************************************************/
void GradShapeFunctionTri(double* dN3, const double* u)
{
	u = u;
	//   d()/dL_1
	dN3[0] = -1.0;
	dN3[1] = 1.0;
	dN3[2] = 0.0;
	//   d()/dL_2
	dN3[3] = -1.0;
	dN3[4] = 0.0;
	dN3[5] = 1.0;
}

/***************************************************************************
   ROCKFLOW - Funktion:  ShapeFunctionTriHQ

   Aufgabe:
        Compute the shape function for numerical intergral of quadratic
   triangle element.
   Formalparameter:
      E:
        const  double *u  : Pointer to a 2 dimension array for
                               the unit coordinates
      R:
   double * N3       : Array of size 6, to store the value of shape
   function

   Programming:
   06/2003     WW        Erste Version

 **************************************************************************/
void ShapeFunctionTriHQ(double* N6, const double* u)
{
	N6[0] = 2. * (1. - u[0] - u[1]) * (0.5 - u[0] - u[1]);
	N6[1] = u[0] * (2. * u[0] - 1.);
	N6[2] = u[1] * (2. * u[1] - 1.);
	N6[3] = 4. * u[0] * (1. - u[0] - u[1]);
	N6[4] = 4. * u[0] * u[1];
	N6[5] = 4. * u[1] * (1. - u[0] - u[1]);
}

/***************************************************************************
   ROCKFLOW - Funktion:  *GradShapeFunctionTriHQ

   Aufgabe:
        Compute the shape function for numerical intergral of quadratic
   triangle element.
   Formalparameter:
      E:
       double * N3       : Array of size 6, gradient of shape function
                          0--5  d()/dL_1
                          6--11 d()/dL_2
   const  double *u  : Pointer to a 2 dimension array for
   the unit coordinates

   Programming:
   06/2006     WW        Erste Version

 **************************************************************************/
void GradShapeFunctionTriHQ(double* dN6, const double* u)
{
	dN6[0] = 4. * (u[0] + u[1]) - 3.; // dN1/dL1
	dN6[6] = 4. * (u[0] + u[1]) - 3.; // dN1/dL2

	dN6[1] = 4. * u[0] - 1.; // dN2/dL1
	dN6[7] = 0.; // dN2/dL2

	dN6[2] = 0.; // dN3/dL1
	dN6[8] = 4. * u[1] - 1.; // dN3/dL2

	dN6[3] = 4. * (1 - 2. * u[0] - u[1]); // dN4/dL1
	dN6[9] = -4. * u[0]; // dN4/dL2

	dN6[4] = 4. * u[1]; // dN5/dL1
	dN6[10] = 4. * u[0]; // dN5/dL2

	dN6[5] = -4. * u[1]; // dN6/dL1
	dN6[11] = 4. * (1 - u[0] - 2. * u[1]); // dN6/dL2
}

/***************************************************************************
   ROCKFLOW - Funktion: realCoordTriHQ

   Aufgabe:
        Mapping to real coordaintes from the local ones of quadratic traingle
   element.
   Formalparameter:
           E:
             double * x         : Array of size 3, real coordiantes
             const double *XY   : Array of size 12, to store the coordinates of
                                  the six  verteces as:
   0-5;  x1, ..., x6,
   6-11; y1, ..., y6
   const double *u    : Array of size 2, unit coordiantes

   Programming:
   06/2003     WW        Erste Version

 **************************************************************************/
void realCoordTriHQ(double* x, const double* XY, const double* u)
{
	x[1] = (1.0 - u[0] - u[1]) * XY[0] + u[0] * XY[1] + u[1] * XY[2];
	x[1] = (1.0 - u[0] - u[1]) * XY[6] + u[0] * XY[7] + u[1] * XY[8];
}

/***************************************************************************
   ROCKFLOW - Funktion:  ShapeFunctionQuad

   Aufgabe:
        Compute the shape function for numerical intergral of linear
   quadralateral element.
   Formalparameter:
           E:
             double * N4     : Array of size 4, to store the value of shape
                                    function
             const  double *u   : Pointer to a 2 dimension array for
   the unit coordinates

   Programming:
   06/2004     WW        Erste Version

 **************************************************************************/
void ShapeFunctionQuad(double* N4, const double* u)
{
	int i;
	N4[0] = (1.0 + u[0]) * (1.0 + u[1]);
	N4[1] = (1.0 - u[0]) * (1.0 + u[1]);
	N4[2] = (1.0 - u[0]) * (1.0 - u[1]);
	N4[3] = (1.0 + u[0]) * (1.0 - u[1]);
	for (i = 0; i < 4; i++)
		N4[i] *= 0.25;
}

/***************************************************************************
   ROCKFLOW - Funktion:  ShapeFunctionQuad

   Aufgabe:
        Compute the shape function for numerical intergral of linear
   quadralateral element.
   Formalparameter:
           E:
             double * dN4     : Array of size 8, to store the gradient
                                of shape function.
             const  double *u   : Pointer to a 2 dimension array for
   the unit coordinates

   Programming:
   06/2004     WW        Erste Version

 **************************************************************************/
void GradShapeFunctionQuad(double* dN4, const double* u)
{
	dN4[0] = +(1.0 + u[1]);
	dN4[1] = -(1.0 + u[1]);
	dN4[2] = -(1.0 - u[1]);
	dN4[3] = +(1.0 - u[1]);
	dN4[4] = +(1.0 + u[0]);
	dN4[5] = +(1.0 - u[0]);
	dN4[6] = -(1.0 - u[0]);
	dN4[7] = -(1.0 + u[0]);
	for (int i = 0; i < 8; i++)
		dN4[i] *= 0.25;
}

/***************************************************************************
   Aufgabe:
        Compute the shape function for numerical intergral of linear
   quadratic quadralateral element.(8 nodes)
   Formalparameter:
      E:
        const  double *u  : Pointer to a 2 dimension array for
                               the unit coordinates
      R:
       double * N8       : Array of size 8, to store the value of shape
                               function

   Programming:
   08/2004     WW              Generalization
 **************************************************************************/
void ShapeFunctionQuadHQ8(double* N8, const double* u)
{
	//
	N8[0] = -0.25 * (1.0 - u[0]) * (1.0 - u[1]) * ((1.0 + u[0] + u[1]));
	N8[1] = 0.25 * (1.0 + u[0]) * (1.0 - u[1]) * ((-1.0 + u[0] - u[1]));
	N8[2] = 0.25 * (1.0 + u[0]) * (1.0 + u[1]) * ((-1.0 + u[0] + u[1]));
	N8[3] = -0.25 * (1.0 - u[0]) * (1.0 + u[1]) * ((1.0 + u[0] - u[1]));
	//
	N8[4] = 0.5 * (1.0 - u[0] * u[0]) * (1.0 - u[1]);
	N8[5] = 0.5 * (1.0 - u[1] * u[1]) * (1.0 + u[0]);
	N8[6] = 0.5 * (1.0 - u[0] * u[0]) * (1.0 + u[1]);
	N8[7] = 0.5 * (1.0 - u[1] * u[1]) * (1.0 - u[0]);
}

/***************************************************************************
   Aufgabe:
        Compute the shape function for numerical intergral of linear
   quadratic quadralateral element.
   Formalparameter:
      E:
        const  double *u  : Pointer to a 2 dimension array for
                               the unit coordinates
      R:
       double * N9       : Array of size 9, to store the value of shape
                               function

   Programming:
   12/1999     R.Kaiser        Erste Version
   06/2004     WW              Generalization
 **************************************************************************/
void ShapeFunctionQuadHQ(double* N9, const double* u)
{
	N9[8] = (1.0 - u[0] * u[0]) * (1.0 - u[1] * u[1]);
	N9[7] = 0.5 * (1.0 - u[1] * u[1]) * (1.0 + u[0]) - 0.5 * N9[8];
	N9[6] = 0.5 * (1.0 - u[0] * u[0]) * (1.0 - u[1]) - 0.5 * N9[8];
	N9[5] = 0.5 * (1.0 - u[1] * u[1]) * (1.0 - u[0]) - 0.5 * N9[8];
	N9[4] = 0.5 * (1.0 - u[0] * u[0]) * (1.0 + u[1]) - 0.5 * N9[8];
	N9[3] = 0.25 * (1.0 + u[0]) * (1.0 - u[1]) - 0.5 * N9[6] - 0.5 * N9[7] - 0.25 * N9[8];
	N9[2] = 0.25 * (1.0 - u[0]) * (1.0 - u[1]) - 0.5 * N9[5] - 0.5 * N9[6] - 0.25 * N9[8];
	N9[1] = 0.25 * (1.0 - u[0]) * (1.0 + u[1]) - 0.5 * N9[4] - 0.5 * N9[5] - 0.25 * N9[8];
	N9[0] = 0.25 * (1.0 + u[0]) * (1.0 + u[1]) - 0.5 * N9[4] - 0.5 * N9[7] - 0.25 * N9[8];
}

/***************************************************************************
   Aufgabe:
        Compute the shape function for numerical intergral of linear
   quadratic quadralateral element.
   Formalparameter:
      E:
       double * N9       : Array of size 9, to store the value of shape
                               function
        const  double *u  : Pointer to a 2 dimension array for
                               the unit coordinates

   Programming:
   12/1999     RK        Erste Version
   06/2004     WW
 **************************************************************************/
void GradShapeFunctionQuadHQ(double* dN9, const double* u)
{
	dN9[8] = -2.0 * u[0] * (1.0 - u[1] * u[1]);
	dN9[7] = +0.5 * (1.0 - u[1] * u[1]) - 0.5 * dN9[8];
	dN9[6] = -1.0 * u[0] * (1.0 - u[1]) - 0.5 * dN9[8];
	dN9[5] = -0.5 * (1.0 - u[1] * u[1]) - 0.5 * dN9[8];
	dN9[4] = -1.0 * u[0] * (1.0 + u[1]) - 0.5 * dN9[8];
	dN9[3] = +0.25 * (1 - u[1]) - 0.5 * dN9[6] - 0.5 * dN9[7] - 0.25 * dN9[8];
	dN9[2] = -0.25 * (1 - u[1]) - 0.5 * dN9[5] - 0.5 * dN9[6] - 0.25 * dN9[8];
	dN9[1] = -0.25 * (1 + u[1]) - 0.5 * dN9[4] - 0.5 * dN9[5] - 0.25 * dN9[8];
	dN9[0] = +0.25 * (1 + u[1]) - 0.5 * dN9[4] - 0.5 * dN9[7] - 0.25 * dN9[8];

	dN9[17] = -2.0 * u[1] * (1.0 - u[0] * u[0]);
	dN9[16] = -1.0 * u[1] * (1.0 + u[0]) - 0.5 * dN9[17];
	dN9[15] = -0.5 * (1.0 - u[0] * u[0]) - 0.5 * dN9[17];
	dN9[14] = -1.0 * u[1] * (1.0 - u[0]) - 0.5 * dN9[17];
	dN9[13] = +0.5 * (1 - u[0] * u[0]) - 0.5 * dN9[17];
	dN9[12] = -0.25 * (1 + u[0]) - 0.5 * dN9[15] - 0.5 * dN9[16] - 0.25 * dN9[17];
	dN9[11] = -0.25 * (1 - u[0]) - 0.5 * dN9[14] - 0.5 * dN9[15] - 0.25 * dN9[17];
	dN9[10] = +0.25 * (1 - u[0]) - 0.5 * dN9[13] - 0.5 * dN9[14] - 0.25 * dN9[17];
	dN9[9] = +0.25 * (1 + u[0]) - 0.5 * dN9[13] - 0.5 * dN9[16] - 0.25 * dN9[17];
}

/***************************************************************************
   Aufgabe:
        Compute the shape function for numerical intergral of linear
   quadratic quadralateral element.
   Formalparameter:
      E:
        double * N8       : Array of size 8, to store the value of shape
                               function
        const  double *u  : Pointer to a 2 dimension array for
                               the unit coordinates

   Programming:
   01/2010     NW Implementation
 **************************************************************************/
void GradShapeFunctionQuadHQ8(double* dN8, const double* u)
{
	double r = u[0];
	double s = u[1];

	// dN/dr
	dN8[0] = (1 - s) * (2 * r + s) * 0.25;
	dN8[1] = (1 - s) * (2 * r - s) * 0.25;
	dN8[2] = (1 + s) * (2 * r + s) * 0.25;
	dN8[3] = (1 + s) * (2 * r - s) * 0.25;
	dN8[4] = -r * (1 - s);
	dN8[5] = (1 - s * s) * 0.5;
	dN8[6] = -r * (1 + s);
	dN8[7] = -(1 - s * s) * 0.5;

	// dN/ds
	dN8[8] = (1 - r) * (r + 2 * s) * 0.25;
	dN8[9] = -(1 + r) * (r - 2 * s) * 0.25;
	dN8[10] = (1 + r) * (r + 2 * s) * 0.25;
	dN8[11] = -(1 - r) * (r - 2 * s) * 0.25;
	dN8[12] = -(1 - r * r) * 0.5;
	dN8[13] = -(1 + r) * s;
	dN8[14] = (1 - r * r) * 0.5;
	dN8[15] = -(1 - r) * s;
}

/***************************************************************************
   Aufgabe:
        Compute the shape function for numerical intergral of
   linear tedrahedra element.
   Formalparameter:
      E:
       double * Nt4      : Array of size 4, to store the value of shape
                            function
        const  double *x  : Pointer to a 3 dimension array for

   Programming:
   09/2004     WW
 **************************************************************************/
void ShapeFunctionTet(double* Nt4, const double* x)
{
	Nt4[0] = 1. - x[0] - x[1] - x[2];
	Nt4[1] = x[0];
	Nt4[2] = x[1];
	Nt4[3] = x[2];
}

/***************************************************************************
   Aufgabe:
        Compute the gradient of shape functions for numerical intergral of
   linear tedrahedra element.
   Formalparameter:
      E:
       double * dNt4      : Array of size 12, to store the value of shape
                             function
                             0--3 /dr
                             4--7 /ds
                             8--12 /dt
   const  double *x   : Pointer to a 3 dimension array for

   Programming:
   09/2004     WW
 **************************************************************************/
void GradShapeFunctionTet(double* dNt4, const double* dummy)
{
	dummy = dummy;
	dNt4[0] = -1.0;
	dNt4[1] = 1.0;
	dNt4[2] = 0.0;
	dNt4[3] = 0.0;

	dNt4[4] = -1.0;
	dNt4[5] = 0.0;
	dNt4[6] = 1.0;
	dNt4[7] = 0.0;

	dNt4[8] = -1.0;
	dNt4[9] = 0.0;
	dNt4[10] = 0.0;
	dNt4[11] = 1.0;
}

/***************************************************************************
   Aufgabe:
        Compute the shape function for numerical intergral of
   quadratic tedrahedra element.
   Formalparameter:
      E:
        const  double *x  : Pointer to a 3 dimension array for
                            the unit coordinates
      R:
       double * N10      : Array of size 10, to store the value of shape
                            function

   Programming:
   09/2004     WW
 **************************************************************************/
void ShapeFunctionTetHQ(double* N10, const double* x)
{
	N10[0] = 2. * (1 - x[0] - x[1] - x[2]) * (0.5 - x[0] - x[1] - x[2]);
	N10[1] = x[0] * (2. * x[0] - 1);
	N10[2] = x[1] * (2. * x[1] - 1);
	N10[3] = x[2] * (2. * x[2] - 1);
	N10[4] = 4.0 * x[0] * (1.0 - x[0] - x[1] - x[2]);
	N10[5] = 4.0 * x[0] * x[1];
	N10[6] = 4.0 * x[1] * (1.0 - x[0] - x[1] - x[2]);
	N10[7] = 4.0 * x[0] * x[2];
	N10[8] = 4.0 * x[1] * x[2];
	N10[9] = 4.0 * x[2] * (1.0 - x[0] - x[1] - x[2]);
}

/***************************************************************************
   Aufgabe:
        Compute the gradient of the shape function for numerical intergral
    of quadratic tedrahedra element with respect to unit coordinates.
   Formalparameter:
      E:
        const  double *x  : Pointer to a 3 dimension array for
                            the unit coordinates
      R:
       double * dN10     : Array of size 30, to store the value
                            0 - 9, d()/dr
   10-19, d()/ds
   20-29, d()/dt
   Programming:
   09/2004     WW
 **************************************************************************/
void GradShapeFunctionTetHQ(double* dN10, const double* x)
{
	dN10[0] = 4.0 * (x[0] + x[1] + x[2]) - 3.0;
	dN10[1] = 4. * x[0] - 1.;
	dN10[2] = 0.0;
	dN10[3] = 0.0;
	dN10[4] = 4.0 * (1.0 - 2.0 * x[0] - x[1] - x[2]);
	dN10[5] = 4.0 * x[1];
	dN10[6] = -4.0 * x[1];
	dN10[7] = 4.0 * x[2];
	dN10[8] = 0.0;
	dN10[9] = -4.0 * x[2];

	dN10[10] = 4. * (x[0] + x[1] + x[2]) - 3.;
	dN10[11] = 0.0;
	dN10[12] = 4. * x[1] - 1.;
	dN10[13] = 0.;
	dN10[14] = -4.0 * x[0];
	dN10[15] = 4.0 * x[0];
	dN10[16] = 4.0 * (1.0 - x[0] - 2.0 * x[1] - x[2]);
	dN10[17] = 0.0;
	dN10[18] = 4.0 * x[2];
	dN10[19] = -4.0 * x[2];

	dN10[20] = 4. * (x[0] + x[1] + x[2]) - 3.;
	dN10[21] = 0.;
	dN10[22] = 0.;
	dN10[23] = 4. * x[2] - 1.;
	dN10[24] = -4.0 * x[0];
	dN10[25] = 0.0;
	dN10[26] = -4.0 * x[1];
	dN10[27] = 4.0 * x[0];
	dN10[28] = 4.0 * x[1];
	dN10[29] = 4.0 * (1.0 - x[0] - x[1] - 2.0 * x[2]);
}

/***************************************************************************
   ROCKFLOW - Funktion: MOmega3D-ShapeFunctionHex
   Aufgabe:
           Berechnet Ansatzfunktion (r,s).
                        / (1+r)(1+s)(1+t) \
 | (1-r)(1+s)(1+t) |
                     1  | (1-r)(1-s)(1+t) |
            N8 = --- | (1+r)(1-s)(1+t) |
                     8  | (1+r)(1+s)(1-t) |
 | (1-r)(1+s)(1-t) |
 | (1-r)(1-s)(1-t) |
 \ (1+r)(1-s)(1-t) /
   Formalparameter:
   Z: *vf - 1x8 Feld
   E: r,s,t (s.o.)
   Ergebnis:
   Vektor
   Aenderungen/Korrekturen:
   08/1995     cb        Erste Version
   09/2004     WW     Generalization
 **************************************************************************/
void ShapeFunctionHex(double* N8, const double* x)
{
	int i;
	N8[0] = (1.0 + x[0]) * (1.0 + x[1]) * (1.0 + x[2]);
	N8[1] = (1.0 - x[0]) * (1.0 + x[1]) * (1.0 + x[2]);
	N8[2] = (1.0 - x[0]) * (1.0 - x[1]) * (1.0 + x[2]);
	N8[3] = (1.0 + x[0]) * (1.0 - x[1]) * (1.0 + x[2]);
	N8[4] = (1.0 + x[0]) * (1.0 + x[1]) * (1.0 - x[2]);
	N8[5] = (1.0 - x[0]) * (1.0 + x[1]) * (1.0 - x[2]);
	N8[6] = (1.0 - x[0]) * (1.0 - x[1]) * (1.0 - x[2]);
	N8[7] = (1.0 + x[0]) * (1.0 - x[1]) * (1.0 - x[2]);
	for (i = 0; i < 8; i++)
		N8[i] *= 0.125;
}

/***************************************************************************
   ROCKFLOW - Funktion: MGradOmega3D
   Aufgabe:
           Berechnet Gradient eines 3D Vektorfeldes,
           dessen Ansatzfunktion (Vektor) bekannt ist (r,s,t).

 |  0   1   2   3   4   5   6   7  |
 |                                 |
      grad (omega)= |  8   9  10  11  12  13  14  15  |
 |                                 |
 | 16  17  18  19  20  21  22  23  |

   Zahlen in der Matrix stehen fuer Positionen im Feld vf.
   Gleichungen s.u..

   Formalparameter:
   Z: *vf - 3x8 Feld
   E: r,s (s.o.)
   Ergebnis:
   3x8 Matrix
   Aenderungen/Korrekturen:
   05/1995     hh        Erste Version
   09/2004     WW        Generalization
 **************************************************************************/

void GradShapeFunctionHex(double* dN8, const double* x)
{
	int i;
	double r = x[0];
	double s = x[1];
	double t = x[2];
	dN8[0] = +(1.0 + s) * (1.0 + t);
	dN8[1] = -(1.0 + s) * (1.0 + t);
	dN8[2] = -(1.0 - s) * (1.0 + t);
	dN8[3] = +(1.0 - s) * (1.0 + t);

	dN8[4] = +(1.0 + s) * (1.0 - t);
	dN8[5] = -(1.0 + s) * (1.0 - t);
	dN8[6] = -(1.0 - s) * (1.0 - t);
	dN8[7] = +(1.0 - s) * (1.0 - t);

	dN8[8] = +(1.0 + r) * (1.0 + t);
	dN8[9] = +(1.0 - r) * (1.0 + t);
	dN8[10] = -(1.0 - r) * (1.0 + t);
	dN8[11] = -(1.0 + r) * (1.0 + t);

	dN8[12] = +(1.0 + r) * (1.0 - t);
	dN8[13] = +(1.0 - r) * (1.0 - t);
	dN8[14] = -(1.0 - r) * (1.0 - t);
	dN8[15] = -(1.0 + r) * (1.0 - t);

	dN8[16] = +(1.0 + r) * (1.0 + s);
	dN8[17] = +(1.0 - r) * (1.0 + s);
	dN8[18] = +(1.0 - r) * (1.0 - s);
	dN8[19] = +(1.0 + r) * (1.0 - s);

	dN8[20] = -(1.0 + r) * (1.0 + s);
	dN8[21] = -(1.0 - r) * (1.0 + s);
	dN8[22] = -(1.0 - r) * (1.0 - s);
	dN8[23] = -(1.0 + r) * (1.0 - s);

	for (i = 0; i < 24; i++)
		dN8[i] *= 0.125;
}

/***************************************************************************
   GEOSYS - Funktion: ShapeFunctionHexHQ
   Task:
     Shape functions for the 20 node hexahedral element
     (Including:
             ShapeFunctionHexHQ_Corner
             ShapeFunctionHexHQ_Middle
             dShapeFunctionHexHQ_Corner
             dShapeFunctionHexHQ_Middle
             )
   Arguments:
   E:
   const  double *x  : Pointer to a 3 dimension array for
   the unit coordinates
   R:
   double * N20      : Array of size 20, to store the value of shape
   function
   Programming:
   09/2004     WW
 **************************************************************************/
double ShapeFunctionHexHQ_Corner(const double r, const double s, const double t)
{
	return 0.125 * (1 + r) * (1 + s) * (1 + t) * (r + s + t - 2.0);
}

double ShapeFunctionHexHQ_Middle(const double r, const double s, const double t)
{
	return 0.25 * (1 - r * r) * (1 + s) * (1 + t);
}

double dShapeFunctionHexHQ_Corner(const double r, const double s, const double t, const int ty)
{
	switch (ty)
	{
		case 0:
			return 0.125 * (1 + s) * (1 + t) * (2.0 * r + s + t - 1.0);
			break;
		case 1:
			return 0.125 * (1 + t) * (1 + r) * (2.0 * s + r + t - 1.0);
			break;
		case 2:
			return 0.125 * (1 + r) * (1 + s) * (2.0 * t + s + r - 1.0);
			break;
	}
	return 0.0;
}

double dShapeFunctionHexHQ_Middle(const double r, const double s, const double t, const int ty)
{
	switch (ty)
	{
		case 0:
			return -0.5 * r * (1 + s) * (1 + t);
			break;
		case 1:
			return 0.25 * (1 - r * r) * (1 + t);
			break;
		case 2:
			return 0.25 * (1 - r * r) * (1 + s);
			break;
	}
	return 0.0;
}

void ShapeFunctionHexHQ(double* N20, const double* x)
{
	double r = x[0];
	double s = x[1];
	double t = x[2];

	N20[0] = ShapeFunctionHexHQ_Corner(r, s, t);
	N20[1] = ShapeFunctionHexHQ_Corner(-r, s, t);
	N20[2] = ShapeFunctionHexHQ_Corner(-r, -s, t);
	N20[3] = ShapeFunctionHexHQ_Corner(r, -s, t);
	N20[4] = ShapeFunctionHexHQ_Corner(r, s, -t);
	N20[5] = ShapeFunctionHexHQ_Corner(-r, s, -t);
	N20[6] = ShapeFunctionHexHQ_Corner(-r, -s, -t);
	N20[7] = ShapeFunctionHexHQ_Corner(r, -s, -t);

	N20[8] = ShapeFunctionHexHQ_Middle(r, s, t);
	N20[10] = ShapeFunctionHexHQ_Middle(r, -s, t);
	N20[14] = ShapeFunctionHexHQ_Middle(r, -s, -t);
	N20[12] = ShapeFunctionHexHQ_Middle(r, s, -t);

	N20[11] = ShapeFunctionHexHQ_Middle(s, t, r);
	N20[15] = ShapeFunctionHexHQ_Middle(s, -t, r);
	N20[13] = ShapeFunctionHexHQ_Middle(s, -t, -r);
	N20[9] = ShapeFunctionHexHQ_Middle(s, t, -r);

	N20[16] = ShapeFunctionHexHQ_Middle(t, r, s);
	N20[17] = ShapeFunctionHexHQ_Middle(t, -r, s);
	N20[18] = ShapeFunctionHexHQ_Middle(t, -r, -s);
	N20[19] = ShapeFunctionHexHQ_Middle(t, r, -s);
}

/***************************************************************************
   GEOSYS - Funktion: ShapeFunctionHexHQ
   Task:
        Compute the gradient of shape functions
        for the 20 node hexahedral element
   Arguments:
      E:
        const  double *x  : Pointer to a 3 dimension array for
                            the unit coordinates
      R:
       double * dN20     : Array of size 60, to store the value of shape
   function
   0 - 9, d()/dr
   10-19, d()/ds
   20-29, d()/dt
   Programming:
   09/2004     WW
 **************************************************************************/
void GradShapeFunctionHexHQ(double* dN20, const double* x)
{
	int co;
	double r = x[0];
	double s = x[1];
	double t = x[2];
	static double sign1[] = {-1.0, 1.0, 1.0};
	static double sign2[] = {1.0, -1.0, 1.0};
	static double sign3[] = {1.0, 1.0, -1.0};
	for (int i = 0; i < 3; i++)
	{
		dN20[20 * i + 0] = dShapeFunctionHexHQ_Corner(r, s, t, i);
		dN20[20 * i + 1] = sign1[i] * dShapeFunctionHexHQ_Corner(-r, s, t, i);
		dN20[20 * i + 2] = sign1[i] * sign2[i] * dShapeFunctionHexHQ_Corner(-r, -s, t, i);
		dN20[20 * i + 3] = sign2[i] * dShapeFunctionHexHQ_Corner(r, -s, t, i);
		dN20[20 * i + 4] = sign3[i] * dShapeFunctionHexHQ_Corner(r, s, -t, i);
		dN20[20 * i + 5] = sign1[i] * sign3[i] * dShapeFunctionHexHQ_Corner(-r, s, -t, i);
		dN20[20 * i + 6] = sign1[i] * sign2[i] * sign3[i] * dShapeFunctionHexHQ_Corner(-r, -s, -t, i);
		dN20[20 * i + 7] = sign2[i] * sign3[i] * dShapeFunctionHexHQ_Corner(r, -s, -t, i);

		dN20[20 * i + 8] = dShapeFunctionHexHQ_Middle(r, s, t, i);
		dN20[20 * i + 10] = sign2[i] * dShapeFunctionHexHQ_Middle(r, -s, t, i);
		dN20[20 * i + 14] = sign2[i] * sign3[i] * dShapeFunctionHexHQ_Middle(r, -s, -t, i);
		dN20[20 * i + 12] = sign3[i] * dShapeFunctionHexHQ_Middle(r, s, -t, i);

		co = (i + 2) % 3;
		dN20[20 * i + 11] = dShapeFunctionHexHQ_Middle(s, t, r, co);
		dN20[20 * i + 15] = sign3[i] * dShapeFunctionHexHQ_Middle(s, -t, r, co);
		dN20[20 * i + 13] = sign1[i] * sign3[i] * dShapeFunctionHexHQ_Middle(s, -t, -r, co);
		dN20[20 * i + 9] = sign1[i] * dShapeFunctionHexHQ_Middle(s, t, -r, co);

		co = (i + 1) % 3;
		dN20[20 * i + 16] = dShapeFunctionHexHQ_Middle(t, r, s, co);
		dN20[20 * i + 17] = sign1[i] * dShapeFunctionHexHQ_Middle(t, -r, s, co);
		dN20[20 * i + 18] = sign1[i] * sign2[i] * dShapeFunctionHexHQ_Middle(t, -r, -s, co);
		dN20[20 * i + 19] = sign2[i] * dShapeFunctionHexHQ_Middle(t, r, -s, co);
	}
}

/***************************************************************************
   GEOSYS/ROCKFLOW - Funktion: ShapeFunctionPri
   Task:
        Compute the gradient of shape functions
        for the 6 node prism element
   Arguments:
      E:
        const  double *x  : Pointer to a 3 dimension array for
                            the unit coordinates
      R:
       double * N     : Array of size 6, to store the value of shape
   function of 6 nodes
   Programming:
   08/2005     WW   (duplicated from MOmegaPrism by MB, 07/2003)
 **************************************************************************/
void ShapeFunctionPri(double* N, const double* x)
{
	double L1 = x[0];
	double L2 = x[1];
	double t = x[2];
	N[0] = 0.5 * (1.0 - L1 - L2) * (1.0 - t);
	N[1] = 0.5 * L1 * (1.0 - t);
	N[2] = 0.5 * L2 * (1.0 - t);
	N[3] = 0.5 * (1.0 - L1 - L2) * (1.0 + t);
	N[4] = 0.5 * L1 * (1.0 + t);
	N[5] = 0.5 * L2 * (1.0 + t);
}

/***************************************************************************
   GEOSYS/ROCKFLOW - Funktion: ShapeFunctionPriHQ
   Task:
        Compute the gradient of shape functions
        for the 6 node prism element
   Arguments:
      E:
        const  double *x  : Pointer to a 3 dimension array for
                            the unit coordinates
      R:
       double * N     : Array of size 15, to store the value of shape
   function of 15 nodes
   Programming:
   08/2005     WW
   09/2014     WW  Efficency improvement.
 **************************************************************************/
void ShapeFunctionPriHQ(double* N, const double* x)
{
	double L1 = x[0];
	double L2 = x[1];
	const double L0 = 1.0 - L1 - L2;
	double t = x[2];
	double tt1 = 1.0 - t * t;

	double v1 = 2.0 * L0 - 1;
	double v2 = 2.0 * L1 - 1;
	double v3 = 2.0 * L2 - 1;
	// Vertex, bottom
	N[0] = 0.5 * L0 * (v1 * (1.0 - t) - tt1);
	N[1] = 0.5 * L1 * (v2 * (1.0 - t) - tt1);
	N[2] = 0.5 * L2 * (v3 * (1.0 - t) - tt1);
	// Vertex, top
	N[3] = 0.5 * L0 * (v1 * (1.0 + t) - tt1);
	N[4] = 0.5 * L1 * (v2 * (1.0 + t) - tt1);
	N[5] = 0.5 * L2 * (v3 * (1.0 + t) - tt1);

	v1 = 2.0 * L0 * L1;
	v2 = 2.0 * L1 * L2;
	v3 = 2.0 * L2 * L0;
	// Middle point, bottom
	N[6] = v1 * (1.0 - t);
	N[7] = v2 * (1.0 - t);
	N[8] = v3 * (1.0 - t);
	// Middle point, top
	N[9] = v1 * (1.0 + t);
	N[10] = v2 * (1.0 + t);
	N[11] = v3 * (1.0 + t);
	// Middle point, center
	N[12] = L0 * tt1;
	N[13] = L1 * tt1;
	N[14] = L2 * tt1;
}

/***************************************************************************
   GEOSYS/ROCKFLOW - Funktion: GradShapeFunctionPri
   Task:
        Compute the gradient of shape functions
        for the 6 node prism element
   Arguments:
      E:
        const  double *x  : Pointer to a 3 dimension array for
                            the unit coordinates
      R:
       double * N     : Array of size 18, to store the value of the grandient
   shape functions
   0--5: dN/dL1
   6--11: dN/dL2
   12--17: dN/dt
   Programming:
   08/2005     WW
 **************************************************************************/
void GradShapeFunctionPri(double* dN, const double* x)
{
	double L1 = x[0];
	double L2 = x[1];
	double t = x[2];
	//  dN/dL1
	dN[0] = -0.5 * (1.0 - t);
	dN[1] = 0.5 * (1.0 - t);
	dN[2] = 0.0;
	dN[3] = -0.5 * (1.0 + t);
	dN[4] = 0.5 * (1.0 + t);
	dN[5] = 0.0;
	//  dN/dL2
	dN[6] = -0.5 * (1.0 - t);
	dN[7] = 0.0;
	dN[8] = 0.5 * (1.0 - t);
	dN[9] = -0.5 * (1.0 + t);
	dN[10] = 0.0;
	dN[11] = 0.5 * (1.0 + t);
	//  dN/dt
	dN[12] = -0.5 * (1.0 - L1 - L2);
	dN[13] = -0.5 * L1;
	dN[14] = -0.5 * L2;
	dN[15] = 0.5 * (1.0 - L1 - L2);
	dN[16] = 0.5 * L1;
	dN[17] = 0.5 * L2;
}

/***************************************************************************
   GEOSYS/ROCKFLOW - Funktion: GradShapeFunctionPriHQ
        Compute the gradient of shape functions
        for the 6 node prism element
   Arguments:
      E:
        const  double *x  : Pointer to a 3 dimension array for
                            the unit coordinates
      R:
       double * N     : Array of size 18, to store the value of the grandient
                       shape functions
   0--15: dN/dL1
   16--29: dN/dL2
   30--44: dN/dt
   Programming:
   08/2005     WW
   09/2014     WW  Correction
 **************************************************************************/
void GradShapeFunctionPriHQ(double* dN, const double* x)
{
	double L1 = x[0];
	double L2 = x[1];
	const double L0 = 1.0 - L1 - L2;
	double t = x[2];
	double tt1 = 1.0 - t * t;

	//---dN/dL1
	double v1 = (4.0 * L0 - 1);
	double v2 = (4.0 * L1 - 1);
	// Vertex, bottom
	dN[0] = -0.5 * (v1 * (1.0 - t) - tt1);
	dN[1] = 0.5 * (v2 * (1.0 - t) - tt1);
	dN[2] = 0.0;
	// Vertex, top
	dN[3] = -0.5 * (v1 * (1.0 + t) - tt1);
	dN[4] = 0.5 * (v2 * (1.0 + t) - tt1);
	dN[5] = 0.0;
	// Middle point, bottom
	dN[6] = 2.0 * (L0 - L1) * (1.0 - t);
	dN[7] = 2.0 * L2 * (1.0 - t);
	dN[8] = -dN[7];
	// Middle point, top
	dN[9] = 2.0 * (L0 - L1) * (1.0 + t);
	dN[10] = 2.0 * L2 * (1.0 + t);
	dN[11] = -dN[10];
	// Middle point, center
	dN[12] = -tt1;
	dN[13] = tt1;
	dN[14] = 0.0;

	//---dN/dL2
	v1 = (4.0 * L2 - 1);
	// Vertex, bottom
	dN[15] = dN[0];
	dN[16] = 0.0;
	dN[17] = 0.5 * (v1 * (1.0 - t) - tt1);
	// Vertex, top
	dN[18] = dN[3];
	dN[19] = 0.0;
	dN[20] = 0.5 * (v1 * (1.0 + t) - tt1);
	// Middle point, bottom
	dN[21] = -2.0 * L1 * (1.0 - t);
	dN[22] = -dN[21];
	v1 = 2.0 * (L0 - L2);
	dN[23] = v1 * (1.0 - t);
	// Middle point, top
	dN[24] = -2.0 * L1 * (1.0 + t);
	dN[25] = -dN[24];
	dN[26] = v1 * (1.0 + t);
	// Middle point, center
	dN[27] = -tt1;
	dN[28] = 0.0;
	dN[29] = tt1;

	//---dN/dt
	v1 = 2.0 * L0 - 1;
	v2 = 2.0 * L1 - 1;
	double v3 = 2.0 * L2 - 1;
	// Vertex, bottom
	dN[30] = 0.5 * L0 * (-v1 + 2.0 * t);
	dN[31] = 0.5 * L1 * (-v2 + 2.0 * t);
	dN[32] = 0.5 * L2 * (-v3 + 2.0 * t);
	// Vertex, top
	dN[33] = 0.5 * L0 * (v1 + 2.0 * t);
	dN[34] = 0.5 * L1 * (v2 + 2.0 * t);
	dN[35] = 0.5 * L2 * (v3 + 2.0 * t);
	// Middle point, bottom
	dN[36] = -2.0 * L0 * L1;
	dN[37] = -2.0 * L1 * L2;
	dN[38] = -2.0 * L2 * L0;
	// Middle point, top
	dN[39] = -dN[36];
	dN[40] = -dN[37];
	dN[41] = -dN[38];
	// Middle point, center
	dN[42] = -2.0 * L0 * t;
	dN[43] = -2.0 * L1 * t;
	dN[44] = -2.0 * L2 * t;
}

/***************************************************************************
   GEOSYS/ROCKFLOW - Funktion: ShapeFunctionPyra
   Task:
        Compute the gradient of shape functions
        for the 5 node pyramid element
   Arguments:
      E:
        const  double *x  : Pointer to a 3 dimension array for
                            the unit coordinates
      R:
        double * N     : Array of size 5, to store the value of shape
                            function of 5 nodes
   Programming:
   01/2010     NW
 **************************************************************************/
void ShapeFunctionPyra(double* N, const double* x)
{
	const double r = x[0];
	const double s = x[1];
	const double t = x[2];

	N[0] = 0.125 * (1 - r) * (1 - s) * (1 - t);
	N[1] = 0.125 * (1 + r) * (1 - s) * (1 - t);
	N[2] = 0.125 * (1 + r) * (1 + s) * (1 - t);
	N[3] = 0.125 * (1 - r) * (1 + s) * (1 - t);
	N[4] = 0.5 * (1 + t);
}

/***************************************************************************
   GEOSYS/ROCKFLOW - Funktion: ShapeFunctionPyraHQ
   Task:
        Compute the gradient of shape functions
        for the 13 node pyramid element
   Arguments:
      E:
        const  double *x  : Pointer to a 3 dimension array for
                            the unit coordinates
      R:
        double * N     : Array of size 13, to store the value of shape
                            function of 13 nodes
   Programming:
   01/2010     NW
 **************************************************************************/
void ShapeFunctionPyraHQ13(double* N, const double* x)
{
	const double r = x[0];
	const double s = x[1];
	const double t = x[2];

	N[0] = -0.0625 * (1.0 - r) * (1.0 - s) * (1.0 - t)
	       * (4.0 + 3.0 * r + 3.0 * s + 2.0 * r * s + 2.0 * t + r * t + s * t + 2.0 * r * s * t);
	N[1] = -0.0625 * (1.0 + r) * (1.0 - s) * (1.0 - t)
	       * (4.0 - 3.0 * r + 3.0 * s - 2.0 * r * s + 2.0 * t - r * t + s * t - 2.0 * r * s * t);
	N[2] = -0.0625 * (1.0 + r) * (1.0 + s) * (1.0 - t)
	       * (4.0 - 3.0 * r - 3.0 * s + 2.0 * r * s + 2.0 * t - r * t - s * t + 2.0 * r * s * t);
	N[3] = -0.0625 * (1.0 - r) * (1.0 + s) * (1.0 - t)
	       * (4.0 + 3.0 * r - 3.0 * s - 2.0 * r * s + 2.0 * t + r * t - s * t - 2.0 * r * s * t);
	N[4] = 0.5 * t * (1.0 + t);
	N[5] = 0.125 * (1.0 - r * r) * (1.0 - s) * (1.0 - t) * (2.0 + s + s * t);
	N[6] = 0.125 * (1.0 + r) * (1.0 - s * s) * (1.0 - t) * (2.0 - r - r * t);
	N[7] = 0.125 * (1.0 - r * r) * (1.0 + s) * (1.0 - t) * (2.0 - s - s * t);
	N[8] = 0.125 * (1.0 - r) * (1.0 - s * s) * (1.0 - t) * (2.0 + r + r * t);
	N[9] = 0.25 * (1.0 - r) * (1.0 - s) * (1.0 - t * t);
	N[10] = 0.25 * (1.0 + r) * (1.0 - s) * (1.0 - t * t);
	N[11] = 0.25 * (1.0 + r) * (1.0 + s) * (1.0 - t * t);
	N[12] = 0.25 * (1.0 - r) * (1.0 + s) * (1.0 - t * t);
}

/***************************************************************************
   GEOSYS/ROCKFLOW - Funktion: GradShapeFunctionPri
   Task:
        Compute the gradient of shape functions
        for the 6 node prism element
   Arguments:
      E:
        const  double *x  : Pointer to a 3 dimension array for
                            the unit coordinates
      R:
        double * N     : Array of size 18, to store the value of the grandient
                         shape functions
                          0--5: dN/dL1
                         6--11: dN/dL2
                        12--17: dN/dt
   Programming:
   08/2005     WW
 **************************************************************************/
void GradShapeFunctionPyra(double* dN, const double* x)
{
	const double r = x[0];
	const double s = x[1];
	const double t = x[2];
	//  dN/dL1
	dN[0] = -0.125 * (1.0 - s) * (1.0 - t);
	dN[1] = 0.125 * (1.0 - s) * (1.0 - t);
	dN[2] = 0.125 * (1.0 + s) * (1.0 - t);
	dN[3] = -0.125 * (1.0 + s) * (1.0 - t);
	dN[4] = 0.0;
	//  dN/dL2
	dN[5] = -0.125 * (1.0 - r) * (1.0 - t);
	dN[6] = -0.125 * (1.0 + r) * (1.0 - t);
	dN[7] = 0.125 * (1.0 + r) * (1.0 - t);
	dN[8] = 0.125 * (1.0 - r) * (1.0 - t);
	dN[9] = 0.0;
	//  dN/dt
	dN[10] = -0.125 * (1.0 - r) * (1.0 - s);
	dN[11] = -0.125 * (1.0 + r) * (1.0 - s);
	dN[12] = -0.125 * (1.0 + r) * (1.0 + s);
	dN[13] = -0.125 * (1.0 - r) * (1.0 + s);
	dN[14] = 0.5;
}

/***************************************************************************
   GEOSYS/ROCKFLOW - Funktion: GradShapeFunctionPyraHQ
        Compute the gradient of shape functions
        for the 13 node pyramid element
   Arguments:
      E:
        const  double *x  : Pointer to a 3 dimension array for
                            the unit coordinates
      R:
        double * N     : Array of size 18, to store the value of the grandient
                         shape functions
                          0--15: dN/dL1
                         16--29: dN/dL2
                         30--44: dN/dt
   Programming:
   08/2005     WW
 **************************************************************************/
void GradShapeFunctionPyraHQ13(double* dN, const double* x)
{
	const double r = x[0];
	const double s = x[1];
	const double t = x[2];
	//---dN/dr
	dN[0] = 0.0625 * (1.0 - s) * (1.0 - t)
	        * (1.0 + 6.0 * r + s + 4.0 * r * s + t + 2.0 * r * t - s * t + 4.0 * r * s * t);
	dN[1] = -0.0625 * (1.0 - s) * (1.0 - t)
	        * (1.0 - 6.0 * r + s - 4.0 * r * s + t - 2.0 * r * t - s * t - 4.0 * r * s * t);
	dN[2] = -0.0625 * (1.0 + s) * (1.0 - t)
	        * (1.0 - 6.0 * r - s + 4.0 * r * s + t - 2.0 * r * t + s * t + 4.0 * r * s * t);
	dN[3] = 0.0625 * (1.0 + s) * (1.0 - t)
	        * (1.0 + 6.0 * r - s - 4.0 * r * s + t + 2.0 * r * t + s * t - 4.0 * r * s * t);
	dN[4] = 0.0;
	dN[5] = -0.25 * r * (1.0 - s) * (1.0 - t) * (2.0 + s + s * t);
	dN[6] = 0.125 * (1.0 - s * s) * (1.0 - t) * (1.0 - 2.0 * r - t - 2 * r * t);
	dN[7] = -0.25 * r * (1.0 + s) * (1.0 - t) * (2.0 - s - s * t);
	dN[8] = -0.125 * (1.0 - s * s) * (1.0 - t) * (1.0 + 2.0 * r - t + 2 * r * t);
	dN[9] = -0.25 * (1.0 - s) * (1.0 - t * t);
	dN[10] = 0.25 * (1.0 - s) * (1.0 - t * t);
	dN[11] = 0.25 * (1.0 + s) * (1.0 - t * t);
	dN[12] = -0.25 * (1.0 + s) * (1.0 - t * t);

	//---dN/ds
	dN[13] = 0.0625 * (1.0 - r) * (1.0 - t)
	         * (1.0 + r + 6.0 * s + 4.0 * r * s + t - r * t + 2.0 * s * t + 4.0 * r * s * t);
	dN[14] = 0.0625 * (1.0 + r) * (1.0 - t)
	         * (1.0 - r + 6.0 * s - 4.0 * r * s + t + r * t + 2.0 * s * t - 4.0 * r * s * t);
	dN[15] = -0.0625 * (1.0 + r) * (1.0 - t)
	         * (1.0 - r - 6.0 * s + 4.0 * r * s + t + r * t - 2.0 * s * t + 4.0 * r * s * t);
	dN[16] = -0.0625 * (1.0 - r) * (1.0 - t)
	         * (1.0 + r - 6.0 * s - 4.0 * r * s + t - r * t - 2.0 * s * t - 4.0 * r * s * t);
	dN[17] = 0.0;
	dN[18] = -0.125 * (1.0 - r * r) * (1.0 - t) * (1.0 + 2.0 * s - t + 2.0 * s * t);
	dN[19] = -0.25 * (1.0 + r) * s * (1.0 - t) * (2.0 - r - r * t);
	dN[20] = 0.125 * (1.0 - r * r) * (1.0 - t) * (1.0 - 2.0 * s - t - 2.0 * s * t);
	dN[21] = -0.25 * (1.0 - r) * s * (1.0 - t) * (2.0 + r + r * t);
	dN[22] = -0.25 * (1.0 - r) * (1.0 - t * t);
	dN[23] = -0.25 * (1.0 + r) * (1.0 - t * t);
	dN[24] = 0.25 * (1.0 + r) * (1.0 - t * t);
	dN[25] = 0.25 * (1.0 - r) * (1.0 - t * t);

	//---dN/dt
	dN[26] = 0.125 * (1.0 - r) * (1.0 - s) * (1.0 + r + s + 2.0 * t + r * t + s * t + 2.0 * r * s * t);
	dN[27] = 0.125 * (1.0 + r) * (1.0 - s) * (1.0 - r + s + 2.0 * t - r * t + s * t - 2.0 * r * s * t);
	dN[28] = 0.125 * (1.0 + r) * (1.0 + s) * (1.0 - r - s + 2.0 * t - r * t - s * t + 2.0 * r * s * t);
	dN[29] = 0.125 * (1.0 - r) * (1.0 + s) * (1.0 + r - s + 2.0 * t + r * t - s * t - 2.0 * r * s * t);
	dN[30] = 0.5 + t;
	dN[31] = -0.25 * (1.0 - r * r) * (1.0 - s) * (1.0 + s * t);
	dN[32] = -0.25 * (1.0 + r) * (1.0 - s * s) * (1.0 - r * t);
	dN[33] = -0.25 * (1.0 - r * r) * (1.0 + s) * (1.0 - s * t);
	dN[34] = -0.25 * (1.0 - r) * (1.0 - s * s) * (1.0 + r * t);
	dN[35] = -0.5 * (1.0 - r) * (1.0 - s) * t;
	dN[36] = -0.5 * (1.0 + r) * (1.0 - s) * t;
	dN[37] = -0.5 * (1.0 + r) * (1.0 + s) * t;
	dN[38] = -0.5 * (1.0 - r) * (1.0 + s) * t;
}
/***************************************************************************
   GeoSys - Funktion:
           CElement::ComputeDetTri(const *double x1, const *double x2,
                                 const *double x3)
   Aufgabe:
         Compute the vulume of a triangle
   Formalparameter:
           E:
             const *double x1    : Vertex 1
             const *double x2    : Vertex 2
             const *double x3    : Vertex 3
   Programming:
   09/2004     WW        Erste Version
 **************************************************************************/
double ComputeDetTri(const double* x1, const double* x2, const double* x3)
{
	static double u[3], v[3], z[3];

	u[0] = x3[0] - x1[0];
	u[1] = x3[1] - x1[1];
	u[2] = x3[2] - x1[2];

	v[0] = x2[0] - x1[0];
	v[1] = x2[1] - x1[1];
	v[2] = x2[2] - x1[2];

	z[0] = u[1] * v[2] - u[2] * v[1];
	z[1] = u[2] * v[0] - u[0] * v[2];
	z[2] = u[0] * v[1] - u[1] * v[0];

	return 0.5 * sqrt(z[0] * z[0] + z[1] * z[1] + z[2] * z[2]);
}

/***************************************************************************
   GeoSys - Funktion:
           CElem::ComputeDetTet(const *double x1, const *double x2,
                                 const *double x3, const *double x4)
   Aufgabe:
         Compute the vulume of a tedrahedra
   Formalparameter:
           E:
             const *double x1    : Vertex 1
             const *double x2    : Vertex 2
             const *double x3    : Vertex 3
   const *double x4    : Vertex 4
   Programming:
   09/2004     WW        Erste Version
 **************************************************************************/
double ComputeDetTex(const double* x1, const double* x2, const double* x3, const double* x4)
{
	return fabs((x1[0] - x4[0]) * ((x2[1] - x4[1]) * (x3[2] - x4[2]) - (x2[2] - x4[2]) * (x3[1] - x4[1]))
	            - (x1[1] - x4[1]) * ((x2[0] - x4[0]) * (x3[2] - x4[2]) - (x2[2] - x4[2]) * (x3[0] - x4[0]))
	            + (x1[2] - x4[2]) * ((x2[0] - x4[0]) * (x3[1] - x4[1]) - (x2[1] - x4[1]) * (x3[0] - x4[0])))
	       / 6.0;
}

/**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   09/2005 WW Implementation
   12/2010 TF
**************************************************************************/
double NormalizeVector(double* x, size_t n)
{
	double nrm(x[0] * x[0]);
	for (size_t i = 1; i < n; i++)
		nrm += x[i] * x[i];
	double sqrt_nrm(sqrt(nrm));
	for (size_t i = 0; i < n; i++)
		x[i] /= sqrt_nrm;
	return sqrt_nrm;
}

/**************************************************************************
   MSHLib-Method:
   Task:         dim == 3
   Programing:
   09/2005 WW Implementation
**************************************************************************/
void CrossProduction(const double* x, const double* y, double* z)
{
	z[0] = x[1] * y[2] - x[2] * y[1];
	z[1] = x[2] * y[0] - x[0] * y[2];
	z[2] = x[0] * y[1] - x[1] * y[0];
}

/**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   01/2006 YD Implementation
**************************************************************************/
double PointProduction(double* x, double* y)
{
	int i;
	double nrm = 0.0;
	for (i = 0; i < 3; i++)
		nrm += x[i] * y[i];
	return nrm;
}

/**************************************************************************
   MSHLib-Method:
   Task: Langevin function: L(x) = coth(x) - 1/x
   Programing:
   12/2009 NW Implementation
**************************************************************************/
double MLangevin(double v)
{
	double s = 0.0;
	if (v < 0.01)
		s = v * (1.0 / 3.0 + v * v * (-1.0 / 45.0 + 18.0 / 8505.0 * v * v));
	else if (0.01 <= v && v < 20)
		s = (exp(v) + exp(-v)) / (exp(v) - exp(-v)) - 1 / v;
	//    s = coth(v)-1/v;
	else if (20 <= v)
		s = 1.0;

	return s;
}

/**************************************************************************
   MSHLib-Method:
   Task: Vector copy routine
   Programing:
   12.03.2010 JT Implementation
**************************************************************************/
void VCopy(double* x, const double* y, const int n)
{
	int i;
	for (i = 0; i < n; i++)
		x[i] = y[i];
}

/**************************************************************************
   MSHLib-Method:
   Task: Flux limiter function: minmod
   Programing:
   04/2010 NW Implementation
**************************************************************************/
double MinMod(double v1, double v2)
{
	if (v1 * v2 < 0.0)
		return 0.0;
	if (fabs(v1) < fabs(v2))
		return v1;
	else
		return v2;
}

/**************************************************************************
   MSHLib-Method:
   Task: Flux limiter function: Superbee
   Programing:
   04/2010 NW Implementation
**************************************************************************/
double SuperBee(double v1, double v2)
{
	if (v1 * v2 < 0.0)
		return 0.0;
	// max{min{2|a|, |b|},min{|a|, 2|b|}}.
	double a1 = std::min(2.0 * fabs(v1), fabs(v2));
	double a2 = std::min(fabs(v1), 2.0 * fabs(v2));
	double ret = std::max(a1, a2);
	if (v1 > 0.0)
		return ret;
	else
		return -ret;
}

/**************************************************************************
   MSHLib-Method:
   Task: Flux limiter function: Superbee
   Programing:
   04/2010 NW Implementation
**************************************************************************/
double GetFCTADiff(double K_ij, double K_ji)
{
	double r = std::min(0.0, -K_ij);
	r = std::min(r, -K_ji);
	return r;
}
/*##########################################################################
 ##########################################################################
   Ende des ROCKFLOW - Moduls: mathlib.c
 ##########################################################################
 ########################################################################## */

/*!
    \brief Binary search a array
     The array must be sorted

     \param arr     an array
     \param target  searched index
     \param start   the start index of the array
     \parm  end     the end index of the array

     By WW. 03.2011

 */
long binarySearch(long* arr, long target, long start, long end)
{
	long middle;
	while (start <= end)
	{
		middle = (start + end) / 2;
		if (arr[middle] == target)
			return middle;
		else if (arr[middle] > target)
			end = middle - 1;
		else
			start = middle + 1;
	}
	return -1;
}
