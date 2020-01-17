/**
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   GeoLib - Object:GEO mathlib
   Task:
   Programing:
   08/2005 CC Implementation
**************************************************************************/

#include "geo_mathlib.h"
#include "mathlib.h"  //WW
#include <cstdio>
#include <stdlib.h>

/**************************************************************************/
/* GEO MathLib - Funktion: EuklVek3dDist
 */
/* Aufgabe:
   Berechnet den Abstand zwischen x und y in der Euklidischen Norm
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E double *x : Zeiger auf Vektor (Punkt x)
   E double *y : Zeiger auf Vektor (Punkt y)
   E long n : Dimension von x und y
 */
/* Ergebnis:
   Abstand
 */
/* Programmaenderungen:
   09/1997     AH         Einbau der Funktion
 */
/**************************************************************************/
double EuklVek3dDist(double* x, double* y)
{
    return EuklVek3dDistCoor(x[0], x[1], x[2], y[0], y[1], y[2]);
}

/**************************************************************************/
/* GEO MathLib - Funktion: EuklVek3dDistCoor
 */
/* Aufgabe:
   Berechnet den Abstand zwischen x und y in der Euklidischen Norm
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E double x1,y1,z1,x2,y2,z2 : Koordinaten des 1. Punktes
                                bzw. des 2. Punktes
   E double *y : Zeiger auf Vektor
   E long n : Dimension von x und y
 */
/* Ergebnis:
   Abstand
 */
/* Programmaenderungen:
   09/1997     AH         Einbau der Funktion
 */
/**************************************************************************/

double EuklVek3dDistCoor(double x1, double y1, double z1, double x2, double y2,
                         double z2)
{
    return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) +
                (z2 - z1) * (z2 - z1));
}
/**************************************************************************/
/* GEO MathLib - Funktion: Vek3dDistCoor
 */
/* Aufgabe:
   Berechnet den Abstand zwischen x und y in der Norm 0,1 und 2
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E double x1,y1,z1,x2,y2,z2 : Koordinaten des 1. Punktes
                                bzw. des 2. Punktes
   E double *y : Zeiger auf Vektor
   E long n : Dimension von x und y
 */
/* Ergebnis:
   Abstand
 */
/* Programmaenderungen:
   08/1998     AH         Einbau der Funktion
 */
/**************************************************************************/
double Vek3dDistCoor(double x1, double y1, double z1, double x2, double y2,
                     double z2, int norm)
{
    double d;

    switch (norm)
    {
        case 0:
            d = fabs(x2 - x1);
            if (fabs(y2 - y1) > d)
                d = fabs(y2 - y1);
            if (fabs(z2 - z1) > d)
                d = fabs(z2 - z1);
            return d;

        case 1:
            return fabs(x2 - x1) + fabs(y2 - y1) + fabs(z2 - z1);

        case 2:
            return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) +
                        (z2 - z1) * (z2 - z1));

        default:
            return -1.;
    }
}
/**************************************************************************
   GEO MathLib - Funktion: M3KreuzProdukt
   Aufgabe:
           Kreuzprodukt zweier 3D-Vektoren
   Formalparameter:
           E: *vec1
           E: *vec2
           A: *ec
   Ergebnis:
           In Vektor vec1 wird das Ergebnis gespeichert.
   Aenderungen/Korrekturen:
   08/1994     hh        Erste Version
**************************************************************************/

int M3KreuzProdukt(double* vec1, double* vec2, double* vec)
{
    vec[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    vec[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    vec[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
    return 1;
}

/***************************************************************************
   GEO MathLib - Funktion: MSkalarprodukt
   Aufgabe:
           Berechnung des Skalarprodukts aus zwei beliebig grossen
           Vektoren
   Formalparameter:
           E: *vecX - Zeiger auf Vektoren
           E: g - Ausdehnung bzw. Grad der Vektoren
   Ergebnis:
           Skalarprodukt
   Aenderungen/Korrekturen:
   07/1994    hh        Erste Version
   11/1995    msr       register
    8/0001    C.Thorenz CBLAS-eingehaengt
 **************************************************************************/

double MSkalarprodukt(double* vec1, double* vec2, long g)
{
#ifndef CBLAS_MSkalarprodukt
    register long i;
    register double sammy = 0.0;
#ifdef SX
#pragma cdir nodep
#endif
    for (i = 0l; i < g; i++)
        sammy += vec1[i] * vec2[i];
    return sammy;
#else
    return cblas_ddot(g, vec1, 1, vec2, 1);
#endif
} /* extern double MSkalarprodukt */
/***************************************************************************
   GEO MathLib - Funktion: M4Determinante
   Aufgabe:
           Berechnung der Determinante einer 4x4-Matrix
           nach der Regel von Sarrus
   Formalparameter:
           E: *matrix
   Ergebnis:
           Determinante
   Aenderungen/Korrekturen:
   10/2001   OK   Erste Version
    08/2005 CC Move
 **************************************************************************/
double M4Determinante(double* m)
{
    double determinante = 0.0;
    double A11[9], A12[9], A13[9], A14[9];

    A11[0] = m[5];
    A11[1] = m[6];
    A11[2] = m[7];
    A11[3] = m[9];
    A11[4] = m[10];
    A11[5] = m[11];
    A11[6] = m[13];
    A11[7] = m[14];
    A11[8] = m[15];

    A12[0] = m[1];
    A12[1] = m[2];
    A12[2] = m[3];
    A12[3] = m[9];
    A12[4] = m[10];
    A12[5] = m[11];
    A12[6] = m[13];
    A12[7] = m[14];
    A12[8] = m[15];

    A13[0] = m[1];
    A13[1] = m[2];
    A13[2] = m[3];
    A13[3] = m[5];
    A13[4] = m[6];
    A13[5] = m[7];
    A13[6] = m[13];
    A13[7] = m[14];
    A13[8] = m[15];

    A14[0] = m[1];
    A14[1] = m[2];
    A14[2] = m[3];
    A14[3] = m[5];
    A14[4] = m[6];
    A14[5] = m[7];
    A14[6] = m[9];
    A14[7] = m[10];
    A14[8] = m[11];

    determinante = +1.0 * M3Determinante(A11) - 1.0 * M3Determinante(A12) +
                   1.0 * M3Determinante(A13) - 1.0 * M3Determinante(A14);

    return determinante;
}
double CalcTetraederVolume(double* x, double* y, double* z)
{
    static double mat4x4[16];
    mat4x4[0] = 1.0;
    mat4x4[1] = x[0];
    mat4x4[2] = y[0];
    mat4x4[3] = z[0];
    mat4x4[4] = 1.0;
    mat4x4[5] = x[1];
    mat4x4[6] = y[1];
    mat4x4[7] = z[1];
    mat4x4[8] = 1.0;
    mat4x4[9] = x[2];
    mat4x4[10] = y[2];
    mat4x4[11] = z[2];
    mat4x4[12] = 1.0;
    mat4x4[13] = x[3];
    mat4x4[14] = y[3];
    mat4x4[15] = z[3];
    return fabs(M4Determinante(mat4x4)) / 6.;
}
double CalcPyramidVolume(double* x, double* y, double* z)
{
    /*
       x[0]=1.0; y[0]=-1.0; z[0]=0.0;
       x[1]=1.0; y[1]=1.0; z[1]=0.0;
       x[2]=-1.0; y[2]=1.0; z[2]=0.0;
       x[3]=-1.0; y[3]=-1.0; z[3]=0.0;
       x[4]=0.0; y[4]=0.0; z[4]=1.0;
     */
    // double hight
    double volume;

    double p1[3], p2[3], p3[3], p5[3], proj[3];
    //--------------------------------------------------------------
    // tet version
    p1[0] = x[0];
    p1[1] = y[0];
    p1[2] = z[0];
    p2[0] = x[1];
    p2[1] = y[1];
    p2[2] = z[1];
    p3[0] = x[2];
    p3[1] = y[2];
    p3[2] = z[2];
    p5[0] = x[4];
    p5[1] = y[4];
    p5[2] = z[4];

    // WW hight =
    MCalcProjectionOfPointOnPlane(p5, p1, p2, p3, proj);
    double xt[4], yt[4], zt[4];
    double volume1, volume2, volume3, volume4;
    xt[0] = x[0];
    xt[1] = x[1];
    xt[2] = x[4];
    xt[3] = proj[0];
    yt[0] = y[0];
    yt[1] = y[1];
    yt[2] = y[4];
    yt[3] = proj[1];
    zt[0] = z[0];
    zt[1] = z[1];
    zt[2] = z[4];
    zt[3] = proj[2];
    volume1 = CalcTetraederVolume(xt, yt, zt);
    xt[0] = x[1];
    xt[1] = x[2];
    xt[2] = x[4];
    xt[3] = proj[0];
    yt[0] = y[1];
    yt[1] = y[2];
    yt[2] = y[4];
    yt[3] = proj[1];
    zt[0] = z[1];
    zt[1] = z[2];
    zt[2] = z[4];
    zt[3] = proj[2];
    volume2 = CalcTetraederVolume(xt, yt, zt);
    xt[0] = x[2];
    xt[1] = x[3];
    xt[2] = x[4];
    xt[3] = proj[0];
    yt[0] = y[2];
    yt[1] = y[3];
    yt[2] = y[4];
    yt[3] = proj[1];
    zt[0] = z[2];
    zt[1] = z[3];
    zt[2] = z[4];
    zt[3] = proj[2];
    volume3 = CalcTetraederVolume(xt, yt, zt);
    xt[0] = x[3];
    xt[1] = x[0];
    xt[2] = x[4];
    xt[3] = proj[0];
    yt[0] = y[3];
    yt[1] = y[0];
    yt[2] = y[4];
    yt[3] = proj[1];
    zt[0] = z[3];
    zt[1] = z[0];
    zt[2] = z[4];
    zt[3] = proj[2];
    volume4 = CalcTetraederVolume(xt, yt, zt);
    volume = fabs(volume1) + fabs(volume2) + fabs(volume3) + fabs(volume4);
    return volume;
}
double CalcPrismVolume(double* x, double* y, double* z)
{
    double tet1, tet2, tet3, prism;
    double mat4x4[16];
    mat4x4[0] = 1.0;
    mat4x4[1] = x[0];
    mat4x4[2] = y[0];
    mat4x4[3] = z[0];
    mat4x4[4] = 1.0;
    mat4x4[5] = x[1];
    mat4x4[6] = y[1];
    mat4x4[7] = z[1];
    mat4x4[8] = 1.0;
    mat4x4[9] = x[2];
    mat4x4[10] = y[2];
    mat4x4[11] = z[2];
    mat4x4[12] = 1.0;
    mat4x4[13] = x[3];
    mat4x4[14] = y[3];
    mat4x4[15] = z[3];
    tet1 = fabs(M4Determinante(mat4x4)) / 6.;
    mat4x4[0] = 1.0;
    mat4x4[1] = x[1];
    mat4x4[2] = y[1];
    mat4x4[3] = z[1];
    mat4x4[4] = 1.0;
    mat4x4[5] = x[2];
    mat4x4[6] = y[2];
    mat4x4[7] = z[2];
    mat4x4[8] = 1.0;
    mat4x4[9] = x[3];
    mat4x4[10] = y[3];
    mat4x4[11] = z[3];
    mat4x4[12] = 1.0;
    mat4x4[13] = x[4];
    mat4x4[14] = y[4];
    mat4x4[15] = z[4];
    tet2 = fabs(M4Determinante(mat4x4)) / 6.;
    mat4x4[0] = 1.0;
    mat4x4[1] = x[2];
    mat4x4[2] = y[2];
    mat4x4[3] = z[2];
    mat4x4[4] = 1.0;
    mat4x4[5] = x[3];
    mat4x4[6] = y[3];
    mat4x4[7] = z[3];
    mat4x4[8] = 1.0;
    mat4x4[9] = x[4];
    mat4x4[10] = y[4];
    mat4x4[11] = z[4];
    mat4x4[12] = 1.0;
    mat4x4[13] = x[5];
    mat4x4[14] = y[5];
    mat4x4[15] = z[5];
    tet3 = fabs(M4Determinante(mat4x4)) / 6.;
    prism = tet1 + tet2 + tet3;
    return prism;
}
/***************************************************************************
   GEO MathLib - Funktion:  MCalcProjectionOfPointOnPlane

   Aufgabe:
           Ermittelt den Projektionspunkt eines Punktes auf einer Flaeche
 (Lotpunkt)

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
double MCalcProjectionOfPointOnPlane(double* pt, double* e1, double* e2,
                                     double* e3, double* proj)
{
    int i;
    double vec1[3], vec2[3], vec3[3], normal[3], abstand, volume, area;

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
    if (area < MKleinsteZahlen)
    {
        printf("\n FEHLER IN CalcDistancePointToPlane: Flaeche <= Null !!!");
        exit(1);
    }
#endif

    /* Abstand zur Ebene */
    abstand = volume / (area + MKleinsteZahlen);

    /* Normierter Normalen-Vektor  */
    for (i = 0; i < 3; i++)
        normal[i] /= (area + MKleinsteZahlen);

    if (MSkalarprodukt(normal, vec3, 3) > 0.)
        for (i = 0; i < 3; i++)
            proj[i] = pt[i] - fabs(abstand) * normal[i];
    else
        for (i = 0; i < 3; i++)
            proj[i] = pt[i] + fabs(abstand) * normal[i];

    return abstand;
}

/***************************************************************************
   GEO MathLib - Funktion: M3Determinante
   Aufgabe:
           Berechnung der Determinante einer 3x3-Matrix
           nach der Regel von Sarrus
   Formalparameter:
           E: *matrix
   Ergebnis:
           Determinante
   Aenderungen/Korrekturen:
   07/1994     hh        Erste Version
 **************************************************************************/

double M3Determinante(double* m)
{
    return m[0] * m[4] * m[8] + m[1] * m[5] * m[6] + m[2] * m[3] * m[7] -
           m[2] * m[4] * m[6] - m[0] * m[5] * m[7] - m[1] * m[3] * m[8];
}

double MCalcDistancePointToPoint(double* pt1, double* pt2)
{
    double vec[3];

    vec[0] = pt1[0] - pt2[0];
    vec[1] = pt1[1] - pt2[1];
    vec[2] = pt1[2] - pt2[2];

    return MBtrgVec(vec, 3);
}

/**************************************************************************
   ROCKFLOW - Function: TOLSortNodes1

   Task:
   Sort nodes descending according to the criterium.

   Parameter: (I: Input; R: Return; X: Both)
           I: long* nodes, double* criterium, int anz

   Return:
  *long nodes (aber sortiert!)

   Programming:
   09/2002   MB   First Version
**************************************************************************/
long* TOLSortNodes1(long* nodes, double* criterium, int anz)
{
    int flag = 1;
    int i;
    int nummer = 0;
    long tempnode;
    double temp;

    do
    {
        flag = 0;
        nummer++;
        for (i = 0; i < (anz - nummer); i++)
            if (criterium[i] < criterium[i + 1])
            {
                flag = 1;
                tempnode = nodes[i];
                temp = criterium[i];
                nodes[i] = nodes[i + 1];
                criterium[i] = criterium[i + 1];
                nodes[i + 1] = tempnode;
                criterium[i + 1] = temp;
            } /* end if */
              /* end for */
    } while (flag == 1);
    return nodes;
}
/***************************************************************************
   ROCKFLOW - Funktion: MPhi2D

   Aufgabe: Berechnet Testfunktion (r,s) ohne Upwind-Parameter
            alpha (im Unterschied zu MPhi2D_SUPG)
            Phi = Omega

   Formalparameter: Z: *vf - 1x4 Feld
                    E: r,s (s.o.)

   Ergebnis: Vektor

   Aenderungen/Korrekturen:
   04/1995     R.Kaiser        Erste Version

 **************************************************************************/

int MPhi2D(double* vf, double r, double s)
{
    int i;
    int ok = 0;
    vf[0] = (1.0 + r) * (1.0 + s);
    vf[1] = (1.0 - r) * (1.0 + s);
    vf[2] = (1.0 - r) * (1.0 - s);
    vf[3] = (1.0 + r) * (1.0 - s);
    for (i = 0; i < 4; i++)
        vf[i] *= 0.25;
    return ok = 1;
}

#ifdef RFW_FRACTURE
/*************************************************************************
   ROCKFLOW - Funktion: LineSegmentIntersection
   Aufgabe:
   Find the intersection of 2 line segments, if it exists.  2D only.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E xline1:  x-coords at the ends of first line segment
   E yline1:  y-coords at the ends of first line segment
   E xline2:  x-coords at the ends of second line segment
   E yline2:  y-coords at the ends of second line segment
   R intercept: the intersection point of the 2 segments

   Ergebnis:
   Returns the value true if an intercept exists
   Programmaenderungen:
   05/2005 RFW Implementierung
 ***************************************************************************/
bool LineSegmentIntersection(vector<double> xline1, vector<double> yline1,
                             vector<double> xline2, vector<double> yline2,
                             vector<double>& intercept)
{
    double determinant, t1;  // WW, t2;
    bool crosses = false;

    determinant = ((yline2[1] - yline2[0]) * (xline1[1] - xline1[0]) -
                   (xline2[1] - xline2[0]) * (yline1[1] - yline1[0]));
    if (determinant == 0)
        crosses = false;
    else
    {
        t1 = ((xline2[1] - xline2[0]) * (yline1[0] - yline2[0]) -
              (yline2[1] - yline2[0]) * (xline1[0] - xline2[0])) /
             determinant;

        if (t1 < 0 || t1 > 1)
            crosses = false;
        else
        {
            crosses = true;
            intercept.push_back(xline1[0] + t1 * (xline1[1] - xline1[0]));
            intercept.push_back(yline1[0] + t1 * (yline1[1] - yline1[0]));
        }
    }

    return crosses;
}
#endif
