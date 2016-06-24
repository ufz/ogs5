/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/***************************************************************************
   ROCKFLOW - Modul: mathlib.h

   Aufgabe:
   ROCKFLOW-Schnittstelle fuer alle mathematischen Operationen, die nicht
   Standard (nicht in math.h enthalten) sind.
   mathlib.h benutzt externe Bibliotheken, die bei Bedarf ausgetauscht
   werden koennen, ohne die ROCKFLOW-Schnittstelle aendern zu muessen.

 **************************************************************************/

#ifndef mathlib_INC

#define mathlib_INC
/* Schutz gegen mehrfaches Einfuegen */

#include <cstddef>

#define noTESTMATH

/* Andere oeffentlich benutzte Module */
//#include "test.h"

/* Die Schnittstellen der Gleichungsloeser und der Speichertechnik werden
   sozusagen durchgeschleift: */
/* #include "matrix_routines.h" */
/* Speichertechnik fuer Matrix des Gesamtgleichungssystems */
#ifndef NEW_EQS // WW. 11.2008
#include "solver.h"
#endif
/* Iterative GLS-Loeser auf Basis der Speichertechnik aus 'matrix.h' (herkoemmliche Verfahren) */

// C++
//#include <vector>

/*##########################################################################
      Mathematische Funktionen
 ########################################################################*/

#ifdef obsolete // 01.2011 WW
extern int MGleichDouble(double zahl1, double zahl2, double tol);
/*   MGleichDouble         - Vergleicht zwei double-Zahlen unter
                             Beruecksichtigung einer Fehlertoleranz */
extern int MOmega1D(double* vf, double r);
/*  Berechnet 1D Ansatzfunktionen */
extern int MOmega2D(double* vf, double r, double s);
/*  Berechnet 2D Ansatzfunktionen */
extern int MOmega3D(double* vf, double r, double s, double t);
/*  Berechnet 3D Ansatzfunktionen */
//#ifdef obsolete //WW. 11.2008
extern int MOmega3DTetrahedron(double* vf, double xx, double yy, double zz, long number);
extern int MOmega2DTriangle(double* vf, double r, double s, long number);
//#endif
/*  Berechnet 2D Ansatzfunktionen fuer Dreiecke */
// extern int MPhi2D(double *vf,double r, double s);
/*  Berechnet 2D Testfunktionen */
extern int MPhi2D_SUPG(double* vf, double r, double s, double* alpha);
/*  Berechnet 2D Testfunktionen (SUPG) */
extern int MPhi3D(double* vf, double r, double s, double t);
/*  Berechnet 3D Testfunktionen */
extern int MPhi3D_SUPG(double* vf, double r, double s, double t, double* alpha);
/*  Berechnet 3D Testfunktionen (SUPG)*/
extern int MGradOmega2D(double* vf, double r, double s);
/*  Berechnet Gradient der 2D Ansatzfunktionen */
extern int MGradOmega3D(double* vf, double r, double s, double t);
/*  Berechnet Gradient der 3D Ansatzfunktionen */
extern int MGradPhi2D(double* vf, double r, double s);
/*  Berechnet Gradient der 2D Testfunktionen */
extern int MGradPhi3D(double* vf, double r, double s, double t);
/* Faktoren fuer die X Punkt Gauss-Integration */
extern void MGetCoor(int typ, long j, double* r, double* s, double* t);
/* lokale Koordinaten der Eckpunkte */

/*   MBtrgVec              - Betrag von Vektor */

/*   MAngleVectors         -  Winkel zwischen Vektoren */ /* MB */
double MAngleVectors(double* v1, double* v2);
/*   MNormiere             - Normiert Vektoren */
void MNormiere(double* vec, long n);
/*   M2Determinante        - Determinante einer 2x2Matrix */
extern double M2Determinante(double* matrix);
/*   M3Determinante        - Determinante einer 3x3Matrix */
// extern double M3Determinante ( double *m );
/*   Mxg2Determinante      - Determinante einer beliebigen Matrix
                             mit goesseren Ausmassen als 2x2  */
/*   M4Determinante        - Determinante einer 4x4Matrix */
// extern double M4Determinante ( double *m );
extern double Mxg2Determinante(double* matrix, long m, long n);
/*   MTranspoVec           - Transponieren beliebiger Vektoren */
extern void MTranspoVec(double* vec, long g);
/*   MTranspoMat           - Transponieren beliebiger Matrizen */
extern void MTranspoMat(double* mat1, long m, long n, double* mat2);
/*   M2Invertiere          - Invertiert 2x2 Matrizen */
extern void M2Invertiere(double* m);
/*   M2InvertiereUndTransponiere - Invertiert und transponiert 2x2 Matrizen */
extern void M2InvertiereUndTransponiere(double* m);
/*   M3Invertiere          - Invertiert 3x3 Matrizen */
extern void M3Invertiere(double* m);
/*   MInvertiere           - Invertieren beliebier regulaerer Matrizen */
extern void MInvertiere(double* matrix, long m, long n);
/*   MAddVektoren          - Addition zweier beliebier Vektoren */
extern int MAddVektoren(double* v1, double* v2, double* vout, long g);
/*   MAddSkalVektoren      - Vektoren mit Skalar multiplizieren und dann addieren */
// WW extern int MAddSkalVektoren(double *v1, double m1, double *v2, double m2, double *vout, long g);
/*   MAddMatrizen          - Addition zweier beliebier Matrizen */
extern int MAddMatrizen(double* m1, double* m2, double* mout, long m, long n);
/*   MMultVecSkalar        - Multiplikation Vektor mit Skalarwert */
extern int MMultVecSkalar(double* vec, double skal, long g);
/*   MMultMatSkalar        - Multiplikation Matrix mit Skalarwert */
extern int MMultMatSkalar(double* matrix, double skal, long m, long n);
/*   MSkalarprodukt        - Skalarprodukt zweier beliebiger Vektoren */
// extern double MSkalarprodukt ( double *vec1, double *vec2, long g );
/*   M3KreuzProdukt        - Kreuzprodukt von 3D-Vektoren */
// extern int M3KreuzProdukt( double *vec1, double *vec2 , double *vec);
/*##########################################################################
   Auf Matrizen und Vektoren aufbauende Funktionen
 ########################################################################*/

int MKTF2Dr2D(double* vec1, double winkel, double* vec2);
int MKTFMat2Dr2D(double winkel, double* tmat);
int MKTF2Dt2D(double* vec1, double dx, double dy, double* vec2);
int MKTF3Dr2D(double* vec1, double* vec2, double winkel, double* vec);
void MKTFMat3Dr2D(double* vec1, double* vec2, double winkel, double* mat);
/* was ist das alles ??? (msr) */

/*##########################################################################
    Prueffunktion fuer Matrix
 ########################################################################*/
/*   MBistDuDiagMat        - Prueft auf Diagonalitaet einer Matrix */
extern int MBistDuDiagMat(double* matrix, long m, long n);

/*##########################################################################
   Sortierfunktionen
 ########################################################################*/
extern void MQSort_LongDouble(void* DataSets, const int NumberOfDataSets, const int SiezeOfDataSet);
extern int MCompare_for_MQSort_LongDouble(const void* arg1, const void* arg2);

extern double* TensorDrehDich(double* d, double* velo);
extern int GetPriMatFromTriMat(double* mat1, double* mat2);
extern int MMultMatMat2(double* mat1, long m1, long n1, double* mat2, double* ergebnis);

/*##########################################################################
   eventuell noch nuetzliche Funktionen,
   die im Moment nicht gebraucht werden
 ########################################################################*/
extern double* JWDMMultMatSkalar(double* matrix, double skal, long m, long n);

// extern long IsNodeInsideTriangle (long n1, long n2, long n3, long node);
// CC
extern double CalcTriangleArea(long n1, long n2, long n3);
// extern long IsNodeInsidePlain (long n1, long n2, long n3, long node);

extern int MOmega2D_9N(double* vf, double r, double s);
extern int MPhi2D_9N(double* vf, double r, double s);
extern int MGradOmega2D_9N(double* vf, double r, double s);
extern int MGradPhi2D_9N(double* vf, double r, double s);

extern int MOmega3D_20N(double* vf, double r, double s, double t);
extern int MPhi3D_20N(double* vf, double r, double s, double t);
extern int MGradOmega3D_20N(double* vf, double r, double s, double t);
extern int MGradPhi3D_20N(double* vf, double r, double s, double t);
#endif // 05.03.2010. WW //#ifdef obsolete

/* Ermittelt min */
extern double MMin(double, double);
/* Ermittelt max */
extern double MMax(double, double);
/* Ermittelt Intervall */
extern double MRange(double a, double b, double c);
double MBtrgVec(double* vec, long n);

#ifndef NEW_EQS // WW. 05.03.2010
/*##########################################################################
    Funktionen fuer Gleichungsloeser (CG)
 ########################################################################*/
extern double MVekNorm1(double* x, long n);
/* Spaltensummennorm */
extern double MVekNorm2(double* x, long n);
/* Euklidische Norm */
// WW extern double MVekNormMax ( double *x, long n );
/* Maximumnorm */
extern void MVekSum(double* x, double alpha, double* y, long n);
/* Fuehrt die Operation x = x + alpha * y durch; n: Vektordimension */
extern void MVekGle(double alpha, double* x, double beta, double* y, double* z, long n);
/* Fuehrt die Operation z = alpha * x + beta * y durch; n: Vektordimension */
extern double MVekDist(double* x, double* y, long n);
/* Abstand zwischen zwei Vektoren */

/*##########################################################################
    Bearbeitungsfunktionen
 ########################################################################*/
/*   MMachVec              - Erzeugt einen Vektor */
extern double* MMachVec(long g);
/*   MNullVec              - Fuellt einen Vektor mit Nullen */
extern void MNullVec(double*, long);
/*   MLoeschVec            - Loescht einen Vektor */
// WW extern void MLoeschVec (double *vec);
/*   MKopierVec            - Kopiert einen Vektor auf einen Anderen */
extern void MKopierVec(double* vecquelle, double* vecziel, long g);
/*   MMachMat              - Erzeugt einer Matrix */
// WW extern double *MMachMat (long m, long n);
/*   MNullMat              - Fuellt eine Matrix mit Nullen*/
// WW extern void MNullMat (double *,long , long);
/*   MLoeschMat            - Loescht einer Matrix */
// WW extern void MLoeschMat (double *mat);
/*   MKopierMat            - Kopiert eine Matrix auf eine Andere */
extern void MKopierMat(double* matq, double* matz, long m, long n);
/* nur fuer mich!!! hh */
// extern void MZeigVec (double *vec, long grad, char *text);
// extern void MZeigMat (double *mat, long m, long n, char *text);
// extern void M2FileVec (double *vec, long grad, char *text);
/* PS: Finger Weg, Michael!!! --> baeeh */
/*   MAddSkalVektoren      - Vektoren mit Skalar multiplizieren und dann addieren */
extern int MAddSkalVektoren(double* v1, double m1, double* v2, double m2, double* vout, long g);
#endif ///////////////////////////

// WW: Only unsed in a member function of CFiniteElementStd implementated by MB
/*   MMultVecVec           - Multiplikation Vektor mit Vektor */
extern int MMultVecVec(double* vec1, long gv1, double* vec2, long gv2, double* mato, long mo, long no);
/*   MMultVecMat           - Multiplikation Vektor mit Matrix */
extern int MMultVecMat(double* vec, long gv, double* mat, long m, long n, double* veco, long go);
/*   MMultMatVec           - Multiplikation Matrix mit Vektor */
extern int MMultMatVec(const double* mat, long m, long n, double* vec, long g, double* veco, long r);
/*   MMultMatMat           - Multiplikation Matrix mit Matrix */
extern int MMultMatMat(double* mat1, long m1, long n1, double* mat2, long m2, long n2, double* mato, long mo, long no);
/*##########################################################################
   Geometrie-Funktionen
 ########################################################################*/
// extern double MCalcDistancePointToPoint(double *pt1,double *pt2);
extern double MCalcDistancePointToLine(double* pt, double* l1, double* l2);
extern double MCalcProjectionOfPointOnLine(double* pt1, double* pt2, double* pt3, double* pt4);
extern double MCalcDistancePointToPlane(double const* const pt, double* e1, double* e2, double* e3);
// extern double MCalcProjectionOfPointOnPlane(double *pt, double *e1, double *e2, double *e3, double *proj);

/*   MNulleVec             - Setze angegebenen Vektor = 0.0 */
extern void MNulleVec(double* vec, long g);
/*   MNulleMat             - Setze angegebene Matrix = 0.0 */
extern void MNulleMat(double* vec, long m, long n);

/*  Berechnet Gradient der 3D Testfunktionen */
extern double MXPGaussPkt(long grd, long pkt);
/* Punkte fuer die X Punkt Gauss-Integration */
extern double MXPGaussFkt(long grd, long pkt);

extern void realCoordTriHQ(double* x, const double* XY, const double* u);

// Family of  element interpolation. WW
extern void ShapeFunctionLine(double* N1, const double* u);
extern void ShapeFunctionLineHQ(double* N1, const double* u);
extern void SamplePointTriHQ(const int nsample, double* SPoints);
extern void SamplePointTet5(const int nsample, double* SPoints);
extern void SamplePointTet15(const int nsample, double* SPoints);
extern void SamplePointPyramid5(const int nsample, double* SPoints);
extern void SamplePointPyramid8(const int i, double* SPoint);
extern void SamplePointPyramid13(const int nsample, double* SPoints);
extern void ShapeFunctionTri(double* N3, const double* u);
extern void ShapeFunctionTriHQ(double* N6, const double* u);
extern void ShapeFunctionQuad(double* N4, const double* u);
extern void ShapeFunctionQuadHQ8(double* N8, const double* u);
extern void ShapeFunctionQuadHQ(double* N9, const double* u);
extern void ShapeFunctionTet(double* Nt4, const double* u);
extern void ShapeFunctionTetHQ(double* N10, const double* u);
extern void ShapeFunctionHex(double* N8, const double* x);
extern void ShapeFunctionHexHQ(double* N9, const double* u);
extern void ShapeFunctionPri(double* N, const double* x);
extern void ShapeFunctionPriHQ(double* N, const double* u);
extern void ShapeFunctionPyra(double* N, const double* x);
extern void ShapeFunctionPyraHQ13(double* N, const double* u);
// Gradient of ...
extern void GradShapeFunctionLine(double* dN1, const double* u);
extern void GradShapeFunctionLineHQ(double* dN1, const double* u);
extern void GradShapeFunctionTri(double* dN3, const double* u);
extern void GradShapeFunctionTriHQ(double* dN3, const double* u);
extern void GradShapeFunctionQuad(double* dN4, const double* u);
extern void GradShapeFunctionQuadHQ(double* dN9, const double* u);
extern void GradShapeFunctionQuadHQ8(double* dN8, const double* u);
extern void GradShapeFunctionTet(double* dNt4, const double* u);
extern void GradShapeFunctionTetHQ(double* dN10, const double* u);
extern void GradShapeFunctionHex(double* N8, const double* x);
extern void GradShapeFunctionHexHQ(double* dN9, const double* u);
extern void GradShapeFunctionPri(double* dN, const double* x);
extern void GradShapeFunctionPriHQ(double* dN, const double* u);
extern void GradShapeFunctionPyra(double* dN, const double* x);
extern void GradShapeFunctionPyraHQ13(double* dN, const double* u);

extern double pai;

extern long binarySearch(long* arr, long target, long start, long end);
// WW Cubic spline
// WW
double ComputeDetTri(const double* x1, const double* x2, const double* x3);
double ComputeDetTex(const double* x1, const double* x2, const double* x3, const double* x4);
void CrossProduction(const double* x, const double* y, double* z);
double NormalizeVector(double* x, size_t n);

// extern double MVectorlength(double dx, double dy, double dz);
extern double PointProduction(double* x, double* y);
extern void VCopy(double* x, const double* y, const int n);

// NW
extern double MLangevin(double v);
extern double MinMod(double v1, double v2);
extern double SuperBee(double v1, double v2);
extern double GetFCTADiff(double v1, double v2);

#endif /* gehoert zum Schutz gegen mehrfaches Einfuegen */

/*##########################################################################
    Ende von: ROCKFLOW - Modul mathlib.h
 ########################################################################*/
