/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   ROCKFLOW - Modul: matrix.h

   Aufgabe:
   Verwaltung von nxn - Matrizen mit verschiedenen Speichertechniken und
   Bereitstellung entsprechender Zugriffsfunktionen sowie einfacher
   mathematischer Operationen auf die entsprechende Matrix.
   Benoetigte Vektoren werden als eindimensionale double-Felder erwartet.

   Speichertechniken:  (1) vollbesetzte Matrix
                           ( param1 = Dimension )
   (2) nur A[i,j]!=0 werden gespeichert (Sparse)
   ( param1 = Dimension )
   (3) symmetrische sparse Matrix
   fuer Preconditioner "incomplete Cholesky"
   nur Aik !=0 mit k>=i werden gespeichert
   ( param1 = Dimension )
   (4) unsymmetrische sparse Matrix
   fuer Preconditioner "incomplete LDU-Zerlegung"
   nur Aik !=0 werden gespeichert
   ( param1 = Dimension )

   Benutzungsbeispiel:

   long dummy =0;             Dummy-Variable mit beliebigem Wert
   void *A = NULL;            Zeiger A fuer Matrix-Wurzel definieren
   void *B = NULL;            Zeiger B fuer Matrix-Wurzel definieren
   A =MXSetupMatrix (10l, 1, dummy);  10x10-Matrix erzeugen, Modell 1
   B =MXSetupMatrix (100l, 2, dummy);  100x100-Matrix erzeugen, Modell 2

   MXSetMatrixPointer(A);     Aktuelle Matrix ist A
   MXInitMatrix();            Matrix A  mit Nullen fuellen
   MXSet(0,0,3.1415926);      Beispiel: A[0,0] = Pi;
   printf("%f",MXGet(0,0));   Beispiel: Pi wieder ausgeben
   A =MXDestroyMatrix();      Speicher wieder freigeben, A =NULL setzen

   MXSetMatrixPointer(B);     Aktuelle Matrix ist B
   MXInitMatrix();            Matrix B  mit Nullen fuellen
   MXSet(2,3,3.1415926);      Beispiel: B[2,3] = Pi;
   printf("%f",MXGet(0,0));   Beispiel: Pi wieder ausgeben
   B =MXDestroyMatrix();      Speicher wieder freigeben, B =NULL setzen

**************************************************************************/

#ifndef matrix_INC

#define matrix_INC
/* Schutz gegen mehrfaches Einfuegen */

/* Andere oeffentlich benutzte Module */

//#include<iostream>
/* Deklarationen */

/* allgemeine Funktionen */
extern void* MXSetupMatrix(long, long, long);
/* Erstellt eine Matrix und initialisiert alle noetigen Pointer */

extern void MXSetMatrixPointer(void* matrix_pointer);
/* Setzt den internen Matrix-Wurzelzeiger auf die angegebene Matrix-Wurzel.
   Alle folgenden Aufrufe der MX-Funktionen beziehen sich dann auf die
   gewaehlte Matrix. Dadurch entfaellt die Parameteruebergabe des
   Matrix-Wurzelzeigers, das Verfahren ermoeglicht aber trotzdem das
   Verwalten verschiedener Matrizen (Gleichungssysteme) gleichzeitig.
   Setzt auch die MX-Funktionszeiger auf die entsprechenden Funktionen des
   angegebenen Speichermodells! */

extern void* MXGetMatrixPointer(void);
/* Liefert einen Zeiger auf die gerade aktive Matrix */

extern long MXGetDim(void);
/* liefert Dimension der zuvor mit MXSetMatrixPointer angegebenen
   Matrix-Struktur. */

extern void MXDumpGLS(const char*, int, double*, double*);
/* Schreibt das durch die zuvor mit MXSetMatrixPointer angegebene
   Matrix und den Vektor gegebene GLeichungssystem eine Datei "GLS" */

extern void MXEstimateStartVector(double*, double*);
/* Versucht mit Hilfe eines Zeilen-Masslumpings eine gute
   Startloseung fuer den Gleichungsloeser zu finden */

extern void MXRandbed(long, double, double*);
/* Einarbeiten einer Randbedingung in Matrix und Rechte Seite */
extern void MXEliminateIrrNode(long index, int anz_nachbarn, long* nachbarn_index, double* rechts);
/* Elimieren irr. Knoten aus dem Gesamtsystem */
extern void MXResiduum(double* x, double* b, double* ergebnis);
/* Berechnet das Residuum "r = b - A x"  mit der zuvor mit
   MXSetMatrixPointer angegebenen Matrix-Struktur. "ergebnis" wird
   komplett ueberschrieben; der Speicher muss bereits allokiert sein. */

/* MX-Funktionszeiger (Typen) */
typedef void* (*MXPCreateMatrix)(long param1, long param2, long param3);
/* Erzeugen einer neuen Matrix-Struktur im zuvor mit MXSetFunctionPointers
   implizit angegebenen Speichermodell. Das Ergebnis ist der
   Matrix-Wurzelzeiger. Die Parameter werden abhaengig vom Speichermodell
   unterschiedlich benutzt (z.B. fuer Bandbreite etc.). */

typedef void* (*MXPDestroyMatrix)(void);
/* Freigeben der zuvor mit MXSetMatrixPointer angegebenen Matrix-Struktur;
   Ergebnis ist immer NULL. */
typedef void (*MXPResizeMatrix)(long dimension);
/* Veraendern der Groesse der zuvor mit MXSetMatrixPointer angegebenen
   Matrix-Struktur. Der Wurzelzeiger bleibt erhalten, die bisher
   eingetragenen Werte nicht. Die Matrix wird nur vergroessert, nie
   verkleinert. */
typedef void (*MXPInitMatrix)(void);
/* Initialisieren (A[i,j] = 0.0) aller Elemente der zuvor mit
   MXSetMatrixPointer angegebenen Matrix-Struktur. */
typedef int (*MXPSet)(long i, long j, double aij);
/* Setzen des Wertes aij als Matrixwert A[i,j] der zuvor mit
   MXSetMatrixPointer angegebenen Matrix-Struktur. */
typedef int (*MXPInc)(long i, long j, double aij_inc);
/* Inkrementieren des Matrixwertes A[i,j] der zuvor mit
   MXSetMatrixPointer angegebenen Matrix-Struktur um den Wert aij_inc. */
typedef int (*MXPDec)(long i, long j, double aij_dec);
/* Dekrementieren des Matrixwertes A[i,j] der zuvor mit
   MXSetMatrixPointer angegebenen Matrix-Struktur um den Wert aij_dec. */
typedef int (*MXPMul)(long i, long j, double aij_mul);
/* Multiplizieren des Matrixwertes A[i,j] der zuvor mit
   MXSetMatrixPointer angegebenen Matrix-Struktur mit dem Wert aij_mul. */
typedef int (*MXPDiv)(long i, long j, double aij_div);
/* Dividieren des Matrixwertes A[i,j] der zuvor mit
   MXSetMatrixPointer angegebenen Matrix-Struktur durch den Wert aij_div. */
typedef double (*MXPGet)(long i, long j);
/* Ergebnis ist der Matrixwert A[i,j] der zuvor mit MXSetMatrixPointer
   angegebenen Matrix-Struktur. */
typedef void (*MXPTrans)(long i, long j, long ii, long jj);
/* Fuehrt in der zuvor mit MXSetMatrixPointer angegebenen Matrix-Struktur
   folgende Operation durch: A[i,j] = A[ii,jj] */
typedef void (*MXPMatVek)(double* vektor, double* ergebnis);
/* Ausf hren des Matrix-Vektor-Produktes "A * vektor = ergebnis" mit der
   zuvor mit MXSetMatrixPointer angegebenen Matrix-Struktur. "ergebnis"
   wird komplett ueberschrieben; der Speicher muss bereits allokiert sein. */
typedef void (*MXPMatTVek)(double* vektor, double* ergebnis);
/* Ausfuehren des Matrix-Vektor-Produktes "A^T * vektor = ergebnis" mit der
   zuvor mit MXSetMatrixPointer angegebenen Matrix-Struktur. "ergebnis"
   wird komplett ueberschrieben; der Speicher muss bereits allokiert sein. */
typedef void (*MXPVorkond)(int aufgabe, double* x, double* b);
/* Vorkonditionierer nach dem mit vorkond gewaehlten Verfahren */
typedef int (*MXPCopyToAMG1R5Structure)(double* A, int* IA, int* JA, int NDA, int NDIA, int NDJA, double*, double*,
                                        double*, double*);
/* Umkopieren einer Matrix auf AMG1R5-Speicherstruktur */

/* zugehoerige MX-Funktionszeiger */
extern MXPCreateMatrix MXCreateMatrix;
extern MXPDestroyMatrix MXDestroyMatrix;
extern MXPResizeMatrix MXResizeMatrix;
extern MXPInitMatrix MXInitMatrix;
extern MXPSet MXSet;
extern MXPInc MXInc;
extern MXPDec MXDec;
extern MXPMul MXMul;
extern MXPDiv MXDiv;
extern MXPGet MXGet;
extern MXPTrans MXTrans;
extern MXPMatVek MXMatVek;
extern MXPMatTVek MXMatTVek;
extern MXPVorkond MXVorkond;
extern MXPCopyToAMG1R5Structure MXCopyToAMG1R5Structure;

#endif
