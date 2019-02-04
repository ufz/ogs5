/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   ROCKFLOW - Modul: matrix.c

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

   Benutzungsbeispiel fuer 2 Matrizen:

   long dummy =0;             Dummy-Variable mit beliebigem Wert
   void *A = NULL;            Zeiger A fuer Matrix-Wurzel definieren
   void *B = NULL;            Zeiger B fuer Matrix-Wurzel definieren
   A =MXSetupMatrix (10l, 1, dummy);  10x10-Matrix erzeugen, Modell 1
   B =MXSetupMatrix (100l, 2, dummy);  100x100-Matrix erzeugen, Modell 2

   MXSetMatrixPointer(A);     Aktuelle Matrix ist A
   MXInitMatrix();            Matrix A  mit Nullen fuellen
   MXSet(0,0,2*3.1415926);    Beispiel: A[0,0] = 2*Pi;
   MXSetMatrixPointer(B);     Aktuelle Matrix ist B
   MXInitMatrix();            Matrix B  mit Nullen fuellen
   MXSet(2,3,3.1415926);      Beispiel: B[2,3] = Pi;
   printf("%f",MXGet(0,0));   Beispiel: Pi wieder ausgeben
   MXSetMatrixPointer(A);     Aktuelle Matrix ist A
   printf("%f",MXGet(0,0));   Beispiel: 2*Pi wieder ausgeben
   A =MXDestroyMatrix();      Speicher wieder freigeben, A =NULL setzen
   B =MXDestroyMatrix();      Speicher wieder freigeben, B =NULL setzen

   Programmaenderungen:
   11/1995     MSR      Erste Version
   04/1999     AH       Das zuletzt bearbeitetes Speichermodell wird
   in der Variable matrix_type gemerkt.
   Damit werden die Funktionszeiger nicht mehr
   veraendert, wenn das gleiche Speichermodell
   angewaehlt wird. Sinnvoll wenn verschiedenen
   Matrizen in staendigem Wechsel bearbeitet wedren.

   11/1999     CT       Neuer Vorkonditionierer, hilfreich bei Systemen
   mit grossem konstantem Anteil in der Loesung

   2/2000      Ra       Speichermodelle 3 und 4
   3/2000      Ra       Alle Modelle: Besser lesbare Fassung mit Makros
   3/2000      Ra       Modelle 1/2: NumDif zugefuegt: Differenz von
   Zeilen/Spaltennummern. Verkuerzt Abfragen beim
   Einarbeiten von Randbedingungen
   3/2000      Ra       Routine MX_Randbed zugefuegt: Einarbeiten
   von Randbedingungen in Matrix und Rechte Seite,
   alle Modelle
   4/2000      CT       Neu: MXGetMatrixPointer
   7/2000      CT       Fehler in M34Randbed beseitigt
   Modell34 optimiert
   Diagonalenvorkonditionierung fuer Modell34
   Neue Zugriffmethode mit impliziter Speicherung des Matrixtyps
   1/2001  C.Thorenz    Diagonale bleibt bei Randbedingungseinbau und
   Irr.Knotenbehandlung auf Wert.
   3/2000  C.Thorenz    CBLAS fuer M2MatVek eingehaengt
   NumDif korrigiert
   Speichermodell 2 vereinfacht
   7/2002  C.Thorenz    Speichermodell 2 merkt sich max/min Spalteneintrag
   -> schnelleres Suchen
**************************************************************************/

#include "makros.h"
#include "solver.h"
#include <cfloat>

#define noTESTMATRIX_PERF
#define noDUMP

/* Header / Andere intern benutzte Module */
#include "memory.h"
#include "display.h"
#include "files0.h"
#include "mathlib.h"
#include "matrix_routines.h"
// MSHLib
#include "msh_mesh.h"
#include "msh_node.h"
extern MeshLib::CFEMesh* FEMGet(const std::string&);

using MeshLib::CNode;
using namespace Display;

#ifdef _OPENMP
#include <omp.h>
#endif
/* Interne (statische) Deklarationen */
static void* wurzel = NULL; /* interner Matrix-Wurzelzeiger */

/* Allgemeine Funktionen ohne Modellspezifisches */
int MXSetFunctionPointers(int matrix_type); /* Setzt Funktionszeiger */
int MXGetMatrixType(void);                  /* Liefert Matrixtyp */
int MX_Dec(long i, long j, double aij_dec);
int MX_Div(long i, long j, double aij_div);
void MX_Trans(long i, long j, long ii, long jj);

/* Struktur und Funktionen zu Speichermodell 1 (vollbesetzte Matrix) */
typedef struct
{
    long info[2]; /* Muss immer zuerst kommen, enthaelt Dimension und Matrixtyp
                   */
    long max_size;
    long NumDif; /* max. Differenz |Spalten-Zeilennummer| */
    double* matrix;
}

Modell1;

void* M1CreateMatrix(long param1, long param2, long param3);
void* M1DestroyMatrix(void);
void M1ResizeMatrix(long dimension);
void M1InitMatrix(void);
int M1Set(long i, long j, double aij);
int M1Inc(long i, long j, double aij_inc);
int M1Mul(long i, long j, double aij_mul);
double M1Get(long i, long j);
void M1MatVek(double* vektor, double* ergebnis);
void M1MatTVek(double* vektor, double* ergebnis);
void M1Vorkond(int aufgabe, double* x, double* b);
int M1CopyToAMG1R5Structure(double* A, int* IA, int* JA, int NDA, int NDIA,
                            int NDJA, double*, double*, double*, double*);

/* Strukturen und Funktionen zu Speichermodell 2 (i, j, Wert) */
typedef struct
{
    int max_anz;  /* Anzahl der allokierten Spalten (ohne Diagonale) */
    int anz;      /* Anzahl der belegten Spalten (ohne Diagonale) */
    int min_col;  /* Kleinste Spaltennummer */
    int max_col;  /* Groesste Spaltennummer */
    long* index;  /* Index fuer Spalteneintraege */
    double* wert; /* Spalteneintraege */
} M2_Zeile;

typedef struct
{
    long info[2]; /* Muss immer zuerst kommen, enthaelt Dimension und Matrixtyp
                   */
    long max_size;
    long NumDif;     /* max. Differenz |Spalten-Zeilennummer| */
    M2_Zeile* zeile; /* Eintraege in den Zeilen; Diagonale in zeile[0] */
} Modell2;

typedef struct
{
    long NumDif;  // dummy
}

Modell5;

void* M2CreateMatrix(long param1, long param2, long param3);
void* M2DestroyMatrix(void);
void M2ResizeMatrix(long dimension);
void M2InitMatrix(void);
int M2Set(long i, long j, double aij);
int M2Inc(long i, long j, double aij_inc);
int M2Mul(long i, long j, double aij_mul);
double M2Get(long i, long j);
void M2MatVek(double* vektor, double* ergebnis);
/*-----------------------------------------------------------------------
 * JAD format Modell 5
 * */
int M5Inc(long i, long j, double aij_inc);
double M5Get(long i, long j);
int M5Set(long i, long j, double e_val);
void M5MatVek(double* b, double* erg);
void M5InitMatrix(void);
void* M5DestroyMatrix(void);
void M5Vorkond(int aufgabe, double* x, double* b);

// void M5CreateMatrix(void);
void* M5CreateMatrix(long param1, long param2, long param3);

void insertionSort1_des(int numbers[], int numbers1[], int array_size);
void transM2toM5(void);
int *jd_ptr1, *jd_ptr2, *jd_ptr, *jd_ptrf, *col_ind;
long* diag5_i;
double *temp_ergebnis, *jdiag;
int m_count, jd_ptr_max, Dim_L;
/*-----------------------------------------------------------------------
 * ITPACKV format Modell 6
 * */
void M6MatVek(double* b, double* erg);
void transM2toM6(void);
double* itpackv;
int* it_col;
int itpack_max;
/*----------------------------------------------------------------------*/
void M2MatTVek(double* vektor, double* ergebnis);
void M2Vorkond(int aufgabe, double* x, double* b);
int M2CopyToAMG1R5Structure(double* A, int* IA, int* JA, int NDA, int NDIA,
                            int NDJA, double*, double*, double*, double*);

/* Strukturen und Funktionen zu Speichermodell 3/4 (i, j, Wert) */
typedef struct
{ /* einzelner Eintrag in der Spalte k mit Index i */
    long Index;
    double Aik[4];
} M34_aik;

typedef struct
{                /* Spalte k der Obermatrix (Zeile k der Untermatrix) */
    int max_anz; /* Anzahl der allokierten Elemente i (ohne Diagonale) */
    int anz;     /* Anzahl eingespeicherter Aik (ohne Diagonale) */
    long rechts; /* rechts[k]: letzte Spalte j, die ein Akj enthaelt */
    M34_aik* Ak; /* Spalten/Zeilenelemente */
} M34_Spalte;

typedef struct
{
    long info[2]; /* Muss immer zuerst kommen, enthaelt Dimension und Matrixtyp
                   */
    long max_size;
    int usym; /* 0:   symmetrisch, Modell 3 */
    /* 1: unsymmetrisch, Modell 4 */
    int stat; /* -1: nicht initialisiert */
    /*  0: initialisiert */
    /*  1: aik eingespeichert oder veraendert, */
    /*  2: ILU-Zerlegung gelaufen */
    int l34_Aik; /* L�nge von Aik */
    /* 2: symmetrisch, 4: unsymmetrisch */
    int i34_Bik; /* Start Bik (Preconditioner) innerhalb aik */
    /* 1: symmetrisch, 2: unsymmetrisch */
    double* Diag;       /* Akk als Vektor */
    double* PreD;       /* Bkk (ILU-Preconditioner) */
    M34_Spalte* Spalte; /* Obermatrixspalten (Untermatrixzeilen) */
} Modell34;

/* Die nachfolgenden Makros dienen der besseren Lesbarkeit (Kuerze!)
   des Codes. Sie werden am Ende dieser Quelle saemtlich undefiniert!
   Allgemeines: */
#define dim ((long*)wurzel)[0] /* Dimension der aktuellen Matrix */
#define matrix_type                                           \
    ((long*)wurzel)[1] /* Speichermodell der aktuellen Matrix \
                        */
#define Maxim(x, v) \
    if (x < (v))    \
    x = v /* x =max(x,v) */
#define Minim(x, v) \
    if (x > (v))    \
    x = v
/* x =min(x,v)
                                          Modell 1: */
#define Matrix1 w->matrix
#define Aik1(i, k) Matrix1[i * dim + k]
/* Modell 2: */
#define Zeil2(i) w->zeile[i]
#define Aik2(i, j) Zeil2(i).wert[j]
#define Ind2(i, j) Zeil2(i).index[j]
#define Diag2(i) Zeil2(i).wert[0]
/* Modelle 3,4: */
#define Sp34(k) w->Spalte[k]
#define Aik34(k, j, teil) Sp34(k).Ak[j].Aik[teil]
#define Bik34(k, j, teil) Sp34(k).Ak[j].Aik[(w->i34_Bik) + teil]
#define Ind34(k, j) Sp34(k).Ak[j].Index

void* M34CreateMatrix(long param1, long param2, long param3);
void* M34DestroyMatrix(void);
void M34ResizeMatrix(long dimension);
void M34InitMatrix(void);
int M34Set(long i, long j, double aij);
int M34Inc(long i, long j, double aij_inc);
int M34Mul(long i, long j, double aij_mul);
double M34Get(long i, long j);
void M34MatVek(double* vektor, double* ergebnis);
void M34MatTVek(double* vektor, double* ergebnis);
void M34Vorkond(int aufgabe, double* x, double* b);
int M34CopyToAMG1R5Structure(double* A, int* IA, int* JA, int NDA, int NDIA,
                             int NDJA, double*, double*, double*, double*);

void MX_Exit(const char* caller, int errcode)
{
    char text[1024];

    strcpy(text, caller);
    strcat(text, ": ");

    switch (errcode)
    {
        case 0:
            strcat(text, "Negative Dimension");
            break;
        case 1:
            strcat(text, "Keine Wurzel eingetragen");
            break;
        case 2:
            strcat(text, "Feldgrenzenueberschreitung");
            break;
        case 3:
            strcat(text, "Nullzeiger");
            break;
        case 4:
            strcat(text, "Unbekanntes Speichermodell");
            break;
        case 5:
            strcat(text, "Argumente falsch");
            break;
    }
    DisplayErrorMsg(strcat(text, " -> Abbruch!"));
    exit(1);
}

/* Definitionen */
MXPCreateMatrix MXCreateMatrix;
MXPDestroyMatrix MXDestroyMatrix;
MXPResizeMatrix MXResizeMatrix;
MXPInitMatrix MXInitMatrix;
MXPSet MXSet;
MXPInc MXInc;
MXPDec MXDec;
MXPMul MXMul;
MXPDiv MXDiv;
MXPGet MXGet;
MXPTrans MXTrans;
MXPMatVek MXMatVek;
MXPMatTVek MXMatTVek;
MXPVorkond MXVorkond;
MXPCopyToAMG1R5Structure MXCopyToAMG1R5Structure;

/* allgemeine Funktionen */

/*************************************************************************
   ROCKFLOW - Funktion: MXSetupMatrix

   Aufgabe:
   Setzt Funktionszeiger, erzeugt neue Matrix-Struktur, initialisiert
   die Matrix

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E long param1: Dimension der Matrix
   E long param2: Speichertechnik der Matrix
   E long param3: nicht benutzt

   Ergebnis:
   Matrix-Wurzelzeiger

   Programmaenderungen:
   7/2000   CT   Erste Version

*************************************************************************/
void* MXSetupMatrix(long param1, long param2, long param3)
{
    MXSetFunctionPointers((int)param2);
    wurzel = MXCreateMatrix(param1, param2, param3);
    MXSetMatrixPointer(wurzel);
    MXInitMatrix();
    return wurzel;
}

/*************************************************************************
   ROCKFLOW - Funktion: MXSetFunctionPointers

   Aufgabe:
   Setzt die MX-Funktionszeiger auf die entsprechenden Funktionen des
   angegebenen Speichermodells

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E int type: Speichermodell

   Ergebnis:
   0 bei unzulaessigem Modelltyp, sonst 1

   Programmaenderungen:
   11/1995     MSR        Erste Version

*************************************************************************/
int MXSetFunctionPointers(int type)
{
    switch (type)
    {
        case 1:
            MXCreateMatrix = M1CreateMatrix;
            MXDestroyMatrix = M1DestroyMatrix;
            MXResizeMatrix = M1ResizeMatrix;
            MXInitMatrix = M1InitMatrix;
            MXSet = M1Set;
            MXInc = M1Inc;
            MXMul = M1Mul;
            MXGet = M1Get;
            MXMatVek = M1MatVek;
            MXMatTVek = M1MatTVek;
            MXVorkond = M1Vorkond;
            MXCopyToAMG1R5Structure = M1CopyToAMG1R5Structure;
            break;
        case 2:
            MXCreateMatrix = M2CreateMatrix;
            MXDestroyMatrix = M2DestroyMatrix;
            MXResizeMatrix = M2ResizeMatrix;
            MXInitMatrix = M2InitMatrix;
            MXSet = M2Set;
            MXInc = M2Inc;
            MXMul = M2Mul;
            MXGet = M2Get;
            MXMatVek = M2MatVek;
            MXMatTVek = M2MatTVek;
            MXVorkond = M2Vorkond;
            MXCopyToAMG1R5Structure = M2CopyToAMG1R5Structure;
            break;
        case 3:
            MXCreateMatrix = M34CreateMatrix;
            MXDestroyMatrix = M34DestroyMatrix;
            MXResizeMatrix = M34ResizeMatrix;
            MXInitMatrix = M34InitMatrix;
            MXSet = M34Set;
            MXInc = M34Inc;
            MXMul = M34Mul;
            MXGet = M34Get;
            MXMatVek = M34MatVek;
            MXMatTVek = M34MatTVek;
            MXVorkond = M34Vorkond;
            MXCopyToAMG1R5Structure = M34CopyToAMG1R5Structure;
            break;
        case 4:
            MXCreateMatrix = M34CreateMatrix;
            MXDestroyMatrix = M34DestroyMatrix;
            MXResizeMatrix = M34ResizeMatrix;
            MXInitMatrix = M34InitMatrix;
            MXSet = M34Set;
            MXInc = M34Inc;
            MXMul = M34Mul;
            MXGet = M34Get;
            MXMatVek = M34MatVek;
            MXMatTVek = M34MatTVek;
            MXVorkond = M34Vorkond;
            MXCopyToAMG1R5Structure = M34CopyToAMG1R5Structure;
            break;
        case 5:
            MXCreateMatrix = M5CreateMatrix;
            MXDestroyMatrix = M5DestroyMatrix;
            MXResizeMatrix = M2ResizeMatrix;
            MXInitMatrix = M5InitMatrix;
            MXSet = M5Set;
            MXInc = M5Inc;  // M2Inc;
            MXMul = M2Mul;
            MXGet = M5Get;
            MXMatVek = M5MatVek;
            MXMatTVek = M2MatTVek;
            MXVorkond = M5Vorkond;
            MXCopyToAMG1R5Structure = M2CopyToAMG1R5Structure;
            break;
        case 6:
            MXCreateMatrix = M2CreateMatrix;
            MXDestroyMatrix = M2DestroyMatrix;
            MXResizeMatrix = M2ResizeMatrix;
            MXInitMatrix = M2InitMatrix;
            MXSet = M2Set;
            MXInc = M2Inc;
            MXMul = M2Mul;
            MXGet = M2Get;
            MXMatVek = M6MatVek;
            MXMatTVek = M2MatTVek;
            MXVorkond = M2Vorkond;
            MXCopyToAMG1R5Structure = M2CopyToAMG1R5Structure;
            break;
        default:
            return 0;
    }

    MXDec = MX_Dec; /* nicht modellspezifisch */
    MXDiv = MX_Div;
    MXTrans = MX_Trans;
    return 1;
}

/*************************************************************************
   ROCKFLOW - Funktion: MXSetMatrixPointer

   Aufgabe:
   Setzt den internen Matrix-Wurzelzeiger auf die angegebene Matrix-Wurzel.
   Alle folgenden Aufrufe der MX-Funktionen beziehen sich dann auf die
   gewaehlte Matrix. Dadurch entfaellt die Parameteruebergabe des
   Matrix-Wurzelzeigers, das Verfahren ermoeglicht aber trotzdem das
   Verwalten verschiedener Matrizen (Gleichungssysteme) gleichzeitig.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E void *matrix_pointer: Matrix-Wurzelzeiger

   Ergebnis:
   - void -

   Programmaenderungen:
   11/1995     MSR        Erste Version
   7/2000     CT         Speichermodell implizit in Matrix-Datenstruktur

*************************************************************************/
void MXSetMatrixPointer(void* matrix_pointer)
{
    wurzel = matrix_pointer;
    MXSetFunctionPointers(matrix_type);
}

/*************************************************************************
   ROCKFLOW - Funktion: MXGetMatrixPointer

   Aufgabe:
   Liefert den internen Matrix-Wurzelzeiger fuer aktuelle Matrix

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)

   Ergebnis:
   void *matrix_pointer: Matrix-Wurzelzeiger-

   Programmaenderungen:
   3/2000     C.Thorenz       Erste Version

*************************************************************************/
void* MXGetMatrixPointer(void)
{
    return wurzel;
}

/*************************************************************************
   ROCKFLOW - Funktion: MXGetDim

   Aufgabe:
   Liefert die Dimension der zuvor mit MXSetMatrixPointer angegebenen
   Matrix-Struktur.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   - void -

   Ergebnis:
   Dimension der Matrix

   Programmaenderungen:
   11/1995     MSR        Erste Version

*************************************************************************/
long MXGetDim(void)
{
    return dim;
}

/*************************************************************************
   ROCKFLOW - Funktion: MXGetMatrixType

   Aufgabe:
   Liefert das Speichermodell der aktuellen Matrix-Struktur.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   - void -

   Ergebnis:
   Speichermodell der Matrix

   Programmaenderungen:
   04/1999     AH        Erste Version
   7/2000     CT        Speichermodell implizit in Matrix-Datenstruktur

*************************************************************************/
int MXGetMatrixType(void)
{
    return matrix_type;
}

/*************************************************************************
   ROCKFLOW - Funktion: M#CreateMatrix

   Aufgabe:
   Erzeugen einer neuen Matrix-Struktur

   speziell Modelle 2, 3, 4:
   Die Struktur muss vor Benutzung unbedingt initialisiert werden!

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E long param1: Dimension der Matrix
   E long param2: Speichertechnik der Matrix
   E long param3: nicht benutzt

   Ergebnis:
   Matrix-Wurzelzeiger

   Programmaenderungen:
   Modell 1:    11/1995     MSR        Erste Version
   Modell 2:    11/1995     MSR        Erste Version
   Modell 3, 4:  2/2000     Ra         Letzte Version
   7/2000     CT        Speichermodell implizit in Matrix-Datenstruktur

*** Modell 1 ************************************************************/
void* M1CreateMatrix(long param1, long param2, long /*param3*/)
{
    static Modell1* w;

#ifdef ERROR_CONTROL
    if (param1 < 0)
        MX_Exit("M1CreateMatrix", 0);
#endif

    w = (Modell1*)Malloc(sizeof(Modell1));
    MXSetMatrixPointer((void*)w);
    Matrix1 = (double*)Malloc(param1 * param1 * sizeof(double));
    dim = w->max_size = param1;
    matrix_type = param2;
    w->NumDif = 0;
    return (void*)w;
}

/**** Modell 2 ************************************************************/
void* M2CreateMatrix(long param1, long param2, long /*param3*/)
{
    static Modell2* w;
    static long i;

#ifdef ERROR_CONTROL
    if (param1 < 0)
        MX_Exit("M2CreateMatrix", 0);
#endif

    w = (Modell2*)Malloc(sizeof(Modell2));
    MXSetMatrixPointer((void*)w);
    w->zeile = (M2_Zeile*)Malloc(param1 * sizeof(M2_Zeile));
    for (i = 0; i < param1; i++)
    {
        Zeil2(i).index = (long*)Malloc(sp2_start * sizeof(long));
        Zeil2(i).wert = (double*)Malloc(sp2_start * sizeof(double));
        Zeil2(i).max_anz = sp2_start;
        Zeil2(i).anz = 0;
        Zeil2(i).min_col = i;
        Zeil2(i).max_col = i;
    }
    dim = w->max_size = param1;
    matrix_type = param2;
    w->NumDif = 0;
    return (void*)w;
}

/**** Modell 3 und 4 ******************************************************/
void* M34CreateMatrix(long param1, long param2, long /*param3*/)
{
    Modell34* w = NULL;
    register long i;

#ifdef ERROR_CONTROL
    if (param1 < 0)
        MX_Exit("M34CreateMatrix", 0);
#endif

    w = (Modell34*)Malloc(sizeof(Modell34));
    MXSetMatrixPointer((void*)w);
    dim = w->max_size = param1;
    matrix_type = param2;
    w->usym = matrix_type - 3;
    w->stat = -1;
    w->l34_Aik = (w->usym + 1) * 2;
    w->i34_Bik = w->usym + 1;
    w->Diag = (double*)Malloc(param1 * sizeof(double));
    /* Diag. Precond. */
    w->PreD = (double*)Malloc(param1 * sizeof(double));
    w->Spalte = (M34_Spalte*)Malloc(param1 * sizeof(M34_Spalte));
    for (i = 0; i < param1; i++)
    { /* Start- und Inkrementlaengen wie Modell 2 */
        Sp34(i).max_anz = sp2_start;
        Sp34(i).anz = 0;
        Sp34(i).rechts = 0;
        Sp34(i).Ak = (M34_aik*)Malloc(sp2_start * sizeof(M34_aik));
    }
    return (void*)w; /* Die ILU-Pointer sind auch schon initialisiert */
}

/*************************************************************************
   ROCKFLOW - Funktion: MxDestroyMatrix

   Aufgabe:
   Freigeben der zuvor mit MXSetMatrixPointer angegebenen Matrix-Struktur.

   Formalparameter: keine

   Ergebnis:
   NULL

   Programmaenderungen:
   Modell 1:    11/1995     MSR        Erste Version
   Modell 2:    11/1995     MSR        Erste Version
   Modell 3, 4:  2/2000     Ra         Letzte Version
*** Modell 1 ************************************************************/
void* M1DestroyMatrix(void)
{
    Modell1* w = (Modell1*)wurzel;
    if (wurzel == NULL)
        return NULL;
    Matrix1 = (double*)Free(Matrix1);
    wurzel = Free(w);
    return NULL;
}

/**** Modell 2 ************************************************************/
void* M2DestroyMatrix(void)
{
    Modell2* w = (Modell2*)wurzel;
    static long i;
    if (wurzel == NULL)
        return NULL;

#ifdef TESTMATRIX_PERF
    /* statistische Auswertung der Speicherstruktur */
    {
        long j, anz_inc = 0, max_inc = 0, zaehler = 0;
        for (j = 0; j < w->max_size; j++)
        {
            Maxim(max_inc, Zeil2(j).max_anz);
            if (Zeil2(j).max_anz > sp2_start)
            {
                anz_inc++;
                zaehler += ((Zeil2(j).max_anz - sp2_start) / sp2_inc);
            }
        }
        DisplayMsgLn("Statistische Auswertung der Speicherstruktur:");
        DisplayMsg(" - max. Dimension des Gleichungssystems: ");
        DisplayLong(w->max_size);
        DisplayMsgLn("");
        DisplayMsg(" - Ausgangsgroesse der Zeileneintraege sp2_start: ");
        DisplayLong((long)sp2_start);
        DisplayMsgLn("");
        DisplayMsg(" - Erhoehung der Zeileneintraege sp2_inc: ");
        DisplayLong((long)sp2_inc);
        DisplayMsgLn("");
        DisplayMsg(" - Anzahl der erhoehten Zeileneintraege: ");
        DisplayLong(anz_inc);
        DisplayMsgLn("");
        DisplayMsg(" - Groesste Groesse eines Zeileneintrages: ");
        DisplayLong(max_inc);
        DisplayMsgLn("");
        DisplayMsg(" - Gesamtzahl der Erhoehungen aller Zeileneintraege: ");
        DisplayLong(zaehler);
        DisplayMsgLn("");
    }
#endif

    for (i = 0; i < w->max_size; i++)
    {
        Zeil2(i).index = (long*)Free(Zeil2(i).index);
        Zeil2(i).wert = (double*)Free(Zeil2(i).wert);
    }
    w->zeile = (M2_Zeile*)Free(w->zeile);
    wurzel = Free(w);
    return NULL;
}

/**** Modell 3, 4 *********************************************************/
void* M34DestroyMatrix(void)
{
    register long i;
    Modell34* w = (Modell34*)wurzel;
    if (w == NULL)
        return NULL;

#ifdef TESTMATRIX_PERF
    { /* statistische Auswertung der Speicherstruktur */
        long anz_inc = 0, max_inc = 0, zaehler = 0;
        for (j = 0; j < w->max_size; j++)
        {
            Maxim(max_inc, Sp34(j).max_anz);
            if (Sp34(j).max_anz > sp2_start)
            {
                anz_inc++;
                zaehler += (Sp34(j).max_anz - sp2_start) / sp2_inc;
            }
        }
        DisplayMsgLn("Statistische Auswertung der Speicherstruktur:");
        DisplayMsg(" - max. Dimension des Gleichungssystems: ");
        DisplayLong(w->max_size);
        DisplayMsgLn("");
        DisplayMsg(" - Ausgangsgroesse der Zeileneintraege sp2_start: ");
        DisplayLong((long)sp2_start);
        DisplayMsgLn("");
        DisplayMsg(" - Erhoehung der Zeileneintraege sp2_inc: ");
        DisplayLong((long)sp2_inc);
        DisplayMsgLn("");
        DisplayMsg(" - Anzahl der erhoehten Zeileneintraege: ");
        DisplayLong(anz_inc);
        DisplayMsgLn("");
        DisplayMsg(" - Groesste Groesse eines Zeileneintrages: ");
        DisplayLong(max_inc);
        DisplayMsgLn("");
        DisplayMsg(" - Gesamtzahl der Erhoehungen aller Zeileneintraege: ");
        DisplayLong(zaehler);
        DisplayMsgLn("");
    }
#endif

    for (i = 0; i < w->max_size; i++)
        Sp34(i).Ak = (M34_aik*)Free(Sp34(i).Ak);

    w->Spalte = (M34_Spalte*)Free(w->Spalte);
    w->Diag = (double*)Free(w->Diag);
    w->PreD = (double*)Free(w->PreD);
    wurzel = Free(w);
    return NULL;
}

/**** Modell 5 ************************************************************/
// WW/PA 08/02/2006
void* M5DestroyMatrix(void)
{
    free(jd_ptr1);
    free(jd_ptr2);
    free(jd_ptr);
    free(jd_ptrf);
    free(temp_ergebnis);
    free(col_ind);
    free(jdiag);
    free(diag5_i);
    return (void*)1;
}

/*************************************************************************
   ROCKFLOW - Funktion: M#ResizeMatrix

   Aufgabe:
   Veraendern der Groesse der zuvor mit MXSetMatrixPointer angegebenen
   Matrix-Struktur. Der Wurzelzeiger bleibt erhalten, die bisher
   eingetragenen Werte nicht. Die Matrix wird nur vergroessert, nie
   verkleinert.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E long dimension: Dimension der Matrix

   Ergebnis:
   - void -

   Programmaenderungen:
   Modell 1:   11/1995     MSR        Erste Version
   Modell 2:   11/1995     MSR        Erste Version
   Modell 3,4   2/2000     Ra         Letzte Version
*** Modell 1 ************************************************************/
void M1ResizeMatrix(long dimension)
{
    Modell1* w = (Modell1*)wurzel;

#ifdef ERROR_CONTROL
    if (w == NULL)
        MX_Exit("M1ResizeMatrix", 1);
    if (dimension < 0)
        MX_Exit("M1ResizeMatrix", 0);
#endif

    dim = dimension;
    if (w->max_size < dimension)
    {
        Matrix1 = (double*)Free(w->matrix);
        Matrix1 = (double*)Malloc(dim * dim * sizeof(double));
        w->max_size = dim;
    }
}

/**** Modell 2 ************************************************************/
void M2ResizeMatrix(long dimension)
{
    Modell2* w = (Modell2*)wurzel;
    static long i;

#ifdef ERROR_CONTROL
    if (w == NULL)
        MX_Exit("M2ResizeMatrix", 1);
    if (dimension < 0)
        MX_Exit("M2ResizeMatrix", 0);
#endif

    dim = dimension;
    if (w->max_size < dim)
    {
        w->zeile = (M2_Zeile*)Realloc(w->zeile, dim * sizeof(M2_Zeile));
        for (i = w->max_size; i < dim; i++)
        {
            Zeil2(i).index = (long*)Malloc(sp2_start * sizeof(long));
            Zeil2(i).wert = (double*)Malloc(sp2_start * sizeof(double));
            Zeil2(i).max_anz = sp2_start;
            Zeil2(i).anz = 0;
            Zeil2(i).min_col = i;
            Zeil2(i).max_col = i;
        }
        w->max_size = dim;
    }
}

/**** Modell 3, 4 *********************************************************/
void M34ResizeMatrix(long dimension)
{
    register long i;
    Modell34* w = (Modell34*)wurzel;

#ifdef ERROR_CONTROL
    if (w == NULL)
        MX_Exit("M34ResizeMatrix", 1);
    if (dimension < 0)
        MX_Exit("M34ResizeMatrix", 0);
#endif

    dim = dimension;
    w->stat = -1;
    if (w->max_size < dim)
    {
        w->Spalte = (M34_Spalte*)Realloc(w->Spalte, dim * sizeof(M34_Spalte));
        w->Diag = (double*)Realloc(w->Diag, dim * sizeof(double));
        w->PreD = (double*)Realloc(w->PreD, dim * sizeof(double));
        for (i = w->max_size; i < dim; i++)
        { /* Laengen wie Modell 2 */
            Sp34(i).max_anz = sp2_start;
            Sp34(i).anz = 0;
            Sp34(i).rechts = 0;
            Sp34(i).Ak = (M34_aik*)Malloc(sp2_start * sizeof(M34_aik));
        }
        w->max_size = dim;
    }
}

/*************************************************************************
   ROCKFLOW - Funktion: M#InitMatrix

   Aufgabe:
   Initialisieren (A[i,j] = 0.0) aller Elemente der zuvor mit
   MXSetMatrixPointer angegebenen Matrix-Struktur.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   - void -

   Ergebnis:  - void -

   Programmaenderungen:
   Modell 1:   11/1995     MSR        Erste Version
   Modell 2:   11/1995     MSR        Erste Version
   Modell 3,4:  2/2000     Ra         Letzte Version
*** Modell 1 ************************************************************/
void M1InitMatrix(void)
{
    Modell1* w = (Modell1*)wurzel;
    register long i, j = dim * dim;
    for (i = 0; i < j; i++)
        Matrix1[i] = 0.0;
    w->NumDif = 0;
}

/**** Modell 2 ************************************************************/
void M2InitMatrix(void)
{
    register long i;
    register int j;
    Modell2* w = (Modell2*)wurzel;
    for (i = 0; i < dim; i++)
    {
        Zeil2(i).anz = 1;
        for (j = 0; j < Zeil2(i).max_anz; j++)
            Aik2(i, j) = 0.0;
        Ind2(i, 0) = i; /* Diagonale immer in Zeilenelement 0! */
        Zeil2(i).min_col = i;
        Zeil2(i).max_col = i;
    }
    w->NumDif = 0;
}

/**** Modell 3,4 **********************************************************/
void M34InitMatrix(void)
{
    register long i;
    Modell34* w = (Modell34*)wurzel;
    for (i = 0; i < dim; i++)
    {
        w->Diag[i] = 0.0;   /* Diagonale */
        Sp34(i).anz = 0;    /* geloescht, keine Nebenelemente */
        Sp34(i).rechts = i; /* fuer RDB */
    }
    w->stat = 0;
}

/**** Modell 5 ************************************************************/
// WW/PA 08/02/2006
void M5InitMatrix(void)
{
    long k;
    for (k = 0; k < m_count; k++)
        jdiag[k] = 0.0;
}

/**************************************************************************/

/* ROCKFLOW - Funktion: M#Set

   Aufgabe:
   Setzen des Wertes aij als Matrixwert A[i,j] der zuvor mit
   MXSetMatrixPointer angegebenen Matrix-Struktur.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E long i, j: Indizes der Matrixposition
   E double aij: Wert A[i,j]

   Ergebnis:
   0

   Programmaenderungen:
   Modell 1:   11/1995     MSR        Erste Version
   Modell 2:   11/1995     MSR        Erste Version
   Modell 3,4   2/2000     Ra         Letzte Version
 *** Modell 1 *********************************************************** */
int M1Set(long i, long j, double aij)
{
    Modell1* w = (Modell1*)wurzel;

#ifdef ERROR_CONTROL
    if ((i >= dim) || (j >= dim) || (i < 0) || (j < 0))
        MX_Exit("M1Set", 2);
#endif

    Aik1(i, j) = aij;
    Maxim(w->NumDif, labs(i - j));
    return 0;
}

/**** Modell 2 ************************************************************/
int M2Set(long i, long j, double aij)
{
    Modell2* w = (Modell2*)wurzel;
    register int k = 0;

#ifdef ERROR_CONTROL
    if ((i >= dim) || (j >= dim) || (i < 0) || (j < 0))
        MX_Exit("M2Set", 2);
#endif

    if (j < Zeil2(i).min_col)
        Zeil2(i).min_col = j;
    if (j > Zeil2(i).max_col)
        Zeil2(i).max_col = j;

    /* Alle Spalteneintraege durchsuchen */
    while (k < Zeil2(i).anz)
    {
        if (j == Ind2(i, k))
        {
            Aik2(i, k) = aij;
            return 0;
        } /* gefunden! */
        k++;
    }
    if (fabs(aij) < MKleinsteZahl)
        return 0; /* Abbruch bei Nullwert */

    (Zeil2(i).anz)++;
    /* Neuen Speicher holen */
    if (Zeil2(i).anz > Zeil2(i).max_anz)
    {
        Zeil2(i).max_anz += sp2_inc; /* Tabelle vergroessern */
        Zeil2(i).index =
            (long*)Realloc(Zeil2(i).index, Zeil2(i).max_anz * sizeof(long));
        Zeil2(i).wert =
            (double*)Realloc(Zeil2(i).wert, Zeil2(i).max_anz * sizeof(double));
    }
    Ind2(i, k) = j;
    Aik2(i, k) = aij;
    Maxim(w->NumDif, labs(i - j));
    return 0;
}

/**** Modell 3,4 **********************************************************/
int M34Set(long i1, long k1, double aik)
{
    register int j = 0;
    register long i = i1, k = k1;
    int u = 0, j1;
    Modell34* w = (Modell34*)wurzel;

#ifdef ERROR_CONTROL
    if ((i >= dim) || (k >= dim) || (i < 0) || (k < 0))
        MX_Exit("M34Set", 2);
#endif

    w->stat = 1;
    if (i == k)
    {
        w->Diag[i] = aik;
        return 0;
    } /* Diagonalelement */
    if (i > k)
    { /* Untermatrixelement */
        if (!(w->usym))
            return 0; /* symmetrische Matrix, ignorieren */
        u = 1;
        i = k;
        k = i1;
    }
    /* Alle Spalteneintraege durchsuchen */
    while (j < Sp34(k).anz && i >= Ind34(k, j))
    {
        if (i == Ind34(k, j))
        {
            Aik34(k, j, u) = aik;
            return 0;
        } /* vorhanden */
        j++;
    }

    /* Das Element ist noch nicht vorhanden */
    // WW  if (fabs(aik)< MKleinsteZahl)
    if (fabs(aik) < DBL_MIN)
        return 0;

    if (Sp34(k).anz == Sp34(k).max_anz)
    { /* Spalte verlaengern */
        j1 = Sp34(k).max_anz;
        Sp34(k).max_anz += sp2_inc;
        Sp34(k).Ak =
            (M34_aik*)Realloc(Sp34(k).Ak, Sp34(k).max_anz * sizeof(M34_aik));
    }
    /* Spalteneintraege sind nach aufsteigender Zeilennummer sortiert! */
    for (j1 = Sp34(k).anz; j1 > j; j1--) /* Rest verschieben */
        Sp34(k).Ak[j1] = Sp34(k).Ak[j1 - 1];

    Sp34(k).anz++;
    Ind34(k, j) = i;
    if (w->usym)
        Aik34(k, j, 1 - u) = 0.0;
    Aik34(k, j, u) = aik;
    Maxim(Sp34(i).rechts, k); /* fuer RDB */
    return 0;
}

/*************************************************************************
   ROCKFLOW - Funktion: M#Inc

   Aufgabe:
   Inkrementieren des Matrixwertes A[i,j] der zuvor mit
   MXSetMatrixPointer angegebenen Matrix-Struktur um den Wert aij_inc.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E long i, j: Indizes der Matrixposition
   E double aij_inc: Inkrement von A[i,j]

   Ergebnis:
   0

   Programmaenderungen:
   Modell 1,2   11/1995     MSR        Erste Version
   Modell 3,4    2/2000     Ra         Letzte Version
   3/2002     CT         NumDif korrigiert
*** Modell 1 ************************************************************/
int M1Inc(long i, long j, double aij_inc)
{
    Modell1* w = (Modell1*)wurzel;

#ifdef ERROR_CONTROL
    if ((i >= dim) || (j >= dim) || (i < 0) || (j < 0))
        MX_Exit("M1Inc", 2);
#endif

    Maxim(w->NumDif, labs(i - j));
    Aik1(i, j) += aij_inc;
    return 0;
}

/**** Modell 2 ************************************************************/
int M2Inc(long i, long j, double aij_inc)
{
    Modell2* w = (Modell2*)wurzel;
    register int k = 0;

#ifdef ERROR_CONTROL
    if ((i >= dim) || (j >= dim) || (i < 0) || (j < 0))
        MX_Exit("M2Inc", 2);
#endif

    if (fabs(aij_inc) < MKleinsteZahl)
        return 0; /* Abbruch bei Nullwert */
    while (k < Zeil2(i).anz)
    { /* alle Eintraege durchsuchen */
        if (j == Ind2(i, k))
        {
            Aik2(i, k) += aij_inc;
            return 0;
        } /* gefunden! */
        k++;
    }
    Maxim(w->NumDif, labs(i - j));
    return M2Set(i, j, aij_inc); /* neuer Eintrag */
}

/**** Modell 3, 4 *********************************************************/
int M34Inc(long i1, long k1, double aik)
{
    register int j = 0;
    register long i = i1, k = k1;
    int u = 0;
    Modell34* w = (Modell34*)wurzel;

#ifdef ERROR_CONTROL
    if ((i >= dim) || (k >= dim) || (i < 0) || (k < 0))
        MX_Exit("M34Inc", 2);
#endif

    w->stat = 1;
    if (i == k)
    {
        w->Diag[i] += aik;
        return 0;
    } /* Diagonalelement */
    if (i > k)
    { /* Untermatrixelement */
        if (!(w->usym))
            return 0; /* symmetrische Matrix, ignorieren */
        u = 1;
        i = k;
        k = i1;
    }
    /* Alle Spalteneintraege durchsuchen */
    while (j < Sp34(k).anz && i >= Ind34(k, j))
    {
        if (i == Ind34(k, j))
        {
            Aik34(k, j, u) += aik;
            return 0;
        }
        j++;
    }
    return M34Set(i1, k1, aik); /* war noch nicht vorhanden */
}

/*************************************************************************
   ROCKFLOW - Funktion: MX#Dec

   Aufgabe:
   Dekrementieren des Matrixwertes A[i,j] der zuvor mit
   MXSetMatrixPointer angegebenen Matrix-Struktur um den Wert aij_dec.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E long i, j: Indizes der Matrixposition
   E double aij_dec: Dekrement von A[i,j]

   Ergebnis:
   0

   Programmaenderungen:
   Modell 1, 2:   11/1995     MSR   Erste Version
   alle Modelle :  2/2000     Ra    Neufassung mit MXInc
*** nicht modellspezifisch! *********************************************/
int MX_Dec(long i, long j, double aij_dec)
{
    return MXInc(i, j, -aij_dec);
}

/*************************************************************************
   ROCKFLOW - Funktion: M#Mul

   Aufgabe:
   Multiplizieren des Matrixwertes A[i,j] der zuvor mit
   MXSetMatrixPointer angegebenen Matrix-Struktur mit dem Wert aij_mul.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E long i, j: Indizes der Matrixposition
   E double aij_mul: Faktor von A[i,j]

   Ergebnis:
   0

   Programmaenderungen:
   Modell 1, 2:  11/1995     MSR        Erste Version
   Modell 3, 4:   2/2000     Ra         Letzte Version
*** Modell 1 ************************************************************/
int M1Mul(long i, long j, double aij_mul)
{
    Modell1* w = (Modell1*)wurzel;

#ifdef ERROR_CONTROL
    if ((i >= dim) || (j >= dim) || (i < 0) || (j < 0))
        MX_Exit("M1Mul", 2);
#endif

    Aik1(i, j) *= aij_mul;
    return 0;
}

/**** Modell 2 ************************************************************/
int M2Mul(long i, long j, double aij_mul)
{
    Modell2* w = (Modell2*)wurzel;
    register int k = 0;

#ifdef ERROR_CONTROL
    if ((i >= dim) || (j >= dim) || (i < 0) || (j < 0))
        MX_Exit("M2Mul", 2);
#endif

    while (k < Zeil2(i).anz)
    { /* Alle Eintraege durchsuchen */
        if (j == Ind2(i, k))
        {
            Aik2(i, k) *= aij_mul;
            return 0;
        } /* gefunden! */
        k++;
    }
    return 0;
}

/**** Modell 3, 4 *********************************************************/
int M34Mul(long i1, long k1, double aik)
{
    Modell34* w = (Modell34*)wurzel;
    register int j = 0;
    register long i = i1, k = k1;
    int u = 0;

#ifdef ERROR_CONTROL
    if ((i >= dim) || (k >= dim) || (i < 0) || (k < 0))
        MX_Exit("M34Mul", 2);
#endif

    w->stat = 1;
    if (i == k)
    {
        w->Diag[i] *= aik;
        return 0;
    } /* Diagonalelement */
    if (i > k)
    { /* Untermatrixelement */
        if (!(w->usym))
            return 0; /* symmetrische Matrix, ignorieren */
        u = 1;
        i = k;
        k = i1;
    }
    /* Alle Spalteneintraege durchsuchen */
    while (j < Sp34(k).anz && i >= Ind34(k, j))
    {
        if (i == Ind34(k, j))
        {
            Aik34(k, j, u) *= aik;
            return 0;
        }
        j++;
    }
    return 0; /* Element nicht vorhanden */
}

/*************************************************************************
   ROCKFLOW - Funktion: M#Div

   Aufgabe:
   Dividieren des Matrixwertes A[i,j] der zuvor mit
   MXSetMatrixPointer angegebenen Matrix-Struktur durch den Wert aij_div.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E long i, j: Indizes der Matrixposition
   E double aij_div: Divisor von A[i,j]

   Ergebnis:
   0

   Programmaenderungen:
   Modell 1, 2:   11/1995     MSR        Erste Version
   alle Modelle :  2/2000     Ra    Neufassung mit MXMul
*** nicht modellspezifisch! *********************************************/
int MX_Div(long i, long j, double aij_div)
{
    return MXMul(i, j, (double)1.0 / aij_div);
}

/*************************************************************************
   ROCKFLOW - Funktion: M#Get

   Aufgabe:
   Ergebnis ist der Matrixwert A[i,j] der zuvor mit MXSetMatrixPointer
   angegebenen Matrix-Struktur.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E long i, j: Indizes der Matrixposition

   Ergebnis:
   A[i,j]

   Programmaenderungen:
   Modell 1, 2:   11/1995     MSR        Erste Version
   Modell 3, 4:    2/2000     Ra         Letzte Version
*** Modell 1 ************************************************************/
double M1Get(long i, long j)
{
    Modell1* w = (Modell1*)wurzel;

#ifdef ERROR_CONTROL
    if ((i >= dim) || (j >= dim) || (i < 0) || (j < 0))
        MX_Exit("M1Get", 2);
#endif

    return Aik1(i, j);
}

/**** Modell 2 ************************************************************/
double M2Get(long i, long j)
{
    Modell2* w = (Modell2*)wurzel;
    register int k = 0;

#ifdef ERROR_CONTROL
    if ((i >= dim) || (j >= dim) || (i < 0) || (j < 0))
        MX_Exit("M2Get", 2);
#endif

    if ((j < Zeil2(i).min_col) || (j > Zeil2(i).max_col))
        return 0.;

    while (k < Zeil2(i).anz)
    { /* Alle Eintraege durchsuchen */
        if (j == Ind2(i, k))
            return Aik2(i, k); /* gefunden! */
        k++;
    }
    return 0.0; /* Kein Eintrag gefunden */
}

/**** Modell 3,4 **********************************************************/
double M34Get(long i1, long k1)
{
    Modell34* w = (Modell34*)wurzel;
    register int j = 0;
    register long i = i1, k = k1;
    int u = 0;

#ifdef ERROR_CONTROL
    if ((i >= dim) || (k >= dim) || (i < 0) || (k < 0))
        MX_Exit("M34Get", 2);
#endif

    if (i == k)
        return w->Diag[i]; /* Diagonalelement */
    if (i > k)
    {
        u = w->usym;
        i = k;
        k = i1;
    } /* Untermatrixelement */
    /* Alle Spalteneintraege durchsuchen */
    while (j < Sp34(k).anz && i >= Ind34(k, j))
    {
        if (i == Ind34(k, j))
            return Aik34(k, j, u);
        j++;
    }
    return 0.0; /* Element ist nicht vorhanden */
}

// WW/PA 08/02/2006
double M5Get(long i, long j)
{
    MeshLib::CFEMesh const* const mesh(FEMGet("GROUNDWATER_FLOW"));
    const long dim1 = mesh->NodesInUsage();

#ifdef ERROR_CONTROL
    if ((i >= dim1) || (j >= dim1) || (i < 0) || (j < 0))
        MX_Exit("M5Get", 2);
#endif

    const size_t jj(mesh->Eqs2Global_NodeIndex[j]);
    CNode const* const nod_i(mesh->nod_vector[mesh->Eqs2Global_NodeIndex[i]]);
    std::vector<size_t> const& connected_nodes(nod_i->getConnectedNodes());
    const size_t n_connected_nodes(connected_nodes.size());

    for (size_t k = 0; k < n_connected_nodes; k++)
        if (connected_nodes[k] == jj)
            // TEST WW		return  jdiag[m_nod_i->m5_index[k]];
            break;

    return 0.0; /* Kein Eintrag gefunden */
}

// WW/PA 08/02/2006
int M5Set(long i, long j, double e_val)
{
    long dim1;
    e_val = e_val;  // OK411
    MeshLib::CFEMesh const* const mesh(FEMGet("GROUNDWATER_FLOW"));
    //  CNode *m_nod_j = NULL;
    dim1 = mesh->NodesInUsage();

#ifdef ERROR_CONTROL
    if ((i >= dim1) || (j >= dim1) || (i < 0) || (j < 0))
        MX_Exit("M5Get", 2);
#endif

    const size_t jj(mesh->Eqs2Global_NodeIndex[j]);
    CNode const* const nod_i(mesh->nod_vector[mesh->Eqs2Global_NodeIndex[i]]);
    std::vector<size_t> const& connected_nodes(nod_i->getConnectedNodes());
    const size_t n_connected_nodes(connected_nodes.size());

    for (size_t k = 0; k < n_connected_nodes; k++)
        if (connected_nodes[k] == jj)
            // TEST WW			jdiag[m_nod_i->m5_index[k]] = e_val;
            break;

    return 0; /* Kein Eintrag gefunden */
}

/*************************************************************************
   ROCKFLOW - Funktion: MxCopyToAMG1R5Structure

   Aufgabe:
   Uebertraegt die RockFlow-Interne Speicherstruktur in die
   Struktur fuer den AMG1R5-Loeser

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)

   Zu komplex es hier darzulegen

   Ergebnis:

   0: korrekte Ausfuehrung

   Programmaenderungen:

   3/2002  C.Thorenz    Erste Version
   1/2004  TK			Compiler Warning: Unreferenced formal parameters deleted

 *****************************************************************************/
int M1CopyToAMG1R5Structure(double* A, int* IA, int* JA, int /*NDA*/,
                            int /*NDIA*/, int /*NDJA*/, double* x, double* U,
                            double* b, double* F)
{
    long i, j, NNA = 0, NIA = 0;
    double a;

    for (i = 0; i < dim; i++)
    {
        /* Vektoren umkopieren */
        U[i] = x[i];
        F[i] = b[i];

        /* Diagonale eintragen */
        j = i;
        a = MXGet(i, j);

        A[NNA] = a;
        JA[NNA] = j + 1;
        NNA++;
        IA[NIA] = NNA;
        NIA++;

        for (j = 0; j < dim; j++)
        {
            /* Restmatrix eintragen */
            a = MXGet(i, j);
            if ((fabs(a) > MKleinsteZahl) && (i != j))
            {
                A[NNA] = a;
                JA[NNA] = j + 1;
                NNA++;
            }
        }
    }

    IA[NIA] = NNA + 1;
    return 0;
}

int M2CopyToAMG1R5Structure(double* A, int* IA, int* JA, int /*NDA*/,
                            int /*NDIA*/, int /*NDJA*/, double* x, double* U,
                            double* b, double* F)
{
    long i, j;
    long NNA = 0, NIA = 0;
    double a;
    Modell2* w = (Modell2*)wurzel;
    register int k;

    for (i = 0; i < dim; i++)
    {
        /* Vektoren umkopieren */
        U[i] = x[i];
        F[i] = b[i];

        j = i;
        a = Aik2(i, 0);

        A[NNA] = a;
        JA[NNA] = j + 1;
        NNA++;
        IA[NIA] = NNA;
        NIA++;

        k = 1; /* Diagonale ist schon beruecksichtigt */
        while (k < Zeil2(i).anz)
        { /* Alle Eintraege durchsuchen */
            j = Ind2(i, k);
            a = Aik2(i, k);
            A[NNA] = a;
            JA[NNA] = j + 1;
            NNA++;
            k++;
        }
    }

    IA[NIA] = NNA + 1;
    return 0;
}

int M34CopyToAMG1R5Structure(double* A, int* IA, int* JA, int /*NDA*/,
                             int /*NDIA*/, int /*NDJA*/, double* x, double* U,
                             double* b, double* F)
{
    long i, j, NNA = 0, NIA = 0;
    double a;

    for (i = 0; i < dim; i++)
    {
        j = i;

        /* Vektoren umkopieren */
        U[i] = x[i];
        F[i] = b[i];

        a = MXGet(i, j);
        A[NNA] = a;
        JA[NNA] = j + 1;
        NNA++;
        IA[NIA] = NNA;
        NIA++;

        for (j = 0; j < dim; j++)
        {
            a = MXGet(i, j);
            if ((fabs(a) > MKleinsteZahl) && (i != j))
            {
                A[NNA] = a;
                JA[NNA] = j + 1;
                NNA++;
            }
        }
    }

    IA[NIA] = NNA + 1;
    return 0;
}

/*************************************************************************
   ROCKFLOW - Funktion: M#Trans

   Aufgabe:
   Fuehrt in der zuvor mit MXSetMatrixPointer angegebenen Matrix-Struktur
   folgende Operation durch: A[i,j] = A[ii,jj]

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E long i, j: Indizes der Matrixposition 1
   E long ii, jj: Indizes der Matrixposition 2

   Ergebnis:  - void -

   Programmaenderungen:
   Modell 1,2:    11/1995     MSR   Erste Version
   alle Modelle :  2/2000     Ra    Neufassung mit MXGet/MXSet
*** nicht modellspezifisch! *********************************************/
void MX_Trans(long i, long j, long ii, long jj)
{
    MXSet(i, j, MXGet(ii, jj));
}

/*************************************************************************
   ROCKFLOW - Funktion: M#MatVek

   Aufgabe:
   Ausfuehren des Matrix-Vektor-Produktes "A * vektor = ergebnis" mit der
   zuvor mit MXSetMatrixPointer angegebenen Matrix-Struktur. "ergebnis"
   wird komplett ueberschrieben; der Speicher muss bereits allokiert sein.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E double *vektor: mit der Matrix zu multiplizierender Vektor
   E double *ergebnis: Ergebnis der Multiplikation

   Ergebnis:
   - void -

   Programmaenderungen:
   Modell 1,2:   11/1995     MSR Erste Version
   Modell 3,4:    2/2000     Ra  Letzte Version
   Modell   2:    3/2002     CT  CBLAS eingehaengt
*** Modell 1 ************************************************************/
void M1MatVek(double* vektor, double* ergebnis)
{
    Modell1* w = (Modell1*)wurzel;
    register long i, k;

#ifdef ERROR_CONTROL
    if ((vektor == NULL) || (ergebnis == NULL))
        MX_Exit("M1MatVek", 3);
#endif
    for (i = 0; i < dim; i++)
    {
        ergebnis[i] = 0.0;
#ifdef SX
#pragma cdir nodep
#endif
        for (k = 0; k < dim; k++)
            ergebnis[i] += Aik1(i, k) * vektor[k];
    }
}

/****** Permutation when transforming from Modell2 to Modell 5 JAD ********/
/**** Modell 2 ************************************************************/
// WW/PA 08/02/2006
int M5Inc(long i, long j, double aij_inc)
{
    MeshLib::CFEMesh const* const mesh(FEMGet("GROUNDWATER_FLOW"));
    const long dim1(mesh->NodesInUsage());

#ifdef ERROR_CONTROL
    if ((i >= dim1) || (j >= dim1) || (i < 0) || (j < 0))
        MX_Exit("M5Inc", 2);
#endif
    if (fabs(aij_inc) < MKleinsteZahl)
        return 0; /* Abbruch bei Nullwert */

    const size_t jj(mesh->Eqs2Global_NodeIndex[j]);
    CNode const* const nod_i(mesh->nod_vector[mesh->Eqs2Global_NodeIndex[i]]);
    std::vector<size_t> const& connected_nodes(nod_i->getConnectedNodes());
    const size_t n_connected_nodes(connected_nodes.size());
    for (size_t k = 0; k < n_connected_nodes; k++)
        if (connected_nodes[k] == jj)
            // TEST WW			 jdiag[m_nod_i->m5_index[k]] += aij_inc;
            break;

    return 1;
}

// WW/PA 08/02/2006
/*
   void Write_Matrix_M5(double *b, ostream& os)
   {
   long i,j , dim1;
   MeshLib::CFEMesh* m_msh = NULL;
   m_msh = FEMGet("GROUNDWATER_FLOW");

   #ifdef SX
   #pragma cdir nodep
   #endif
   dim1 = m_msh->NodesInUsage();
   for (i = 0; i < dim1; i++)
   {
   for (j = 0; j < dim1; j++)
   {
   if(fabs(MXGet(i,j))>MKleinsteZahl)
   os<<i<<"  "<<j<<"  "<<MXGet(i,j)<<"\n";
   }
   if(fabs(b[i])>MKleinsteZahl)
   os<<i<<"  "<<dim1+1<<"  "<<b[i]<<"\n";   // os<<"\n";
   }
   }
 */

// void transM2toM5(void)
// void M5CreateMatrix(void)
// PA/WW 08/02/2006
void* M5CreateMatrix(long param1, long param2, long param3)
{
    param1 = param1;
    param3 = param3;
    //
    Modell5* w = (Modell5*)wurzel;
    w = (Modell5*)Malloc(sizeof(Modell5));
    MXSetMatrixPointer((void*)w);

    matrix_type = param2;

    int i, ii, jj, count1;
    long k, index, Size, dim1;  // OK
    /*------------------------------------------------------------*/
    jd_ptr_max = 0;
    MeshLib::CFEMesh* m_msh = NULL;
    m_msh = FEMGet("GROUNDWATER_FLOW");
#ifdef SX
#pragma cdir nodep
#endif
    dim1 = m_msh->NodesInUsage();
    Dim_L = dim1;
    for (k = 0; k < dim1; k++)
    {
        index = m_msh->Eqs2Global_NodeIndex[k];  //
        // WW
        Size = (int)m_msh->nod_vector[index]->getConnectedNodes().size();
        if (Size > jd_ptr_max)
            jd_ptr_max = Size;
#ifdef SX
#pragma cdir nodep
#endif
        for (i = 0; i < Size; i++)
            m_count++;
    }
    /*------------------------------------------------------------*/
    jd_ptr1 = (int*)malloc(dim1 * sizeof(int));
    jd_ptr2 = (int*)malloc(dim1 * sizeof(int));
    /*------------------------------------------------------------*/
    jd_ptr = (int*)malloc(jd_ptr_max * sizeof(int));
    jd_ptrf = (int*)malloc((jd_ptr_max + 1) * sizeof(int));
    temp_ergebnis = (double*)malloc(dim1 * sizeof(double));
    col_ind = (int*)malloc(m_count * sizeof(int));
    jdiag = (double*)malloc(m_count * sizeof(double));
    diag5_i = (long*)malloc(dim1 * sizeof(long));

    /*------------------------------------------------------------*/
    for (k = 0; k < jd_ptr_max; k++)
        jd_ptr[k] = 0;

    for (k = 0; k < dim1; k++)
    {
        index = m_msh->Eqs2Global_NodeIndex[k];  //
        // WW
        Size = (int)m_msh->nod_vector[index]->getConnectedNodes().size();
        for (i = 0; i < Size; i++)
            jd_ptr[i]++;
    }
    //	  printf("In transM2toM5 dim=%ld\n",dim1);
    for (k = 0; k < dim1; k++)
    {
        //	  printf("Zeil2(%d).anz=%d\n",k,Zeil2(k).anz);
        index = m_msh->Eqs2Global_NodeIndex[k];  //
        //       Zeil2(k).anz;
        jd_ptr1[k] = (int)m_msh->nod_vector[index]->getConnectedNodes().size();
        jd_ptr2[k] = k;
    }

    /*	  printf("Before insertionSort1_des\n");
         for(k=0;k<dim;k++)
         printf("jd_ptr1[%d]=%d\n",k,jd_ptr1[k]);
         for(k=0;k<dim;k++)
         printf("jd_ptr2[%d]=%d\n",k,jd_ptr2[k]);
     */
    insertionSort1_des(jd_ptr1, jd_ptr2, dim1);
    /*
         printf("After insertionSort1_des\n");
         for(k=0;k<dim;k++)
         printf("jd_ptr1[%d]=%d\n",k,jd_ptr1[k]);
         for(k=0;k<dim;k++)
         printf("jd_ptr2[%d]=%d\n",k,jd_ptr2[k]);
     */
    for (k = 0; k < jd_ptr_max + 1; k++)
        jd_ptrf[0] = 0;

    for (k = 0; k < jd_ptr_max; k++)
        jd_ptrf[k + 1] = jd_ptrf[k] + jd_ptr[k];

    //
    /*
          count1 = 0;
          for (k = 0; k < jd_ptr_max; k++)
          {
          for (i = 0; i < jd_ptr[k]; i++)
            {
            ii = jd_ptr2[i];
              jdiag[count1] = Aik2(ii, k);
       //	     col_ind[count1] = Ind2(ii, k);

       //TEST WW OUT
       cout<<ii<<"  "<< col_ind[count1]<<"    ";
       count1++;
       }
       cout<<"\n";
       }

       //TEST WW OUT
       cout<<"----------------------------"<<"\n";

     */

    count1 = 0;
    long col_i = 0;
    for (k = 0; k < jd_ptr_max; k++)
        for (i = 0; i < jd_ptr[k]; i++)
        {
            // Row of equation system
            ii = jd_ptr2[i];
            index = m_msh->Eqs2Global_NodeIndex[ii];
            //			 col_ind[count1] =
            // m_msh->go[m_msh->nod_vector[index]->connected_nodes[k]];
            col_i = m_msh->nod_vector[index]->getConnectedNodes()[k];
            col_ind[count1] = m_msh->nod_vector[col_i]->GetEquationIndex();
            jj = col_ind[count1];
            // TEST WW m_msh->nod_vector[index]->m5_index[k]=count1;
            if (ii == jj)
                diag5_i[ii] = count1;
            // TEST WW OUT
            //			 cout<<" Row: "<< ii<<"  "<< col_ind[count1]<<"    ";
            //             cout<<ii<<"    ";
            //             cout<< col_ind[count1]<<"    ";
            //////////////////////////////////

            count1++;
        }
    // TEST WW
    //          cout<<"\n";
    ////////////////////
    return (void*)w;
}

/*-------------------------------------------------------------*/
#ifdef SX
void insertionSort1_des(int* restrict numbers, int* restrict numbers1,
                        int array_size)
#else
void insertionSort1_des(int* numbers, int* numbers1, int array_size)
#endif
{
    int i, j, index, index1;
#ifdef SX
#pragma cdir nodep
#endif
    for (i = 1; i < array_size; i++)
    {
        index = numbers[i];
        index1 = numbers1[i];
        j = i;
        while ((j > 0) && (numbers[j - 1] < index))
        {
            numbers[j] = numbers[j - 1];
            numbers1[j] = numbers1[j - 1];
            j = j - 1;
        }
        numbers[j] = index;
        numbers1[j] = index1;
    }
} /*End of insertionSort1_des.*/

/****** Transforming from Modell2 to Modell 6 ITPACKV ********/
void transM2toM6(void)
{
    Modell2* w = (Modell2*)wurzel;
    int i, k, count1, dim1;

    dim1 = ((long*)wurzel)[0];
    /*------------------------------------  TEST ----------------------
    *
       for (k = 0; k < dim1; k++)
       {
          printf("Zeil2(%d).anz=%d\n",k, Zeil2(k).anz);
          for (i = 0; i < Zeil2(k).anz; i++)
          {

          printf("Ind2(%d,%d) = %d\n",k,i,Ind2(k, i));
         }
       }
       for (k = 0; k < dim1; k++)
       {
       printf("Zeil2(%d).anz=%d\n",k, Zeil2(k).anz);
       for (i = 0; i < Zeil2(k).anz; i++)
       {

       printf("Aik2(%d,%d) = %le\n",k,i,Aik2(k, i));
       }
       }

       ------------------------------------  TEST ----------------------*/

    /*------------------------------------------------------------*/
    itpack_max = 0;

#ifdef SX
#pragma cdir nodep
#endif
    for (k = 0; k < dim1; k++)
        if (Zeil2(k).anz > itpack_max)
            itpack_max = Zeil2(k).anz;

    itpackv = (double*)malloc(itpack_max * dim1 * sizeof(double));
    it_col = (int*)malloc(itpack_max * dim1 * sizeof(int));
    //	  temp_ergebnis = (double *)malloc(dim*sizeof(double));

    count1 = 0;
    for (i = 0; i < itpack_max; i++)
        for (k = 0; k < dim1; k++)
        {
            if (i < Zeil2(k).anz)
            {
                itpackv[count1] = Aik2(k, i);
                it_col[count1] = Ind2(k, i);
                count1++;
            }
            else
            {
                itpackv[count1] = 0.0;
                it_col[count1] = Ind2(0, 0);
                count1++;
            }
        }
    printf("AT End of transM2toM6\n");
    /*---------------------------------------------------------
         for(k=0; k<count1; k++)
         {
          printf("itpackv[%d]=%le\n",k,itpackv[k]);
         }
         for(k=0; k<count1; k++)
         {
          printf("it_col[%d]=%d\n",k,it_col[k]);
         }
     */
}

/**** Modell 2 ************************************************************/
void M2MatVek(double* vektor, double* ergebnis)
{
    Modell2* w = (Modell2*)wurzel;
    register long k;
    register int i;
#ifdef CBLAS_M2MatVek
    double* help;

    help = (double*)Malloc(dim * sizeof(double));
#endif

#ifdef ERROR_CONTROL
    if ((vektor == NULL) || (ergebnis == NULL))
        MX_Exit("M2MatVek", 3);
#endif

    for (k = 0; k < dim; k++)
    {
#ifndef CBLAS_M2MatVek
        ergebnis[k] = 0.0;
#ifdef SX
#pragma cdir nodep
#endif
        for (i = 0; i < Zeil2(k).anz; i++)
            ergebnis[k] += Aik2(k, i) * vektor[Ind2(k, i)];

#else
        for (i = 0; i < Zeil2(k).anz; i++)
            help[i] = vektor[Ind2(k, i)];

        ergebnis[k] = cblas_ddot(Zeil2(k).anz, Zeil2(k).wert, 1, help, 1);
#endif
    }

#ifdef CBLAS_M2MatVek
    help = (double*)Free(help);
#endif
}

/*--------------------------------------------------------------------
 * Matrxi-Vektor Multiply with Jagged Diagonal Format (Modell 5)
   PA Initial version
   ---------------------------------------------------------------------*/
#ifdef SX
void M5MatVek(double* restrict b, double* restrict erg)
#else
void M5MatVek(double* b, double* erg)
#endif
{
    //  Modell2 *w = (Modell2 *) wurzel;
    int i, dim1;  // k,
    int j, num;
    long col_len;
    MeshLib::CFEMesh* m_msh = NULL;  // WW
    m_msh = FEMGet("GROUNDWATER_FLOW");
    dim1 = m_msh->NodesInUsage();

#ifdef SX
#pragma cdir nodep
#endif
    for (i = 0; i < dim1; i++)
        erg[i] = 0.0;

    for (i = 0; i < dim1; i++)
        temp_ergebnis[i] = 0.0;

    for (i = 0; i < jd_ptr_max; i++)
    {
        col_len = jd_ptrf[i + 1] - jd_ptrf[i];
        num = jd_ptrf[i];
        //	  printf("num=%d\n",num);
        //	  printf("col_len=%d\n",col_len);

#ifdef _OPENMP
#pragma omp parallel for private(j) \
    shared(temp_ergebnis, jdiag, b, num, col_ind)
#endif
#ifdef SX
#pragma cdir nodep
#endif
        for (j = 0; j < col_len; j++)
            temp_ergebnis[j] += jdiag[num + j] * b[col_ind[num + j]];
    }

#ifdef SX
#pragma cdir nodep
#endif
    for (i = 0; i < dim1; i++)
        erg[jd_ptr2[i]] = temp_ergebnis[i];
}

/*--------------------------------------------------------------------
* Matrxi-Vektor Multiply with ITPACKV
   --------------------------------------------------------------------*/
#ifdef SX
void M6MatVek(double* restrict b, double* restrict erg)
#else
void M6MatVek(double* b, double* erg)
#endif
{
    //  Modell2 *w = (Modell2 *) wurzel;
    int i, j, num, dim1;

    dim1 = ((long*)wurzel)[0];

#ifdef SX
#pragma cdir nodep
#endif
    for (i = 0; i < dim1; i++)
        erg[i] = 0.0;

    num = 0;
    //  printf("In M6MatVek itpack_max=%d\n",itpack_max);
    for (i = 0; i < itpack_max; i++)
    {
//#pragma omp parallel for private(j) shared(erg,itpackv,b,num,it_col)
#ifdef SX
#pragma cdir nodep
#endif
        for (j = 0; j < dim1; j++)
            erg[j] = erg[j] + itpackv[num + j] * b[it_col[num + j]];
        num = num + dim1;
    }
}

/**** Modell 3,4 **********************************************************/
void H_MX34_mul(double* x, double* r, int o)
/* Hilfsroutine zu Matrix-Vektorprodukt Modelle 3 und 4
   Ergebnis: r = A*x  bei o=0
   bzw.      r = At*x bei o=1
   A = Matrix Modell 3/4
   Die Vektoren x und r muessen allokiert sein!
 */
{
    Modell34* w = (Modell34*)wurzel;
    register long k, i;
    register int j;
    int u = 1 - o; /* o=0: normal, o=1: mit Transponierter */

#ifdef ERROR_CONTROL
    if ((x == NULL) || (r == NULL))
        MX_Exit("M34Mat(T)Vek", 3);
#endif

    if (!(w->usym))
    {
        o = 0;
        u = 0;
    } /* symmetrische Matrix: Oik=Uki */
    for (k = 0; k < dim; k++)
    {
        r[k] = w->Diag[k] * x[k];
#ifdef SX
#pragma cdir nodep
#endif
        for (j = 0; j < Sp34(k).anz; j++)
        {
            i = Ind34(k, j);
            r[i] += Aik34(k, j, o) * x[k];
            r[k] += Aik34(k, j, u) * x[i];
        }
    }
}

void M34MatVek(double* vektor, double* ergebnis)
{
    H_MX34_mul(vektor, ergebnis, (int)0);
}

/*************************************************************************
   ROCKFLOW - Funktion: M#MatTVek

   Aufgabe:
   Ausfuehren des Matrix-Vektor-Produktes "A^T * vektor = ergebnis" mit der
   zuvor mit MXSetMatrixPointer angegebenen Matrix-Struktur. "ergebnis"
   wird komplett ueberschrieben; der Speicher muss bereits allokiert sein.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E double *vektor: mit der transponierten Matrix zu multiplizierender
   Vektor
   E double *ergebnis: Ergebnis der Multiplikation

   Ergebnis:
   - void -

   Programmaenderungen:
   Modell 1,2:   11/1995     MSR        Erste Version
   Modell 3,4:    2/2000     Ra         Letzte Version
 *** Modell 1 ***********************************************************/
void M1MatTVek(double* vektor, double* ergebnis)
{
    Modell1* w = (Modell1*)wurzel;
    register long i, k;

#ifdef ERROR_CONTROL
    if ((vektor == NULL) || (ergebnis == NULL))
        MX_Exit("M1MatTVek", 3);
#endif

    for (i = 0; i < dim; i++)
    {
        ergebnis[i] = 0.0;
        for (k = 0; k < dim; k++)
            ergebnis[i] += Aik1(k, i) * vektor[k];
    }
}

/**** Modell 2 ***********************************************************/
void M2MatTVek(double* vektor, double* ergebnis)
{
    Modell2* w = (Modell2*)wurzel;
    register long k;
    register int i;

#ifdef ERROR_CONTROL
    if ((vektor == NULL) || (ergebnis == NULL))
        MX_Exit("M2MatTVek", 3);
#endif

    for (k = 0; k < dim; k++)
        ergebnis[k] = 0.0;
    for (k = 0; k < dim; k++)
        for (i = 0; i < Zeil2(k).anz; i++)
            ergebnis[Ind2(k, i)] += Aik2(k, i);
}

/**************************************************************************/
void M34MatTVek(double* vektor, double* ergebnis)
{
    H_MX34_mul(vektor, ergebnis, (int)1); /* siehe M34MatVek */
}

/*************************************************************************
   ROCKFLOW - Funktion: MXResiduum

   Aufgabe:
   Berechnet das Residuum "r = b - A x"  mit der zuvor mit
   MXSetMatrixPointer angegebenen Matrix-Struktur. "ergebnis" wird
   komplett ueberschrieben; der Speicher muss bereits allokiert sein.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E double *x: Vektor x
   E double *b: Vektor b
   E double *ergebnis: Residuum
   Ueberschreiben: b und ergebnis duerfen _nicht_ gleich sein (CT)       */

/* Ergebnis:
   - void -

   Programmaenderungen:
   Modell 1,2:  11/1995   MSR        Erste Version
                10/1998   C.Thorenz  Hilfsvariable gegen kleine Differenzen
                                     grosser Zahlen
   Modell 3,4:   2/2000   Ra         Letzte Version
   Alle Modelle: 3/2002   CT         Vereinfacht und alle Modelle
   zusammengefasst
 */

void MXResiduum(double* x, double* b, double* ergebnis)
{
#ifdef ERROR_CONTROL
    if ((x == NULL) || (b == NULL) || (ergebnis == NULL))
        MX_Exit("MXResiduum", 3);
    if (b == ergebnis)
        MX_Exit("MXResiduum", 5);
#endif

    MXMatVek(x, ergebnis);

    // WW
    long Dimension = 0;
    if (matrix_type == 5)
        Dimension = Dim_L;
    else
        Dimension = dim;
    //
    // WW
    MAddSkalVektoren(ergebnis, -1., b, 1., ergebnis, Dimension);
}

/*************************************************************************
   ROCKFLOW - Funktion: MXRandbed

   Aufgabe:
   Vorgegebene Randbedingung Xi= Ri in Matrix und Rechte Seite einbauen

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E long ir  : Gleichungsnummer, fuer die Randbedingung vorgegeben wird
   E double Ri: vorgegebener Wert der Randbedingung
   E double *ReSei: Zeiger auf Rechte Seite

   Ergebnis:  - void -

   Programmaenderungen:
   alle Modelle:    3/2000     Ra        Erste Version
   6/2000     C.Thorenz Fehler in Modell 34 beseitigt
   1/2001     C.Thorenz Diagonale bleibt auf Wert, stattdessen wird
   der rechte Seite Eintrag skaliert
   => keine boese Veraenderung an einer Stelle der Diagonale
   12/2001     C.Thorenz Diagonalwert = 0 abgefangen
   02/2006      WW/PA M5 storage
*** Unterscheidung nach Modell innerhalb der Prozedur! ******************/
void MXRandbed(long ir, double Ri, double* ReSei)
{
    register long i = 0, ip, im;
    long p, q;
    int k;
    double diag;

    /* Diagonalelement holen */
    diag = MXGet(ir, ir);

    /* Wenn das Diagonalelement ~= 0 ist, ein anderes suchen.
       Sollte eigentlich nicht vorkommen. */

    if (fabs(diag) < DBL_MIN)  // TODO This probably won't work as intended
        for (i = 0; i < dim; i++)
        {
            ip = ir + i;
            im = ir - i;
            if (ip < dim)
                if (fabs(diag = MXGet(ip, ip)) > DBL_MIN)
                    break;
            if (im >= 0)
                if (fabs(diag = MXGet(im, im)) > DBL_MIN)
                    break;
        }

    /* Dieser Fall sollte _nie_ auftreten */
    if (fabs(diag) < DBL_MIN)
        diag = MKleinsteZahl;

    // TEST WW
    /*
       #ifdef ERROR_CONTROL
       if ((ir >= dim) || (ir < 0))
        MX_Exit("MXRandbed", 2);
       if (ReSei == NULL)
        MX_Exit("MXRandbed", 3);
       #endif
     */
    switch (matrix_type)
    {
        case 1:
        {
            Modell1* w = (Modell1*)wurzel;
            p = ir - w->NumDif;
            if (p < 0)
                p = 0;
            q = ir + w->NumDif + 1;
            if (q > dim)
                q = dim;
            for (i = p; i < q; i++)
            {
                ReSei[i] -= Aik1(i, ir) * Ri; /* Spalte ausmultiplizieren */
                Aik1(i, ir) = Aik1(ir, i) =
                    0.0; /* Spalte und Zeile Null setzen */
            }
            ReSei[ir] = Ri * diag; /* Randwert einsetzen */
            MXSet(ir, ir, diag);
        }
        break;

        case 2:
        {
            Modell2* w = (Modell2*)wurzel;
            p = ir - w->NumDif;
            if (p < 0)
                p = 0;
            q = ir + w->NumDif + 1;
            if (q > dim)
                q = dim;
            for (i = p; i < q; i++)
                if (i != ir)
                { /* alle anderen Zeilen mit Spalte iR */
                    k = 0;
                    while (++k <
                           Zeil2(i).anz) /* Alle Eintraege (ausser Diag.) */
                        if (ir == Ind2(i, k))
                        {
                            ReSei[i] -= Aik2(i, k) * Ri;
                            Aik2(i, k) = 0.0;
                            goto End_Zeil;
                        } /*[i,iR] */
                End_Zeil:;
                }

            for (k = 0; k < Zeil2(ir).anz; k++)
                Aik2(ir, k) = 0.0; /* Rest der Zeile Null */
            ReSei[ir] = Ri * diag; /* Randwert einsetzen */
            MXSet(ir, ir, diag);
        }
        break;

        case 3:
        case 4:
        {
            Modell34* w = (Modell34*)wurzel;
            int u = w->usym;
            for (k = 0; k < Sp34(ir).anz; k++)
            { /* ueber (links von) Diag. */
                /* Obermatrix ausmultipl. */
                ReSei[Ind34(ir, k)] -= Aik34(ir, k, 0) * Ri;
                /* Zeile / Spalte Null */
                Aik34(ir, k, 0) = Aik34(ir, k, u) = 0.0;
            }

            MXSet(ir, ir, diag); /* Randwert einsetzen */
            ReSei[ir] = Ri * diag;

            for (i = ir + 1; i <= Sp34(ir).rechts; i++)
                for (k = 0; k < Sp34(i).anz; k++)
                {
                    if (Ind34(i, k) > ir)
                        break;
                    if (Ind34(i, k) == ir)
                    {
                        ReSei[i] -=
                            Aik34(i, k, u) * Ri; /* Untermatrix ausmult. */
                        /* Zeile / Spalte Null */
                        Aik34(i, k, 0) = Aik34(i, k, u) = 0.0;
                    }
                }
            w->stat = 1; /* Matrix veraendert! */
        }
        break;
        case 5:  // WW/PA  08/02/2006
        {
            // WW long dim1;
            long jr;
            MeshLib::CFEMesh const* const mesh(FEMGet("GROUNDWATER_FLOW"));
            CNode* m_nod_j = NULL;
            // WW dim1 = m_msh->NodesInUsage();
            ReSei[ir] = Ri * MXGet(ir, ir);

            CNode const* const nod_i(
                mesh->nod_vector[mesh->Eqs2Global_NodeIndex
                                     [i]]);  // TODO check this. i could be 0.
            std::vector<size_t> const& connected_nodes(
                nod_i->getConnectedNodes());
            const size_t n_connected_nodes(connected_nodes.size());

            for (size_t k = 0; k < n_connected_nodes; k++)
            {
                m_nod_j = mesh->nod_vector[connected_nodes[k]];
                jr = m_nod_j->GetEquationIndex();
                if (ir == jr)
                    continue;
                MXSet(ir, jr, 0.0);
                ReSei[jr] -= MXGet(jr, ir) * Ri;
                MXSet(jr, ir, 0.0);
            }
        }
        break;
        case 6:
        {
            Modell2* w = (Modell2*)wurzel;
            p = ir - w->NumDif;
            if (p < 0)
                p = 0;
            q = ir + w->NumDif + 1;
            if (q > dim)
                q = dim;
            for (i = p; i < q; i++)
                if (i != ir)
                { /* alle anderen Zeilen mit Spalte iR */
                    k = 0;
                    while (++k <
                           Zeil2(i).anz) /* Alle Eintraege (ausser Diag.) */
                        if (ir == Ind2(i, k))
                        {
                            ReSei[i] -= Aik2(i, k) * Ri;
                            Aik2(i, k) = 0.0;
                            goto End_Zeil6;
                        } /*[i,iR] */
                End_Zeil6:;
                }

            for (k = 0; k < Zeil2(ir).anz; k++)
                Aik2(ir, k) = 0.0; /* Rest der Zeile Null */
            ReSei[ir] = Ri * diag; /* Randwert einsetzen */
            MXSet(ir, ir, diag);
        }
        break;

#ifdef ERROR_CONTROL
        default:
            MX_Exit("MXRandbed", 4);
#endif
    }
} /* MXRandbed */

/*************************************************************************
   ROCKFLOW - Funktion: MXEliminateIrrNode

   Aufgabe:
   Eliminieren eines irr. Knotens aus dem fertig aufgestellten
   Gesamtgleichungssystem und Eintragen der Zwangsbedingungen

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E long index: Index des irr. Knotens
   E int anz_nachbarn: Anzahl (2 oder 4) der reg. Nachbarknoten
   E long *nachbarn_index: Feld mit den Indizes der Nachbarn
   E double *rechts: Rechte Seite Vektor

   Ergebnis:  - void -

   Programmaenderungen:
   9/2000  C.Thorenz  Erste Version
   1/2001  C.Thorenz  Diagonale bleibt auf Wert, stattdessen werden
   die Zwangsbedingungen skaliert
   => keine boese Veraenderung an einer Stelle der Diagonale

 ************************************************************************/
void MXEliminateIrrNode(long index, int anz_nachbarn, long* nachbarn_index,
                        double* rechts)
{
    long zeile, spalte, j, k, l;
    double inc, gewicht = 1. / (double)anz_nachbarn, diag = MXGet(index, index);

#ifdef ERROR_CONTROL
    if ((index >= dim) || (index < 0))
        MX_Exit("MXEliminateIrrNode", 2);
    if (rechts == NULL)
        MX_Exit("MXEliminateIrrNode", 3);
#endif

    switch (matrix_type)
    {
        default:
        {
            /* Spalte des irr. Knotens eliminieren : */
            /* Ueber alle Zeilen gehen, um nach Eintragen fuer den Knoten zu
             * suchen */
            for (j = 0l; j < dim; j++)
            {
                /* Der Knoteneintrag wird ausmultipliziert und den Nachbarn
                 * zugeschlagen */
                /* index ist der irregulaere Knoten, j ist die Zeile,
                 * nachbarn_index die Spalte */
                inc = gewicht * MXGet(j, index);
                if (fabs(inc) > 0.)
                    for (k = 0; k < anz_nachbarn; k++)
                        MXInc(j, nachbarn_index[k], inc);

                /* Spalte des Knotens zu Null setzen, da der Wert oben verteilt
                 * wurde */
                if (fabs(MXGet(j, index)) > 0.)
                    MXSet(j, index, 0.);
            } /* endfor j */

#ifdef DUMP
            MXDumpGLS("a_spalte", 1, rechts, rechts);
#endif

            /* Zeile des irr. Knotens eliminieren : */
            /* Ueber Zeilen der regulaeren Nachbarn gehen  */
            /* Der Knoteneintrag wird ausmultipliziert und den Nachbarn
             * zugeschlagen */
            for (k = 0; k < dim; k++)
            {
                inc = gewicht * MXGet(index, k);
                if (fabs(inc) > 0.)
                {
                    for (j = 0l; j < anz_nachbarn; j++)
                        MXInc(nachbarn_index[j], k, inc);
                    /* Zeile des Knotens zu Null setzen, da hier die
                     * Interpolation eingesetzt wird */
                    MXSet(index, k, 0.);
                } /* endif */
            }

            /* Der rechte Seite-Eintrag wird den Nachbarn zugeschlagen */
            for (j = 0l; j < anz_nachbarn; j++)
                rechts[nachbarn_index[j]] += gewicht * rechts[index];

#ifdef DUMP
            MXDumpGLS("a_zeile", 1, rechts, rechts);
#endif

            /* Neuer Zeileneintrag fuer den irr. Knoten, dies ist
               die Zwangsbedingung aus der Interpolation */

            MXSet(index, index, diag);
            rechts[index] = 0.;
            for (k = 0; k < anz_nachbarn; k++)
                /* i ist die Zeile des irregulaere Knoten, nachbarn_index die
                 * Spalte */
                MXSet(index, nachbarn_index[k], -gewicht * diag);

#ifdef DUMP
            MXDumpGLS("a_zwang", 1, rechts, rechts);
#endif
        }
        break;

        case 2:
        {
            Modell2* w = (Modell2*)wurzel;

            /* Spalte des irr. Knotens eliminieren : */
            /* Ueber alle Zeilen gehen, um nach Eintragen fuer den Knoten zu
             * suchen */
            for (zeile = 0l; zeile < dim; zeile++)
            {
                /* Der Knoteneintrag wird ausmultipliziert und den Nachbarn
                 * zugeschlagen */
                /* index ist der irregulaere Knoten, nachbarn_index[] die Spalte
                 */

                spalte = index;
                k = 0;

                while (k < Zeil2(zeile).anz)
                { /* Alle Eintraege durchsuchen */
                    if (spalte == Ind2(zeile, k))
                    { /* gefunden! */
                        inc =
                            gewicht *
                            Aik2(zeile, k); /* Mit dem Gewicht multiplizieren */
                        Aik2(zeile, k) = 0.; /* Eigenen Eintrag loeschen */
                        for (l = 0; l < anz_nachbarn; l++)
                            /* Den reg. Nachbarn zuschlagen */
                            MXInc(zeile, nachbarn_index[l], inc);
                        break;
                    }
                    else
                        k++;
                }
            } /* endfor j */

            /* Zeile des irr. Knotens eliminieren : */
            /* Ueber Zeilen der regulaeren Nachbarn gehen  */
            for (spalte = 0; spalte < dim; spalte++)
            {
                zeile = index;
                k = 0;

                while (k < Zeil2(zeile).anz)
                { /* Alle Eintraege durchsuchen */
                    if (spalte == Ind2(zeile, k))
                    { /* gefunden! */
                        inc =
                            gewicht *
                            Aik2(zeile, k); /* Mit dem Gewicht multiplizieren */
                        Aik2(zeile, k) = 0.; /* Eigenen Eintrag loeschen */
                        for (l = 0; l < anz_nachbarn; l++)
                            /* Den reg. Nachbarn zuschlagen */
                            MXInc(nachbarn_index[l], spalte, inc);
                        break;
                    }
                    else
                        k++;
                } /* endwhile */
            }     /* endfor */

            for (j = 0l; j < anz_nachbarn; j++)
                /* Der rechte Seite-Eintrag wird den Nachbarn zugeschlagen */
                rechts[nachbarn_index[j]] += gewicht * rechts[index];
            /* endfor j */

            /* Neuer Zeileneintrag fuer den irr. Knoten, dies ist
               die Zwangsbedingung aus der Interpolation */
            MXSet(index, index, diag);

            for (k = 0; k < anz_nachbarn; k++)
                /* i ist die Zeile des irregulaere Knoten, nachbarn_index die
                 * Spalte */
                MXSet(index, nachbarn_index[k], -gewicht * diag);

            rechts[index] = 0.;
        }
        break;

        case 3:
        case 4:
        {
            Modell34* w = (Modell34*)wurzel;
            long i34, j34, k34;
            int u;

            /* Spalte des irr. Knotens eliminieren : */
            /* Ueber alle Zeilen gehen, um nach Eintragen fuer den Knoten zu
             * suchen */
            for (j = 0l; j < dim; j++)
            {
                /* Der Knoteneintrag wird ausmultipliziert und den Nachbarn
                 * zugeschlagen */
                inc = 0.;

                if (j == index)
                {
                    inc = gewicht * (w->Diag[j]); /* Diagonalelement */
                    w->Diag[j] = 0.;
                }
                else
                {
                    /* Variablen setzen */
                    j34 = 0;

                    /* Untermatrixelement? */
                    if (j > index)
                    {
                        u = w->usym;
                        i34 = index;
                        k34 = j;
                    }
                    else
                    {
                        u = 0;
                        i34 = j;
                        k34 = index;
                    }
                    /* Alle Spalteneintraege durchsuchen */
                    while (j34 < Sp34(k34).anz && i34 >= Ind34(k34, j34))
                    {
                        if (i34 == Ind34(k34, j34))
                        {
                            inc = gewicht * Aik34(k34, j34, u);
                            /* Spalte des Knotens zu Null setzen, da der Wert
                             * mit "inc" verteilt wird */
                            Aik34(k34, j34, u) = 0.;
                            break;
                        }
                        j34++;
                    }
                } /* endif */

                if (fabs(inc) > 0.)
                    for (k = 0; k < anz_nachbarn; k++)
                        MXInc(j, nachbarn_index[k], inc);
            } /* endfor j */

#ifdef DUMP
            MXDumpGLS("b_spalte", 1, rechts, rechts);
#endif

            /* Zeile des irr. Knotens eliminieren : */
            /* Ueber Zeilen der regulaeren Nachbarn gehen  */
            /* Der Knoteneintrag wird ausmultipliziert und den Nachbarn
             * zugeschlagen */
            for (k = 0; k < dim; k++)
            {
                inc = 0.;

                if (index == k)
                {
                    inc = gewicht * (w->Diag[k]); /* Diagonalelement */
                    w->Diag[k] = 0.;
                }
                else
                {
                    /* Variablen setzen */
                    j34 = 0;

                    /* Untermatrixelement */
                    if (index > k)
                    {
                        u = w->usym;
                        i34 = k;
                        k34 = index;
                    }
                    else
                    {
                        u = 0;
                        i34 = index;
                        k34 = k;
                    }

                    /* Alle Spalteneintraege durchsuchen */
                    while (j34 < Sp34(k34).anz && i34 >= Ind34(k34, j34))
                    {
                        if (i34 == Ind34(k34, j34))
                        {
                            inc = gewicht * Aik34(k34, j34, u);
                            /* Spalte des Knotens zu Null setzen, da der Wert
                             * mit "inc" verteilt wird */
                            Aik34(k34, j34, u) = 0.;
                            for (j = 0l; j < anz_nachbarn; j++)
                                MXInc(nachbarn_index[j], k, inc);
                            break;
                        }
                        j34++;
                    } /* endwhile */
                }     /* endif */
            }         /* endfor */

            /* Der rechte Seite-Eintrag wird den Nachbarn zugeschlagen */
            for (j = 0l; j < anz_nachbarn; j++)
                rechts[nachbarn_index[j]] += gewicht * rechts[index];
                /* endfor j */

#ifdef DUMP
            MXDumpGLS("b_zeile", 1, rechts, rechts);
#endif

            /* Neuer Zeileneintrag fuer den irr. Knoten, dies ist
               die Zwangsbedingung aus der Interpolation */

            MXSet(index, index, diag);
            rechts[index] = 0.;
            for (k = 0; k < anz_nachbarn; k++)
                /* i ist die Zeile des irregulaere Knoten, nachbarn_index die
                 * Spalte */
                MXSet(index, nachbarn_index[k], -gewicht * diag);

#ifdef DUMP
            MXDumpGLS("b_zwang", 1, rechts, rechts);
#endif

            break;
        }
        case 5:
        {
            Modell2* w = (Modell2*)wurzel;

            /* Spalte des irr. Knotens eliminieren : */
            /* Ueber alle Zeilen gehen, um nach Eintragen fuer den Knoten zu
             * suchen */
            for (zeile = 0l; zeile < dim; zeile++)
            {
                /* Der Knoteneintrag wird ausmultipliziert und den Nachbarn
                 * zugeschlagen */
                /* index ist der irregulaere Knoten, nachbarn_index[] die Spalte
                 */

                spalte = index;
                k = 0;

                while (k < Zeil2(zeile).anz)
                { /* Alle Eintraege durchsuchen */
                    if (spalte == Ind2(zeile, k))
                    { /* gefunden! */
                        inc =
                            gewicht *
                            Aik2(zeile, k); /* Mit dem Gewicht multiplizieren */
                        Aik2(zeile, k) = 0.; /* Eigenen Eintrag loeschen */
                        for (l = 0; l < anz_nachbarn; l++)
                            /* Den reg. Nachbarn zuschlagen */
                            MXInc(zeile, nachbarn_index[l], inc);
                        break;
                    }
                    else
                        k++;
                }
            } /* endfor j */

            /* Zeile des irr. Knotens eliminieren : */
            /* Ueber Zeilen der regulaeren Nachbarn gehen  */
            for (spalte = 0; spalte < dim; spalte++)
            {
                zeile = index;
                k = 0;

                while (k < Zeil2(zeile).anz)
                { /* Alle Eintraege durchsuchen */
                    if (spalte == Ind2(zeile, k))
                    { /* gefunden! */
                        inc =
                            gewicht *
                            Aik2(zeile, k); /* Mit dem Gewicht multiplizieren */
                        Aik2(zeile, k) = 0.; /* Eigenen Eintrag loeschen */
                        for (l = 0; l < anz_nachbarn; l++)
                            /* Den reg. Nachbarn zuschlagen */
                            MXInc(nachbarn_index[l], spalte, inc);
                        break;
                    }
                    else
                        k++;
                } /* endwhile */
            }     /* endfor */

            for (j = 0l; j < anz_nachbarn; j++)
                /* Der rechte Seite-Eintrag wird den Nachbarn zugeschlagen */
                rechts[nachbarn_index[j]] += gewicht * rechts[index];
            /* endfor j */

            /* Neuer Zeileneintrag fuer den irr. Knoten, dies ist
               die Zwangsbedingung aus der Interpolation */
            MXSet(index, index, diag);

            for (k = 0; k < anz_nachbarn; k++)
                /* i ist die Zeile des irregulaere Knoten, nachbarn_index die
                 * Spalte */
                MXSet(index, nachbarn_index[k], -gewicht * diag);

            rechts[index] = 0.;
        }
        break;
        case 6:
        {
            Modell2* w = (Modell2*)wurzel;

            /* Spalte des irr. Knotens eliminieren : */
            /* Ueber alle Zeilen gehen, um nach Eintragen fuer den Knoten zu
             * suchen */
            for (zeile = 0l; zeile < dim; zeile++)
            {
                /* Der Knoteneintrag wird ausmultipliziert und den Nachbarn
                 * zugeschlagen */
                /* index ist der irregulaere Knoten, nachbarn_index[] die Spalte
                 */

                spalte = index;
                k = 0;

                while (k < Zeil2(zeile).anz)
                { /* Alle Eintraege durchsuchen */
                    if (spalte == Ind2(zeile, k))
                    { /* gefunden! */
                        inc =
                            gewicht *
                            Aik2(zeile, k); /* Mit dem Gewicht multiplizieren */
                        Aik2(zeile, k) = 0.; /* Eigenen Eintrag loeschen */
                        for (l = 0; l < anz_nachbarn; l++)
                            /* Den reg. Nachbarn zuschlagen */
                            MXInc(zeile, nachbarn_index[l], inc);
                        break;
                    }
                    else
                        k++;
                }
            } /* endfor j */

            /* Zeile des irr. Knotens eliminieren : */
            /* Ueber Zeilen der regulaeren Nachbarn gehen  */
            for (spalte = 0; spalte < dim; spalte++)
            {
                zeile = index;
                k = 0;

                while (k < Zeil2(zeile).anz)
                { /* Alle Eintraege durchsuchen */
                    if (spalte == Ind2(zeile, k))
                    { /* gefunden! */
                        inc =
                            gewicht *
                            Aik2(zeile, k); /* Mit dem Gewicht multiplizieren */
                        Aik2(zeile, k) = 0.; /* Eigenen Eintrag loeschen */
                        for (l = 0; l < anz_nachbarn; l++)
                            /* Den reg. Nachbarn zuschlagen */
                            MXInc(nachbarn_index[l], spalte, inc);
                        break;
                    }
                    else
                        k++;
                } /* endwhile */
            }     /* endfor */

            for (j = 0l; j < anz_nachbarn; j++)
                /* Der rechte Seite-Eintrag wird den Nachbarn zugeschlagen */
                rechts[nachbarn_index[j]] += gewicht * rechts[index];
            /* endfor j */

            /* Neuer Zeileneintrag fuer den irr. Knoten, dies ist
               die Zwangsbedingung aus der Interpolation */
            MXSet(index, index, diag);

            for (k = 0; k < anz_nachbarn; k++)
                /* i ist die Zeile des irregulaere Knoten, nachbarn_index die
                 * Spalte */
                MXSet(index, nachbarn_index[k], -gewicht * diag);

            rechts[index] = 0.;
        }
        break;
    }
}

/*************************************************************************
   ROCKFLOW - Funktion: M#Vorkond

   Aufgabe:
   Vorkonditionierer nach dem mit vorkond gewaehleten Verfahren.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E int aufgabe: 0:Start des Vork., (vor Loesen des GLS)
   Matrix und b (bei Extraktion auch x) werden veraendert
   1:Ende des Vork.   (nach Loesen des GLS)
   Ruecktransformation, nur x wird bearbeitet, free
   (neu 3/2000)   2:Linkstransformation mit der genaeherten Inversen L:
   x <-- L*b, x darf gleich b sein
   (mehrfach innerhalb des Loesers aufrufbar)
   3:Linkstransformation mit der genaeherten Inversen L(trans)
   x <-- L(t)*b, x darf gleich b sein
   (mehrfach innerhalb des Loesers aufrufbar)

   Es ist L ~= A**(-1). Das GLS A*x =b wird transformiert zu
   L*A *x ~= I*x = L*b
   zu transformieren sind die Residuen  r = L*(b-Ax)           -> Aufgabe 2
   oder bei quadriertem Problem         r = A(t)*L(t)*L*(b-Ax) -> Aufgabe 3

   E double *b : Vektor b
   E double *x : Vektor x

   Ergebnis:  - void -

   Programmaenderungen:
   Modell 1:   11/1995     MSR        Erste Version
   3/2000     Ra         Gauss-Precond. (direkte Loesung)
   6/2000     C.Thorenz  Hilfsvariable Diagonalenskalierung
   Modell 2:   11/1995     MSR        Erste Version
   10/1999     C.Thorenz  Zweite Version
   6/2000     C.Thorenz  Hilfsvariable Diagonalenskalierung
   Modell 3,4: 2/2000     Ra         1. Version (nur Extraktion und ILU)
   6/2000     C.Thorenz  Diagonalenskalierung eingebaut (Nur Mod.4)
*** Definitionen fuer Preconditioner ************************************/
#define fastNull MKleinsteZahl
#define VK_Skalierung 1
#define VK_Extraktion 10
#define VK_iLDU 100
#define VK_Modus(vk) ((int)(vorkond % (vk * 10) / vk))
/* werden am Ende von matrix.c undefiniert
   Hilfsfunktion Skalarprodukt a*b, (Elementabstand da in a, 1 in db) */
double H_M1skprod(double* a, long da, double* b, long n)
{
    double s = 0.0;
    register long j, i = 0;
    for (j = 0; j < n; j++)
    {
        s += a[i] * b[j];
        i += da;
    }
    return s;
}

/**** Modell 1 ************************************************************/
void M1Vorkond(int aufgabe, double* x, double* b)
{
    Modell1* w = (Modell1*)wurzel;
    register long i, k;
    register double h;
    static double *x0, *r0, *Gik = NULL;

#ifdef ERROR_CONTROL
    if (b == NULL || x == NULL)
        MX_Exit("M1Vorkond", 3);
#endif

    switch (aufgabe)
    {
        case 0: /* Start des Vorkonditionierers */
            if
                VK_Modus(VK_Extraktion)
                { /* immer zuerst! */
                    x0 = (double*)Malloc(sizeof(double) * dim);
                    r0 = (double*)Malloc(sizeof(double) * dim);
                    MXResiduum(x, b,
                               r0); /* Rechte Seite ex A*(Startloesung x) */
                    for (i = 0; i < dim; i++)
                    {
                        b[i] = r0[i];
                        x0[i] = x[i];
                        x[i] = 0.0;
                    }
                    r0 = (double*)Free(r0);
                }
            if
                VK_Modus(VK_Skalierung)
                { /* Diagonal-Skalierung */
                    for (i = 0; i < dim; i++)
                    {
                        if (fabs(Aik1(i, i)) > DBL_MIN)
                        {
                            h = 1. / Aik1(i, i);
                            b[i] *= h;
                            for (k = 0; k < dim; k++)
                                Aik1(i, k) *= h;
                        }
                        else
                        {
                            DisplayMsg("!!! Equation system: Line: ");
                            DisplayLong(i);
                            DisplayMsg(" Value: ");
                            DisplayDouble(Aik1(i, i), 0, 0);
                            DisplayMsgLn(
                                "!!! Diagonal near zero! Disable diagonal "
                                "preconditioner!");
                            exit(1);
                        }
                    }
                }
            break;

        case 1: /* Ende des Vorkonditionierers */
            if
                VK_Modus(VK_iLDU) Gik = (double*)Free(Gik);

            if
                VK_Modus(VK_Extraktion)
                { /* immer zuletzt: Startloesung addieren */
                    for (i = 0; i < dim; i++)
                        x[i] += x0[i];
                    x0 = (double*)Free(x0);
                }
            break;

        case 2: /* Linkstransformation des Gesamtsystems x <= L*b */
            if
                VK_Modus(VK_iLDU)
                { /*  Gauss anstelle L(D)U-Zerlegung */
                    long ia, diag, ka = 0;
                    if (Gik == NULL)
                    { /* Matrix kopieren (transponiert!) und zerlegen */
                        long ig = 0;
                        Gik = (double*)Malloc(dim * dim * sizeof(double));
                        for (k = 0; k < dim; k++)
                        {
                            ia = k;
                            for (i = 0; i < dim; i++)
                            {
                                Gik[ig++] = w->matrix[ia];
                                ia += dim;
                            }
                        }

                        for (k = 0; k < dim; k++)
                        { /* alle Spalten/Zeilen */
                            diag = 0;
                            for (i = 0; i < k; i++)
                            { /* Spalte ohne Okk */
                                Gik[ka + i] -=
                                    H_M1skprod(&Gik[i], dim, &Gik[ka], i);
                                Gik[ka + i] *= Gik[diag];
                                diag += dim + 1; /* /Uii (Kehrwert) */
                            }                    /* i */

                            ia = 0; /* Spaltenanfang Spalte i */
                            for (i = 0; i <= k; i++)
                            { /* Zeile incl. Ukk */
                                Gik[ia + k] -=
                                    H_M1skprod(&Gik[k], dim, &Gik[ia], i);
                                ia += dim;               /* Uik */
                            }                            /* i */
                            Gik[diag] = 1.0 / Gik[diag]; /* 1.0/Ukk */
                            ka += dim;                   /* naechste Spalte k */
                        }
                    } /* k - Matrix zerlegt */
                    /* Gleichungssystem aufloesen */
                    diag = 0;
                    for (i = 0; i < dim; i++)
                    { /*  vorwaerts einsetzen mit Uik */
                        x[i] =
                            (b[i] - H_M1skprod(&Gik[i], dim, x, i)) * Gik[diag];
                        diag += dim + 1; /* naechstes Diagonalelement */
                    }
                    diag--; /* rechts neben letztem Diagonalelement */
                    for (i = dim; i > 0; i--)
                    { /* rueckwaerts einsetzen mit Oik */
                        x[i - 1] -= H_M1skprod(&Gik[diag], dim, &x[i], dim - i);
                        diag -= dim + 1;
                    } /* i */
                }     /* Modus iLDU, (zerlegen und) aufloesen */
        case 3:       /* Linkstransformation des Gesamtsystems x <= L(t)*b */
            if
                VK_Modus(VK_iLDU)
                { /*  Gauss anstelle L(D)U-Zerlegung */
                    long ia = 0, diag = dim * dim - 1;
                    /* Nur aufloesen, Matrix ist immer schon zerlegt! */
                    for (i = 0; i < dim; i++)
                    { /*  vorwaerts einsetzen mit Oik */
                        x[i] = b[i] - H_M1skprod(&Gik[ia], 1l, x, i);
                        ia += dim;
                    }

                    ia--; /* rechts neben letztem Diagonalelement */
                    for (i = dim; i > 0; i--)
                    { /* rueckwaerts einsetzen mit Uik */
                        x[i - 1] -= H_M1skprod(&Gik[ia], 1l, &x[i], dim - i);
                        x[i - 1] *= Gik[diag];
                        diag -= dim + 1;
                    } /* i */
                }     /* Modus iLDU, aufloesen mit L(t) */
    }                 /* aufgabe */
} /* M1Vorkond */

/**** Modell 2 ************************************************************/
void M2Vorkond(int aufgabe, double* x, double* b)
{
    Modell2* w = (Modell2*)wurzel;
    register long k;
    register int i;
    static double *x0, *r0, h;

#ifdef ERROR_CONTROL
    if (b == NULL || x == NULL)
        MX_Exit("M2Vorkond", 3);
#endif
    //======================================================================
    switch (aufgabe)
    {
        //--------------------------------------------------------------------
        case 0: /* Start des Vorkonditionierers */
            if
                VK_Modus(VK_Extraktion)
                { /* immer zuerst! */
                    x0 = (double*)Malloc(sizeof(double) * dim);
                    r0 = (double*)Malloc(sizeof(double) * dim);
                    MXResiduum(x, b,
                               r0); /* Rechte Seite ex A*(Startloesung x) */
                    for (i = 0; i < dim; i++)
                    {
                        b[i] = r0[i];
                        x0[i] = x[i];
                        x[i] = 0.0;
                    }
                    r0 = (double*)Free(r0);
                }
            if
                VK_Modus(VK_Skalierung)
                { /* Diagonal-Skalierung */
                    for (k = 0; k < dim; k++)
                    {
                        if (fabs(Diag2(k)) > DBL_MIN)
                        {
                            h = 1. / Diag2(k);
                            b[k] *= h;
                            for (i = 0; i < Zeil2(k).anz; i++)
                                Aik2(k, i) *= h;
                        }
                        else
                        {
                            DisplayMsg("!!! Equation system: Line: ");
                            DisplayLong(k);
                            DisplayMsg(" Value: ");
                            DisplayDouble(Diag2(k), 0, 0);
                            DisplayMsgLn("");
                            DisplayMsgLn(
                                "!!! Diagonal near zero! Disable diagonal "
                                "preconditioner!");
                            exit(1);
                        }
                    }
                }
            break;
        //--------------------------------------------------------------------
        case 1: /* Ende des Vorkonditionierers */
            if
                VK_Modus(VK_Extraktion)
                { /* immer zuletzt: Startloesung addieren */
                    for (i = 0; i < dim; i++)
                        x[i] += x0[i];
                    x0 = (double*)Free(x0);
                }
            break;
        //--------------------------------------------------------------------
        case 2:
        //--------------------------------------------------------------------
        case 3: /* Linkstransformationen */
            if
                VK_Modus(VK_iLDU) /*  incomplete L(D)U-Zerlegung geht nicht! */
                    DisplayMsgLn("Modell 2: kein ILU-Vorkonditionierer!");
            //--------------------------------------------------------------------
    }
    //======================================================================
}

/**** Modell 3,4 **********************************************************/
void M34Vorkond(int aufgabe, double* x, double* b)
{
    Modell34* w = (Modell34*)wurzel;
    long i, k, l;
    register long zk;
    register int ji, jk;
    int j, u = w->usym, o = 0;
    double Oik, Uki = 0.0, Dkk, r;
    static double *x0, *r0, h;

#ifdef ERROR_CONTROL
    if (b == NULL || x == NULL)
        MX_Exit("M34Vorkond", 3);
#endif

    switch (aufgabe)
    {
        case 0: /* Start des Vorkonditionierers - ILU! */
            if (VK_Modus(VK_Extraktion)) /* immer zuerst! */
            {
                x0 = (double*)Malloc(sizeof(double) * dim);
                r0 = (double*)Malloc(sizeof(double) * dim);
                MXResiduum(x, b, r0); /* Rechte Seite ex A*(Startloesung x) */
                for (i = 0; i < dim; i++)
                {
                    b[i] = r0[i];
                    x0[i] = x[i];
                    x[i] = 0.0;
                }
                r0 = (double*)Free(r0);
            }
            if (VK_Modus(VK_Skalierung))
            {
                for (k = 0; k < dim; k++)
                    if (fabs(w->Diag[k]) > DBL_MIN)
                    {
                        h = 1. / w->Diag[k];
                        b[k] *= h;
                        for (l = 0; l < dim; l++)
                            M34Mul(k, l, h);
                    }
                    /*
                       else {
                       DisplayMsg("!!! Equation system: Line: ");
                       DisplayLong(k);
                       DisplayMsg(" Value: ");
                       DisplayDouble(w -> Diag[k], 0, 0);
                       DisplayMsgLn("");
                       DisplayMsgLn("!!! Diagonal near zero! Disable diagonal
                       preconditioner!"); exit(1);
                       }*/

#ifdef geht_nicht
                for (k = 0; k < dim; k++)
                    for (j = 0; j < Sp34(k).anz; j++)
                    {
                        i = Ind34(k, j);
                        Aik34(k, j, u) /= w->Diag[k]; /* Untermatrix */
                        Aik34(k, j, 0) /= w->Diag[i]; /* Obermatrix  o=0! */
                    }
                for (k = 0; k < dim; k++)
                    w->Diag[k] = 1.;
#endif
            }
            break;

        case 1: /* Ende des Vorkonditionierers */
            if (VK_Modus(
                    VK_Extraktion)) /* immer zuletzt: Startloesung addieren */
            {
                for (i = 0; i < dim; i++)
                    x[i] += x0[i];
                x0 = (double*)Free(x0);
            }
            break; /* Matrix fuer ILU-Zerlegung wird nicht freigegeben! */

        case 3: /* Linkstransformation des Gesamtsystems x <= L(t)*b */
            u = 0;
            o = w->usym;
            /* kein break! "Fall through" nach Aufgabe 2 */

        case 2:                    /* Linkstransformation  x <= L*b */
            if (VK_Modus(VK_iLDU)) /* incomplete L(D)U-Zerlegung */
            {
                if (w->stat < 2) /* Matrix ist noch nicht zerlegt */
                {
                    for (k = 0; k < dim; k++) /* Spalten reduzieren */
                    {
                        for (j = 0; j < Sp34(k).anz; j++)
                        {
                            Oik = Aik34(k, j, o);
                            if (u)
                                Uki = Aik34(k, j, u);
                            i = Ind34(k, j);

                            /* Skalarprodukte abziehen */
                            jk = 0;
                            zk = Ind34(k, jk); /* oberstes Element Spalte k */
                            for (ji = 0; ji < Sp34(i).anz; ji++)
                            {
                                while (Ind34(i, ji) > zk)
                                    zk = Ind34(k, ++jk); /* zk existiert! */
                                if (zk == Ind34(i, ji))
                                {
                                    Oik -= Bik34(i, ji, u) * Bik34(k, jk, o);
                                    if (u)
                                        Uki -=
                                            Bik34(i, ji, o) * Bik34(k, jk, u);
                                }
                            } /* Ende ji Skalarprodukt */
                            Bik34(k, j, o) = Oik;
                            if (u)
                                Bik34(k, j, u) = Uki;
                        } /* j, Spaltenelemente staffeln */

                        /* Diagonale extrahieren mit Sk.prod. fuer
                         * Diagonalelement */
                        Dkk = w->Diag[k];
                        for (jk = 0; jk < Sp34(k).anz; jk++)
                        {
                            Oik = Bik34(k, jk, o);
                            zk = Ind34(k, jk);
                            Bik34(k, jk, o) *= w->PreD[zk];
                            if (u)
                                Bik34(k, jk, u) *= w->PreD[zk];
                            Dkk -= Bik34(k, jk, u) * Oik;
                        }
                        if (fabs(Dkk) < fastNull)
                            /*sign */
                            Dkk = (Dkk < 0.0) ? -fastNull : fastNull;
                        w->PreD[k] = 1.0 / Dkk; /* Kehrwert speichern */
                    }
                    w->stat = 2; /* merken! */
                }                /* Ende der Zerlegung */
                /* genaeherte Loesung von  A*x = b. Ergebnis: x */
                for (k = 0; k < dim;
                     k++) /* vorwaerts einsetzen mit Untermatrix */
                {
                    r = b[k];
                    for (j = 0; j < Sp34(k).anz; j++)
                        r -= Bik34(k, j, u) * x[Ind34(k, j)];
                    x[k] = r;
                }
                for (k = 0; k < dim; k++)
                    x[k] *= w->PreD[k]; /* Diagonal-Normierung */
                for (k = dim - 1; k > 0;
                     k--) /* rueckwaerts einsetzen mit Obermatrix */
                {
                    r = x[k];
                    for (j = 0; j < Sp34(k).anz; j++)
                        x[Ind34(k, j)] -= Bik34(k, j, o) * r;
                }
            } /* Ende ILU-Vorkonditionierer */
    }         /* switch aufgabe */
} /*M34Precond */

/**** Modell 5 ************************************************************/
// WW/PA 16/02/2006
void M5Vorkond(int aufgabe, double* x, double* b)
{
    Modell2* w = (Modell2*)wurzel;
    register long k;
    register int i;
    static double *x0, *r0, h;
    long ii;  //, count1;
    double v_diag = 0.0;

#ifdef ERROR_CONTROL
    if (b == NULL || x == NULL)
        MX_Exit("M2Vorkond", 3);
#endif
    //======================================================================
    switch (aufgabe)
    {
        //--------------------------------------------------------------------
        case 0: /* Start des Vorkonditionierers */
            if
                VK_Modus(VK_Extraktion)
                { /* immer zuerst! */
                    x0 = (double*)Malloc(sizeof(double) * dim);
                    r0 = (double*)Malloc(sizeof(double) * dim);
                    MXResiduum(x, b,
                               r0); /* Rechte Seite ex A*(Startloesung x) */
                    for (i = 0; i < dim; i++)
                    {
                        b[i] = r0[i];
                        x0[i] = x[i];
                        x[i] = 0.0;
                    }
                    r0 = (double*)Free(r0);
                }
            if
                VK_Modus(VK_Skalierung)
                {
                    for (k = 0; k < Dim_L; k++)
                    {
                        v_diag = jdiag[diag5_i[k]];
                        if (fabs(v_diag) > DBL_MIN)
                        {
                            h = 1. / v_diag;
                            b[k] *= h;
                        }
                        else
                        {
                            DisplayMsg("!!! Equation system: Line: ");
                            DisplayLong(k);
                            DisplayMsg(" Value: ");
                            DisplayDouble(Diag2(k), 0, 0);
                            DisplayMsgLn("");
                            DisplayMsgLn(
                                "!!! Diagonal near zero! Disable diagonal "
                                "preconditioner!");
                            exit(1);
                        }
                    }

                    // For test, must be improved

                    double val;
                    for (k = 0; k < Dim_L; k++)
                    {
                        v_diag = MXGet(k, k);
                        for (ii = 0; ii < Dim_L; ii++)
                        {
                            val = MXGet(k, ii);
                            val /= v_diag;
                            MXSet(k, ii, val);
                        }
                    }

                    /*
                       long count1 = 0;
                       for (k = 0; k < jd_ptr_max; k++)
                       {
                           for (i = 0; i < jd_ptr[k]; i++)
                           {
                              // Row of equation system
                              ii = jd_ptr2[i];
                              v_diag = jdiag[diag5_i[ii]];
                              if(fabs(v_diag) > DBL_MIN)
                                 jdiag[count1] /= v_diag;  //
                       count1++;
                       }
                       }
                     */
                    /////////////////////////////////////////////////

                    /*
                       // Diagonal-Skalierung
                       for (k = 0; k < dim; k++)
                       {
                        if(fabs(Diag2(k)) > DBL_MIN) {
                          h = 1. / Diag2(k);
                          b[k] *= h;
                          for (i = 0; i < Zeil2(k).anz; i++)
                            Aik2(k, i) *= h;
                        } else {
                          DisplayMsg("!!! Equation system: Line: ");
                       DisplayLong(k);
                       DisplayMsg(" Value: ");
                       DisplayDouble(Diag2(k), 0, 0);
                       DisplayMsgLn("");
                       DisplayMsgLn("!!! Diagonal near zero! Disable diagonal
                       preconditioner!"); exit(1);
                       }
                       }
                     */
                }
            break;
        //--------------------------------------------------------------------
        case 1: /* Ende des Vorkonditionierers */
            if
                VK_Modus(VK_Extraktion)
                { /* immer zuletzt: Startloesung addieren */
                    for (i = 0; i < dim; i++)
                        x[i] += x0[i];
                    x0 = (double*)Free(x0);
                }
            break;
        //--------------------------------------------------------------------
        case 2:
        //--------------------------------------------------------------------
        case 3: /* Linkstransformationen */
            if
                VK_Modus(VK_iLDU) /*  incomplete L(D)U-Zerlegung geht nicht! */
                    DisplayMsgLn("Modell 2: kein ILU-Vorkonditionierer!");
            //--------------------------------------------------------------------
    }
    //======================================================================
}

/*************************************************************************
   ROCKFLOW - Funktion: MXDumpGLS

   Aufgabe:
   Erzeugt einen Dump des Gleichungssystems

   Die Matrix muss mit MXSetMatrixPointer aktiviert worden sein.
   Der Vektor muss die zugehoerige Laenge (dim) haben, die von MXGetDim
   geliefert wird.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E char name       : Filename fuer die Ausgabedatei
   E int modus       : modus=0 : Komplettes GLS wird ausgegebeen
   : modus=1 : Eine Map der ungleich-Null Werte wird erzeugt
   : modus=2 : Komplette Matrix wird ausgegebeen
   : modus=3 : Eine Map der ungleich-Null Werte der Matrix wird erzeugt
   E double *rechts  : rechte Seite des GLS
   E double *ergebnis: Ergebnis des GLS

   Ergebnis:
   - void -

   Programmaenderungen:
   3/1998       C. Thorenz        Erste Version
   5/1998       C. Thorenz        Verallgemeinert
   09/2002   OK   besser > 0 als > MKleinsteZahl

*** benutzt keine modellspezifischen Dinge ******************************/
void MXDumpGLS(const char* name, int modus, double* rechts, double* ergebnis)
{
    register long i, j;
    long NonZeroMatrix = 0, NonZeroRHS = 0;
    FILE* dumpfile;

    dumpfile = fopen(name, "w");
    if (dumpfile == NULL)
        return;

    for (i = 0; i < dim; i++)
    {                             /* Nicht-Null-Elemente zaehlen */
        for (j = 0; j < dim; j++) /* !quadratisches Problem! (Ra) */
            if (fabs(MXGet(i, j)) > MKleinsteZahl)
                NonZeroMatrix++;
        if (fabs(rechts[i]) > MKleinsteZahl)
            NonZeroRHS++;
    }

    if (modus == 0)
    { /* Map des gesamten GLS */
        FilePrintString(dumpfile, "Variables=i,j,X\nZone, i=");
        FilePrintLong(dumpfile, dim * (dim + 2));
        FilePrintString(dumpfile, ", f=point\n");

        for (i = 0; i < dim; i++)
        {
            for (j = 0; j < dim; j++)
                FilePrintDouble(dumpfile, MXGet(i, j));
            FilePrintDouble(dumpfile, rechts[i]);
            FilePrintDouble(dumpfile, ergebnis[i]);
            FilePrintString(dumpfile, "\n");
        }
    }
    if (modus == 1)
    { /* Map der ungleich-null Werte des GLS */
        FilePrintString(dumpfile, "Variables=i,j,X\n");
        for (i = 0; i < dim; i++)
        {
            for (j = 0; j < dim; j++)
                if (fabs(MXGet(i, j)) > 0.0)
                {
                    FilePrintDouble(dumpfile, (double)i);
                    FilePrintDouble(dumpfile, (double)j);
                    FilePrintDouble(dumpfile, MXGet(i, j));
                    FilePrintString(dumpfile, "\n");
                }
            /*
               if (fabs(rechts[i]) > 0.0)
               {
                FilePrintDouble(dumpfile, (double) i);
                FilePrintDouble(dumpfile, (double) (dim + 1));
                FilePrintDouble(dumpfile, rechts[i]);
                FilePrintString(dumpfile, "\n");
               }
             */
            if (ergebnis)
                if (fabs(ergebnis[i]) > MKleinsteZahl)
                {
                    FilePrintDouble(dumpfile, (double)i);
                    FilePrintDouble(dumpfile, (double)(dim + 2));
                    FilePrintDouble(dumpfile, ergebnis[i]);
                    FilePrintString(dumpfile, "\n");
                }
        }

        // TEST  WW
        FilePrintString(dumpfile, "RHS\n");
        for (i = 0; i < dim; i++)
        {
            FilePrintDouble(dumpfile, (double)i);
            FilePrintDouble(dumpfile, rechts[i]);
            FilePrintString(dumpfile, "\n");
        }
    }
    if (modus == 2)
    {
        FilePrintString(dumpfile, "Variables=i,j,X\nZone, i=");
        FilePrintLong(dumpfile, dim * dim);
        FilePrintString(dumpfile, ", f=point\n");
        for (i = 0; i < dim; i++)
        {
            for (j = 0; j < dim; j++)
                FilePrintDouble(dumpfile, MXGet(i, j));
            FilePrintString(dumpfile, "\n");
        }
    }
    if (modus == 3)
    {
        FilePrintString(dumpfile, "Variables=i,j,X\n");
        for (i = 0; i < dim; i++)
        {
            for (j = 0; j < dim; j++)
                if (fabs(MXGet(i, j)) > MKleinsteZahl)
                {
                    FilePrintDouble(dumpfile, (double)i);
                    FilePrintDouble(dumpfile, (double)j);
                    FilePrintDouble(dumpfile, MXGet(i, j));
                    FilePrintString(dumpfile, "\n");
                }
        }
    }
    fclose(dumpfile);
}

/*************************************************************************
   ROCKFLOW - Funktion: MXEstimateStartVector

   Aufgabe:
   Versucht, eine gute Startloesung fuer das GLS zu finden
   Modell 1 und 2: mit Zeilen-Masslumping (Thorenz)
   Modell 3 und 4: ueber Preconditioner   (Ratke)

   Die Matrix muss mit MXSetMatrixPointer aktiviert worden sein und
   die Vektoren muessen die zugehoerige Laenge (dim) haben!

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E double *ReSei  : Rechte Seite des GLS
   E double *x0     : Ermittelter Startvektor

   Ergebnis:
   - void -

   Programmaenderungen:
   3/1998   C. Thorenz  Erste Version
   10/1998  C. Thorenz  Erste Version
   3/2000   R. Ratke    Benutzung von MXMatVek (modellunabhaeng,
   vermeidet quadratischen Aufwand!)
*************************************************************************/
void MXEstimateStartVector(double* ReSei, double* x0)
{
    register long i;

    if (VK_Modus(VK_iLDU) && matrix_type != 2)
        MXVorkond(2, x0, ReSei); /* Startvektor aus ILU-Zerlegung */

    else
    { /* Startvektor aus Zeilen-Lumping */
        double* sum = (double*)Malloc(dim * sizeof(double));
        for (i = 0; i < dim; i++)
            x0[i] = 1.0;
        MXMatVek(x0, sum); /* sum: Zeilensummen von A */
        for (i = 0; i < dim; i++)
            x0[i] = ReSei[i] / sum[i];
        sum = (double*)Free(sum);
    }
}

/* lokale Makros entfernen ******************************* */
#undef dim
#undef matrix_type
#undef Maxim
#undef Minim
/* Modell 1: */
#undef Matrix1
#undef Aik1
/* Modell 2: */
#undef Zeil2
#undef Aik2
#undef Ind2
#undef Diag2
/* Modelle 3,4: */
#undef Sp34
#undef Aik34
#undef Bik34
#undef Ind34
/* Makros fuer Preconditioner entfernen */
#undef fastNull
#undef VK_Skalierung
#undef VK_Extraktion
#undef VK_iLDU
#undef VK_Modus
/* Ende Matrix.c */
