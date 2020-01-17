/**
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/*************************************************************************
   ROCKFLOW - Modul: solver.c

   Aufgabe:
   Solver zum Loesen von linearen und nichtlinearen Gleichungssystemen
   (Zusammenfassung der Dateien loeserx.c)

   Programmaenderungen:

   06/1999   OK   erweiterte Lesefunktion
   10/1999   AH   Systemzeit
   11/1999   CT   Fehlerausgabe "SOLVER_SHOW_ERROR" korrigiert
   11/1999   CT   Neuer Vorkonditionierer eingehaengt
   03/2000   RR   1. Speichertechniken 3 und 4 sowie
   2. ILDU-Vorkonditionierer zugefuegt
   3. einfache Vektoroperationen direkt mit Register!
   09/2000   CT   Neues Fehlerkriterium, Funktionen zu alten Schluesselworten
entfernt 11/2000   OK   alte Lesefunktionen raus 11/2000   CT   0/0 Problem in
BICGStab 01/2001   CT   - Neues Fehlerkriterium: max(||Ax||, ||b||, ||r0||)
   - Restart in BiCGSTAB
   - Zweites Abbruchkriterium in BiCGSTAB
   03/2002   CT   Interface zu AMG1R5 und UMF-Pack
   03/2003   RK   Quellcode bereinigt, Globalvariablen entfernt

*************************************************************************/
#include <cfloat>
#include <iostream>
using namespace std;
#include "files0.h"
#include "makros.h"
#include "mathlib.h"
#include "matrix_routines.h"
#include "rf_pcs.h"  //OK_MOD"
#include "rf_tim_new.h"
#include "solver.h"
#include "tools.h"
#include "mathlib.h"
#include "memory.h"
#include "display.h"

/* AMG-Solver */
#ifdef AMG1R5
#include "amg1r5.h"
#endif

/* UMFPACK-Solver */
#ifdef UMFPACK
#include <umfpack.h>
#endif

using namespace Display;

/**** Definitionen fuer Preconditioner (Ra, 3/2000) */
#define VK_Skalierung 1
#define VK_Extraktion 10
#define VK_iLDU 100
#define VK_Modus(vk) ((int)(vorkond % (vk * 10) / vk))
/* werden am Ende dieser Quelle wieder undefiniert! */

#define noTESTLOES
#define noTESTLOES1
#define noTESTLOES4
#define noTESTLOES6
#define noSOLVER_SHOW_ERROR

/* Interne (statische) Deklarationen
   ToDo -> InitSolverParameter */
/* #LINEAR_SOLVER_PRESSURE, #LINEAR_SOLVER_SATURATION, #LINEAR_SOLVER_TRANSPORT
 */
int loeser_flow = 2, loeser_tran = 2, loeser_temp = 2, loeser_satu = 2;
double gls_iter_theta;
int vorkond = 1, vorkond_flow = 100, vorkond_tran = 100, vorkond_temp = 100,
    vorkond_satu = 100;
int speichertechnik_flow = 4, speichertechnik_tran = 4,
    speichertechnik_temp = 4, speichertechnik_satu = 4;
int linear_error_type = 6, linear_error_type_flow = 6,
    linear_error_type_tran = 6, linear_error_type_temp = 6,
    linear_error_type_satu = 6;
double eps_flow = 1.e-10, eps_tran = 1.e-10, eps_temp = 1.e-10,
       eps_satu = 1.e-10, cg_eps = 1.e-10;
int maxiter_flow = -1, maxiter_tran = -1, maxiter_temp = -1, maxiter_satu = -1,
    cg_maxiter = -1;
int repeat_flow = 1, repeat_tran = 1, repeat_temp = 1, repeat_satu = 1,
    cg_repeat = 1;

int sp2_start, sp2_inc;

double rel_eps;

/* #ITERATION_FLOW, #ITERATION_TRANSPORT, #ITERATION_CONTROL */
int nonlinear_method;
int nonlinear_maxiter;
int nonlinear_convergence_type;
double nonlinear_abs_eps;
double nonlinear_rel_eps;
double nonlinear_rel_cg_eps;
int nonlinear_assemble;

int nonlinear_method_flow;
int nonlinear_maxiter_flow;
int nonlinear_convergence_type_flow;
double nonlinear_abs_eps_flow;
double nonlinear_rel_eps_flow;
double nonlinear_rel_cg_eps_flow;
int nonlinear_assemble_flow;

int nonlinear_method_tran;
int nonlinear_maxiter_tran;
int nonlinear_convergence_type_tran;
double nonlinear_abs_eps_tran;
double nonlinear_rel_eps_tran;
double nonlinear_rel_cg_eps_tran;
int nonlinear_assemble_tran;

int iteration_min_iter;
int iteration_max_iter;
double iteration_weight_plus;
double iteration_weight_up;
double iteration_weight_down;
double iteration_min_dt;
double iteration_max_dt;

IntFuncDXDXL LinearSolver;
IntFuncDXDXLVXL NonlinearSolver;
#ifndef USE_MPI
IntFuncDXDXL LoeserFlow = SpBICGSTAB, LoeserTran = SpBICGSTAB,
             LoeserTemp = SpBICGSTAB, LoeserSatu = SpBICGSTAB;
#endif

/* LU Dekomposition */
void ludcmp_3(double* a, long n, long* indx);
void lubksb_3(double* a, long n, long* indx, double* b);

/*************************************************************************
   ROCKFLOW - Funktion: InitSolverParameter

   Aufgabe:
   Initialisierung von Solver-Parametern

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   void

   Ergebnis:
   - void -

   Programmaenderungen:
   07/1999   OK   Implementierung
   11/2000   CT   Repeat eingefuegt, Genauigkeiten hoeher

*************************************************************************/
void InitSolverParameter(void)
{
    /* #LINEAR_SOLVER_PRESSURE */
    loeser_flow = 2;
    maxiter_flow = -1;
    eps_flow = 1.e-10;
    repeat_flow = -1;
    vorkond_flow = VK_Extraktion + VK_iLDU; /* Ra, 3/2000 */
    speichertechnik_flow = 4;               /* Ra, 3/2000 */
    linear_error_type_flow = 6;

    /* #LINEAR_SOLVER_SATURATION */
    loeser_satu = 2;
    maxiter_satu = -1;
    eps_satu = 1.e-10;
    repeat_satu = -1;
    vorkond_tran = VK_Extraktion + VK_iLDU; /* Ra, 3/2000 */
    speichertechnik_tran = 4;
    linear_error_type_satu = 6;

    /* #LINEAR_SOLVER_TRANSPORT */
    loeser_tran = 2;
    maxiter_tran = -1;
    eps_tran = 1.e-10;
    repeat_tran = -1;
    vorkond_tran = VK_Extraktion + VK_iLDU; /* Ra, 3/2000 */
    speichertechnik_tran = 4;
    linear_error_type_tran = 6;

    /* #LINEAR_SOLVER_TEMPERATURE */
    loeser_temp = 2;
    maxiter_temp = -1;
    eps_temp = 1.e-10;
    repeat_temp = -1;
    vorkond_temp = VK_Extraktion + VK_iLDU; /* Ra, 3/2000 */
    speichertechnik_temp = 4;
    linear_error_type_temp = 6;

    rel_eps = 1.e-4;

    /* #ITERATION_FLOW */
    nonlinear_method_flow = 1;
    nonlinear_maxiter_flow = 1000;
    nonlinear_convergence_type_flow = 1;
    nonlinear_abs_eps_flow = 1.;
    nonlinear_rel_eps_flow = 1.e-9;
    nonlinear_rel_cg_eps_flow = 0.;
    nonlinear_assemble_flow = 1;

    /* #ITERATION_TRANSPORT */
    nonlinear_method_tran = 1;
    nonlinear_maxiter_tran = 1000;
    nonlinear_convergence_type_tran = 1;
    nonlinear_abs_eps_tran = 1.;
    nonlinear_rel_eps_tran = 1.e-9;
    nonlinear_rel_cg_eps_tran = 0.;
    nonlinear_assemble_tran = 1;

    /* #ITERATION_TIME_CONTROL */
    iteration_min_iter = 5;
    iteration_max_iter = 100;
    iteration_weight_plus = 10.;
    iteration_weight_up = 1.4;
    iteration_weight_down = 0.7;
    iteration_min_dt = 1.e-3;
    iteration_max_dt = 1.e+6;
}

/*************************************************************************
   ROCKFLOW - Modul: loeser1.c

   Aufgabe:
   Iterative-Gleichungsloeser mit Speichertechnik aus 'matrix.h'
   Herk¸mmliche Verfahren wie Richardson, Jacobi, Gauss-Seidel und SOR.
   Konvergieren zwar sehr langsam im Vergleich zu den neuen Verfahren
   wie PCG oder BICGSTAB, dennoch haben sie gute Eigenschaften, die spaeter
   benoetigt werden k¸nnen (Stichwort: Multigrid, Decomposition domain,
   Vorkonditionierung, etc.)

   Programmaenderungen:
   02/1998     A. Habbar  Richardson-Loeser, JOR-Loeser, SOR-Loeser (Erste
Version)

*************************************************************************/

/*************************************************************************
   ROCKFLOW - Funktion: Richardson

   Aufgabe:
   Gleichungsloeser (Richardson-Iteration)
   Algorithmus laueft ueber den Defekt, Implementierung ueber das Residuum.
   Die Matrix wird durch die Speichertechnik implizit bereitgestellt.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E double *b: Rechte Seite des Gleichungssystems
   X double *x: Ergebnisvektor, Speicher muss bereits reserviert sein
   E long n   : Dimension des Gleichungssystems

   Ergebnis:
   Anzahl der Iterationen

   Programmaenderungen:
   02/1998     AH     Erste Version
   9/2000      C.Thorenz  Neues Fehlerkriterium, Funktionen zu alten
   Schluesselworten entfernt
   01/2001     CT   Neues Fehlerkriterium: max(||Ax||, ||b||, ||r0||)

*************************************************************************/
int SpRichardson(double* b, double* x, long n)
{
    double *r, *s;
    double eps = cg_eps;
    int k = 0;  // WW, max_iter = 0;
    double r0norm = 1., b0norm = 1., x0norm = 1.;

    if (linear_error_type == 5)
    {
        b0norm = VEKNORM_BICGSTAB(b, n);
        x0norm = VEKNORM_BICGSTAB(x, n);
    }
    /* Ggf. starten der Vorkonditionierung */
    if (vorkond)
        MXVorkond(0, x, b);

        // WW if (cg_maxiter > 0)
        // WW   max_iter = cg_maxiter;

        // OK411    if (cg_maxiter == -1)
        // OK411        max_iter = NodeListLength;

#ifdef TESTLOES1
    DisplayMsgLn("SpRichard");
#endif

    r = (double*)Malloc(n * sizeof(double)); /* MNulleVec(r,n); */
    s = (double*)Malloc(n * sizeof(double)); /* MNulleVec(r,n); */

    MXResiduum(x, b, r);
    if (linear_error_type == 0)
        eps = cg_eps;
    else if (linear_error_type == 1)
        eps = cg_eps * VEKNORM_BICGSTAB(b, n);
    else if (linear_error_type == 2)
        eps = cg_eps * VEKNORM_BICGSTAB(r, n);
    else if (linear_error_type == 3)
    {
        if ((r0norm = VEKNORM_BICGSTAB(r, n)) > 1.)
            eps = cg_eps * r0norm;
        else
            eps = cg_eps;
    }
    else if (linear_error_type == 4)
        eps = cg_eps * (VEKNORM_BICGSTAB(x, n));
    else if (linear_error_type == 5)
    {
        b0norm = VEKNORM_BICGSTAB(b, n);
        x0norm = VEKNORM_BICGSTAB(x, n);
        r0norm = VEKNORM_BICGSTAB(r, n);
        eps = cg_eps * MMax(MMax(x0norm, b0norm), r0norm);
    }
    else if (linear_error_type == 6)
    {
        MXMatVek(x, s);
        b0norm = VEKNORM_BICGSTAB(b, n);
        x0norm = VEKNORM_BICGSTAB(s, n);
        r0norm = VEKNORM_BICGSTAB(r, n);
        eps = cg_eps * MMax(MMax(x0norm, b0norm), r0norm);
    }
    s = (double*)Free(s);

    if (VEKNORM_BICGSTAB(r, n) <= eps)
    {
        r = (double*)Free(r);
        /* Ggf. abschliessen der Vorkonditionierung */
        if (vorkond)
            MXVorkond(1, x, b);
        return 0;
    }
    for (;;)
    {
        if
            VK_Modus(VK_iLDU) MXVorkond(2, r, r); /* Ra, 3/2000 */
        MVekSum(x, gls_iter_theta, r, n);
#ifdef SOLVER_SHOW_RESULTS
        printf("\n%ld %f %f %f %f %f", (long)k, x[(long)(n * .1)],
               x[(long)(n * .3)], x[(long)(n * .5)], x[(long)(n * .7)],
               x[(long)(n * .9)]);
#endif

        MXResiduum(x, b, r);
        k++;
        if (linear_error_type == 4)
            eps = cg_eps * (VEKNORM_BICGSTAB(x, n));

        if (VEKNORM_BICGSTAB(r, n) <= eps)
            break;

        if (k >= cg_maxiter)
            break;
#ifdef SOLVER_SHOW_ERROR
        if (k % (int)MMax((cg_maxiter / 10), 1.) == 1)
        {
            DisplayMsg("      Iteration-Nr.: ");
            DisplayLong((long)k);
            DisplayMsg(", Fehler = ");
            DisplayDouble(VEKNORM_BICGSTAB(r, n) / eps, 4, 1);
            DisplayMsgLn(" ");
        }
#endif
    }

    r = (double*)Free(r);
    /* Ggf. abschliessen der Vorkonditionierung */
    if (vorkond)
        MXVorkond(1, x, b);
    return k;
}

/*************************************************************************
   ROCKFLOW - Funktion: JOR

   Aufgabe:
   Gleichungsloeser (Jacobi over relaxation, gedaempfte Jacobi-Verfahren fuer
   Faktor gls_iter_theta zwischen 0 und 1, Jacobi oder Gesamtschrittverfahren
   fuer gls_iter_theta gleich 1).
   Sehr gut geeignet fuer Mehrgitterverfahren. Laesst sich auch gut
   parallelisieren. Konvergenz ist gesichert nur in dem Fall wenn die
   Koeffizientenmatrix diagonaldominant ist.
   Algorithmus laueft ueber den Defekt, Implementierung ueber das Residuum.
   Die Matrix wird durch die Speichertechnik implizit bereitgestellt.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E double *b: Rechte Seite des Gleichungssystems
   X double *x: Ergebnisvektor, Speicher muss bereits reserviert sein
   E long n   : Dimension des Gleichungssystems

   Ergebnis:
   Anzahl der Iterationen

   Programmaenderungen:
   02/1998   AH     Erste Version
   9/2000    C.Thorenz  Neues Fehlerkriterium
   01/2001   CT   Neues Fehlerkriterium: max(||Ax||, ||b||, ||r0||)

*************************************************************************/
int SpJOR(double* b, double* x, long n)
{
    double *r, *s;
    static double eps;
    int k = 0, max_iter = 0;
    static double r0norm, b0norm, x0norm;
    register long i;
    static double h;

    /* Ggf. starten der Vorkonditionierung */
    if (vorkond)
        MXVorkond(0, x, b);

    if (cg_maxiter > 0)
        max_iter = cg_maxiter;
        // OK411    if (cg_maxiter == -1)
        // OK411        max_iter = NodeListLength;

#ifdef TESTLOES1
    DisplayMsgLn("SpJacobi");
#endif

    r = (double*)Malloc(n * sizeof(double)); /* MNulleVec(r,n); */
    s = (double*)Malloc(n * sizeof(double)); /* MNulleVec(r,n); */

    MXResiduum(x, b, r);
    if (linear_error_type == 0)
        eps = cg_eps;
    else if (linear_error_type == 1)
        eps = cg_eps * VEKNORM_BICGSTAB(b, n);
    else if (linear_error_type == 2)
        eps = cg_eps * VEKNORM_BICGSTAB(r, n);
    else if (linear_error_type == 3)
    {
        if ((r0norm = VEKNORM_BICGSTAB(r, n)) > 1.)
            eps = cg_eps * r0norm;
        else
            eps = cg_eps;
    }
    else if (linear_error_type == 4)
        eps = cg_eps * (VEKNORM_BICGSTAB(x, n));
    else if (linear_error_type == 5)
    {
        b0norm = VEKNORM_BICGSTAB(b, n);
        x0norm = VEKNORM_BICGSTAB(x, n);
        r0norm = VEKNORM_BICGSTAB(r, n);
        eps = cg_eps * MMax(MMax(x0norm, b0norm), r0norm);
    }
    else if (linear_error_type == 6)
    {
        MXMatVek(x, s);
        b0norm = VEKNORM_BICGSTAB(b, n);
        x0norm = VEKNORM_BICGSTAB(s, n);
        r0norm = VEKNORM_BICGSTAB(r, n);
        eps = cg_eps * MMax(MMax(x0norm, b0norm), r0norm);
    }
    s = (double*)Free(s);

    if (VEKNORM_BICGSTAB(r, n) <= eps)
    {
        r = (double*)Free(r);
        /* Ggf. abschliessen der Vorkonditionierung */
        if (vorkond)
            MXVorkond(1, x, b);
        return 0;
    }
    for (;;)
    {
        for (i = 0; i <= n - 1; i++)
        {
            // WW           if ((h = MXGet(i, i)) != 0.0)
            h = MXGet(i, i);

            if (fabs(h) > MKleinsteZahl)
                r[i] = r[i] / h;
            else
            {
/* Eventuell Zeilenvertauschen */
#ifdef TESTLOES1
                DisplayMsg("Fehler im Jacobi-Loser: Diagonalelement = 0.0 !!");
                DisplayMsgLn("");
#endif
                k = -1;
                break;
            }
        }
        if (k == -1)
            break;

        if
            VK_Modus(VK_iLDU) MXVorkond(2, r, r); /* Ra, 3/2000 */
        MVekSum(x, gls_iter_theta, r, n);
#ifdef SOLVER_SHOW_RESULTS
        printf("\n%ld %f %f %f %f %f", (long)k, x[(long)(n * .1)],
               x[(long)(n * .3)], x[(long)(n * .5)], x[(long)(n * .7)],
               x[(long)(n * .9)]);
#endif
        MXResiduum(x, b, r);
        k++;
        if (linear_error_type == 4)
            eps = cg_eps * (VEKNORM_BICGSTAB(x, n));

        if (VEKNORM_BICGSTAB(r, n) <= eps)
            break;

        if (k >= max_iter)
            break;

#ifdef SOLVER_SHOW_ERROR
        if (k % (int)MMax((cg_maxiter / 10), 1.) == 1)
        {
            DisplayMsg("      Iteration-Nr.: ");
            DisplayLong((long)k);
            DisplayMsg(", Fehler = ");
            DisplayDouble(VEKNORM_BICGSTAB(r, n) / eps, 4, 1);
            DisplayMsgLn(" ");
        }
#endif
    }

    r = (double*)Free(r);
    /* Ggf. abschliessen der Vorkonditionierung */
    if (vorkond)
        MXVorkond(1, x, b);

    return k;
}

/*************************************************************************
   ROCKFLOW - Funktion: SOR

   Aufgabe:
   Gleichungsloeser (SOR: successive over relaxation, gedaempfte Gauss-Seidel
   fuer Faktor gls_iter_theta zwischen 0 und 1, Gauss-Seidel oder
   Einzelschrittverfahren fuer gls_iter_theta gleich 1).

   Sehr gut geeignet fuer Mehrgitterverfahren (Interessant die Wahl des
Faktors). Laesst sich aber nicht (schlecht) parallelisieren.  Konvergenz ist
gesichert nur in dem Fall wenn die Koeffizientenmatrix diagonaldominant ist,
oder im Falle von symmetrischen und positiv definiten Matrizen fuer theta-Faktor
   zwischen 0 und 2.
   Algorithmus laueft ueber den Defekt, Implementierung ueber das Residuum.
   Die Matrix wird durch die Speichertechnik implizit bereitgestellt.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E double *b: Rechte Seite des Gleichungssystems
   X double *x: Ergebnisvektor, Speicher muss bereits reserviert sein
   E long n   : Dimension des Gleichungssystems

   Ergebnis:
   Anzahl der Iterationen

   Programmaenderungen:
   02/1998   AH     Erste Version
   9/2000    C.Thorenz  Neues Fehlerkriterium
   01/2001   CT   Neues Fehlerkriterium: max(||Ax||, ||b||, ||r0||)

*************************************************************************/
int SpSOR(double* b, double* x, long n)
{
    double *r, *s;
    static double eps;
    int k = 0, max_iter = 0;
    static double r0norm, b0norm, x0norm;
    register long i, j;
    static double h, sum;

    /* Ggf. starten der Vorkonditionierung */
    if (vorkond)
        MXVorkond(0, x, b);

    if (cg_maxiter > 0)
        max_iter = cg_maxiter;
        // OK411    if (cg_maxiter == -1)
        // OK411        max_iter = NodeListLength;

#ifdef TESTLOES1
    DisplayMsgLn("SpGaussSeidel");
#endif

    r = (double*)Malloc(n * sizeof(double)); /* MNulleVec(r,n); */
    s = (double*)Malloc(n * sizeof(double)); /* MNulleVec(r,n); */

    MXResiduum(x, b, r);
    if (linear_error_type == 0)
        eps = cg_eps;
    else if (linear_error_type == 1)
        eps = cg_eps * VEKNORM_BICGSTAB(b, n);
    else if (linear_error_type == 2)
        eps = cg_eps * VEKNORM_BICGSTAB(r, n);
    else if (linear_error_type == 3)
    {
        if ((r0norm = VEKNORM_BICGSTAB(r, n)) > 1.)
            eps = cg_eps * r0norm;
        else
            eps = cg_eps;
    }
    else if (linear_error_type == 4)
        eps = cg_eps * (VEKNORM_BICGSTAB(x, n));
    else if (linear_error_type == 5)
    {
        b0norm = VEKNORM_BICGSTAB(b, n);
        x0norm = VEKNORM_BICGSTAB(x, n);
        r0norm = VEKNORM_BICGSTAB(r, n);
        eps = cg_eps * MMax(MMax(x0norm, b0norm), r0norm);
    }
    else if (linear_error_type == 6)
    {
        MXMatVek(x, s);
        b0norm = VEKNORM_BICGSTAB(b, n);
        x0norm = VEKNORM_BICGSTAB(s, n);
        r0norm = VEKNORM_BICGSTAB(r, n);
        eps = cg_eps * MMax(MMax(x0norm, b0norm), r0norm);
    }

    s = (double*)Free(s);

    if (VEKNORM_BICGSTAB(r, n) <= eps)
    {
        r = (double*)Free(r);
        /* Ggf. abschliessen der Vorkonditionierung */
        if (vorkond)
            MXVorkond(1, x, b);
        return 0;
    }
    for (;;)
    {
        for (i = 0; i <= n - 1; i++)
        {
            sum = 0.0;
            for (j = 0; j <= n - 1; j++)
                sum = sum + MXGet(i, j) * x[j];
            /* nicht mit MXResiduum machbar, weil Einzelschrittverfahren! (Ra)
               ?? Preconditioner ?? */
            h = MXGet(i, i);

            if (fabs(h) > MKleinsteZahl)
                x[i] = x[i] - gls_iter_theta * (sum - b[i]) / h;
            else
            {
/* Eventuell Zeilenvertauschen */
#ifdef TESTLOES1
                DisplayMsg(
                    "Fehler im Gauss-Seidel-Loeser: Diagonalelement = 0.0 !!");
                DisplayMsgLn("");
#endif
                k = -1;
                break;
            }
        }
#ifdef SOLVER_SHOW_RESULTS
        printf("\n%ld %f %f %f %f %f", (long)k, x[(long)(n * .1)],
               x[(long)(n * .3)], x[(long)(n * .5)], x[(long)(n * .7)],
               x[(long)(n * .9)]);
#endif
        if (k == -1)
            break;
        MXResiduum(x, b, r);
        k++;
        if (linear_error_type == 4)
            eps = cg_eps * (VEKNORM_BICGSTAB(x, n));

        if (VEKNORM_BICGSTAB(r, n) <= eps)
            break;

        if (k >= max_iter)
            break;

#ifdef SOLVER_SHOW_ERROR
        if (k % (int)MMax((cg_maxiter / 10), 1.) == 1)
        {
            DisplayMsg("      Iteration-Nr.: ");
            DisplayLong((long)k);
            DisplayMsg(", Fehler = ");
            DisplayDouble(VEKNORM_BICGSTAB(r, n) / eps, 4, 1);
            DisplayMsgLn(" ");
        }
#endif
    }

    r = (double*)Free(r);
    /* Ggf. abschliessen der Vorkonditionierung */
    if (vorkond)
        MXVorkond(1, x, b);

    return k;
}

/*************************************************************************
   ROCKFLOW - Funktion: Gauss

   Aufgabe:
   Gleichungsloeser:
   LR-Faktorisierung von matrix nach dem Gausschen Algorithmus mit
   partieller Pivotisierung; Berechnung der Loesung aus LR-Faktorisierung
   Der Ergebnisvektor wird erst auf die rechte Seite gespeichert, und
   hinterher umkopiert. Hier kann noch Speicherplatz gespart werden.
   Die alten Inhalte von matrix und vecb werden zerstoert!
   Quelle: Schwetlick/Kretzschmar, Numerische Verfahren 1991

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E double *matrix: Linke Seite Gleichungssystem
   E double *vecb: Rechte Seite Gleichungssystem
   X double *vecx: Ergebnisvektor, Speicher muss bereits reserviert sein
   E int g: Dimension des Gleichungssystems

   Ergebnis:
   - void -

   Programmaenderungen:
   12/1994     MSR        Erste Version

*************************************************************************/
void Gauss(double* matrix, double* vecb, double* vecx, int g)
{
    /* Matrizen sind in C zeilenweise abgespeichert:
       matrix[i][j] -> matrix[i*g+j] */
    static int* s;
    static double z, hilf;
    register int k, i, j, sk;
    s = (int*)Malloc(sizeof(int) * (g - 1));
    /* LR-Faktorisierung */
    for (k = 0; k < (g - 1); k++)
    {
        /* Pivotsuche */
        z = 0.0;
        sk = 0;
        for (i = k; i < g; i++)
        {
            hilf = fabs(matrix[i * g + k]); /* matrix[i][k] */
            if (hilf > z)
            {
                z = hilf;
                sk = i;
            }
        }
        s[k] = sk;
        /* evtl. Zeilen vertauschen */
        if (sk > k)
            for (j = 0; j < g; j++)
            {
                z = matrix[k * g + j]; /* matrix[k][j] */
                /* matrix[k][j], matrix[sk][j] */
                matrix[k * g + j] = matrix[sk * g + j];
                matrix[sk * g + j] = z; /* matrix[sk][j] */
            }
        /* Berechnung der Eliminationskoeffizienten */
        for (i = (k + 1); i < g; i++)
            matrix[i * g + k] /=
                matrix[k * g + k]; /* matrix[i][k], matrix[k][k] */
        /* Spaltenweise Berechnung der neuen Restmatrix */
        for (j = (k + 1); j < g; j++)
            for (i = (k + 1); i < g; i++)
                matrix[i * g + j] -= (matrix[i * g + k] * matrix[k * g + j]);
        /* matrix[i][j], matrix[i][k], matrix[k][j] */
    }
    /* Loesung berechnen */
    /* vecb transformieren */
    for (k = 0; k < (g - 1); k++)
    {
        sk = s[k];
        if (sk > k)
        {
            z = vecb[k];
            vecb[k] = vecb[sk];
            vecb[sk] = z;
        }
    }
    for (k = 1; k < g; k++)
        for (j = 0; j < k; j++)
            vecb[k] -= (matrix[k * g + j] * vecb[j]);
    /* matrix[k][j] */
    /* vecx berechnen */
    for (k = (g - 1); k >= 0; k--)
    {
        for (j = (k + 1); j < g; j++)
            vecb[k] -= (matrix[k * g + j] * vecb[j]); /* matrix[k][j] */
        vecb[k] /= matrix[k * g + k];                 /* matrix[k][k] */
    }
    /* Umspeichern des Ergebnisses von vecb nach vecx */
    for (k = 0; k < g; k++)
        vecx[k] = vecb[k];
    /* Speicher freigeben */
    s = (int*)Free(s);
}

/*************************************************************************
   ROCKFLOW - Funktion: SpAMG1R5

   Aufgabe:
   Interface zu "Algebraischem Multigrid-Loeser" (AMG1R5)

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E double *b : Vektor b
   R double *x : Ergebnisvektor x;
   der Speicher muss bereits allokiert sein
   E long g : Dimension der Vektoren

   Ergebnis:
   Keins

   Programmaenderungen:
   10/2001     C.Thorenz  Erste Version

*************************************************************************/
int SpAMG1R5(double* b, double* x, long n)
{
#ifdef AMG1R5
    /* Variablen */
    static double *r, *t;
    long max_iter;
    int k = 0;
    static double r0norm, b0norm, x0norm;

    double *A, *U, *F;
    int *IA, *JA, *IG;
    int NDA, NDIA, NDJA, NDU, NDF, NDIG, NNU, NNA;
    int MATRIX, IFIRST, ISWTCH, IOUT, IPRINT, IERR, MADAPT;
    double EPS;
    int LEVELX, NCYC, NRD, NSOLCO, NRU, NWT, NTR;
    double ECG1, ECG2, EWT2;

    static double a_size_mult = 1.;
    static double u_size_mult = 1.;
    static double f_size_mult = 1.;
    static double ia_size_mult = 1.;
    static double ig_size_mult = 1.;
    static double ja_size_mult = 1.;

    /* Ermitteln des Residuums */
    r = (double*)Malloc(n * sizeof(double));
    t = (double*)Malloc(n * sizeof(double));

    /* Ggf. starten der Vorkonditionierung */
    if (vorkond)
        MXVorkond(0, x, b);

    /* Fehler anpassen */
    MXResiduum(x, b, r);

    switch (linear_error_type)
    {
        default:
        case 0:
            EPS = cg_eps;
            break;
        case 1:
            EPS = cg_eps * VEKNORM_BICGSTAB(b, n);
            break;
        case 2:
            EPS = cg_eps * VEKNORM_BICGSTAB(r, n);
            break;
        case 3:
            if ((r0norm = VEKNORM_BICGSTAB(r, n)) > 1.)
                EPS = cg_eps * r0norm;
            else
                EPS = cg_eps;
            break;
        case 4:
            EPS = cg_eps * (VEKNORM_BICGSTAB(x, n));
            break;
        case 5:
            b0norm = VEKNORM_BICGSTAB(b, n);
            x0norm = VEKNORM_BICGSTAB(x, n);
            r0norm = VEKNORM_BICGSTAB(r, n);
            EPS = cg_eps * MMax(MMax(x0norm, b0norm), r0norm);
            break;
        case 6:
            MXMatVek(x, t);
            b0norm = VEKNORM_BICGSTAB(b, n);
            x0norm = VEKNORM_BICGSTAB(t, n);
            r0norm = VEKNORM_BICGSTAB(r, n);
            EPS = cg_eps * MMax(MMax(x0norm, b0norm), r0norm);
            break;
    }

    switch (cg_maxiter)
    {
        case -1:
            max_iter = NodeListLength;
            break;
        default:
            max_iter = cg_maxiter;
            break;
    }

    MATRIX = 22;
    IFIRST = 10;
    ISWTCH = 4;
    IOUT = 12;
    IPRINT = 10606;
    LEVELX = 0;
    MADAPT = 0;
    NCYC = 102 * (long)pow(10., (double)(1 + (int)log10((double)max_iter))) +
           max_iter;
    NRD = 0;
    NSOLCO = 10;
    NRU = 0;
    NWT = 2;
    NTR = 0;
    ECG1 = 0.;
    ECG2 = 0.25;
    EWT2 = 0.35;

restart:

    NNU = n;
    NNA = 20 * n * a_size_mult;

    NDA = (5 * NNU + 10 * NNA) * a_size_mult;
    NDU = 5 * NNU * u_size_mult;
    NDF = 5 * NNU * f_size_mult;
    NDIA = 5 * NNU * ia_size_mult;
    NDJA = (5 * NNU + 10 * NNA) * ja_size_mult;
    NDIG = 10 * NNU * ig_size_mult;

    A = (double*)Malloc(NDA * sizeof(double));
    U = (double*)Malloc(NDU * sizeof(double));
    F = (double*)Malloc(NDF * sizeof(double));

    IA = (int*)Malloc(NDIA * sizeof(int));
    JA = (int*)Malloc(NDJA * sizeof(int));
    IG = (int*)Malloc(NDIG * sizeof(int));

    MXCopyToAMG1R5Structure(A, IA, JA, NDA, NDIA, NDJA, x, U, b, F);

    amg1r5_(A, IA, JA, U, F, IG, &NDA, &NDIA, &NDJA, &NDU, &NDF, &NDIG, &NNU,
            &MATRIX, &ISWTCH, &IOUT, &IPRINT, &LEVELX, &IFIRST, &NCYC, &EPS,
            &MADAPT, &NRD, &NSOLCO, &NRU, &ECG1, &ECG2, &EWT2, &NWT, &NTR,
            &IERR);

    /*
        aux1r5_(A,IA,JA,U,F,IG,&NDA,&NDIA,&NDJA,&NDU,&NDF,&NDIG,&NNU,&MATRIX,
               &EPS,&IFIRST,&ISWTCH,&IOUT,&IPRINT,&IERR);
     */

    MKopierVec(U, x, n);

    A = (double*)Free(A);
    U = (double*)Free(U);
    F = (double*)Free(F);

    IA = (int*)Free(IA);
    JA = (int*)Free(JA);
    IG = (int*)Free(IG);

    t = (double*)Free(t);
    r = (double*)Free(r);

    switch (IERR)
    {
        case 1:
            a_size_mult *= 1.2;
            goto restart;
        case 2:
            ia_size_mult *= 1.2;
            goto restart;
        case 3:
            ja_size_mult *= 1.2;
            goto restart;
        case 4:
            u_size_mult *= 1.2;
            goto restart;
        case 5:
            f_size_mult *= 1.2;
            goto restart;
        case 6:
            ig_size_mult *= 1.2;
            goto restart;
        case -11:
            printf(" \n AMG-Solver: A-entry missing!\n ");
            break;
        case -12:
            printf(" \n AMG-Solver: Parameters erroneous!\n ");
            break;
        case 13:
            printf(" \n AMG-Solver: Diagonal not stored first!\n ");
            break;
        case 14:
            printf(" \n AMG-Solver: Diagonal not positiv!\n ");
            break;
        case 15:
            printf(" \n AMG-Solver: IA erroneous!\n ");
            break;
        case 16:
            printf(" \n AMG-Solver: JA erroneous!\n ");
            break;
        case 17:
            printf(" \n AMG-Solver: ISWITCH erroneous!\n ");
            break;
        case 31:
            printf(" \n AMG-Solver: CG not defined!\n ");
            break;
        case -32:
            printf(" \n AMG-Solver: YALE-SMP not possible!\n ");
            break;
    }

    /* Ggf. abschliessen der Vorkonditionierung */
    if (vorkond)
        MXVorkond(1, x, b);

    return k;
#else  // ifdef AMG1R5
    n = n;
    *x = *x;
    *b = *b;
    DisplayMsg("!!!! Error: AMG1R5-Solver not included in this version!");
    DisplayMsgLn(" ");
    exit(1);
#endif
}

/*************************************************************************
   ROCKFLOW - Funktion: SpUMF

   Aufgabe:
   Interface zu UMFPACK

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E double *b : Vektor b
   R double *x : Ergebnisvektor x;
   der Speicher muss bereits allokiert sein
   E long g : Dimension der Vektoren

   Ergebnis:
   Keins

   Programmaenderungen:
   10/2001     C.Thorenz  Erste Version fuer UMFPack 3.1
   4/2002     C.Thorenz  Anpassung an UMFPack 4.0

*************************************************************************/
int SpUMF(double* b, double* x, long n)
{
#ifdef UMFPACK
    /* Variablen */
    static double *r, *t;
    long i, j, k = 0;
    double r0norm, b0norm, x0norm;
    double* Ax = NULL;
    int *Ap = NULL, *Ai = NULL;
    int NA, NNA, NP, NNP, status;
    double EPS, a;

    double *Control = (double*)NULL, *Info = (double*)NULL;
    void *Symbolic = NULL, *Numeric = NULL;

    static double a_size_mult = 1.;

    /* Ggf. starten der Vorkonditionierung */
    if (vorkond)
        MXVorkond(0, x, b);

    r = (double*)Malloc(n * sizeof(double));
    t = (double*)Malloc(n * sizeof(double));

    /* Ermitteln des Residuums */
    MXResiduum(x, b, r);

    switch (linear_error_type)
    {
        default:
        case 0:
            EPS = cg_eps;
            break;
        case 1:
            EPS = cg_eps * VEKNORM_BICGSTAB(b, n);
            break;
        case 2:
            EPS = cg_eps * VEKNORM_BICGSTAB(r, n);
            break;
        case 3:
            if ((r0norm = VEKNORM_BICGSTAB(r, n)) > 1.)
                EPS = cg_eps * r0norm;
            else
                EPS = cg_eps;
            break;
        case 4:
            EPS = cg_eps * (VEKNORM_BICGSTAB(x, n));
            break;
        case 5:
            b0norm = VEKNORM_BICGSTAB(b, n);
            x0norm = VEKNORM_BICGSTAB(x, n);
            r0norm = VEKNORM_BICGSTAB(r, n);
            EPS = cg_eps * MMax(MMax(x0norm, b0norm), r0norm);
            break;
        case 6:
            MXMatVek(x, t);
            b0norm = VEKNORM_BICGSTAB(b, n);
            x0norm = VEKNORM_BICGSTAB(t, n);
            r0norm = VEKNORM_BICGSTAB(r, n);
            EPS = cg_eps * MMax(MMax(x0norm, b0norm), r0norm);
            break;
    }

    /* Fehler ausgeben */
    DisplayMsg("   Relativer UMF-Anfangsfehler = ");
    DisplayDouble(VEKNORM_BICGSTAB(r, n) / EPS, 8, 3);
    DisplayMsgLn("");

    if (VEKNORM_BICGSTAB(r, n) <= EPS)
    {
        r = (double*)Free(r);
        t = (double*)Free(t);

        /* Ggf. abschliessen der Vorkonditionierung */
        if (vorkond)
            MXVorkond(1, x, b);
        return 0;
    }

    k++;

restart:

    NA = 20 * n * a_size_mult;
    NP = (n + 1) * a_size_mult;

    Ax = (double*)Realloc(Ax, NA * sizeof(double));
    Ai = (int*)Realloc(Ai, NA * sizeof(long));
    Ap = (int*)Realloc(Ap, NP * sizeof(long));

    NNA = 0;
    NNP = 0;

    for (j = 0; j < n; j++)
    {
        /* Spaltenanfang speichern */
        Ap[NNP] = NNA;
        NNP++;
        for (i = 0; i < n; i++)
        {
            /* Ineffektiv!! Besser auf Sparse-Strukturen zugreifen */
            a = MXGet(i, j);
            if (fabs(a) > MFastNull)
            {
                /* Spalteneintraege speichern */
                Ax[NNA] = a;
                Ai[NNA] = i;
                NNA++;
                if (NNA == NA)
                {
                    a_size_mult *= 1.2;
                    Ax = (double*)Free(Ax);
                    Ai = (int*)Free(Ai);
                    Ap = (int*)Free(Ap);
                    goto restart;
                }
            }
        }
    }

    Ap[NNP] = NNA;
    NNP++;

#ifdef UMFPACK31
    umfpack_i_symbolic(n, Ap, Ai, &Symbolic, Control, Info);
    if (Symbolic == NULL)
        DisplayMsgLn("UMFPACK: Symbolic ist NULL!");
    umfpack_i_numeric(Ap, Ai, Ax, Symbolic, &Numeric, Control, Info);
    if (Numeric == NULL)
        DisplayMsgLn("UMFPACK: Numeric ist NULL!");
    status = umfpack_i_solve("Ax=b", Ap, Ai, Ax, x, b, Numeric, Control, Info);
    umfpack_i_report_info(Control, Info);
    umfpack_i_report_status(Control, status);
#endif
#ifdef UMFPACK40
    status = umfpack_di_symbolic(n, n, Ap, Ai, &Symbolic, Control, Info);
    umfpack_di_report_status(Control, status);
    if (status == UMFPACK_ERROR_out_of_memory)
    {
        a_size_mult *= 1.2;
        Ax = (double*)Free(Ax);
        Ai = (int*)Free(Ai);
        Ap = (int*)Free(Ap);
        goto restart;
    }
    if (status != UMFPACK_OK)
    {
        DisplayMsgLn("UMFPACK wird abgebrochen!");
        DisplayMsgLn("GLS wird in ERROR.GLS gespeichert.");
        MXDumpGLS("ERROR.GLS", 1, b, NULL);
        exit(1);
    }
    if (Symbolic == NULL)
        DisplayMsgLn("UMFPACK: Symbolic ist NULL!");

    status = umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, Control, Info);
    umfpack_di_report_status(Control, status);
    if (status != UMFPACK_OK)
    {
        DisplayMsgLn("UMFPACK wird abgebrochen!");
        DisplayMsgLn("GLS wird in ERROR.GLS gespeichert.");
        MXDumpGLS("ERROR.GLS", 1, b, NULL);
        exit(1);
    }
    if (Numeric == NULL)
        DisplayMsgLn("UMFPACK: Numeric ist NULL!");

    status =
        umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, Control, Info);
    umfpack_di_report_info(Control, Info);
    umfpack_di_report_status(Control, status);
    if (status != UMFPACK_OK)
    {
        DisplayMsgLn("UMFPACK wird abgebrochen!");
        DisplayMsgLn("GLS wird in ERROR.GLS gespeichert.");
        MXDumpGLS("ERROR.GLS", 1, b, NULL);
        exit(1);
    }
#endif

#ifdef UMFPACK31
    umfpack_i_free_symbolic(&Symbolic);
    umfpack_i_free_numeric(&Numeric);
#endif
#ifdef UMFPACK40
    umfpack_di_free_symbolic(&Symbolic);
    umfpack_di_free_numeric(&Numeric);
#endif

    MXResiduum(x, b, r);
    DisplayMsg("   Relativer UMF-Endfehler = ");
    DisplayDouble(VEKNORM_BICGSTAB(r, n) / EPS, 8, 3);
    DisplayMsgLn("");

    /* Ggf. abschliessen der Vorkonditionierung */
    if (vorkond)
        MXVorkond(1, x, b);

    Ax = (double*)Free(Ax);
    Ai = (int*)Free(Ai);
    Ap = (int*)Free(Ap);

    r = (double*)Free(r);
    t = (double*)Free(t);

    return k;
#else  // ifdef UMFPACK
    n = n;
    *x = *x;
    *b = *b;
    DisplayMsg("!!!! Error: UMF-Solver not included in this version!");
    DisplayMsgLn(" ");
    exit(1);
#endif
}

/*************************************************************************
   ROCKFLOW - Funktion: SpBICG

   Aufgabe:
   BICG-Loeser fuer LGS der Form A*x=b

   Die Matrix A steht ueber die Speichertechnik implizit zur Verfuegung.
   Die benutzte Vektornorm wird in makros.h gesetzt.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E double *b : Vektor b
   R double *x : Ergebnisvektor x;
   der Speicher muss bereits allokiert sein
   E long g : Dimension der Vektoren

   Ergebnis:
   Anzahl der Iterationen

   Programmaenderungen:
   ??/????     DAHMANE         Erste Version
   11/1995     MSR/DAHMANE     an ROCKFLOW und Speichertechnik angepasst
   6/1997     C.Thorenz  Berechnung des Abbruchfehlers EPS geaendert
   Speicherloch bei vorzeitigem Abbruch beseitigt
   10/1997    C.Thorenz  Neues Abbruchkriterium
   VekNorm(Residuum)/VekNorm(x) < EPS
   9/2000     C.Thorenz  Neues Fehlerkriterium
   01/2001    CT   Neues Fehlerkriterium: max(||Ax||, ||b||, ||r0||)

*************************************************************************/
int SpBICG(double* b, double* x, long n)
{
    /* Variablen */
    static double *r, *rs, *p, *ps, *v, *tmp;
    static double alpha, beta, rho, rho1, eps;
    register int i; /* schnellere Vektoroperationen, Ra, 3/2000 */
    int k = 0, max_iter = 0;
    static double r0norm, b0norm, x0norm;

    /* Ggf. starten der Vorkonditionierung */
    if (vorkond)
        MXVorkond(0, x, b);

    rho1 = rho = 1.0;

    if (cg_maxiter > 0)
        max_iter = cg_maxiter;
        // OK411    if (cg_maxiter == -1)
        // OK411        max_iter = NodeListLength;

#ifdef TESTLOES4
    DisplayMsgLn("SpBICG");
#endif

    r = (double*)Malloc(n * sizeof(double)); /* MNulleVec(r,n); */
    rs = (double*)Malloc(n * sizeof(double));

    MXResiduum(x, b, r);

    if (linear_error_type == 0)
        eps = cg_eps;
    else if (linear_error_type == 1)
        eps = cg_eps * VEKNORM_BICG(b, n);
    else if (linear_error_type == 2)
        eps = cg_eps * VEKNORM_BICG(r, n);
    else if (linear_error_type == 3)
    {
        if ((r0norm = VEKNORM_BICG(r, n)) > 1.)
            eps = cg_eps * r0norm;
        else
            eps = cg_eps;
    }
    else if (linear_error_type == 4)
        eps = cg_eps * (VEKNORM_BICG(x, n));
    else if (linear_error_type == 5)
    {
        b0norm = VEKNORM_BICGSTAB(b, n);
        x0norm = VEKNORM_BICGSTAB(x, n);
        r0norm = VEKNORM_BICGSTAB(r, n);
        eps = cg_eps * MMax(MMax(x0norm, b0norm), r0norm);
    }
    else if (linear_error_type == 6)
    {
        MXMatVek(x, rs);
        b0norm = VEKNORM_BICGSTAB(b, n);
        x0norm = VEKNORM_BICGSTAB(rs, n);
        r0norm = VEKNORM_BICGSTAB(r, n);
        eps = cg_eps * MMax(MMax(x0norm, b0norm), r0norm);
    }

#ifdef TESTLOES4
    DisplayMsg("eps = ");
    DisplayDouble(eps, 22, 20);
    DisplayMsgLn("");
#endif

    if (VEKNORM_BICG(r, n) <= eps)
    {
        r = (double*)Free(r);
        rs = (double*)Free(rs);

        /* Ggf. abschliessen der Vorkonditionierung */
        if (vorkond)
            MXVorkond(1, x, b);
        return 0;
    }
    p = (double*)Malloc(n * sizeof(double));
    ps = (double*)Malloc(n * sizeof(double));
    v = (double*)Malloc(n * sizeof(double));
    tmp = (double*)Malloc(n * sizeof(double));

    if
        VK_Modus(VK_iLDU) MXVorkond(2, r, r); /* Ra, 3/2000 */
    for (i = 0; i < n; i++)
    {
        p[i] = ps[i] = 0.0;
        rs[i] = r[i];
    }

    for (;;)
    {
        rho = MSkalarprodukt(rs, r, n);
        beta = rho / rho1;
        for (i = 0; i < n; i++)
        {
            p[i] = p[i] * beta + r[i];
            ps[i] = ps[i] * beta + rs[i];
        }

        MXMatVek(p, v);
        if
            VK_Modus(VK_iLDU) MXVorkond(2, v, v); /* Ra, 3/2000 */
        alpha = rho / MSkalarprodukt(ps, v, n);
        for (i = 0; i < n; i++)
        {
            x[i] += alpha * p[i];
            r[i] -= alpha * v[i];
        }

#ifdef SOLVER_SHOW_RESULTS
        printf("\n%ld %f %f %f %f %f", (long)k, x[(long)(n * .1)],
               x[(long)(n * .3)], x[(long)(n * .5)], x[(long)(n * .7)],
               x[(long)(n * .9)]);
#endif

        k++;
        if (linear_error_type == 4)
            eps = cg_eps * VEKNORM_BICG(x, n);
        if (VEKNORM_BICG(r, n) <= eps)
            break;

        if
            VK_Modus(VK_iLDU)
            {
                MXVorkond(3, v, ps);
                MXMatTVek(v, tmp);
            }
        else
            MXMatTVek(ps, tmp); /* Ra, 3/2000 */
        for (i = 0; i < n; i++)
            rs[i] -= alpha * tmp[i];

        if (k >= max_iter)
            break;

#ifdef SOLVER_SHOW_ERROR
        if (k % (int)MMax((cg_maxiter / 10), 1.) == 1)
        {
            DisplayMsg("      Iteration-Nr.: ");
            DisplayLong((long)k);
            DisplayMsg(", Fehler = ");
            DisplayDouble(VEKNORM_BICGSTAB(r, n) / eps, 4, 1);
            DisplayMsgLn(" ");
        }
#endif

        rho1 = rho;
    }

    r = (double*)Free(r);
    rs = (double*)Free(rs);
    p = (double*)Free(p);
    ps = (double*)Free(ps);
    v = (double*)Free(v);
    tmp = (double*)Free(tmp);
    /* msr: evtl. Speicher immer nur vergroessern, nie ganz freigeben
       (wie bei masscont geplant) */

    /* Ggf. abschliessen der Vorkonditionierung */
    if (vorkond)
        MXVorkond(1, x, b);

    return k;
}

/*************************************************************************
   ROCKFLOW - Funktion: SpBICGSTAB

   Aufgabe:
   BICGSTAB-Loeser fuer LGS der Form A*x=b

   Die Matrix A steht ueber die Speichertechnik implizit zur Verfuegung.
   Die benutzte Vektornorm wird in makros.h gesetzt.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E double *b : Vektor b
   R double *x : Ergebnisvektor x;
   der Speicher muss bereits allokiert sein
   E long g : Dimension der Vektoren

   Ergebnis:
   Anzahl der Iterationen

   Programmaenderungen:

   ??/????   DAHMANE       Erste Version
   11/1995   MSR/DAHMANE   an ROCKFLOW und Speichertechnik angepasst
   6/1997   C.Thorenz  Berechnung des Abbruchfehlers EPS geaendert
   Speicherloch bei vorzeitigem Abbruch beseitigt
   10/1997   C.Thorenz  Neues Abbruchkriterium
   VekNorm(Residuum)/VekNorm(x) < EPS
   9/2000   C.Thorenz  Neues Fehlerkriterium
   11/2000   C.Thorenz  0/0-Problem hoffentlich besiegt
   1/2001   C.Thorenz  Neues Fehlerkriterium: max(||Ax||, ||b||, ||r0||)
   1/2001   C.Thorenz  Ueberarbeitet nach Barrett, Berry, Chan et al.
   (Basierend auf bicgstab.h C++-Template)

*************************************************************************/
#ifdef SX
int SpBICGSTAB(double* restrict b, double* restrict x, long n)
#else
int SpBICGSTAB(double* b, double* x, long n)
#endif
{
/* Variablen */
#ifdef SX
    double *restrict r, *restrict r2, *restrict rs, *restrict p, *restrict s,
        *restrict t, *restrict v;
#else
    double *r, *r2, *rs, *p, *s, *t, *v;
#endif
    double alpha, beta, omega, rho, rho1, eps = 0.;
    register long i; /* schnellere Vektoroperationen, Ra, 3/2000 */
    int k = 0, max_iter = 0, repeat = 0;
    double r0norm = 0., b0norm = 0., x0norm = 0., tt, ts, rsv;
    // WW double error_rel;
    // MXDumpGLS("rf_pcs.txt",1,b,x); abort();
    /* Ggf. starten der Vorkonditionierung */
    if (vorkond)
        // WW      cout << "        Preconditioning" << "\n";
        MXVorkond(0, x, b);
    r = (double*)Malloc(n * sizeof(double));
    r2 = (double*)Malloc(n * sizeof(double));
    rs = (double*)Malloc(n * sizeof(double));
    p = (double*)Malloc(n * sizeof(double));
    s = (double*)Malloc(n * sizeof(double));
    t = (double*)Malloc(n * sizeof(double));
    v = (double*)Malloc(n * sizeof(double));

    if (cg_maxiter > 0)
        max_iter = cg_maxiter;
    // OK411    if (cg_maxiter == -1)
    // OK411        max_iter = NodeListLength;

restart:

    /* Zaehler fuer Wiederholungen des gesamten Loesers */
    repeat++;

    alpha = omega = rho1 = rho = 1.0;

#ifdef TESTLOES4
    DisplayMsgLn("SpBICGSTAB");
#endif

    MXResiduum(x, b, r);

    if (linear_error_type == 0)
        eps = cg_eps;
    else if (linear_error_type == 1)
        eps = cg_eps * VEKNORM_BICGSTAB(b, n);
    else if (linear_error_type == 2)
        eps = cg_eps * VEKNORM_BICGSTAB(r, n);
    else if (linear_error_type == 3)
    {
        if ((r0norm = VEKNORM_BICGSTAB(r, n)) > 1.)
            eps = cg_eps * r0norm;
        else
            eps = cg_eps;
    }
    else if (linear_error_type == 4)
        eps = cg_eps * (VEKNORM_BICGSTAB(x, n));
    else if (linear_error_type == 5)
    {
        b0norm = VEKNORM_BICGSTAB(b, n);
        x0norm = VEKNORM_BICGSTAB(x, n);
        r0norm = VEKNORM_BICGSTAB(r, n);
        eps = cg_eps * MMax(MMax(x0norm, b0norm), r0norm);
    }
    else if (linear_error_type == 6)
    {
        MXMatVek(x, t);
        b0norm = VEKNORM_BICGSTAB(b, n);
        x0norm = VEKNORM_BICGSTAB(t, n);
        r0norm = VEKNORM_BICGSTAB(r, n);
        eps = cg_eps * MMax(MMax(x0norm, b0norm), r0norm);
    }
    // OK
    r0norm = VEKNORM_BICGSTAB(r, n);
// WW error_rel = r0norm/eps;
// WW    cout << "\n  SpBICGSTAB iteration: 0/" << max_iter << " Error: " <<
// error_rel << "\n";
#ifdef TESTLOES4
    DisplayMsg("eps = ");
    DisplayDouble(eps, 22, 20);
    DisplayMsgLn("");
#endif
    /* Fehlerkriterium oder max_iter erreicht?  */
    if ((k >= max_iter) || (VEKNORM_BICGSTAB(r, n) <= eps))
        goto end;

    if
        VK_Modus(VK_iLDU) MXVorkond(2, r, r);

    MKopierVec(r, rs, n);
    MKopierVec(r, p, n);
    MNulleVec(v, n);

    for (;;)
    {
        /* Zaehler fuer Wiederholungen innerhalb des Loesers */
        k++;
        if (k > max_iter)
            goto end;

        rho = MSkalarprodukt(rs, r, n);

        if (fabs(rho) < DBL_MIN)
            /* Bei rho==0 ist eigentlich alles hoffnungslos, wir versuchen einen
             * Restart ... */
            goto restart;

        if (k > 1)
        {
            beta = (rho / rho1) * (alpha / omega);
#ifdef SX
#pragma cdir nodep
#endif
            for (i = 0; i < n; i++)
                p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }

        MXMatVek(p, v);
        if
            VK_Modus(VK_iLDU) MXVorkond(2, v, v); /* Ra, 3/2000 */

        rsv = MSkalarprodukt(rs, v, n);

        alpha = rho / rsv;
#ifdef SX
#pragma cdir nodep
#endif
        for (i = 0; i < n; i++)
            s[i] = r[i] - alpha * v[i];

        /* Erste Abbruchmoeglichkeit */
        if (VEKNORM_BICGSTAB(s, n) <= eps)
        {
/* Die Loesung nochmal updaten */
#ifdef SX
#pragma cdir nodep
#endif
            for (i = 0; i < n; i++)
                x[i] += alpha * p[i];
            break;
        }

        MXMatVek(s, t);

        if
            VK_Modus(VK_iLDU) MXVorkond(2, t, t); /* Ra, 3/2000 */

        ts = MSkalarprodukt(t, s, n);
        tt = MSkalarprodukt(t, t, n);

        /* tt ueberpruefen */
        if (fabs(tt) > DBL_MIN)
        {
            if ((log(fabs(ts)) - log(tt)) < log(DBL_MAX))
                /* Alles o.k. */
                omega = ts / tt;
            else
                /* Division unzulaessig, nur Vorzeichen behalten ist sinnvoll */
                omega = (double)Signum(ts) / (double)Signum(tt) / MKleinsteZahl;
        }
        else
            /* tt ist Null, darf eigentlich nicht vorkommen */
            omega = 1.;

#ifdef SX
#pragma cdir nodep
#endif
        for (i = 0; i < n; i++)
            x[i] += alpha * p[i] + omega * s[i];

        MXResiduum(x, b, r2);
#ifdef SX
#pragma cdir nodep
#endif
        for (i = 0; i < n; i++)
            r[i] = s[i] - omega * t[i];

#ifdef SOLVER_SHOW_RESULTS
        printf("\n%ld %f %f %f %f %f", (long)k, x[(long)(n * .1)],
               x[(long)(n * .3)], x[(long)(n * .5)], x[(long)(n * .7)],
               x[(long)(n * .9)]);
#endif
        // OK
        // WW error_rel = VEKNORM_BICGSTAB(r2,n)/eps;
        // WW    printf("\r        SpBICGSTAB iteration: %i/%i Error:
        // %f",k,max_iter,error_rel);
        if (linear_error_type == 4)
            eps = cg_eps * (VEKNORM_BICGSTAB(x, n));
        if (VEKNORM_BICGSTAB(r2, n) <= eps)
            break;

#ifdef SOLVER_SHOW_ERROR
        if (k % (int)MMax((max_iter / 10), 1.) == 1)
        {
            DisplayMsg("      Iteration-Nr.: ");
            DisplayLong((long)k);
            DisplayMsg(", Fehler = ");
            DisplayDouble(VEKNORM_BICGSTAB(r2, n) / eps, 4, 1);
            DisplayMsgLn(" ");
        }
#endif

        rho1 = rho;
    }

    /* Den Loeser nochmal ausfuehren? */
    if ((repeat <= cg_repeat) || (cg_repeat == -1))
        goto restart;
    else
        goto end;

end:
    r = (double*)Free(r);
    r2 = (double*)Free(r2);
    rs = (double*)Free(rs);
    p = (double*)Free(p);
    s = (double*)Free(s);
    t = (double*)Free(t);
    v = (double*)Free(v);

    /* Ggf. abschliessen der Vorkonditionierung */
    if (vorkond)
        MXVorkond(1, x, b);

    // WW
    printf("\t  SpBICGSTAB iteration: %i/%i \n", k, max_iter);

    return k;
}

/*************************************************************************
   ROCKFLOW - Funktion: SpQMRCGSTAB

   Aufgabe:
   QMRCGSTAB-Loeser fuer LGS der Form A*x=b

   Die Matrix A steht ueber die Speichertechnik implizit zur Verfuegung.
   Die benutzte Vektornorm wird in makros.h gesetzt.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E double *b : Vektor b
   R double *x : Ergebnisvektor x;
   der Speicher muss bereits allokiert sein
   E long g : Dimension der Vektoren

   Ergebnis:
   Anzahl der Iterationen

   Programmaenderungen:
   ??/????     DAHMANE         Erste Version
   11/1995     MSR/DAHMANE     an ROCKFLOW und Speichertechnik angepasst
   6/1997      C.Thorenz  Berechnung des Abbruchfehlers EPS geaendert
   Speicherloch bei vorzeitigem Abbruch beseitigt
   10/1997     C.Thorenz  Neues Abbruchkriterium
   VekNorm(Residuum)/VekNorm(x) < EPS
   9/2000      C.Thorenz  Neues Fehlerkriterium
   1/2001      CT   Neues Fehlerkriterium: max(||Ax||, ||b||, ||r0||)

*************************************************************************/
int SpQMRCGSTAB(double* b, double* x, long n)
{
    /* Variablen */
    static double *r, *rs, *d, *ds, *p, *s, *t, *v, *xs;
    static double alpha, beta, omega, rho, rho1, eps;
    static double tau, teta, eta, taus, tetas, etas, c;
    register int i; /* schnellere Vektoroperationen, Ra, 3/2000 */
    int k = 0, max_iter = 0;
    static double r0norm, b0norm, x0norm;

    /* Ggf. starten der Vorkonditionierung */
    if (vorkond)
        MXVorkond(0, x, b);

    alpha = 1.0;
    omega = 1.0;
    rho1 = rho = 1.0;
    teta = eta = 0.0;

    if (cg_maxiter > 0)
        max_iter = cg_maxiter;
        // OK411    if (cg_maxiter == -1)
        // OK411        max_iter = NodeListLength;

#ifdef TESTLOES4
    DisplayMsgLn("SpQMRCGSTAB");
#endif

    r = (double*)Malloc(n * sizeof(double)); /* MNulleVec(r,n); */
    s = (double*)Malloc(n * sizeof(double));
    MXResiduum(x, b, r);

    if (linear_error_type == 0)
        eps = cg_eps;
    else if (linear_error_type == 1)
        eps = cg_eps * VEKNORM_QMRCGSTAB(b, n);
    else if (linear_error_type == 2)
        eps = cg_eps * VEKNORM_QMRCGSTAB(r, n);
    else if (linear_error_type == 3)
    {
        if ((r0norm = VEKNORM_QMRCGSTAB(r, n)) > 1.)
            eps = cg_eps * r0norm;
        else
            eps = cg_eps;
    }
    else if (linear_error_type == 4)
        eps = cg_eps * (VEKNORM_QMRCGSTAB(x, n));
    else if (linear_error_type == 5)
        eps = cg_eps * MMax(MMax(x0norm, b0norm), r0norm);
    else if (linear_error_type == 6)
    {
        MXMatVek(x, s);
        b0norm = VEKNORM_BICGSTAB(b, n);
        x0norm = VEKNORM_BICGSTAB(s, n);
        r0norm = VEKNORM_BICGSTAB(r, n);
        eps = cg_eps * MMax(MMax(x0norm, b0norm), r0norm);
    }

#ifdef TESTLOES4
    DisplayMsg("eps = ");
    DisplayDouble(eps, 22, 20);
    DisplayMsgLn("");
#endif

    if (VEKNORM_QMRCGSTAB(r, n) <= eps)
    {
        r = (double*)Free(r);
        s = (double*)Free(s);
        /* Ggf. abschliessen der Vorkonditionierung */
        if (vorkond)
            MXVorkond(1, x, b);
        return 0;
    }
    rs = (double*)Malloc(n * sizeof(double));
    d = (double*)Malloc(n * sizeof(double));
    ds = (double*)Malloc(n * sizeof(double));
    p = (double*)Malloc(n * sizeof(double));
    s = (double*)Malloc(n * sizeof(double));
    t = (double*)Malloc(n * sizeof(double));
    v = (double*)Malloc(n * sizeof(double));
    xs = (double*)Malloc(n * sizeof(double));
    for (i = 0; i < n; i++)
        d[i] = p[i] = v[i] = 0.0;

    if
        VK_Modus(VK_iLDU) MXVorkond(2, r, r); /* Ra, 3/2000 */
    for (i = 0; i < n; i++)
        rs[i] = r[i];

    for (;;)
    {
        rho = MSkalarprodukt(rs, r, n);
        beta = (rho / rho1) * (alpha / omega);
        for (i = 0; i < n; i++)
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
        MXMatVek(p, v);
        if
            VK_Modus(VK_iLDU) MXVorkond(2, v, v); /* Ra, 3/2000 */
        alpha = rho / MSkalarprodukt(rs, v, n);
        for (i = 0; i < n; i++)
            s[i] = r[i] - alpha * v[i];
        /* first quasi-minimization */
        tetas = VEKNORM_QMRCGSTAB(s, n) / tau;
        c = 1. / sqrt(1. + tetas * tetas);
        taus = tau * tetas * c;
        etas = c * c * alpha;
        c = teta * teta * eta / alpha;
        for (i = 0; i < n; i++)
        {
            ds[i] = p[i] + c * d[i];
            xs[i] = x[i] + etas * ds[i];
        }
        /* update r */
        MXMatVek(s, t);
        if
            VK_Modus(VK_iLDU) MXVorkond(2, t, t); /* Ra, 3/2000 */
        omega = MSkalarprodukt(t, s, n) / MSkalarprodukt(t, t, n);
        for (i = 0; i < n; i++)
            r[i] = s[i] - omega * t[i];

        /* second quasi-minimization */
        teta = VEKNORM_QMRCGSTAB(r, n) / taus;
        c = 1. / sqrt(1 + teta * teta);
        tau = taus * teta * c;
        eta = c * c * omega;
        c = tetas * tetas * etas / omega;
        for (i = 0; i < n; i++)
        {
            d[i] = s[i] + c * ds[i];
            x[i] = xs[i] + eta * d[i];
        }

#ifdef SOLVER_SHOW_RESULTS
        printf("\n%ld %f %f %f %f %f", (long)k, x[(long)(n * .1)],
               x[(long)(n * .3)], x[(long)(n * .5)], x[(long)(n * .7)],
               x[(long)(n * .9)]);
#endif
        k++;

        if (linear_error_type == 4)
            eps = cg_eps * (VEKNORM_QMRCGSTAB(x, n));
        if (VEKNORM_QMRCGSTAB(r, n) <= eps)
            break;

        if (k >= max_iter)
            break;

#ifdef SOLVER_SHOW_ERROR
        if (k % (int)MMax((cg_maxiter / 10), 1.) == 1)
        {
            DisplayMsg("      Iteration-Nr.: ");
            DisplayLong((long)k);
            DisplayMsg(", Fehler = ");
            DisplayDouble(VEKNORM_BICGSTAB(r, n) / eps, 4, 1);
            DisplayMsgLn(" ");
        }
#endif
        rho1 = rho;
    }

    r = (double*)Free(r);
    rs = (double*)Free(rs);
    d = (double*)Free(d);
    ds = (double*)Free(ds);
    p = (double*)Free(p);
    s = (double*)Free(s);
    t = (double*)Free(t);
    v = (double*)Free(v);
    xs = (double*)Free(xs);
    /* msr: evtl. Speicher immer nur vergroessern, nie ganz freigeben
       (wie bei masscont geplant) */

    /* Ggf. abschliessen der Vorkonditionierung */
    if (vorkond)
        MXVorkond(1, x, b);

    return k;
}

/*************************************************************************
   ROCKFLOW - Funktion: SpMGMRES

   Aufgabe:
   MGMRES-Loeser fuer LGS der Form A*x=b

   Die Matrix A steht ueber die Speichertechnik implizit zur Verfuegung.

   ... weitere Angaben

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E double *b : Vektor b
   R double *x : Ergebnisvektor x;
   der Speicher muss bereits allokiert sein
   E long g : Dimension der Vektoren

   Ergebnis:
   Anzahl der Iterationen

   Programmaenderungen:
   3/2002    C.Thorenz   Erste Version, funktioniert noch nicht

*************************************************************************/
int SpMGMRES(double* b, double* x, long n)
{
#ifdef GMRES_WORKS
    double *V, *U, *r, *x0, *y, *c, *cs, *sn, *s, **v, *w, *help;
    double beta, tol = cg_eps, resid, normb;
    long max_iter, i, i1, i2, j, k;
    const long m = 30;
    double H[31][31];

    if (n <= 0)
        return -1;

    V = MMachVec(n * (m + 1));
    U = MMachVec(m * (m + 1) / 2);
    r = MMachVec(n);
    x0 = MMachVec(n);
    w = MMachVec(n);
    help = MMachVec(n);
    y = MMachVec(m + 1);
    c = MMachVec(m);
    s = MMachVec(m + 1);
    cs = MMachVec(m + 1);
    sn = MMachVec(m + 1);
    v = (double**)Malloc((m + 1) * sizeof(double*)); /* Pointerfeld auf V */

    for (i = 0; i <= m; ++i)
        v[i] = V + i * n; /* Zuweisen */

    if (cg_maxiter >= 0)
        max_iter = cg_maxiter;
    if (cg_maxiter == -1)
        max_iter = NodeListLength;

    if
        VK_Modus(VK_iLDU) MXVorkond(2, b, help);
    normb = VEKNORM_BICGSTAB(help, n);

    if (normb == 0.0)
        normb = 1;

    MKopierVec(x, x0, n);

    j = 1;
    while (j < max_iter)
    {
        /*
           "aussere Iteration
           Muss hier x0 oder x verwendet werden???
           r = M.solve(b - A * x);
         */
        MXResiduum(x0, b, r);
        if
            VK_Modus(VK_iLDU) MXVorkond(2, r, r);

        beta = VEKNORM_BICGSTAB(r, n);

        if ((resid = beta / normb) < tol)
        {
            tol = resid;
            max_iter = j;
            goto end;
        }

        MVekGle(0., r, 1. / beta, r, v[0], n);
        MNulleVec(s, m + 1);
        s[0] = beta;

        for (i = 0; i < m && j <= max_iter; i++, j++)
        {
            /*      w = M.solve(A * v[i]); */
            MXMatVek(v[i], w);
            if
                VK_Modus(VK_iLDU) MXVorkond(2, w, w);

            for (k = 0; k <= i; k++)
            {
                H[k][i] = MSkalarprodukt(w, v[k], n);
                MVekGle(1., w, -H[k][i], v[k], w, n);
            }

            H[i + 1][i] = VEKNORM_BICGSTAB(w, n);
            MVekGle(0., w, 1. / H[i + 1][i], w, v[i + 1], n);

            for (k = 0; k < i; k++)
                ApplyPlaneRotation(&H[k][i], &H[k + 1][i], &cs[k], &sn[k]);

            GeneratePlaneRotation(&H[i][i], &H[i + 1][i], &cs[i], &sn[i]);
            ApplyPlaneRotation(&H[i][i], &H[i + 1][i], &cs[i], &sn[i]);
            ApplyPlaneRotation(&s[i], &s[i + 1], &cs[i], &sn[i]);

            if ((resid = fabs(s[i + 1]) / normb) < tol)
            {
                /*        Update(x, i, H, s, v); */
                tol = resid;
                max_iter = j;
                goto end;
            }
        }
        /* Backsolve:  Update(x, m - 1, H, s, v); */
        for (i1 = m - 1; i1 >= 0; i1--)
        {
            y[i1] /= H[i1][i1];
            for (i2 = i1 - 1; i2 >= 0; i2--)
                y[i2] -= H[i2][i1] * y[i1];
        }

        for (i2 = 0; i2 <= m - 1; i2++)
        {
            MKopierVec(v[i2], help, n);
            MMultVecSkalar(help, y[i2], n);
            MAddVektoren(x, help, x, n);
        }
    }

end:
    MLoeschVec(V);
    MLoeschVec(U);
    MLoeschVec(r);
    MLoeschVec(y);
    MLoeschVec(c);
    MLoeschVec(s);
    v = (double**)Free(v);

    return max_iter;
#else  // ifdef GMRES_WORKS
    n = n;
    *x = *x;
    *b = *b;
    DisplayMsg("!!!! Error: GMRES-Solver not yet finished!");
    DisplayMsgLn(" ");
    exit(1);
#endif
}

#ifdef GMRES_WORKS
void GeneratePlaneRotation(double* dx, double* dy, double* cs, double* sn)
{
    static double temp;
    if (*dy == 0.0)
    {
        *cs = 1.0;
        *sn = 0.0;
    }
    else if (abs(*dy) > abs(*dx))
    {
        temp = *dx / *dy;
        *sn = 1.0 / sqrt(1.0 + temp * temp);
        *cs = temp * *sn;
    }
    else
    {
        temp = *dy / *dx;
        *cs = 1.0 / sqrt(1.0 + temp * temp);
        *sn = temp * *cs;
    }
}

void ApplyPlaneRotation(double* dx, double* dy, double* cs, double* sn)
{
    static double temp;
    temp = *cs * *dx + *sn * *dy;
    *dy = -*sn * *dx + *cs * *dy;
    *dx = temp;
}
#endif

/*************************************************************************
   ROCKFLOW - Funktion: SpCG

   Aufgabe:
   CG-Loeser fuer LGS der Form A*x=b wobei A sym. und pos. def.

   Die Matrix A steht ueber die Speichertechnik implizit zur Verfuegung.
   Die benutzte Vektornorm wird in makros.h gesetzt.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E double *b : Vektor b
   R double *x : Ergebnisvektor x;
   der Speicher muss bereits allokiert sein
   E long g : Dimension der Vektoren

   Ergebnis:
   Anzahl der Iterationen

   Programmaenderungen:
   ??/????     DAHMANE         Erste Version
   09/1997     AH              Einbau der symmetrischen Loeser
   9/2000     C.Thorenz  Neues Fehlerkriterium
*************************************************************************/
int SpCG(double* b, double* x, long n)
{
    /* Variablen */
    static double *d, *r, *s, *tmp;
    static double alpha, beta, tmpr, eps;
    int k = 0, max_iter = 0;
    register int i; /* schnellere Vektoroperationen, Ra, 3/2000 */
    static double r0norm, b0norm, x0norm;

    /* Ggf. starten der Vorkonditionierung */
    if (vorkond)
        MXVorkond(0, x, b);

    if (cg_maxiter > 0)
        max_iter = cg_maxiter;
        // OK411    if (cg_maxiter == -1)
        // OK411        max_iter = NodeListLength;

#ifdef TESTLOES4
    DisplayMsgLn("SpCG");
#endif

    r = (double*)Malloc(n * sizeof(double));
    s = (double*)Malloc(n * sizeof(double));

    MXResiduum(x, b, r);
    if (linear_error_type == 0)
        eps = cg_eps;
    else if (linear_error_type == 1)
        eps = cg_eps * VEKNORM_CG(b, n);
    else if (linear_error_type == 2)
        eps = cg_eps * VEKNORM_CG(r, n);
    else if (linear_error_type == 3)
    {
        if ((r0norm = VEKNORM_CG(r, n)) > 1.)
            eps = cg_eps * r0norm;
        else
            eps = cg_eps;
    }
    else if (linear_error_type == 4)
        eps = cg_eps * (VEKNORM_CG(x, n));
    else if (linear_error_type == 5)
        eps = cg_eps * MMax(MMax(x0norm, b0norm), r0norm);
    else if (linear_error_type == 6)
    {
        MXMatVek(x, s);
        b0norm = VEKNORM_BICGSTAB(b, n);
        x0norm = VEKNORM_BICGSTAB(s, n);
        r0norm = VEKNORM_BICGSTAB(r, n);
        eps = cg_eps * MMax(MMax(x0norm, b0norm), r0norm);
    }

    s = (double*)Free(s);

#ifdef TESTLOES4
    DisplayMsg("eps = ");
    DisplayDouble(eps, 22, 20);
    DisplayMsgLn("");
#endif

    if (VEKNORM_CG(r, n) <= eps)
    {
        r = (double*)Free(r);
        /* Ggf. abschliessen der Vorkonditionierung */
        if (vorkond)
            MXVorkond(1, x, b);
        return 0;
    }
    tmp = (double*)Malloc(n * sizeof(double));
    d = (double*)Malloc(n * sizeof(double));

    if
        VK_Modus(VK_iLDU) MXVorkond(2, r, r); /* Ra, 3/2000 */
    for (i = 0; i < n; i++)
        d[i] = r[i];

    for (;;)
    {
        MXMatVek(d, tmp);
        if
            VK_Modus(VK_iLDU) MXVorkond(2, tmp, tmp); /* Ra, 3/2000 */

        alpha = (tmpr = MSkalarprodukt(r, r, n)) / MSkalarprodukt(d, tmp, n);
        for (i = 0; i < n; i++)
        {
            x[i] += alpha * d[i];
            r[i] -= alpha * tmp[i];
        }

#ifdef SOLVER_SHOW_RESULTS
        printf("\n%ld %f %f %f %f %f", (long)k, x[(long)(n * .1)],
               x[(long)(n * .3)], x[(long)(n * .5)], x[(long)(n * .7)],
               x[(long)(n * .9)]);
#endif

        k++;
        if (linear_error_type == 4)
            eps = cg_eps * VEKNORM_CG(x, n);
        if (VEKNORM_CG(r, n) <= eps)
            break;
        if (k >= max_iter)
            break;

#ifdef SOLVER_SHOW_ERROR
        if (k % (int)MMax((cg_maxiter / 10), 1.) == 1)
        {
            DisplayMsg("      Iteration-Nr.: ");
            DisplayLong((long)k);
            DisplayMsg(", Fehler = ");
            DisplayDouble(VEKNORM_BICGSTAB(r, n) / eps, 4, 1);
            DisplayMsgLn(" ");
        }
#endif

        beta = MSkalarprodukt(r, r, n) / tmpr;
        for (i = 0; i < n; i++)
            d[i] = r[i] + beta * d[i];
    }

    r = (double*)Free(r);
    tmp = (double*)Free(tmp);
    d = (double*)Free(d);

    /* Ggf. abschliessen der Vorkonditionierung */
    if (vorkond)
        MXVorkond(1, x, b);

    printf("\t  SpCG iteration: %i/%i \n", k, max_iter);

    return k;
}

/*************************************************************************
   ROCKFLOW - Funktion: SpCGNR

   Aufgabe:
   CG-Loeser fuer LGS der Form A*x=b wobei A beliebige reelle unsymmetrische
   Matrix ist.
   Das Verfahren beruht auf die Anwendung auf die normale Gleichung:
   A'A x = A'b
   Die Idee ist ganz simple: Ist A regulaer so ist A'A symmetrish
   und positiv definit.
   Der Algorithmus wurde so weit modifiziert, dass die Matrixmultiplikation
   A'A vermieden wird. Stichwort: Rundungsfehler.

   Die Matrix A steht ueber die Speichertechnik implizit zur Verfuegung.
   Die benutzte Vektornorm wird in makros.h gesetzt.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E double *b : Vektor b
   R double *x : Ergebnisvektor x;
   der Speicher muss bereits allokiert sein
   E long g : Dimension der Vektoren

   Ergebnis:
   Anzahl der Iterationen

   Programmaenderungen:
   ??/????     DAHMANE         Erste Version
   02/1998     A. Habbar       Einbau in Rockflow

*************************************************************************/
int SpCGNR(double* b, double* x, long n)
{
    /* Variablen */
    static double *d, *r, *s, *tmp, *h = NULL;
    static double alpha, beta, tmpr, tmpr1, eps;
    int k = 0, max_iter = 0;
    register int i; /* schnellere Vektoroperationen, Ra, 3/2000 */
    static double r0norm, b0norm, x0norm;

    /* Ggf. starten der Vorkonditionierung */
    if (vorkond)
        MXVorkond(0, x, b);

    if (cg_maxiter > 0)
        max_iter = cg_maxiter;
        // OK411    if (cg_maxiter == -1)
        // OK411        max_iter = NodeListLength;

#ifdef TESTLOES4
    DisplayMsgLn("SpCGNE");
#endif

    r = (double*)Malloc(n * sizeof(double));
    s = (double*)Malloc(n * sizeof(double));

    MXResiduum(x, b, r);
    if (linear_error_type == 0)
        eps = cg_eps;
    else if (linear_error_type == 1)
        eps = cg_eps * VEKNORM_CG(b, n);
    else if (linear_error_type == 2)
        eps = cg_eps * VEKNORM_CG(r, n);
    else if (linear_error_type == 3)
    {
        if ((r0norm = VEKNORM_CG(r, n)) > 1.)
            eps = cg_eps * r0norm;
        else
            eps = cg_eps;
    }
    else if (linear_error_type == 4)
        eps = cg_eps * (VEKNORM_CG(x, n));
    else if (linear_error_type == 5)
        eps = cg_eps * MMax(MMax(x0norm, b0norm), r0norm);
    else if (linear_error_type == 6)
    {
        MXMatVek(x, s);
        b0norm = VEKNORM_BICGSTAB(b, n);
        x0norm = VEKNORM_BICGSTAB(s, n);
        r0norm = VEKNORM_BICGSTAB(r, n);
        eps = cg_eps * MMax(MMax(x0norm, b0norm), r0norm);
    }

    s = (double*)Free(s);

    if (VEKNORM_CG(r, n) <= eps)
    {
        r = (double*)Free(r);
        /* Ggf. abschliessen der Vorkonditionierung */
        if (vorkond)
            MXVorkond(1, x, b);
        return 0;
    }
    tmp = (double*)Malloc(n * sizeof(double));
    d = (double*)Malloc(n * sizeof(double));
    if
        VK_Modus(VK_iLDU) h = (double*)Malloc(n * sizeof(double));

    if
        VK_Modus(VK_iLDU)
        {
            MXVorkond(2, r, r);
            MXVorkond(3, h, r); /* Ra, 3/2000 */
            MXMatTVek(h, tmp);
        }
    else
        MXMatTVek(r, tmp);

    tmpr = MSkalarprodukt(tmp, tmp, n);
    for (i = 0; i < n; i++)
        d[i] = tmp[i];

    for (;;)
    {
        MXMatVek(d, tmp);
        if
            VK_Modus(VK_iLDU) MXVorkond(2, tmp, tmp); /* Ra, 3/2000 */

        alpha = tmpr / MSkalarprodukt(tmp, tmp, n);
        for (i = 0; i < n; i++)
        {
            x[i] += alpha * d[i];
            r[i] -= alpha * tmp[i];
        }

#ifdef SOLVER_SHOW_RESULTS
        printf("\n%ld %f %f %f %f %f", (long)k, x[(long)(n * .1)],
               x[(long)(n * .3)], x[(long)(n * .5)], x[(long)(n * .7)],
               x[(long)(n * .9)]);
#endif
        k++;
        if (linear_error_type == 4)
            eps = cg_eps * VEKNORM_BICG(x, n);
        if (VEKNORM_CG(r, n) <= eps)
            break;
        if (k >= max_iter)
            break;

#ifdef SOLVER_SHOW_ERROR
        if (k % (int)MMax((cg_maxiter / 10), 1.) == 1)
        {
            DisplayMsg("      Iteration-Nr.: ");
            DisplayLong((long)k);
            DisplayMsg(", Fehler = ");
            DisplayDouble(VEKNORM_BICGSTAB(r, n) / eps, 4, 1);
            DisplayMsgLn(" ");
        }
#endif

        if
            VK_Modus(VK_iLDU)
            {
                MXVorkond(3, h, r); /* Ra, 3/2000 */
                MXMatTVek(h, tmp);
            }
        else
            MXMatTVek(r, tmp);
        beta = (tmpr1 = MSkalarprodukt(tmp, tmp, n)) / tmpr;
        for (i = 0; i < n; i++)
            d[i] = tmp[i] + beta * d[i];
        tmpr = tmpr1;
    }

    r = (double*)Free(r);
    tmp = (double*)Free(tmp);
    d = (double*)Free(d);
    h = (double*)Free(h);

    /* Ggf. abschliessen der Vorkonditionierung */
    if (vorkond)
        MXVorkond(1, x, b);

    return k;
}

/*************************************************************************
   ROCKFLOW - Funktion: SpCGS

   Aufgabe:
   CGS-Loeser fuer LGS der Form A*x=b wobei A beliebige reelle unsymmetrische
   Matrix ist.
   Das Verfahren hat den Vorteil im Vergleich zu BICG, dass die
   transponierte Matrix nicht mehr benoetigt wird.

   Die Matrix A steht ueber die Speichertechnik implizit zur Verfuegung.
   Die benutzte Vektornorm wird in makros.h gesetzt.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E double *b : Vektor b
   R double *x : Ergebnisvektor x;
   der Speicher muss bereits allokiert sein
   E long g : Dimension der Vektoren

   Ergebnis:
   Anzahl der Iterationen

   Programmaenderungen:
   ??/????     DAHMANE         Erste Version
   02/1998     A. Habbar       Einbau in Rockflow

*************************************************************************/
int SpCGS(double* b, double* x, long n)
{
    /* Variablen */
    static double *r, *s, *rs, *p, *q, *u, *v, *w, *us, *tmp;
    static double alpha, beta, rho, rho1;
    static double eps;
    static double r0norm, b0norm, x0norm;
    register int i; /* schnellere Vektoroperationen, Ra, 3/2000 */
    int k = 0, max_iter = 0;

    /* Ggf. starten der Vorkonditionierung */
    if (vorkond)
        MXVorkond(0, x, b);

    rho1 = rho = 1.0;

    if (cg_maxiter > 0)
        max_iter = cg_maxiter;
        // OK411    if (cg_maxiter == -1)
        // OK411        max_iter = NodeListLength;

#ifdef TESTLOES4
    DisplayMsgLn("SpCGS");
#endif

    r = (double*)Malloc(n * sizeof(double));
    s = (double*)Malloc(n * sizeof(double));

    MXResiduum(x, b, r);
    if (linear_error_type == 0)
        eps = cg_eps;
    else if (linear_error_type == 1)
        eps = cg_eps * VEKNORM_CG(b, n);
    else if (linear_error_type == 2)
        eps = cg_eps * VEKNORM_CG(r, n);
    else if (linear_error_type == 3)
    {
        if ((r0norm = VEKNORM_CG(r, n)) > 1.)
            eps = cg_eps * r0norm;
        else
            eps = cg_eps;
    }
    else if (linear_error_type == 4)
        eps = cg_eps * (VEKNORM_CG(x, n));
    else if (linear_error_type == 5)
        eps = cg_eps * MMax(MMax(x0norm, b0norm), r0norm);
    else if (linear_error_type == 6)
    {
        MXMatVek(x, s);
        b0norm = VEKNORM_BICGSTAB(b, n);
        x0norm = VEKNORM_BICGSTAB(s, n);
        r0norm = VEKNORM_BICGSTAB(r, n);
        eps = cg_eps * MMax(MMax(x0norm, b0norm), r0norm);
    }
    s = (double*)Free(s);

    if (VEKNORM_CG(r, n) <= eps)
    {
        r = (double*)Free(r);
        /* Ggf. abschliessen der Vorkonditionierung */
        if (vorkond)
            MXVorkond(1, x, b);
        return 0;
    }
    rs = (double*)Malloc(n * sizeof(double));
    p = (double*)Malloc(n * sizeof(double));
    q = (double*)Malloc(n * sizeof(double));
    u = (double*)Malloc(n * sizeof(double));
    v = (double*)Malloc(n * sizeof(double));
    w = (double*)Malloc(n * sizeof(double));
    us = (double*)Malloc(n * sizeof(double));
    tmp = (double*)Malloc(n * sizeof(double));

    if
        VK_Modus(VK_iLDU) MXVorkond(2, r, r); /* Ra, 3/2000 */
    for (i = 0; i < n; i++)
        rs[i] = r[i];

    for (;;)
    {
        rho = MSkalarprodukt(rs, r, n);
        beta = rho / rho1;
        alpha = beta * beta;
        for (i = 0; i < n; i++)
        {
            u[i] = r[i] + beta * q[i];
            p[i] = u[i] + alpha * p[i] + beta * q[i];
        }
        MXMatVek(p, v);
        if
            VK_Modus(VK_iLDU) MXVorkond(2, v, v); /* Ra, 3/2000 */
        alpha = rho / MSkalarprodukt(rs, v, n);
        for (i = 0; i < n; i++)
        {
            q[i] = u[i] - alpha * v[i];
            us[i] = u[i] + q[i];
            x[i] += alpha * us[i];
        }
#ifdef SOLVER_SHOW_RESULTS
        printf("\n%ld %f %f %f %f %f", (long)k, x[(long)(n * .1)],
               x[(long)(n * .3)], x[(long)(n * .5)], x[(long)(n * .7)],
               x[(long)(n * .9)]);
#endif
        MXMatVek(us, tmp);
        if
            VK_Modus(VK_iLDU) MXVorkond(2, tmp, tmp); /* Ra, 3/2000 */
        for (i = 0; i < n; i++)
            r[i] -= alpha * tmp[i];

        k++;

        if (linear_error_type == 4)
            eps = cg_eps * (VEKNORM_QMRCGSTAB(x, n));
        if (VEKNORM_QMRCGSTAB(r, n) <= eps)
            break;

#ifdef SOLVER_SHOW_ERROR
        if (k % (int)MMax((cg_maxiter / 10), 1.) == 1)
        {
            DisplayMsg("      Iteration-Nr.: ");
            DisplayLong((long)k);
            DisplayMsg(", Fehler = ");
            DisplayDouble(VEKNORM_QMRCGSTAB(r, n) / eps, 4, 1);
            DisplayMsgLn(" ");
        }
#endif

        if (k >= max_iter)
            break;
        rho1 = rho;
    }

    r = (double*)Free(r);
    rs = (double*)Free(rs);
    p = (double*)Free(p);
    q = (double*)Free(q);
    u = (double*)Free(u);
    v = (double*)Free(v);
    w = (double*)Free(w);
    us = (double*)Free(us);
    tmp = (double*)Free(tmp);

    /* Ggf. abschliessen der Vorkonditionierung */
    if (vorkond)
        MXVorkond(1, x, b);

    return k;
}

/*************************************************************************
   ROCKFLOW - Modul: loeser5.c

   Aufgabe:
   Direkter Gleichungsloeser mit Speichertechnik aus 'matrix.h';
   durch Speichertechnik als direkter Loeser sehr langsam.

   Programmaenderungen:
   11/1995     MSR        Erste Version

*************************************************************************/

/*************************************************************************
   ROCKFLOW - Funktion: Gauss

   Aufgabe:
   Gleichungsloeser:
   LR-Faktorisierung von matrix nach dem Gausschen Algorithmus mit
   partieller Pivotisierung; Berechnung der Loesung aus LR-Faktorisierung
   Der Ergebnisvektor wird erst auf die rechte Seite gespeichert, und
   hinterher umkopiert. Hier kann noch Speicherplatz gespart werden.
   Die alten Inhalte von matrix und vecb werden zerstoert!
   Quelle: Schwetlick/Kretzschmar, Numerische Verfahren 1991

   Die Matrix wird durch die Speichertechnik implizit bereitgestellt.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E double *vecb: Rechte Seite Gleichungssystem
   X double *vecx: Ergebnisvektor, Speicher muss bereits reserviert sein
   E long g: Dimension des Gleichungssystems

   Ergebnis:
   0 (Anzahl der Iterationen)

   Programmaenderungen:
   11/1995     MSR        umgeschrieben auf Speichertechnik

*************************************************************************/
int SpGauss(double* vecb, double* vecx, long g)
{
    /* Matrizen sind in C zeilenweise abgespeichert:
       matrix[i][j] -> matrix[i*g+j] */
    static int* s;
    static double z, hilf;
    register int k, i, j, sk;
#ifdef TESTLOES
    DisplayMsgLn("SpGAUSS");
#endif
    s = (int*)Malloc(sizeof(int) * (g - 1));
    /* LR-Faktorisierung */
    for (k = 0; k < (g - 1); k++)
    {
        /* Pivotsuche */
        z = 0.0;
        sk = 0;
        for (i = k; i < g; i++)
        {
            hilf = fabs(MXGet(i, k)); /* matrix[i][k] */
            if (hilf > z)
            {
                z = hilf;
                sk = i;
            }
        }
        s[k] = sk;
        /* evtl. Zeilen vertauschen */
        if (sk > k)
            for (j = 0; j < g; j++)
            {
                z = MXGet(k, j);      /* matrix[k][j] */
                MXTrans(k, j, sk, j); /* matrix[k][j] = matrix[sk][j] */
                MXSet(sk, j, z);      /* matrix[sk][j] */
            }
        /* Berechnung der Eliminationskoeffizienten */
        for (i = (k + 1); i < g; i++)
            MXDiv(i, k, MXGet(k, k)); /* matrix[i][k], matrix[k][k] */
        /* Spaltenweise Berechnung der neuen Restmatrix */
        for (j = (k + 1); j < g; j++)
            for (i = (k + 1); i < g; i++)
                MXDec(i, j, (MXGet(i, k) * MXGet(k, j)));
        /* matrix[i][j], matrix[i][k], matrix[k][j] */
    }
    /* Loesung berechnen */
    /* vecb transformieren */
    for (k = 0; k < (g - 1); k++)
    {
        sk = s[k];
        if (sk > k)
        {
            z = vecb[k];
            vecb[k] = vecb[sk];
            vecb[sk] = z;
        }
    }
    for (k = 1; k < g; k++)
        for (j = 0; j < k; j++)
            vecb[k] -= (MXGet(k, j) * vecb[j]);
    /* matrix[k][j] */
    /* vecx berechnen */
    for (k = (g - 1); k >= 0; k--)
    {
        for (j = (k + 1); j < g; j++)
            vecb[k] -= (MXGet(k, j) * vecb[j]); /* matrix[k][j] */
        vecb[k] /= MXGet(k, k);                 /* matrix[k][k] */
    }
    /* Umspeichern des Ergebnisses von vecb nach vecx */
    for (k = 0; k < g; k++)
        vecx[k] = vecb[k];
    /* Speicher freigeben */
    s = (int*)Free(s);
    /* Ende */
    return 0;
}

/*************************************************************************
   ROCKFLOW - Modul: loeser6.c

   Aufgabe:
   Initialisierung und Aufruf der nichtlinearen Gleichungsloeser

   Programmaenderungen:
   11/1995     MSR        Erste Version
   6/1997     C.Thorenz  Berechnung des Abbruchfehlers EPS geaendert
   Speicherloch bei vorzeitigem Abbruch beseitigt
   01.07.1997  R.Kaiser   Korrekturen und Aenderungen aus dem aTM
   uebertragen
   08/1997     O.Kolditz  Einbau der nichtli. Loeser (PICARD) im Gas-Modell
   29.08.1997  A.Habbar   Eigenstaendige Funktion PICARD.
   09/1997     O.Kolditz  Optimierung des nichtlinearen Gleichungsloesers
   29.08.1997  A.Habbar   Weitere Optimierungen
   10/1997     C.Thorenz  Variables Abbruchkriterium fuer CG-Loeser

   letzte Aenderung: OK 29.11.97

*************************************************************************/

/*************************************************************************
   ROCKFLOW - Funktion: NonLinearSolve

   Aufgabe:
   Initilisierung und Aufruf der nichtlinearen Gleichungsloeser der
   Form :  A(x) * x = b

   Die Matrix A steht ueber die Speichertechnik implizit zur Verfuegung.
   Die benutzte Vektornorm wird in makros.h gesetzt.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E long   cas: 1 fuer Stroemung, 2 fuer Transport
   E double *b : Vektor b
   R double *x : Ergebnisvektor x;
   der Speicher muss bereits allokiert sein
   E long n    : Dimension des Gleichungssystems
   E void(*f)(): Zeiger auf eine Funktion, die das globale System wieder
   aufbaut.
   E long ind  : Index der Unbekannte in der internen Datenstruktur

   Ergebnis:
   Anzahl der Iterationen

   Programmaenderungen:
   ??/????     DAHMANE         Erste Version
   11/1995     MSR/DAHMANE     an ROCKFLOW und Speichertechnik angepasst
   6/1997     C.Thorenz  Berechnung des Abbruchfehlers EPS geaendert
   Speicherloch bei vorzeitigem Abbruch beseitigt
   10/1997     C.Thorenz  Variables Abbruchkriterium fuer CG-Loeser

*************************************************************************/
int NonLinearSolve(long cas, double* b, double* x, long n,
                   void (*f)(double*, double*, double), long ind)
{
    /* Variablen */
    static long iter;

#ifdef TESTLOES
    static char text[20];
#endif

    if (cas == 1)
    {
        /* Gleichungssystem-Parameter setzen */
        switch (loeser_flow)
        {
            case 1:
                LinearSolver = SpGauss;
                break;
            case 2:
#ifndef USE_MPI
                LinearSolver = SpBICGSTAB;
#endif
                break;
            case 3:
                LinearSolver = SpBICG;
                break;
            case 4:
                LinearSolver = SpQMRCGSTAB;
                break;
            case 5:
                LinearSolver = SpCG;
                break;
        }
        cg_maxiter = maxiter_flow;
        cg_eps = eps_flow;
        vorkond = vorkond_flow;
        linear_error_type = linear_error_type_flow;

        /* Nichtlinearer-Loeser-Parameter setzen */
        switch (nonlinear_method_flow)
        {
            case 1:
                NonlinearSolver = SpPICARD;
                break;
            case 2:
                NonlinearSolver = SpNEWTON;
                break;
        }
        nonlinear_maxiter = nonlinear_maxiter_flow;
        nonlinear_convergence_type = nonlinear_convergence_type_flow;
        nonlinear_abs_eps = nonlinear_abs_eps_flow;
        nonlinear_rel_eps = nonlinear_rel_eps_flow;
        nonlinear_assemble = nonlinear_assemble_flow;
        nonlinear_rel_cg_eps = nonlinear_rel_cg_eps_flow;
#ifdef TESTLOES6
        switch (nonlinear_method_flow)
        {
            case 1:
                strcpy(text, "SpPICARD");
                break;
            case 2:
                strcpy(text, "SpNEWTON");
                break;
        }
        DisplayMsgLn(text);
#endif
    }
    else if (cas == 2)
    {
        /* Gleichungssystem-Parameter setzen */
        switch (loeser_tran)
        {
            case 1:
                LinearSolver = SpGauss;
                break;
            case 2:
#ifndef USE_MPI
                LinearSolver = SpBICGSTAB;
#endif
                break;
            case 3:
                LinearSolver = SpBICG;
                break;
            case 4:
                LinearSolver = SpQMRCGSTAB;
                break;
            case 5:
                LinearSolver = SpCG;
                break;
        }
        cg_maxiter = maxiter_tran;
        cg_eps = eps_tran;
        vorkond = vorkond_tran;
        linear_error_type = linear_error_type_tran;

        /* Nichtlinearer-Loeser-Parameter setzen */
        switch (nonlinear_method_tran)
        {
            case 1:
                NonlinearSolver = SpPICARD;
                break;
            case 2:
                NonlinearSolver = SpNEWTON;
                break;
        }
        nonlinear_maxiter = nonlinear_maxiter_tran;
        nonlinear_convergence_type = nonlinear_convergence_type_tran;
        nonlinear_abs_eps = nonlinear_abs_eps_tran;
        nonlinear_rel_eps = nonlinear_rel_eps_tran;
        nonlinear_assemble = nonlinear_assemble_tran;
        nonlinear_rel_cg_eps = nonlinear_rel_cg_eps_tran;

#ifdef TESTLOES
        switch (nonlinear_method_tran)
        {
            case 1:
                strcpy(text, "SpPICARD");
                break;
            case 2:
                strcpy(text, "SpNEWTON");
                break;
        }
        DisplayMsgLn(text);
#endif
    }
    iter = NonlinearSolver(b, x, n, f, ind);

    if (iter == nonlinear_maxiter)
    {
        DisplayErrorMsg("Abbruch! Maximale Anzahl an Iterationen erreicht.");
        return 0;
    }
    return iter;
}

/*************************************************************************
   ROCKFLOW - Funktion: SpPICARD

   Aufgabe:
   Loesung des nichtlinearen Gleichungssystems der Form :  A(x) * x = b
   mit Hilfe der Picard-Iteration.

   Die Matrix A steht ueber die Speichertechnik implizit zur Verfuegung.
   Die benutzte Vektornorm wird in makros.h gesetzt.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E double *b : Vektor b
   R double *x : Ergebnisvektor x;
   der Speicher muss bereits allokiert sein
   E long n    : Dimension des Gleichungssystems
   E void(*f)(): Zeiger auf eine Funktion, die das globale System wieder
   aufbaut.
   E long ind  : Index der Unbekannte in der internen Datenstruktur

   Ergebnis:
   Anzahl der Iterationen

   Programmaenderungen:
   08/1997     O.KOLDITZ    Erste Version in Gas-Modell
   09/1997     MSR/AH       - an ROCKFLOW und Speichertechnik angepasst
   - Eigenstaendige Funktion PICARD.
   10/1997     C.Thorenz  Variables Abbruchkriterium fuer CG-Loeser

   letzte Aenderung: OK 29.11.97
 */

/**************************************************************************/
int SpPICARD(double* b, double* x, long n,
             void (*f)(double* b, double* x, double aktuelle_zeit), long ind)
/*int SpPICARD ( double *b, double *x, long n, void (*f)(), long ind ) */
{
    ind = ind;  // OK411
    /* Variablen */
    static double *r, *xs, *rs;
    int k = 0;
    static long iter;
    static double error;
    static double nonlinear_eps;

    static int iterations_old_timestep = -1;
    static double cg_eps_original; /* Hilfsspeicher fuer die urspruengliche
                                      Loesungsgenauigkeit */
    static int dirty; /* Kennzeichnet, dass die Iteration trotz unterschreiten
                         der Fehlerschranke nicht abgebrochen werden darf */
    static int rebuild_matrix; /* Kennzeichnet, dass die Systemmatrix neu
                                  aufgebaut werden muss */

#define GERMAN

#ifdef TESTLOES6
    DisplayMsgLn("SpPICARD");
#endif

    cg_eps_original = cg_eps;
    /* Fuer ersten Aufruf geringere Genauigkeit */
    if ((nonlinear_rel_cg_eps > 0) && (iterations_old_timestep == -1))
        cg_eps = nonlinear_rel_cg_eps;
    if ((nonlinear_rel_cg_eps > 0) && (iterations_old_timestep > 2))
        cg_eps = nonlinear_rel_eps;

    /****************************************************************************/
    /* Sichern der Loesung des letzten Iterationsschrittes
       Speicherreservierung fuer temporaere Felder */
    if (nonlinear_convergence_type == 1)
    {
        xs = (double*)Malloc(n * sizeof(double));
        MKopierVec(x, xs, n);
    }
    if (nonlinear_convergence_type == 2)
    {
        r = (double*)Malloc(n * sizeof(double));
        rs = (double*)Malloc(n * sizeof(double));
        MXResiduum(x, b, r);
        MKopierVec(r, rs, n);
    }
    if (nonlinear_convergence_type == 3)
    {
        xs = (double*)Malloc(n * sizeof(double));
        MKopierVec(x, xs, n);
    }
    /*  bs  = (double *) Malloc(n*sizeof(double));
       MKopierVec(b,bs,n);
     */
    rebuild_matrix = 0;

    /****************************************************************************/
    /* Loesen des Gleichungssystems und Fehlerberechnung */
    for (;;)
    {
        k++;
        dirty = 0;

        /*  M2ZeigMat("GLOB-MATRIX");
           MZeigVec(b,n,"GLOB-R.S"); */

        /* Neuaufbau des Gleichungssystems */
        if ((rebuild_matrix) || (!(k % nonlinear_assemble_flow)))
        {
#ifdef TESTLOES6
            DisplayMsgLn("Neuaufbau der Systemmatrix");
#endif
            f(b, x, aktuelle_zeit);
            rebuild_matrix = 0;
        }
        else
            dirty = 1;
        /* Loesen des Gleichungssystems */

        iter = LinearSolver(b, x, n);

        /* Speichern des Ergebnisses des Iterationsschritts in nval[ind] */
        // OK411        TransferNodeVals(x, ind);

        /* Fehlerberechnung */
        switch (nonlinear_convergence_type)
        {
            case 1:
                error = MVekDist(x, xs, n);
                nonlinear_eps =
                    nonlinear_abs_eps + nonlinear_rel_eps * VEKNORM_BICG(xs, n);
                break;
            case 2:
                MXResiduum(x, b, r);
                error = VEKNORM_BICG(r, n);
                nonlinear_eps =
                    nonlinear_abs_eps + nonlinear_rel_eps * VEKNORM_BICG(rs, n);
                break;
            case 3:
                error = MVekDist(x, xs, n);
                nonlinear_eps = nonlinear_rel_eps * VEKNORM_BICG(x, n);
                break;
        }

#ifdef GERMAN
        DisplayMsg("      Iteration-Nr.: ");
        DisplayLong((long)k);
        DisplayMsg(", GLS-Iter. = ");
        DisplayLong((long)iter);
        DisplayMsg(", Fehler/Abbruchf. = ");
        DisplayDouble(error / nonlinear_eps, 4, 1);
        DisplayMsgLn(" ");
#endif
#ifdef ENGLISH
        DisplayMsg("      Iteration-Nr.: ");
        DisplayLong((long)k);
        DisplayMsg(", Linear solver iterations = ");
        DisplayLong((long)iter);
        DisplayMsg(", error/errorcrit. = ");
        DisplayDouble(error / nonlinear_eps, 4, 1);
        DisplayMsgLn("");
#endif

        if (nonlinear_rel_cg_eps > 0)
        {
            if (iter == 0)
                nonlinear_rel_cg_eps = nonlinear_rel_cg_eps / 10.;
            if ((error / nonlinear_eps < 10.) && (cg_eps > cg_eps_original))
            {
                dirty = 1;
                cg_eps = cg_eps_original;
            }
            else
                cg_eps = (error / nonlinear_eps) * nonlinear_rel_eps *
                         nonlinear_rel_cg_eps;
            if (cg_eps < cg_eps_original)
                cg_eps = cg_eps_original;
        }
        /**************************************************************************/
        /* Fehlerabfrage */
        if (error <= nonlinear_eps || k >= nonlinear_maxiter)
        {
            if (dirty)
                rebuild_matrix = 1;
            else
                break;
        }
        if (nonlinear_convergence_type == 1)
            MKopierVec(x, xs, n);
        else if (nonlinear_convergence_type == 2)
            MKopierVec(r, rs, n);
        /**************************************************************************/

        /*    MKopierVec(bs,b,n); */
        if (nonlinear_convergence_type == 1)
            MKopierVec(xs, x, n);
    }
    /* Ende der Iterationsschleife
     **************************************************************************
     */

    iterations_old_timestep = k;

    /* Speicherfreigabe */
    if ((nonlinear_convergence_type == 1) || (nonlinear_convergence_type == 3))
        xs = (double*)Free(xs);
    else if (nonlinear_convergence_type == 2)
    {
        r = (double*)Free(r);
        rs = (double*)Free(rs);
    }
    /*  bs= Free(bs); */

    return k;
}

/*************************************************************************
   ROCKFLOW - Funktion: SpNEWTON

   Aufgabe:
   Loesung des nichtlinearen Gleichungssystems der Form :  A(x) * x = b
   mit Hilfe des modifizierten Newton-Raphson-Verfahrens.

   Die Matrix A steht ueber die Speichertechnik implizit zur Verfuegung.
   Die benutzte Vektornorm wird in makros.h gesetzt.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E double *b : Vektor b
   R double *x : Ergebnisvektor x;
   der Speicher muss bereits allokiert sein
   E long n    : Dimension des Gleichungssystems
   E void(*f)(): Zeiger auf eine Funktion, die das globale System wieder
   aufbaut.
   E long ind  : Index der Unbekannte in der internen Datenstruktur

   Ergebnis:
   Anzahl der Iterationen

   Programmaenderungen:
   ??/????     DAHMANE         Erste Version
   11/1995     MSR/DAHMANE     an ROCKFLOW und Speichertechnik angepasst
   6/1997     C.Thorenz  Berechnung des Abbruchfehlers EPS geaendert
   Speicherloch bei vorzeitigem Abbruch beseitigt
   10/1997      C.Thorenz  Einfuehrung eines relativen Fehlers fuer
   den CG-Loeser

   last modified: OK 28.06.1999 Dummy-Zuweisungen

*************************************************************************/
#ifdef ksdjfjdshashdhahsfgkfdg
int SpNEWTON(double* b, double* x, long n,
             void (*f)(double* b, double* x, double dummy), long ind)
#else
int SpNEWTON(double*, double*, long, void (*)(double*, double*, double), long)
#endif
/*int SpNEWTON ( double *b, double *x, long n, void (*f)(), long ind ) */
{
    // WW double ddummy;
    // WW long ldummy;
    // WW void (*g) (double *b, double *x, double dummy);

    DisplayMsgLn("Newton noch nicht implementiert !!! ");

    // WW ddummy = b[0];
    // WW ddummy = x[0];
    // WW ldummy = n;
    // WW ldummy = ind;
    // WW g = f;

    return nonlinear_maxiter;

#ifdef ksdjfjdshashdhahsfgkfdg
    /* Variablen */
    static double *r, *dx, *xs;
    static double eps;
    int k = 0, max_iter = 0;
    static float omega = 0.5;
    static double error;

#ifdef TESTLOES6
    DisplayMsgLn("SpNEWTON");
#endif

    k = 0;
    r = (double*)Malloc(n * sizeof(double));
    dx = (double*)Malloc(n * sizeof(double));
    xs = (double*)Malloc(n * sizeof(double));

#ifdef RELATIVE_EPS
    eps = VEKNORM_BICG(b, n);
    eps *= cg_eps;
#else
    eps = cg_eps;
#endif

#ifdef TESTLOES6
    DisplayMsg("eps = ");
    DisplayDouble(eps, 22, 20);
    DisplayMsgLn("");
#endif

    MXResiduum(x, b, r);
    if (VEKNORM_BICG(r, n) <= eps)
    {
        r = (double*)Free(r);
        dx = (double*)Free(r);
        xs = (double*)Free(xs);
        return k;
    }
    for (;;)
    {
        k++;
        MKopierVec(x, xs, n);
        LinearSolver(r, dx, n);
        MVekSum(x, omega, dx, n);
        TransferNodeVals(x, ind);
        MXResiduum(x, b, r);
        switch (nonlinear_convergence_type)
        {
            case 1:
                error = MVekDist(x, xs, n);
                break;
            case 2:
                error = VEKNORM_BICG(r, n);
                break;
        }
        if (error <= nonlinear_abs_eps)
            break;
        if (k >= nonlinear_maxiter)
            break;
        if ((k % nonlinear_assemble) == 0)
            f(b, x, ind);
    }

    r = (double*)Free(r);
    dx = (double*)Free(dx);
    xs = (double*)Free(xs);

    return k;
#endif
}

/**** Definitionen fuer Preconditioner entfernen */
#undef VK_Skalierung
#undef VK_Extraktion
#undef VK_iLDU
#undef VK_Modus
