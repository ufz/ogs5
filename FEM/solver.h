/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

 /***************************************************************************
   ROCKFLOW - Modul: solver.h

   Aufgabe:
   ROCKFLOW-Schnittstelle fuer alle Loeser

 **************************************************************************/

#ifndef solver_INC
#define solver_INC

#include "rf_num_new.h"
/*-----------------------------------------------------------------------*/
/* Objekt-Eigenschaften*/

/* Deklarationen fuer Schluesselwoerter
 #LINEAR_SOLVER_PRESSURE, #LINEAR_SOLVER_SATURATION, #LINEAR_SOLVER_TRANSPORT */
extern int loeser_flow, loeser_tran, loeser_temp, loeser_satu;
extern double gls_iter_theta;
extern int vorkond, vorkond_flow, vorkond_tran, vorkond_temp, vorkond_satu;
extern int speichertechnik_flow, speichertechnik_tran, speichertechnik_temp, speichertechnik_satu;
extern int linear_error_type,linear_error_type_flow,linear_error_type_tran,linear_error_type_temp,
           linear_error_type_satu;
extern double eps_flow, eps_tran, eps_temp, eps_satu, cg_eps;
extern int maxiter_flow, maxiter_tran, maxiter_temp, maxiter_satu, cg_maxiter;
extern int repeat_flow, repeat_tran, repeat_temp, repeat_satu, cg_repeat;

/* Funktionszeiger auf Loeser */
extern IntFuncDXDXL LoeserFlow, LoeserTran, LoeserTemp, LoeserSatu;
extern int sp2_start, sp2_inc;                    /* Speichertechnik-Feineinstellung */
extern double rel_eps;

/* Deklarationen fuer Schluesselwoerter #ITERATION_FLOW und #ITERATION_TRANSPORT */
/* Nichtlinear-Gleichungsloeser-Nummern */
extern int nonlinear_method_flow, nonlinear_method_tran;
/* max. Iterationen fuer Nichtlinear-Loeser */
extern int nonlinear_maxiter_flow, nonlinear_maxiter_tran;
/* Konvergenztyp fuer Nichtlinear-Loeser */
extern int nonlinear_convergence_type_flow, nonlinear_convergence_type_tran;
/* abs. Abbruchkriterium */
extern double nonlinear_abs_eps_flow, nonlinear_abs_eps_tran;
/* rel. Abbruchkriterium */
extern double nonlinear_rel_eps_flow, nonlinear_rel_eps_tran;
/* rel. Genauigkeit fuer CG-Loeser */
extern double nonlinear_rel_cg_eps,nonlinear_rel_cg_eps_flow, nonlinear_rel_cg_eps_tran;
/* Wiederaufbau des globalen Systems */
extern int nonlinear_assemble_flow, nonlinear_assemble_tran;

/* andere benoetigte globale variablen */
extern int nonlinear_method;                      /* Nichtlinear-Gleichungsloeser-Nummern */
extern int nonlinear_maxiter;                     /* max. Iterationen fuer Nichtlinear-Loeser */
extern int nonlinear_convergence_type;            /* Konvergenztyp fuer Nichtlinear-Loeser */
extern double nonlinear_abs_eps;                  /* abs. Abbruchkriterium */
extern double nonlinear_rel_eps;                  /* rel. Abbruchkriterium */
extern double nonlinear_rel_cg_eps;               /* rel. Genauigkeit fuer CG-Loeser */
extern int nonlinear_assemble;                    /* Wiederaufbau des globalen Systems */

/* Deklarationen fuer Schluesselwoerter #ITERATION_TIME_CONTROL */
extern int iteration_min_iter;
extern int iteration_max_iter;
extern double iteration_weight_plus;
extern double iteration_weight_up;
extern double iteration_weight_down;
extern double iteration_min_dt;
extern double iteration_max_dt;

/*-----------------------------------------------------------------------*/
/* Objekt-Methoden */
extern void InitSolverParameter(void);            /* ToDo raus */
extern void ConfigRenumberProperties(void);

/* Lineare Solver */
extern IntFuncDXDXL LinearSolver;                 /* Funktionszeiger auf linearen Loeser */
extern int SpRichardson ( double* b, double* x, long n );
/* Richardson-Loeser */
extern int SpJOR ( double* b, double* x, long n );
/* Jacobi-Loeser */
extern int SpSOR ( double* b, double* x, long n );
/* SOR-Loeser (Successive Over Relaxation) */
extern void Gauss ( double* matrix, double* vecb, double* vecx, int g );
/* Gauss-Loeser */
extern void LU_Decomposition (double* matrix, double* vecb, double* vecx, int g);
/* LU Dekomoposition */
extern int SpBICG ( double* b, double* x, long n );
/* BICG-Loeser */
extern int SpBICGSTAB_old ( double* b, double* x, long n );
#ifdef SX
int SpBICGSTAB(double* restrict b, double* restrict x, long n);
#else
int SpBICGSTAB(double* b, double* x,  long n);
#endif
/* BICGSTAB-Loeser */
extern int SpQMRCGSTAB ( double* b, double* x, long n );
/* QMRCGSTAB-Loeser */
extern int SpMGMRES ( double* b, double* x, long n );
/* MGMRES-Loeser */
extern int SpCG ( double* b, double* x, long n );
/* CG-Loeser */
extern int SpCGNR ( double* b, double* x, long n );
/* CGNR-Loeser */
extern int SpCGS ( double* b, double* x, long n );
/* CGS-Loeser */
extern int SpGauss ( double* vecb, double* vecx, long g );
/* Gauss-Loeser */
extern int SpAMG1R5( double* b, double* x, long n );
/* Algebraischer-Multigrid-Loeser AMG1R5 */
extern int SpUMF ( double* b, double* x, long n );
/* UMF-Loeser aus UMFPack */

/* Nicht-Lineare Solver */
extern IntFuncDXDXLVXL NonlinearSolver;           /* Funktionszeiger auf nichtlinearen Loeser */
extern int NonLinearSolve( long cas, double* b, double* x, long n, void (*)(double*,
                                                                            double*,
                                                                            double), long ind );
/* Initialisierung und Aufruf des NLGS */
extern int SpPICARD ( double* b, double* x, long n, void (*)(double*,double*,double), long );
/* PICARD-Loeser */
extern int SpNEWTON ( double* b, double* x, long n, void (*)(double*,double*,double), long );
/* mod-NEWTON-Loeser */
#endif
