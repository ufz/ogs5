/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Stiff_Bulirsch-Stoer.h"
#include "stdlib.h"
#include <iostream>


/* interne Deklarationen */
double *dvector(long nl, long nh);
void free_dvector(double *v, long nl, long nh);
void stifbs(double y[], double dydx[], int nv, double *xx, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
	void (*derivs)(double, double [], double [], int, long, double));
void odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1,
	double hmin, double *hnext, int *nok, int *nbad,
	void (*derivs)(double, double [], double [], int, long, double),
	bool (*stifbs)(double [], double [], int, double *, double, double, double [],
	double *, double *, void (*)(double, double [], double [], int, long, double)),
	bool (*rkqs)(double [], double [], int, double *, double , double , double [],
	double *, double *, void (*)(double , double [], double [], int, long, double), long),
	long, int);



double *dvector(long nl, long nh)
{
	// Removed
	std::cout << "Error: dvector not implemented!" << std::endl;
	return NULL;
}


void free_dvector(double *v, long nl, long nh)
{
	// Removed
	std::cout << "Error: free_dvector not implemented!" << std::endl;
}


/* Input: y=current_conc, dydx=their_derivs, xx=current_time, htry=suggested_stepsize     */
/* Ouput: y=updated_conc,  xx=end_time, hdid=achieved_stepsize , hnext= estimated_next_ss */
bool stifbs(double y[], double dydx[], int nv, double *xx, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
	void (*derivs)(double, double [], double [], int, long, double), long node)
{
	// Removed
	std::cout << "Error: stifbs not implemented!" << std::endl;
	return false;
}

bool rkqs(double y[], double dydx[], int n, double *x, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
	void (*derivs)(double , double [], double [], int, long, double), long node)
{
	// Removed
	std::cout << "Error: rkqs not implemented!" << std::endl;
	return false;
}

bool odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1,
	double hmin, double *nexth, int *nok, int *nbad,
	void (*derivs)(double, double [], double [], int, long, double),
	bool (*stifbs)(double [], double [], int, double *, double, double, double [],
	  double *, double *, void (*)(double, double [], double [], int, long, double), long),
  bool (*rkqs)(double [], double [], int, double *, double , double , double [],
	  double *, double *, void (*)(double , double [], double [], int, long, double),long),
  long node, int SolverType)
{
	// Removed
	std::cout << "Error: odeint not implemented!" << std::endl;
	return false;
}
