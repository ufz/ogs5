/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef STIFF_BULIRSCH_STOER
#define STIFF_BULIRSCH_STOER

/* Deklarationen */

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#define MAXSTEP 20000
#define TINY 1.0e-30
#define NR_END 0
#define FREE_ARG char*

/* externe Funktionen */
extern void derivs(double x, double y[], double dydx[], int n, long node);
extern void jacobn(double x, double y[], double dfdx[], double** dfdy, int n, long node);

#endif
