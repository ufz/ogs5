/**
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************/
/* ROCKFLOW - Modul: prototyp.h
 */
/* Aufgabe:
   Funktionsprototypen
 */
/* Programmaenderungen:
   05/99   OK   Implementierung
 */
/**************************************************************************/

#ifndef prototyp_INC
#define prototyp_INC

#include <stdio.h>

/* Deklarationen der Funktionszeigerprototypen */
typedef void (*VoidFuncVoid)(void);
typedef void* (*VoidXFuncVoid)(void);
typedef void* (*VoidXFuncVoidX)(void*);
typedef void (*VoidFuncCharX)(char*);
typedef void (*VoidFuncInt)(int);
typedef void (*VoidFuncIntFileX)(int, FILE*);
typedef void (*VoidFuncIntFileXLX)(int, FILE*, long*);
typedef void (*VoidFuncILD)(int, long, double);
typedef void (*VoidFuncILID)(int, long, int, double);
typedef void (*VoidFuncILDX)(int, long, double*);
typedef void (*VoidFuncILDXDX)(int, long, double*, double*);
typedef void (*VoidFuncILDDDDDX)(int, long, double, double, double, double,
                                 double*);
typedef void (*VoidFuncLong)(long);
typedef void (*VoidFuncLD)(long, double);
typedef void (*VoidFuncLDI)(long, double, int);
typedef void (*VoidFuncLDX)(long, double*);
typedef void (*VoidFuncLDXD)(long, double*, double);
typedef void (*VoidFuncLDXDX)(long, double*, double*);
typedef void (*VoidFuncLDXDXD)(long, double*, double*, double);
typedef void (*VoidFuncLI)(long, int);
typedef void (*VoidFuncLID)(long, int, double);
typedef void (*VoidFuncLIDX)(long, int, double*);
typedef void (*VoidFuncLLL)(long, long, long);
typedef void (*VoidFuncLLX)(long, long*);
typedef void (*VoidFuncLX)(long*);
typedef void (*VoidFuncDouble)(double);
typedef void (*VoidFuncLXDXDXD)(long, double*, double*, double*);
typedef void (*VoidFuncLIID)(long, int, int, double);
typedef void (*VoidFuncLIIDX)(long, int, int, double*);
typedef int (*IntFuncLIDX)(long, int, double*);
typedef int (*IntFuncVoid)(void);
typedef int (*IntFuncLong)(long);
typedef int (*IntFuncFileX)(FILE*);
typedef int (*IntFuncDXDXL)(double*, double*, long);
typedef int (*IntFuncDXDXLVXL)(double*, double*, long,
                               void (*)(double*, double*, double), long);
typedef int (*IntFuncII)(int, int);
typedef int (*IntFuncIII)(int, int, int);
typedef long* (*LongXFuncLIX)(long, int*);
typedef long (*LongFuncVoid)(void);
typedef double (*DoubleFuncVoid)(void);
typedef double (*DoubleFuncLong)(long);
typedef double (*DoubleFuncLongDouble)(long, double);
typedef double (*DoubleFuncLDDD)(long, double, double, double);
typedef double (*DoubleFuncLDDDD)(long, double, double, double, double);
typedef double (*DoubleFuncLI)(long, int);
typedef double (*DoubleFuncLII)(long, int, int);
typedef double (*DoubleFuncLIII)(long, int, int, int);
typedef double (*DoubleFuncLIIII)(long, int, int, int, int);
typedef double (*DoubleFuncILD)(int, long, double);
/* UJ rf3217 */
typedef double (*DoubleFuncILDDD)(int, long, double, double, double);
typedef double (*DoubleFuncILDDDD)(int, long, double, double, double, double);
typedef double (*DoubleFuncILDDID)(int, long, double, double, double, double);
typedef double (*DoubleFuncILLD)(int, long, long, double);
typedef double (*DoubleFuncILLI)(int, long, long, int);
typedef double (*DoubleFuncIV)(int, ...);
typedef double* (*DoubleXFuncLong)(long);
typedef double* (*DoubleXFuncLI)(long, int);
typedef double* (*DoubleXFuncLDX)(long, double*);
typedef double* (*DoubleXFuncLII)(long, int, int);
typedef void (*VoidFuncDXCDX)(double*, const double*);
typedef void* (*VoidXFuncIVoidX)(int, void*);
typedef double (*DoubleFuncIL)(int, long);
#endif
