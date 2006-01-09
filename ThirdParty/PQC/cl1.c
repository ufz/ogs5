#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "phqalloc.h"
#include "output.h"
#include "phrqtype.h"
void 
cl1_space (int check, int n2d, int klm, int nklmd);
extern void zero_double(LDBLE *target, int n);
LDBLE *x_arg = NULL, *res_arg = NULL, *scratch = NULL;
int    x_arg_max = 0, res_arg_max = 0, scratch_max = 0;

int cl1 (int k, int l, int m, int n,
	 int nklmd, int n2d,
	 LDBLE * q,
	 int *kode, LDBLE toler,
	 int *iter, LDBLE * x, LDBLE * res, LDBLE * error,
	 LDBLE * cu, int *iu, int *s, int check);
static char const svnid[] = "$Id: cl1.c 4 2009-04-21 17:29:29Z delucia $";

extern void *free_check_null (void *ptr);
extern void malloc_error (void);

/* debug
#define DEBUG_CL1
#define CHECK_ERRORS
 */

int
cl1 (int k, int l, int m, int n,
     int nklmd, int n2d,
     LDBLE * q,
     int *kode, LDBLE toler,
     int *iter, LDBLE * x, LDBLE * res, LDBLE * error,
     LDBLE * cu, int *iu, int *s, int check)
{
  /* System generated locals */
  union double_or_int
  {
    int ival;
    LDBLE dval;
  } *q2;

  /* Local variables */
  static int nklm;
  static LDBLE xmin, xmax;
  static int iout, i, j;
  static LDBLE z;
  static int maxit, n1, n2;
  static LDBLE pivot;
  static int ia, ii, kk, in, nk, js;
  static LDBLE sn;
  static int iphase, kforce;
  static LDBLE zu, zv;
  static LDBLE tpivot;
  static int klm, jmn, nkl, jpn;
  static LDBLE cuv, sum;
  static int klm1;
  int q_dim, cu_dim;
  int kode_arg;
  LDBLE check_toler;
#ifdef CHECK_ERRORS
  extern char **col_name, **row_name;
  extern int *row_back, *col_back;
#endif
/* THIS SUBROUTINE USES A MODIFICATION OF THE SIMPLEX */
/* METHOD OF LINEAR PROGRAMMING TO CALCULATE AN L1 SOLUTION */
/* TO A K BY N SYSTEM OF LINEAR EQUATIONS */
/*             AX=B */
/* SUBJECT TO L LINEAR EQUALITY CONSTRAINTS */
/*             CX=D */
/* AND M LINEAR INEQUALITY CONSTRAINTS */
/*             EX.LE.F. */
/* DESCRIPTION OF PARAMETERS */
/* K      NUMBER OF ROWS OF THE MATRIX A (K.GE.1). */
/* L      NUMBER OF ROWS OF THE MATRIX C (L.GE.0). */
/* M      NUMBER OF ROWS OF THE MATRIX E (M.GE.0). */
/* N      NUMBER OF COLUMNS OF THE MATRICES A,C,E (N.GE.1). */
/* KLMD   SET TO AT LEAST K+L+M FOR ADJUSTABLE DIMENSIONS. */
/* KLM2D  SET TO AT LEAST K+L+M+2 FOR ADJUSTABLE DIMENSIONS. */
/* NKLMD  SET TO AT LEAST N+K+L+M FOR ADJUSTABLE DIMENSIONS. */
/* N2D    SET TO AT LEAST N+2 FOR ADJUSTABLE DIMENSIONS */
/* Q      TWO DIMENSIONAL REAL ARRAY WITH KLM2D ROWS AND */
/*        AT LEAST N2D COLUMNS. */
/*        ON ENTRY THE MATRICES A,C AND E, AND THE VECTORS */
/*        B,D AND F MUST BE STORED IN THE FIRST K+L+M ROWS */
/*        AND N+1 COLUMNS OF Q AS FOLLOWS */
/*             A B */
/*         Q = C D */
/*             E F */
/*        THESE VALUES ARE DESTROYED BY THE SUBROUTINE. */
/* KODE   A CODE USED ON ENTRY TO, AND EXIT */
/*        FROM, THE SUBROUTINE. */
/*        ON ENTRY, THIS SHOULD NORMALLY BE SET TO 0. */
/*        HOWEVER, IF CERTAIN NONNEGATIVITY CONSTRAINTS */
/*        ARE TO BE INCLUDED IMPLICITLY, RATHER THAN */
/*        EXPLICITLY IN THE CONSTRAINTS EX.LE.F, THEN KODE */
/*        SHOULD BE SET TO 1, AND THE NONNEGATIVITY */
/*        CONSTRAINTS INCLUDED IN THE ARRAYS X AND */
/*        RES (SEE BELOW). */
/*        ON EXIT, KODE HAS ONE OF THE */
/*        FOLLOWING VALUES */
/*             0- OPTIMAL SOLUTION FOUND, */
/*             1- NO FEASIBLE SOLUTION TO THE */
/*                CONSTRAINTS, */
/*             2- CALCULATIONS TERMINATED */
/*                PREMATURELY DUE TO ROUNDING ERRORS, */
/*             3- MAXIMUM NUMBER OF ITERATIONS REACHED. */
/* TOLER  A SMALL POSITIVE TOLERANCE. EMPIRICAL */
/*        EVIDENCE SUGGESTS TOLER = 10**(-D*2/3), */
/*        WHERE D REPRESENTS THE NUMBER OF DECIMAL */
/*        DIGITS OF ACCURACY AVAILABLE. ESSENTIALLY, */
/*        THE SUBROUTINE CANNOT DISTINGUISH BETWEEN ZERO */
/*        AND ANY QUANTITY WHOSE MAGNITUDE DOES NOT EXCEED */
/*        TOLER. IN PARTICULAR, IT WILL NOT PIVOT ON ANY */
/*        NUMBER WHOSE MAGNITUDE DOES NOT EXCEED TOLER. */
/* ITER   ON ENTRY ITER MUST CONTAIN AN UPPER BOUND ON */
/*        THE MAXIMUM NUMBER OF ITERATIONS ALLOWED. */
/*        A SUGGESTED VALUE IS 10*(K+L+M). ON EXIT ITER */
/*        GIVES THE NUMBER OF SIMPLEX ITERATIONS. */
/* X      ONE DIMENSIONAL REAL ARRAY OF SIZE AT LEAST N2D. */
/*        ON EXIT THIS ARRAY CONTAINS A */
/*        SOLUTION TO THE L1 PROBLEM. IF KODE=1 */
/*        ON ENTRY, THIS ARRAY IS ALSO USED TO INCLUDE */
/*        SIMPLE NONNEGATIVITY CONSTRAINTS ON THE */
/*        VARIABLES. THE VALUES -1, 0, OR 1 */
/*        FOR X(J) INDICATE THAT THE J-TH VARIABLE */
/*        IS RESTRICTED TO BE .LE.0, UNRESTRICTED, */
/*        OR .GE.0 RESPECTIVELY. */
/* RES    ONE DIMENSIONAL REAL ARRAY OF SIZE AT LEAST KLMD. */
/*        ON EXIT THIS CONTAINS THE RESIDUALS B-AX */
/*        IN THE FIRST K COMPONENTS, D-CX IN THE */
/*        NEXT L COMPONENTS (THESE WILL BE =0),AND */
/*        F-EX IN THE NEXT M COMPONENTS. IF KODE=1 ON */
/*        ENTRY, THIS ARRAY IS ALSO USED TO INCLUDE SIMPLE */
/*        NONNEGATIVITY CONSTRAINTS ON THE RESIDUALS */
/*        B-AX. THE VALUES -1, 0, OR 1 FOR RES(I) */
/*        INDICATE THAT THE I-TH RESIDUAL (1.LE.I.LE.K) IS */
/*        RESTRICTED TO BE .LE.0, UNRESTRICTED, OR .GE.0 */
/*        RESPECTIVELY. */
/* ERROR  ON EXIT, THIS GIVES THE MINIMUM SUM OF */
/*        ABSOLUTE VALUES OF THE RESIDUALS. */
/* CU     A TWO DIMENSIONAL REAL ARRAY WITH TWO ROWS AND */
/*        AT LEAST NKLMD COLUMNS USED FOR WORKSPACE. */
/* IU     A TWO DIMENSIONAL INTEGER ARRAY WITH TWO ROWS AND */
/*        AT LEAST NKLMD COLUMNS USED FOR WORKSPACE. */
/* S      INTEGER ARRAY OF SIZE AT LEAST KLMD, USED FOR */
/*        WORKSPACE. */
/*      DOUBLE PRECISION DBLE */
/*      REAL */

/* INITIALIZATION. */
  if (svnid == NULL)
    fprintf (stderr, " ");


  zv = 0;
  kode_arg = *kode;
  cl1_space(check, n2d, k+l+m, nklmd);

/* Parameter adjustments */
  q_dim = n2d;
  q2 = (union double_or_int *) q;
  cu_dim = nklmd;

/* Function Body */
  maxit = *iter;
  n1 = n + 1;
  n2 = n + 2;
  nk = n + k;
  nkl = nk + l;
  klm = k + l + m;
  klm1 = klm + 1;
  nklm = n + klm;
  kforce = 1;
  *iter = 0;
  js = 0;
  ia = -1;

/* SET UP LABELS IN Q. */
  for (j = 0; j < n; ++j)
  {
    q2[klm1 * q_dim + j].ival = j + 1;
  }
/* L10: */
  for (i = 0; i < klm; ++i)
  {
    q2[i * q_dim + n1].ival = n + i + 1;
    if (q2[i * q_dim + n].dval < 0.)
    {
      for (j = 0; j < n1; ++j)
      {
	q2[i * q_dim + j].dval = -q2[i * q_dim + j].dval;
      }
      q2[i * q_dim + n1].ival = -q2[i * q_dim + n1].ival;
/* L20: */
    }
  }
/* L30: */
/* SET UP PHASE 1 COSTS. */
  iphase = 2;
#ifdef DEBUG_CL1
  output_msg (OUTPUT_MESSAGE, "Set up phase 1 costs\n");
#endif
/* Zero first row of cu and iu */
  memcpy ((void *) &(cu[0]), (void *) &(scratch[0]),
	  (size_t) nklm * sizeof (LDBLE));
  for (j = 0; j < nklm; ++j)
  {
    iu[j] = 0;
  }
/* L40: */
#ifdef DEBUG_CL1
  output_msg (OUTPUT_MESSAGE, "L40\n");
#endif
  if (l != 0)
  {
    for (j = nk; j < nkl; ++j)
    {
      cu[j] = 1.;
      iu[j] = 1;
    }
/* L50: */
    iphase = 1;
  }

/* Copy first row of cu and iu to second row */
  memcpy ((void *) &(cu[cu_dim]), (void *) &(cu[0]),
	  (size_t) nklm * sizeof (LDBLE));
  memcpy ((void *) &(iu[cu_dim]), (void *) &(iu[0]),
	  (size_t) nklm * sizeof (int));

/* L60: */
#ifdef DEBUG_CL1
  output_msg (OUTPUT_MESSAGE, "L60\n");
#endif
  if (m != 0)
  {
    for (j = nkl; j < nklm; ++j)
    {
      cu[cu_dim + j] = 1.;
      iu[cu_dim + j] = 1;
      jmn = j - n;
      if (q2[jmn * q_dim + n1].ival < 0)
      {
	iphase = 1;
      }
    }
/* L70: */
  }
/* L80: */
#ifdef DEBUG_CL1
  output_msg (OUTPUT_MESSAGE, "L80\n");
#endif
  if (*kode != 0)
  {
    for (j = 0; j < n; ++j)
    {
      if (x[j] < 0.)
      {
/* L90: */
	cu[j] = 1.;
	iu[j] = 1;
      }
      else if (x[j] > 0.)
      {
	cu[cu_dim + j] = 1.;
	iu[cu_dim + j] = 1;
      }
    }
/* L110: */
#ifdef DEBUG_CL1
    output_msg (OUTPUT_MESSAGE, "L110\n");
#endif
    for (j = 0; j < k; ++j)
    {
      jpn = j + n;
      if (res[j] < 0.)
      {
/* L120: */
	cu[jpn] = 1.;
	iu[jpn] = 1;
	if (q2[j * q_dim + n1].ival > 0)
	{
	  iphase = 1;
	}
      }
      else if (res[j] > 0.)
      {
/* L130: */
	cu[cu_dim + jpn] = 1.;
	iu[cu_dim + jpn] = 1;
	if (q2[j * q_dim + n1].ival < 0)
	{
	  iphase = 1;
	}
      }
    }
/* L140: */
  }
/* L150: */
#ifdef DEBUG_CL1
  output_msg (OUTPUT_MESSAGE, "L150\n");
#endif
  if (iphase == 2)
  {
    goto L500;
  }
/* COMPUTE THE MARGINAL COSTS. */
L160:
#ifdef DEBUG_CL1
  output_msg (OUTPUT_MESSAGE, "L160\n");
#endif
  for (j = js; j < n1; ++j)
  {
    sum = 0.;
    for (i = 0; i < klm; ++i)
    {
      ii = q2[i * q_dim + n1].ival;
      if (ii < 0)
      {
	z = cu[cu_dim - ii - 1];
      }
      else
      {
	z = cu[ii - 1];
      }
      sum += q2[i * q_dim + j].dval * z;
    }
    q2[klm * q_dim + j].dval = sum;
  }
  for (j = js; j < n; ++j)
  {
    ii = q2[klm1 * q_dim + j].ival;
    if (ii < 0)
    {
      z = cu[cu_dim - ii - 1];
    }
    else
    {
      z = cu[ii - 1];
    }
    q2[klm * q_dim + j].dval -= z;
  }
/* DETERMINE THE VECTOR TO ENTER THE BASIS. */
L240:
#ifdef DEBUG_CL1
  output_msg (OUTPUT_MESSAGE, "L240, xmax %e\n", xmax);
#endif
  xmax = 0.;
  if (js >= n)
  {
    goto L490;			/* test for optimality */
  }
  for (j = js; j < n; ++j)
  {
    zu = q2[klm * q_dim + j].dval;
    ii = q2[klm1 * q_dim + j].ival;
    if (ii > 0)
    {
      zv = -zu - cu[ii - 1] - cu[cu_dim + ii - 1];
    }
    else
    {
      ii = -ii;
      zv = zu;
      zu = -zu - cu[ii - 1] - cu[cu_dim + ii - 1];
    }
/* L260 */
    if (kforce == 1 && ii > n)
    {
      continue;
    }
    if (iu[ii - 1] != 1 && zu > xmax)
    {
      xmax = zu;
      in = j;
    }
/* L270 */
    if (iu[cu_dim + ii - 1] != 1 && zv > xmax)
    {
      xmax = zv;
      in = j;
    }
  }
/* L280 */
#ifdef DEBUG_CL1
  output_msg (OUTPUT_MESSAGE, "L280 xmax %e, toler %e\n", xmax, toler);
#endif
  if (xmax <= toler)
  {
#ifdef DEBUG_CL1
    output_msg (OUTPUT_MESSAGE, "xmax before optimality test %e\n", xmax);
#endif
    goto L490;			/* test for optimality */
  }
  if (q2[klm * q_dim + in].dval != xmax)
  {
    for (i = 0; i < klm1; ++i)
    {
      q2[i * q_dim + in].dval = -q2[i * q_dim + in].dval;
    }
    q2[klm1 * q_dim + in].ival = -q2[klm1 * q_dim + in].ival;
/* L290: */
    q2[klm * q_dim + in].dval = xmax;
  }
/* DETERMINE THE VECTOR TO LEAVE THE BASIS. */
  if (iphase != 1 && ia != -1)
  {
    xmax = 0.;
/* find maximum absolute value in column "in" */
    for (i = 0; i <= ia; ++i)
    {
      z = fabs (q2[i * q_dim + in].dval);
      if (z > xmax)
      {
	xmax = z;
	iout = i;
      }
    }
/* L310: */
#ifdef DEBUG_CL1
    output_msg (OUTPUT_MESSAGE, "L310, xmax %e\n", xmax);
#endif
/* switch row ia with row iout, use memcpy */
    if (xmax > toler)
    {
      memcpy ((void *) &(scratch[0]), (void *) &(q2[ia * q_dim]),
	      (size_t) n2 * sizeof (LDBLE));
      memcpy ((void *) &(q2[ia * q_dim]), (void *) &(q2[iout * q_dim]),
	      (size_t) n2 * sizeof (LDBLE));
      memcpy ((void *) &(q2[iout * q_dim]), (void *) &(scratch[0]),
	      (size_t) n2 * sizeof (LDBLE));
/* L320: */
/* set pivot to row ia, column in */
      iout = ia;
      --ia;
      pivot = q2[iout * q_dim + in].dval;
      goto L420;		/* Gauss Jordan */
    }
  }
/* L330: */
#ifdef DEBUG_CL1
  output_msg (OUTPUT_MESSAGE, "L330, xmax %e\n", xmax);
#endif
  kk = -1;
/* divide column n1 by positive value in column "in" greater than toler */
  for (i = 0; i < klm; ++i)
  {
    z = q2[i * q_dim + in].dval;
    if (z > toler)
    {
      ++kk;
      res[kk] = q2[i * q_dim + n].dval / z;
      s[kk] = i;
    }
  }
/* L340: */
#ifdef DEBUG_CL1
  if (kk < 0)
  {
    output_msg (OUTPUT_MESSAGE, "kode = 2 in loop 340.\n");
  }
#endif
L350:
#ifdef DEBUG_CL1
  output_msg (OUTPUT_MESSAGE, "L350, xmax %e\n", xmax);
#endif
  if (kk < 0)
  {
/* no positive value found in L340 or bypass intermediate verticies */
    *kode = 2;
    goto L590;
  }
/* L360: */
#ifdef DEBUG_CL1
  output_msg (OUTPUT_MESSAGE, "L360, xmax %e\n", xmax);
#endif
/* find minimum residual */
  xmin = res[0];
  iout = s[0];
  j = 0;
  if (kk != 0)
  {
    for (i = 1; i <= kk; ++i)
    {
      if (res[i] < xmin)
      {
	j = i;
	xmin = res[i];
	iout = s[i];
      }
    }
/* L370: */
/* put kk in position j */
    res[j] = res[kk];
    s[j] = s[kk];
  }
/* L380: */
#ifdef DEBUG_CL1
  output_msg (OUTPUT_MESSAGE, "L380 iout %d, xmin %e, xmax %e\n", iout, xmin,
	      xmax);
#endif
  --kk;
  pivot = q2[iout * q_dim + in].dval;
  ii = q2[iout * q_dim + n1].ival;
  if (iphase != 1)
  {
    if (ii < 0)
    {
/* L390: */
      if (iu[-ii - 1] == 1)
      {
	goto L420;
      }
    }
    else
    {
      if (iu[cu_dim + ii - 1] == 1)
      {
	goto L420;
      }
    }
  }
/* L400: */
#ifdef DEBUG_CL1
  output_msg (OUTPUT_MESSAGE, "L400\n");
#endif
  ii = abs (ii);
  cuv = cu[ii - 1] + cu[cu_dim + ii - 1];
  if (q2[klm * q_dim + in].dval - pivot * cuv > toler)
  {

/* BYPASS INTERMEDIATE VERTICES. */
    for (j = js; j < n1; ++j)
    {
      z = q2[iout * q_dim + j].dval;
      q2[klm * q_dim + j].dval -= z * cuv;
      q2[iout * q_dim + j].dval = -z;
    }
/* L410: */
    q2[iout * q_dim + n1].ival = -q2[iout * q_dim + n1].ival;
    goto L350;
  }
/* GAUSS-JORDAN ELIMINATION. */
L420:
#ifdef DEBUG_CL1
  output_msg (OUTPUT_MESSAGE, "Gauss Jordon %d\n", *iter);
#endif
  if (*iter >= maxit)
  {
    *kode = 3;
    goto L590;
  }
/* L430: */
#ifdef DEBUG_CL1
  output_msg (OUTPUT_MESSAGE, "L430\n");
#endif
  ++(*iter);
  for (j = js; j < n1; ++j)
  {
    if (j != in)
    {
      q2[iout * q_dim + j].dval /= pivot;
    }
  }
/* L440: */
  for (j = js; j < n1; ++j)
  {
    if (j != in)
    {
      z = -q2[iout * q_dim + j].dval;
      for (i = 0; i < klm1; ++i)
      {
	if (i != iout)
	{
	  q2[i * q_dim + j].dval += z * q2[i * q_dim + in].dval;
	}
      }
/* L450: */
    }
  }
/* L460: */
#ifdef DEBUG_CL1
  output_msg (OUTPUT_MESSAGE, "L460\n");
#endif
  tpivot = -pivot;
  for (i = 0; i < klm1; ++i)
  {
    if (i != iout)
    {
      q2[i * q_dim + in].dval /= tpivot;
    }
  }
/* L470: */
  q2[iout * q_dim + in].dval = 1. / pivot;
  ii = q2[iout * q_dim + n1].ival;
  q2[iout * q_dim + n1].ival = q2[klm1 * q_dim + in].ival;
  q2[klm1 * q_dim + in].ival = ii;
  ii = abs (ii);
  if (iu[ii - 1] == 0 || iu[cu_dim + ii - 1] == 0)
  {
    goto L240;
  }
/* switch column */
  for (i = 0; i < klm1; ++i)
  {
    z = q2[i * q_dim + in].dval;
    q2[i * q_dim + in].dval = q2[i * q_dim + js].dval;
    q2[i * q_dim + js].dval = z;
  }
  i = q2[klm1 * q_dim + in].ival;
  q2[klm1 * q_dim + in].ival = q2[klm1 * q_dim + js].ival;
  q2[klm1 * q_dim + js].ival = i;
/* L480: */
  ++js;
  goto L240;
/* TEST FOR OPTIMALITY. */
L490:
#ifdef DEBUG_CL1
  output_msg (OUTPUT_MESSAGE, "L490\n");
#endif
  if (kforce == 0)
  {
    if (iphase == 1)
    {
      if (q2[klm * q_dim + n].dval <= toler)
      {
	goto L500;
      }
#ifdef DEBUG_CL1
      output_msg (OUTPUT_MESSAGE, "q2[klm1-1, n1-1] > *toler. %e\n",
		  q2[(klm1 - 1) * q_dim + n1 - 1].dval);
#endif
      *kode = 1;
      goto L590;
    }
    *kode = 0;
    goto L590;
  }
  if (iphase != 1 || q2[klm * q_dim + n].dval > toler)
  {
    kforce = 0;
    goto L240;
  }
/* SET UP PHASE 2 COSTS. */
L500:
#ifdef DEBUG_CL1
  output_msg (OUTPUT_MESSAGE, "Set up phase 2 costs %d\n", *iter);
#endif
  iphase = 2;
  for (j = 0; j < nklm; ++j)
  {
    cu[j] = 0.;
  }
/* L510: */
  for (j = n; j < nk; ++j)
  {
    cu[j] = 1.;
  }
  memcpy ((void *) &(cu[cu_dim]), (void *) &(cu[0]),
	  (size_t) nklm * sizeof (LDBLE));
/* L520: */
  for (i = 0; i < klm; ++i)
  {
    ii = q2[i * q_dim + n1].ival;
    if (ii <= 0)
    {
      if (iu[cu_dim - ii - 1] == 0)
      {
	continue;
      }
      cu[cu_dim - ii - 1] = 0.;
    }
    else
    {
/* L530: */
      if (iu[ii - 1] == 0)
      {
	continue;
      }
      cu[ii - 1] = 0.;
    }
/* L540: */
    ++ia;
/* switch row */
    memcpy ((void *) &(scratch[0]), (void *) &(q2[ia * q_dim]),
	    (size_t) n2 * sizeof (LDBLE));
    memcpy ((void *) &(q2[ia * q_dim]), (void *) &(q2[i * q_dim]),
	    (size_t) n2 * sizeof (LDBLE));
    memcpy ((void *) &(q2[i * q_dim]), (void *) &(scratch[0]),
	    (size_t) n2 * sizeof (LDBLE));
/* L550: */
  }
/* L560: */
  goto L160;


/* PREPARE OUTPUT. */
L590:
#ifdef DEBUG_CL1
  output_msg (OUTPUT_MESSAGE, "L590\n");
#endif
  sum = 0.;
  for (j = 0; j < n; ++j)
  {
    x[j] = 0.;
  }
/* L600: */
  for (i = 0; i < klm; ++i)
  {
    res[i] = 0.;
  }
/* L610: */
  for (i = 0; i < klm; ++i)
  {
    ii = q2[i * q_dim + n1].ival;
    sn = 1.;
    if (ii < 0)
    {
      ii = -ii;
      sn = -1.;
    }
    if (ii <= n)
    {
/* L620: */
      x[ii - 1] = sn * q2[i * q_dim + n].dval;
    }
    else
    {
/* L630: */
      res[ii - n - 1] = sn * q2[i * q_dim + n].dval;
      if (ii >= n1 && ii <= nk)
      {
/*     *    DBLE(Q(I,N1)) */
	sum += q2[i * q_dim + n].dval;
      }
    }
  }
/* L640: */
#ifdef DEBUG_CL1
  output_msg (OUTPUT_MESSAGE, "L640\n");
#endif
  *error = sum;
  /*
   *  Check calculation
   */
  if ((check == 1) && (*kode == 0))
  {
    check_toler = 10. * toler;
    /*
     *  Check optimization constraints
     */
    if (kode_arg == 1)
    {
      for (i = 0; i < k; i++)
      {
	if (res_arg[i] < 0.0)
	{
	  if (res[i] > check_toler)
	  {
#ifdef CHECK_ERRORS
	    output_msg (OUTPUT_MESSAGE,
			"\tCL1: optimization constraint not satisfied row %d, res %s, constraint %f.\n",
			row_name[row_back[i]], res[i], res_arg[i]);
#endif
	    *kode = 1;
	  }
	}
	else if (res_arg[i] > 0.0)
	{
	  if (res[i] < -check_toler)
	  {
#ifdef CHECK_ERRORS
	    output_msg (OUTPUT_MESSAGE,
			"\tCL1: optimization constraint not satisfied row %s, res %e, constraint %f.\n",
			row_name[row_back[i]], res[i], res_arg[i]);
#endif
	    *kode = 1;
	  }
	}
      }
    }
    /*
     *  Check equalities
     */
    for (i = k; i < k + l; i++)
    {
      if (fabs (res[i]) > check_toler)
      {
#ifdef CHECK_ERRORS
	output_msg (OUTPUT_MESSAGE,
		    "\tCL1: equality constraint not satisfied row %s, res %e, tolerance %e.\n",
		    row_name[row_back[i]], res[i], check_toler);
#endif
	*kode = 1;
      }
    }
    /*
     *  Check inequalities
     */
    for (i = k + l; i < k + l + m; i++)
    {
      if (res[i] < -check_toler)
      {
#ifdef CHECK_ERRORS
	output_msg (OUTPUT_MESSAGE,
		    "\tCL1: inequality constraint not satisfied row %s, res %e, tolerance %e.\n",
		    row_name[row_back[i]], res[i], check_toler);
#endif
	*kode = 1;
      }
    }
    /*
     *   Check dissolution/precipitation constraints
     */
    if (kode_arg == 1)
    {
      for (i = 0; i < n; i++)
      {
	if (x_arg[i] < 0.0)
	{
	  if (x[i] > check_toler)
	  {
#ifdef CHECK_ERRORS
	    output_msg (OUTPUT_MESSAGE,
			"\tCL1: dis/pre constraint not satisfied column %s, x %e, constraint %f.\n",
			col_name[col_back[i]], x[i], x_arg[i]);
#endif
	    *kode = 1;
	  }
	}
	else if (x_arg[i] > 0.0)
	{
	  if (x[i] < -check_toler)
	  {
#ifdef CHECK_ERRORS
	    output_msg (OUTPUT_MESSAGE,
			"\tCL1: dis/pre constraint not satisfied column %s, x %e, constraint %f.\n",
			col_name[col_back[i]], x[i], x_arg[i]);
#endif
	    *kode = 1;
	  }
	}
      }
    }
    if (*kode == 1)
    {
      output_msg (OUTPUT_MESSAGE,
		  "\n\tCL1: Roundoff errors in optimization.\n\t     Try using -multiple_precision in INVERSE_MODELING\n");
    }
  }
  return 0;
}

void
cl1_space (int check, int n2d, int klm, int nklmd)
{
  if (check == 1)
  {
    if (x_arg == NULL)
    {
      x_arg = (LDBLE *) PHRQ_malloc ((size_t) (n2d * sizeof (LDBLE)));
    }
    else if (n2d > x_arg_max)
    {
      x_arg = (LDBLE *) PHRQ_realloc (x_arg, (size_t) (n2d * sizeof (LDBLE)));
      x_arg_max = n2d;
    }
    if (x_arg == NULL)
      malloc_error ();
    zero_double(x_arg, n2d);

    if (res_arg == NULL) 
    {
      res_arg = (LDBLE *) PHRQ_malloc ((size_t) ((klm) * sizeof (LDBLE)));
    }
    else if (klm > res_arg_max)
    {
      res_arg = (LDBLE *) PHRQ_realloc (res_arg, (size_t) ((klm) * sizeof (LDBLE)));
      res_arg_max = klm;
    }
    if (res_arg == NULL)
      malloc_error ();
    zero_double(res_arg, klm);
  }

/* Make scratch space */
  if (scratch == NULL) 
  {
    scratch = (LDBLE *) PHRQ_malloc ((size_t) nklmd * sizeof (LDBLE));
  }
  else if (nklmd > scratch_max)
  {
    scratch = (LDBLE *) PHRQ_realloc (scratch, (size_t) nklmd * sizeof (LDBLE));
    scratch_max = nklmd;
  }
  if (scratch == NULL)
    malloc_error ();
  zero_double(scratch, nklmd);
}
