#define EXTERNAL extern
#include "global.h"
#include "phqalloc.h"
#include "output.h"
#include "phrqproto.h"

static char const svnid[] =
  "$Id: integrate.c 4 2009-04-21 17:29:29Z delucia $";

#define MAX_QUAD 20
#define K_POLY 5
static LDBLE g_function (LDBLE x_value);
static LDBLE midpnt (LDBLE x1, LDBLE x2, int n);
static void polint (LDBLE * xa, LDBLE * ya, int n, LDBLE xv, LDBLE * yv,
		    LDBLE * dy);
static LDBLE qromb_midpnt (LDBLE x1, LDBLE x2);
static LDBLE z, xd, alpha;
static struct surface_charge *surface_charge_ptr;

static LDBLE calc_psi_avg (LDBLE surf_chrg_eq);
static int calc_all_donnan_music (void);
static int calc_init_donnan_music (void);

/* ---------------------------------------------------------------------- */
int
calc_all_g (void)
/* ---------------------------------------------------------------------- */
{
  int i, j, k;
  int converge, converge1;
  int count_g, count_charge;
  LDBLE new_g, xd1;
  LDBLE epsilon;

  if (svnid == NULL)
    fprintf (stderr, " ");
  if (use.surface_ptr == NULL)
    return (OK);
/*
 *   calculate g for each surface
 */
#ifdef SKIP
  if (punch.high_precision == FALSE)
  {
    epsilon = 1e-8;
    G_TOL = 1e-9;
  }
  else
  {
    epsilon = 1.e-12;
    G_TOL = 1e-10;
  }
#endif
  epsilon = convergence_tolerance;
  if (convergence_tolerance >= 1e-8)
  {
    G_TOL = 1e-9;
  }
  else
  {
    G_TOL = 1e-10;
  }

  converge = TRUE;
  count_charge = 0;
  for (j = 0; j < count_unknowns; j++)
  {
    if (x[j]->type != SURFACE_CB)
      continue;
    if (debug_diffuse_layer == TRUE)
      output_msg (OUTPUT_MESSAGE, "Calc_all_g, X[%d]\n", j);
    surface_charge_ptr = x[j]->surface_charge;
    count_g = 1;
    x[j]->surface_charge->g[0].charge = 0.0;
    x[j]->surface_charge->g[0].g = 0.0;
    x[j]->surface_charge->g[0].dg = 0.0;
    xd = exp (-2 * x[j]->master[0]->s->la * LOG_10);
    /* alpha = 0.02935 @ 25;                (ee0RT/2)**1/2, (L/mol)**1/2 C / m**2 */
    /* 1000 J/kJ and 1000 L/m**3 */
    alpha =
      sqrt (EPSILON * EPSILON_ZERO * (R_KJ_DEG_MOL * 1000.0) * 1000.0 * tk_x *
	    0.5);
/*
 *   calculate g for given surface for each species
 */
    for (i = 0; i < count_s_x; i++)
    {
      if (s_x[i]->type > HPLUS)
	continue;
      for (k = 0; k < count_g; k++)
      {
	if (equal (x[j]->surface_charge->g[k].charge, s_x[i]->z, G_TOL) ==
	    TRUE)
	{
	  s_x[i]->diff_layer[count_charge].charge = x[j]->surface_charge;
	  s_x[i]->diff_layer[count_charge].count_g = k;
	  break;
	}
      }
      if (k < count_g)
	continue;

      if (x[j]->surface_charge->grams > 0.0)
      {
	z = s_x[i]->z;
	if ((use.surface_ptr->only_counter_ions == FALSE) ||
	    (((x[j]->master[0]->s->la > 0) && (z < 0))
	     || ((x[j]->master[0]->s->la < 0) && (z > 0))))
	{
	  if (xd > 0.1)
	  {
	    new_g = qromb_midpnt (1.0, xd);
	  }
	  else if (xd > 0.01)
	  {
	    new_g = qromb_midpnt (1.0, 0.1);
	    new_g += qromb_midpnt (0.1, xd);
	  }
	  else if (xd > 0.001)
	  {
	    new_g = qromb_midpnt (1.0, 0.1);
	    new_g += qromb_midpnt (0.1, 0.01);
	    new_g += qromb_midpnt (0.01, xd);
	  }
	  else if (xd > 0.0001)
	  {
	    new_g = qromb_midpnt (1.0, 0.1);
	    new_g += qromb_midpnt (0.1, 0.01);
	    new_g += qromb_midpnt (0.01, .001);
	    new_g += qromb_midpnt (0.001, xd);
	  }
	  else if (xd > 0.00001)
	  {
	    new_g = qromb_midpnt (1.0, 0.1);
	    new_g += qromb_midpnt (0.1, 0.01);
	    new_g += qromb_midpnt (0.01, .001);
	    new_g += qromb_midpnt (0.001, .0001);
	    new_g += qromb_midpnt (0.0001, xd);
	  }
	  else if (xd > 0.000001)
	  {
	    new_g = qromb_midpnt (1.0, 0.1);
	    new_g += qromb_midpnt (0.1, 0.01);
	    new_g += qromb_midpnt (0.01, .001);
	    new_g += qromb_midpnt (0.001, .0001);
	    new_g += qromb_midpnt (0.0001, .00001);
	    new_g += qromb_midpnt (0.00001, xd);
	  }
	  else if (xd > 0.0000001)
	  {
	    new_g = qromb_midpnt (1.0, 0.1);
	    new_g += qromb_midpnt (0.1, 0.01);
	    new_g += qromb_midpnt (0.01, .001);
	    new_g += qromb_midpnt (0.001, .0001);
	    new_g += qromb_midpnt (0.0001, .00001);
	    new_g += qromb_midpnt (0.00001, .000001);
	    new_g += qromb_midpnt (0.000001, xd);
	  }
	  else if (xd > 0.00000001)
	  {
	    new_g = qromb_midpnt (1.0, 0.1);
	    new_g += qromb_midpnt (0.1, 0.01);
	    new_g += qromb_midpnt (0.01, .001);
	    new_g += qromb_midpnt (0.001, .0001);
	    new_g += qromb_midpnt (0.0001, .00001);
	    new_g += qromb_midpnt (0.00001, .000001);
	    new_g += qromb_midpnt (0.000001, .0000001);
	    new_g += qromb_midpnt (0.0000001, xd);
	  }
	  else
	  {
	    new_g = qromb_midpnt (1.0, 0.1);
	    new_g += qromb_midpnt (0.1, 0.01);
	    new_g += qromb_midpnt (0.01, .001);
	    new_g += qromb_midpnt (0.001, .0001);
	    new_g += qromb_midpnt (0.0001, .00001);
	    new_g += qromb_midpnt (0.00001, .000001);
	    new_g += qromb_midpnt (0.000001, .0000001);
	    new_g += qromb_midpnt (0.0000001, .00000001);
	    new_g += qromb_midpnt (0.00000001, xd);
	  }
	}
	else
	{
	  new_g = 0;
	}
      }
      else
      {
	new_g = 0.0;
      }
      if ((use.surface_ptr->only_counter_ions == TRUE) && new_g < 0)
	new_g = 0;
      x[j]->surface_charge->g[count_g].charge = s_x[i]->z;
      converge1 = TRUE;
      if (fabs (new_g) >= 1.)
      {
	if (fabs ((new_g - x[j]->surface_charge->g[count_g].g) / new_g) >
	    epsilon)
	{
	  converge1 = FALSE;
	}
      }
      else
      {
	if (fabs (new_g - x[j]->surface_charge->g[count_g].g) > epsilon)
	{
	  converge1 = FALSE;
	}
      }
      if (converge1 == FALSE)
      {
	converge = FALSE;
	if (debug_diffuse_layer == TRUE)
	{
	  output_msg (OUTPUT_MESSAGE, "\t%12f\t%12.4e\t%12.4e\t%12.4e\n",
		      (double) x[j]->surface_charge->g[count_g].charge,
		      (double) x[j]->surface_charge->g[count_g].g,
		      (double) new_g,
		      (double) (new_g - x[j]->surface_charge->g[count_g].g));
	}
      }
      x[j]->surface_charge->g[count_g].g = new_g;
      if (new_g == 0)
      {
	x[j]->surface_charge->g[count_g].dg = 0;
      }
      else
      {
	if (x[j]->surface_charge->grams > 0.0)
	{
	  x[j]->surface_charge->g[count_g].dg =
	    surface_charge_ptr->grams * surface_charge_ptr->specific_area *
	    alpha * g_function (xd) / F_C_MOL;
	  x[j]->surface_charge->g[count_g].dg *=
	    -2. / (exp (x[j]->master[0]->s->la * LOG_10) *
		   exp (x[j]->master[0]->s->la * LOG_10));
	  if ((xd - 1) < 0.0)
	  {
	    x[j]->surface_charge->g[count_g].dg *= -1.0;
	  }
	  if (fabs (x[j]->surface_charge->g[count_g].dg) < 1e-8)
	  {
	    xd1 = exp (-2 * 1e-3 * LOG_10);


	    new_g = qromb_midpnt (1.0, xd1);
	    x[j]->surface_charge->g[count_g].dg = new_g / .001;
	  }
	}
	else
	{
	  x[j]->surface_charge->g[count_g].dg = 0.0;
	}
      }
      s_x[i]->diff_layer[count_charge].charge = x[j]->surface_charge;
      s_x[i]->diff_layer[count_charge].count_g = count_g;
      count_g++;

    }
    if (debug_diffuse_layer == TRUE)
    {
      output_msg (OUTPUT_MESSAGE,
		  "\nSurface component %d: charge,\tg,\tdg/dlny,\txd\n",
		  count_charge);
      for (i = 0; i < count_g; i++)
      {
	output_msg (OUTPUT_MESSAGE, "\t%12f\t%12.4e\t%12.4e\t%12.4e\n",
		    (double) x[j]->surface_charge->g[i].charge,
		    (double) x[j]->surface_charge->g[i].g,
		    (double) x[j]->surface_charge->g[i].dg, (double) xd);
      }
    }
    count_charge++;
  }
  return (converge);
}

/* ---------------------------------------------------------------------- */
LDBLE
g_function (LDBLE x_value)
/* ---------------------------------------------------------------------- */
{
  LDBLE sum, return_value, sum1;
  int i, j;
  LDBLE ln_x_value;

  if (equal (x_value, 1.0, G_TOL * 100) == TRUE)
    return (0.0);
  sum = 0.0;
  ln_x_value = log (x_value);
  for (j = 0; j < use.surface_ptr->charge[0].count_g; j++)
  {
    use.surface_ptr->charge[0].g[j].psi_to_z =
      exp (ln_x_value * use.surface_ptr->charge[0].g[j].charge) - 1.0;
  }
  for (i = 0; i < count_s_x; i++)
  {
    if (s_x[i]->type < H2O && s_x[i]->z != 0.0)
    {
      for (j = 0; j < use.surface_ptr->charge[0].count_g; j++)
      {
	if (use.surface_ptr->charge[0].g[j].charge == s_x[i]->z)
	{
	  sum += s_x[i]->moles * use.surface_ptr->charge[0].g[j].psi_to_z;
	  break;
	}
      }
    }
  }
  if (sum < 0.0)
  {
    sum = 0.0;
    sum1 = 0.0;
    output_msg (OUTPUT_MESSAGE, "Species\tmoles\tX**z-1\tsum\tsum charge\n");
    for (i = 0; i < count_s_x; i++)
    {
      if (s_x[i]->type < H2O && s_x[i]->z != 0.0)
      {
	sum += s_x[i]->moles * (pow (x_value, s_x[i]->z) - 1.0);
	sum1 += s_x[i]->moles * s_x[i]->z;
	output_msg (OUTPUT_MESSAGE, "%s\t%e\t%e\t%e\t%e\n", s_x[i]->name,
		    (double) s_x[i]->moles,
		    (double) (pow (x_value, (double) s_x[i]->z) - 1.0),
		    (double) sum, (double) sum1);
      }
    }
    sprintf (error_string, "Negative sum in g_function, %e\t%e.",
	     (double) sum, (double) x_value);
    error_msg (error_string, CONTINUE);
    sprintf (error_string,
	     "Solutions must be charge balanced, charge imbalance is %e\n",
	     (double) sum1);
    error_msg (error_string, STOP);
  }

  return_value =
    (exp (ln_x_value * z) -
     1) / sqrt ((x_value * x_value * mass_water_aq_x * sum));
  return (return_value);
}

/* ---------------------------------------------------------------------- */
void
polint (LDBLE * xa, LDBLE * ya, int n, LDBLE xv, LDBLE * yv, LDBLE * dy)
/* ---------------------------------------------------------------------- */
{
  int i, m, ns;
  LDBLE den, dif, dift, ho, hp, w;
  LDBLE *c, *d;

  ns = 1;
  dif = fabs (xv - xa[1]);
/*
 *   Malloc work space
 */
  c = (double *) PHRQ_malloc ((size_t) (n + 1) * sizeof (LDBLE));
  if (c == NULL)
    malloc_error ();
  d = (double *) PHRQ_malloc ((size_t) (n + 1) * sizeof (LDBLE));
  if (d == NULL)
    malloc_error ();



  for (i = 1; i <= n; i++)
  {
    dift = fabs (xv - xa[i]);
    if (dift < dif)
    {
      ns = i;
      dif = dift;
    }
    c[i] = ya[i];
    d[i] = ya[i];
  }

  *yv = ya[ns--];
  for (m = 1; m < n; m++)
  {
    for (i = 1; i <= n - m; i++)
    {
      ho = xa[i] - xv;
      hp = xa[i + m] - xv;
      w = c[i + 1] - d[i];
      if ((den = ho - hp) == 0.0)
      {
	error_msg ("In subroutine polint.", STOP);
      }
      den = w / den;
      d[i] = hp * den;
      c[i] = ho * den;
    }
    if (2 * ns < (n - m))
    {
      *dy = c[ns + 1];
    }
    else
    {
      *dy = d[ns--];
    }
    *yv += *dy;

/*		*yv += (*dy = (2 * ns < (n-m) ? c[ns+1] : d[ns--])); */
  }
  c = (double *) free_check_null (c);
  d = (double *) free_check_null (d);
  return;
}

/* ---------------------------------------------------------------------- */
LDBLE
midpnt (LDBLE x1, LDBLE x2, int n)
/* ---------------------------------------------------------------------- */
{
  LDBLE xv, tnm, sum, del, ddel;
  static LDBLE sv;
  int it, j;

  if (n == 1)
  {
    sv = (x2 - x1) * g_function (0.5 * (x1 + x2));
    return (sv);
  }
  else
  {
    for (it = 1, j = 1; j < n - 1; j++)
      it *= 3;
    tnm = (LDBLE) it;
    del = (x2 - x1) / (3 * tnm);
    ddel = del + del;
    xv = x1 + 0.5 * del;
    sum = 0.0;
    for (j = 1; j <= it; j++)
    {
#if defined(PHREEQCI_GUI)
      if (WaitForSingleObject (g_hKill /*g_eventKill */ , 0) == WAIT_OBJECT_0)
      {
	error_msg ("Execution canceled by user.", CONTINUE);
	RaiseException (USER_CANCELED_RUN, 0, 0, NULL);
      }
#endif
      sum += g_function (xv);
      xv += ddel;
      sum += g_function (xv);
      xv += del;
    }
    sv = (sv + (x2 - x1) * sum / tnm) / 3.0;
    return sv;
  }
}

/* ---------------------------------------------------------------------- */
LDBLE
qromb_midpnt (LDBLE x1, LDBLE x2)
/* ---------------------------------------------------------------------- */
{
  LDBLE ss, dss;
  LDBLE sv[MAX_QUAD + 2], h[MAX_QUAD + 2];
  int j;

  h[0] = 1.0;
  sv[0] = midpnt (x1, x2, 1);
  for (j = 1; j < MAX_QUAD; j++)
  {
    sv[j] = midpnt (x1, x2, j + 1);
    h[j] = h[j - 1] / 9.0;

    if (fabs (sv[j] - sv[j - 1]) <= G_TOL * fabs (sv[j]))
    {
      sv[j] *= surface_charge_ptr->grams * surface_charge_ptr->specific_area * alpha / F_C_MOL;	/* (ee0RT/2)**1/2, (L/mol)**1/2 C / m**2 */
      if ((x2 - 1) < 0.0)
	sv[j] *= -1.0;
      if (debug_diffuse_layer == TRUE)
      {
	output_msg (OUTPUT_MESSAGE, "Iterations in qromb_midpnt: %d\n", j);
      }
      return (sv[j]);
    }

    if (j >= K_POLY - 1)
    {
      polint (&h[j - K_POLY], &sv[j - K_POLY], K_POLY, 0.0, &ss, &dss);
      if (fabs (dss) <= G_TOL * fabs (ss) || fabs (dss) < G_TOL)
      {
	ss *= surface_charge_ptr->grams * surface_charge_ptr->specific_area * alpha / F_C_MOL;	/* (ee0RT/2)**1/2, (L/mol)**1/2 C / m**2 */
	if ((x2 - 1) < 0.0)
	  ss *= -1.0;
	if (debug_diffuse_layer == TRUE)
	{
	  output_msg (OUTPUT_MESSAGE, "Iterations in qromb_midpnt: %d\n", j);
	}
	return (ss);
      }
    }

  }
  sprintf (error_string,
	   "\nToo many iterations integrating diffuse layer.\n");
  error_msg (error_string, STOP);
  return (-999.9);
}

/* ---------------------------------------------------------------------- */
int
calc_init_g (void)
/* ---------------------------------------------------------------------- */
{
  int i, j, k;
  int count_g, count_charge;

  if (use.surface_ptr == NULL)
    return (OK);

/*
 *   calculate g for each surface
 */
  count_charge = 0;
  for (j = 0; j < count_unknowns; j++)
  {
    if (x[j]->type != SURFACE_CB)
      continue;
    surface_charge_ptr = x[j]->surface_charge;
    count_g = 0;
    if (x[j]->surface_charge->g != NULL)
    {
      count_g = x[j]->surface_charge->count_g;
    }
    if (count_g == 0)
    {
      x[j]->surface_charge->g =
	(struct surface_diff_layer *) PHRQ_malloc ((size_t)
						   sizeof (struct
							   surface_diff_layer));
      if (x[j]->surface_charge->g == NULL)
	malloc_error ();
      x[j]->surface_charge->g[0].charge = 0.0;
      x[j]->surface_charge->g[0].g = 0.0;
      x[j]->surface_charge->g[0].dg = 0.0;
      xd = exp (-2 * x[j]->master[0]->s->la * LOG_10);
      /* alpha = 0.02935 @ 25;                (ee0RT/2)**1/2, (L/mol)**1/2 C / m**2 */
      /*  second 1000 is liters/m**3 */
      alpha =
	sqrt (EPSILON * EPSILON_ZERO * (R_KJ_DEG_MOL * 1000.0) * 1000.0 *
	      tk_x * 0.5);
    }
/*
 *   calculate g for given surface for each species
 */
    count_g = 1;
    for (i = 0; i < count_s_x; i++)
    {
      if (s_x[i]->type > HPLUS)
	continue;
      for (k = 0; k < count_g; k++)
      {
	if (equal (x[j]->surface_charge->g[k].charge, s_x[i]->z, G_TOL) ==
	    TRUE)
	{
	  s_x[i]->diff_layer[count_charge].charge = x[j]->surface_charge;
	  s_x[i]->diff_layer[count_charge].count_g = k;
	  s_x[i]->diff_layer[count_charge].g_moles = 0.0;
	  s_x[i]->diff_layer[count_charge].dg_g_moles = 0.0;
	  break;
	}
      }
      if (k >= count_g)
      {

	/* malloc space to save g for charge */
	x[j]->surface_charge->g =
	  (struct surface_diff_layer *) PHRQ_realloc (x[j]->surface_charge->g,
						      (size_t) (count_g +
								1) *
						      sizeof (struct
							      surface_diff_layer));
	if (x[j]->surface_charge->g == NULL)
	  malloc_error ();

	/* save g for charge */
	x[j]->surface_charge->g[count_g].charge = s_x[i]->z;
	if (x[j]->surface_charge->grams > 0.0)
	{
	  x[j]->surface_charge->g[count_g].g =
	    2 * alpha * sqrt (mu_x) * (pow (xd, s_x[i]->z / 2.0) -
				       1) * surface_charge_ptr->grams *
	    surface_charge_ptr->specific_area / F_C_MOL;
	  x[j]->surface_charge->g[count_g].dg = -s_x[i]->z;
	  if ((use.surface_ptr->only_counter_ions == TRUE) &&
	      x[j]->surface_charge->g[count_g].g < 0)
	  {
	    x[j]->surface_charge->g[count_g].g = 0;
	    x[j]->surface_charge->g[count_g].dg = 0;
	  }
	}
	else
	{
	  x[j]->surface_charge->g[count_g].g = 0.0;
	  x[j]->surface_charge->g[count_g].dg = -s_x[i]->z;
	}
	/* save g for species */
	s_x[i]->diff_layer[count_charge].charge = x[j]->surface_charge;
	s_x[i]->diff_layer[count_charge].count_g = count_g;
	s_x[i]->diff_layer[count_charge].g_moles = 0.0;
	s_x[i]->diff_layer[count_charge].dg_g_moles = 0.0;
	count_g++;
      }
    }
    if (debug_diffuse_layer == TRUE)
    {
      output_msg (OUTPUT_MESSAGE, "\nSurface component %d: charge,\tg,\tdg\n",
		  count_charge);
      for (i = 0; i < count_g; i++)
      {
	output_msg (OUTPUT_MESSAGE, "\t%12f\t%12.4e\t%12.4e\n",
		    (double) x[j]->surface_charge->g[i].charge,
		    (double) x[j]->surface_charge->g[i].g,
		    (double) x[j]->surface_charge->g[i].dg);
      }
    }
    count_charge++;
    x[j]->surface_charge->count_g = count_g;
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
initial_surface_water (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   In initial surface calculation, need to calculate
 *   mass of water in diffuse layer.
 *   diffuse layer water + aqueous solution water = bulk water.
 *   Ionic strength is fixed, so diffuse-layer water will not change
 */
  int i;
  LDBLE debye_length, b, r, rd, ddl_limit, rd_limit, fraction, sum_surfs, s;
  double damp_aq;

/*
 *   Debye  length = 1/k = sqrt[eta*eta_zero*R*T/(2*F**2*mu_x*1000)], Dzombak and Morel, p 36
 *
 *   1000 converts kJ to J; 1000 converts Liters to meter**3; debye_length is in meters.
 */
  debye_length = (EPSILON * EPSILON_ZERO * R_KJ_DEG_MOL * 1000.0 * tk_x)
    / (2. * F_C_MOL * F_C_MOL * mu_x * 1000.);
  debye_length = sqrt (debye_length);

  /*  ddl is at most the fraction ddl_limit of bulk water */
  ddl_limit = use.surface_ptr->DDL_limit;

/*
 *   Loop through all surface components, calculate each H2O surface (diffuse layer),
 *   H2O aq, and H2O bulk (diffuse layers plus aqueous).
 */

  if (use.surface_ptr->debye_lengths > 0)
  {
    sum_surfs = 0.0;
    for (i = 0; i < count_unknowns; i++)
    {
      if (x[i]->type != SURFACE_CB)
        continue;
      sum_surfs += x[i]->surface_charge->specific_area * x[i]->surface_charge->grams;
    }

    rd = debye_length * use.surface_ptr->debye_lengths;
    use.surface_ptr->thickness = rd;

    if (state == INITIAL_SURFACE)
    {
      /* distribute water over DDL (rd) and free pore (r - rd) */
      /* find r: free pore (m3) = pi * (r - rd)^2 * L, where L = A / (2*pi*r),
         A = sum_surfs = sum of the surface areas */
      b = -2 * (rd + use.solution_ptr->mass_water / (1000 * sum_surfs));
      r = 0.5 * (-b + sqrt (b * b - 4 * rd * rd));
      /* DDL (m3) = pi * (r^2 - (r - rd)^2) * L */
      rd_limit = (1 - sqrt (1 - ddl_limit)) * r;
      /* rd should be smaller than r and the ddl limit */
      if (rd > rd_limit)
      {
	mass_water_surfaces_x = use.solution_ptr->mass_water * ddl_limit / (1 - ddl_limit);
	r = 0.002 * (mass_water_surfaces_x + use.solution_ptr->mass_water) / sum_surfs;
	rd_limit = (1 - sqrt (1 - ddl_limit)) * r;
	use.surface_ptr->thickness = rd = rd_limit;
      }
      else
	mass_water_surfaces_x =
	  (r * r / pow (r - rd, 2) - 1) * use.solution_ptr->mass_water;
      for (i = 0; i < count_unknowns; i++)
      {
	if (x[i]->type != SURFACE_CB)
	  continue;
        s = x[i]->surface_charge->specific_area * x[i]->surface_charge->grams;
	x[i]->surface_charge->mass_water = mass_water_surfaces_x * s / sum_surfs;
      }
    }
    else
    {
      r = 0.002 * mass_water_bulk_x / sum_surfs;
      rd_limit = (1 - sqrt (1 - ddl_limit)) * r;
      if (rd > rd_limit)
      {
	use.surface_ptr->thickness = rd = rd_limit;
	fraction = ddl_limit;
      }
      else
	fraction = 1 - pow (r - rd, 2) / (r * r);
      damp_aq = 1.0;
      if (g_iterations > 10)
	damp_aq = 0.2;
      else if (g_iterations > 5)
	damp_aq = 0.5;
      mass_water_surfaces_x = damp_aq * fraction * mass_water_bulk_x +
	(1 - damp_aq) * mass_water_surfaces_x;
      for (i = 0; i < count_unknowns; i++)
      {
	if (x[i]->type != SURFACE_CB)
	  continue;
        s = x[i]->surface_charge->specific_area * x[i]->surface_charge->grams;
	x[i]->surface_charge->mass_water = mass_water_surfaces_x * s / sum_surfs;
      }
    }
  }
  else
  {
    /* take constant thickness of, default 1e-8 m (100 Angstroms) */
    mass_water_surfaces_x = 0.0;
    for (i = 0; i < count_unknowns; i++)
    {
      if (x[i]->type != SURFACE_CB)
	continue;
      x[i]->surface_charge->mass_water = x[i]->surface_charge->specific_area *
	x[i]->surface_charge->grams * use.surface_ptr->thickness * 1000;
      mass_water_surfaces_x += x[i]->surface_charge->mass_water;
    }
  }

  if (use.surface_ptr->type == CD_MUSIC)
    mass_water_bulk_x = mass_water_aq_x + mass_water_surfaces_x;
  else
  {
    /*  for variable distribution of water over DDL and bulk... */
    if (state > INITIAL_SURFACE)
      mass_water_aq_x = mass_water_bulk_x - mass_water_surfaces_x;
    else
      mass_water_bulk_x = mass_water_aq_x + mass_water_surfaces_x;
  }

  return (OK);
}

/* ---------------------------------------------------------------------- */
int
sum_diffuse_layer (struct surface_charge *surface_charge_ptr1)
/* ---------------------------------------------------------------------- */
{
  int i, j, count_g;
  LDBLE mass_water_surface;
  LDBLE molality, moles_excess, moles_surface;

  if (use.surface_ptr == NULL)
    return (OK);
/*
 *   Find position of component in list of components
 */
  i = 0;

  for (j = 0; j < use.surface_ptr->count_charge; j++)
  {
    if (&(use.surface_ptr->charge[j]) == surface_charge_ptr1)
    {
      i = j;
      break;
    }
  }
  if (j >= use.surface_ptr->count_charge)
  {
    sprintf (error_string, "In sum_diffuse_layer, component not found, %s.",
	     surface_charge_ptr1->name);
    error_msg (error_string, STOP);
  }
/*
 *   Loop through all surface components, calculate each H2O surface (diffuse layer),
 *   H2O aq, and H2O bulk (diffuse layers plus aqueous).
 */
  count_elts = 0;
  paren_count = 0;
  mass_water_surface = surface_charge_ptr1->mass_water;
  for (j = 0; j < count_s_x; j++)
  {
    if (s_x[j]->type > HPLUS)
      continue;
    molality = under (s_x[j]->lm);
    count_g = s_x[j]->diff_layer[i].count_g;
#ifdef SKIP
    moles_excess = mass_water_bulk_x *
/*			s_x[j]->diff_layer[i].charge->g[count_g].g * molality; */
      surface_charge_ptr1->g[count_g].g * molality;
#endif
    moles_excess =
      mass_water_aq_x * molality * surface_charge_ptr1->g[count_g].g;

    moles_surface = mass_water_surface * molality + moles_excess;
/*
 *   Accumulate elements in diffuse layer
 */
    add_elt_list (s_x[j]->next_elt, moles_surface);
  }
  add_elt_list (s_h2o->next_elt, mass_water_surface / gfw_water);

  if (count_elts > 0)
  {
    qsort (elt_list, (size_t) count_elts,
	   (size_t) sizeof (struct elt_list), elt_list_compare);
    elt_list_combine ();
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
calc_all_donnan (void)
/* ---------------------------------------------------------------------- */
{
  int i, j, k;
  int count_g, count_charge, converge;
  char name[MAX_LENGTH];
  LDBLE new_g, f_psi, surf_chrg_eq, psi_avg, f_sinh, A_surf, ratio_aq;
  LDBLE new_g2, f_psi2, surf_chrg_eq2, psi_avg2, dif;
  LDBLE cz, cm, cp;

  if (use.surface_ptr == NULL)
    return (OK);
  if (use.surface_ptr->type == CD_MUSIC)
    return (calc_all_donnan_music ());
  f_sinh =
    sqrt (8000.0 * EPSILON * EPSILON_ZERO * (R_KJ_DEG_MOL * 1000.0) * tk_x *
	  mu_x);
  cz = cm = 1.0;
  cp = 1.0;
/*
 *   calculate g for each surface...
 */
  converge = TRUE;
  count_charge = 0;
  for (j = 0; j < count_unknowns; j++)
  {
    if (x[j]->type != SURFACE_CB)
      continue;
    surface_charge_ptr = x[j]->surface_charge;

    if (debug_diffuse_layer == TRUE)
      output_msg (OUTPUT_MESSAGE, "Calc_all_g, X[%d]\n", j);
/*
 *  sum eq of each charge number in solution...
 */
    count_g = x[j]->surface_charge->count_g;
    for (i = 0; i < count_g; i++)
    {
      charge_group[i].eq = 0.0;
    }
    for (i = 0; i < count_s_x; i++)
    {
      if (s_x[i]->type > HPLUS)
	continue;
      for (k = 0; k < count_g; k++)
      {
	if (equal (charge_group[k].z, s_x[i]->z, G_TOL) == TRUE)
	{
#ifdef SKIP
	  if (s_x[i]->z > 0)
	    cz = s_x[i]->z * pow (cp, s_x[i]->z);
	  else
	    cz = s_x[i]->z * pow (cm, -s_x[i]->z);
	  charge_group[k].eq += cz * s_x[i]->moles;
#endif
	  charge_group[k].eq += s_x[i]->z * s_x[i]->moles;
	  break;
	}
      }
    }
    /* find surface charge from potential... */
    A_surf =
      x[j]->surface_charge->specific_area * x[j]->surface_charge->grams;
    f_psi = x[j]->master[0]->s->la * LOG_10;
    surf_chrg_eq = A_surf * f_sinh * sinh (f_psi) / F_C_MOL;
    /* also for the derivative... */
    dif = 1e-5;
    f_psi2 = f_psi + dif;
    surf_chrg_eq2 = A_surf * f_sinh * sinh (f_psi2) / F_C_MOL;


    /* find psi_avg that matches surface charge... */
    psi_avg = calc_psi_avg (surf_chrg_eq);
    psi_avg2 = calc_psi_avg (surf_chrg_eq2);

    /*output_msg(OUTPUT_MESSAGE, "psi's  %e %e %e\n", f_psi, psi_avg, surf_chrg_eq); */

    /* fill in g's */
    ratio_aq = surface_charge_ptr->mass_water / mass_water_aq_x;

    for (k = 0; k < count_g; k++)
    {
      x[j]->surface_charge->g[k].charge = charge_group[k].z;
#ifdef SKIP
      if (charge_group[k].z > 0)
	cz = pow (cp, charge_group[k].z);
      else
	cz = pow (cm, -charge_group[k].z);
      new_g = cz * (exp (-charge_group[k].z * psi_avg)) - 1;
      if (new_g < -ratio_aq)
	new_g = -ratio_aq + G_TOL * 1e-5;
      new_g2 = cz * (exp (-charge_group[k].z * psi_avg2)) - 1;
      if (new_g2 < -ratio_aq)
	new_g2 = -ratio_aq + G_TOL * 1e-5;
#endif
      new_g = ratio_aq * (exp (-charge_group[k].z * psi_avg) - 1);
      if (use.surface_ptr->only_counter_ions &&
	  ((surf_chrg_eq < 0 && charge_group[k].z < 0)
	   || (surf_chrg_eq > 0 && charge_group[k].z > 0)))
	new_g = -ratio_aq;
      if (new_g <= -ratio_aq)
	new_g = -ratio_aq + G_TOL * 1e-3;
      new_g2 = ratio_aq * (exp (-charge_group[k].z * psi_avg2) - 1);
      if (use.surface_ptr->only_counter_ions &&
	  ((surf_chrg_eq < 0 && charge_group[k].z < 0)
	   || (surf_chrg_eq > 0 && charge_group[k].z > 0)))
	new_g2 = -ratio_aq;
      if (new_g2 <= -ratio_aq)
	new_g2 = -ratio_aq + G_TOL * 1e-3;
      if (fabs (new_g) >= 1)
      {
	if (fabs ((new_g - x[j]->surface_charge->g[k].g) / new_g) >
	    convergence_tolerance)
	{
	  converge = FALSE;
	}
      }
      else
      {
	if (fabs (new_g - x[j]->surface_charge->g[k].g) >
	    convergence_tolerance)
	{
	  converge = FALSE;
	}
      }
      x[j]->surface_charge->g[k].g = new_g;
      if (new_g != 0)
      {
	x[j]->surface_charge->g[k].dg = (new_g2 - new_g) / dif;
      }
      else
      {
	x[j]->surface_charge->g[k].dg = -charge_group[k].z;
      }
      /* save g for species */
      for (i = 0; i < count_s_x; i++)
      {
	if (equal (charge_group[k].z, s_x[i]->z, G_TOL) == TRUE)
	{
	  s_x[i]->diff_layer[count_charge].charge = x[j]->surface_charge;
	  s_x[i]->diff_layer[count_charge].count_g = k;
	}
      }
    }
    if (debug_diffuse_layer == TRUE)
    {
      strcpy (name, x[j]->master[0]->elt->name);
      replace ("_psi", "", name);
/*			surf_chrg_eq = calc_surface_charge(name);
 */
      output_msg (OUTPUT_MESSAGE,
		  "\nDonnan all on %s (%d): charge, \tg, \tdg, Psi_surface = %8f V. \n",
		  name, count_charge,
		  x[j]->master[0]->s->la * 2 * LOG_10 * R_KJ_DEG_MOL * tk_x /
		  F_KJ_V_EQ);
      for (i = 0; i < count_g; i++)
      {
	output_msg (OUTPUT_MESSAGE, "\t%12f\t%12.4e\t%12.4e\n",
		    (double) x[j]->surface_charge->g[i].charge,
		    (double) x[j]->surface_charge->g[i].g,
		    (double) x[j]->surface_charge->g[i].dg);
      }
    }
    count_charge++;
  }
  return (converge);
}

/* ---------------------------------------------------------------------- */
int
calc_init_donnan (void)
/* ---------------------------------------------------------------------- */
{
  int i, j, k;
  int count_g, count_charge;
  char name[MAX_LENGTH];
  LDBLE f_psi, surf_chrg_eq, psi_avg, f_sinh, A_surf, ratio_aq;

  if (use.surface_ptr == NULL)
    return (OK);
  if (use.surface_ptr->type == CD_MUSIC)
    return (calc_init_donnan_music ());
  f_sinh =
    sqrt (8000.0 * EPSILON * EPSILON_ZERO * (R_KJ_DEG_MOL * 1000.0) * tk_x *
	  mu_x);
  if (convergence_tolerance >= 1e-8)
  {
    G_TOL = 1e-9;
  }
  else
  {
    G_TOL = 1e-13;
  }
/*
 *  sum eq of each charge number in solution...
 */
  charge_group = (struct Charge_Group *) free_check_null (charge_group);
  charge_group =
    (struct Charge_Group *) PHRQ_malloc ((size_t)
					 sizeof (struct Charge_Group));
  if (charge_group == NULL)
    malloc_error ();
  charge_group[0].z = 0.0;
  charge_group[0].eq = 0.0;

  count_g = 1;
  for (i = 0; i < count_s_x; i++)
  {
    if (s_x[i]->type > HPLUS)
      continue;
    for (k = 0; k < count_g; k++)
    {
      if (equal (charge_group[k].z, s_x[i]->z, G_TOL) == TRUE)
      {
	charge_group[k].eq += s_x[i]->z * s_x[i]->moles;
	break;
      }
    }
    if (k >= count_g)
    {
      charge_group =
	(struct Charge_Group *) PHRQ_realloc (charge_group,
					      (size_t) (count_g +
							1) *
					      sizeof (struct Charge_Group));
      if (charge_group == NULL)
	malloc_error ();
      charge_group[count_g].z = s_x[i]->z;
      charge_group[count_g].eq = s_x[i]->z * s_x[i]->moles;

      count_g++;
    }
  }
/*
 *   calculate g for each surface...
 */
  count_charge = 0;
  for (j = 0; j < count_unknowns; j++)
  {
    if (x[j]->type != SURFACE_CB)
      continue;
    surface_charge_ptr = x[j]->surface_charge;

    x[j]->surface_charge->g =
      (struct surface_diff_layer *) PHRQ_malloc ((size_t) count_g *
						 sizeof (struct
							 surface_diff_layer));
    if (x[j]->surface_charge->g == NULL)
      malloc_error ();
    x[j]->surface_charge->count_g = count_g;

    /* find surface charge from potential... */
    A_surf =
      x[j]->surface_charge->specific_area * x[j]->surface_charge->grams;
    f_psi = x[j]->master[0]->s->la * LOG_10;
    surf_chrg_eq = A_surf * f_sinh * sinh (f_psi) / F_C_MOL;

    /* find psi_avg that matches surface charge... */
    psi_avg = calc_psi_avg (0);	/*(surf_chrg_eq); */

    /* fill in g's */
    ratio_aq = surface_charge_ptr->mass_water / mass_water_aq_x;

    for (k = 0; k < count_g; k++)
    {
      x[j]->surface_charge->g[k].charge = charge_group[k].z;
      x[j]->surface_charge->g[k].g =
	ratio_aq * (exp (-charge_group[k].z * psi_avg) - 1);
      if (use.surface_ptr->only_counter_ions
	  && ((surf_chrg_eq < 0 && charge_group[k].z < 0)
	      || (surf_chrg_eq > 0 && charge_group[k].z > 0)))
	x[j]->surface_charge->g[k].g = -ratio_aq;
      if (x[j]->surface_charge->g[k].g != 0)
      {
	x[j]->surface_charge->g[k].dg = -A_surf * f_sinh * cosh (f_psi) /
	  (charge_group[k].eq * F_C_MOL);
      }
      else
      {
	x[j]->surface_charge->g[k].dg = -charge_group[k].z;
      }
      /* save g for species */
      for (i = 0; i < count_s_x; i++)
      {
	if (equal (charge_group[k].z, s_x[i]->z, G_TOL) == TRUE)
	{
	  s_x[i]->diff_layer[count_charge].charge = x[j]->surface_charge;
	  s_x[i]->diff_layer[count_charge].count_g = k;
	  s_x[i]->diff_layer[count_charge].g_moles = 0.0;
	  s_x[i]->diff_layer[count_charge].dg_g_moles = 0.0;
	}
      }
    }
    if (debug_diffuse_layer == TRUE)
    {
      strcpy (name, x[j]->master[0]->elt->name);
      replace ("_psi", "", name);
/*			surf_chrg_eq = calc_surface_charge(name);
 */
      output_msg (OUTPUT_MESSAGE,
		  "\nDonnan init on %s : charge, \tg, \tdg, Psi_surface = %8f V. \n",
		  name,
		  x[j]->master[0]->s->la * 2 * LOG_10 * R_KJ_DEG_MOL * tk_x /
		  F_KJ_V_EQ);
      for (i = 0; i < count_g; i++)
      {
	output_msg (OUTPUT_MESSAGE, "\t%12f\t%12.4e\t%12.4e\n",
		    (double) x[j]->surface_charge->g[i].charge,
		    (double) x[j]->surface_charge->g[i].g,
		    (double) x[j]->surface_charge->g[i].dg);
      }
    }
    count_charge++;
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
LDBLE
calc_psi_avg (LDBLE surf_chrg_eq)
/* ---------------------------------------------------------------------- */
{
/*
 * calculate the average (F * Psi / RT) that lets the DL charge counter the surface charge
 */
  int i, iter, count_g;
  LDBLE fd, fd1, p, temp, ratio_aq;
/*	LDBLE dif;
 */
  count_g = surface_charge_ptr->count_g;
  ratio_aq = surface_charge_ptr->mass_water / mass_water_aq_x;
  p = 0;
  if (surf_chrg_eq == 0)
    return (0.0);
  else if (surf_chrg_eq < 0)
    p = -0.5 * log (-surf_chrg_eq * ratio_aq / mu_x + 1);
  else if (surf_chrg_eq > 0)
    p = 0.5 * log (surf_chrg_eq * ratio_aq / mu_x + 1);
/*
 * Optimize p in SS{s_x[i]->moles * z_i * g(p)} = -surf_chrg_eq
 *  g(p) = exp(-p * z_i) * ratio_aq
 * Elsewhere in PHREEQC, g is the excess, after subtraction of conc's for p = 0:
 *                      g(p) = (exp(-p *z_i) - 1) * ratio_aq
 */
  iter = 0;
  do
  {
    fd = surf_chrg_eq;
    fd1 = 0.0;
    for (i = 1; i < count_g; i++)
    {
      if (use.surface_ptr->type == CD_MUSIC)
	temp = exp (-charge_group[i].z * p);
      else
	/*  multiply with ratio_aq for multiplier options cp and cm
	   in calc_all_donnan (not used now)...  */
	temp = exp (-charge_group[i].z * p) * ratio_aq;

      if (use.surface_ptr->only_counter_ions &&
	  ((surf_chrg_eq < 0 && charge_group[i].z < 0)
	   || (surf_chrg_eq > 0 && charge_group[i].z > 0)))
	temp = 0.0;
      fd += charge_group[i].eq * temp;
      fd1 -= charge_group[i].z * charge_group[i].eq * temp;
    }
    fd /= -fd1;
    p += (fd > 1) ? 1 : ((fd < -1) ? -1 : fd);
    if (fabs (p) < G_TOL)
      p = 0.0;
    iter++;
    if (iter > 50)
    {
      sprintf (error_string,
	       "\nToo many iterations for surface in subroutine calc_psi_avg.\n");
      error_msg (error_string, STOP);
    }
  }
  while (fabs (fd) > 1e-12 && p != 0.0);
  if (debug_diffuse_layer == TRUE)
    output_msg (OUTPUT_MESSAGE,
		"iter in calc_psi_avg = %d. g(+1) = %8f. surface charge = %8f.\n",
		iter, (double) (exp (-p) - 1), (double) surf_chrg_eq);

  return (p);
}

/* ---------------------------------------------------------------------- */
int
calc_all_donnan_music (void)
/* ---------------------------------------------------------------------- */
{
  int i, j, k;
  int count_g, count_charge, converge;
  char name[MAX_LENGTH];
  LDBLE new_g, f_psi, surf_chrg_eq, psi_avg, f_sinh, A_surf, ratio_aq;
  LDBLE cz, cm, cp;

  if (use.surface_ptr == NULL)
    return (OK);
  f_sinh =
    sqrt (8000.0 * EPSILON * EPSILON_ZERO * (R_KJ_DEG_MOL * 1000.0) * tk_x *
	  mu_x);
  cz = cm = 1.0;
  cp = 1.0;
/*
 *   calculate g for each surface...
 */
  converge = TRUE;
  count_charge = 0;
  for (j = 0; j < count_unknowns; j++)
  {
    if (x[j]->type != SURFACE_CB)
      continue;
    surface_charge_ptr = x[j]->surface_charge;

    if (debug_diffuse_layer == TRUE)
      output_msg (OUTPUT_MESSAGE, "Calc_all_g, X[%d]\n", j);
/*
 *  sum eq of each charge number in solution...
 */
    count_g = x[j]->surface_charge->count_g;
    for (i = 0; i < count_g; i++)
    {
      charge_group[i].eq = 0.0;
    }
    for (i = 0; i < count_s_x; i++)
    {
      if (s_x[i]->type > HPLUS)
	continue;
      for (k = 0; k < count_g; k++)
      {
	if (equal (charge_group[k].z, s_x[i]->z, G_TOL) == TRUE)
	{
	  charge_group[k].eq += s_x[i]->z * s_x[i]->moles;
	  break;
	}
      }
    }
    /* find surface charge from potential... */
    A_surf =
      x[j]->surface_charge->specific_area * x[j]->surface_charge->grams;
    f_psi = x[j + 2]->master[0]->s->la * LOG_10;	/* -FPsi/RT */
    f_psi = f_psi / 2;
    surf_chrg_eq = A_surf * f_sinh * sinh (f_psi) / F_C_MOL;
    psi_avg = calc_psi_avg (surf_chrg_eq);

    /*output_msg(OUTPUT_MESSAGE, "psi's  %e %e %e\n", f_psi, psi_avg, surf_chrg_eq); */

    /* fill in g's */
    ratio_aq = surface_charge_ptr->mass_water / mass_water_aq_x;

    for (k = 0; k < count_g; k++)
    {
      x[j]->surface_charge->g[k].charge = charge_group[k].z;
      new_g = exp (charge_group[k].z * psi_avg) - 1;
      if (new_g < -ratio_aq)
	new_g = -ratio_aq + G_TOL * 1e-5;
      if (fabs (new_g) >= 1)
      {
	if (fabs ((new_g - x[j]->surface_charge->g[k].g) / new_g) >
	    convergence_tolerance)
	{
	  converge = FALSE;
	}
      }
      else
      {
	if (fabs (new_g - x[j]->surface_charge->g[k].g) >
	    convergence_tolerance)
	{
	  converge = FALSE;
	}
      }
      x[j]->surface_charge->g[k].g = new_g;
      /* save g for species */
      for (i = 0; i < count_s_x; i++)
      {
	if (equal (charge_group[k].z, s_x[i]->z, G_TOL) == TRUE)
	{
	  s_x[i]->diff_layer[count_charge].charge = x[j]->surface_charge;
	  s_x[i]->diff_layer[count_charge].count_g = k;
	}
      }
    }
    if (debug_diffuse_layer == TRUE)
    {
      strcpy (name, x[j + 2]->master[0]->elt->name);
      replace ("_psi", "", name);
/*			surf_chrg_eq = calc_surface_charge(name);
 */
      output_msg (OUTPUT_MESSAGE,
		  "\nDonnan all on %s (%d): charge, \tg, \tdg, Psi_surface = %8f V. \n",
		  name, count_charge,
		  x[j]->master[0]->s->la * LOG_10 * R_KJ_DEG_MOL * tk_x /
		  F_KJ_V_EQ);
      for (i = 0; i < count_g; i++)
      {
	output_msg (OUTPUT_MESSAGE, "\t%12f\t%12.4e\t%12.4e\n",
		    (double) x[j]->surface_charge->g[i].charge,
		    (double) x[j]->surface_charge->g[i].g,
		    (double) x[j]->surface_charge->g[i].dg);
      }
    }
    count_charge++;
  }
  return (converge);
}

/* ---------------------------------------------------------------------- */
int
calc_init_donnan_music (void)
/* ---------------------------------------------------------------------- */
{
  int i, j, k;
  int count_g, count_charge;
  char name[MAX_LENGTH];
  LDBLE psi_avg, f_sinh, ratio_aq;

  if (use.surface_ptr == NULL)
    return (OK);
  f_sinh =
    sqrt (8000.0 * EPSILON * EPSILON_ZERO * (R_KJ_DEG_MOL * 1000.0) * tk_x);
  if (convergence_tolerance >= 1e-8)
  {
    G_TOL = 1e-9;
  }
  else
  {
    G_TOL = 1e-10;
  }
/*
 *  sum eq of each charge number in solution...
 */
  charge_group = (struct Charge_Group *) free_check_null (charge_group);
  charge_group =
    (struct Charge_Group *) PHRQ_malloc ((size_t)
					 sizeof (struct Charge_Group));
  if (charge_group == NULL)
    malloc_error ();
  charge_group[0].z = 0.0;
  charge_group[0].eq = 0.0;

  count_g = 1;
  for (i = 0; i < count_s_x; i++)
  {
    if (s_x[i]->type > HPLUS)
      continue;
    for (k = 0; k < count_g; k++)
    {
      if (equal (charge_group[k].z, s_x[i]->z, G_TOL) == TRUE)
      {
	charge_group[k].eq += s_x[i]->z * s_x[i]->moles;
	break;
      }
    }
    if (k >= count_g)
    {
      charge_group =
	(struct Charge_Group *) PHRQ_realloc (charge_group,
					      (size_t) (count_g +
							1) *
					      sizeof (struct Charge_Group));
      if (charge_group == NULL)
	malloc_error ();
      charge_group[count_g].z = s_x[i]->z;
      charge_group[count_g].eq = s_x[i]->z * s_x[i]->moles;

      count_g++;
    }
  }
/*
 *   calculate g for each surface...
 */
  count_charge = 0;
  for (j = 0; j < count_unknowns; j++)
  {
    if (x[j]->type != SURFACE_CB)
      continue;
    surface_charge_ptr = x[j]->surface_charge;

    x[j]->surface_charge->g =
      (struct surface_diff_layer *) PHRQ_malloc ((size_t) count_g *
						 sizeof (struct
							 surface_diff_layer));
    if (x[j]->surface_charge->g == NULL)
      malloc_error ();
    x[j]->surface_charge->count_g = count_g;

    /* find psi_avg that matches surface charge... */

    psi_avg = calc_psi_avg (0);

    /* fill in g's */
    ratio_aq = surface_charge_ptr->mass_water / mass_water_aq_x;

    for (k = 0; k < count_g; k++)
    {
      x[j]->surface_charge->g[k].charge = charge_group[k].z;
      x[j]->surface_charge->g[k].g = exp (charge_group[k].z * psi_avg) - 1;
      /* save g for species */
      for (i = 0; i < count_s_x; i++)
      {
	if (equal (charge_group[k].z, s_x[i]->z, G_TOL) == TRUE)
	{
	  s_x[i]->diff_layer[count_charge].charge = x[j]->surface_charge;
	  s_x[i]->diff_layer[count_charge].count_g = k;
	  s_x[i]->diff_layer[count_charge].g_moles = 0.0;
	  s_x[i]->diff_layer[count_charge].dg_g_moles = 0.0;
	}
      }
    }
    if (debug_diffuse_layer == TRUE)
    {
      strcpy (name, x[j + 2]->master[0]->elt->name);
      replace ("_psi", "", name);
/*			surf_chrg_eq = calc_surface_charge(name);
 */
      output_msg (OUTPUT_MESSAGE,
		  "\nDonnan all on %s (%d): charge, \tg, \tdg, Psi_surface = %8f V. \n",
		  name, count_charge,
		  x[j]->master[0]->s->la * LOG_10 * R_KJ_DEG_MOL * tk_x /
		  F_KJ_V_EQ);
      for (i = 0; i < count_g; i++)
      {
	output_msg (OUTPUT_MESSAGE, "\t%12f\t%12.4e\t%12.4e\n",
		    (double) x[j]->surface_charge->g[i].charge,
		    (double) x[j]->surface_charge->g[i].g,
		    (double) x[j]->surface_charge->g[i].dg);
      }
    }
    count_charge++;
  }
  return (OK);
}
