#define  EXTERNAL extern
#include "global.h"
#include "phqalloc.h"
#include "output.h"
#include "phrqproto.h"

#include "sundialstypes.h"	/* definitions of types realtype and                        */
							 /* integertype, and the constant FALSE            */
#include "cvode.h"		/* prototypes for CVodeMalloc, CVode, and            */
							 /* CVodeFree, constants OPT_SIZE, BDF, NEWTON,   */
							 /* SV, SUCCESS, NST,NFE,NSETUPS, NNI, NCFN, NETF */
#include "cvdense.h"		/* prototype for CVDense, constant DENSE_NJE    */
#include "nvector_serial.h"	/* definitions of type N_Vector and macro          */
							 /* NV_Ith_S, prototypes for N_VNew, N_VFree      */
#include "dense.h"		/* definitions of type DenseMat, macro DENSE_ELEM */
#define KINETICS_EXTERNAL extern
#include "kinetics.h"

// MDL: to see the renamed function P_step()
#include "step.h"

/* These macros are defined in order to write code which exactly matches
   the mathematical problem description given above.

   Ith(v,i) references the ith component of the vector v, where i is in
   the range [1..NEQ] and NEQ is defined below. The Ith macro is defined
   using the N_VIth macro in nvector.h. N_VIth numbers the components of
   a vector starting from 0.

   IJth(A,i,j) references the (i,j)th element of the dense matrix A, where
   i and j are in the range [1..NEQ]. The IJth macro is defined using the
   DENSE_ELEM macro in dense.h. DENSE_ELEM numbers rows and columns of a
   dense matrix starting from 0. */

#define Ith(v,i)	NV_Ith_S(v,i-1)	/* Ith numbers components 1..NEQ */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1)	/* IJth numbers rows,cols 1..NEQ */

static void f (integertype N, realtype t, N_Vector y, N_Vector ydot,
	       void *f_data);

static void Jac (integertype N, DenseMat J, RhsFn f, void *f_data, realtype t,
		 N_Vector y, N_Vector fy, N_Vector ewt, realtype h,
		 realtype uround, void *jac_data, long int *nfePtr,
		 N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);

static char const svnid[] =
  "$Id: kinetics.c 4 2009-04-21 17:29:29Z delucia $";
static int calc_final_kinetic_reaction (struct kinetics *kinetics_ptr);
static int calc_kinetic_reaction (struct kinetics *kinetics_ptr,
				  LDBLE time_step);
static int rk_kinetics (int i, LDBLE kin_time, int use_mix, int nsaver,
			LDBLE step_fraction);
static int set_reaction (int i, int use_mix, int use_kinetics);
static int set_transport (int i, int use_mix, int use_kinetics, int nsaver);
static int store_get_equi_reactants (int k, int kin_end);

#define MAX_DIVIDE 2
#define KINETICS_TOL 1e-8;

LDBLE *m_original;
LDBLE *m_temp;

extern LDBLE min_value;
extern int count_total_steps;
#ifdef SKIP
/* appt */
static int change_surf_flag;
#endif
/* ---------------------------------------------------------------------- */
int
calc_kinetic_reaction (struct kinetics *kinetics_ptr, LDBLE time_step)
/* ---------------------------------------------------------------------- */
{
/*
 *	Go through kinetic components to
 *	determine rates and
 *	a list of elements and amounts in
 *	the reaction.
 */
  int i, j, return_value;
  LDBLE coef;
/*	char token[MAX_LENGTH];
	char *ptr;
	struct phase *phase_ptr;
 */ char command[] = "run";
  struct rate *rate_ptr;
/*	LDBLE t1, t2; */
  if (svnid == NULL)
    fprintf (stderr, " ");
/*
 *   Go through list and generate list of elements and
 *   coefficient of elements in reaction
 */
  return_value = OK;
  count_elts = 0;
  paren_count = 0;

  rate_time = time_step;

/*	t1 = clock(); */
  for (i = 0; i < kinetics_ptr->count_comps; i++)
  {
    coef = 0.0;
/*
 *   Send command to basic interpreter
 */
    rate_ptr = rate_search (kinetics_ptr->comps[i].rate_name, &j);
    if (rate_ptr == NULL)
    {
      sprintf (error_string, "Rate not found for %s",
	       kinetics_ptr->comps[i].rate_name);
      error_msg (error_string, STOP);
    }
    else
    {
      rate_moles = NAN;
      rate_m = kinetics_ptr->comps[i].m;
      rate_m0 = kinetics_ptr->comps[i].m0;
      rate_p = kinetics_ptr->comps[i].d_params;
      count_rate_p = kinetics_ptr->comps[i].count_d_params;
      if (rate_ptr->new_def == TRUE)
      {
	if (basic_compile
	    (rates[j].commands, &rates[j].linebase, &rates[j].varbase,
	     &rates[j].loopbase) != 0)
	{
	  sprintf (error_string, "Fatal Basic error in rate %s.",
		   kinetics_ptr->comps[i].rate_name);
	  error_msg (error_string, STOP);
	}

	rate_ptr->new_def = FALSE;
      }
      if (basic_run
	  (command, rates[j].linebase, rates[j].varbase,
	   rates[j].loopbase) != 0)
      {
	sprintf (error_string, "Fatal Basic error in rate %s.",
		 kinetics_ptr->comps[i].rate_name);
	error_msg (error_string, STOP);
      }
      if (rate_moles == NAN)
      {
	sprintf (error_string, "Moles of reaction not SAVE'd for %s.",
		 kinetics_ptr->comps[i].rate_name);
	error_msg (error_string, STOP);
      }
      else
      {

	coef = rate_moles;
      }
    }
/*
 *   Accumulate moles of reaction for component
 */
    kinetics_ptr->comps[i].moles += coef;
    if (coef == 0.0)
      continue;
  }
/*	t2=clock();
	printf("secs in reac %e, t2 %e\n", t2-t1, t1);
 */
  return (return_value);
}

/* ---------------------------------------------------------------------- */
int
calc_final_kinetic_reaction (struct kinetics *kinetics_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *	Go through kinetic components to
 *	using extrapolated values, which were
 *	stored in moles in run_kinetics
 */
  int i, j, k;
  LDBLE coef;
  char token[MAX_LENGTH];
  char *ptr;
  struct phase *phase_ptr;
  struct master *master_ptr;
/*
 *   Go through list and generate list of elements and
 *   coefficient of elements in reaction
 */
  kinetics_ptr->totals =
    (struct elt_list *) free_check_null (kinetics_ptr->totals);
  count_elts = 0;
  paren_count = 0;
  for (i = 0; i < kinetics_ptr->count_comps; i++)
  {
    if (kinetics_ptr->comps[i].moles > m_temp[i])
    {
      kinetics_ptr->comps[i].moles = m_temp[i];
      kinetics_ptr->comps[i].m = 0;
    }
    coef = kinetics_ptr->comps[i].moles;
    if (coef == 0.0)
      continue;
/*
 *   Reactant is a pure phase, copy formula into token
 */
    for (j = 0; j < kinetics_ptr->comps[i].count_list; j++)
    {
      phase_ptr = NULL;
      strcpy (token, kinetics_ptr->comps[i].list[j].name);
      phase_ptr = phase_bsearch (token, &k, FALSE);
      if (phase_ptr != NULL)
      {
	add_elt_list (phase_ptr->next_elt,
		      coef * kinetics_ptr->comps[i].list[j].coef);
      }
      else
      {
	ptr = kinetics_ptr->comps[i].list[j].name;
	get_elts_in_species (&ptr,
			     coef * kinetics_ptr->comps[i].list[j].coef);
      }
    }
#ifdef SKIP
    phase_ptr = NULL;
    if (kinetics_ptr->comps[i].count_list == 1)
    {
      strcpy (token, kinetics_ptr->comps[i].list[0].name);
      phase_ptr = phase_bsearch (token, &j, FALSE);
    }
    if (phase_ptr != NULL)
    {
      add_elt_list (phase_ptr->next_elt,
		    coef * kinetics_ptr->comps[i].list[0].coef);
    }
    else
    {
      for (j = 0; j < kinetics_ptr->comps[i].count_list; j++)
      {
	ptr = kinetics_ptr->comps[i].list[j].name;
	get_elts_in_species (&ptr,
			     coef * kinetics_ptr->comps[i].list[j].coef);
      }
    }
#endif
    if (use.exchange_ptr != NULL && use.exchange_ptr->related_rate == TRUE)
    {
      for (j = 0; j < use.exchange_ptr->count_comps; j++)
      {
	if (use.exchange_ptr->comps[j].rate_name != NULL)
	{
	  if (strcmp_nocase
	      (kinetics_ptr->comps[i].rate_name,
	       use.exchange_ptr->comps[j].rate_name) == 0)
	  {
	    /* found kinetics component */
	    add_elt_list (use.exchange_ptr->comps[j].formula_totals,
			  -coef *
			  use.exchange_ptr->comps[j].phase_proportion);
	  }
	}
      }
    }
    if (use.surface_ptr != NULL && use.surface_ptr->related_rate == TRUE)
    {
      for (j = 0; j < use.surface_ptr->count_comps; j++)
      {
	if (use.surface_ptr->comps[j].rate_name != NULL)
	{
	  if (strcmp_nocase
	      (kinetics_ptr->comps[i].rate_name,
	       use.surface_ptr->comps[j].rate_name) == 0)
	  {
	    /* found kinetics component */
	    ptr = use.surface_ptr->comps[j].formula;
/* Surface = 0 when m becomes low ...
 */
	    if (0.9 * use.surface_ptr->comps[j].phase_proportion *
		(kinetics_ptr->comps[i].m) < MIN_RELATED_SURFACE)
	    {
	      master_ptr = master_bsearch (ptr);
	      master_ptr->total = 0.0;
	    }
	    else
	    {
	      get_elts_in_species (&ptr,
				   -coef *
				   use.surface_ptr->comps[j].
				   phase_proportion);
	    }
	  }
	}
      }
    }
  }
  kinetics_ptr->totals = elt_list_save ();
  /*
     output_msg(OUTPUT_MESSAGE, "Calc_final_kinetic_reaction                                                 \n");
     elt_list_print(kinetics_ptr->totals);
   */
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
rk_kinetics (int i, LDBLE kin_time, int use_mix, int nsaver,
	     LDBLE step_fraction)
/* ---------------------------------------------------------------------- */
{
/*
 * Runge-Kutta-Fehlberg method; 6 evaluations of the derivative
 *	give O(h^5) global error and error estimate
 *	calc_kinetic_reaction(.., ..) calculates moles of intermediate derivatives;
 *	these are calc'd for the whole step h.
 *	calc_final_kinetic reaction(..) translates moles to PHREEQC reaction.
 */
  int j, k, m, save_old;
  int bad, step_bad, step_ok;
  int n_reactions;
  LDBLE h, h_old, h_sum;
  LDBLE *rk_moles;
  LDBLE error, error_max, safety, moles_max, moles_reduction;
  struct kinetics *kinetics_ptr;
  int equal_rate, zero_rate;

  struct pp_assemblage *pp_assemblage_save = NULL;
  struct s_s_assemblage *s_s_assemblage_save = NULL;

#ifdef SKIP
  struct gas_phase *gas_phase_save = NULL;
  struct solution *solution_save = NULL;
  struct exchange *exchange_save = NULL;
  struct surface *surface_save = NULL;
#endif

  static LDBLE b31 = 3. / 40., b32 = 9. / 40.,
    b51 = -11. / 54., b53 = -70. / 27., b54 = 35. / 27.,
    b61 = 1631. / 55296., b62 = 175. / 512., b63 = 575. / 13824., b64 =
    44275. / 110592., b65 = 253. / 4096., c1 = 37. / 378., c3 =
    250. / 621., c4 = 125. / 594., c6 = 512. / 1771., dc5 = -277. / 14336.;
  LDBLE dc1 = c1 - 2825. / 27648., dc3 = c3 - 18575. / 48384., dc4 =
    c4 - 13525. / 55296., dc6 = c6 - 0.25;
#ifdef SKIP
/* appt */
  change_surf_flag = FALSE;
#endif
/*
 *  Save kinetics i and solution i, if necessary
 */
  save_old = -2 - (count_cells * (1 + stag_data->count_stag) + 2);
  kinetics_duplicate (i, save_old);
  if (nsaver != i)
  {
    solution_duplicate (i, save_old);
  }

/*
 *   Malloc some space
 */
  if (kinetics_bsearch (i, &m) == NULL)
    return (OK);
  n_reactions = kinetics[m].count_comps;
  rk_moles =
    (LDBLE *) PHRQ_malloc ((size_t) 6 * n_reactions * sizeof (LDBLE));
  if (rk_moles == NULL)
    malloc_error ();

  /*if (use_mix != NOMIX) last_model.force_prep = TRUE; */
  set_and_run_wrapper (i, use_mix, FALSE, i, step_fraction);
  run_reactions_iterations += iterations;

  saver ();
  if (state == TRANSPORT || state == PHAST)
  {
    set_transport (i, NOMIX, TRUE, i);
  }
  else if (state == ADVECTION)
  {
    set_advection (i, NOMIX, TRUE, i);
  }
  else if (state == REACTION)
  {
    set_reaction (i, NOMIX, TRUE);
  }
  if (use.pp_assemblage_ptr != NULL)
  {
    pp_assemblage_save =
      (struct pp_assemblage *) PHRQ_malloc (sizeof (struct pp_assemblage));
    if (pp_assemblage_save == NULL)
      malloc_error ();
  }
  if (use.s_s_assemblage_ptr != NULL)
  {
    s_s_assemblage_save =
      (struct s_s_assemblage *) PHRQ_malloc (sizeof (struct s_s_assemblage));
    if (s_s_assemblage_save == NULL)
      malloc_error ();
  }

  kinetics_ptr = kinetics_bsearch (i, &m);

  step_bad = step_ok = 0;
  bad = FALSE;
  h_sum = 0.;
  h = h_old = kin_time;
  moles_max = 0.1;
  moles_reduction = 1.0;
  safety = 0.7;
  if (kinetics_ptr->rk < 1)
    kinetics_ptr->rk = 1;
  else if (kinetics_ptr->rk > 3)
    kinetics_ptr->rk = 6;

  if (kinetics_ptr->rk == 6)
    equal_rate = FALSE;
  else
    equal_rate = TRUE;
/*
 * if step_divide > 1, initial timestep is divided
 * if			 < 1, step_divide indicates maximal reaction...
 */
  if (kinetics_ptr->step_divide > 1.0)
  {
    h = h_old = kin_time / kinetics_ptr->step_divide;
    equal_rate = FALSE;
  }
  else if (kinetics_ptr->step_divide < 1.0)
    moles_max = kinetics_ptr->step_divide;

  rate_sim_time = rate_sim_time_start + h_sum;

  status (0, NULL);
  while (h_sum < kin_time)
  {

    if (step_bad > kinetics_ptr->bad_step_max)
    {
      sprintf (error_string,
	       "Bad RK steps > %d. Please decrease (time)step or increase -bad_step_max.",
	       kinetics_ptr->bad_step_max);
      error_msg (error_string, STOP);
    }

  MOLES_TOO_LARGE:
    if (moles_reduction > 1.0)
    {
      h_old = h;
      h = safety * h / (1.0 + moles_reduction);
      moles_reduction = 1.0;
      equal_rate = FALSE;
      bad = TRUE;
    }
/*
 *   find k1
 */
    if (bad == TRUE)
    {
      for (j = 0; j < n_reactions; j++)
      {
	rk_moles[j] *= (h / h_old);
	kinetics_ptr->comps[j].moles = rk_moles[j] * 0.2;
	kinetics_ptr->comps[j].m = m_temp[j];
      }
      bad = FALSE;
    }
    else
    {
/*
 *   define pointers for calc_kinetic_, they are lost after saver()...
 */
      if (state == TRANSPORT || state == PHAST)
      {
	set_transport (i, NOMIX, TRUE, i);
      }
      else if (state == ADVECTION)
      {
	set_advection (i, NOMIX, TRUE, i);
      }
      else if (state == REACTION)
      {
	set_reaction (i, NOMIX, TRUE);
      }
      /*
       *   Moles of minerals and solid solutions may change to make positive
       *   concentrations. Reactions may take out more than is present in
       *   solution.
       */
      if (pp_assemblage_save != NULL)
      {
	pp_assemblage_copy (use.pp_assemblage_ptr, pp_assemblage_save,
			    use.pp_assemblage_ptr->n_user);
      }
      if (s_s_assemblage_save != NULL)
      {
	s_s_assemblage_copy (use.s_s_assemblage_ptr, s_s_assemblage_save,
			     use.s_s_assemblage_ptr->n_user);
      }
      for (j = 0; j < n_reactions; j++)
      {
	kinetics_ptr->comps[j].moles = 0.;
	m_temp[j] = kinetics_ptr->comps[j].m;
      }

      rate_sim_time = rate_sim_time_start + h_sum;
      calc_kinetic_reaction (kinetics_ptr, h);

      /* store k1 in rk_moles ... */
      for (j = 0; j < n_reactions; j++)
      {
	if (moles_reduction * moles_max < fabs (kinetics_ptr->comps[j].moles))
	{
	  moles_reduction = fabs (kinetics_ptr->comps[j].moles) / moles_max;
	}
	/*  define reaction for calculating k2 ... */
	rk_moles[j] = kinetics_ptr->comps[j].moles;
	kinetics_ptr->comps[j].moles *= 0.2;
      }
      if (moles_reduction > 1.0)
	goto MOLES_TOO_LARGE;
    }
/*
 * Quit rk with rk = 1 and equal rates ...
 */
    if (kinetics_ptr->rk == 1 && equal_rate)
    {
      zero_rate = TRUE;
      for (j = 0; j < n_reactions; j++)
      {
	kinetics_ptr->comps[j].moles = rk_moles[j];
	if (fabs (kinetics_ptr->comps[j].moles) > MIN_TOTAL)
	  zero_rate = FALSE;
      }

      if (zero_rate == FALSE)
      {
	calc_final_kinetic_reaction (kinetics_ptr);
	for (j = 0; j < n_reactions; j++)
	{
	  kinetics_ptr->comps[j].m = m_temp[j] - kinetics_ptr->comps[j].moles;
	  if (kinetics_ptr->comps[j].m < 1.e-30)
	    kinetics_ptr->comps[j].m = 0;
	  kinetics_ptr->comps[j].moles = 0.;
	}
	if (set_and_run_wrapper (i, NOMIX, TRUE, i, 0.) == MASS_BALANCE)
	{
	  run_reactions_iterations += iterations;
	  moles_reduction = 9;
	  goto MOLES_TOO_LARGE;
	}
	run_reactions_iterations += iterations;
	calc_kinetic_reaction (kinetics_ptr, h);
	for (j = 0; j < n_reactions; j++)
	{
	  if (fabs (rk_moles[j] - kinetics_ptr->comps[j].moles) >
	      kinetics_ptr->comps[j].tol)
	  {
	    equal_rate = FALSE;
	    break;
	  }
	}
      }
      if (zero_rate || equal_rate)
      {
	/* removing the following line causes different results for 
	   example 6 distributed with the program */
	saver ();

	/*  Free space */

	if (pp_assemblage_save != NULL)
	{
	  pp_assemblage_free (pp_assemblage_save);
	}
	if (s_s_assemblage_save != NULL)
	{
	  s_s_assemblage_free (s_s_assemblage_save);
	}
	goto EQUAL_RATE_OUT;
      }
      else
      {
	kinetics_ptr->rk = 3;
	for (j = 0; j < n_reactions; j++)
	{
	  kinetics_ptr->comps[j].moles = 0.2 * rk_moles[j];
	}
      }
    }
/*
 * Continue with rk ...
 */
    calc_final_kinetic_reaction (kinetics_ptr);
    if (set_and_run_wrapper (i, NOMIX, TRUE, i, 0.) == MASS_BALANCE)
    {
      run_reactions_iterations += iterations;
      moles_reduction = 9;
      goto MOLES_TOO_LARGE;
    }
    run_reactions_iterations += iterations;

/*
 *   find k2
 */
    for (j = 0; j < n_reactions; j++)
    {
      kinetics_ptr->comps[j].m = m_temp[j] - kinetics_ptr->comps[j].moles;
      kinetics_ptr->comps[j].moles = 0.;
    }
    rate_sim_time = rate_sim_time_start + h_sum + 0.2 * h;
    calc_kinetic_reaction (kinetics_ptr, h);

    /*   Reset to values of last saver() */
    if (pp_assemblage_save != NULL)
    {
      pp_assemblage_free (use.pp_assemblage_ptr);
      pp_assemblage_copy (pp_assemblage_save, use.pp_assemblage_ptr,
			  pp_assemblage_save->n_user);
    }
    if (s_s_assemblage_save != NULL)
    {
      s_s_assemblage_free (use.s_s_assemblage_ptr);
      s_s_assemblage_copy (s_s_assemblage_save, use.s_s_assemblage_ptr,
			   s_s_assemblage_save->n_user);
    }

    /* store k2 in rk_moles */
    k = n_reactions;
    for (j = 0; j < n_reactions; j++)
    {
      if (moles_reduction * moles_max < fabs (kinetics_ptr->comps[j].moles))
      {
	moles_reduction = fabs (kinetics_ptr->comps[j].moles) / moles_max;
      }
      /*  define reaction for calculating k3 */
      rk_moles[k + j] = kinetics_ptr->comps[j].moles;

      kinetics_ptr->comps[j].moles = b31 * rk_moles[j]
	+ b32 * rk_moles[k + j];
/*
 * check for equal_rate ...
 */
      if (equal_rate
	  && fabs (rk_moles[j] - rk_moles[k + j]) >
	  kinetics_ptr->comps[j].tol)
      {
	equal_rate = FALSE;
      }
    }
    if (moles_reduction > 1.0)
      goto MOLES_TOO_LARGE;
/*
 * Quit rk with rk = 2 and equal rates ...
 */
    if (kinetics_ptr->rk == 2 && equal_rate)
    {
      for (j = 0; j < n_reactions; j++)
      {
	kinetics_ptr->comps[j].moles =
	  0.3 * rk_moles[j] + 0.7 * rk_moles[k + j];
      }
      calc_final_kinetic_reaction (kinetics_ptr);
      for (j = 0; j < n_reactions; j++)
      {
	kinetics_ptr->comps[j].m = m_temp[j] - kinetics_ptr->comps[j].moles;
	if (kinetics_ptr->comps[j].m < 1.e-30)
	  kinetics_ptr->comps[j].m = 0;
	kinetics_ptr->comps[j].moles = 0.;
      }
      if (set_and_run_wrapper (i, NOMIX, TRUE, i, 0.) == MASS_BALANCE)
      {
	run_reactions_iterations += iterations;
	moles_reduction = 9;
	goto MOLES_TOO_LARGE;
      }
      run_reactions_iterations += iterations;
/*
 * Move next calc'n to rk = 1 when initial rate equals final rate ...
 */
      calc_kinetic_reaction (kinetics_ptr, h);
      for (j = 0; j < n_reactions; j++)
      {
	if (fabs (rk_moles[j] - kinetics_ptr->comps[j].moles) >
	    kinetics_ptr->comps[j].tol)
	{
	  equal_rate = FALSE;
	  break;
	}
      }
      if (equal_rate)
	kinetics_ptr->rk = 1;

      saver ();

      /*  Free space */

      if (pp_assemblage_save != NULL)
      {
	pp_assemblage_free (pp_assemblage_save);
      }
      if (s_s_assemblage_save != NULL)
      {
	s_s_assemblage_free (s_s_assemblage_save);
      }
      goto EQUAL_RATE_OUT;
    }
/*
 * Continue runge_kutta..
 */
    calc_final_kinetic_reaction (kinetics_ptr);
    if (set_and_run_wrapper (i, NOMIX, TRUE, i, 0.) == MASS_BALANCE)
    {
      run_reactions_iterations += iterations;
      moles_reduction = 9;
      goto MOLES_TOO_LARGE;
    }
    run_reactions_iterations += iterations;
/*
 *   find k3
 */
    for (j = 0; j < n_reactions; j++)
    {
      kinetics_ptr->comps[j].m = m_temp[j] - kinetics_ptr->comps[j].moles;
      kinetics_ptr->comps[j].moles = 0.;
    }
    rate_sim_time = rate_sim_time_start + h_sum + 0.3 * h;
    calc_kinetic_reaction (kinetics_ptr, h);

    /*   Reset to values of last saver() */
    if (pp_assemblage_save != NULL)
    {
      pp_assemblage_free (use.pp_assemblage_ptr);
      pp_assemblage_copy (pp_assemblage_save, use.pp_assemblage_ptr,
			  pp_assemblage_save->n_user);
    }
    if (s_s_assemblage_save != NULL)
    {
      s_s_assemblage_free (use.s_s_assemblage_ptr);
      s_s_assemblage_copy (s_s_assemblage_save, use.s_s_assemblage_ptr,
			   s_s_assemblage_save->n_user);
    }

    /* store k3 in rk_moles */
    k = 2 * n_reactions;
    for (j = 0; j < n_reactions; j++)
    {
      if (moles_reduction * moles_max < fabs (kinetics_ptr->comps[j].moles))
      {
	moles_reduction = fabs (kinetics_ptr->comps[j].moles) / moles_max;
      }
      /*  define reaction for calculating k4 ... */
      rk_moles[k + j] = kinetics_ptr->comps[j].moles;

      kinetics_ptr->comps[j].moles = 0.3 * rk_moles[j]
	- 0.9 * rk_moles[n_reactions + j] + 1.2 * rk_moles[k + j];
/*
 * check for equal_rate ...
 */
      if (equal_rate
	  && fabs (rk_moles[j] - rk_moles[k + j]) >
	  kinetics_ptr->comps[j].tol)
	equal_rate = FALSE;
    }
    if (moles_reduction > 1.0)
      goto MOLES_TOO_LARGE;
/*
 * Quit rk with rk = 3 and equal rates ...
 */
    if (kinetics_ptr->rk == 3 && equal_rate)
    {
      for (j = 0; j < n_reactions; j++)
      {
	kinetics_ptr->comps[j].moles = 0.5 * rk_moles[j]
	  - 1.5 * rk_moles[n_reactions + j] + 2 * rk_moles[k + j];
      }
      calc_final_kinetic_reaction (kinetics_ptr);
      for (j = 0; j < n_reactions; j++)
      {
	kinetics_ptr->comps[j].m = m_temp[j] - kinetics_ptr->comps[j].moles;
	if (kinetics_ptr->comps[j].m < 1.e-30)
	  kinetics_ptr->comps[j].m = 0;
	kinetics_ptr->comps[j].moles = 0.;
      }

      if (set_and_run_wrapper (i, NOMIX, TRUE, i, 0.) == MASS_BALANCE)
      {
	run_reactions_iterations += iterations;
	moles_reduction = 9;
	goto MOLES_TOO_LARGE;
      }
      run_reactions_iterations += iterations;
/*
 * Move next calc'n to rk = 1 when initial rate equals final rate ...
 */
      calc_kinetic_reaction (kinetics_ptr, h);
      for (j = 0; j < n_reactions; j++)
      {
	if (fabs (rk_moles[j] - kinetics_ptr->comps[j].moles) >
	    kinetics_ptr->comps[j].tol)
	{
	  equal_rate = FALSE;
	  break;
	}
      }
      if (equal_rate)
	kinetics_ptr->rk = 1;

      saver ();

      /*  Free space */

      if (pp_assemblage_save != NULL)
      {
	pp_assemblage_free (pp_assemblage_save);
      }
      if (s_s_assemblage_save != NULL)
      {
	s_s_assemblage_free (s_s_assemblage_save);
      }
      goto EQUAL_RATE_OUT;
    }
/*
 * Continue runge_kutta..
 */

    calc_final_kinetic_reaction (kinetics_ptr);
    if (set_and_run_wrapper (i, NOMIX, TRUE, i, 0.) == MASS_BALANCE)
    {
      run_reactions_iterations += iterations;
      moles_reduction = 9;
      goto MOLES_TOO_LARGE;
    }
    run_reactions_iterations += iterations;
/*
 *   find k4
 */
    for (j = 0; j < n_reactions; j++)
    {
      kinetics_ptr->comps[j].m = m_temp[j] - kinetics_ptr->comps[j].moles;
      kinetics_ptr->comps[j].moles = 0.;
    }
    rate_sim_time = rate_sim_time_start + h_sum + 0.6 * h;
    calc_kinetic_reaction (kinetics_ptr, h);

    /*   Reset to values of last saver() */
    if (pp_assemblage_save != NULL)
    {
      pp_assemblage_free (use.pp_assemblage_ptr);
      pp_assemblage_copy (pp_assemblage_save, use.pp_assemblage_ptr,
			  pp_assemblage_save->n_user);
    }
    if (s_s_assemblage_save != NULL)
    {
      s_s_assemblage_free (use.s_s_assemblage_ptr);
      s_s_assemblage_copy (s_s_assemblage_save, use.s_s_assemblage_ptr,
			   s_s_assemblage_save->n_user);
    }

    /* store k4 in rk_moles */
    k = 3 * n_reactions;
    for (j = 0; j < n_reactions; j++)
    {
      if (moles_reduction * moles_max < fabs (kinetics_ptr->comps[j].moles))
      {
	moles_reduction = fabs (kinetics_ptr->comps[j].moles) / moles_max;
      }

      /*  define reaction for calculating k5 */
      rk_moles[k + j] = kinetics_ptr->comps[j].moles;
      kinetics_ptr->comps[j].moles = b51 * rk_moles[j]
	+ 2.5 * rk_moles[n_reactions + j]
	+ b53 * rk_moles[2 * n_reactions + j] + b54 * rk_moles[k + j];
    }
    if (moles_reduction > 1.0)
      goto MOLES_TOO_LARGE;
    calc_final_kinetic_reaction (kinetics_ptr);
    if (set_and_run_wrapper (i, NOMIX, TRUE, i, 0.) == MASS_BALANCE)
    {
      run_reactions_iterations += iterations;
      moles_reduction = 9;
      goto MOLES_TOO_LARGE;
    }
    run_reactions_iterations += iterations;
/*
 *   find k5
 */
    for (j = 0; j < n_reactions; j++)
    {
      kinetics_ptr->comps[j].m = m_temp[j] - kinetics_ptr->comps[j].moles;
      kinetics_ptr->comps[j].moles = 0.;
    }
    rate_sim_time = rate_sim_time_start + h_sum + h;
    calc_kinetic_reaction (kinetics_ptr, h);

    /*   Reset to values of last saver() */
    if (pp_assemblage_save != NULL)
    {
      pp_assemblage_free (use.pp_assemblage_ptr);
      pp_assemblage_copy (pp_assemblage_save, use.pp_assemblage_ptr,
			  pp_assemblage_save->n_user);
    }
    if (s_s_assemblage_save != NULL)
    {
      s_s_assemblage_free (use.s_s_assemblage_ptr);
      s_s_assemblage_copy (s_s_assemblage_save, use.s_s_assemblage_ptr,
			   s_s_assemblage_save->n_user);
    }

    /* store k5 in rk_moles */
    k = 4 * n_reactions;
    for (j = 0; j < n_reactions; j++)
    {
      if (moles_reduction * moles_max < fabs (kinetics_ptr->comps[j].moles))
      {
	moles_reduction = fabs (kinetics_ptr->comps[j].moles) / moles_max;
      }

      /*  define reaction for calculating k6 */
      rk_moles[k + j] = kinetics_ptr->comps[j].moles;
      kinetics_ptr->comps[j].moles = b61 * rk_moles[j]
	+ b62 * rk_moles[n_reactions + j]
	+ b63 * rk_moles[2 * n_reactions + j]
	+ b64 * rk_moles[3 * n_reactions + j] + b65 * rk_moles[k + j];
    }
    if (moles_reduction > 1.0)
      goto MOLES_TOO_LARGE;
    calc_final_kinetic_reaction (kinetics_ptr);
    if (set_and_run_wrapper (i, NOMIX, TRUE, i, 0.) == MASS_BALANCE)
    {
      run_reactions_iterations += iterations;
      moles_reduction = 9;
      goto MOLES_TOO_LARGE;
    }
    run_reactions_iterations += iterations;
/*
 *   find k6
 */
    for (j = 0; j < n_reactions; j++)
    {
      kinetics_ptr->comps[j].m = m_temp[j] - kinetics_ptr->comps[j].moles;
      kinetics_ptr->comps[j].moles = 0.;
    }
    rate_sim_time = rate_sim_time_start + h_sum + 0.875 * h;
    calc_kinetic_reaction (kinetics_ptr, h);

    /*   Reset to values of last saver() */
    if (pp_assemblage_save != NULL)
    {
      pp_assemblage_free (use.pp_assemblage_ptr);
      pp_assemblage_copy (pp_assemblage_save, use.pp_assemblage_ptr,
			  pp_assemblage_save->n_user);
    }
    if (s_s_assemblage_save != NULL)
    {
      s_s_assemblage_free (use.s_s_assemblage_ptr);
      s_s_assemblage_copy (s_s_assemblage_save, use.s_s_assemblage_ptr,
			   s_s_assemblage_save->n_user);
    }

    /* store k6 in rk_moles */
    k = 5 * n_reactions;
    for (j = 0; j < n_reactions; j++)
    {
      rk_moles[k + j] = kinetics_ptr->comps[j].moles;
    }

/*
 *   Evaluate error
 */
    error_max = 0.;
    for (j = 0; j < n_reactions; j++)
    {
      error = fabs (dc1 * rk_moles[j]
		    + dc3 * rk_moles[2 * n_reactions + j]
		    + dc4 * rk_moles[3 * n_reactions + j]
		    + dc5 * rk_moles[4 * n_reactions + j]
		    + dc6 * rk_moles[5 * n_reactions + j]);

      /* tol is in moles/l */
      error /= kinetics_ptr->comps[j].tol;
      if (error > error_max)
	error_max = error;
    }

/*
 *   repeat with smaller step
 */
/* printf("timest %g ; error_max %g\n", h, error_max); */
    if (error_max > 1)
    {
      h_old = h;
      if (step_ok == 0)
	h = h * safety / error_max;
      else
	h = h * safety * pow (error_max, -0.25);
      bad = TRUE;
      step_bad++;
    }
    else
    {
/*
 *   OK, calculate result
 */
      for (j = 0; j < n_reactions; j++)
      {
	kinetics_ptr->comps[j].moles = c1 * rk_moles[j]
	  + c3 * rk_moles[2 * n_reactions + j]
	  + c4 * rk_moles[3 * n_reactions + j]
	  + c6 * rk_moles[5 * n_reactions + j];
      }
      calc_final_kinetic_reaction (kinetics_ptr);
      for (j = 0; j < n_reactions; j++)
      {
	kinetics_ptr->comps[j].m = m_temp[j] - kinetics_ptr->comps[j].moles;
	if (kinetics_ptr->comps[j].m < 1.e-30)
	  kinetics_ptr->comps[j].m = 0;
	kinetics_ptr->comps[j].moles = 0.;
      }

      if (set_and_run_wrapper (i, NOMIX, TRUE, i, 0.) == MASS_BALANCE)
      {
	run_reactions_iterations += iterations;
	moles_reduction = 9;
	goto MOLES_TOO_LARGE;
      }
      run_reactions_iterations += iterations;
/*
 * Move next calc'n to rk = 1 when initial rate equals final rate ...
 */
      calc_kinetic_reaction (kinetics_ptr, h);
      for (j = 0; j < n_reactions; j++)
      {
	if (fabs (rk_moles[j] - kinetics_ptr->comps[j].moles) >
	    kinetics_ptr->comps[j].tol)
	{
	  equal_rate = FALSE;
	  break;
	}
      }
      if (equal_rate && kinetics_ptr->rk < 6)
	kinetics_ptr->rk = 1;

      saver ();

      step_ok++;
      h_sum += h;
      /*  Free space */

      if (pp_assemblage_save != NULL)
      {
	pp_assemblage_free (pp_assemblage_save);
      }
      if (s_s_assemblage_save != NULL)
      {
	s_s_assemblage_free (s_s_assemblage_save);
      }
/*
 *   and increase step size ...
 */
      if (h_sum < kin_time)
      {
	if (error_max > 0.000577)
	{
	  h = h * safety * pow (error_max, -0.2);
	}
	else
	{
	  h *= 4;
	}
	if (h > (kin_time - h_sum))
	  h = (kin_time - h_sum);
      }
    }
#if !defined(PHREEQCI_GUI)
#ifndef PHREEQ98
    if (pr.status == TRUE && status_on == TRUE)
    {
      char str[MAX_LENGTH];
      backspace_screen (37);
      sprintf (str, "RK-steps: Bad%4d. OK%5d. Time %3d%%", step_bad, step_ok,
	       (int) (100 * h_sum / kin_time));
      output_msg (OUTPUT_SCREEN, "%-37s", str);
    }
#endif
#endif
  }

EQUAL_RATE_OUT:

/*
 *   Run one more time to get distribution of species
 */
  if (state >= REACTION || nsaver != i)
  {
#ifdef SKIP
/* appt */
    saver ();
    change_surf_flag = TRUE;
#endif
    set_and_run_wrapper (i, NOMIX, FALSE, nsaver, 0.);
    run_reactions_iterations += iterations;
  }
/*	saver();
   *//* reset for printing */
  if (use_mix == DISP)
  {
    use.mix_ptr = &mix[count_mix - count_cells + i - 1];
    use.mix_in = TRUE;
    use.n_mix_user = i;
  }
  else if ((use_mix == STAG || use_mix == TRUE) && state == TRANSPORT)
  {
    use.mix_ptr = mix_search (i, &use.n_mix, FALSE);
    if (use.mix_ptr != NULL)
    {
      use.mix_in = TRUE;
      use.n_mix_user = i;
    }
  }
/*
 *  Restore solution i, if necessary
 */
  if (nsaver != i)
  {
    solution_duplicate (save_old, i);
  }
  rk_moles = (LDBLE *) free_check_null (rk_moles);

#ifdef SKIP
  if (state != TRANSPORT)
  {
#ifdef DOS
    output_msg (OUTPUT_SCREEN, "\n");
#else
    output_msg (OUTPUT_SCREEN, "\n%-80s", " ");
#endif
  }
#endif
  rate_sim_time = rate_sim_time_start + kin_time;
  use.kinetics_in = TRUE;

  /*  Free space */

  if (pp_assemblage_save != NULL)
  {
    pp_assemblage_save =
      (struct pp_assemblage *) free_check_null (pp_assemblage_save);
  }
  if (s_s_assemblage_save != NULL)
  {
    s_s_assemblage_save =
      (struct s_s_assemblage *) free_check_null (s_s_assemblage_save);
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
set_and_run_wrapper (int i, int use_mix, int use_kinetics, int nsaver,
		     LDBLE step_fraction)
/* ---------------------------------------------------------------------- */
{
  int j, converge, max_try;
  int old_diag, old_itmax;
  LDBLE old_tol, old_min_value, old_step, old_pe, old_pp_column_scale;
  LDBLE small_pe_step, small_step;
  struct pp_assemblage *pp_assemblage_save = NULL;
  struct s_s_assemblage *s_s_assemblage_save = NULL;
  struct kinetics *kinetics_save = NULL;

  /* 0 -- normal */
  /* 1 -- try smaller step size, more iterations */
  /* 2 -- try diagonal scaling */
  /* 3 -- try smaller tolerance */
  /* 4 -- try alternate scaling */
  small_pe_step = 5.;
  small_step = 10.;
  converge = FALSE;

  old_diag = diagonal_scale;
  old_itmax = itmax;
  old_tol = ineq_tol;
  old_step = step_size;
  old_pe = pe_step_size;
  old_min_value = min_value;
  old_pp_column_scale = pp_column_scale;

  if (state == TRANSPORT || state == PHAST)
  {
    set_transport (i, use_mix, use_kinetics, i);
  }
  else if (state == ADVECTION)
  {
    set_advection (i, use_mix, use_kinetics, i);
  }
  else if (state == REACTION)
  {
    set_reaction (i, use_mix, use_kinetics);
  }
  if (use.pp_assemblage_ptr != NULL)
  {
    pp_assemblage_save =
      (struct pp_assemblage *) PHRQ_malloc (sizeof (struct pp_assemblage));
    if (pp_assemblage_save == NULL)
      malloc_error ();
    pp_assemblage_copy (use.pp_assemblage_ptr, pp_assemblage_save,
			use.pp_assemblage_ptr->n_user);
  }
  if (use.s_s_assemblage_ptr != NULL)
  {
    s_s_assemblage_save =
      (struct s_s_assemblage *) PHRQ_malloc (sizeof (struct s_s_assemblage));
    if (s_s_assemblage_save == NULL)
      malloc_error ();
    s_s_assemblage_copy (use.s_s_assemblage_ptr, s_s_assemblage_save,
			 use.s_s_assemblage_ptr->n_user);
  }
  if (use.kinetics_ptr != NULL)
  {
    kinetics_save =
      (struct kinetics *) PHRQ_malloc (sizeof (struct kinetics));
    if (kinetics_save == NULL)
      malloc_error ();
    kinetics_copy (use.kinetics_ptr, kinetics_save, use.kinetics_ptr->n_user);
  }

  if (pitzer_model == TRUE)
  {
    diagonal_scale = TRUE;
    always_full_pitzer = FALSE;
    max_try = 2;
  }
  else
  {
    max_try = 11;
  }
  /*max_try = 1; */
  for (j = 0; j < max_try; j++)
  {
    if (j == 1)
    {
      always_full_pitzer = TRUE;
      if (pe_step_size <= small_pe_step && step_size <= small_step)
	continue;
      itmax *= 2;
      step_size = small_step;
      pe_step_size = small_pe_step;
      sprintf (error_string,
	       "Trying smaller step size, pe step size %g, %g ... \n",
	       (double) step_size, (double) pe_step_size);
      warning_msg (error_string);
    }
    else if (j == 2)
    {
      itmax *= 2;
      if (diagonal_scale == TRUE)
      {
	diagonal_scale = FALSE;
      }
      else
      {
	diagonal_scale = TRUE;
      }
      sprintf (error_string, "Trying diagonal scaling ...\n");
      warning_msg (error_string);
    }
    else if (j == 3)
    {
      itmax *= 2;
      ineq_tol /= 10.;
      sprintf (error_string, "Trying reduced tolerance %g ...\n",
	       (double) ineq_tol);
      warning_msg (error_string);
    }
    else if (j == 4)
    {
      itmax *= 2;
      ineq_tol *= 10.;
      sprintf (error_string, "Trying increased tolerance %g ...\n",
	       (double) ineq_tol);
      warning_msg (error_string);
    }
    else if (j == 5)
    {
      itmax *= 2;
      if (diagonal_scale == TRUE)
      {
	diagonal_scale = FALSE;
      }
      else
      {
	diagonal_scale = TRUE;
      }
      ineq_tol /= 10.;
      sprintf (error_string,
	       "Trying diagonal scaling and reduced tolerance %g ...\n",
	       (double) ineq_tol);
      warning_msg (error_string);
    }
    else if (j == 6)
    {
      itmax *= 2;
      pp_column_scale = 1e-10;
      sprintf (error_string, "Trying scaling pure_phase columns %g ...\n",
	       (double) pp_column_scale);
      warning_msg (error_string);
    }
    else if (j == 7)
    {
      itmax *= 2;
      pp_column_scale = 1e-10;
      if (diagonal_scale == TRUE)
      {
	diagonal_scale = FALSE;
      }
      else
      {
	diagonal_scale = TRUE;
      }
      sprintf (error_string,
	       "Trying scaling pure_phase columns and diagonal scale %g ...\n",
	       (double) pp_column_scale);
      warning_msg (error_string);
    }
    else if (j == 8)
    {
      itmax *= 2;
      min_value *= 10;
      sprintf (error_string, "Trying increased scaling %g ...\n",
	       (double) min_value);
      warning_msg (error_string);
    }
    else if (j == 9)
    {
      aqueous_only = 5;
      sprintf (error_string,
	       "Skipping optimize equations for first %d iterations ...\n",
	       aqueous_only);
      warning_msg (error_string);
    }
    else if (j == 10)
    {
      negative_concentrations = TRUE;
      sprintf (error_string,
	       "Adding inequality to make concentrations greater than zero.\n");
      warning_msg (error_string);
    }
    if (j > 0)
    {
      if (pp_assemblage_save != NULL)
      {
	pp_assemblage_free (use.pp_assemblage_ptr);
	pp_assemblage_copy (pp_assemblage_save, use.pp_assemblage_ptr,
			    pp_assemblage_save->n_user);
      }
      if (s_s_assemblage_save != NULL)
      {
	s_s_assemblage_free (use.s_s_assemblage_ptr);
	s_s_assemblage_copy (s_s_assemblage_save, use.s_s_assemblage_ptr,
			     s_s_assemblage_save->n_user);
      }
      if (kinetics_save != NULL)
      {
	kinetics_free (use.kinetics_ptr);
	kinetics_copy (kinetics_save, use.kinetics_ptr,
		       kinetics_save->n_user);
      }
    }

    converge = set_and_run (i, use_mix, use_kinetics, nsaver, step_fraction);
    /* reset values */
    diagonal_scale = old_diag;
    itmax = old_itmax;
    ineq_tol = old_tol;
    step_size = old_step;
    pe_step_size = old_pe;
    min_value = old_min_value;
    pp_column_scale = old_pp_column_scale;
    aqueous_only = 0;
    negative_concentrations = FALSE;
    if (converge == TRUE)
    {
      break;
    }
    else if (converge == MASS_BALANCE)
    {
      break;
    }
    warning_msg
      ("Numerical method failed with this set of convergence parameters.\n");
  }
  if (converge == FALSE && use.kinetics_ptr != NULL
      && use.kinetics_ptr->use_cvode == TRUE)
  {
    sprintf (error_string,
	     "Numerical method failed on all parameter combinations, retrying integration");
    /*output_msg(OUTPUT_MESSAGE, "Numerical method failed on all parameter combinations, retrying integration\n"); */
    warning_msg (error_string);
    converge = MASS_BALANCE;
  }
  if (converge == FALSE)
  {
/*
 *   write to error.inp what failed to converge.
 */
    if (output_open (OUTPUT_DUMP, "error.inp") != OK)
    {
      sprintf (error_string, "Can't open file, %s.", "error.inp");
      error_msg (error_string, CONTINUE);
      input_error++;
    }
    else
    {
      if (use.irrev_in == TRUE)
      {
	dump_reaction (use.n_mix_user);
      }
      if (use.kinetics_ptr != NULL)
      {
	dump_kinetics (use.kinetics_ptr->n_user);
      }
      output_msg (OUTPUT_DUMP, "END\n");
      if (use.solution_ptr != NULL)
      {
	dump_solution (use.n_solution_user);
      }
      else if (use.mix_ptr != NULL)
      {
	dump_mix (use.n_mix_user);
      }
      if (use.pp_assemblage_in == TRUE)
      {
	dump_pp_assemblage (use.n_pp_assemblage_user);
      }
      if (use.exchange_in == TRUE)
      {
	dump_exchange (use.n_exchange_user);
      }
      if (use.surface_in == TRUE)
      {
	dump_surface (use.n_surface_user);
      }
      if (use.gas_phase_in == TRUE)
      {
	dump_gas_phase (use.n_gas_phase_user);
      }
      if (use.s_s_assemblage_in == TRUE)
      {
	dump_s_s_assemblage (use.n_s_s_assemblage_user);
      }
      output_msg (OUTPUT_DUMP, "END\n");
    }
    /* if (state == TRANSPORT && dump_modulus == 0) dump(); */
    check_residuals ();
    pr.all = TRUE;
    pr.gas_phase = use.gas_phase_in;
    pr.pp_assemblage = use.pp_assemblage_in;
    pr.s_s_assemblage = use.s_s_assemblage_in;
    pr.surface = use.surface_in;
    pr.exchange = use.exchange_in;
    pr.totals = TRUE;
    pr.species = TRUE;
    pr.saturation_indices = TRUE;
    pr.irrev = use.irrev_in;
    pr.mix = use.mix_in;
    pr.reaction = TRUE;
    pr.use = TRUE;
    sum_species ();
    print_all ();
    sprintf (error_string,
	     "Numerical method failed on all combinations of convergence parameters");
    error_msg (error_string, STOP);
  }
  if (pp_assemblage_save != NULL)
  {
    pp_assemblage_free (pp_assemblage_save);
    pp_assemblage_save =
      (struct pp_assemblage *) free_check_null (pp_assemblage_save);
  }
  if (s_s_assemblage_save != NULL)
  {
    s_s_assemblage_free (s_s_assemblage_save);
    s_s_assemblage_save =
      (struct s_s_assemblage *) free_check_null (s_s_assemblage_save);
  }
  if (kinetics_save != NULL)
  {
    kinetics_free (kinetics_save);
    kinetics_save = (struct kinetics *) free_check_null (kinetics_save);
  }
  if (converge == MASS_BALANCE)
  {
    return (MASS_BALANCE);
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
set_and_run (int i, int use_mix, int use_kinetics, int nsaver,
	     LDBLE step_fraction)
/* ---------------------------------------------------------------------- */
{
/*
 *   i			--user number for soln, reaction, etc.
 *   use_mix	--integer flag
					state == TRANSPORT: DISP, STAG, NOMIX
					state == REACTION: TRUE, FALSE
 *   use_kinetics --true or false flag to calculate kinetic reactions
 *   nsaver	   --user number to store solution
 *   step_fraction--fraction of irreversible reaction to add
 */
  int n, n1, n2, converge;
  if (state == TRANSPORT || state == PHAST)
  {
    set_transport (i, use_mix, use_kinetics, nsaver);
  }
  else if (state == ADVECTION)
  {
    set_advection (i, use_mix, use_kinetics, nsaver);
  }
  else if (state == REACTION)
  {
    set_reaction (i, use_mix, use_kinetics);
  }
  cell = i;
/*
 *   Take step
 */
  if (state >= REACTION)
  {
    // MDL: function step redefined to P_step
    if (P_step (step_fraction) == MASS_BALANCE)
    {
      return (MASS_BALANCE);
    }
/*
 *  Always use solution, exchange, and surface -1
 */
    use.solution_ptr = solution_bsearch (-1, &n, TRUE);

    /* new */
    if (use.exchange_ptr != NULL)
    {
      use.exchange_ptr = exchange_bsearch (-1, &n1);
    }
    if (use.surface_ptr != NULL)
    {
      use.surface_ptr = surface_bsearch (-1, &n2);
    }
  }
  /* end new */
#ifdef SKIP
/* a quick one ...*/
  if (state == TRANSPORT /*&& ishift == 0 */  &&
      (cell == 1 || cell == count_cells))
    last_model.force_prep = TRUE;
#endif
  if (use.surface_ptr != NULL)
  {
    dl_type_x = use.surface_ptr->dl_type;
  }
  if (use.surface_ptr != NULL && dl_type_x != NO_DL)
  {
    converge = surface_model ();
  }
  else
  {
    prep ();
#ifdef SKIP
    k_temp (solution[n]->tc);
#endif
    k_temp (use.solution_ptr->tc);
    set (FALSE);
    converge = model ();
  }
  sum_species ();
  return (converge);
}

/* ---------------------------------------------------------------------- */
int
set_transport (int i, int use_mix, int use_kinetics, int nsaver)
/* ---------------------------------------------------------------------- */
{
/*
 *   i			--user number for soln, reaction, etc.
 *   use_mix	  --integer flag
					state == TRANSPORT: DISP, STAG, NOMIX
					state == REACTION: TRUE, FALSE
 *   use_kinetics --true or false flag to calculate kinetic reactions
 *   nsaver	   --user number to store solution
 */
  int n;

  cell = i;
  reaction_step = 1;
#ifdef SKIP
  if (pr.use == TRUE && pr.all == TRUE)
  {
    output_msg (OUTPUT_MESSAGE, "\nCell %d\n", i);
  }
#endif
/*
 *   Find mixture or solution
 */

  use.mix_ptr = NULL;
  use.mix_in = FALSE;
  if (use_mix == DISP)
  {
    use.mix_ptr = &mix[count_mix - count_cells + i - 1];
    use.mix_in = TRUE;
    use.n_mix_user = i;
    use.n_mix_user_orig = i;
  }
  else if (use_mix == STAG && multi_Dflag != TRUE)
  {
    use.mix_ptr = mix_search (i, &use.n_mix, FALSE);
    if (use.mix_ptr != NULL)
    {
      use.mix_in = TRUE;
      use.n_mix_user = i;
      use.n_mix_user_orig = i;
    }
    else
    {
      use.solution_ptr = solution_bsearch (i, &use.n_solution, FALSE);
      if (use.solution_ptr == NULL)
      {
	sprintf (error_string, "Solution %d not found.", use.n_solution_user);
	error_msg (error_string, STOP);
      }
      use.n_solution_user = i;
      use.solution_in = TRUE;
    }
  }
  else
  {
    use.solution_ptr = solution_bsearch (i, &use.n_solution, FALSE);
    if (use.solution_ptr == NULL)
    {
      sprintf (error_string, "Solution %d not found.", use.n_solution_user);
      error_msg (error_string, STOP);
    }
    use.n_solution_user = i;
    use.solution_in = TRUE;
  }
  save.solution = TRUE;
  save.n_solution_user = nsaver;
  save.n_solution_user_end = nsaver;
/*
 *   Find pure phase assemblage
 */

  use.pp_assemblage_ptr = pp_assemblage_bsearch (i, &use.n_pp_assemblage);
  if (use.pp_assemblage_ptr != NULL)
  {
    use.pp_assemblage_in = TRUE;
    use.n_pp_assemblage_user = i;
    save.pp_assemblage = TRUE;
    save.n_pp_assemblage_user = i;
    save.n_pp_assemblage_user_end = i;
  }
  else
  {
    use.pp_assemblage_in = FALSE;
    save.pp_assemblage = FALSE;
  }
/*
 *   Find irreversible reaction
 */
  use.irrev_ptr = irrev_bsearch (i, &use.n_irrev);
  if (use.irrev_ptr != NULL)
  {
    use.irrev_in = TRUE;
    use.n_irrev_user = i;
  }
  else
  {
    use.irrev_in = FALSE;
  }
/*
 *   Find exchange
 */
  use.exchange_ptr = exchange_bsearch (i, &use.n_exchange);
  if (use.exchange_ptr != NULL)
  {
    use.exchange_in = TRUE;
    use.n_exchange_user = i;
    save.exchange = TRUE;
    save.n_exchange_user = i;
    save.n_exchange_user_end = i;
  }
  else
  {
    use.exchange_in = FALSE;
    save.exchange = FALSE;
  }

/*
 *   Find surface
 */
  use.surface_ptr = surface_bsearch (i, &use.n_surface);
  if (use.surface_ptr != NULL)
  {
    use.surface_in = TRUE;
    use.n_surface_user = i;
    save.surface = TRUE;
    save.n_surface_user = i;
    save.n_surface_user_end = i;
#ifdef SKIP
/* appt */
    if (change_surf_flag)
    {
      for (n = 0; n < change_surf_count; n++)
      {
	if (change_surf[n].cell_no != i)
	  break;
	reformat_surf (change_surf[n].comp_name, change_surf[n].fraction,
		       change_surf[n].new_comp_name, change_surf[n].new_Dw,
		       change_surf[n].cell_no);
      }
      change_surf_count = 0;
      change_surf_flag = FALSE;
    }
#endif
  }
  else
  {
    use.surface_in = FALSE;
    save.surface = FALSE;
    dl_type_x = NO_DL;
  }
/*
 *   Find temperature;  temp retardation is done in step
 */
  use.temperature_ptr = temperature_bsearch (i, &use.n_temperature);
  if (use.temperature_ptr != NULL)
  {
    use.temperature_in = TRUE;
    use.n_temperature_user = i;
  }
  else
  {
    use.temperature_in = FALSE;
  }
/*
 *   Find gas
 */
  use.gas_phase_ptr = gas_phase_bsearch (i, &use.n_gas_phase);
  if (use.gas_phase_ptr != NULL)
  {
    use.gas_phase_in = TRUE;
    use.n_gas_phase_user = i;
    save.gas_phase = TRUE;
    save.n_gas_phase_user = i;
    save.n_gas_phase_user_end = i;
  }
  else
  {
    use.gas_phase_in = FALSE;
    save.gas_phase = FALSE;
  }
/*
 *   Find s_s_assemblage
 */
  use.s_s_assemblage_ptr = s_s_assemblage_bsearch (i, &use.n_s_s_assemblage);
  if (use.s_s_assemblage_ptr != NULL)
  {
    use.s_s_assemblage_in = TRUE;
    use.n_s_s_assemblage_user = i;
    save.s_s_assemblage = TRUE;
    save.n_s_s_assemblage_user = i;
    save.n_s_s_assemblage_user_end = i;
  }
  else
  {
    use.s_s_assemblage_in = FALSE;
    save.s_s_assemblage = FALSE;
  }
/*
 *   Find kinetics
 */
  if (use_kinetics == TRUE
      && (use.kinetics_ptr = kinetics_bsearch (i, &n)) != NULL)
  {
    use.n_kinetics_user = i;
    use.n_kinetics = n;
    use.kinetics_in = TRUE;
    save.kinetics = TRUE;
    save.n_kinetics_user = i;
    save.n_kinetics_user_end = i;
  }
  else
  {
    use.kinetics_ptr = NULL;
    use.kinetics_in = FALSE;
    save.kinetics = FALSE;
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
set_reaction (int i, int use_mix, int use_kinetics)
/* ---------------------------------------------------------------------- */
{
/*
 *   i			--user number for soln, reaction, etc.
 *   use_mix	  --integer flag
					state == TRANSPORT: DISP, STAG, NOMIX
					state == REACTION: TRUE, FALSE
 *   use_kinetics --true or false flag to calculate kinetic reactions
 */
/*
 *   Find mixture or solution
 */
  use.mix_ptr = NULL;
  use.solution_ptr = NULL;
  if (use_mix == TRUE && use.mix_in == TRUE)
  {
    use.mix_ptr = mix_bsearch (i, &use.n_mix);
    if (use.mix_ptr == NULL)
    {
      sprintf (error_string, "MIX %d not found.", i);
      error_msg (error_string, STOP);
    }
  }
  else
  {
    use.solution_ptr = solution_bsearch (i, &use.n_solution, FALSE);
    if (use.solution_ptr == NULL)
    {
      sprintf (error_string, "Solution %d not found.", i);
      error_msg (error_string, STOP);
    }
  }
/*
 *   Find pure phase assemblage
 */
  if (use.pp_assemblage_in == TRUE)
  {
    use.pp_assemblage_ptr = pp_assemblage_bsearch (i, &use.n_pp_assemblage);
    if (use.pp_assemblage_ptr == NULL)
    {
      sprintf (error_string, "PP_ASSEMBLAGE %d not found.", i);
      error_msg (error_string, STOP);
    }
  }

/*
 *   Find irreversible reaction
 */
  if (use.irrev_in == TRUE)
  {
    use.irrev_ptr = irrev_bsearch (i, &use.n_irrev);
    if (use.irrev_ptr == NULL)
    {
      sprintf (error_string, "REACTION %d not found.", i);
      error_msg (error_string, STOP);
    }
  }
/*
 *   Find exchange
 */
  if (use.exchange_in == TRUE)
  {
    use.exchange_ptr = exchange_bsearch (i, &use.n_exchange);
    if (use.exchange_ptr == NULL)
    {
      sprintf (error_string, "EXCHANGE %d not found.", i);
      error_msg (error_string, STOP);
    }
  }
/*
 *   Find surface
 */
  dl_type_x = NO_DL;
  if (use.surface_in == TRUE)
  {
    use.surface_ptr = surface_bsearch (i, &use.n_surface);
    if (use.surface_ptr == NULL)
    {
      sprintf (error_string, "SURFACE %d not found.", i);
      error_msg (error_string, STOP);
    }
  }
/*
 *   Find temperature;  temp retardation is done in step
 */
  if (use.temperature_in == TRUE)
  {
    use.temperature_ptr = temperature_bsearch (i, &use.n_temperature);
    if (use.temperature_ptr == NULL)
    {
      sprintf (error_string, "TEMPERATURE %d not found.", i);
      error_msg (error_string, STOP);
    }
  }
/*
 *   Find gas
 */
  if (use.gas_phase_in == TRUE)
  {
    use.gas_phase_ptr = gas_phase_bsearch (i, &use.n_gas_phase);
    if (use.gas_phase_ptr == NULL)
    {
      sprintf (error_string, "GAS_PHASE %d not found.", i);
      error_msg (error_string, STOP);
    }
  }
/*
 *   Find s_s_assemblage
 */
  if (use.s_s_assemblage_in == TRUE)
  {
    use.s_s_assemblage_ptr =
      s_s_assemblage_bsearch (i, &use.n_s_s_assemblage);
    if (use.s_s_assemblage_ptr == NULL)
    {
      sprintf (error_string, "Solid-solution Assemblage %d not found.", i);
      error_msg (error_string, STOP);
    }
  }
/*
 *   Find kinetics
 */
  if (use_kinetics == TRUE && use.kinetics_in == TRUE)
  {
    use.kinetics_ptr = kinetics_bsearch (i, &use.n_kinetics);
    if (use.kinetics_ptr == NULL)
    {
      sprintf (error_string, "KINETICS %d not found.", i);
      error_msg (error_string, STOP);
    }
  }
  else
  {
    use.kinetics_ptr = NULL;
  }
  return (OK);
}

#ifdef DEBUG_KINETICS
/* ---------------------------------------------------------------------- */
int
dump_kinetics_stderr (int k)
/* ---------------------------------------------------------------------- */
{
/*
 *      Dumps kinetics data
 */
  int i, j, n;
  struct kinetics *kinetics_ptr;
  struct kinetics_comp *kinetics_comp_ptr;

  if ((kinetics_ptr = kinetics_search (k, &n, FALSE)) == NULL)
    return (OK);

  output_msg (OUTPUT_MESSAGE, "KINETICS XXX %d\n", k);
  output_msg (OUTPUT_MESSAGE, "\t-step_divide   %15.6e\n",
	      (double) kinetics[n].step_divide);
  output_msg (OUTPUT_MESSAGE, "\t-cvode_steps   %15.6e\n",
	      (double) kinetics[n].cvode_steps);
  output_msg (OUTPUT_MESSAGE, "\t-cvode_order   %15.6e\n",
	      (double) kinetics[n].cvode_order);

  output_msg (OUTPUT_MESSAGE, "totals\n");
  for (i = 0; kinetics_ptr->totals[i].elt != NULL; i++)
  {
    output_msg (OUTPUT_MESSAGE, "\t%s\t%e\tSolution_total: %e\n", kinetics_ptr->totals[i].elt->name,
		(double) kinetics_ptr->totals[i].coef, total(kinetics_ptr->totals[i].elt->name));
  }

  for (i = 0; i < kinetics[n].count_comps; i++)
  {
    output_msg (OUTPUT_MESSAGE, "%-15s\n", kinetics_ptr->comps[i].rate_name);

    kinetics_comp_ptr = &kinetics_ptr->comps[i];
    output_msg (OUTPUT_MESSAGE, "\t-formula ");
    for (j = 0; j < kinetics_comp_ptr->count_list; j++)
    {
      output_msg (OUTPUT_MESSAGE, "   %s  %12.3e",
		  kinetics_comp_ptr->list[j].name,
		  (double) kinetics_comp_ptr->list[j].coef);
    }
    output_msg (OUTPUT_MESSAGE, "\n");

    output_msg (OUTPUT_MESSAGE, "\t-tol %15.2e\n",
		(double) kinetics_ptr->comps[i].tol);
    output_msg (OUTPUT_MESSAGE, "\t-m0  %15.6e\n",
		(double) kinetics_ptr->comps[i].m0);
    output_msg (OUTPUT_MESSAGE, "\t-m   %15.6e\n",
		(double) kinetics_ptr->comps[i].m);

    if (kinetics_comp_ptr->count_d_params != 0)
    {
      output_msg (OUTPUT_MESSAGE, "\t-parm");
      for (j = 0; j < kinetics_comp_ptr->count_d_params; j++)
      {
	output_msg (OUTPUT_MESSAGE, "%15.6e",
		    (double) kinetics_comp_ptr->d_params[j]);
      }
      output_msg (OUTPUT_MESSAGE, "\n");
    }

/* not dumped:		kinetics_comp_ptr->count_c_params = 0 */
  }
  output_msg (OUTPUT_MESSAGE, "\n");

  output_msg (OUTPUT_SCREEN, "KINETICS XXX %d\n", k);
  output_msg (OUTPUT_SCREEN, "\t-step_divide   %15.6e\n",
	      (double) kinetics[n].step_divide);
  output_msg (OUTPUT_SCREEN, "\t-cvode_steps   %15.6e\n",
	      (double) kinetics[n].cvode_steps);
  output_msg (OUTPUT_SCREEN, "\t-cvode_order   %15.6e\n",
	      (double) kinetics[n].cvode_order);

  output_msg (OUTPUT_SCREEN, "totals\n");
  for (i = 0; kinetics_ptr->totals[i].elt != NULL; i++)
  {
    output_msg (OUTPUT_SCREEN, "\t%s\t%e\tSolution_total: %e\n", kinetics_ptr->totals[i].elt->name,
		(double) kinetics_ptr->totals[i].coef, total(kinetics_ptr->totals[i].elt->name));
  }

  for (i = 0; i < kinetics[n].count_comps; i++)
  {
    output_msg (OUTPUT_SCREEN, "%-15s\n", kinetics_ptr->comps[i].rate_name);

    kinetics_comp_ptr = &kinetics_ptr->comps[i];
    output_msg (OUTPUT_SCREEN, "\t-formula ");
    for (j = 0; j < kinetics_comp_ptr->count_list; j++)
    {
      output_msg (OUTPUT_SCREEN, "   %s  %12.3e",
		  kinetics_comp_ptr->list[j].name,
		  (double) kinetics_comp_ptr->list[j].coef);
    }
    output_msg (OUTPUT_SCREEN, "\n");

    output_msg (OUTPUT_SCREEN, "\t-tol %15.2e\n",
		(double) kinetics_ptr->comps[i].tol);
    output_msg (OUTPUT_SCREEN, "\t-m0  %15.6e\n",
		(double) kinetics_ptr->comps[i].m0);
    output_msg (OUTPUT_SCREEN, "\t-m   %15.6e\n",
		(double) kinetics_ptr->comps[i].m);

    if (kinetics_comp_ptr->count_d_params != 0)
    {
      output_msg (OUTPUT_SCREEN, "\t-parm");
      for (j = 0; j < kinetics_comp_ptr->count_d_params; j++)
      {
	output_msg (OUTPUT_SCREEN, "%15.6e",
		    (double) kinetics_comp_ptr->d_params[j]);
      }
      output_msg (OUTPUT_SCREEN, "\n");
    }

/* not dumped:		kinetics_comp_ptr->count_c_params = 0 */
  }
  output_msg (OUTPUT_SCREEN, "\n");
  return (OK);
}
#endif
/* ---------------------------------------------------------------------- */
int
run_reactions (int i, LDBLE kin_time, int use_mix, LDBLE step_fraction)
/* ---------------------------------------------------------------------- */
{
/*
 * Kinetics calculations
 * Rates and moles of each reaction are calculated in calc_kinetic_reaction
 * Total number of moles in reaction is stored in kinetics[i].totals
 */

  int j, n, converge, iter;
  int pr_all_save;
  int nsaver;
  struct kinetics *kinetics_ptr;
  struct pp_assemblage *pp_assemblage_ptr;
  struct s_s_assemblage *s_s_assemblage_ptr;
  struct Use use_save;
  int save_old, m, n_reactions /*, nok, nbad */ ;

  /* CVODE definitions */
  realtype ropt[OPT_SIZE], reltol, t, tout, tout1, sum_t;
  long int iopt[OPT_SIZE];
  int flag;
/*
 *   Set nsaver
 */
  run_reactions_iterations = 0;
  kin_time_x = kin_time;
  nsaver = i;
  if (state == TRANSPORT || state == PHAST)
  {
    if (use_mix == DISP)
    {
      nsaver = -2;
    }
    else if (use_mix == STAG)
    {
      nsaver = -2 - i;
    }
  }
  if (state == ADVECTION)
  {
    nsaver = -2;
  }
/*
 * Check that reaction exists for this cell ..
 */
  if (kin_time <= 0 ||
      (state == REACTION && use.kinetics_in == FALSE) ||
      (state == TRANSPORT && kinetics_bsearch (i, &n) == NULL) ||
      (state == PHAST && kinetics_bsearch (i, &n) == NULL) ||
      (state == ADVECTION && kinetics_bsearch (i, &n) == NULL))
  {
    converge = set_and_run_wrapper (i, use_mix, FALSE, nsaver, step_fraction);
    if (converge == MASS_BALANCE)
      error_msg ("Negative concentration in system. Stopping calculation.",
		 STOP);
    run_reactions_iterations += iterations;
  }
  else
  {
/*
 *   Save moles of kinetic reactants for printout...
 */
    kinetics_ptr = kinetics_bsearch (i, &n);

    m_temp =
      (LDBLE *) PHRQ_malloc ((size_t) kinetics_ptr->count_comps *
			     sizeof (LDBLE));
    if (m_temp == NULL)
      malloc_error ();

    m_original =
      (LDBLE *) PHRQ_malloc ((size_t) kinetics_ptr->count_comps *
			     sizeof (LDBLE));
    if (m_original == NULL)
      malloc_error ();

    for (j = 0; j < kinetics_ptr->count_comps; j++)
    {
      m_original[j] = kinetics_ptr->comps[j].m;
      m_temp[j] = kinetics_ptr->comps[j].m;
    }
/*
*   Start the loop for timestepping ...
 *   Use either Runge-Kutta-Fehlberg, or final result extrapolation
 */
    pr_all_save = pr.all;
    pr.all = FALSE;
/*
 *   This condition makes output equal for incremental_reactions TRUE/FALSE...
 *		(if (incremental_reactions == FALSE || reaction_step == 1)
 */
    store_get_equi_reactants (i, FALSE);
    if (kinetics_ptr->use_cvode == FALSE)
    {
/* in case dispersivity is not wanted..
			if (multi_Dflag)
				rk_kinetics(i, kin_time, NOMIX, nsaver, step_fraction);
			else
 */ rk_kinetics (i, kin_time, use_mix, nsaver, step_fraction);
    }
    else
    {
      save_old = -2 - (count_cells * (1 + stag_data->count_stag) + 2);
      if (nsaver != i)
      {
	solution_duplicate (i, save_old);
      }
      for (j = 0; j < OPT_SIZE; j++)
      {
	iopt[j] = 0;
	ropt[j] = 0;
      }

/*
 *	Do mix first
 */
      kinetics_ptr = kinetics_bsearch (i, &m);
      n_reactions = kinetics_ptr->count_comps;
      cvode_n_user = i;
      cvode_kinetics_ptr = (void *) kinetics_ptr;
      cvode_n_reactions = n_reactions;
      cvode_rate_sim_time_start = rate_sim_time_start;
      cvode_rate_sim_time = rate_sim_time;

      if (multi_Dflag)
	converge = set_and_run_wrapper (i, NOMIX, FALSE, i, 0.0);
      else
	converge = set_and_run_wrapper (i, use_mix, FALSE, i, 0.0);
      if (converge == MASS_BALANCE)
	error_msg ("Negative concentration in system. Stopping calculation.",
		   STOP);
      saver ();
      pp_assemblage_ptr = pp_assemblage_bsearch (i, &n);
      s_s_assemblage_ptr = s_s_assemblage_bsearch (i, &n);
      if (pp_assemblage_ptr != NULL)
      {
	cvode_pp_assemblage_save =
	  (struct pp_assemblage *)
	  PHRQ_malloc (sizeof (struct pp_assemblage));
	if (cvode_pp_assemblage_save == NULL)
	  malloc_error ();
	pp_assemblage_copy (pp_assemblage_ptr, cvode_pp_assemblage_save,
			    pp_assemblage_ptr->n_user);
      }
      if (s_s_assemblage_ptr != NULL)
      {
	cvode_s_s_assemblage_save =
	  (struct s_s_assemblage *)
	  PHRQ_malloc (sizeof (struct s_s_assemblage));
	if (cvode_s_s_assemblage_save == NULL)
	  malloc_error ();
	s_s_assemblage_copy (s_s_assemblage_ptr, cvode_s_s_assemblage_save,
			     s_s_assemblage_ptr->n_user);
      }

      /* allocate space for CVODE */
      kinetics_machEnv = M_EnvInit_Serial (n_reactions);
      kinetics_y = N_VNew (n_reactions, kinetics_machEnv);	/* Allocate y, abstol vectors */
      if (kinetics_y == NULL)
	malloc_error ();
      cvode_last_good_y = N_VNew (n_reactions, kinetics_machEnv);	/* Allocate y, abstol vectors */
      if (cvode_last_good_y == NULL)
	malloc_error ();
      cvode_prev_good_y = N_VNew (n_reactions, kinetics_machEnv);	/* Allocate y, abstol vectors */
      if (cvode_prev_good_y == NULL)
	malloc_error ();
      kinetics_abstol = N_VNew (n_reactions, kinetics_machEnv);
      if (kinetics_abstol == NULL)
	malloc_error ();

/*
 *	Set y to 0.0
 */
      for (j = 0; j < n_reactions; j++)
      {
	kinetics_ptr->comps[j].moles = 0.0;
	Ith (kinetics_y, j + 1) = 0.0;
	Ith (kinetics_abstol, j + 1) = kinetics_ptr->comps[j].tol;
	/*Ith(abstol,j+1) = 1e-8; */
	/* m_temp[j] = kinetics_ptr->comps[j].m; */
      }
      reltol = 0.0;

      /* Call CVodeMalloc to initialize CVODE:

         NEQ   is the problem size = number of equations
         f       is the user's right hand side function in y'=f(t,y)
         T0     is the initial time
         y       is the initial dependent variable vector
         BDF   specifies the Backward Differentiation Formula
         NEWTON  specifies a Newton iteration
         SV     specifies scalar relative and vector absolute tolerances
         &reltol is a pointer to the scalar relative tolerance
         abstol  is the absolute tolerance vector
         FALSE   indicates there are no optional inputs in iopt and ropt
         iopt is an array used to communicate optional integer input and output
         ropt is an array used to communicate optional real input and output

         A pointer to CVODE problem memory is returned and stored in cvode_mem. */
      /* Don't know what this does */
      /*
         iopt[SLDET] = TRUE;
         cvode_mem = CVodeMalloc(n_reactions, f, 0.0, y, BDF, NEWTON, SV, &reltol, abstol, NULL, NULL, TRUE, iopt, ropt, machEnv);
         cvode_mem = CVodeMalloc(n_reactions, f, 0.0, y, ADAMS, FUNCTIONAL, SV, &reltol, abstol, NULL, NULL, FALSE, iopt, ropt, machEnv);
	 iopt[MXSTEP] is maximum number of steps that CVODE tries.
       */
      iopt[MXSTEP] = kinetics_ptr->cvode_steps;
      iopt[MAXORD] = kinetics_ptr->cvode_order;
      kinetics_cvode_mem =
	CVodeMalloc (n_reactions, f, 0.0, kinetics_y, BDF, NEWTON, SV,
		     &reltol, kinetics_abstol, NULL, NULL, TRUE, iopt, ropt,
		     kinetics_machEnv);
      if (kinetics_cvode_mem == NULL)
	malloc_error ();

      /* Call CVDense to specify the CVODE dense linear solver with the
         user-supplied Jacobian routine Jac. */

      flag = CVDense (kinetics_cvode_mem, Jac, NULL);
      if (flag != SUCCESS)
      {
	error_msg ("CVDense failed.", STOP);
      }
      t = 0;
      tout = kin_time;
      /*ropt[HMAX] = tout/10.; */
      /*ropt[HMIN] = 1e-17; */
      use_save = use;
      flag = CVode (kinetics_cvode_mem, tout, kinetics_y, &t, NORMAL);
      rate_sim_time = rate_sim_time_start + t;
      /*
         printf("At t = %0.4e   y =%14.6e  %14.6e  %14.6e\n",
         t, Ith(y,1), Ith(y,2), Ith(y,3));
       */
      iter = 0;
      sum_t = 0;
    RESTART:
      while (flag != SUCCESS)
      {
	sum_t += cvode_last_good_time;
	sprintf(error_string, "CVode incomplete at cvode_steps %d. Cell: %d\tTime: %e\tCvode calls: %d, continuing...", (int) iopt[NST], cell_no, sum_t, iter + 1);
	output_msg(OUTPUT_STDERR, "%s\n", error_string);
	if (state == PHAST) output_msg(OUTPUT_SEND_MESSAGE + 1, "%s\n", error_string);

#ifdef DEBUG_KINETICS
	 if (iter > 5) dump_kinetics_stderr(cell_no);
#endif

	cvode_last_good_time = 0;
	if (++iter >= kinetics_ptr->bad_step_max)
	{
	  m_temp = (LDBLE *) free_check_null (m_temp);
	  m_original = (LDBLE *) free_check_null (m_original);
	  error_msg ("Repeated restart of integration.", STOP);
	}
	tout1 = tout - sum_t;
	t = 0;
	N_VScale (1.0, cvode_last_good_y, kinetics_y);
	for (j = 0; j < OPT_SIZE; j++)
	{
	  iopt[j] = 0;
	  ropt[j] = 0;
	}
	CVodeFree (kinetics_cvode_mem);	/* Free the CVODE problem memory */
	iopt[MXSTEP] = kinetics_ptr->cvode_steps;
	iopt[MAXORD] = kinetics_ptr->cvode_order;
	kinetics_cvode_mem =
	  CVodeMalloc (n_reactions, f, 0.0, kinetics_y, BDF, NEWTON, SV,
		       &reltol, kinetics_abstol, NULL, NULL, TRUE, iopt,
		       ropt, kinetics_machEnv);
	if (kinetics_cvode_mem == NULL)
	  malloc_error ();

	/* Call CVDense to specify the CVODE dense linear solver with the
	   user-supplied Jacobian routine Jac. */

	flag = CVDense (kinetics_cvode_mem, Jac, NULL);
	if (flag != SUCCESS)
	{
	  error_msg ("CVDense failed.", STOP);
	}
	flag = CVode (kinetics_cvode_mem, tout1, kinetics_y, &t, NORMAL);
	/*
	   sprintf(error_string, "CVode failed, flag=%d.\n", flag);
	   error_msg(error_string, STOP);
	 */
      }
      /*
         odeint(&ystart[-1], n_reactions, 0, kin_time, kinetics_ptr->comps[0].tol, kin_time/kinetics_ptr->step_divide, 0.0, &nok, &nbad, i, nsaver );
       */
      for (j = 0; j < n_reactions; j++)
      {
	kinetics_ptr->comps[j].moles = Ith (kinetics_y, j + 1);
	kinetics_ptr->comps[j].m =
	  m_original[j] - kinetics_ptr->comps[j].moles;
	if (kinetics_ptr->comps[j].m < 0)
	{
	  kinetics_ptr->comps[j].moles = m_original[j];
	  kinetics_ptr->comps[j].m = 0.0;
	}
	/* output_msg(OUTPUT_MESSAGE,"%d y[%d] %g\n", i, j, ystart[j]); */
      }
      if (use.pp_assemblage_ptr != NULL)
      {
	pp_assemblage_free (use.pp_assemblage_ptr);
	pp_assemblage_copy (cvode_pp_assemblage_save, use.pp_assemblage_ptr,
			    i);
      }
      if (use.s_s_assemblage_ptr != NULL)
      {
	s_s_assemblage_free (use.s_s_assemblage_ptr);
	s_s_assemblage_copy (cvode_s_s_assemblage_save,
			     use.s_s_assemblage_ptr, i);
      }
      calc_final_kinetic_reaction (kinetics_ptr);
      if (set_and_run_wrapper (i, NOMIX, TRUE, nsaver, 1.0) == MASS_BALANCE)
      {
	/*error_msg("FAIL 2 after successful integration in CVode", CONTINUE); */
	warning_msg ("FAIL 2 after successful integration in CVode");
	flag = -1;
	goto RESTART;
      }
      for (j = 0; j < kinetics_ptr->count_comps; j++)
      {
	kinetics_ptr->comps[j].m =
	  m_original[j] - kinetics_ptr->comps[j].moles;
      }
/*
 *  Restore solution i, if necessary
 */
      if (nsaver != i)
      {
	solution_duplicate (save_old, i);
      }
      free_cvode ();
      use.mix_in = use_save.mix_in;
      use.mix_ptr = use_save.mix_ptr;
    }

    rate_sim_time = rate_sim_time_start + kin_time;
    store_get_equi_reactants (i, TRUE);
    pr.all = pr_all_save;

    kinetics_ptr = kinetics_bsearch (i, &n);
    for (j = 0; j < kinetics_ptr->count_comps; j++)
    {
      kinetics_ptr->comps[j].moles = m_original[j] - kinetics_ptr->comps[j].m;
/*						if (kinetics_ptr->comps[j].moles < 1.e-15) kinetics_ptr->comps[j].moles = 0.0;
 */ }
    m_temp = (LDBLE *) free_check_null (m_temp);
    m_original = (LDBLE *) free_check_null (m_original);
  }
  iterations = run_reactions_iterations;
  if (cvode_pp_assemblage_save != NULL)
  {
    pp_assemblage_free (cvode_pp_assemblage_save);
    cvode_pp_assemblage_save =
      (struct pp_assemblage *) free_check_null (cvode_pp_assemblage_save);
  }
  if (cvode_s_s_assemblage_save != NULL)
  {
    s_s_assemblage_free (cvode_s_s_assemblage_save);
    cvode_s_s_assemblage_save =
      (struct s_s_assemblage *) free_check_null (cvode_s_s_assemblage_save);
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
free_cvode (void)
/* ---------------------------------------------------------------------- */
{
  if (kinetics_y != NULL)
    N_VFree (kinetics_y);	/* Free vector */
  kinetics_y = NULL;
  if (cvode_last_good_y != NULL)
    N_VFree (cvode_last_good_y);	/* Free vector */
  cvode_last_good_y = NULL;
  if (cvode_prev_good_y != NULL)
    N_VFree (cvode_prev_good_y);	/* Free vector */
  cvode_prev_good_y = NULL;
  if (kinetics_abstol != NULL)
    N_VFree (kinetics_abstol);	/* Free vector */
  kinetics_abstol = NULL;
  if (kinetics_cvode_mem != NULL)
    CVodeFree (kinetics_cvode_mem);	/* Free the CVODE problem memory */
  kinetics_cvode_mem = NULL;
  if (kinetics_machEnv != NULL)
    M_EnvFree_Serial (kinetics_machEnv);	/* Free the machine environment memory */
  kinetics_machEnv = NULL;
  if (cvode_pp_assemblage_save != NULL)
  {
    pp_assemblage_free (cvode_pp_assemblage_save);
    cvode_pp_assemblage_save =
      (struct pp_assemblage *) free_check_null (cvode_pp_assemblage_save);
  }
  if (cvode_s_s_assemblage_save != NULL)
  {
    s_s_assemblage_free (cvode_s_s_assemblage_save);
    cvode_s_s_assemblage_save =
      (struct s_s_assemblage *) free_check_null (cvode_s_s_assemblage_save);
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
set_advection (int i, int use_mix, int use_kinetics, int nsaver)
/* ---------------------------------------------------------------------- */
{
/*
 *   i			--user number for soln, reaction, etc.
 *   use_mix	  --integer flag
					state == TRANSPORT: DISP, STAG, NOMIX
					state == REACTION: TRUE, FALSE
			state == ADVECTION: TRUE, FALSE
 *   use_kinetics --true or false flag to calculate kinetic reactions
 *   nsaver	   --user number to store solution
 */
  int n;

  cell = i;
  reaction_step = 1;
/*
 *   Find mixture or solution
 */

  use.mix_ptr = NULL;
  use.mix_in = FALSE;
  if (use_mix == TRUE &&
      (use.mix_ptr = mix_search (i, &use.n_mix, FALSE)) != NULL)
  {
    use.mix_in = TRUE;
    use.n_mix_user = i;
    use.n_mix_user_orig = i;
    use.n_solution_user = i;
  }
  else
  {
    use.solution_ptr = solution_bsearch (i, &use.n_solution, FALSE);
    if (use.solution_ptr == NULL)
    {
      sprintf (error_string, "Solution %d not found.", use.n_solution_user);
      error_msg (error_string, STOP);
    }
    use.n_solution_user = i;
    use.solution_in = TRUE;
  }
  save.solution = TRUE;
  save.n_solution_user = nsaver;
  save.n_solution_user_end = nsaver;
/*
 *   Find pure phase assemblage
 */

  use.pp_assemblage_ptr = pp_assemblage_bsearch (i, &use.n_pp_assemblage);
  if (use.pp_assemblage_ptr != NULL)
  {
    use.pp_assemblage_in = TRUE;
    use.n_pp_assemblage_user = i;
    save.pp_assemblage = TRUE;
    save.n_pp_assemblage_user = i;
    save.n_pp_assemblage_user_end = i;
  }
  else
  {
    use.pp_assemblage_in = FALSE;
    save.pp_assemblage = FALSE;
  }
/*
 *   Find irreversible reaction
 */
  use.irrev_ptr = irrev_bsearch (i, &use.n_irrev);
  if (use.irrev_ptr != NULL)
  {
    use.irrev_in = TRUE;
    use.n_irrev_user = i;
  }
  else
  {
    use.irrev_in = FALSE;
  }
/*
 *   Find exchange
 */
  use.exchange_ptr = exchange_bsearch (i, &use.n_exchange);
  if (use.exchange_ptr != NULL)
  {
    use.exchange_in = TRUE;
    use.n_exchange_user = i;
    save.exchange = TRUE;
    save.n_exchange_user = i;
    save.n_exchange_user_end = i;
  }
  else
  {
    use.exchange_in = FALSE;
    save.exchange = FALSE;
  }

/*
 *   Find surface
 */
  use.surface_ptr = surface_bsearch (i, &use.n_surface);
  if (use.surface_ptr != NULL)
  {
    use.surface_in = TRUE;
    use.n_surface_user = i;
    save.surface = TRUE;
    save.n_surface_user = i;
    save.n_surface_user_end = i;
  }
  else
  {
    use.surface_in = FALSE;
    save.surface = FALSE;
    dl_type_x = NO_DL;
  }
/*
 *   Find temperature;  temp retardation is done in step
 */
  use.temperature_ptr = temperature_bsearch (i, &use.n_temperature);
  if (use.temperature_ptr != NULL)
  {
    use.temperature_in = TRUE;
    use.n_temperature_user = i;
  }
  else
  {
    use.temperature_in = FALSE;
  }
/*
 *   Find gas
 */
  use.gas_phase_ptr = gas_phase_bsearch (i, &use.n_gas_phase);
  if (use.gas_phase_ptr != NULL)
  {
    use.gas_phase_in = TRUE;
    use.n_gas_phase_user = i;
    save.gas_phase = TRUE;
    save.n_gas_phase_user = i;
    save.n_gas_phase_user_end = i;
  }
  else
  {
    use.gas_phase_in = FALSE;
    save.gas_phase = FALSE;
  }
/*
 *   Find s_s_assemblage
 */
  use.s_s_assemblage_ptr = s_s_assemblage_bsearch (i, &use.n_s_s_assemblage);
  if (use.s_s_assemblage_ptr != NULL)
  {
    use.s_s_assemblage_in = TRUE;
    use.n_s_s_assemblage_user = i;
    save.s_s_assemblage = TRUE;
    save.n_s_s_assemblage_user = i;
    save.n_s_s_assemblage_user_end = i;
  }
  else
  {
    use.s_s_assemblage_in = FALSE;
    save.s_s_assemblage = FALSE;
  }
/*
 *   Find kinetics
 */
  if (use_kinetics == TRUE
      && (use.kinetics_ptr = kinetics_bsearch (i, &n)) != NULL)
  {
    use.n_kinetics_user = i;
    use.n_kinetics = n;
    use.kinetics_in = TRUE;
    save.kinetics = TRUE;
    save.n_kinetics_user = i;
    save.n_kinetics_user_end = i;
  }
  else
  {
    use.kinetics_ptr = NULL;
    use.kinetics_in = FALSE;
    save.kinetics = FALSE;
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
store_get_equi_reactants (int l, int kin_end)
/* ---------------------------------------------------------------------- */
{
  int i, j, k;
  static int count_pp, count_pg, count_s_s;
  static LDBLE *x0_moles;

  if (use.pp_assemblage_in == TRUE)
  {
    use.pp_assemblage_ptr = pp_assemblage_bsearch (l, &use.n_pp_assemblage);
  }
  else
    use.pp_assemblage_ptr = NULL;
  if (use.gas_phase_in == TRUE)
  {
    use.gas_phase_ptr = gas_phase_bsearch (l, &use.n_gas_phase);
  }
  else
    use.gas_phase_ptr = NULL;
  if (use.s_s_assemblage_in == TRUE)
  {
    use.s_s_assemblage_ptr =
      s_s_assemblage_bsearch (l, &use.n_s_s_assemblage);
  }
  else
    use.s_s_assemblage_ptr = NULL;

  if (kin_end == FALSE)
  {
    count_pp = count_s_s = count_pg = 0;
    if (use.pp_assemblage_ptr != NULL)
      count_pp = use.pp_assemblage_ptr->count_comps;
    if (use.gas_phase_ptr != NULL)
      count_pg = use.gas_phase_ptr->count_comps;
    if (use.s_s_assemblage_ptr != NULL)
    {
      for (j = 0; j < use.s_s_assemblage_ptr->count_s_s; j++)
      {
	for (i = 0; i < use.s_s_assemblage_ptr->s_s[j].count_comps; i++)
	{
	  count_s_s++;
	}
      }
    }
    k = count_pp + count_s_s + count_pg;
    x0_moles = NULL;
    if (k == 0)
      return (OK);
    x0_moles = (LDBLE *) PHRQ_malloc ((size_t) k * sizeof (LDBLE));
    if (x0_moles == NULL)
      malloc_error ();

    k = -1;
    for (j = 0; j < count_pp; j++)
    {
      x0_moles[++k] = use.pp_assemblage_ptr->pure_phases[j].moles;
    }
    for (j = 0; j < count_pg; j++)
    {
      x0_moles[++k] = use.gas_phase_ptr->comps[j].moles;
    }
    if (count_s_s != 0)
    {
      for (j = 0; j < use.s_s_assemblage_ptr->count_s_s; j++)
      {
	for (i = 0; i < use.s_s_assemblage_ptr->s_s[j].count_comps; i++)
	{
	  x0_moles[++k] = use.s_s_assemblage_ptr->s_s[j].comps[i].moles;
	}
/*!!!! also miscibility gap comps ??
 */
      }
    }
  }
  else
  {
    k = -1;
    for (j = 0; j < count_pp; j++)
    {
      use.pp_assemblage_ptr->pure_phases[j].moles = x0_moles[++k];
      use.pp_assemblage_ptr->pure_phases[j].delta = 0.0;
    }
    for (j = 0; j < count_pg; j++)
    {
      use.gas_phase_ptr->comps[j].moles = x0_moles[++k];
    }
    if (count_s_s != 0)
    {
      for (j = 0; j < use.s_s_assemblage_ptr->count_s_s; j++)
      {
	for (i = 0; i < use.s_s_assemblage_ptr->s_s[j].count_comps; i++)
	{
	  use.s_s_assemblage_ptr->s_s[j].comps[i].initial_moles =
	    x0_moles[++k];
	}
/*!!!! also miscibility gap comps ??
 */
      }
    }
/*
 *   This condition makes output equal for incremental_reactions TRUE/FALSE...
 *		if (incremental_reactions == FALSE || reaction_step == count_total_steps)
 */
    x0_moles = (LDBLE *) free_check_null (x0_moles);
  }
  return (OK);
}

#ifdef SKIP
/* ---------------------------------------------------------------------- */
static void
derivs (LDBLE x, LDBLE y[], LDBLE dydx[], int n_reactions, int n_user,
	struct kinetics *kinetics_ptr, LDBLE step_fraction)
/* ---------------------------------------------------------------------- */
#endif
     static void f (integertype N, realtype t, N_Vector y, N_Vector ydot,
		    void *f_data)
{
  int i, n_reactions, n_user;
  LDBLE step_fraction;
  struct kinetics *kinetics_ptr;

  cvode_error = FALSE;
  n_reactions = cvode_n_reactions;
  n_user = cvode_n_user;
  kinetics_ptr = (struct kinetics *) cvode_kinetics_ptr;
  step_fraction = cvode_step_fraction;
  rate_sim_time = cvode_rate_sim_time;

  for (i = 0; i < n_reactions; i++)
  {
    /*
       kinetics_ptr->comps[i].moles = y[i + 1];
       kinetics_ptr->comps[i].m = m_original[i] - y[i + 1];
     */
    kinetics_ptr->comps[i].moles = Ith (y, i + 1);
    kinetics_ptr->comps[i].m = m_original[i] - Ith (y, i + 1);
    if (kinetics_ptr->comps[i].m < 0)
    {
      /*
         NOTE: y is not correct if it is greater than m_original
         However, it seems to work to let y wander off, but use
         .moles as the correct integral.
         It does not work to reset Y to m_original, presumably
         because the rational extrapolation gets screwed up.
       */

      /*
         Ith(y,i + 1) = m_original[i];
       */
      kinetics_ptr->comps[i].moles = m_original[i];
      kinetics_ptr->comps[i].m = 0.0;
    }
  }
  calc_final_kinetic_reaction (kinetics_ptr);
  /*      if (set_and_run(n_user, FALSE, TRUE, n_user, step_fraction) == MASS_BALANCE) { */
  if (use.pp_assemblage_ptr != NULL)
  {
    pp_assemblage_free (use.pp_assemblage_ptr);
    pp_assemblage_copy (cvode_pp_assemblage_save, use.pp_assemblage_ptr,
			n_user);
  }
  if (use.s_s_assemblage_ptr != NULL)
  {
    s_s_assemblage_free (use.s_s_assemblage_ptr);
    s_s_assemblage_copy (cvode_s_s_assemblage_save, use.s_s_assemblage_ptr,
			 n_user);
  }

  if (set_and_run_wrapper (n_user, FALSE, TRUE, n_user, 0.0) == MASS_BALANCE)
  {
    run_reactions_iterations += iterations;
    cvode_error = TRUE;
    /*
       error_msg("Mass balance error in f", CONTINUE);
     */
    return;
  }
  if (cvode_test == TRUE)
  {
    return;
  }
  run_reactions_iterations += iterations;
  for (i = 0; i < n_reactions; i++)
  {
    kinetics_ptr->comps[i].moles = 0.0;
  }
  calc_kinetic_reaction (kinetics_ptr, 1.0);
  for (i = 0; i < n_reactions; i++)
  {
    /*
       dydx[i + 1] = kinetics_ptr->comps[i].moles;
     */
    Ith (ydot, i + 1) = kinetics_ptr->comps[i].moles;
  }
  return;
}

/*
static void Jac(integertype N, DenseMat J, RhsFn f, void *f_data, realtype t,
				N_Vector y, N_Vector fy, N_Vector ewt, realtype h,
				realtype uround, void *jac_data, long int *nfePtr,
				N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
*/
#ifdef SKIP
/* ---------------------------------------------------------------------- */
void
jacobn (LDBLE x, LDBLE y[], LDBLE dfdx[], LDBLE ** dfdy, int n_reactions,
	int n_user, struct kinetics *kinetics_ptr, LDBLE step_fraction)
/* ---------------------------------------------------------------------- */
#endif
     static void Jac (integertype N, DenseMat J, RhsFn f, void *f_data,
		      realtype t, N_Vector y, N_Vector fy, N_Vector ewt,
		      realtype h, realtype uround, void *jac_data,
		      long int *nfePtr, N_Vector vtemp1, N_Vector vtemp2,
		      N_Vector vtemp3)
{
  int count_cvode_errors;
  int i, j, n_reactions, n_user;
  LDBLE *initial_rates, del;
  struct kinetics *kinetics_ptr;
  LDBLE step_fraction;

  cvode_error = FALSE;
  n_reactions = cvode_n_reactions;
  n_user = cvode_n_user;
  kinetics_ptr = (struct kinetics *) cvode_kinetics_ptr;
  step_fraction = cvode_step_fraction;
  rate_sim_time = cvode_rate_sim_time;

  initial_rates =
    (LDBLE *) PHRQ_malloc ((size_t) n_reactions * sizeof (LDBLE));
  if (initial_rates == NULL)
    malloc_error ();

  for (i = 0; i < n_reactions; i++)
  {
    /*
       kinetics_ptr->comps[i].moles = y[i + 1];
       kinetics_ptr->comps[i].m = m_original[i] - y[i + 1];
     */
    kinetics_ptr->comps[i].moles = Ith (y, i + 1);
    kinetics_ptr->comps[i].m = m_original[i] - Ith (y, i + 1);
    if (kinetics_ptr->comps[i].m < 0)
    {
      /*
         NOTE: y is not correct if it is greater than m_original
         However, it seems to work to let y wander off, but use
         .moles as the correct integral.
         It does not work to reset Y to m_original, presumably
         because the rational extrapolation gets screwed up.
       */

      /*
         Ith(y,i + 1) = m_original[i];
       */
      kinetics_ptr->comps[i].moles = m_original[i];
      kinetics_ptr->comps[i].m = 0.0;
    }
  }
  calc_final_kinetic_reaction (kinetics_ptr);
  /* if (set_and_run(n_user, FALSE, TRUE, n_user, step_fraction) == MASS_BALANCE) { */
  if (use.pp_assemblage_ptr != NULL)
  {
    pp_assemblage_free (use.pp_assemblage_ptr);
    pp_assemblage_copy (cvode_pp_assemblage_save, use.pp_assemblage_ptr,
			n_user);
  }
  if (set_and_run_wrapper (n_user, FALSE, TRUE, n_user, 0.0) == MASS_BALANCE)
  {
    run_reactions_iterations += iterations;
    cvode_error = TRUE;
    /*
       error_msg("Mass balance error in jacobian", CONTINUE);
     */
    initial_rates = (LDBLE *) free_check_null (initial_rates);
    return;
  }
  run_reactions_iterations += iterations;
  for (i = 0; i < n_reactions; i++)
    kinetics_ptr->comps[i].moles = 0.0;
  calc_kinetic_reaction (kinetics_ptr, 1.0);
  for (i = 0; i < n_reactions; i++)
  {
    initial_rates[i] = kinetics_ptr->comps[i].moles;
  }
  for (i = 0; i < n_reactions; i++)
  {
    /* calculate reaction up to current time */
    del = 1e-12;
    cvode_error = TRUE;
    count_cvode_errors = 0;
    while (cvode_error == TRUE)
    {
      del /= 10.;
      for (j = 0; j < n_reactions; j++)
      {
	/*
	   kinetics_ptr->comps[j].moles = y[j + 1];
	   kinetics_ptr->comps[j].m = m_original[j] - y[j + 1];
	 */
	kinetics_ptr->comps[j].moles = Ith (y, j + 1);
	kinetics_ptr->comps[j].m = m_original[j] - Ith (y, j + 1);
	if (kinetics_ptr->comps[i].m < 0)
	{
	  /*
	     NOTE: y is not correct if it is greater than m_original
	     However, it seems to work to let y wander off, but use
	     .moles as the correct integral.
	     It does not work to reset Y to m_original, presumably
	     because the rational extrapolation gets screwed up.
	   */

	  /*
	     Ith(y,i + 1) = m_original[i];
	   */
	  kinetics_ptr->comps[i].moles = m_original[i];
	  kinetics_ptr->comps[i].m = 0.0;
	}
      }

      /* Add small amount of ith reaction */
      kinetics_ptr->comps[i].m -= del;
      if (kinetics_ptr->comps[i].m < 0)
      {
	kinetics_ptr->comps[i].m = 0;
      }
      kinetics_ptr->comps[i].moles += del;
      calc_final_kinetic_reaction (kinetics_ptr);
      if (use.pp_assemblage_ptr != NULL)
      {
	pp_assemblage_free (use.pp_assemblage_ptr);
	pp_assemblage_copy (cvode_pp_assemblage_save, use.pp_assemblage_ptr,
			    n_user);
      }
#ifdef SKIP
      if (set_and_run_wrapper (n_user, FALSE, TRUE, n_user, step_fraction) ==
	  MASS_BALANCE)
      {
	run_reactions_iterations += iterations;
	/*
	   error_msg("Mass balance error in jacobian 2", CONTINUE);
	 */
	cvode_error = TRUE;
	initial_rates = (LDBLE *) free_check_null (initial_rates);
	return;
      }
#endif
      if (set_and_run_wrapper (n_user, FALSE, TRUE, n_user, step_fraction) ==
	  MASS_BALANCE)
      {
	count_cvode_errors++;
	cvode_error = TRUE;
	if (count_cvode_errors > 30)
	{
	  initial_rates = (LDBLE *) free_check_null (initial_rates);
	  return;
	}
	run_reactions_iterations += iterations;
	continue;
      }
      cvode_error = FALSE;
      run_reactions_iterations += iterations;
      /*kinetics_ptr->comps[i].moles -= del; */
      for (j = 0; j < n_reactions; j++)
	kinetics_ptr->comps[j].moles = 0.0;
      calc_kinetic_reaction (kinetics_ptr, 1.0);

      /* calculate new rates for df/dy[i] */
      /* dfdx[i + 1] = 0.0; */
      for (j = 0; j < n_reactions; j++)
      {
	IJth (J, j + 1, i + 1) =
	  (kinetics_ptr->comps[j].moles - initial_rates[j]) / del;
      }
    }
  }
  for (i = 0; i < n_reactions; i++)
  {
    kinetics_ptr->comps[i].moles = 0;
  }
  initial_rates = (LDBLE *) free_check_null (initial_rates);
  return;
}

void
cvode_init (void)
{
  cvode_kinetics_ptr = NULL;
  cvode_test = 0;
  cvode_error = 0;
  cvode_n_user = -99;
  cvode_n_reactions = -99;
  cvode_step_fraction = 0.0;
  cvode_rate_sim_time = 0.0;
  cvode_rate_sim_time_start = 0.0;
  cvode_last_good_time = 0.0;
  cvode_prev_good_time = 0.0;
  cvode_last_good_y = NULL;
  cvode_prev_good_y = NULL;
  kinetics_machEnv = NULL;
  kinetics_y = kinetics_abstol = NULL;
  kinetics_cvode_mem = NULL;
  cvode_pp_assemblage_save = NULL;
  cvode_s_s_assemblage_save = NULL;
  return;
}
