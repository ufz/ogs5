#define EXTERNAL extern
#include "global.h"
#include "phqalloc.h"
#include "output.h"
#include "phrqproto.h"

static char const svnid[] =
  "$Id: inverse.c 4 2009-04-21 17:29:29Z delucia $";

#define MAX_MODELS 20
#define MIN_TOTAL_INVERSE 1e-14

#ifdef INVERSE_CL1MP
/* cl1mp.c */
extern int cl1mp (int k, int l, int m, int n,
		  int nklmd, int n2d,
		  LDBLE * q_arg,
		  int *kode, LDBLE toler,
		  int *iter, LDBLE * x_arg, LDBLE * res_arg, LDBLE * error,
		  LDBLE * cu_arg, int *iu, int *s, int check,
		  LDBLE censor_arg);
#endif

/* variables local to module */

int max_row_count, max_column_count;
static int carbon;
char **col_name, **row_name;
static int count_rows, count_optimize;
static int col_phases, col_redox, col_epsilon, col_ph, col_water,
  col_isotopes, col_phase_isotopes;
static int row_mb, row_fract, row_charge, row_carbon, row_isotopes,
  row_epsilon, row_isotope_epsilon, row_water;
static LDBLE *zero, *array1, *res, *delta1, *delta2, *delta3, *cu,
  *delta_save;
static LDBLE *min_delta, *max_delta;
static int *iu, *is;
static int klmd, nklmd, n2d, kode, iter;
static LDBLE toler, error, max_pct, scaled_error;
static struct master *master_alk;
int *row_back, *col_back;

static unsigned long *good, *bad, *minimal;
static int max_good, max_bad, max_minimal;
static int count_good, count_bad, count_minimal, count_calls;
static unsigned long soln_bits, phase_bits, current_bits, temp_bits;

/* subroutines */

static int add_to_file(char *filename, char *string);
static int bit_print (unsigned long bits, int l);
static int carbon_derivs (struct inverse *inv_ptr);
static int check_isotopes (struct inverse *inv_ptr);
static int check_solns (struct inverse *inv_ptr);
static int count_isotope_unknowns (struct inverse *inv_ptr,
				   struct isotope **isotope_unknowns);
struct isotope *get_isotope (struct solution *solution_ptr, char *elt);
static struct conc *get_total (struct solution *solution_ptr, char *elt);
static int isotope_balance_equation (struct inverse *inv_ptr, int row, int n);
static int post_mortem (void);
static unsigned long get_bits (unsigned long bits, int position, int number);
static unsigned long minimal_solve (struct inverse *inv_ptr,
				    unsigned long minimal_bits);
static void dump_netpath(struct inverse *inv_ptr);
static int dump_netpath_pat (struct inverse *inv_ptr);
static int next_set_phases (struct inverse *inv_ptr, int first_of_model_size,
			    int model_size);
static int phase_isotope_inequalities (struct inverse *inv_ptr);
static int print_model (struct inverse *inv_ptr);
static int punch_model_heading (struct inverse *inv_ptr);
static int punch_model (struct inverse *inv_ptr);
void print_isotope (FILE *netpath_file, struct solution *solution_ptr, char *elt, char *string);
void print_total (FILE *netpath_file, struct solution *solution_ptr, char *elt, char *string);
void print_total_multi (FILE *netpath_file, struct solution *solution_ptr, char *string, char *elt0, char *elt1, char *elt2, char *elt3, char *elt4 );
void print_total_pat (FILE *netpath_file, char *elt, char *string);
static int range (struct inverse *inv_ptr, unsigned long cur_bits);
static int save_bad (unsigned long bits);
static int save_good (unsigned long bits);
static int save_minimal (unsigned long bits);
static unsigned long set_bit (unsigned long bits, int position, int value);
static int setup_inverse (struct inverse *inv_ptr);
int set_initial_solution (int n_user_old, int n_user_new);
static int set_ph_c (struct inverse *inv_ptr,
		     int i,
		     struct solution *soln_ptr_orig,
		     int n_user_new,
		     LDBLE d_alk, LDBLE ph_factor, LDBLE alk_factor);
static int shrink (struct inverse *inv_ptr, LDBLE * array_in,
		   LDBLE * array_out, int *k, int *l, int *m, int *n,
		   unsigned long cur_bits, LDBLE * delta_l, int *col_back_l,
		   int *row_back_l);
static int solve_inverse (struct inverse *inv_ptr);
static int solve_with_mask (struct inverse *inv_ptr, unsigned long cur_bits);
static int subset_bad (unsigned long bits);
static int subset_minimal (unsigned long bits);
static int superset_minimal (unsigned long bits);
static int write_optimize_names (struct inverse *inv_ptr);

#ifdef SKIP
#define SCALE_EPSILON .0009765625
#define SCALE_WATER   .0009765625
#define SCALE_ALL     16
#endif
#define SCALE_EPSILON .0009765625
#define SCALE_WATER   1.
#define SCALE_ALL     1.

static FILE *netpath_file;
static int count_inverse_models, count_pat_solutions;
/* ---------------------------------------------------------------------- */
int
inverse_models (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Go through list of inverse models, make calculations
 *   for any marked "new".
 */
  int n, print1;
  char string[MAX_LENGTH];

  if (svnid == NULL)
    fprintf (stderr, " ");
  array1 = NULL;
  zero = NULL;
  res = NULL;
  delta1 = NULL;
  delta2 = NULL;
  delta3 = NULL;
  delta_save = NULL;
  cu = NULL;
  iu = NULL;
  is = NULL;
  col_name = NULL;
  row_name = NULL;
  min_delta = NULL;
  max_delta = NULL;
  good = NULL;
  bad = NULL;
  minimal = NULL;

  print1 = TRUE;
  state = INVERSE;
  dl_type_x = NO_DL;

  for (n = 0; n < count_inverse; n++)
  {
    if (inverse[n].new_def == TRUE)
    {
/*
 * dump .lon file
 */
      if (inverse[n].netpath != NULL) dump_netpath(&inverse[n]);

/*
 * open .pat file
 */
      if (inverse[n].pat != NULL) 
      {
	strcpy(string, inverse[n].pat);
	if (replace(".pat", ".pat", string) != TRUE)
	{
	  strcat(string, ".pat");
	}
	netpath_file = fopen(string, "w");
	if (netpath_file == NULL) {
	  sprintf (error_string, "Can't open file, %s.", string);
	  error_msg (error_string, STOP);
	}
	count_inverse_models = 0;
	count_pat_solutions = 0;
	/* Header */
	fprintf(netpath_file, "2.14               # File format\n");
      }
/*
 *  Fill in stucture "use".  
 */
      use.inverse_in = TRUE;
      use.inverse_ptr = &inverse[n];
      use.n_inverse_user = inverse[n].n_user;
/*
 *  Initial prints
 */
      sprintf (error_string, "Beginning of inverse modeling %d calculations.",
	       inverse[n].n_user);
      dup_print (error_string, TRUE);

      if (inverse[n].mp == TRUE)
      {
	output_msg (OUTPUT_MESSAGE,
		    "Using Cl1MP multiprecision optimization routine.\n");

      }
      else
      {
	output_msg (OUTPUT_MESSAGE,
		    "Using Cl1 standard precision optimization routine.\n");
      }
      status (0, NULL);

/*
 *  Setup and solve
 */
      count_calls = 0;
      setup_inverse (&(inverse[n]));
      punch_model_heading (&inverse[n]);
      solve_inverse (&(inverse[n]));
      if (inverse[n].count_isotope_unknowns > 0)
      {
	inverse[n].isotope_unknowns =
	  (struct isotope *) free_check_null (inverse[n].isotope_unknowns);
      }
      inverse[n].new_def = FALSE;
      if (inverse[n].pat != NULL) 
      {
	fclose(netpath_file);
	netpath_file = NULL;
      }
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
setup_inverse (struct inverse *inv_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Fill in array for an inverse problem
 */
  int i, j, k, i_alk, i_carb;
  int max;
  int count_rows_t;
  int column, row;
  int temp, temppun;
  LDBLE isotope_number;
  LDBLE f, coef, cb, conc;
  char token[MAX_LENGTH];
  struct phase *phase_ptr;
  struct solution *solution_ptr;
  struct reaction *rxn_ptr;
  struct master *master_ptr;
/*
 *   Determine array sizes, row and column positions
 */
  toler = inv_ptr->tolerance;
  if (inv_ptr->mp == TRUE)
  {
    toler = inv_ptr->mp_tolerance;
  }
/*
 *   Alkalinity derivatives with pH and carbon
 */
  carbon = 1;
  temp = pr.status;
  pr.status = FALSE;
  temppun = punch.in;
  punch.in = FALSE;
  carbon_derivs (inv_ptr);
  pr.status = temp;
  punch.in = temppun;
  state = INVERSE;
/* 
 *   tidy isotopes if necessary
 */
  inv_ptr->count_isotope_unknowns = 0;
  if (inv_ptr->count_isotopes > 0)
  {
    inv_ptr->count_isotope_unknowns =
      count_isotope_unknowns (inv_ptr, &inv_ptr->isotope_unknowns);
    if (input_error > 0)
    {
      error_msg ("Stopping because of input errors.", STOP);
    }
    check_isotopes (inv_ptr);
    if (input_error > 0)
    {
      error_msg ("Stopping because of input errors.", STOP);
    }
  }

/*
 *    count unknowns
 */
  max_column_count = inv_ptr->count_elts * inv_ptr->count_solns +	/* epsilons */
    inv_ptr->count_solns +	/* solutions */
    inv_ptr->count_phases +	/* phases */
    inv_ptr->count_redox_rxns +	/* redox reactions */
    carbon * inv_ptr->count_solns +	/* pH */
    1 +				/* water */
    inv_ptr->count_isotope_unknowns * inv_ptr->count_solns +	/* isotopes in solution */
    inv_ptr->count_isotopes * inv_ptr->count_phases +	/* isotopes in phases */
    1 + 1;			/* rhs, ineq */
  count_unknowns = max_column_count - 2;
  col_phases = inv_ptr->count_solns;
  col_redox = col_phases + inv_ptr->count_phases;
  col_epsilon = col_redox + inv_ptr->count_redox_rxns;
  col_ph = col_epsilon + inv_ptr->count_elts * inv_ptr->count_solns;
  col_water = col_ph + carbon * inv_ptr->count_solns;
  col_isotopes = col_water + 1;
  col_phase_isotopes =
    col_isotopes + inv_ptr->count_isotope_unknowns * inv_ptr->count_solns;
  max_row_count = inv_ptr->count_solns * inv_ptr->count_elts +	/* optimize */
    carbon * inv_ptr->count_solns +	/* optimize ph */
    1 +				/* optimize water */
    inv_ptr->count_solns * inv_ptr->count_isotope_unknowns +	/* optimize isotopes */
    inv_ptr->count_isotopes * inv_ptr->count_phases +	/* optimize phase isotopes */
    inv_ptr->count_elts +	/* mass balances */
    1 + 1 +			/* fractions, init and final */
    inv_ptr->count_solns +	/* charge balances */
    carbon * inv_ptr->count_solns +	/* dAlk = dC + dph */
    inv_ptr->count_isotopes +	/* isotopes */
    2 * inv_ptr->count_solns * inv_ptr->count_elts +	/* epsilon constraints */
    2 * carbon * inv_ptr->count_solns +	/* epsilon on ph */
    2 +				/* epsilon for water */
    2 * inv_ptr->count_isotope_unknowns * inv_ptr->count_solns +	/* epsilon for isotopes */
    2 * inv_ptr->count_isotopes * inv_ptr->count_phases +	/* epsilon for isotopes in phases */
    2;				/* work space */

  row_mb = inv_ptr->count_solns * inv_ptr->count_elts +
    carbon * inv_ptr->count_solns + 1 +
    inv_ptr->count_solns * inv_ptr->count_isotope_unknowns +
    inv_ptr->count_isotopes * inv_ptr->count_phases;
  row_fract = row_mb + inv_ptr->count_elts;
  row_charge = row_fract + 2;
  row_carbon = row_charge + inv_ptr->count_solns;
  row_isotopes = row_carbon + carbon * inv_ptr->count_solns;
  row_epsilon = row_isotopes + inv_ptr->count_isotopes;
/*   The next three are not right, some rows of epsilon are deleted */
/*
	row_ph_epsilon = row_epsilon + 2 * inv_ptr->count_solns * inv_ptr->count_elts; 
	row_water_epsilon = row_ph + 2 * carbon * inv_ptr->count_solns;
	row_isotope_epsilon
	row_isotope_phase_epsilon
 */
/*
 *   Malloc space for arrays
 */
  array = (LDBLE *) free_check_null (array);
  array =
    (LDBLE *) PHRQ_malloc ((size_t) max_column_count * max_row_count *
			   sizeof (LDBLE));
  if (array == NULL)
    malloc_error ();

  array1 =
    (LDBLE *) PHRQ_malloc ((size_t) max_column_count * max_row_count *
			   sizeof (LDBLE));
  if (array1 == NULL)
    malloc_error ();

  col_name =
    (char **) PHRQ_malloc ((size_t) max_column_count * sizeof (char *));
  if (col_name == NULL)
    malloc_error ();

  row_name = (char **) PHRQ_malloc ((size_t) max_row_count * sizeof (char *));
  if (row_name == NULL)
    malloc_error ();

  delta = (LDBLE *) free_check_null (delta);
  delta = (LDBLE *) PHRQ_malloc ((size_t) max_column_count * sizeof (LDBLE));
  if (delta == NULL)
    malloc_error ();

  delta1 = (LDBLE *) PHRQ_malloc ((size_t) max_column_count * sizeof (LDBLE));
  if (delta1 == NULL)
    malloc_error ();

  delta2 = (LDBLE *) PHRQ_malloc ((size_t) max_column_count * sizeof (LDBLE));
  if (delta2 == NULL)
    malloc_error ();

  delta3 = (LDBLE *) PHRQ_malloc ((size_t) max_column_count * sizeof (LDBLE));
  if (delta3 == NULL)
    malloc_error ();

  delta_save =
    (LDBLE *) PHRQ_malloc ((size_t) max_column_count * sizeof (LDBLE));
  if (delta_save == NULL)
    malloc_error ();

  min_delta =
    (LDBLE *) PHRQ_malloc ((size_t) max_column_count * sizeof (LDBLE));
  if (min_delta == NULL)
    malloc_error ();

  max_delta =
    (LDBLE *) PHRQ_malloc ((size_t) max_column_count * sizeof (LDBLE));
  if (max_delta == NULL)
    malloc_error ();

  res = (LDBLE *) PHRQ_malloc ((size_t) max_row_count * sizeof (LDBLE));
  if (res == NULL)
    malloc_error ();

  if (max_column_count < max_row_count)
  {
    max = max_row_count;
  }
  else
  {
    max = max_column_count;
  }
  zero = (LDBLE *) PHRQ_malloc ((size_t) max * sizeof (LDBLE));
  if (zero == NULL)
    malloc_error ();
/*
 *   Define zero and zero array, delta
 */
  for (i = 0; i < max; i++)
    zero[i] = 0.0;

  memcpy ((void *) &(delta[0]), (void *) &(zero[0]),
	  (size_t) max_column_count * sizeof (LDBLE));
  memcpy ((void *) &(min_delta[0]), (void *) &(zero[0]),
	  (size_t) max_column_count * sizeof (LDBLE));
  memcpy ((void *) &(max_delta[0]), (void *) &(zero[0]),
	  (size_t) max_column_count * sizeof (LDBLE));
  for (i = 0; i < max_row_count; i++)
  {
    memcpy ((void *) &(array[i * max_column_count]), (void *) &(zero[0]),
	    (size_t) max_column_count * sizeof (LDBLE));
  }
/*
 *   begin filling array
 */
  count_rows = 0;
/*
 *   optimization
 */
  count_optimize = inv_ptr->count_solns * inv_ptr->count_elts +	/* optimize */
    carbon * inv_ptr->count_solns +	/* optimize ph */
    1 +				/* optimize water */
    inv_ptr->count_solns * inv_ptr->count_isotope_unknowns +	/* optimize isotopes */
    inv_ptr->count_isotopes * inv_ptr->count_phases;	/* optimize phase isotopes */

  for (i = 0; i < count_optimize; i++)
  {
    row_name[count_rows] = string_hsave ("optimize");
    count_rows++;
  }
  write_optimize_names (inv_ptr);
/*
 *   equalities
 */

/*
 *   Mass_balance: solution data
 */

  /* initialize master species */
  for (i = 0; i < count_master; i++)
  {
    master[i]->in = -1;
    if (strstr (master[i]->elt->name, "Alk") == master[i]->elt->name)
    {
      master_alk = master[i];
    }
  }
  /* mark master species included in model, write row names */
  count_rows_t = count_rows;
  i_alk = -1;
  i_carb = -1;
  for (i = 0; i < inv_ptr->count_elts; i++)
  {
    master_ptr = inv_ptr->elts[i].master;
    if (master_ptr == master_alk)
      i_alk = i;
    if (strcmp (master_ptr->elt->name, "C(4)") == 0)
      i_carb = i;
    inv_ptr->elts[i].master->in = count_rows_t;
    row_name[count_rows_t] = inv_ptr->elts[i].master->elt->name;
    count_rows_t++;
  }
  /* put concentrations in array */
  for (i = 0; i < inv_ptr->count_solns; i++)
  {
    xsolution_zero ();
    solution_ptr = solution_bsearch (inv_ptr->solns[i], &j, TRUE);
    if (solution_ptr == NULL)
    {
      sprintf (error_string, "Solution number %d not found.",
	       inv_ptr->solns[i]);
      error_msg (error_string, STOP);
    }
    /* write master species concentrations */
    for (j = 0; solution_ptr->totals[j].description != NULL; j++)
    {
      master_ptr = master_bsearch (solution_ptr->totals[j].description);
      master_ptr->total += solution_ptr->totals[j].moles;
      /*   List elements not included in model */
      if (master_ptr->in < 0)
      {
	sprintf (error_string,
		 "%s is included in solution %d, but is not included as a mass-balance constraint.",
		 solution_ptr->totals[j].description, inv_ptr->solns[i]);
	warning_msg (error_string);
      }
    }
    master_alk->total = solution_ptr->total_alkalinity;
    f = 1.0;
    if (i == (inv_ptr->count_solns - 1))
    {
      f = -1.0;
    }
    column = i;
    sprintf (token, "soln %d", i);
    col_name[column] = string_hsave (token);
    for (j = 0; j < count_master; j++)
    {
      if (master[j]->in >= 0)
      {
	array[master[j]->in * max_column_count + i] = f * master[j]->total;
	if (master[j]->s == s_eminus)
	{
	  array[master[j]->in * max_column_count + i] = 0.0;
	}
      }
    }
    /* calculate charge balance for elements in model */
    cb = 0;
    for (j = 0; j < count_master; j++)
    {
      if (master[j]->in >= 0)
      {
	if (master[j]->s == s_eminus)
	{
	  coef = 0.0;
	}
	else if (master[j] == master_alk)
	{
	  coef = -1.0;
	}
	else
	{
	  coef = master[j]->s->z + master[j]->s->alk;
	}
	cb += coef * master[j]->total;
      }
    }
    if (fabs (cb) < toler)
      cb = 0.0;
    array[(row_charge + i) * max_column_count + i] = cb;
  }

/*   mass_balance: phase data */

  for (i = 0; i < inv_ptr->count_phases; i++)
  {
    phase_ptr = inv_ptr->phases[i].phase;
    rxn_ptr = phase_ptr->rxn_s;
    column = col_phases + i;
    col_name[column] = phase_ptr->name;
    for (j = 1; rxn_ptr->token[j].s != NULL; j++)
    {
      if (rxn_ptr->token[j].s->secondary != NULL)
      {
	master_ptr = rxn_ptr->token[j].s->secondary;
      }
      else
      {
	master_ptr = rxn_ptr->token[j].s->primary;
      }
      if (master_ptr == NULL)
      {
	sprintf (error_string, "Setup_inverse, reaction for phase, %s.",
		 phase_ptr->name);
	error_msg (error_string, STOP);
      }
      if (master_ptr->s == s_hplus)
	continue;
      if (master_ptr->s == s_h2o)
      {
	row = row_fract;
/* turn off h2o from minerals in water mass balance */
	if (inv_ptr->mineral_water != TRUE)
	  continue;

      }
      else
      {
	row = master_ptr->in;
      }
      /* e- has coef of 0 for some reason */
      coef = master_ptr->coef;
      if (coef <= 0)
	coef = 1.0;
      array[row * max_column_count + column] = rxn_ptr->token[j].coef * coef;
    }
    row = master_alk->in;	/* include alkalinity for phase */
    array[row * max_column_count + column] = calc_alk (rxn_ptr);
  }

/*   mass balance: redox reaction data */

  k = 0;
  for (i = 0; i < inv_ptr->count_elts; i++)
  {
    if (inv_ptr->elts[i].master->s->primary == NULL)
    {
      coef = inv_ptr->elts[i].master->coef;
      rxn_ptr = inv_ptr->elts[i].master->rxn_primary;
      column = col_redox + k;
      col_name[column] = inv_ptr->elts[i].master->elt->name;
      k++;
      for (j = 0; rxn_ptr->token[j].s != NULL; j++)
      {
	if (rxn_ptr->token[j].s->secondary != NULL)
	{
	  master_ptr = rxn_ptr->token[j].s->secondary;
	}
	else
	{
	  master_ptr = rxn_ptr->token[j].s->primary;
	}
	if (master_ptr == NULL)
	{
	  sprintf (error_string,
		   "Subroutine setup_inverse, element not found, %s.",
		   rxn_ptr->token[j].s->name);
	  error_msg (error_string, STOP);
	}
	if (master_ptr->s == s_hplus)
	  continue;
	if (master_ptr->s == s_h2o)
	{
	  row = row_fract;
/* turn off h2o from minerals in water mass balance */
	  if (inv_ptr->mineral_water != TRUE)
	    continue;
	}
	else
	{
	  row = master_ptr->in;
	}
	array[row * max_column_count + column] = rxn_ptr->token[j].coef;
	/* if coefficient of element is not 1.0 in master species */
	if (j != 0)
	  array[row * max_column_count + column] /= coef;
      }
      row = master_alk->in;	/* include alkalinity for redox reaction */
      array[row * max_column_count + column] =
	(calc_alk (rxn_ptr) - inv_ptr->elts[i].master->s->alk) / coef;
    }
  }

/*   mass-balance: epsilons */

  column = col_epsilon;
  for (i = 0; i < inv_ptr->count_elts; i++)
  {
    row = inv_ptr->elts[i].master->in;
    for (j = 0; j < inv_ptr->count_solns; j++)
    {
      if (j < (inv_ptr->count_solns - 1))
      {
	array[row * max_column_count + column] = 1.0;
      }
      else
      {
	array[row * max_column_count + column] = -1.0;
      }
      if (inv_ptr->elts[i].master->s == s_eminus)
      {
	array[row * max_column_count + column] = 0.0;
      }
      sprintf (token, "%s %d", row_name[row], j);
      col_name[column] = string_hsave (token);
      column++;
    }
  }
  count_rows += inv_ptr->count_elts;

/*   put names in col_name for ph */

  for (i = 0; i < inv_ptr->count_solns; i++)
  {
    sprintf (token, "ph %d", i);
    col_name[column] = string_hsave (token);
    column++;
  }
/*   put names in col_name for water */

  sprintf (token, "water");
  col_name[column] = string_hsave (token);
  column++;

/*   put names of isotopes in col_name */
  for (i = 0; i < inv_ptr->count_solns; i++)
  {
    for (j = 0; j < inv_ptr->count_isotope_unknowns; j++)
    {
      sprintf (token, "%d%s %d",
	       (int) inv_ptr->isotope_unknowns[j].isotope_number,
	       inv_ptr->isotope_unknowns[j].elt_name, i);
      col_name[column] = string_hsave (token);
      column++;
    }
  }

/*   put phase isotopes in col_name */

  if (inv_ptr->count_isotopes > 0)
  {
    /* isotopes of phases phases */
    for (i = 0; i < inv_ptr->count_phases; i++)
    {
      for (j = 0; j < inv_ptr->count_isotopes; j++)
      {
	sprintf (token, "%d%s %s",
		 (int) inv_ptr->isotopes[j].isotope_number,
		 inv_ptr->isotopes[j].elt_name,
		 inv_ptr->phases[i].phase->name);
	col_name[column] = string_hsave (token);
	column++;
      }
    }
  }
/*   
 *   Initial solution mixing fractions or water mass balance
 */
  for (i = 0; i < inv_ptr->count_solns; i++)
  {
    solution_ptr = solution_bsearch (inv_ptr->solns[i], &j, TRUE);
    if (i < inv_ptr->count_solns - 1)
    {
      array[count_rows * max_column_count + i] =
	1.0 / gfw_water * solution_ptr->mass_water;
    }
    else
    {
      array[count_rows * max_column_count + inv_ptr->count_solns - 1] =
	-1.0 / gfw_water * solution_ptr->mass_water;
    }
  }
  /* coefficient for water uncertainty */
  if (inv_ptr->water_uncertainty > 0)
  {
    array[count_rows * max_column_count + col_water] = 1.0;
  }
  row_name[count_rows] = string_hsave ("H2O");
  row_water = count_rows;
  count_rows++;

/*
 *   Final solution fraction equals 1.0
 */

  array[count_rows * max_column_count + inv_ptr->count_solns - 1] = 1.0;
  array[count_rows * max_column_count + count_unknowns] = 1.0;
  row_name[count_rows] = string_hsave ("fract, final");
  count_rows++;

/*
 *   Charge balance:
 */

  for (i = 0; i < inv_ptr->count_solns; i++)
  {
/*		solution_ptr = solution_bsearch(inv_ptr->solns[i], &j, TRUE); */
/*		array[count_rows * max_column_count + i] = solution_ptr->cb; */
    for (j = 0; j < inv_ptr->count_elts; j++)
    {
      column = col_epsilon + j * inv_ptr->count_solns + i;
      coef = inv_ptr->elts[j].master->s->z + inv_ptr->elts[j].master->s->alk;
      if (inv_ptr->elts[j].master == master_alk)
      {
	coef = -1.0;
      }
      array[count_rows * max_column_count + column] = coef;
      if (inv_ptr->elts[j].master->s == s_eminus)
      {
	array[count_rows * max_column_count + column] = 0.0;
      }
    }
    sprintf (token, "%s %d", "charge", i);
    row_name[count_rows] = string_hsave (token);
    count_rows++;
  }
/*
 *   dC = (dC/dph)*dph + (dC/dAlk)*dAlk for each solution
 */
/*
 *   dAlk = (dAlk/dC)*dC + (dAlk/dpH)*dpH for each solution
 */
  for (i = 0; i < inv_ptr->count_solns; i++)
  {
    if (inv_ptr->dalk_dph[i] != 0 || inv_ptr->dalk_dc[i] != 0)
    {
      column = col_ph + i;
      array[count_rows * max_column_count + column] = inv_ptr->dalk_dph[i];
      column = col_epsilon + i_alk * inv_ptr->count_solns + i;
      array[count_rows * max_column_count + column] = -1.0;
      column = col_epsilon + i_carb * inv_ptr->count_solns + i;
      array[count_rows * max_column_count + column] = inv_ptr->dalk_dc[i];
    }
    sprintf (token, "%s %d", "dAlk", i);
    row_name[count_rows] = string_hsave (token);
    count_rows++;
  }
/*
 *   Isotope mass balances
 */
  if (input_error > 0)
  {
    error_msg ("Stopping because of input errors.", STOP);
  }
  if (inv_ptr->count_isotopes != 0)
  {
    for (j = 0; j < inv_ptr->count_isotopes; j++)
    {
      isotope_balance_equation (inv_ptr, count_rows, j);
      sprintf (token, "%d%s", (int) inv_ptr->isotopes[j].isotope_number,
	       inv_ptr->isotopes[j].elt_name);
      row_name[count_rows] = string_hsave (token);
      count_rows++;
    }
  }
/*
 *   inequalities
 */
  row_epsilon = count_rows;
  for (i = 0; i < inv_ptr->count_solns; i++)
  {
    for (j = 0; j < inv_ptr->count_elts; j++)
    {
      if (inv_ptr->elts[j].master->s == s_eminus)
	continue;
      column = col_epsilon + j * inv_ptr->count_solns + i;

/* calculate magnitude of bound */

      coef = inv_ptr->elts[j].uncertainties[i];
      if (coef <= 0.0)
      {
	coef = -coef;
      }
      else
      {
	coef =
	  array[inv_ptr->elts[j].master->in * max_column_count + i] * coef;
	coef = fabs (coef);
      }

      if (coef < toler)
	coef = 0;

/* zero column if uncertainty is zero */
      if (coef == 0.0)
      {
	for (k = 0; k < count_rows; k++)
	{
	  array[k * max_column_count + column] = 0.0;
	}
	continue;
      }

/*  this statement probably obviates some of the following logic. */
/*			coef += toler; */

/* scale epsilon optimization equation */

      if (coef < toler)
      {
	array[(column - col_epsilon) * max_column_count + column] =
	  SCALE_EPSILON / toler;
      }
      else
      {
	array[(column - col_epsilon) * max_column_count + column] =
	  SCALE_EPSILON / coef;
      }

/* set upper limit of change in positive direction */
      if (coef < toler)
      {
	coef = toler;
	f = 10;
      }
      else
      {
	f = 1.0;
      }
      array[count_rows * max_column_count + column] = 1.0 * f;
      array[count_rows * max_column_count + i] = -coef * f;
      sprintf (token, "%s %s", inv_ptr->elts[j].master->elt->name, "eps+");
      row_name[count_rows] = string_hsave (token);
      count_rows++;

/* set lower limit of change in negative direction */
      conc = array[inv_ptr->elts[j].master->in * max_column_count + i];

      /* if concentration is zero, only positive direction allowed */
      if (conc == 0.0)
      {
	delta[column] = 1.0;
	continue;
      }
      /* if uncertainty is less than tolerance, set uncertainty to toler */
      if (coef <= toler)
      {
/*				f = 10 * toler / coef; */
	coef = toler;
	f = 10;
      }
      else
      {
	f = 1.0;
      }
      /* if uncertainty is greater than concentration,
         maximum negative is equal to concentrations,
         except alkalinity */
      if (coef > fabs (conc) &&
	  (strstr (inv_ptr->elts[j].master->elt->name, "Alkalinity") !=
	   inv_ptr->elts[j].master->elt->name))
	coef = fabs (conc) + toler;

      array[count_rows * max_column_count + i] = -coef * f;
      array[count_rows * max_column_count + column] = -1.0 * f;
      sprintf (token, "%s %s", inv_ptr->elts[j].master->elt->name, "eps-");
      row_name[count_rows] = string_hsave (token);
      count_rows++;
    }
  }
/*
 *   inequalities for pH
 */
  /* row_ph_epsilon = count_rows; */
  if (inv_ptr->carbon == TRUE)
  {
    for (i = 0; i < inv_ptr->count_solns; i++)
    {
      column = col_ph + i;
      coef = inv_ptr->ph_uncertainties[i];

/* scale epsilon in optimization equation */

      array[(column - col_epsilon) * max_column_count + column] =
	SCALE_EPSILON / coef;

/* set upper limit of change in positive direction */

      array[count_rows * max_column_count + column] = 1.0;
      array[count_rows * max_column_count + i] = -coef;
      sprintf (token, "%s %s", "pH", "eps+");
      row_name[count_rows] = string_hsave (token);
      count_rows++;

/* set lower limit of change in negative direction */

      array[count_rows * max_column_count + column] = -1.0;
      array[count_rows * max_column_count + i] = -coef;
      sprintf (token, "%s %s", "pH", "eps-");
      row_name[count_rows] = string_hsave (token);
      count_rows++;
    }
  }
/*
 *   inequalities for water
 */
  column = col_water;
  coef = inv_ptr->water_uncertainty;
  /* row_water_epsilon = count_rows; */
  if (coef > 0.0)
  {
/* set upper limit of change in positive direction */
    array[count_rows * max_column_count + column] = 1.0;
    array[count_rows * max_column_count + count_unknowns] = coef;
    sprintf (token, "%s %s", "water", "eps+");
    row_name[count_rows] = string_hsave (token);
    count_rows++;

/* set lower limit of change in negative direction */

    array[count_rows * max_column_count + column] = -1.0;
    array[count_rows * max_column_count + count_unknowns] = coef;
    sprintf (token, "%s %s", "water", "eps-");
    row_name[count_rows] = string_hsave (token);
    count_rows++;
  }
/*
 *   inequalities for isotopes
 */
  row_isotope_epsilon = count_rows;
  if (inv_ptr->count_isotopes > 0)
  {
    for (i = 0; i < inv_ptr->count_solns; i++)
    {
      solution_ptr = solution_bsearch (inv_ptr->solns[i], &k, TRUE);
      for (j = 0; j < inv_ptr->count_isotope_unknowns; j++)
      {
	column = col_isotopes + (i * inv_ptr->count_isotope_unknowns) + j;
	master_ptr = inv_ptr->isotope_unknowns[j].master;
	isotope_number = inv_ptr->isotope_unknowns[j].isotope_number;
	for (k = 0; k < solution_ptr->count_isotopes; k++)
	{
	  if (solution_ptr->isotopes[k].master == master_ptr &&
	      solution_ptr->isotopes[k].isotope_number == isotope_number)
	  {
	    coef = solution_ptr->isotopes[k].x_ratio_uncertainty;

/* scale epsilon in optimization equation */

	    array[(column - col_epsilon) * max_column_count + column] =
	      SCALE_EPSILON / coef;

/* set upper limit of change in positive direction */
	    array[count_rows * max_column_count + column] = 1.0;
	    array[count_rows * max_column_count + i] = -coef;
	    sprintf (token, "%d%s %s",
		     (int) solution_ptr->isotopes[k].isotope_number,
		     solution_ptr->isotopes[k].elt_name, "eps+");
	    row_name[count_rows] = string_hsave (token);
	    count_rows++;

/* set lower limit of change in negative direction */

	    array[count_rows * max_column_count + column] = -1.0;
	    array[count_rows * max_column_count + i] = -coef;
	    sprintf (token, "%d%s %s",
		     (int) solution_ptr->isotopes[k].isotope_number,
		     solution_ptr->isotopes[k].elt_name, "eps-");
	    row_name[count_rows] = string_hsave (token);
	    count_rows++;
	    break;
	  }
	}
      }
    }
  }
/*
 *   inequalities for isotopes in phases
 */
  /*      row_isotope_phase_epsilon = count_rows; */
  phase_isotope_inequalities (inv_ptr);
  if (input_error > 0)
  {
    error_msg ("Stopping because of input errors.", STOP);
  }
/*
 *   Set non-negativity constraints
 */

  for (i = 0; i < inv_ptr->count_phases; i++)
  {
    if (inv_ptr->phases[i].constraint == PRECIPITATE)
    {
      delta[col_phases + i] = -1.0;
    }
    else if (inv_ptr->phases[i].constraint == DISSOLVE)
    {
      delta[col_phases + i] = 1.0;
    }
  }
  for (i = 0; i < (inv_ptr->count_solns - 1); i++)
  {
    delta[i] = 1.0;
  }
/*
 *   Scale water equation
 */
  for (i = 0; i < max_column_count; i++)
  {
    array[row_water * max_column_count + i] *= SCALE_WATER;
  }
/*
 *   Arrays are complete
 */
  if (debug_inverse == TRUE)
  {
    for (i = 0; i < count_unknowns; i++)
    {
      output_msg (OUTPUT_MESSAGE, "%d\t%s\n", i, col_name[i]);
    }
    for (i = 0; i < count_rows; i++)
    {
      k = 0;
      output_msg (OUTPUT_MESSAGE, "%d\t%s\n", i, row_name[i]);
      for (j = 0; j < count_unknowns + 1; j++)
      {
	if (k > 7)
	{
	  output_msg (OUTPUT_MESSAGE, "\n");
	  k = 0;
	}
	output_msg (OUTPUT_MESSAGE, "%11.2e",
		    (double) array[i * max_column_count + j]);
	k++;
      }
      if (k != 0)
      {
	output_msg (OUTPUT_MESSAGE, "\n");
      }
      output_msg (OUTPUT_MESSAGE, "\n");
    }
    output_msg (OUTPUT_MESSAGE, "row_mb %d\n", row_mb);
    output_msg (OUTPUT_MESSAGE, "row_fract %d\n", row_fract);
    output_msg (OUTPUT_MESSAGE, "row_charge %d\n", row_charge);
    output_msg (OUTPUT_MESSAGE, "row_carbon %d\n", row_carbon);
    output_msg (OUTPUT_MESSAGE, "row_isotopes %d\n", row_isotopes);
    output_msg (OUTPUT_MESSAGE, "row_epsilon %d\n", row_epsilon);

    output_msg (OUTPUT_MESSAGE, "col_phases %d\n", col_phases);
    output_msg (OUTPUT_MESSAGE, "col_redox %d\n", col_redox);
    output_msg (OUTPUT_MESSAGE, "col_epsilon %d\n", col_epsilon);
    output_msg (OUTPUT_MESSAGE, "col_ph %d\n", col_ph);
    output_msg (OUTPUT_MESSAGE, "col_water %d\n", col_water);
    output_msg (OUTPUT_MESSAGE, "col_isotopes %d\n", col_isotopes);
    output_msg (OUTPUT_MESSAGE, "col_phase_isotopes %d\n",
		col_phase_isotopes);
    output_msg (OUTPUT_MESSAGE, "count_unknowns %d\n", count_unknowns);
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
solve_inverse (struct inverse *inv_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Exhaustively search for mass-balance models with two options
 *      -minimal on or off
 *      -range   on or off
 *      
 */
  int i, j, n;
  int quit, print, first;
  int first_of_model_size, model_size;
  unsigned long minimal_bits, good_bits;
  char token[MAX_LENGTH];

  n = count_unknowns;		/* columns in A, C, E */
  klmd = max_row_count - 2;
  nklmd = n + klmd;
  n2d = n + 2;

  max_good = MAX_MODELS;
  max_bad = MAX_MODELS;
  max_minimal = MAX_MODELS;

  good =
    (unsigned long *) PHRQ_malloc ((size_t) max_good *
				   sizeof (unsigned long));
  if (good == NULL)
    malloc_error ();
  count_good = 0;

  bad =
    (unsigned long *) PHRQ_malloc ((size_t) max_bad * sizeof (unsigned long));
  if (bad == NULL)
    malloc_error ();
  count_bad = 0;

  minimal =
    (unsigned long *) PHRQ_malloc ((size_t) max_minimal *
				   sizeof (unsigned long));
  if (minimal == NULL)
    malloc_error ();
  count_minimal = 0;

  col_back = (int *) PHRQ_malloc ((size_t) max_column_count * sizeof (int));
  if (col_back == NULL)
    malloc_error ();

  row_back = (int *) PHRQ_malloc ((size_t) max_row_count * sizeof (int));
  if (row_back == NULL)
    malloc_error ();

/*
 *   Allocate space for arrays
 */
  cu = (LDBLE *) PHRQ_malloc ((size_t) 2 * nklmd * sizeof (LDBLE));
  if (cu == NULL)
    malloc_error ();
  memset(cu, 0, ((size_t) (2 * nklmd * sizeof(LDBLE))));
  iu = (int *) PHRQ_malloc ((size_t) 2 * nklmd * sizeof (int));
  if (iu == NULL)
    malloc_error ();
  is = (int *) PHRQ_malloc ((size_t) klmd * sizeof (int));
  if (is == NULL)
    malloc_error ();

  for (i = 0; i < 79; i++)
    token[i] = '=';
  token[79] = '\0';
/*
 *   Set solutions, largest bit is final solution, smallest bit is initial solution 1 
 *   Set phases, largest bit is last phase, smallest bit is first phase
 *   Set current bits to complete list.
 */
  soln_bits = 0;
  if (inv_ptr->count_solns + inv_ptr->count_phases > 32)
  {
    error_msg
      ("For inverse modeling, sum of initial solutions and phases must be <= 32.\n\tFor all reasonable calculations, the sum should be much less than 32.",
       STOP);
  }
  for (i = inv_ptr->count_solns; i > 0; i--)
  {
    temp_bits = 1 << (i - 1);
    soln_bits += temp_bits;
  }
  if (check_solns (inv_ptr) == ERROR)
  {
    error_msg ("Calculations terminating.", STOP);
  }
/*
 *   solutions are in highest bits, phases are in lower bits;
 */
/*
 *   All combinations of solutions 
 */
  first = TRUE;
  for (;
       get_bits (soln_bits, inv_ptr->count_solns - 2,
		 inv_ptr->count_solns - 1) > 0; soln_bits--)
  {
/*
 *   Loop through all models of of descending size
 */
    for (model_size = inv_ptr->count_phases; model_size >= 0; model_size--)
    {
      first_of_model_size = TRUE;
      quit = TRUE;
      while (next_set_phases (inv_ptr, first_of_model_size, model_size) ==
	     TRUE)
      {
	first_of_model_size = FALSE;
	current_bits = (soln_bits << inv_ptr->count_phases) + phase_bits;

	if (subset_bad (current_bits) == TRUE
	    || subset_minimal (current_bits) == TRUE)
	  continue;
	quit = FALSE;
/*  
 *   Switch for finding minimal models only
 */
	if (inv_ptr->minimal == TRUE
	    && superset_minimal (current_bits) == TRUE)
	  continue;
/*
 *   Solve for minimum epsilons, continue if no solution found.  
 */
	if (solve_with_mask (inv_ptr, current_bits) == ERROR)
	{
	  save_bad (current_bits);
	  if (first == TRUE)
	  {
	    post_mortem ();
	    quit = TRUE;
	    break;
	  }
	  else
	  {
	    continue;
	  }
	}
	first = FALSE;
/*
 *   Model has been found, set bits 
 */
	good_bits = current_bits;
	for (i = 0; i < inv_ptr->count_phases; i++)
	{
	  if (equal (delta1[i + inv_ptr->count_solns], 0.0, TOL) == TRUE)
	  {
	    good_bits = set_bit (good_bits, i, 0);
	  }
	}
	for (i = 0; i < inv_ptr->count_solns; i++)
	{
	  if (equal (delta1[i], 0.0, TOL) == TRUE)
	  {
	    good_bits = set_bit (good_bits, i + inv_ptr->count_phases, 0);
	  }
	}
/*
 *   Determine if model is new
 */
	for (j = 0; j < count_good; j++)
	{
	  if (good_bits == good[j])
	    break;
	}
/*
 *  Calculate ranges and print model only if NOT looking for minimal models
 */
	print = FALSE;
	if (j >= count_good && inv_ptr->minimal == FALSE)
	{
	  print = TRUE;
	  save_good (good_bits);
	  if (inv_ptr->range == TRUE)
	  {
	    range (inv_ptr, good_bits);
	  }
	  print_model (inv_ptr);
	  punch_model (inv_ptr);
	  dump_netpath_pat (inv_ptr);
	}
/*
 *   If superset of a minimal model continue
 */
	minimal_bits = good_bits;
	if (superset_minimal (minimal_bits) == TRUE)
	{
	  if (print == TRUE)
	  {
	    if (pr.inverse == TRUE && pr.all == TRUE)
	    {
	      output_msg (OUTPUT_MESSAGE, "%s\n\n", token);
	    }
	  }
	  continue;
	}
/*
 *   If not superset of minimal model, find minimal model
 */
	minimal_bits = minimal_solve (inv_ptr, minimal_bits);
	if (minimal_bits == good_bits && print == TRUE)
	{
	  if (pr.inverse == TRUE && pr.all == TRUE)
	  {
	    output_msg (OUTPUT_MESSAGE,
			"\nModel contains minimum number of phases.\n");
	  }
	}
	if (print == TRUE)
	{
	  if (pr.inverse == TRUE && pr.all == TRUE)
	  {
	    output_msg (OUTPUT_MESSAGE, "%s\n\n", token);
	  }
	}
	for (j = 0; j < count_good; j++)
	{
	  if (minimal_bits == good[j])
	    break;
	}
	if (j >= count_good)
	{
	  save_good (minimal_bits);
	  if (inv_ptr->range == TRUE)
	  {
	    range (inv_ptr, minimal_bits);
	  }
	  print_model (inv_ptr);
	  if (pr.inverse == TRUE && pr.all == TRUE)
	  {
	    output_msg (OUTPUT_MESSAGE,
			"\nModel contains minimum number of phases.\n");
	    output_msg (OUTPUT_MESSAGE, "%s\n\n", token);
	  }
	  punch_model (inv_ptr);
	  dump_netpath_pat (inv_ptr);
	}
	save_minimal (minimal_bits);
      }
      if (quit == TRUE)
	break;
    }
  }
/*
 *   Summary print
 */
  if (pr.inverse == TRUE && pr.all == TRUE)
  {
    output_msg (OUTPUT_MESSAGE, "\nSummary of inverse modeling:\n\n");
    output_msg (OUTPUT_MESSAGE, "\tNumber of models found: %d\n", count_good);
    output_msg (OUTPUT_MESSAGE, "\tNumber of minimal models found: %d\n",
		count_minimal);
    output_msg (OUTPUT_MESSAGE,
		"\tNumber of infeasible sets of phases saved: %d\n",
		count_bad);
    output_msg (OUTPUT_MESSAGE, "\tNumber of calls to cl1: %d\n",
		count_calls);
  }
  array = (LDBLE *) free_check_null (array);
  delta = (LDBLE *) free_check_null (delta);
  array1 = (LDBLE *) free_check_null (array1);
  zero = (LDBLE *) free_check_null (zero);
  res = (LDBLE *) free_check_null (res);
  delta1 = (LDBLE *) free_check_null (delta1);
  delta2 = (LDBLE *) free_check_null (delta2);
  delta3 = (LDBLE *) free_check_null (delta3);
  delta_save = (LDBLE *) free_check_null (delta_save);
  cu = (LDBLE *) free_check_null (cu);
  iu = (int *) free_check_null (iu);
  is = (int *) free_check_null (is);
  col_name = (char **) free_check_null (col_name);
  row_name = (char **) free_check_null (row_name);
  col_back = (int *) free_check_null (col_back);
  row_back = (int *) free_check_null (row_back);
  min_delta = (LDBLE *) free_check_null (min_delta);
  max_delta = (LDBLE *) free_check_null (max_delta);
  good = (unsigned long *) free_check_null (good);
  bad = (unsigned long *) free_check_null (bad);
  minimal = (unsigned long *) free_check_null (minimal);

  return (OK);
}

/* ---------------------------------------------------------------------- */
unsigned long
minimal_solve (struct inverse *inv_ptr, unsigned long minimal_bits)
/* ---------------------------------------------------------------------- */
{
/*
 *   Starting with phases indicated in minimal bits, sequentially
 *   remove phases to find minimal solution
 */
  int i;
  unsigned long temp_bits_l;
  if (debug_inverse == TRUE)
  {
    output_msg (OUTPUT_MESSAGE, "Beginning minimal solve: \n");
    bit_print (minimal_bits, inv_ptr->count_phases + inv_ptr->count_solns);
  }
  for (i = 0; i < inv_ptr->count_phases + inv_ptr->count_solns - 1; i++)
  {
    if (get_bits (minimal_bits, i, 1) == 0)
      continue;
    temp_bits_l = 1 << i;	/* 0's and one 1 */
    temp_bits_l = ~temp_bits_l;	/* 1's and one 0 */
    minimal_bits = minimal_bits & temp_bits_l;
    if (debug_inverse == TRUE)
    {
      output_msg (OUTPUT_MESSAGE, "Solving for minimal\n");
      bit_print (minimal_bits, inv_ptr->count_phases + inv_ptr->count_solns);
    }

/*
 *   minimal_bits can not be superset of a minimal model, but
 *   could be subset of one of the sets of minerals with no feasible solution
 *   If it is a subset, then replace mineral and go on to next 
 */
    if (subset_bad (minimal_bits) == TRUE)
    {
      /* put bit back */
      minimal_bits = minimal_bits | ~temp_bits_l;	/* 0's and one 1 */
      continue;
    }
    if (solve_with_mask (inv_ptr, minimal_bits) == ERROR)
    {
      save_bad (minimal_bits);
      /* put bit back */
      minimal_bits = minimal_bits | ~temp_bits_l;	/* 0's and one 1 */
    }

  }
  if (debug_inverse == TRUE)
  {
    output_msg (OUTPUT_MESSAGE, "\n\nMINIMAL MODEL\n\n");
    bit_print (minimal_bits, inv_ptr->count_phases + inv_ptr->count_solns);
  }
  solve_with_mask (inv_ptr, minimal_bits);
  return (minimal_bits);
}

/* ---------------------------------------------------------------------- */
int
solve_with_mask (struct inverse *inv_ptr, unsigned long cur_bits)
/* ---------------------------------------------------------------------- */
{
/*
 *   Uses cur_bits to zero out columns of the array and then solves.
 */
  int i, k, l, m, n;

/*
 *   Calculate dimensions
 */
  k = row_mb;			/* rows in A */
  l = row_epsilon - row_mb;	/* rows in C */
  m = count_rows - row_epsilon;	/* rows in E */
  n = count_unknowns;



  memcpy ((void *) &(res[0]), (void *) &(zero[0]),
	  (size_t) max_row_count * sizeof (LDBLE));
  memcpy ((void *) &(delta2[0]), (void *) &(delta[0]),
	  (size_t) max_column_count * sizeof (LDBLE));
  memcpy ((void *) &(delta_save[0]), (void *) &(zero[0]),
	  (size_t) max_column_count * sizeof (LDBLE));

  shrink (inv_ptr, array, array1,
	  &k, &l, &m, &n, cur_bits, delta2, col_back, row_back);
  /*
   *  Save delta constraints
   */
  for (i = 0; i < n; i++)
  {
    delta_save[col_back[i]] = delta2[i];
  }


  if (debug_inverse == TRUE)
  {
    output_msg (OUTPUT_MESSAGE, "\nColumns\n");
    for (i = 0; i < n; i++)
    {
      output_msg (OUTPUT_MESSAGE, "\t%d\t%s\n", i, col_name[col_back[i]]);
    }

    output_msg (OUTPUT_MESSAGE, "\nRows\n");
    for (i = 0; i < k + l + m; i++)
    {
      output_msg (OUTPUT_MESSAGE, "\t%d\t%s\n", i, row_name[row_back[i]]);
    }

    output_msg (OUTPUT_MESSAGE, "\nA and B arrays:\n\n");
    array_print (array1, k + l + m, n + 1, max_column_count);

    output_msg (OUTPUT_MESSAGE, "\nInput delta vector:\n");
    for (i = 0; i < n; i++)
    {
      output_msg (OUTPUT_MESSAGE, "%6d  %-12.12s %10.2e", i,
		  col_name[col_back[i]], (double) delta2[i]);
      output_msg (OUTPUT_MESSAGE, "\n");
    }

    for (i = 0; i < k + l + m; i++)
    {
      if (res[i] == 0)
	continue;
      output_msg (OUTPUT_MESSAGE, "\nInput res is non zero:\n");
      output_msg (OUTPUT_MESSAGE, "%6d  %-12.12s %10.2e", i,
		  row_name[row_back[i]], (double) res[i]);
      output_msg (OUTPUT_MESSAGE, "\n");
    }
  }
/*
 *   Call CL1
 */

  if (debug_inverse == TRUE)
  {
    output_msg (OUTPUT_MESSAGE,
		"k, l, m, n, max_col, max_row\t%d\t%d\t%d\t%d\t%d\t%d\n", k,
		l, m, n, max_column_count, max_row_count);
  }

  kode = 1;
  iter = 1000;
  count_calls++;

#ifdef INVERSE_CL1MP
  if (inv_ptr->mp == TRUE)
  {
    cl1mp (k, l, m, n,
	   nklmd, n2d, array1,
	   &kode, inv_ptr->mp_tolerance, &iter,
	   delta2, res, &error, cu, iu, is, TRUE, inv_ptr->mp_censor);
  }
  else
  {
    cl1 (k, l, m, n,
	 nklmd, n2d, array1,
	 &kode, toler, &iter, delta2, res, &error, cu, iu, is, TRUE);
  }
#else
  cl1 (k, l, m, n,
       nklmd, n2d, array1,
       &kode, toler, &iter, delta2, res, &error, cu, iu, is, TRUE);
#endif
  if (kode == 3)
  {
    sprintf (error_string,
	     "Exceeded maximum iterations in inverse modeling: %d.\n"
	     "Recompile program with larger limit.", iter);
    error_msg (error_string, STOP);
  }
  memcpy ((void *) &(delta1[0]), (void *) &(zero[0]),
	  (size_t) max_column_count * sizeof (LDBLE));
  for (i = 0; i < n; i++)
  {
    delta1[col_back[i]] = delta2[i];
  }

/*
 *   Debug, write results 
 */

  if (debug_inverse == TRUE)
  {
    output_msg (OUTPUT_MESSAGE, "kode: %d\titer: %d\terror: %e\n", kode, iter,
		(double) error);
    output_msg (OUTPUT_MESSAGE, "\nsolution vector:\n");
    for (i = 0; i < n; i++)
    {
      output_msg (OUTPUT_MESSAGE, "%6d  %-12.12s %10.2e", i,
		  col_name[col_back[i]], (double) delta2[i]);
      output_msg (OUTPUT_MESSAGE, "\n");
    }

    output_msg (OUTPUT_MESSAGE, "\nresidual vector:\n");
    for (i = 0; i < (k + l + m); i++)
    {
      output_msg (OUTPUT_MESSAGE, "%6d  %-12.12s %10.2e\n", i,
		  row_name[row_back[i]], (double) res[i]);
    }
  }

  if (kode != 0)
  {
    return (ERROR);
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
unsigned long
get_bits (unsigned long bits, int position, int number)
/* ---------------------------------------------------------------------- */
{
/*
 *   Returns number of bits from position and below.
 *   position begins at 0.
 */
  return ((bits >> (position + 1 - number)) & ~(~0 << number));
}

/* ---------------------------------------------------------------------- */
int
save_minimal (unsigned long bits)
/* ---------------------------------------------------------------------- */
{
/*
 *   Keeps list of minimal models
 */
  minimal[count_minimal] = bits;
  count_minimal++;
  if (count_minimal >= max_minimal)
  {
    max_minimal *= 2;
    minimal =
      (unsigned long *) PHRQ_realloc (minimal,
				      (size_t) max_minimal *
				      sizeof (unsigned long));
    if (minimal == NULL)
      malloc_error ();
  }
  return (TRUE);
}

/* ---------------------------------------------------------------------- */
int
save_good (unsigned long bits)
/* ---------------------------------------------------------------------- */
{
/*
 *   Keeps list of good models, not necessarily minimal
 */
  good[count_good] = bits;
  count_good++;
  if (count_good >= max_good)
  {
    max_good *= 2;
    good =
      (unsigned long *) PHRQ_realloc (good,
				      (size_t) max_good *
				      sizeof (unsigned long));
    if (good == NULL)
      malloc_error ();
  }
  return (TRUE);
}

/* ---------------------------------------------------------------------- */
int
save_bad (unsigned long bits)
/* ---------------------------------------------------------------------- */
{
/*
 *   Keeps list of sets of phases with no feasible solution
 */
  bad[count_bad] = bits;
  count_bad++;
  if (count_bad >= max_bad)
  {
    max_bad *= 2;
    bad =
      (unsigned long *) PHRQ_realloc (bad,
				      (size_t) max_bad *
				      sizeof (unsigned long));
    if (bad == NULL)
      malloc_error ();
  }
  return (TRUE);
}

/* ---------------------------------------------------------------------- */
int
superset_minimal (unsigned long bits)
/* ---------------------------------------------------------------------- */
{
/*
 *   Checks whether bits is a superset of any of the minimal models
 */
  int i;
  unsigned long temp_bits_l;
  for (i = 0; i < count_minimal; i++)
  {
    temp_bits_l = bits | minimal[i];
    if (temp_bits_l == bits)
    {
      return (TRUE);
    }
  }
  return (FALSE);
}

/* ---------------------------------------------------------------------- */
int
subset_bad (unsigned long bits)
/* ---------------------------------------------------------------------- */
{
/*
 *   Checks whether bits is a superset of any of the bad models
 */
  int i;
  unsigned long temp_bits_l;
  for (i = 0; i < count_bad; i++)
  {
    temp_bits_l = bits | bad[i];
    if (temp_bits_l == bad[i])
    {
      return (TRUE);
    }
  }
  return (FALSE);
}

/* ---------------------------------------------------------------------- */
int
subset_minimal (unsigned long bits)
/* ---------------------------------------------------------------------- */
{
/*
 *   Checks whether bits is a subset of any of the minimal models
 */
  int i;
  unsigned long temp_bits_l;
  for (i = 0; i < count_minimal; i++)
  {
    temp_bits_l = bits | minimal[i];
    if (temp_bits_l == minimal[i])
    {
      return (TRUE);
    }
  }
  return (FALSE);
}

/* ---------------------------------------------------------------------- */
int
bit_print (unsigned long bits, int l)
/* ---------------------------------------------------------------------- */
{
/*
 *   Prints l bits of an unsigned long
 */
  int i;

  for (i = l - 1; i >= 0; i--)
  {
    output_msg (OUTPUT_MESSAGE, "%lu  ", get_bits (bits, i, 1));
  }
  output_msg (OUTPUT_MESSAGE, "\n");
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
print_model (struct inverse *inv_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Prints model
 */
  int i, j, k;
  int column;
  int print_msg;
  struct solution *solution_ptr;
  struct master *master_ptr;
  struct isotope *isotope_ptr;
  LDBLE d1, d2, d3, d4;
  char token[MAX_LENGTH];
/*
 *   Update screen
 */
  status (count_good, NULL);
/*
 *   print solution data, epsilons, and revised data
 */
  if (pr.inverse == FALSE || pr.all == FALSE)
    return (OK);
  max_pct = 0;
  scaled_error = 0;
  for (i = 0; i < inv_ptr->count_solns; i++)
  {
    if (equal (delta1[i], 0.0, toler) == TRUE)
      continue;
    solution_ptr = solution_bsearch (inv_ptr->solns[i], &j, TRUE);
    xsolution_zero ();
    for (j = 0; solution_ptr->totals[j].description != NULL; j++)
    {
      master_ptr = master_bsearch (solution_ptr->totals[j].description);
      master_ptr->total = solution_ptr->totals[j].moles;
    }

    output_msg (OUTPUT_MESSAGE, "\nSolution %d: %s\n", inv_ptr->solns[i],
		solution_ptr->description);
    output_msg (OUTPUT_MESSAGE, "\n%15.15s   %12.12s   %12.12s   %12.12s\n",
		"  ", "Input", "Delta", "Input+Delta");
    master_alk->total = solution_ptr->total_alkalinity;
    if (inv_ptr->carbon == TRUE)
    {
      d1 = solution_ptr->ph;
      d2 = delta1[col_ph + i] / delta1[i];
      d3 = d1 + d2;
      if (equal (d1, 0.0, MIN_TOTAL_INVERSE) == TRUE)
	d1 = 0.0;
      if (equal (d2, 0.0, MIN_TOTAL_INVERSE) == TRUE)
	d2 = 0.0;
      if (equal (d3, 0.0, MIN_TOTAL_INVERSE) == TRUE)
	d3 = 0.0;
      output_msg (OUTPUT_MESSAGE, "%15.15s   %12.3e  +%12.3e  =%12.3e\n",
		  "pH", (double) d1, (double) d2, (double) d3);
      if (inv_ptr->ph_uncertainties[i] > 0)
      {
	scaled_error += fabs (d2) / inv_ptr->ph_uncertainties[i];
/* debug 
				output_msg(OUTPUT_MESSAGE, "%e\t%e\t%e\n", fabs(d2) / inv_ptr->ph_uncertainties[i], fabs(d2), inv_ptr->ph_uncertainties[i]);
 */
      }
      else if (d2 != 0.0)
      {
	error_msg ("Computing delta pH/uncertainty", CONTINUE);
      }
    }
    for (j = 0; j < inv_ptr->count_elts; j++)
    {
      if (inv_ptr->elts[j].master->s == s_eminus)
	continue;
      d1 = inv_ptr->elts[j].master->total;
      d2 = delta1[col_epsilon + j * inv_ptr->count_solns + i] / delta1[i];
      d3 = d1 + d2;

      if (equal (d1, 0.0, MIN_TOTAL_INVERSE) == TRUE)
	d1 = 0.0;
      if (equal (d2, 0.0, MIN_TOTAL_INVERSE) == TRUE)
	d2 = 0.0;
      if (equal (d3, 0.0, MIN_TOTAL_INVERSE) == TRUE)
	d3 = 0.0;

      output_msg (OUTPUT_MESSAGE, "%15.15s   %12.3e  +%12.3e  =%12.3e\n",
		  inv_ptr->elts[j].master->elt->name, (double) d1,
		  (double) d2, (double) d3);
      if (equal (d1, 0.0, MIN_TOTAL_INVERSE) == FALSE)
      {
	d3 = fabs (d2 / d1);
	if (d3 > max_pct)
	  max_pct = d3;
      }
      d4 = 0;
      if (inv_ptr->elts[j].uncertainties[i] > 0)
      {
	d4 = fabs (inv_ptr->elts[j].uncertainties[i] * d1);
      }
      else if (inv_ptr->elts[j].uncertainties[i] < 0)
      {
	d4 = -inv_ptr->elts[j].uncertainties[i];
      }
      if (d4 > 0)
      {
	scaled_error += fabs (d2) / d4;
/* debug 
				output_msg(OUTPUT_MESSAGE, "%e\t%e\t%e\n", fabs(d2) / d4, fabs(d2), d4);
 */
      }
      else if (d2 != 0.0)
      {
	error_msg ("Computing delta element/uncertainty", CONTINUE);
      }
    }
    if (inv_ptr->count_isotopes > 0)
    {
      /* adjustments to solution isotope composition */
      for (j = 0; j < inv_ptr->count_isotope_unknowns; j++)
      {
	for (k = 0; k < solution_ptr->count_isotopes; k++)
	{
	  if (inv_ptr->isotope_unknowns[j].elt_name !=
	      solution_ptr->isotopes[k].elt_name ||
	      inv_ptr->isotope_unknowns[j].isotope_number !=
	      solution_ptr->isotopes[k].isotope_number)
	    continue;
	  d1 = solution_ptr->isotopes[k].ratio;
	  d2 =
	    delta1[col_isotopes + i * inv_ptr->count_isotope_unknowns +
		   j] / delta1[i];
	  d3 = d1 + d2;

	  if (equal (d1, 0.0, MIN_TOTAL_INVERSE) == TRUE)
	    d1 = 0.0;
	  if (equal (d2, 0.0, MIN_TOTAL_INVERSE) == TRUE)
	    d2 = 0.0;
	  if (equal (d3, 0.0, MIN_TOTAL_INVERSE) == TRUE)
	    d3 = 0.0;
	  sprintf (token, "%d%s",
		   (int) inv_ptr->isotope_unknowns[j].isotope_number,
		   inv_ptr->isotope_unknowns[j].elt_name);
	  output_msg (OUTPUT_MESSAGE, "%15.15s   %12g  +%12g  =%12g\n", token,
		      (double) d1, (double) d2, (double) d3);
/*
					if (equal(d1, 0.0, MIN_TOTAL_INVERSE) == FALSE ) {
						d3 = fabs(d2/d1);
						if (d3 > max_pct) max_pct = d3;
					}
 */
	  if (solution_ptr->isotopes[k].x_ratio_uncertainty > 0)
	  {
	    scaled_error +=
	      fabs (d2) / solution_ptr->isotopes[k].x_ratio_uncertainty;
/* debug
						output_msg(OUTPUT_MESSAGE, "%e\t%e\t%e\n", fabs(d2) / solution_ptr->isotopes[k].x_ratio_uncertainty , fabs(d2), solution_ptr->isotopes[k].x_ratio_uncertainty);
 */
	  }
	  else if (d2 != 0.0)
	  {
	    error_msg ("Computing delta solution isotope/uncertainty",
		       CONTINUE);
	  }
	}
      }
    }
  }

/*
 *    Adjustments to phases
 */
  print_msg = FALSE;
  if (inv_ptr->count_isotopes > 0)
  {
    output_msg (OUTPUT_MESSAGE, "\nIsotopic composition of phases:\n");
    for (i = 0; i < inv_ptr->count_phases; i++)
    {
      if (inv_ptr->phases[i].count_isotopes == 0)
	continue;
      j = col_phases + i;
      if (equal (delta1[j], 0.0, toler) == TRUE &&
	  equal (min_delta[j], 0.0, toler) == TRUE &&
	  equal (max_delta[j], 0.0, toler) == TRUE)
	continue;
      isotope_ptr = inv_ptr->phases[i].isotopes;
      for (j = 0; j < inv_ptr->count_isotopes; j++)
      {
	for (k = 0; k < inv_ptr->phases[i].count_isotopes; k++)
	{
	  if (inv_ptr->isotopes[j].elt_name !=
	      isotope_ptr[k].elt_name ||
	      inv_ptr->isotopes[j].isotope_number !=
	      isotope_ptr[k].isotope_number)
	    continue;
	  d1 = isotope_ptr[k].ratio;
	  column = col_phase_isotopes + i * inv_ptr->count_isotopes + j;
	  if (delta1[col_phases + i] != 0.0)
	  {
	    d2 = delta1[column] / delta1[col_phases + i];
	  }
	  else
	  {
	    continue;
	  }
	  d3 = d1 + d2;
	  if (equal (d1, 0.0, 1e-7) == TRUE)
	    d1 = 0.0;
	  if (equal (d2, 0.0, 1e-7) == TRUE)
	    d2 = 0.0;
	  if (equal (d3, 0.0, 1e-7) == TRUE)
	    d3 = 0.0;
	  sprintf (token, "%d%s %s",
		   (int) inv_ptr->isotopes[j].isotope_number,
		   inv_ptr->isotopes[j].elt_name,
		   inv_ptr->phases[i].phase->name);
	  output_msg (OUTPUT_MESSAGE, "%15.15s   %12g  +%12g  =%12g", token,
		      (double) d1, (double) d2, (double) d3);
	  if (fabs (d2) > (isotope_ptr[k].ratio_uncertainty + toler))
	  {
	    output_msg (OUTPUT_MESSAGE, " **");
	    print_msg = TRUE;
	  }
	  output_msg (OUTPUT_MESSAGE, "\n");
	  if (isotope_ptr[k].ratio_uncertainty > 0)
	  {
	    scaled_error += fabs (d2) / isotope_ptr[k].ratio_uncertainty;
/* debug
						output_msg(OUTPUT_MESSAGE, "%e\t%e\t%e\n", fabs(d2) / isotope_ptr[k].ratio_uncertainty, fabs(d2), isotope_ptr[k].ratio_uncertainty);
 */
	  }
	  else if (d2 != 0.0)
	  {
	    error_msg ("Computing delta phase isotope/uncertainty", CONTINUE);
	  }
	}
      }
    }
  }
  if (print_msg == TRUE)
  {
    output_msg (OUTPUT_MESSAGE,
		"\n**\tWARNING: The adjustment to at least one isotopic"
		"\n\tcomposition of a phase exceeded the specified uncertainty"
		"\n\tfor the phase.  If the phase is not constrained to dissolve"
		"\n\tor precipitate, then the isotopic composition of the phase"
		"\n\tis also unconstrained.\n");
  }
  output_msg (OUTPUT_MESSAGE, "\n%-20.20s   %7s   %12.12s   %12.12s\n",
	      "Solution fractions:", " ", "Minimum", "Maximum");
  for (i = 0; i < inv_ptr->count_solns; i++)
  {
    d1 = delta1[i];
    d2 = min_delta[i];
    d3 = max_delta[i];
    if (equal (d1, 0.0, MIN_TOTAL_INVERSE) == TRUE)
      d1 = 0.0;
    if (equal (d2, 0.0, MIN_TOTAL_INVERSE) == TRUE)
      d2 = 0.0;
    if (equal (d3, 0.0, MIN_TOTAL_INVERSE) == TRUE)
      d3 = 0.0;
    output_msg (OUTPUT_MESSAGE, "%11s%4d   %12.3e   %12.3e   %12.3e\n",
		"Solution", inv_ptr->solns[i], (double) d1, (double) d2,
		(double) d3);
  }

  output_msg (OUTPUT_MESSAGE, "\n%-25.25s   %2s   %12.12s   %12.12s\n",
	      "Phase mole transfers:", " ", "Minimum", "Maximum");
  for (i = col_phases; i < col_redox; i++)
  {
    if (equal (delta1[i], 0.0, toler) == TRUE &&
	equal (min_delta[i], 0.0, toler) == TRUE &&
	equal (max_delta[i], 0.0, toler) == TRUE)
      continue;
    d1 = delta1[i];
    d2 = min_delta[i];
    d3 = max_delta[i];
    if (equal (d1, 0.0, MIN_TOTAL_INVERSE) == TRUE)
      d1 = 0.0;
    if (equal (d2, 0.0, MIN_TOTAL_INVERSE) == TRUE)
      d2 = 0.0;
    if (equal (d3, 0.0, MIN_TOTAL_INVERSE) == TRUE)
      d3 = 0.0;
    output_msg (OUTPUT_MESSAGE, "%15.15s   %12.3e   %12.3e   %12.3e   %s\n",
		col_name[i], (double) d1, (double) d2, (double) d3,
		inv_ptr->phases[i - col_phases].phase->formula);
  }

  output_msg (OUTPUT_MESSAGE, "\n%-25.25s\n", "Redox mole transfers:");
  for (i = col_redox; i < col_epsilon; i++)
  {
    if (equal (delta1[i], 0.0, toler) == TRUE)
      continue;
    output_msg (OUTPUT_MESSAGE, "%15.15s   %12.3e\n", col_name[i],
		(double) delta1[i]);
  }

  output_msg (OUTPUT_MESSAGE,
	      "\nSum of residuals (epsilons in documentation):      %12.3e\n",
	      ((double) error / SCALE_EPSILON));
  output_msg (OUTPUT_MESSAGE,
	      "Sum of delta/uncertainty limit:                    %12.3e\n",
	      (double) scaled_error);
  output_msg (OUTPUT_MESSAGE,
	      "Maximum fractional error in element concentration: %12.3e\n",
	      (double) max_pct);
/*
 *   Flush buffer after each model
 */
  output_fflush (OUTPUT_MESSAGE);
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
punch_model_heading (struct inverse *inv_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Prints model headings to selected output file
 */
  int i;
  char token[MAX_LENGTH];
  if (punch.in == FALSE || pr.punch == FALSE || punch.inverse == FALSE)
    return (OK);
/*
 *  Print sum of residuals and maximum fractional error
 */
  if (punch.high_precision == FALSE)
  {
    output_msg (OUTPUT_PUNCH, "%12s\t%12s\t%12s\t", "Sum_resid",
		"Sum_Delta/U", "MaxFracErr");
  }
  else
  {
    output_msg (OUTPUT_PUNCH, "%20s\t%20s\t%20s\t", "Sum_resid",
		"Sum_Delta/U", "MaxFracErr");
  }
/*
 *   Print solution numbers
 */
  for (i = 0; i < inv_ptr->count_solns; i++)
  {
    sprintf (token, "Soln %d", inv_ptr->solns[i]);
    if (punch.high_precision == FALSE)
    {
      output_msg (OUTPUT_PUNCH, "%12s\t%12s\t%12s\t", token, "min", "max");
    }
    else
    {
      output_msg (OUTPUT_PUNCH, "%20s\t%20s\t%20s\t", token, "min", "max");
    }
  }
/*
 *   Print phase names
 */
  for (i = col_phases; i < col_redox; i++)
  {
    if (punch.high_precision == FALSE)
    {
      output_msg (OUTPUT_PUNCH, "%12s\t%12s\t%12s\t", col_name[i], "     min",
		  "     max");
    }
    else
    {
      output_msg (OUTPUT_PUNCH, "%20s\t%20s\t%20s\t", col_name[i], "     min",
		  "     max");
    }
  }
  output_msg (OUTPUT_PUNCH, "\n");
/*
 *   Flush buffer after each model
 */
  output_fflush (OUTPUT_PUNCH);
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
punch_model (struct inverse *inv_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Prints model to selected output file
 */
  int i;
  LDBLE d1, d2, d3;
  if (punch.in == FALSE || pr.punch == FALSE || punch.inverse == FALSE)
    return (OK);
/*
 *   write residual info
 */
  if (punch.high_precision == FALSE)
  {
    output_msg (OUTPUT_PUNCH, "%12.4e\t", ((double) error / SCALE_EPSILON));
    output_msg (OUTPUT_PUNCH, "%12.4e\t", (double) scaled_error);
    output_msg (OUTPUT_PUNCH, "%12.4e\t", (double) max_pct);
  }
  else
  {
    output_msg (OUTPUT_PUNCH, "%20.12e\t", ((double) error / SCALE_EPSILON));
    output_msg (OUTPUT_PUNCH, "%20.12e\t", (double) scaled_error);
    output_msg (OUTPUT_PUNCH, "%20.12e\t", (double) max_pct);
  }
/*
 *   write solution fractions
 */
  for (i = 0; i < inv_ptr->count_solns; i++)
  {
    d1 = delta1[i];
    d2 = min_delta[i];
    d3 = max_delta[i];
    if (equal (d1, 0.0, MIN_TOTAL_INVERSE) == TRUE)
      d1 = 0.0;
    if (equal (d2, 0.0, MIN_TOTAL_INVERSE) == TRUE)
      d2 = 0.0;
    if (equal (d3, 0.0, MIN_TOTAL_INVERSE) == TRUE)
      d3 = 0.0;
    if (punch.high_precision == FALSE)
    {
      output_msg (OUTPUT_PUNCH, "%12.4e\t%12.4e\t%12.4e\t", (double) d1,
		  (double) d2, (double) d3);
    }
    else
    {
      output_msg (OUTPUT_PUNCH, "%20.12e\t%20.12e\t%20.12e\t", (double) d1,
		  (double) d2, (double) d3);
    }
  }
/*
 *   write phase transfers
 */
  for (i = col_phases; i < col_redox; i++)
  {
    d1 = delta1[i];
    d2 = min_delta[i];
    d3 = max_delta[i];
    if (equal (d1, 0.0, MIN_TOTAL_INVERSE) == TRUE)
      d1 = 0.0;
    if (equal (d2, 0.0, MIN_TOTAL_INVERSE) == TRUE)
      d2 = 0.0;
    if (equal (d3, 0.0, MIN_TOTAL_INVERSE) == TRUE)
      d3 = 0.0;
    if (punch.high_precision == FALSE)
    {
      output_msg (OUTPUT_PUNCH, "%12.4e\t%12.4e\t%12.4e\t", (double) d1,
		  (double) d2, (double) d3);
    }
    else
    {
      output_msg (OUTPUT_PUNCH, "%20.12e\t%20.12e\t%20.12e\t", (double) d1,
		  (double) d2, (double) d3);
    }
  }
  output_msg (OUTPUT_PUNCH, "\n");
/*
 *   Flush buffer after each model
 */
  output_fflush (OUTPUT_PUNCH);
  return (OK);
}

/* ---------------------------------------------------------------------- */
unsigned long
set_bit (unsigned long bits, int position, int value)
/* ---------------------------------------------------------------------- */
{
/*
 *   Sets a single bit
 */
  unsigned long temp_bits_l;

  temp_bits_l = 1 << position;
  if (value == 0)
  {
    temp_bits_l = ~temp_bits_l;
    temp_bits_l = bits & temp_bits_l;
  }
  else
  {
    temp_bits_l = bits | temp_bits_l;
  }
  return (temp_bits_l);
}

/* ---------------------------------------------------------------------- */
int
next_set_phases (struct inverse *inv_ptr,
		 int first_of_model_size, int model_size)
/* ---------------------------------------------------------------------- */
{
  int i, j, k;
  unsigned long temp_bits_l;
  static int min_position[32], max_position[32], now[32];

/*
 *   min_ and max_position are arrays, logically with length
 *      of model_size, that contain minimum and maximum
 *      phase numbers that can be in that model position.
 *
 *   now contains a list of phase numbers to mark as in for this
 *      model
 */

/*
 *   Initialize for a given model_size
 */
  if (first_of_model_size == TRUE)
  {
    for (i = 0; i < model_size; i++)
    {
      min_position[i] = i;
      now[i] = i;
      max_position[i] = inv_ptr->count_phases - model_size + i;
    }
  }
  else
  {
/*
 *   Determine next combination of phases for fixed model_size
 */
    for (i = (model_size - 1); i >= 0; i--)
    {
      if (now[i] < max_position[i])
      {
	now[i]++;
	if (i < (model_size - 1))
	{
	  k = now[i];
	  for (j = (i + 1); j < model_size; j++)
	  {
	    k++;
	    now[j] = k;
	  }
	}
	break;
      }
    }
    if (i < 0)
      return (FALSE);
  }
/*
 *   Set bits which switch in phases
 */
  temp_bits_l = 0;
  for (j = 0; j < model_size; j++)
  {
    temp_bits_l += (1 << now[j]);
  }
  phase_bits = temp_bits_l;
  return (TRUE);
}

/* ---------------------------------------------------------------------- */
int
range (struct inverse *inv_ptr, unsigned long cur_bits)
/* ---------------------------------------------------------------------- */
{
/*
 *   Takes the model from cur_bits and sequentially determines the 
 *   minimum and maximum values for each solution fraction and
 *   each phase mass transfer.
 */
  int i, j;
  int k, l, m, n;
  int f;
  unsigned long bits;
  LDBLE error2;
/*
 *   Include forced solutions and phases in range calculation
 */
  for (i = 0; i < inv_ptr->count_solns + inv_ptr->count_phases; i++)
  {
    if (i < inv_ptr->count_phases)
    {
      if (inv_ptr->phases[i].force == TRUE)
      {
	cur_bits = set_bit (cur_bits, i, 1);
      }
    }
    else
    {
      if (inv_ptr->force_solns[i - inv_ptr->count_phases] == TRUE)
      {
	cur_bits = set_bit (cur_bits, i, 1);
      }
    }
  }

  memcpy ((void *) &(min_delta[0]), (void *) &(zero[0]),
	  (size_t) max_column_count * sizeof (LDBLE));
  memcpy ((void *) &(max_delta[0]), (void *) &(zero[0]),
	  (size_t) max_column_count * sizeof (LDBLE));
/*
 *   Switch bits so that phases are high and solutions are low
 */
  bits =
    get_bits (cur_bits, inv_ptr->count_phases + inv_ptr->count_solns - 1,
	      inv_ptr->count_solns);
  bits +=
    (get_bits (cur_bits, inv_ptr->count_phases - 1, inv_ptr->count_phases) <<
     inv_ptr->count_solns);
/*
 *   Do range calculation
 */
  for (i = 0; i < inv_ptr->count_solns + inv_ptr->count_phases; i++)
  {
    if (inv_ptr->count_solns == i + 1)
    {
      min_delta[i] = 1.0;
      max_delta[i] = 1.0;
      continue;
    }
    if (get_bits (bits, i, 1) == 0)
      continue;
/*
 *   Calculate min and max
 */
    for (f = -1; f < 2; f += 2)
    {
      k = row_mb;		/* rows in A */
      l = row_epsilon - row_mb;	/* rows in C */
      m = count_rows - row_epsilon;	/* rows in E */
      n = count_unknowns;	/* number of variables */
/*
 *   Copy equations
 */
      memcpy ((void *) &(array1[0]), (void *) &(array[0]),
	      (size_t) max_column_count * max_row_count * sizeof (LDBLE));
      memcpy ((void *) &(delta2[0]), (void *) &(delta[0]),
	      (size_t) max_column_count * sizeof (LDBLE));
      memcpy ((void *) &(delta3[0]), (void *) &(zero[0]),
	      (size_t) max_column_count * sizeof (LDBLE));
      memcpy ((void *) &(delta_save[0]), (void *) &(zero[0]),
	      (size_t) max_column_count * sizeof (LDBLE));
      memcpy ((void *) &(res[0]), (void *) &(zero[0]),
	      (size_t) max_row_count * sizeof (LDBLE));

/*
 *   Change optimization
 */
      for (j = 0; j < k; j++)
      {
	memcpy ((void *) &(array1[j * max_column_count]), (void *) &(zero[0]),
		(size_t) max_column_count * sizeof (LDBLE));
      }
      array1[i] = 1.0;
      if (f < 1)
      {
	array1[n] = -fabs (inv_ptr->range_max);
      }
      else
      {
	array1[n] = fabs (inv_ptr->range_max);
      }
      shrink (inv_ptr, array1, array1,
	      &k, &l, &m, &n, cur_bits, delta2, col_back, row_back);
      /*
       *  Save delta constraints
       */
      for (j = 0; j < n; j++)
      {
	delta_save[col_back[j]] = delta2[j];
      }
      if (debug_inverse == TRUE)
      {
	output_msg (OUTPUT_MESSAGE, "\nInput delta:\n\n");
	for (j = 0; j < n; j++)
	{
	  output_msg (OUTPUT_MESSAGE, "\t%d %s\t%g\n", j,
		      col_name[col_back[j]], (double) delta2[j]);
	}
	output_msg (OUTPUT_MESSAGE, "\nA and B arrays:\n\n");
	array_print (array1, k + l + m, n + 1, max_column_count);
      }
      kode = 1;
      iter = 200;
      count_calls++;
#ifdef INVERSE_CL1MP
      if (inv_ptr->mp == TRUE)
      {
	cl1mp (k, l, m, n,
	       nklmd, n2d, array1,
	       &kode, inv_ptr->mp_tolerance, &iter,
	       delta2, res, &error2, cu, iu, is, TRUE, inv_ptr->mp_censor);
      }
      else
      {
	cl1 (k, l, m, n,
	     nklmd, n2d, array1,
	     &kode, toler, &iter, delta2, res, &error2, cu, iu, is, TRUE);
      }
#else
      cl1 (k, l, m, n,
	   nklmd, n2d, array1,
	   &kode, toler, &iter, delta2, res, &error2, cu, iu, is, TRUE);
#endif
      if (kode != 0)
      {
	output_msg (OUTPUT_MESSAGE, "Error in subroutine range. Kode = %d\n",
		    kode);
      }

      if (debug_inverse == TRUE)
      {
	output_msg (OUTPUT_MESSAGE, "kode: %d\titer: %d\terror: %e\n", kode,
		    iter, (double) error2);
	output_msg (OUTPUT_MESSAGE, "k, l, m, n: %d\t%d\t%d\t%d\n", k, l, m,
		    n);
	output_msg (OUTPUT_MESSAGE, "\nsolution vector %s\n", col_name[i]);
	for (j = 0; j < n; j++)
	{
	  output_msg (OUTPUT_MESSAGE, "%6d  %-12.12s %10.2e", j,
		      col_name[col_back[j]], (double) delta2[j]);
	  output_msg (OUTPUT_MESSAGE, "\n");
	}

	output_msg (OUTPUT_MESSAGE, "\nresidual vector:\n");
	for (j = 0; j < (k + l + m); j++)
	{
	  output_msg (OUTPUT_MESSAGE, "%6d  %-12.12s %10.2e\n", j,
		      row_name[row_back[j]], (double) res[j]);
	}
      }
      for (j = 0; j < n; j++)
      {
	if (col_back[j] == i)
	  break;
      }
      if (f < 0)
      {
	min_delta[i] = delta2[j];
      }
      else
      {
	max_delta[i] = delta2[j];
      }
      for (j = 0; j < n; j++)
      {
	delta3[col_back[j]] = delta2[j];
      }
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
shrink (struct inverse *inv_ptr, LDBLE * array_in, LDBLE * array_out,
	int *k, int *l, int *m, int *n,
	unsigned long cur_bits,
	LDBLE * delta_l, int *col_back_l, int *row_back_l)
/* ---------------------------------------------------------------------- */
{
/*
 *   Shrink eliminates any rows that are all zeros and any columns
 *   that are not in cur_bits, result is put in array_out
 *
 *   k, l, m, n return the new sizes of the array.
 *   delta is remapped to retain any non-negativity constraints
 *   Col_back maps columns that remain back to original columns
 *   Row_back maps rows that remain back to original rows
 */
  int i, j, row;
  int k1, l1, m1;
  int cur_col, column;
  int nonzero;
/*
 *   Copy array_in to array_out
 */
  if (array_in != array_out)
  {
    for (i = 0; i < (*k + *l + *m); i++)
    {
      memcpy (&(array_out[i * max_column_count]),
	      &(array_in[i * max_column_count]),
	      (size_t) max_column_count * sizeof (LDBLE));
    }
  }
/*
 *   Determine columns to eliminate
 */
  for (i = 0; i < (*n + 1); i++)
    col_back_l[i] = i;

/*   
 *   Drop phases not in model
 */
  for (i = 0; i < inv_ptr->count_phases; i++)
  {
    if (get_bits (cur_bits, i, 1) == 0)
    {
      col_back_l[col_phases + i] = -1;
      /* drop isotopes */
      if (inv_ptr->count_isotopes > 0)
      {
	for (j = 0; j < inv_ptr->count_isotopes; j++)
	{
	  column = col_phase_isotopes + i * inv_ptr->count_isotopes + j;
	  col_back_l[column] = -1;
	}
      }
    }
  }
/*   
 *   Drop solutions not in model
 */
  for (i = 0; i < (inv_ptr->count_solns - 1); i++)
  {
    if (get_bits (cur_bits, inv_ptr->count_phases + i, 1) == 0)
    {
      col_back_l[i] = -1;
      /* drop all epsilons for the solution */
      for (j = 0; j < inv_ptr->count_elts; j++)
      {
	column = col_epsilon + j * inv_ptr->count_solns + i;
	col_back_l[column] = -1;
      }
      /* drop pH for the solution */
      if (inv_ptr->carbon == TRUE)
      {
	column = col_ph + i;
	col_back_l[column] = -1;
      }
      /* drop isotopes */
      if (inv_ptr->count_isotopes > 0)
      {
	for (j = 0; j < inv_ptr->count_isotope_unknowns; j++)
	{
	  column = col_isotopes + i * inv_ptr->count_isotope_unknowns + j;
	  col_back_l[column] = -1;
	}
      }
    }
  }

/*   
 *   Drop epsilons not used
 */
  for (i = col_epsilon; i < *n; i++)
  {
    if (col_back_l[i] < 0)
      continue;
    for (j = 0; j < (*k + *l + *m); j++)
    {
      if (array_out[j * max_column_count + i] != 0)
	break;
    }
    if (j == (*k + *l + *m))
    {
      col_back_l[i] = -1;
    }
  }
/*
 *   rewrite array_out
 */
  cur_col = 0;
  for (i = 0; i < (*n + 1); i++)
  {
    if (col_back_l[i] < 0)
      continue;
    if (cur_col == col_back_l[i])
    {
      cur_col++;
      continue;
    }
    for (j = 0; j < (*k + *l + *m); j++)
    {
      array_out[j * max_column_count + cur_col] =
	array_out[j * max_column_count + i];
    }
    col_back_l[cur_col] = col_back_l[i];
    delta_l[cur_col] = delta_l[i];
    cur_col++;
  }
  *n = cur_col - 1;
/* 
 *   Eliminate unnecessary optimization eqns
 */
  row = 0;
  k1 = 0;
  for (i = 0; i < *k; i++)
  {
    if (memcmp (&(array_out[i * max_column_count]), &(zero[0]),
		(size_t) (*n) * sizeof (LDBLE)) == 0)
    {
      continue;
    }
/*
		memcpy(&(array_out[row * max_column_count]), &(array_out[i * max_column_count]),
		       (size_t) max_column_count * sizeof(LDBLE));
 */
    memcpy (&(array_out[row * max_column_count]),
	    &(array_out[i * max_column_count]),
	    (size_t) (*n + 1) * sizeof (LDBLE));
    row_back_l[row] = i;
    row++;
    k1++;
  }

/* 
 *   Eliminate unnecessary equality eqns
 */
  l1 = 0;
  for (i = *k; i < (*k + *l); i++)
  {
    nonzero = FALSE;
    for (j = 0; j < *n; j++)
    {
      if (equal (array_out[i * max_column_count + j], 0.0, toler) == FALSE)
      {
	nonzero = TRUE;
	break;
      }
    }
    if (nonzero == FALSE)
      continue;
/*
		if (memcmp(&(array_out[i * max_column_count]), &(zero[0]), 
			   (size_t) (*n) * sizeof(LDBLE)) == 0) {
			continue;
		}
 */
/*
		memcpy(&(array_out[row * max_column_count]), &(array_out[i * max_column_count]),
		       (size_t) max_column_count * sizeof(LDBLE));
 */
    memcpy (&(array_out[row * max_column_count]),
	    &(array_out[i * max_column_count]),
	    (size_t) (*n + 1) * sizeof (LDBLE));
    row_back_l[row] = i;
    row++;
    l1++;
  }
/* 
 *   Eliminate unnecessary inequality eqns
 */
  m1 = 0;
  for (i = (*k + *l); i < (*k + *l + *m); i++)
  {
    nonzero = FALSE;
    for (j = 0; j < *n; j++)
    {
      if (equal (array_out[i * max_column_count + j], 0.0, toler) == FALSE)
      {
	nonzero = TRUE;
	break;
      }
    }
    if (nonzero == FALSE)
      continue;
/*
		if (memcmp(&(array_out[i * max_column_count]), &(zero[0]), 
			   (size_t) (*n) * sizeof(LDBLE)) == 0) {
			continue;
		}
 */
/*
		memcpy(&(array_out[row * max_column_count]), &(array_out[i * max_column_count]),
		       (size_t) max_column_count * sizeof(LDBLE));
 */
    memcpy (&(array_out[row * max_column_count]),
	    &(array_out[i * max_column_count]),
	    (size_t) (*n + 1) * sizeof (LDBLE));
    row_back_l[row] = i;
    row++;
    m1++;
  }

  *k = k1;
  *l = l1;
  *m = m1;
/*
 *  Scale all inequality rows
 */

  for (i = *k + *l; i < *k + *l + *m; i++)
  {
    for (j = 0; j < *n + 1; j++)
    {
      array_out[i * max_column_count + j] *= SCALE_ALL;
    }
  }

  return (OK);
}

/* ---------------------------------------------------------------------- */
int
check_solns (struct inverse *inv_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Check_solns checks that each solution can be charge balanced within
 *   the given constraints. If not, it is an error and the program will
 *   terminate.
 */
  int i, j;
  int k, l, m, n;
  int return_value;
  unsigned long bits;
  LDBLE error2;

  memcpy ((void *) &(min_delta[0]), (void *) &(zero[0]),
	  (size_t) max_column_count * sizeof (LDBLE));
  memcpy ((void *) &(max_delta[0]), (void *) &(zero[0]),
	  (size_t) max_column_count * sizeof (LDBLE));

/*
 *   Switch bits so that phases are high and solutions are low
 */
  return_value = OK;
  for (i = 0; i < inv_ptr->count_solns; i++)
  {
    bits = 0;
    bits += 1 << (inv_ptr->count_phases + i);
/*
 *   Check for feasibility of charge balance with given uncertainties
 */
    k = row_mb;			/* rows in A */
    l = row_epsilon - row_mb;	/* rows in C */
    m = count_rows - row_epsilon;	/* rows in E */
    n = count_unknowns;		/* number of variables */
/* debug
	output_msg(OUTPUT_MESSAGE, "\nColumns\n");
	for (j = 0; j < n; j++) {
		output_msg(OUTPUT_MESSAGE, "\t%d\t%s\n", j, col_name[j]);
	}

	output_msg(OUTPUT_MESSAGE, "\nRows\n");
	for (j = 0; j < k + l + m; j++) {
		output_msg(OUTPUT_MESSAGE, "\t%d\t%s\n", j, row_name[j]);
	}

 	output_msg(OUTPUT_MESSAGE, "\nA and B arrays:\n\n");
 	array_print(array, k + l + m,
 		    n + 1, max_column_count);
 */
/*
 *   Copy equations
 */
    memcpy ((void *) &(array1[0]), (void *) &(array[0]),
	    (size_t) max_column_count * max_row_count * sizeof (LDBLE));
    memcpy ((void *) &(delta2[0]), (void *) &(delta[0]),
	    (size_t) max_column_count * sizeof (LDBLE));
    memcpy ((void *) &(res[0]), (void *) &(zero[0]),
	    (size_t) max_row_count * sizeof (LDBLE));

/*
 *   Keep optimization
 */
/*
 *   Zero out mass balance rows and fraction rows
 */
    for (j = row_mb; j < row_charge; j++)
    {
      memcpy ((void *) &(array1[j * max_column_count]), (void *) &(zero[0]),
	      (size_t) max_column_count * sizeof (LDBLE));
    }
/*
 *   Set fraction of solution to 1.0
 */
    array1[(row_charge - 1) * max_column_count + i] = 1.0;
    array1[(row_charge - 1) * max_column_count + n] = 1.0;

/*
 *   Zero out charge balance rows for other solutions
 */
    for (j = 0; j < inv_ptr->count_solns; j++)
    {
      if (j == i)
	continue;
      memcpy ((void *) &(array1[(row_charge + j) * max_column_count]),
	      (void *) &(zero[0]),
	      (size_t) max_column_count * sizeof (LDBLE));
    }

/*
 *   Zero out isotope mole balance
 */
    for (j = row_isotopes; j < row_epsilon; j++)
    {
      memcpy ((void *) &(array1[j * max_column_count]), (void *) &(zero[0]),
	      (size_t) max_column_count * sizeof (LDBLE));
    }

/*
 *   Zero out isotope uncertainties
 */
    for (j = row_isotope_epsilon; j < count_rows; j++)
    {
      memcpy ((void *) &(array1[j * max_column_count]), (void *) &(zero[0]),
	      (size_t) max_column_count * sizeof (LDBLE));
    }
/*
 *   Can't Zero out epsilon constraint rows for other solutions because not sure which
 *   are which
 */

    shrink (inv_ptr, array1, array1,
	    &k, &l, &m, &n, bits, delta2, col_back, row_back);
/* Debug

		output_msg(OUTPUT_MESSAGE, "\nColumns\n");
		for (j = 0; j < n; j++) {
			output_msg(OUTPUT_MESSAGE, "\t%d\t%s\n", j, col_name[col_back[j]]);
		}
		
		output_msg(OUTPUT_MESSAGE, "\nRows\n");
		for (j = 0; j < k + l + m; j++) {
			output_msg(OUTPUT_MESSAGE, "\t%d\t%s\n", j, row_name[row_back[j]]);
		}
		
		output_msg(OUTPUT_MESSAGE, "\nA and B arrays:\n\n");
		array_print(array1, k + l + m,
			    n + 1, max_column_count);

		output_msg(OUTPUT_MESSAGE, "\nInput delta vector:\n");
		for (j=0; j < n; j++) {
			output_msg(OUTPUT_MESSAGE, "%6d  %-12.12s %10.2e", j, col_name[col_back[j]], delta2[j]);
			output_msg(OUTPUT_MESSAGE, "\n");
		}
 */

    kode = 1;
    iter = 200;
    count_calls++;
    cl1 (k, l, m, n,
	 nklmd, n2d, array1,
	 &kode, toler, &iter, delta2, res, &error2, cu, iu, is, TRUE);

    if (kode != 0)
    {
      sprintf (error_string,
	       "Not possible to balance solution %d with input uncertainties.",
	       inv_ptr->solns[i]);
      error_msg (error_string, CONTINUE);
      return_value = ERROR;
    }

/* Debug
		output_msg(OUTPUT_MESSAGE, "kode: %d\titer: %d\terror: %e\n", kode, iter, error);
		output_msg(OUTPUT_MESSAGE, "k, l, m, n: %d\t%d\t%d\t%d\n", k, l, m, n);

		output_msg(OUTPUT_MESSAGE, "\nsolution vector %s\n", col_name[i]);
		for (j = 0; j < n; j++) {
			output_msg(OUTPUT_MESSAGE, "%6d  %-12.12s %10.2e", j, col_name[col_back[j]], delta2[j]);
			output_msg(OUTPUT_MESSAGE, "\n");
		}

		output_msg(OUTPUT_MESSAGE, "\nresidual vector:\n");
		for (j = 0; j < (k + l + m); j++) {
			output_msg(OUTPUT_MESSAGE, "%6d  %-12.12s %10.2e\n", j, row_name[row_back[j]], res[j]);
		}
 */
  }
  return (return_value);
}

/* ---------------------------------------------------------------------- */
int
post_mortem (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Post_mortem simply identifies which equality and inequality of the
 *   array have not been satisfied. 
 *
 */
  int i, j;
  LDBLE sum;
/*
 *   Check equalities
 */
  output_msg (OUTPUT_MESSAGE,
	      "\nPost_mortem examination of inverse modeling:\n\n");
  for (i = row_mb; i < row_epsilon; i++)
  {
    sum = 0;
    for (j = 0; j < count_unknowns; j++)
    {
      sum += delta1[j] * array[i * max_column_count + j];
    }

    if (equal (sum, array[(i * max_column_count) + count_unknowns], toler) ==
	FALSE)
    {
      output_msg (OUTPUT_MESSAGE,
		  "\tERROR: equality not satisfied for %s, %e.\n",
		  row_name[i],
		  sum - array[(i * max_column_count) + count_unknowns]);
    }
  }
/*
 *   Check inequalities
 */
  for (i = row_epsilon; i < count_rows; i++)
  {
    sum = 0;
    for (j = 0; j < count_unknowns; j++)
    {
      sum += delta1[j] * array[i * max_column_count + j];
    }

    if (sum > array[(i * max_column_count) + count_unknowns] + toler)
    {
      output_msg (OUTPUT_MESSAGE,
		  "\tERROR: inequality not satisfied for %s, %e\n",
		  row_name[i],
		  sum - array[(i * max_column_count) + count_unknowns]);
    }
  }
/*
 *   Check dissolution/precipitation constraints
 */
  for (i = 0; i < count_unknowns; i++)
  {
    if (delta_save[i] > 0.5 && delta1[i] < -toler)
    {
      output_msg (OUTPUT_MESSAGE,
		  "\tERROR: Dissolution/precipitation constraint not satisfied for column %d, %s, %e.\n",
		  i, col_name[i], delta1[i]);
    }
    else if (delta_save[i] < -0.5 && delta1[i] > toler)
    {
      output_msg (OUTPUT_MESSAGE,
		  "\tERROR: Dissolution/precipitation constraint not satisfied for column %d, %s, %e.\n",
		  i, col_name[i], delta1[i]);
    }
  }

  return (OK);
}

/* ---------------------------------------------------------------------- */
int
carbon_derivs (struct inverse *inv_ptr)
/* ---------------------------------------------------------------------- */
{
  int i, j, k, n, temp;
  LDBLE c_uncertainty, d_carbon, alk_plus, alk_minus;
  struct solution *solution_ptr_orig, *solution_ptr;

  inv_ptr->dalk_dph = (LDBLE *) free_check_null (inv_ptr->dalk_dph);
  inv_ptr->dalk_dph =
    (LDBLE *) PHRQ_malloc ((size_t) inv_ptr->count_solns * sizeof (LDBLE));
  if (inv_ptr->dalk_dph == NULL)
    malloc_error ();

  inv_ptr->dalk_dc = (LDBLE *) free_check_null (inv_ptr->dalk_dc);
  inv_ptr->dalk_dc =
    (LDBLE *) PHRQ_malloc ((size_t) inv_ptr->count_solns * sizeof (LDBLE));
  if (inv_ptr->dalk_dc == NULL)
    malloc_error ();

  for (i = 0; i < inv_ptr->count_solns; i++)
  {
    solution_ptr_orig = solution_bsearch (inv_ptr->solns[i], &n, TRUE);
    if (solution_ptr_orig == NULL)
    {
      sprintf (error_string, "Solution %d for inverse "
	       "modeling not found.", inv_ptr->solns[i]);
      error_msg (error_string, STOP);
    }
/* 
 *   Find carbon uncertainty
 */
    c_uncertainty = 0;
    d_carbon = 0;
    for (j = 0; j < inv_ptr->count_elts; j++)
    {
      if (inv_ptr->elts[j].master == s_co3->secondary)
      {
	c_uncertainty = inv_ptr->elts[j].uncertainties[i];
	break;
      }
    }
    if (c_uncertainty < 0.0)
    {
      d_carbon = -c_uncertainty;
    }
    else if (c_uncertainty > 0.0)
    {
      for (k = 0; solution_ptr_orig->totals[k].description != NULL; k++)
      {
	if (strcmp (solution_ptr_orig->totals[k].description, "C(4)") == 0)
	{
	  d_carbon = solution_ptr_orig->totals[k].moles /
	    solution_ptr_orig->mass_water * c_uncertainty;
	  break;
	}

      }
    }

/*
 *   Make four copies of solution
 *   Modify ph and carbon in solutions
 */
    set_ph_c (inv_ptr, i, solution_ptr_orig, -5, 0.0, 1.0, 0.0);
    set_ph_c (inv_ptr, i, solution_ptr_orig, -4, 0.0, -1.0, 0.0);
    if (c_uncertainty != 0)
    {
      set_ph_c (inv_ptr, i, solution_ptr_orig, -3, d_carbon, 0.0, 1.0);
      set_ph_c (inv_ptr, i, solution_ptr_orig, -2, d_carbon, 0.0, -1.0);
    }
/* */
    temp = pr.all;
    pr.all = FALSE;
    initial_solutions (FALSE);
    pr.all = temp;
/*
 *   dAlk/dpH
 */
    solution_ptr = solution_bsearch (-5, &n, TRUE);
    alk_plus = solution_ptr->total_alkalinity;
    solution_ptr = solution_bsearch (-4, &n, TRUE);
    alk_minus = solution_ptr->total_alkalinity;
    inv_ptr->dalk_dph[i] = (alk_plus - alk_minus) /
      (2.0 * inv_ptr->ph_uncertainties[i]);
/*
 *   dAlk/dC
 */
    if (d_carbon != 0)
    {
      solution_ptr = solution_bsearch (-3, &n, TRUE);
      alk_plus = solution_ptr->total_alkalinity;
      solution_ptr = solution_bsearch (-2, &n, TRUE);
      alk_minus = solution_ptr->total_alkalinity;
      inv_ptr->dalk_dc[i] = (alk_plus - alk_minus) / (2.0 * d_carbon);
    }
    else
    {
      inv_ptr->dalk_dc[i] = 0.0;
    }
    if (debug_inverse == TRUE)
    {
      output_msg (OUTPUT_MESSAGE, "dAlk/dph = %e\tdAlk/dC = %e\n",
		  (double) inv_ptr->dalk_dph[i],
		  (double) inv_ptr->dalk_dc[i]);
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
set_ph_c (struct inverse *inv_ptr,
	  int i,
	  struct solution *solution_ptr_orig,
	  int n_user_new, LDBLE d_carbon, LDBLE ph_factor, LDBLE c_factor)
/* ---------------------------------------------------------------------- */
{
  int j, n_user_orig;
  struct solution *solution_ptr;
  struct conc *conc_ptr;

  n_user_orig = inv_ptr->solns[i];
  solution_duplicate (n_user_orig, n_user_new);
  solution_ptr = solution_bsearch (n_user_new, &j, TRUE);
  solution_ptr->new_def = TRUE;
  solution_ptr->n_user_end = n_user_new;
  solution_ptr->ph += inv_ptr->ph_uncertainties[i] * ph_factor;
  for (j = 0; solution_ptr->totals[j].description != NULL; j++)
  {
    conc_ptr = &solution_ptr->totals[j];
    conc_ptr->input_conc = conc_ptr->moles / solution_ptr_orig->mass_water;
    conc_ptr->units = string_hsave ("Mol/kgw");
    if (strcmp (conc_ptr->description, "C(4)") == 0)
    {
      conc_ptr->input_conc += d_carbon * c_factor;
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
isotope_balance_equation (struct inverse *inv_ptr, int row, int n)
/* ---------------------------------------------------------------------- */
/*
 *   routine fills in an isotope balance equation
 *
 *   row is the row in array that needs to be filled
 *   n is the isotope number in inv_ptr
 */
{
  int i, j, k;
  LDBLE isotope_number;
  int column;
  LDBLE f;
  struct master *primary_ptr;
  struct solution *solution_ptr;
  struct isotope *isotope_ptr;
/*
 *   Determine primary master species and isotope number for
 *   isotope mass-balance equation
 */
  column = 0;
  primary_ptr = master_bsearch_primary (inv_ptr->isotopes[n].elt_name);
  isotope_number = inv_ptr->isotopes[n].isotope_number;
  /* isotope element must be defined */
  if (primary_ptr == NULL)
  {
    sprintf (error_string, "In isotope calculation: element not defined: %s.",
	     inv_ptr->isotopes[n].elt_name);
    error_msg (error_string, CONTINUE);
    input_error++;
  }

  /* isotope element must be primary */
  if (primary_ptr->primary != TRUE)
  {
    sprintf (error_string, "Isotope mass-balance may only be used"
	     " for total element concentrations.\n"
	     "Secondary species not allowed: %s.",
	     inv_ptr->isotopes[n].elt_name);
    error_msg (error_string, CONTINUE);
    input_error++;
  }

/*
 *   Fill in terms for each solution
 */
  for (i = 0; i < inv_ptr->count_solns; i++)
  {
    if (i == (inv_ptr->count_solns - 1))
    {
      f = -1.0;
    }
    else
    {
      f = 1.0;
    }

    /* mixing fraction term */
    solution_ptr = solution_bsearch (inv_ptr->solns[i], &j, TRUE);
    isotope_ptr = solution_ptr->isotopes;
    for (j = 0; j < solution_ptr->count_isotopes; j++)
    {
      if (isotope_ptr[j].primary == primary_ptr &&
	  isotope_ptr[j].isotope_number == isotope_number)
      {
	array[row * max_column_count + i] +=
	  f * isotope_ptr[j].total * isotope_ptr[j].ratio;
      }
    }

    /* epsilon of total moles of element valence * ratio */
    for (j = 0; j < solution_ptr->count_isotopes; j++)
    {

      /* What to do with H and O, skip for now ??? */
      if (primary_ptr == s_hplus->primary || primary_ptr == s_h2o->primary)
	continue;

      if (isotope_ptr[j].primary == primary_ptr &&
	  isotope_ptr[j].isotope_number == isotope_number)
      {

	/* find column of master for solution i */
	for (k = 0; k < inv_ptr->count_elts; k++)
	{
	  if (isotope_ptr[j].master == inv_ptr->elts[k].master)
	    break;
	}
	column = col_epsilon + (k * inv_ptr->count_solns) + i;
	array[row * max_column_count + column] += f * isotope_ptr[j].ratio;
      }
    }

    /* epsilon of ratio * total of element valence */
    for (j = 0; j < solution_ptr->count_isotopes; j++)
    {

      if (isotope_ptr[j].primary == primary_ptr &&
	  isotope_ptr[j].isotope_number == isotope_number)
      {

	/* find column of epsilon for ratio of valence */
	for (k = 0; k < inv_ptr->count_isotope_unknowns; k++)
	{
	  if (isotope_ptr[j].master == inv_ptr->isotope_unknowns[k].master &&
	      isotope_ptr[j].isotope_number ==
	      inv_ptr->isotope_unknowns[k].isotope_number)
	  {
	    column = col_isotopes + (i * inv_ptr->count_isotope_unknowns) + k;
	  }
	}
	array[row * max_column_count + column] += f * isotope_ptr[j].total;
      }
    }
  }
/*
 *   Fill in terms for each phase
 */
  for (i = 0; i < inv_ptr->count_phases; i++)
  {
    if (inv_ptr->phases[i].count_isotopes <= 0)
      continue;
    isotope_ptr = inv_ptr->phases[i].isotopes;
    for (j = 0; j < inv_ptr->phases[i].count_isotopes; j++)
    {
      if (isotope_ptr[j].primary == primary_ptr &&
	  isotope_ptr[j].isotope_number == isotope_number)
      {
	/* term for alpha phase unknowns */
	column = col_phases + i;
	array[row * max_column_count + column] =
	  isotope_ptr[j].ratio * isotope_ptr[j].coef;
	/* term for phase isotope uncertainty unknown */
	column = col_phase_isotopes + i * inv_ptr->count_isotopes + n;
	array[row * max_column_count + column] = isotope_ptr[j].coef;
	break;
      }
    }

  }
  return OK;
}

/* ---------------------------------------------------------------------- */
int
count_isotope_unknowns (struct inverse *inv_ptr,
			struct isotope **isotope_unknowns)
/* ---------------------------------------------------------------------- */
{
/*
 *  Go through elements for which isotope balances are requested
 *  and make a array of isotope structures
 *  return total number of isotope unknowns and structure array
 */
  int i, k;
  LDBLE isotope_number;
  struct master *primary_ptr;
  int count_isotopes;
  struct isotope *isotopes;

  if (inv_ptr->count_isotopes == 0)
  {
    *isotope_unknowns = NULL;
    return (0);
  }
  isotopes =
    (struct isotope *) PHRQ_malloc ((size_t) sizeof (struct isotope));
  if (isotopes == NULL)
    malloc_error ();
  count_isotopes = 0;

  for (i = 0; i < inv_ptr->count_isotopes; i++)
  {
    primary_ptr = master_bsearch (inv_ptr->isotopes[i].elt_name);
    isotope_number = inv_ptr->isotopes[i].isotope_number;
    if (primary_ptr == NULL)
    {
      sprintf (error_string, "Element not found for isotope calculation: %s.",
	       inv_ptr->isotopes[i].elt_name);
      error_msg (error_string, CONTINUE);
      input_error++;
      break;
    }
    if (primary_ptr->primary != TRUE)
    {
      sprintf (error_string, "Isotope mass-balance may only be used"
	       " for total element concentrations.\n"
	       "Secondary species not allowed: %s.",
	       inv_ptr->isotopes[i].elt_name);
      error_msg (error_string, CONTINUE);
      input_error++;
      break;
    }

    /* nonredox element */
    if (primary_ptr->s->secondary == NULL)
    {
      isotopes =
	(struct isotope *) PHRQ_realloc (isotopes,
					 (size_t) (count_isotopes +
						   1) *
					 sizeof (struct isotope));
      if (isotopes == NULL)
	malloc_error ();
      isotopes[count_isotopes].primary = primary_ptr;
      isotopes[count_isotopes].master = primary_ptr;
      isotopes[count_isotopes].isotope_number = isotope_number;
      isotopes[count_isotopes].elt_name = primary_ptr->elt->name;
      count_isotopes++;

      /* redox element */
    }
    else
    {

      /* find master */
      for (k = 0; k < count_master; k++)
      {
	if (master[k] == primary_ptr)
	  break;
      }

      /* sum all secondary for master */
      k++;
      for (; k < count_master; k++)
      {
	if (master[k]->elt->primary != primary_ptr)
	  break;
	isotopes =
	  (struct isotope *) PHRQ_realloc (isotopes,
					   (size_t) (count_isotopes +
						     1) *
					   sizeof (struct isotope));
	if (isotopes == NULL)
	  malloc_error ();
	isotopes[count_isotopes].primary = primary_ptr;
	isotopes[count_isotopes].master = master[k];
	isotopes[count_isotopes].isotope_number = isotope_number;
	isotopes[count_isotopes].elt_name = master[k]->elt->name;
	count_isotopes++;
      }
    }
  }
  *isotope_unknowns = isotopes;
  return (count_isotopes);
}

/* ---------------------------------------------------------------------- */
int
check_isotopes (struct inverse *inv_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *  Go through elements for which isotope balances are requested
 *  and make sure each solution has isotope ratios defined
 */
  int i, ii, j, k, l;
  int err, found_isotope;
  LDBLE isotope_number;
  struct master *master_ptr, *primary_ptr;
  struct solution *solution_ptr;
  struct phase *phase_ptr;
  char token[MAX_LENGTH];

/*
 *  Check solutions for necessary isotope data
 */
  for (j = 0; j < inv_ptr->count_solns; j++)
  {
    solution_ptr = solution_bsearch (inv_ptr->solns[j], &i, TRUE);
    xsolution_zero ();
    add_solution (solution_ptr, 1.0, 1.0);
/*
 *   Go through inverse isotopes and make sure isotope data for each solution
 *   inv_ptr->isotopes has elements; inv_ptr->i_u has redox states and uncertainties
 */
    for (i = 0; i < inv_ptr->count_isotopes; i++)
    {
      err = FALSE;
      primary_ptr = master_bsearch (inv_ptr->isotopes[i].elt_name);
      isotope_number = inv_ptr->isotopes[i].isotope_number;
      found_isotope = FALSE;
      for (k = 0; k < solution_ptr->count_isotopes; k++)
      {
	if (solution_ptr->isotopes[k].primary == primary_ptr &&
	    solution_ptr->isotopes[k].isotope_number == isotope_number)
	{
	  found_isotope = TRUE;
	  break;
	}
      }
      if (found_isotope == TRUE)
	continue;

      /* did not find isotope, which is ok if element not in solution */
      if (primary_ptr == s_h2o->primary || primary_ptr == s_hplus->primary)
      {
	err = TRUE;
      }
      else if (primary_ptr->total > 0)
      {
	err = TRUE;
      }
      if (err == TRUE)
      {
	sprintf (error_string,
		 "In solution %d, isotope ratio(s) are needed for element: %g%s.",
		 solution_ptr->n_user, (double) isotope_number,
		 primary_ptr->elt->name);
	error_msg (error_string, CONTINUE);
	input_error++;
	continue;
      }
    }
/*
 *   Go through solution isotopes and set uncertainties
 */
    for (k = 0; k < solution_ptr->count_isotopes; k++)
    {
      solution_ptr->isotopes[k].x_ratio_uncertainty = NAN;
/*
 *  Search for secondary or primary master in inverse uncertainties
 */
      ii = -1;
      for (i = 0; i < inv_ptr->count_i_u; i++)
      {
	master_ptr = master_bsearch (inv_ptr->i_u[i].elt_name);
	if (master_ptr == solution_ptr->isotopes[k].master)
	{
	  ii = i;
	  break;
	}
	if (master_ptr == solution_ptr->isotopes[k].primary)
	{
	  ii = i;
	}
      }
      /* solution isotope data not being used in inverse */
      if (ii == -1)
	continue;

      i = ii;
      /* use inverse-defined uncertainties first */
      if (j < inv_ptr->i_u[i].count_uncertainties
	  && inv_ptr->i_u[i].uncertainties[j] != NAN)
      {
	solution_ptr->isotopes[k].x_ratio_uncertainty =
	  inv_ptr->i_u[i].uncertainties[j];

	/* use solution-defined uncertainties second */
      }
      else if (solution_ptr->isotopes[k].ratio_uncertainty != NAN)
      {
	solution_ptr->isotopes[k].x_ratio_uncertainty =
	  solution_ptr->isotopes[k].ratio_uncertainty;
	/* use isotope defaults third */
      }
      else
      {
	sprintf (token, "%g%s",
		 (double) solution_ptr->isotopes[k].isotope_number,
		 solution_ptr->isotopes[k].elt_name);
	for (l = 0; l < count_iso_defaults; l++)
	{
	  if (strcmp (token, iso_defaults[l].name) == 0)
	  {
	    solution_ptr->isotopes[k].x_ratio_uncertainty =
	      iso_defaults[l].uncertainty;
	    sprintf (error_string,
		     "Solution %d,  element %g%s: default isotope ratio uncertainty is used, %g.",
		     solution_ptr->n_user,
		     (double) solution_ptr->isotopes[k].isotope_number,
		     solution_ptr->isotopes[k].elt_name,
		     (double) solution_ptr->isotopes[k].x_ratio_uncertainty);
	    warning_msg (error_string);
	    break;
	  }
	}
      }
      if (solution_ptr->isotopes[k].x_ratio_uncertainty == NAN)
      {
	sprintf (error_string,
		 "In solution %d, isotope ratio uncertainty is needed for element: %g%s.",
		 solution_ptr->n_user,
		 (double) solution_ptr->isotopes[k].isotope_number,
		 solution_ptr->isotopes[k].elt_name);
	error_msg (error_string, CONTINUE);
	input_error++;
      }
    }
  }
/*
 *  Check phases for necessary isotope data
 */
  for (j = 0; j < inv_ptr->count_phases; j++)
  {
    for (i = 0; i < inv_ptr->count_isotopes; i++)
    {
      primary_ptr = master_bsearch (inv_ptr->isotopes[i].elt_name);
      isotope_number = inv_ptr->isotopes[i].isotope_number;
      found_isotope = FALSE;
      for (k = 0; k < inv_ptr->phases[j].count_isotopes; k++)
      {
	if (inv_ptr->phases[j].isotopes[k].primary == primary_ptr &&
	    inv_ptr->phases[j].isotopes[k].isotope_number == isotope_number)
	{
	  found_isotope = TRUE;
	  break;
	}
      }
      if (found_isotope == TRUE)
	continue;

      /* did not find isotope, which is ok if element not in solution */
      phase_ptr = inv_ptr->phases[j].phase;
      k = 0;
      while (phase_ptr->next_elt[k].elt != NULL)
      {
	if (phase_ptr->next_elt[k].elt->primary == primary_ptr)
	{
	  if (s_hplus->primary == primary_ptr ||
	      s_h2o->primary == primary_ptr)
	  {
	    k++;
	    continue;
	  }
	  else
	  {
	    sprintf (error_string,
		     "In phase %s, isotope ratio(s) are needed for element: %g%s.",
		     phase_ptr->name, (double) isotope_number,
		     primary_ptr->elt->name);
	    error_msg (error_string, CONTINUE);
	    input_error++;
	    break;
	  }
	}
	k++;
      }
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
phase_isotope_inequalities (struct inverse *inv_ptr)
/* ---------------------------------------------------------------------- */
{
  int i, j, k;
  int column;
  char token[MAX_LENGTH];
  if (inv_ptr->count_isotopes <= 0)
    return OK;
  for (i = 0; i < inv_ptr->count_phases; i++)
  {
    if (inv_ptr->phases[i].count_isotopes <= 0)
      continue;

    for (j = 0; j < inv_ptr->phases[i].count_isotopes; j++)
    {
      /* find index number */
      for (k = 0; k < inv_ptr->count_isotopes; k++)
      {
	if (inv_ptr->phases[i].isotopes[j].elt_name ==
	    inv_ptr->isotopes[k].elt_name
	    && inv_ptr->phases[i].isotopes[j].isotope_number ==
	    inv_ptr->isotopes[k].isotope_number)
	{
	  break;
	}
      }
      if (k >= inv_ptr->count_isotopes)
	break;
      column = col_phase_isotopes + i * inv_ptr->count_isotopes + k;
/*
 *   zero column if uncertainty is zero
 */
      if (inv_ptr->phases[i].isotopes[j].ratio_uncertainty == 0)
      {
	for (k = 0; k < count_rows; k++)
	{
	  array[k * max_column_count + column] = 0.0;
	}
	continue;
      }

/* 
 *   optimization
 */
      array[(column - col_epsilon) * max_column_count + column] =
	SCALE_EPSILON / inv_ptr->phases[i].isotopes[j].ratio_uncertainty;
/*
 *   two inequalities to account for absolute value
 */
      /* for phases constrained to precipitate */
      if (inv_ptr->phases[i].constraint == PRECIPITATE)
      {
	array[count_rows * max_column_count + col_phases + i] =
	  inv_ptr->phases[i].isotopes[j].ratio_uncertainty;
	array[count_rows * max_column_count + column] = 1.0;
	sprintf (token, "%s %s", inv_ptr->phases[i].phase->name, "iso pos");
	row_name[count_rows] = string_hsave (token);
	count_rows++;

	array[count_rows * max_column_count + col_phases + i] =
	  inv_ptr->phases[i].isotopes[j].ratio_uncertainty;
	array[count_rows * max_column_count + column] = -1.0;
	sprintf (token, "%s %s", inv_ptr->phases[i].phase->name, "iso neg");
	row_name[count_rows] = string_hsave (token);
	count_rows++;

	/* for phases constrained to dissolve */
      }
      else if (inv_ptr->phases[i].constraint == DISSOLVE)
      {
	array[count_rows * max_column_count + col_phases + i] =
	  -inv_ptr->phases[i].isotopes[j].ratio_uncertainty;
	array[count_rows * max_column_count + column] = -1.0;
	sprintf (token, "%s %s", inv_ptr->phases[i].phase->name, "iso pos");
	row_name[count_rows] = string_hsave (token);
	count_rows++;

	array[count_rows * max_column_count + col_phases + i] =
	  -inv_ptr->phases[i].isotopes[j].ratio_uncertainty;
	array[count_rows * max_column_count + column] = 1.0;
	sprintf (token, "%s %s", inv_ptr->phases[i].phase->name, "iso neg");
	row_name[count_rows] = string_hsave (token);
	count_rows++;

	/* Error if phase is not constrained */
      }
      else
      {
	sprintf (error_string,
		 "In isotope calculations, all phases containing isotopes must be"
		 " constrained.\nPhase %s is not constrained.\n",
		 inv_ptr->phases[i].phase->name);
	error_msg (error_string, CONTINUE);
	input_error++;
	continue;
      }
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
write_optimize_names (struct inverse *inv_ptr)
/* ---------------------------------------------------------------------- */
{
  int i, j, row;
  char token[MAX_LENGTH];
  row = 0;
/*
 *   epsilons for analytical data
 */
  for (j = 0; j < inv_ptr->count_elts; j++)
  {
    for (i = 0; i < inv_ptr->count_solns; i++)
    {
      sprintf (token, "%s %s %d", "optimize",
	       inv_ptr->elts[j].master->elt->name, inv_ptr->solns[i]);
      row_name[row] = string_hsave (token);
      row++;
    }
  }
/*
 *   pH
 */
  if (carbon > 0)
  {
    for (i = 0; i < inv_ptr->count_solns; i++)
    {
      sprintf (token, "%s %s %d", "optimize", "pH", inv_ptr->solns[i]);
      row_name[row] = string_hsave (token);
      row++;
    }
  }
/*
 *   water
 */
  sprintf (token, "%s %s", "optimize", "water");
  row_name[row] = string_hsave (token);
  row++;
/*
 *   solution isotopes
 */
  for (i = 0; i < inv_ptr->count_solns; i++)
  {
    for (j = 0; j < inv_ptr->count_isotope_unknowns; j++)
    {
      sprintf (token, "%s %d%s %d", "optimize",
	       (int) inv_ptr->isotope_unknowns[j].isotope_number,
	       inv_ptr->isotope_unknowns[j].elt_name, inv_ptr->solns[i]);
      row_name[row] = string_hsave (token);
      row++;
    }
  }
/*
 *   phase isotopes
 */

  for (i = 0; i < inv_ptr->count_phases; i++)
  {
    for (j = 0; j < inv_ptr->count_isotopes; j++)
    {
      sprintf (token, "%s %s %d%s", "optimize",
	       inv_ptr->phases[i].phase->name,
	       (int) inv_ptr->isotopes[j].isotope_number,
	       inv_ptr->isotopes[j].elt_name);
      row_name[row] = string_hsave (token);
      row++;
    }
  }
  return OK;
}
/* ---------------------------------------------------------------------- */
void
dump_netpath (struct inverse *inverse_ptr)
/* ---------------------------------------------------------------------- */
{
  int i, j, l;
  char string[MAX_LENGTH];
  char *ptr;

  /*if (netpath_dumped) return;*/
  if (inverse_ptr->netpath == NULL) return;

  /* open file */
  strcpy(string, inverse_ptr->netpath);
  if (replace(".lon", ".lon", string) != TRUE)
  {
    strcat(string, ".lon");
  }
  netpath_file = fopen(string, "w");
  if (netpath_file == NULL) {
    sprintf (error_string, "Can't open file, %s.", inverse_ptr->netpath);
    error_msg (error_string, STOP);
  }
  add_to_file("netpath.fil", inverse_ptr->netpath);

  /* Header */
  fprintf(netpath_file, "2.14                                                       # File format\n");

  /* write out each solution */
  for (i = 0; i < count_solution; i++) {
    if (solution[i]->n_user < 0) continue;

    /* flags and description */
    ptr = solution[i]->description;
    j = copy_token(string, &ptr, &l);
    if (j != EMPTY) 
    {
      sprintf(string, "%s", solution[i]->description);
    } else {
      sprintf(string, "Solution %d", solution[i]->n_user);
    }
    fprintf(netpath_file, "4020%s\n", string);

    /* lat/lon */
    fprintf(netpath_file, "                                                           # Lat/lon\n");
    
    /* well number */
    fprintf(netpath_file, "%15d                                            # Well number\n", solution[i]->n_user);

    /* total number of wells */
    fprintf(netpath_file, "%15d                                            # Total wells\n", count_solution);

    /* address */
    fprintf(netpath_file, "                                                           # Address1\n");
    fprintf(netpath_file, "                                                           # Address2\n");
    fprintf(netpath_file, "                                                           # Address3\n");
    fprintf(netpath_file, "                                                           # Address4\n");
    fprintf(netpath_file, "                                                           # Address5\n");

    /* temperature */
    fprintf(netpath_file, "%15g                                            # Temperature\n", solution[i]->tc);

    /* pH */
    fprintf(netpath_file, "%15g                                            # pH\n", solution[i]->ph);

    /* DO */
    print_total(netpath_file, solution[i], "O(0)", "Dissolved Oxygen");

    /* TDIC */
    print_total(netpath_file, solution[i], "C(4)", "TDIC");

    /* Tritium */
    print_isotope(netpath_file, solution[i], "3H(1)", "Tritium");

    /* H2S */
    print_total(netpath_file, solution[i], "S(-2)", "H2S");

    /* Calcium */
    print_total(netpath_file, solution[i], "Ca", "Calcium");

    /* Eh */
    fprintf(netpath_file, "%15g                                            # Eh\n", 0.059 * solution[i]->solution_pe);

    /* Magnesium */
    print_total(netpath_file, solution[i], "Mg", "Magnesium");

    /* Sodium */
    print_total(netpath_file, solution[i], "Na", "Sodium");

    /* Potassium */
    print_total(netpath_file, solution[i], "K", "Potassium");

    /* Chloride */
    print_total(netpath_file, solution[i], "Cl", "Chloride");

    /* Sulfate */
    print_total(netpath_file, solution[i], "S(6)", "Sulfate");

    /* Fluoride */
    print_total(netpath_file, solution[i], "F", "Fluoride");

    /* Silica */
    print_total(netpath_file, solution[i], "Si", "Silica");

    /* Bromide */
    print_total(netpath_file, solution[i], "Br", "Bromide");

    /* Boron */
    print_total(netpath_file, solution[i], "B", "Boron");

    /* Barium */
    print_total(netpath_file, solution[i], "Ba", "Barium");

    /* Lithium */
    print_total(netpath_file, solution[i], "Li", "Lithium");

    /* Strontium */
    print_total(netpath_file, solution[i], "Sr", "Strontium");

    /* Iron */
    print_total_multi(netpath_file, solution[i], "Iron", "Fe", "Fe(2)", "Fe(3)", "", "");


    /* Manganese */
    print_total_multi(netpath_file, solution[i], "Manganese", "Mn", "Mn(2)", "Mn(3)", "Mn(6)", "Mn(7)");

    /* Nitrate */
    print_total(netpath_file, solution[i], "N(5)", "Nitrate");

    /* Ammonium */
    print_total_multi(netpath_file, solution[i], "Ammonium", "N(-3)", "Amm", "", "", "");

    /* Phosphate */
    print_total(netpath_file, solution[i], "P", "Phosphate");

    /* DOC */
    print_total_multi(netpath_file, solution[i], "DOC", "Fulvate", "Humate", "", "", "");

    /* Sp. Cond. */
    fprintf(netpath_file, "                                                           # Sp. Cond.\n");

    /* Density */
    fprintf(netpath_file, "                                                           # Density\n");

    /* Delta C-13 TDIC */
    print_isotope(netpath_file, solution[i], "13C(4)", "Delta C-13 TDIC");

    /* C-14 TDIC */
    print_isotope(netpath_file, solution[i], "14C(4)", "C-14 TDIC");

    /* Delta S-34 (SO4) */
    print_isotope(netpath_file, solution[i], "34S(6)", "Delta S-34 (SO4)");

    /* Delta S-34 (H2S) */
    print_isotope(netpath_file, solution[i], "34S(-2)", "Delta S-34 (H2S)");

    /* Delta Deuterium */
    print_isotope(netpath_file, solution[i], "2H(1)", "Delta Deuterium");

    /* Delta O-18 */
    print_isotope(netpath_file, solution[i], "18O(-2)", "Delta O-18");

    /* CH4 (aq) */
    print_total(netpath_file, solution[i], "C(-4)", "CH4 (aq)");

    /* Sr 87/86 */
    print_isotope(netpath_file, solution[i], "87Sr", "Sr 87/86");

    /* Al */
    print_total(netpath_file, solution[i], "Al", "Alumninum");

    /* N2 (aq) */
    print_total(netpath_file, solution[i], "N(0)", "N2 (aq)");

    /* N-15 of N2 (aq) */
    print_isotope(netpath_file, solution[i], "15N(0)", "N-15 of N2 (aq)");

    /* N-15 of Nitrate */
    print_isotope(netpath_file, solution[i], "15N(5)", "N-15 of Nitrate");

    /* N-15 of Ammonium */
    print_isotope(netpath_file, solution[i], "15N(-3)", "N-15 of Ammonium");

    /* Formation */
    fprintf(netpath_file, "                                                           # Formation\n");

  }
  if (netpath_file != NULL) 
  {
    fclose(netpath_file);
    netpath_file = NULL;
  }
  return;
}
/* ---------------------------------------------------------------------- */
struct conc *
get_total (struct solution *solution_ptr, char *elt)
/* ---------------------------------------------------------------------- */
{
  int i;

  for (i = 0; solution_ptr->totals[i].description != NULL; i++)
  {
    if (strcmp(elt, solution_ptr->totals[i].description) == 0) return(&solution_ptr->totals[i]);
  }
  return(NULL);

}
/* ---------------------------------------------------------------------- */
struct isotope *
get_isotope (struct solution *solution_ptr, char *elt)
/* ---------------------------------------------------------------------- */
{
  int i;

  for (i = 0; i < solution_ptr->count_isotopes; i++)
  {
    if (strcmp(elt, solution_ptr->isotopes[i].isotope_name) == 0) return(&solution_ptr->isotopes[i]);
  }
  return(NULL);
}

/* ---------------------------------------------------------------------- */
void
print_total (FILE *netpath_file, struct solution *solution_ptr, char *elt, char *string)
/* ---------------------------------------------------------------------- */
{
  struct conc *tot_ptr;
  tot_ptr = get_total(solution_ptr, elt);
  if (tot_ptr == NULL) {
      fprintf(netpath_file, "                                                           # %s\n", string);
    } else {
      fprintf(netpath_file, "%15g                                            # %s\n", 
	      1000 * tot_ptr->moles / solution_ptr->mass_water, string);
    }
}

/* ---------------------------------------------------------------------- */
void
print_isotope (FILE *netpath_file, struct solution *solution_ptr, char *elt, char *string)
/* ---------------------------------------------------------------------- */
{
  struct isotope *iso_ptr;
  iso_ptr = get_isotope(solution_ptr, elt);
  if (iso_ptr == NULL) {
      fprintf(netpath_file, "                                                           # %s\n", string);
    } else {
      fprintf(netpath_file, "%15g                                            # %s\n", 
	      iso_ptr->ratio, string);
    }
}
/* ---------------------------------------------------------------------- */
void
print_total_multi (FILE *netpath_file, struct solution *solution_ptr, char *string, char *elt0, char *elt1, char *elt2, char *elt3, char *elt4 )
/* ---------------------------------------------------------------------- */
{
  char elts[5][MAX_LENGTH];
  struct conc *tot_ptr;
  double sum;
  int i, found;

  strcpy(elts[0], elt0);
  strcpy(elts[1], elt1);
  strcpy(elts[2], elt2);
  strcpy(elts[3], elt3);
  strcpy(elts[4], elt4);


  sum = 0;
  found = FALSE;
  for (i = 0; i < 5; i++)
  {
    tot_ptr = get_total(solution_ptr, elts[i]);
    if (tot_ptr == NULL) {
      continue;
    } else {
      sum += tot_ptr->moles;
      found = TRUE;
    }
  }
  if (found != TRUE)
  {
    fprintf(netpath_file, "                                                           # %s\n", string);
  } else 
  {
    fprintf(netpath_file, "%15g                                            # %s\n", 
	    1000 * sum / solution_ptr->mass_water, string);
  }
  return;
}

/* ---------------------------------------------------------------------- */
int
dump_netpath_pat (struct inverse *inv_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Prints model
 */
  int i, j, k;
  struct solution *solution_ptr, *solution_ptr_orig;
  struct master *master_ptr;
  LDBLE d1, d2, d3;
  char string[MAX_LENGTH], string1[MAX_LENGTH], token[MAX_LENGTH];
  char *ptr, *str_ptr;
  int l;
  double sum, sum1, sum_iso, d;
  double *array_save, *delta_save;
  int count_unknowns_save, max_row_count_save, max_column_count_save, temp, count_current_solutions;
  int solnmap[10][2];
  struct isotope *isotope_ptr;
  FILE *model_file;
  struct elt_list *next_elt;
  int exch, column;
  double f;
  struct rxn_token *rxn_ptr;
/*
 *   print solution data, epsilons, and revised data
 */
  if (inv_ptr->pat == NULL) return(OK);

  array_save = array;
  delta_save = delta;
  count_unknowns_save = count_unknowns;
  max_row_count_save = max_row_count;
  max_column_count_save = max_column_count;

  array = NULL;
  delta = NULL;
  count_unknowns = 0;
  max_row_count = 0;
  max_column_count = 0;

  count_current_solutions = 0;
  count_inverse_models++;

  for (i = 0; i < inv_ptr->count_solns; i++)
  {
    if (equal (delta1[i], 0.0, toler) == TRUE)
      continue;
    solution_ptr_orig = solution_bsearch (inv_ptr->solns[i], &j, TRUE);

    solution_duplicate(solution_ptr_orig->n_user, -6);
    solution_ptr = solution_bsearch (-6, &j, TRUE);
    xsolution_zero ();

    /* Adjust pH */
    if (inv_ptr->carbon == TRUE)
    {
      d1 = solution_ptr->ph;
      d2 = delta1[col_ph + i] / delta1[i];
      d3 = d1 + d2;
      solution_ptr->ph = d3;
    }

    /* put original totals in master */
    for (j = 0; solution_ptr->totals[j].description != NULL; j++)
    {
      master_ptr = master_bsearch (solution_ptr->totals[j].description);
      master_ptr->total = solution_ptr->totals[j].moles;
    }

    /* ignore alkalinity */
    /*master_alk->total = solution_ptr->total_alkalinity;*/

    /* update total in master */
    for (j = 0; j < inv_ptr->count_elts; j++)
    {
      if (inv_ptr->elts[j].master->s == s_eminus)
	continue;
      d1 = inv_ptr->elts[j].master->total;
      d2 = delta1[col_epsilon + j * inv_ptr->count_solns + i] / delta1[i];
      d3 = d1 + d2;
      inv_ptr->elts[j].master->total = d3;
    }

    /* put updated total back in solution */
    for (j = 0; solution_ptr->totals[j].description != NULL; j++)
    {
      master_ptr = master_bsearch (solution_ptr->totals[j].description);
      solution_ptr->totals[j].moles = master_ptr->total;
    }


    /* update isotopes in solution */
    if (inv_ptr->count_isotopes > 0)
    {
      /* adjustments to solution isotope composition */
      for (j = 0; j < inv_ptr->count_isotope_unknowns; j++)
      {
	for (k = 0; k < solution_ptr->count_isotopes; k++)
	{
	  if (inv_ptr->isotope_unknowns[j].elt_name !=
	      solution_ptr->isotopes[k].elt_name ||
	      inv_ptr->isotope_unknowns[j].isotope_number !=
	      solution_ptr->isotopes[k].isotope_number)
	    continue;
	  d1 = solution_ptr->isotopes[k].ratio;
	  d2 =
	    delta1[col_isotopes + i * inv_ptr->count_isotope_unknowns +
		   j] / delta1[i];
	  d3 = d1 + d2;
	  solution_ptr->isotopes[k].ratio = d3;
	}
      }
    }

    set_initial_solution(-6, -7);
    temp = pr.all;
    pr.all = FALSE;
    initial_solutions(FALSE);
    pr.all = temp;
    solution_ptr = solution_bsearch(-7, &j, TRUE);

    /* Header */
    ptr = solution_ptr_orig->description;
    if (copy_token(string, &ptr, &l) != EMPTY)
    { 
      fprintf(netpath_file, "%d. %s\n", count_inverse_models, solution_ptr_orig->description);
    } else {
      fprintf(netpath_file, "%d. Solution %d\n", count_inverse_models, solution_ptr_orig->n_user);
    }

    /* bookkeeping */
    count_pat_solutions++;
    solnmap[count_current_solutions][0] = solution_ptr_orig->n_user;
    solnmap[count_current_solutions][1] = count_pat_solutions;
    count_current_solutions++;

    /* Dump info to .pat file */
    print_total_pat(netpath_file, "C", "C");
    print_total_pat(netpath_file, "S", "S");
    print_total_pat(netpath_file, "Ca", "CA");
    print_total_pat(netpath_file, "Al", "AL");
    print_total_pat(netpath_file, "Mg", "MG");
    print_total_pat(netpath_file, "Na", "NA");
    print_total_pat(netpath_file, "K", "K");
    print_total_pat(netpath_file, "Cl", "CL");
    print_total_pat(netpath_file, "F", "F");
    print_total_pat(netpath_file, "Si", "SI");
    print_total_pat(netpath_file, "Br", "BR");
    print_total_pat(netpath_file, "B", "B");
    print_total_pat(netpath_file, "Ba", "BA");
    print_total_pat(netpath_file, "Li", "LI");
    print_total_pat(netpath_file, "Sr", "SR");
    print_total_pat(netpath_file, "Fe", "FE");
    print_total_pat(netpath_file, "Mn", "MN");
    print_total_pat(netpath_file, "N", "N");
    print_total_pat(netpath_file, "P", "P");
    fprintf(netpath_file, "%14g     # TEMP\n", solution_ptr->tc);
    print_total_pat(netpath_file, "S(-2)", "H2S");
    print_total_pat(netpath_file, "S(6)", "SO4");

    /* N15 */
    sum_iso = 0;
    sum = 0;
    for (k = 0; k < solution_ptr->count_isotopes; k++)
    {
      if (strstr(solution_ptr->isotopes[k].isotope_name, "15N") != NULL)
      {
	d = total(solution_ptr->isotopes[k].elt_name);
	sum_iso += solution_ptr->isotopes[k].ratio * d;
	sum += d;
      }
    }
    if (sum == 0)
    {
      fprintf(netpath_file, "%14g*    # N15\n", sum);
    } else
    {
      fprintf(netpath_file, "%14g     # N15\n", sum_iso/sum);
    }    

    /* RS of N */
    sum = 0;
    sum = total("N(-3)")*-3 + total("N(0)")*0 + total("N(3)")*3 + total("N(5)")*5 ;
    sum1 = total("N(-3)") + total("N(0)") + total("N(3)") + total("N(5)");
    if (sum1 == 0) 
    {
      fprintf(netpath_file, "%14g*    # RS of N\n", sum1);
    } else 
    {
      fprintf(netpath_file, "%14g     # RS of N\n", sum / sum1);
    }

    /* DOX */
    print_total_pat(netpath_file, "O(0)", "DOX");

    /*HCO3*/
    d = 1000 * sum_match_species("*HCO3*", "C");
    if (d == 0.0) {
      fprintf(netpath_file, "%14g*    # HCO3\n", d);
    } else 
    {
      fprintf(netpath_file, "%14g     # HCO3\n", d);
    }

    /* pH */
    fprintf(netpath_file, "%14g     # PH\n", solution_ptr->ph);

    /*H2CO3**/
    d = 1000 * (molality("H2CO3") + molality("CO2"));
    if (d == 0.0) {
      fprintf(netpath_file, "%14g*    # H2CO3\n", d);
    } else 
    {
      fprintf(netpath_file, "%14g     # H2CO3\n", d);
    }

    /*CO3*/
    d = sum_match_species("*CO3*", "C");
    d -= sum_match_species("*HCO3*", "C");
    d *= 1000.0;
    if (d == 0.0) {
      fprintf(netpath_file, "%14g*     # CO3\n", d);
    } else 
    {
      fprintf(netpath_file, "%14g     # CO3\n", d);
    }

    /* CARBONATES */
    print_total_pat(netpath_file, "C(4)", "CARBONATES");
    print_total_pat(netpath_file, "Fe(2)", "FE2+");
    print_total_pat(netpath_file, "Fe(3)", "FE3+");
    print_total_pat(netpath_file, "Mn(2)", "MN2+");
    print_total_pat(netpath_file, "Mn(3)", "MN3+");
    print_total_pat(netpath_file, "Mn(6)", "MN6+");
    print_total_pat(netpath_file, "Mn(7)", "MN7+");
    print_total_pat(netpath_file, "C(-4)", "CH4");
    print_total_pat(netpath_file, "Doc", "DOC");

    /*RS OF DOC*/
    fprintf(netpath_file, "%14g*    # RS OF DOC\n", 0.0);

    /* Blank */
    print_total_pat(netpath_file, "Blank", "BLANK");

    /*C13*/
    sum_iso = 0;
    sum = 0;
    for (k = 0; k < solution_ptr->count_isotopes; k++)
    {
      if (strstr(solution_ptr->isotopes[k].isotope_name, "13C") != NULL)
      {
	d = total(solution_ptr->isotopes[k].elt_name);
	sum_iso += solution_ptr->isotopes[k].ratio * d;
	sum += d;
      }
    }
    if (sum == 0)
    {
      fprintf(netpath_file, "%14g*    # C13\n", sum);
    } else
    {
      fprintf(netpath_file, "%14g     # C13\n", sum_iso/sum);
    }    

    /*C14*/    
    sum_iso = 0;
    sum = 0;
    for (k = 0; k < solution_ptr->count_isotopes; k++)
    {
      if (strstr(solution_ptr->isotopes[k].isotope_name, "14C") != NULL)
      {
	d = total(solution_ptr->isotopes[k].elt_name);
	sum_iso += solution_ptr->isotopes[k].ratio * d;
	sum += d;
      }
    }
    if (sum == 0)
    {
      fprintf(netpath_file, "%14g*    # C14\n", sum);
    } else
    {
      fprintf(netpath_file, "%14g     # C14\n", sum_iso/sum);
    }    

    /*SR87*/
    sum_iso = 0;
    sum = 0;
    for (k = 0; k < solution_ptr->count_isotopes; k++)
    {
      if (strstr(solution_ptr->isotopes[k].isotope_name, "87Sr") != NULL)
      {
	d = total(solution_ptr->isotopes[k].elt_name);
	sum_iso += solution_ptr->isotopes[k].ratio * d;
	sum += d;
      }
    }
    if (sum == 0)
    {
      fprintf(netpath_file, "%14g*    # SR87\n", sum);
    } else
    {
      fprintf(netpath_file, "%14g     # SR87\n", sum_iso/sum);
    }    
    
    /*D*/
    sum_iso = 0;
    sum = 0;
    for (k = 0; k < solution_ptr->count_isotopes; k++)
    {
      if (strstr(solution_ptr->isotopes[k].isotope_name, "2H") != NULL)
      {
	d = total(solution_ptr->isotopes[k].elt_name);
	sum_iso += solution_ptr->isotopes[k].ratio * d;
	sum += d;
      }
    }
    if (sum == 0)
    {
      fprintf(netpath_file, "%14g*    # D\n", sum);
    } else
    {
      fprintf(netpath_file, "%14g     # D\n", sum_iso/sum);
    }    

    /*O-18*/
    sum_iso = 0;
    sum = 0;
    for (k = 0; k < solution_ptr->count_isotopes; k++)
    {
      if (strstr(solution_ptr->isotopes[k].isotope_name, "18O") != NULL)
      {
	if (strcmp(solution_ptr->isotopes[k].elt_name, "O(-2)") == 0) 
	{
	  d = solution_ptr->total_o - total("O(0)");
	} else if (strcmp(solution_ptr->isotopes[k].elt_name, "H(1)") == 0) 
	{
	  d = solution_ptr->total_h - total("H(0)");
	} else
	{
	  d = total(solution_ptr->isotopes[k].elt_name);
	}
	sum_iso += solution_ptr->isotopes[k].ratio * d;
	sum += d;
      }
    }
    if (sum == 0)
    {
      fprintf(netpath_file, "%14g*    # O-18\n", sum);
    } else
    {
      fprintf(netpath_file, "%14g     # O-18\n", sum_iso/sum);
    }    

    /*TRITIUM*/
    sum_iso = 0;
    sum = 0;
    for (k = 0; k < solution_ptr->count_isotopes; k++)
    {
      if (strstr(solution_ptr->isotopes[k].isotope_name, "3H") != NULL)
      {
	d = total(solution_ptr->isotopes[k].elt_name);
	sum_iso += solution_ptr->isotopes[k].ratio * d;
	sum += d;
      }
    }
    if (sum == 0)
    {
      fprintf(netpath_file, "%14g*    # TRITIUM\n", sum);
    } else
    {
      fprintf(netpath_file, "%14g     # TRITIUM\n", sum_iso/sum);
    }    

    /*34SSO4*/
    sum_iso = 0;
    sum = 0;
    for (k = 0; k < solution_ptr->count_isotopes; k++)
    {
      if (strstr(solution_ptr->isotopes[k].isotope_name, "34S(6)") != NULL)
      {
	d = total(solution_ptr->isotopes[k].elt_name);
	sum_iso += solution_ptr->isotopes[k].ratio * d;
	sum += d;
      }
    }
    if (sum == 0)
    {
      fprintf(netpath_file, "%14g*    # 34SSO4\n", sum);
    } else
    {
      fprintf(netpath_file, "%14g     # 34SSO4\n", sum_iso/sum);
    }    

    /*34SH2S*/
    sum_iso = 0;
    sum = 0;
    for (k = 0; k < solution_ptr->count_isotopes; k++)
    {
      if (strstr(solution_ptr->isotopes[k].isotope_name, "34S(-2)") != NULL)
      {
	d = total(solution_ptr->isotopes[k].elt_name);
	sum_iso += solution_ptr->isotopes[k].ratio * d;
	sum += d;
      }
    }
    if (sum == 0)
    {
      fprintf(netpath_file, "%14g*    # 34SH2S\n", sum);
    } else
    {
      fprintf(netpath_file, "%14g     # 34SH2S\n", sum_iso/sum);
    }    

    /* Well number */
    fprintf(netpath_file, "%14d     # Well number\n", count_pat_solutions);
  }
  free_model_allocs();
  array = array_save;
  delta = delta_save;
  count_unknowns = count_unknowns_save;
  max_row_count = max_row_count_save;
  max_column_count = max_column_count_save;

/*
 * Open model file
 */
  strcpy(string, inv_ptr->pat);
  replace(".pat","", string);
  string_trim(string);
  sprintf(string1, "%s-%d.mod",string, count_inverse_models);
  model_file = fopen(string1, "w");
  if (model_file == NULL) {
    sprintf (error_string, "Can't open file, %s.", string);
    error_msg (error_string, STOP);
  }
  add_to_file("model.fil", string1);


/*
 * Write header
 */
  fprintf(model_file,"%s\n", string);

/*
 * Write well numbers
 */
  for (i = 0; i < count_current_solutions; i++) 
  {
    fprintf(model_file,"%3d", solnmap[i][1]);
  }
  fprintf(model_file,"\n");
/*
 * Write elements
 */
  xsolution_zero();
  for (j = 0; j < count_master; j++)
  {
    master[j]->in = FALSE;
  }
  for (j = 0; j < inv_ptr->count_elts; j++)
  {
    master_ptr = inv_ptr->elts[j].master;
    master_ptr = master_ptr->elt->primary;
    if (strcmp(master_ptr->elt->name, "Alkalinity") == 0) continue;
    if (strcmp(master_ptr->elt->name, "H") == 0)  continue;
    if (strcmp(master_ptr->elt->name, "O") == 0) continue;
    if (strcmp(master_ptr->elt->name, "X") == 0) continue;
    if (strcmp(master_ptr->elt->name, "E") == 0) continue;
    master_ptr->in = TRUE;
  }
  for (j = 0; j < count_master; j++)
  {
    if (master[j]->in == TRUE) 
    {
      strcpy(string, master[j]->elt->name);
      str_toupper(string);
      fprintf(model_file," %-2s", string);
    }
  }
  fprintf(model_file," %-2s", "RS");
/*
 * Add isotope mole balance
 */
  for (j = 0; j < inv_ptr->count_isotopes; j++)
  {
    sprintf(string, "%d%s", (int) inv_ptr->isotopes[j].isotope_number, inv_ptr->isotopes[j].elt_name);
    if (strcmp(string,"13C") == 0) fprintf(model_file," %-2s", "I1");
    if (strcmp(string,"14C") == 0) fprintf(model_file," %-2s", "I2");
    if (strcmp(string,"34S") == 0) fprintf(model_file," %-2s", "I3");
    if (strcmp(string,"87Sr") == 0) fprintf(model_file," %-2s", "I4");
    if (strcmp(string,"15N") == 0) fprintf(model_file," %-2s", "I9");
    if (strcmp(string,"2H") == 0) fprintf(model_file," %-2s", "D ");
    if (strcmp(string,"3H") == 0) fprintf(model_file," %-2s", "TR");
    if (strcmp(string,"18O") == 0) fprintf(model_file," %-2s", "18");
  }    

  /* end of element line */
  fprintf(model_file,"\n");

/*
 * Write phase information
 */
  for (i = 0; i < inv_ptr->count_phases; i++)
  {
    j = col_phases + i;
    /* skip if not in model */
/*    if (equal (delta1[j], 0.0, toler) == TRUE) continue;*/

    /* Do not include Na exchange phase */
    if (strcmp_nocase(inv_ptr->phases[i].name, "NaX") == 0) continue;
/*
 * Determine if exchange reaction
 */
    exch = FALSE;
    for (next_elt = inv_ptr->phases[i].phase->next_elt; next_elt->elt != NULL; next_elt++)
    {
      if (strcmp(next_elt->elt->name, "X") == 0) {
	exch = TRUE;
	break;
      }
    }
/*
 * Write phase name and constraints
 */
    strncpy(string, inv_ptr->phases[i].name, (size_t) 8);
    str_ptr = string_pad(string, 10);
    str_ptr[10] = '\0';
    if (inv_ptr->phases[i].force == TRUE)
    {
      str_ptr[8] = 'F';
    } else 
    {
      str_ptr[8] = ' ';
    }
    switch (inv_ptr->phases[i].constraint)
    {
    case EITHER:
      str_ptr[9] = ' ';
      break;
    case PRECIPITATE:
      if (exch == TRUE)
      {
	str_ptr[9] = '+';
      } else
      {
	str_ptr[9] = '-';
      }
      break;
    case DISSOLVE:
      if (exch == TRUE)
      {
	str_ptr[9] = '-';
      } else
      {
	str_ptr[9] = '+';
      }
      break;
    }
    fprintf(model_file,"%-10s", str_ptr);
    str_ptr = (char *) free_check_null(str_ptr);
/*
 *  Write stoichiometry
 */
    for (next_elt = inv_ptr->phases[i].phase->next_elt; next_elt->elt != NULL; next_elt++)
    {
      f = 1.0;
      if (exch == TRUE) f = -1.0;
      master_ptr = next_elt->elt->primary;
      if (strcmp(master_ptr->elt->name, "Alkalinity") == 0) continue;
      if (strcmp(master_ptr->elt->name, "H") == 0)  continue;
      if (strcmp(master_ptr->elt->name, "O") == 0) continue;
      if (strcmp(master_ptr->elt->name, "E") == 0) continue;
      strcpy(string, master_ptr->elt->name);
      if (strcmp(master_ptr->elt->name, "X") == 0)
      {
	strcpy(string, "Na");
	f  = 1.0;
      }
      str_toupper(string);
      fprintf(model_file," %-2s%12.7f", string, next_elt->coef*f);
    }
/*
 * Calculate RS
 */
    sum = 0;
    for (rxn_ptr = inv_ptr->phases[i].phase->rxn_s->token + 1; rxn_ptr->s != NULL; rxn_ptr++)
    {
      if (rxn_ptr->s == s_hplus) continue;
      if (rxn_ptr->s == s_h2o) continue;
      if (rxn_ptr->s->secondary == NULL && rxn_ptr->s != s_eminus) continue;
      if (rxn_ptr->s == s_o2)
      {
	sum += 4*rxn_ptr->coef;
      } else if (rxn_ptr->s == s_h2)
      {
	sum += -2*rxn_ptr->coef;
      } else if (rxn_ptr->s == s_eminus)
      {
	sum += -1*rxn_ptr->coef;
      } else
      {
	strcpy(string, rxn_ptr->s->secondary->elt->name);
	replace("("," ", string);
	replace(")"," ", string);
	ptr = string;
	copy_token(token, &ptr, &l);
	copy_token(string1, &ptr, &l);
	sscanf(string1,"%lf", &f);
	sum += f*rxn_ptr->coef;
      }
    }
    if (sum != 0.0) fprintf(model_file," %-2s%12.7f", "RS", sum);
/*
 * Add isotopes
 */

    for (k = 0; k < inv_ptr->phases[i].count_isotopes; k++)
    {
      isotope_ptr = inv_ptr->phases[i].isotopes;
      d1 = isotope_ptr[k].ratio;
      for (j = 0; j < inv_ptr->count_isotopes; j++)
      {
	if ((inv_ptr->isotopes[j].elt_name != isotope_ptr[k].elt_name) ||
	    (inv_ptr->isotopes[j].isotope_number != isotope_ptr[k].isotope_number))
	    continue;
	break;
      }
      d2 = 0.0;
      if (j < inv_ptr->count_isotopes) 
      {
	column = col_phase_isotopes + i * inv_ptr->count_isotopes + j;
	if (delta1[col_phases + i] != 0.0)
	{
	  d2 = delta1[column] / delta1[col_phases + i];
	}
      }
      d3 = d1 + d2;
      sprintf(string, "%d%s", (int) isotope_ptr[k].isotope_number, isotope_ptr[k].elt_name);
      if (strcmp(string,"13C") == 0) fprintf(model_file," %-2s%12.7f", "I1", d3);
      if (strcmp(string,"14C") == 0) fprintf(model_file," %-2s%12.7f", "I2", d3);
      if (strcmp(string,"34S") == 0) 
      {
	fprintf(model_file," %-2s%12.7f", "I3", d3);
	fprintf(model_file," %-2s%12.7f", "I7", 0.0);
      }
      if (strcmp(string,"87Sr") == 0) fprintf(model_file," %-2s%12.7f", "I4", d3);
      if (strcmp(string,"15N") == 0) fprintf(model_file," %-2s%12.7f", "I9", d3);
      if (strcmp(string,"2H") == 0) fprintf(model_file," %-2s%12.7f", "D ", d3);
      if (strcmp(string,"3H") == 0) fprintf(model_file," %-2s%12.7f", "TR", d3);
      if (strcmp(string,"18O") == 0) fprintf(model_file," %-2s%12.7f", "18", d3);
    }
/*
 *  Done with stoichiometry
 */
    fprintf(model_file,"\n");
  }
  fprintf(model_file,"\n");
/*
 * Write extra stuff at bottom
 */
  /* 
     (Iflag(i),i=2,6), (P(i),i=1,2), (Isdocrs(i),i=0,5), Disalong, 
     (C14dat(i),i=1,5), (Usera(i),i=1,2), 
     (C14dat(i),i=8,9), i10, i11, (C14dat(i),i=12,13), 
     (Dbdata(Well(0),i),i=44,47), Dbdata(Well(0),49),
     ((Dbdata(Well(iwell),i),i=44,47),Dbdata(Well(iwell),49),Usera(iwell),iwell=1,Iflag(1)+1)
9030 FORMAT (5(I2),2(F8.4),6(I1),F6.3,/,
            7(F8.3),/,
            2(F8.3),2(F8.0), 2(F8.3),/,
            5(F8.3),/,
	    7(6(F8.3),/))
  */
  /* iflags */
  /*fprintf(model_file,"%2d", i);*/  /* not written, 1, mixing, number of mixing wells -1 */
  fprintf(model_file,"%2d", 3);  /* 2, exchange */
  i = 0;
  if (inv_ptr->count_isotopes > 0) i = 1;
  fprintf(model_file,"%2d", i);  /* 3, Rayleigh */
  fprintf(model_file,"%2d", 1);  /* 4, A0 model */
  fprintf(model_file,"%2d", 0);  /* 5, Mook/Deines */
  fprintf(model_file,"%2d", 0);  /* 6, Evaporation/Dilution */
  /* p */
  fprintf(model_file,"%8.4f%8.4f", 1.0, 0.0);  /* P(1),(2) fraction of CO2, fraction of Ca in exch */
  fprintf(model_file,"000000");  /* isdocrs(0,5) doc, rs? */
  fprintf(model_file,"%6.3f\n", 1.0);  /* disalong */
  fprintf(model_file, "   0.000 100.000   0.000   0.000 -25.000 100.000 100.000\n");
  fprintf(model_file, "   0.000   0.000      0.      0.   0.000   0.000  \n");
  /* Final well data */
  fprintf(model_file, " -40.000 -25.000   0.000   0.000   0.000\n");
  /* other wells */
  for (i = 0; i < count_current_solutions - 1; i++)
  {
    fprintf(model_file, " -40.000 -25.000   0.000   0.000   0.000 100.000\n");
  }
  fprintf(model_file, "\n");
  fclose (model_file);

/*
 *    Adjustments to phases
 */
#ifdef SKIP

  if (inv_ptr->count_isotopes > 0)
  {
    /*output_msg (OUTPUT_MESSAGE, "\nIsotopic composition of phases:\n");*/
    for (i = 0; i < inv_ptr->count_phases; i++)
    {
      if (equal (delta1[j], 0.0, toler) == TRUE) continue;


      if (inv_ptr->phases[i].count_isotopes == 0)
	continue;
      j = col_phases + i;
      if (equal (delta1[j], 0.0, toler) == TRUE &&
	  equal (min_delta[j], 0.0, toler) == TRUE &&
	  equal (max_delta[j], 0.0, toler) == TRUE)
	continue;
      isotope_ptr = inv_ptr->phases[i].isotopes;
      for (j = 0; j < inv_ptr->count_isotopes; j++)
      {
	for (k = 0; k < inv_ptr->phases[i].count_isotopes; k++)
	{
	  if (inv_ptr->isotopes[j].elt_name !=
	      isotope_ptr[k].elt_name ||
	      inv_ptr->isotopes[j].isotope_number !=
	      isotope_ptr[k].isotope_number)
	    continue;
	  d1 = isotope_ptr[k].ratio;
	  column = col_phase_isotopes + i * inv_ptr->count_isotopes + j;
	  if (delta1[col_phases + i] != 0.0)
	  {
	    d2 = delta1[column] / delta1[col_phases + i];
	  }
	  else
	  {
	    continue;
	  }
	  d3 = d1 + d2;
	  if (equal (d1, 0.0, 1e-7) == TRUE)
	    d1 = 0.0;
	  if (equal (d2, 0.0, 1e-7) == TRUE)
	    d2 = 0.0;
	  if (equal (d3, 0.0, 1e-7) == TRUE)
	    d3 = 0.0;
	  sprintf (token, "%d%s %s",
		   (int) inv_ptr->isotopes[j].isotope_number,
		   inv_ptr->isotopes[j].elt_name,
		   inv_ptr->phases[i].phase->name);
	  output_msg (OUTPUT_MESSAGE, "%15.15s   %12g  +%12g  =%12g", token,
		      (double) d1, (double) d2, (double) d3);
	  if (fabs (d2) > (isotope_ptr[k].ratio_uncertainty + toler))
	  {
	    output_msg (OUTPUT_MESSAGE, " **");
	    print_msg = TRUE;
	  }
	  output_msg (OUTPUT_MESSAGE, "\n");
	  if (isotope_ptr[k].ratio_uncertainty > 0)
	  {
	    scaled_error += fabs (d2) / isotope_ptr[k].ratio_uncertainty;
/* debug
						output_msg(OUTPUT_MESSAGE, "%e\t%e\t%e\n", fabs(d2) / isotope_ptr[k].ratio_uncertainty, fabs(d2), isotope_ptr[k].ratio_uncertainty);
 */
	  }
	  else if (d2 != 0.0)
	  {
	    error_msg ("Computing delta phase isotope/uncertainty", CONTINUE);
	  }
	}
      }
    }
  }
#endif
#ifdef SKIP
  if (print_msg == TRUE)
  {
    output_msg (OUTPUT_MESSAGE,
		"\n**\tWARNING: The adjustment to at least one isotopic"
		"\n\tcomposition of a phase exceeded the specified uncertainty"
		"\n\tfor the phase.  If the phase is not constrained to dissolve"
		"\n\tor precipitate, then the isotopic composition of the phase"
		"\n\tis also unconstrained.\n");
  }
  output_msg (OUTPUT_MESSAGE, "\n%-20.20s   %7s   %12.12s   %12.12s\n",
	      "Solution fractions:", " ", "Minimum", "Maximum");
  for (i = 0; i < inv_ptr->count_solns; i++)
  {
    d1 = delta1[i];
    d2 = min_delta[i];
    d3 = max_delta[i];
    if (equal (d1, 0.0, MIN_TOTAL_INVERSE) == TRUE)
      d1 = 0.0;
    if (equal (d2, 0.0, MIN_TOTAL_INVERSE) == TRUE)
      d2 = 0.0;
    if (equal (d3, 0.0, MIN_TOTAL_INVERSE) == TRUE)
      d3 = 0.0;
    output_msg (OUTPUT_MESSAGE, "%11s%4d   %12.3e   %12.3e   %12.3e\n",
		"Solution", inv_ptr->solns[i], (double) d1, (double) d2,
		(double) d3);
  }

  output_msg (OUTPUT_MESSAGE, "\n%-25.25s   %2s   %12.12s   %12.12s\n",
	      "Phase mole transfers:", " ", "Minimum", "Maximum");
  for (i = col_phases; i < col_redox; i++)
  {
    if (equal (delta1[i], 0.0, toler) == TRUE &&
	equal (min_delta[i], 0.0, toler) == TRUE &&
	equal (max_delta[i], 0.0, toler) == TRUE)
      continue;
    d1 = delta1[i];
    d2 = min_delta[i];
    d3 = max_delta[i];
    if (equal (d1, 0.0, MIN_TOTAL_INVERSE) == TRUE)
      d1 = 0.0;
    if (equal (d2, 0.0, MIN_TOTAL_INVERSE) == TRUE)
      d2 = 0.0;
    if (equal (d3, 0.0, MIN_TOTAL_INVERSE) == TRUE)
      d3 = 0.0;
    output_msg (OUTPUT_MESSAGE, "%15.15s   %12.3e   %12.3e   %12.3e   %s\n",
		col_name[i], (double) d1, (double) d2, (double) d3,
		inv_ptr->phases[i - col_phases].phase->formula);
  }

  output_msg (OUTPUT_MESSAGE, "\n%-25.25s\n", "Redox mole transfers:");
  for (i = col_redox; i < col_epsilon; i++)
  {
    if (equal (delta1[i], 0.0, toler) == TRUE)
      continue;
    output_msg (OUTPUT_MESSAGE, "%15.15s   %12.3e\n", col_name[i],
		(double) delta1[i]);
  }

  output_msg (OUTPUT_MESSAGE,
	      "\nSum of residuals (epsilons in documentation):      %12.3e\n",
	      ((double) error / SCALE_EPSILON));
  output_msg (OUTPUT_MESSAGE,
	      "Sum of delta/uncertainty limit:                    %12.3e\n",
	      (double) scaled_error);
  output_msg (OUTPUT_MESSAGE,
	      "Maximum fractional error in element concentration: %12.3e\n",
	      (double) max_pct);
/*
 *   Flush buffer after each model
 */
  output_fflush (OUTPUT_MESSAGE);
#endif
  return (OK);
}
/* ---------------------------------------------------------------------- */
void
print_total_pat (FILE *netpath_file, char *elt, char *string)
/* ---------------------------------------------------------------------- */
{
  double d;
  d  = 1000.0 * total(elt);
  if (d == 0) {
    fprintf(netpath_file, "%14g%1s    # %s\n", d, "*", string);
  } else {
    fprintf(netpath_file, "%14g%1s    # %s\n", d, " ", string);
  }
}
/* ---------------------------------------------------------------------- */
int
set_initial_solution (int n_user_old, int n_user_new)
/* ---------------------------------------------------------------------- */
{
  int j;
  struct solution *solution_ptr;
  struct conc *conc_ptr;

  solution_duplicate (n_user_old, n_user_new);
  solution_ptr = solution_bsearch (n_user_new, &j, TRUE);
  solution_ptr->new_def = TRUE;
  solution_ptr->n_user_end = n_user_new;
  for (j = 0; solution_ptr->totals[j].description != NULL; j++)
  {
    conc_ptr = &solution_ptr->totals[j];
    conc_ptr->input_conc = conc_ptr->moles / solution_ptr->mass_water;
    conc_ptr->units = string_hsave ("Mol/kgw");
  }
  return (OK);
}
/* ---------------------------------------------------------------------- */
int
add_to_file(char *filename, char *string)
/* ---------------------------------------------------------------------- */
{
  FILE *model_file;
  char c;
  int i; 
  char string_line[MAX_LENGTH];

  model_file = fopen(filename, "r");
  if (model_file == NULL) {
    model_file = fopen(filename, "w");
    if (model_file == NULL)
    {
      sprintf (error_string, "Can't open file, %s.", filename);
      error_msg (error_string, STOP);
    }
  }
  i = 0;
  for (;;)
  {
    c = getc(model_file);
    if (c != EOF && c != '\n')
    {
      string_line[i] = c;
      i++;
      continue;
    }
    string_line[i] = '\0';
    /* new line or eof */
    string_trim(string_line);
    if (strcmp(string_line, string) == 0)
    {
      fclose(model_file);
      return(OK);
    }
    /* eof, add line */
    if (c == EOF) 
    {
      fclose(model_file);
      model_file = fopen(filename, "a");
      fprintf(model_file, "%s\n", string);
      fclose(model_file);
      return(OK);
    }
    /* new line */
    i = 0;
  }
}
