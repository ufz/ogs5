#define EXTERNAL extern
#include "global.h"
#include "phqalloc.h"
#include "output.h"
#include "phrqproto.h"

static char const svnid[] =
  "$Id: structures.c 4 2009-04-21 17:29:29Z delucia $";

static int exchange_compare_int (const void *ptr1, const void *ptr2);
static int gas_phase_compare_int (const void *ptr1, const void *ptr2);

static int inverse_compare (const void *ptr1, const void *ptr2);
static int inverse_free (struct inverse *inverse_ptr);

static int irrev_compare (const void *ptr1, const void *ptr2);
static int irrev_compare_int (const void *ptr1, const void *ptr2);

static int kinetics_compare_int (const void *ptr1, const void *ptr2);

static int logk_init (struct logk *logk_ptr);

static int master_compare_string (const void *ptr1, const void *ptr2);
static int master_free (struct master *master_ptr);

#if defined(PHREEQCI_GUI)
int mix_compare (const void *ptr1, const void *ptr2);
#else
static int mix_compare (const void *ptr1, const void *ptr2);
#endif
static int mix_compare_int (const void *ptr1, const void *ptr2);

static struct phase *phase_alloc (void);
static int phase_compare_string (const void *ptr1, const void *ptr2);
static int phase_free (struct phase *phase_ptr);
static int phase_init (struct phase *phase_ptr);

static int pp_assemblage_compare_int (const void *ptr1, const void *ptr2);

static int rate_compare (const void *ptr1, const void *ptr2);
static int rate_compare_string (const void *ptr1, const void *ptr2);

static struct species *s_alloc (void);
static int s_free (struct species *s_ptr);
static int s_init (struct species *s_ptr);

static int s_s_assemblage_compare_int (const void *ptr1, const void *ptr2);

static int solution_compare (const void *ptr1, const void *ptr2);
static int solution_compare_int (const void *ptr1, const void *ptr2);

static int species_list_compare (const void *ptr1, const void *ptr2);

static int surface_compare_int (const void *ptr1, const void *ptr2);

static int temperature_compare (const void *ptr1, const void *ptr2);
static int temperature_compare_int (const void *ptr1, const void *ptr2);

static int rxn_token_temp_compare (const void *ptr1, const void *ptr2);
static int trxn_multiply (LDBLE coef);

#ifdef PHREEQCI_GUI
extern void free_spread (void);
#endif
#if defined(USE_MPI) && defined(HDF5_CREATE) && defined(MERGE_FILES)
extern void MergeFinalize (void);
#endif
extern LDBLE *scratch, *x_arg, *res_arg;
extern LDBLE *normal, *ineq_array, *zero, *res, *delta1, *cu;
extern int *iu, *is, *back_eq;
extern int x_arg_max, res_arg_max, scratch_max;

/* ---------------------------------------------------------------------- */
int
clean_up (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Free all allocated memory, except strings
 */
  int i, j;
  if (svnid == NULL)
    fprintf (stderr, " ");

  description_x = (char *) free_check_null (description_x);
  isotopes_x = (struct isotope *) free_check_null (isotopes_x);
  moles_per_kilogram_string =
    (char *) free_check_null (moles_per_kilogram_string);
  pe_string = (char *) free_check_null (pe_string);
/* model */
  last_model.exchange =
    (struct master **) free_check_null (last_model.exchange);
  last_model.gas_phase =
    (struct phase **) free_check_null (last_model.gas_phase);
  last_model.pp_assemblage =
    (struct phase **) free_check_null (last_model.pp_assemblage);
  last_model.s_s_assemblage =
    (char **) free_check_null (last_model.s_s_assemblage);
  last_model.add_formula = (char **) free_check_null (last_model.add_formula);
  last_model.si = (double *) free_check_null (last_model.si);
  last_model.surface_comp =
    (char **) free_check_null (last_model.surface_comp);
  last_model.surface_charge =
    (char **) free_check_null (last_model.surface_charge);

  /* model */
  free_model_allocs ();

/* species */

  for (j = 0; j < count_s; j++)
  {
    s_free (s[j]);
    s[j] = (struct species *) free_check_null (s[j]);
  }
  s = (struct species **) free_check_null (s);

/* master species */

  for (j = 0; j < count_master; j++)
  {
    master_free (master[j]);
  }
  master = (struct master **) free_check_null (master);

/* elements */

  for (j = 0; j < count_elements; j++)
  {
    elements[j] = (struct element *) free_check_null (elements[j]);
  }
  elements = (struct element **) free_check_null (elements);

/* solutions */

  for (j = 0; j < count_solution; j++)
  {
    solution_free (solution[j]);
  }
  solution = (struct solution **) free_check_null (solution);

/* surfaces */

  for (j = 0; j < count_surface; j++)
  {
    surface_free (&surface[j]);
  }
  surface = (struct surface *) free_check_null (surface);

/* exchange */

  for (j = 0; j < count_exchange; j++)
  {
    exchange_free (&exchange[j]);
  }
  exchange = (struct exchange *) free_check_null (exchange);

/* pp assemblages */

  for (j = 0; j < count_pp_assemblage; j++)
  {
    pp_assemblage_free (&pp_assemblage[j]);
  }
  pp_assemblage = (struct pp_assemblage *) free_check_null (pp_assemblage);

/* s_s assemblages */

  for (j = 0; j < count_s_s_assemblage; j++)
  {
    s_s_assemblage_free (&s_s_assemblage[j]);
  }
  s_s_assemblage = (struct s_s_assemblage *) free_check_null (s_s_assemblage);

/* irreversible reactions */

  for (j = 0; j < count_irrev; j++)
  {
    irrev_free (&irrev[j]);
  }
  irrev = (struct irrev *) free_check_null (irrev);

/* temperature */

  for (j = 0; j < count_temperature; j++)
  {
    temperature_free (&temperature[j]);
  }
  temperature = (struct temperature *) free_check_null (temperature);

/* unknowns */

  for (j = 0; j < max_unknowns; j++)
  {
    unknown_free (x[j]);
  }
  x = (struct unknown **) free_check_null (x);

/* mixtures */

  for (j = 0; j < count_mix; j++)
  {
    mix_free (&mix[j]);
  }
  mix = (struct mix *) free_check_null (mix);

/* phases */

  for (j = 0; j < count_phases; j++)
  {
    phase_free (phases[j]);
    phases[j] = (struct phase *) free_check_null (phases[j]);
  }
  phases = (struct phase **) free_check_null (phases);

/* inverse */
  for (j = 0; j < count_inverse; j++)
  {
    inverse_free (&(inverse[j]));
  }
  inverse = (struct inverse *) free_check_null (inverse);

/* gases */
  for (j = 0; j < count_gas_phase; j++)
  {
    gas_phase_free (&gas_phase[j]);
  };
  gas_phase = (struct gas_phase *) free_check_null (gas_phase);

/* kinetics */
  for (j = 0; j < count_kinetics; j++)
  {
    kinetics_free (&kinetics[j]);
  }
  kinetics = (struct kinetics *) free_check_null (kinetics);

/* rates */
  for (j = 0; j < count_rates; j++)
  {
    rate_free (&rates[j]);
  }
  rates = (struct rate *) free_check_null (rates);

/* logk hash table */
  for (j = 0; j < count_logk; j++)
  {
    free_check_null (logk[j]->add_logk);
    logk[j] = (struct logk *) free_check_null (logk[j]);
  }
  logk = (struct logk **) free_check_null (logk);

/* save_values */
  for (j = 0; j < count_save_values; j++)
  {
    save_values[j].subscripts =
      (int *) free_check_null (save_values[j].subscripts);
  }
  save_values = (struct save_values *) free_check_null (save_values);

/* model */

/* global solution */

  pe_data_free (pe_x);
  pe_x = NULL;

/* species_list */

  species_list = (struct species_list *) free_check_null (species_list);

/* transport data */

  stag_data = (struct stag_data *) free_check_null (stag_data);
  cell_data = (struct cell_data *) free_check_null (cell_data);

/* punch */

  punch.totals = (struct name_master *) free_check_null (punch.totals);
  punch.molalities =
    (struct name_species *) free_check_null (punch.molalities);
  punch.activities =
    (struct name_species *) free_check_null (punch.activities);
  punch.pure_phases =
    (struct name_phase *) free_check_null (punch.pure_phases);
  punch.si = (struct name_phase *) free_check_null (punch.si);
  punch.gases = (struct name_phase *) free_check_null (punch.gases);
  punch.s_s = (struct name_phase *) free_check_null (punch.s_s);
  punch.kinetics = (struct name_phase *) free_check_null (punch.kinetics);
  advection_punch = (int *) free_check_null (advection_punch);
  advection_print = (int *) free_check_null (advection_print);
  punch.isotopes = (struct name_master *) free_check_null (punch.isotopes);
  punch.calculate_values =
    (struct name_master *) free_check_null (punch.calculate_values);

/*  user_print and user_punch */
  rate_free (user_print);
  rate_free (user_punch);
  user_print = (struct rate *) free_check_null (user_print);

  user_punch = (struct rate *) free_check_null (user_punch);
  user_punch_headings = (char **) free_check_null (user_punch_headings);

  /*
     Free llnl aqueous model parameters
   */
  llnl_temp = (LDBLE *) free_check_null (llnl_temp);
  llnl_adh = (LDBLE *) free_check_null (llnl_adh);
  llnl_bdh = (LDBLE *) free_check_null (llnl_bdh);
  llnl_bdot = (LDBLE *) free_check_null (llnl_bdot);
  llnl_co2_coefs = (LDBLE *) free_check_null (llnl_co2_coefs);
  /*
   * Copier space
   */
  copier_free (&copy_solution);
  copier_free (&copy_pp_assemblage);
  copier_free (&copy_exchange);
  copier_free (&copy_surface);
  copier_free (&copy_s_s_assemblage);
  copier_free (&copy_gas_phase);
  copier_free (&copy_kinetics);
  copier_free (&copy_mix);
  copier_free (&copy_irrev);
  copier_free (&copy_temperature);

#ifdef PHREEQ98
  rate_free (user_graph);
  user_graph = free_check_null (user_graph);
  user_graph_headings = free_check_null (user_graph_headings);
#endif
  /* master_isotope */
  for (i = 0; i < count_master_isotope; i++)
  {
    master_isotope[i] =
      (struct master_isotope *) free_check_null (master_isotope[i]);
  }
  master_isotope =
    (struct master_isotope **) free_check_null (master_isotope);
  hdestroy_multi (master_isotope_hash_table);
  master_isotope_hash_table = NULL;

  /* calculate_value */
  for (i = 0; i < count_calculate_value; i++)
  {
    calculate_value_free (calculate_value[i]);
    calculate_value[i] =
      (struct calculate_value *) free_check_null (calculate_value[i]);
  }
  calculate_value =
    (struct calculate_value **) free_check_null (calculate_value);
  hdestroy_multi (calculate_value_hash_table);
  calculate_value_hash_table = NULL;

  /* isotope_ratio */
  for (i = 0; i < count_isotope_ratio; i++)
  {
    isotope_ratio[i] =
      (struct isotope_ratio *) free_check_null (isotope_ratio[i]);
  }
  isotope_ratio = (struct isotope_ratio **) free_check_null (isotope_ratio);
  hdestroy_multi (isotope_ratio_hash_table);
  isotope_ratio_hash_table = NULL;

  /* isotope_alpha */
  for (i = 0; i < count_isotope_alpha; i++)
  {
    isotope_alpha[i] =
      (struct isotope_alpha *) free_check_null (isotope_alpha[i]);
  }
  isotope_alpha = (struct isotope_alpha **) free_check_null (isotope_alpha);
  hdestroy_multi (isotope_alpha_hash_table);
  isotope_alpha_hash_table = NULL;

  free_tally_table ();

  /* CVODE memory */
  free_cvode ();

  pitzer_clean_up ();


/* hash tables */

  hdestroy_multi (elements_hash_table);
  hdestroy_multi (species_hash_table);
  hdestroy_multi (logk_hash_table);
  hdestroy_multi (phases_hash_table);
  hdestroy_multi (keyword_hash_table);
  keyword_hash = (struct key *) free_check_null (keyword_hash);

  elements_hash_table = NULL;
  species_hash_table = NULL;
  logk_hash_table = NULL;
  phases_hash_table = NULL;
  keyword_hash_table = NULL;

/* strings */
  free_hash_strings (strings_hash_table);
  hdestroy_multi (strings_hash_table);
  strings_hash_table = NULL;

/* basic commands hash table */
  cmd_free ();
  change_surf = (struct Change_Surf *) free_check_null (change_surf);

/* miscellaneous work space */

  elt_list = (struct elt_list *) free_check_null (elt_list);
  trxn.token = (struct rxn_token_temp *) free_check_null (trxn.token);
  mb_unknowns = (struct unknown_list *) free_check_null (mb_unknowns);

  // MDL: this was the original
  line = (char *) free_check_null (line);
  line_save = (char *) free_check_null (line_save);

  zeros = (LDBLE *) free_check_null (zeros);
  scratch = (LDBLE *) free_check_null (scratch);
  x_arg = (LDBLE *) free_check_null (x_arg);
  res_arg = (LDBLE *) free_check_null (res_arg);

  normal = (LDBLE *) free_check_null (normal);
  ineq_array = (LDBLE *) free_check_null (ineq_array);
  back_eq = (int *) free_check_null (back_eq);
  zero = (LDBLE *) free_check_null (zero);
  res = (LDBLE *) free_check_null (res);
  delta1 = (LDBLE *) free_check_null (delta1);
  cu = (LDBLE *) free_check_null (cu);
  iu = (int *) free_check_null (iu);
  is = (int *) free_check_null (is);

  //x_arg = res_arg = scratch = NULL;
  x_arg_max = res_arg_max = scratch_max = 0;


/* free user database name if defined */
  user_database = (char *) free_check_null (user_database);
  selected_output_file_name =
    (char *) free_check_null (selected_output_file_name);
  dump_file_name = (char *) free_check_null (dump_file_name);
#ifdef PHREEQCI_GUI
  free_spread ();
  title_x = free_check_null (title_x);
#endif
#if defined(USE_MPI) && defined(HDF5_CREATE) && defined(MERGE_FILES)
  MergeFinalize ();
#endif

  PHRQ_free_all ();

  clean_up_output_callbacks ();

  count_solution = 0;
  count_pp_assemblage = 0;
  count_exchange = 0;
  count_surface = 0;
  count_gas_phase = 0;
  count_kinetics = 0;
  count_s_s_assemblage = 0;

  count_elements = 0;
  count_irrev = 0;
  count_master = 0;
  count_mix = 0;
  count_phases = 0;
  count_s = 0;
  count_temperature = 0;
  count_logk = 0;
  count_master_isotope = 0;

  count_rates = 0;
  count_inverse = 0;
  count_save_values = 0;

  llnl_count_temp = 0;
  llnl_count_adh = 0;
  llnl_count_bdh = 0;
  llnl_count_bdot = 0;
  llnl_count_co2_coefs = 0;

  count_calculate_value = 0;
  count_isotope_ratio = 0;
  count_isotope_alpha = 0;

  return (OK);
}

/* ---------------------------------------------------------------------- */
int
reinitialize (void)
/* ---------------------------------------------------------------------- */
{
  int j;
/* solutions */
  for (j = 0; j < count_solution; j++)
  {
    solution_free (solution[j]);
  }
  count_solution = 0;

/* surfaces */

  for (j = 0; j < count_surface; j++)
  {
    surface_free (&surface[j]);
  }
  count_surface = 0;

/* exchange */

  for (j = 0; j < count_exchange; j++)
  {
    exchange_free (&exchange[j]);
  }
  count_exchange = 0;

/* pp assemblages */

  for (j = 0; j < count_pp_assemblage; j++)
  {
    pp_assemblage_free (&pp_assemblage[j]);
  }
  count_pp_assemblage = 0;

/* s_s assemblages */

  for (j = 0; j < count_s_s_assemblage; j++)
  {
    s_s_assemblage_free (&s_s_assemblage[j]);
  }
  count_s_s_assemblage = 0;

/* gases */
  for (j = 0; j < count_gas_phase; j++)
  {
    gas_phase_free (&gas_phase[j]);
  };
  count_gas_phase = 0;

/* kinetics */
  for (j = 0; j < count_kinetics; j++)
  {
    kinetics_free (&kinetics[j]);
  }
  count_kinetics = 0;

/* irreversible reactions */

  for (j = 0; j < count_irrev; j++)
  {
    irrev_free (&irrev[j]);
  }
  count_irrev = 0;

/* temperature */

  for (j = 0; j < count_temperature; j++)
  {
    temperature_free (&temperature[j]);
  }
  count_temperature = 0;
  return (OK);
}

/* **********************************************************************
 *
 *   Routines related to structure "element"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
int
element_compare (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
  const struct element *element_ptr1, *element_ptr2;
  element_ptr1 = *(const struct element **) ptr1;
  element_ptr2 = *(const struct element **) ptr2;
/*      return(strcmp_nocase(element_ptr1->name, element_ptr2->name)); */
  return (strcmp (element_ptr1->name, element_ptr2->name));

}

/* ---------------------------------------------------------------------- */
struct element *
element_store ( char *element)
/* ---------------------------------------------------------------------- */
{
/*
 *   Function locates the string "element" in the hash table for elements.
 *
 *   If found, pointer to the appropriate element structure is returned.
 *
 *   If the string is not found, a new entry is made at the end of
 *   the elements array (position count_elements) and count_elements is
 *   incremented. A new entry is made in the hash table. Pointer to
 *   the new structure is returned.
 *
 *   Arguments:
 *      element    input, character string to be located or stored.
 *
 *   Returns:
 *      The address of an elt structure that contains the element data.
 */
  int n;
  struct element *elts_ptr;
  ENTRY item, *found_item;
/*
 *   Search list
 */
  item.key = element;
  item.data = NULL;
  found_item = hsearch_multi (elements_hash_table, item, FIND);
  if (found_item != NULL)
  {
    elts_ptr = (struct element *) (found_item->data);
    return (elts_ptr);
  }
/*
 *   Save new elt structure and return pointer to it
 */
  /* make sure there is space in elements */
  elements[count_elements] =
    (struct element *) PHRQ_malloc ((size_t) sizeof (struct element));
  if (elements[count_elements] == NULL)
    malloc_error ();
  /* set name pointer in elements structure */
  elements[count_elements]->name = string_hsave (element);
  /* set return value */
  elements[count_elements]->master = NULL;
  elements[count_elements]->primary = NULL;
  elements[count_elements]->gfw = 0.0;
  n = count_elements++;
  if (count_elements >= max_elements)
  {
    space ((void **) ((void *) &elements), count_elements, &max_elements,
	   sizeof (struct element *));
  }
/*
 *   Update hash table
 */
  item.key = elements[n]->name;
  item.data = (void *) elements[n];
  found_item = hsearch_multi (elements_hash_table, item, ENTER);
  if (found_item == NULL)
  {
    sprintf (error_string, "Hash table error in element_store.");
    error_msg (error_string, CONTINUE);
  }
  return (elements[n]);
}
/* **********************************************************************
 *
 *   Routines related to structure "elt_list"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
int
elt_list_combine (void)
/* ---------------------------------------------------------------------- */
/*
 *      Function goes through the list of elements pointed to by elt_list
 *      and combines the coefficients of elements that are the same.
 *      Assumes elt_list has been sorted by element name.
 */
{
  int i, j;

  if (count_elts < 1)
  {
    output_msg (OUTPUT_MESSAGE, "elt_list_combine: How did this happen?\n");
    return (ERROR);
  }
  if (count_elts == 1)
  {
    return (OK);
  }
  j = 0;
  for (i = 1; i < count_elts; i++)
  {
    if (elt_list[i].elt == elt_list[j].elt)
    {
      elt_list[j].coef += elt_list[i].coef;
    }
    else
    {
      j++;
      if (i != j)
      {
	elt_list[j].elt = elt_list[i].elt;
	elt_list[j].coef = elt_list[i].coef;
      }
    }
  }
  count_elts = j + 1;
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
elt_list_compare (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
  const struct elt_list *a, *b;

  a = (const struct elt_list *) ptr1;
  b = (const struct elt_list *) ptr2;
  return (strncmp (a->elt->name, b->elt->name, MAX_LENGTH));
}

/* ---------------------------------------------------------------------- */
struct elt_list *
elt_list_dup (struct elt_list *elt_list_ptr_old)
/* ---------------------------------------------------------------------- */
{
/*
 *  Duplicates the elt_list structure pointed to by elt_list_ptr_old.
 */
  int i, count_totals;
  struct elt_list *elt_list_ptr_new;
/*
 *   Count totals data and copy
 */
  if (elt_list_ptr_old == NULL)
    return (NULL);
  for (i = 0; elt_list_ptr_old[i].elt != NULL; i++);
  count_totals = i;
/*
 *   Malloc space and store element data
 */
  elt_list_ptr_new =
    (struct elt_list *) PHRQ_malloc ((size_t) (count_totals + 1) *
				     sizeof (struct elt_list));
  if (elt_list_ptr_new == NULL)
    malloc_error ();
  memcpy (elt_list_ptr_new, elt_list_ptr_old,
	  (size_t) (count_totals + 1) * sizeof (struct elt_list));
  return (elt_list_ptr_new);
}

/* ---------------------------------------------------------------------- */
int
elt_list_print (struct elt_list *elt_list_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *  Duplicates the elt_list structure pointed to by elt_list_ptr_old.
 */
  int i;
/*
 *   Debug print for element list
 */
  if (elt_list_ptr == NULL)
    return (ERROR);
  output_msg (OUTPUT_MESSAGE, "Elt_list\n");
  for (i = 0; elt_list_ptr[i].elt != NULL; i++)
  {
    output_msg (OUTPUT_MESSAGE, "\t%s\t%e\n", elt_list_ptr[i].elt->name,
		(double) elt_list_ptr[i].coef);
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
elt_list_print0 (void)
/* ---------------------------------------------------------------------- */
{
/*
 *  Duplicates the elt_list structure pointed to by elt_list_ptr_old.
 */
  int i;
/*
 *   Debug print for element list
 */
  output_msg (OUTPUT_STDERR, "Elt_list\n");
  for (i = 0; i < count_elts; i++)
  {
    output_msg (OUTPUT_STDERR, "\t%s\t%e\n", elt_list[i].elt->name,
		(double) elt_list[i].coef);
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
struct elt_list *
elt_list_save (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Takes data from work space elt_list, allocates a new elt_list structure,
 *   copies data from work space to new structure, and returns pointer to
 *   new structure.
 */
  int j;
  struct elt_list *elt_list_ptr;
/*
 *   Sort elements in reaction and combine
 */
  if (count_elts > 0)
  {
    qsort (elt_list, (size_t) count_elts,
	   (size_t) sizeof (struct elt_list), elt_list_compare);
    elt_list_combine ();
  }
/*
 *   Malloc space and store element data
 */
  elt_list_ptr =
    (struct elt_list *) PHRQ_malloc ((size_t) (count_elts + 1) *
				     sizeof (struct elt_list));
  if (elt_list_ptr == NULL)
    malloc_error ();
  for (j = 0; j < count_elts; j++)
  {
    elt_list_ptr[j].elt = elt_list[j].elt;
    elt_list_ptr[j].coef = elt_list[j].coef;
  }
  elt_list_ptr[count_elts].elt = NULL;
  return (elt_list_ptr);
}

/* **********************************************************************
 *
 *   Routines related to structure "exchange"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
struct exchange *
exchange_alloc (void)
/* ---------------------------------------------------------------------- */
{
  struct exchange *exchange_ptr;
  exchange_ptr = (struct exchange *) PHRQ_malloc (sizeof (struct exchange));
  if (exchange_ptr == NULL)
    malloc_error ();

  exchange_ptr->n_user = -1;
  exchange_ptr->n_user_end = -1;
  exchange_ptr->new_def = 0;
  exchange_ptr->description = NULL;
  exchange_ptr->solution_equilibria = 0;
  exchange_ptr->n_solution = 0;
  exchange_ptr->count_comps = 0;
  exchange_ptr->comps = NULL;
  exchange_ptr->related_phases = 0;
  exchange_ptr->related_rate = 0;
  exchange_ptr->pitzer_exchange_gammas = 0;

  return (exchange_ptr);
}

/* ---------------------------------------------------------------------- */
struct exchange *
exchange_bsearch (int k, int *n)
/* ---------------------------------------------------------------------- */
{
/*
 *   Binary search exchange to find if user exchange number k exists.
 *   Exchange is assumed to be in sort order by user number.
 */
  void *void_ptr;
  if (count_exchange == 0)
  {
    *n = -999;
    return (NULL);
  }
  void_ptr = (void *)
    bsearch ((char *) &k,
	     (char *) exchange,
	     (size_t) count_exchange,
	     (size_t) sizeof (struct exchange), exchange_compare_int);
  if (void_ptr == NULL)
  {
    *n = -999;
    return (NULL);
  }
  *n = (int) ((struct exchange *) void_ptr - exchange);
  return ((struct exchange *) void_ptr);
}

/* ---------------------------------------------------------------------- */
int
exchange_comp_compare (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
  const struct exch_comp *exch_comp_ptr1, *exch_comp_ptr2;
  exch_comp_ptr1 = (const struct exch_comp *) ptr1;
  exch_comp_ptr2 = (const struct exch_comp *) ptr2;
  return (strcmp_nocase
	  (exch_comp_ptr1->master->elt->name,
	   exch_comp_ptr2->master->elt->name));
}

/* ---------------------------------------------------------------------- */
int
exchange_compare (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
  const struct exchange *exchange_ptr1, *exchange_ptr2;
  exchange_ptr1 = (const struct exchange *) ptr1;
  exchange_ptr2 = (const struct exchange *) ptr2;
  if (exchange_ptr1->n_user > exchange_ptr2->n_user)
    return (1);
  if (exchange_ptr1->n_user < exchange_ptr2->n_user)
    return (-1);
  return (0);

}

/* ---------------------------------------------------------------------- */
static int
exchange_compare_int (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
/*
 *   Compare exchange user numbers
 */
  const int *nptr1;
  const struct exchange *nptr2;

  nptr1 = (const int *) ptr1;
  nptr2 = (const struct exchange *) ptr2;
  if (*nptr1 > nptr2->n_user)
    return (1);
  if (*nptr1 < nptr2->n_user)
    return (-1);
  return (0);
}

/* ---------------------------------------------------------------------- */
int
exchange_copy (struct exchange *exchange_old_ptr,
	       struct exchange *exchange_new_ptr, int n_user_new)
/* ---------------------------------------------------------------------- */
{
/*
 *   Copies exchange data from exchange_old_ptr to new location, exchange_new_ptr.
 *   Space for the exchange_new_ptr structure must already be malloced.
 *   Space for exchange components is malloced here.
 */
  int count_comps, i;
  char token[MAX_LENGTH];
/*
 *   Store data for structure exchange
 */
  memcpy (exchange_new_ptr, exchange_old_ptr, sizeof (struct exchange));
  exchange_new_ptr->n_user = n_user_new;
  exchange_new_ptr->n_user_end = n_user_new;
  count_comps = exchange_old_ptr->count_comps;
  exchange_new_ptr->new_def = exchange_old_ptr->new_def;
  exchange_new_ptr->solution_equilibria = FALSE;
  exchange_new_ptr->n_solution = -2;
  exchange_new_ptr->count_comps = count_comps;
  sprintf (token, "Exchanger defined in simulation %d.", simulation);
  exchange_new_ptr->description = string_duplicate (token);
/*
 *   Count exchange components and allocate space
 */
  exchange_new_ptr->comps =
    (struct exch_comp *) PHRQ_malloc ((size_t) (count_comps) *
				      sizeof (struct exch_comp));
  if (exchange_new_ptr->comps == NULL)
    malloc_error ();
/*
 *   Write exch_comp structure for each exchange component
 */
  memcpy (exchange_new_ptr->comps, exchange_old_ptr->comps,
	  (size_t) (count_comps) * sizeof (struct exch_comp));
  for (i = 0; i < count_comps; i++)
  {
    exchange_new_ptr->comps[i].totals =
      elt_list_dup (exchange_old_ptr->comps[i].totals);
    exchange_new_ptr->comps[i].formula_totals =
      elt_list_dup (exchange_old_ptr->comps[i].formula_totals);
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
exchange_copy_to_last (int n, int n_user)
/* ---------------------------------------------------------------------- */
{
/*
 *   Copies an exchange definition from position n to position count_exchange.
 *   New exchange structure is given user number n_user.
 */
  space ((void **) ((void *) &exchange), count_exchange, &max_exchange,
	 sizeof (struct exchange));
  exchange_copy (&exchange[n], &exchange[count_exchange], n_user);
  count_exchange++;
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
exchange_delete (int n_user_old)
/* ---------------------------------------------------------------------- */
/*
 *   Frees space for user number n_user_old, removes structure from
 *   array exchange.
 */
{
  int i;
  int n_old;
  struct exchange *exchange_ptr_old;
/*
 *   Find n_user_old in structure array
 */
  exchange_ptr_old = exchange_bsearch (n_user_old, &n_old);
  if (exchange_ptr_old != NULL)
  {
    /*
     *   Delete exchange
     */
    exchange_free (&exchange[n_old]);

    for (i = n_old + 1; i < count_exchange; i++)
    {
      memcpy ((void *) &exchange[i - 1], (void *) &exchange[i],
	      (size_t) sizeof (struct exchange));
    }
    count_exchange--;
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
exchange_duplicate (int n_user_old, int n_user_new)
/* ---------------------------------------------------------------------- */
{
/*
 *   Checks if n_user_new exists, and frees space if it does
 *   Copies exchange[n_user_old] to old n_user_new space if
 *   found or to exchange[count_exchange] if not found.
 *   Exchange array may not be in sort order after the copy.
 */
  int sort;
  int n_old, n_new;
  struct exchange *exchange_ptr_old, *exchange_ptr_new;
/*
 *   Find n_user_old in structure array exchange
 */
  if (n_user_old == n_user_new)
    return (OK);
  exchange_ptr_old = exchange_bsearch (n_user_old, &n_old);
  if (exchange_ptr_old == NULL)
  {
    sprintf (error_string, "Exchange %d not found.", n_user_old);
    error_msg (error_string, CONTINUE);
    input_error++;
    return (ERROR);
  }
/*
 *   Find n_user_new in structure array exchange or make new space
 */
  sort = FALSE;
  exchange_ptr_new = exchange_bsearch (n_user_new, &n_new);
  if (exchange_ptr_new != NULL)
  {
    exchange_free (exchange_ptr_new);
  }
  else
  {
    space ((void **) ((void *) &exchange), count_exchange, &max_exchange,
	   sizeof (struct exchange));
    if (n_user_new < exchange[count_exchange - 1].n_user)
      sort = TRUE;
    n_new = count_exchange++;
  }
/*
 *   Copy data
 */
  exchange_copy (&exchange[n_old], &exchange[n_new], n_user_new);
  if (sort == TRUE)
    exchange_sort ();
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
exchange_free (struct exchange *exchange_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Frees all data associated with exchange structure.
 */
  int i;
  if (exchange_ptr == NULL)
    return (ERROR);
/*
 *   Free space allocated for exchange structure
 */
  exchange_ptr->description =
    (char *) free_check_null (exchange_ptr->description);
  for (i = 0; i < exchange_ptr->count_comps; i++)
  {
    exchange_ptr->comps[i].totals =
      (struct elt_list *) free_check_null (exchange_ptr->comps[i].totals);
    exchange_ptr->comps[i].formula_totals =
      (struct elt_list *) free_check_null (exchange_ptr->comps[i].
					   formula_totals);
  }
  exchange_ptr->comps =
    (struct exch_comp *) free_check_null (exchange_ptr->comps);
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
exchange_init (struct exchange *exchange_ptr, int n_user, int n_user_end,
	       char *description)
/* ---------------------------------------------------------------------- */
{
/*
 *   Frees all data associated with exchange structure.
 */

  if (exchange_ptr == NULL)
    return (ERROR);
  exchange_ptr->n_user = n_user;
  exchange_ptr->n_user_end = n_user_end;
  exchange_ptr->new_def = TRUE;
  exchange_ptr->description = string_duplicate (description);
  exchange_ptr->solution_equilibria = FALSE;
  exchange_ptr->n_solution = 0;
  exchange_ptr->count_comps = 0;
  exchange_ptr->comps =
    (struct exch_comp *) PHRQ_malloc (sizeof (struct exch_comp));
  if (exchange_ptr->comps == NULL)
    malloc_error ();
  exchange_ptr->related_phases = FALSE;
  exchange_ptr->related_rate = FALSE;
  exchange_ptr->pitzer_exchange_gammas = TRUE;
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
exchange_ptr_to_user (struct exchange *exchange_ptr_old, int n_user_new)
/* ---------------------------------------------------------------------- */
{
/*
 *   Checks if n_user_new exists, and frees space if it does
 *   Copies exchange_old_ptr to old n_user_new space if
 *   found or to exchange[count_exchange] if not found.
 *   Exchange array may not be in sort order after the copy.
 */
  int sort;
  int n_new;
  struct exchange *exchange_ptr_new;
/*
 *   Find n_user_old in structure array exchange
 */
  if (exchange_ptr_old == NULL)
  {
    sprintf (error_string, "Exchange ptr is NULL");
    error_msg (error_string, CONTINUE);
    input_error++;
    return (ERROR);
  }
/*
 *   Find n_user_new in structure array exchange or make new space
 */
  sort = FALSE;
  exchange_ptr_new = exchange_bsearch (n_user_new, &n_new);
  if (exchange_ptr_new == exchange_ptr_old)
    return (OK);
  if (exchange_ptr_new != NULL)
  {
    exchange_free (exchange_ptr_new);
  }
  else
  {
    space ((void **) ((void *) &exchange), count_exchange, &max_exchange,
	   sizeof (struct exchange));
    if (n_user_new < exchange[count_exchange - 1].n_user)
      sort = TRUE;
    n_new = count_exchange++;
  }
/*
 *   Copy data
 */
  exchange_copy (exchange_ptr_old, &exchange[n_new], n_user_new);
  if (sort == TRUE)
    exchange_sort ();
  return (OK);
}

/* ---------------------------------------------------------------------- */
struct exchange *
exchange_replicate (struct exchange *exchange_old_ptr, int n_user_new)
/* ---------------------------------------------------------------------- */
{
  struct exchange *exchange_ptr;
  exchange_ptr = exchange_alloc ();
  exchange_copy (exchange_old_ptr, exchange_ptr, n_user_new);
  return (exchange_ptr);
}

/* ---------------------------------------------------------------------- */
struct exchange *
exchange_search (int n_user, int *n, int print)
/* ---------------------------------------------------------------------- */
{
/*
 *  Does linear search for user number n_user.
 *  Returns pointer to exchange structure and position number, n,
 *     in exchange array if found.
 *  Returns NULL if not found.
 */
  int i;
  for (i = 0; i < count_exchange; i++)
  {
    if (n_user == exchange[i].n_user)
    {
      break;
    }
  }
  if (i >= count_exchange)
  {
    if (print == TRUE)
    {
      sprintf (error_string, "Exchange %d not found.", n_user);
      error_msg (error_string, CONTINUE);
    }
    *n = -999;
    return (NULL);
  }
  *n = i;
  return (&exchange[i]);
}

/* ---------------------------------------------------------------------- */
int
exchange_sort (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Sort array of exchange structures
 */
  if (count_exchange > 0)
  {
    qsort (exchange, (size_t) count_exchange,
	   (size_t) sizeof (struct exchange), exchange_compare);
  }
  return (OK);
}

/* **********************************************************************
 *
 *   Routines related to structure "gas"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
int
gas_comp_compare (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
  const struct gas_comp *gas_comp_ptr1, *gas_comp_ptr2;
  gas_comp_ptr1 = (const struct gas_comp *) ptr1;
  gas_comp_ptr2 = (const struct gas_comp *) ptr2;
  return (strcmp_nocase (gas_comp_ptr1->name, gas_comp_ptr2->name));

}

/* ---------------------------------------------------------------------- */
struct gas_phase *
gas_phase_alloc (void)
/* ---------------------------------------------------------------------- */
{
  struct gas_phase *gas_phase_ptr;
  gas_phase_ptr =
    (struct gas_phase *) PHRQ_malloc (sizeof (struct gas_phase));
  if (gas_phase_ptr == NULL)
    malloc_error ();
  gas_phase_ptr->n_user = -1;
  gas_phase_ptr->n_user_end = -1;
  gas_phase_ptr->description = NULL;
  gas_phase_ptr->new_def = 0;
  gas_phase_ptr->solution_equilibria = 0;
  gas_phase_ptr->n_solution = 0;
  gas_phase_ptr->type = 0;
  gas_phase_ptr->total_p = 0;
  gas_phase_ptr->total_moles = 0;
  gas_phase_ptr->volume = 0;
  gas_phase_ptr->temperature = 0;
  gas_phase_ptr->count_comps = 0;
  gas_phase_ptr->comps = NULL;
  return (gas_phase_ptr);
}

/* ---------------------------------------------------------------------- */
struct gas_phase *
gas_phase_bsearch (int k, int *n)
/* ---------------------------------------------------------------------- */
{
/*
 *   Binary search gas_phase to find if user gas_phase number k exists.
 *   Assumes array gas_phase is in sort order.
 */
  void *void_ptr;
  if (count_gas_phase == 0)
  {
    *n = -999;
    return (NULL);
  }
  void_ptr = (void *)
    bsearch ((char *) &k,
	     (char *) gas_phase,
	     (size_t) count_gas_phase,
	     (size_t) sizeof (struct gas_phase), gas_phase_compare_int);
  if (void_ptr == NULL)
  {
    *n = -999;
    return (NULL);
  }
  *n = (int) ((struct gas_phase *) void_ptr - gas_phase);
  return ((struct gas_phase *) void_ptr);
}

/* ---------------------------------------------------------------------- */
int
gas_phase_compare (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
  const struct gas_phase *gas_phase_ptr1, *gas_phase_ptr2;
  gas_phase_ptr1 = (const struct gas_phase *) ptr1;
  gas_phase_ptr2 = (const struct gas_phase *) ptr2;
  if (gas_phase_ptr1->n_user > gas_phase_ptr2->n_user)
    return (1);
  if (gas_phase_ptr1->n_user < gas_phase_ptr2->n_user)
    return (-1);
  return (0);

}

/* ---------------------------------------------------------------------- */
static int
gas_phase_compare_int (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
/*
 *   Compare gas phase user numbers
 */
  const int *nptr1;
  const struct gas_phase *nptr2;

  nptr1 = (const int *) ptr1;
  nptr2 = (const struct gas_phase *) ptr2;
  if (*nptr1 > nptr2->n_user)
    return (1);
  if (*nptr1 < nptr2->n_user)
    return (-1);
  return (0);
}

/* ---------------------------------------------------------------------- */
int
gas_phase_copy (struct gas_phase *gas_phase_old_ptr,
		struct gas_phase *gas_phase_new_ptr, int n_user_new)
/* ---------------------------------------------------------------------- */
{
/*
 *   Copies gas_phase structure from gas_phase_old_ptr to new location gas_phase_new_ptr.
 *   Space for a gas_phase structure must already be malloced
 *   Space for gas_phase components is alloced here.
 */
  int count_comps;
  char token[MAX_LENGTH];
/*
 *   Store data for structure gas
 */
  memcpy (gas_phase_new_ptr, gas_phase_old_ptr,
	  (size_t) sizeof (struct gas_phase));
  count_comps = gas_phase_old_ptr->count_comps;
  gas_phase_new_ptr->n_user = n_user_new;
  gas_phase_new_ptr->n_user_end = n_user_new;
  gas_phase_new_ptr->new_def = FALSE;
  sprintf (token, "Gas defined in simulation %d.", simulation);
  gas_phase_new_ptr->description = string_duplicate (token);
/*
 *    Allocate space for gas components
 */
  gas_phase_new_ptr->comps =
    (struct gas_comp *) PHRQ_malloc ((size_t) (count_comps) *
				     sizeof (struct gas_comp));
  if (gas_phase_new_ptr->comps == NULL)
    malloc_error ();
/*
 *   Write gas_comp structure for each gas component
 */
  memcpy (gas_phase_new_ptr->comps, gas_phase_old_ptr->comps,
	  (size_t) (count_comps) * sizeof (struct gas_comp));
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
gas_phase_copy_to_last (int n, int n_user)
/* ---------------------------------------------------------------------- */
{
/*
 *   Copies an gas_phase definition to position count_gas_phase.
 */
  space ((void **) ((void *) &gas_phase), count_gas_phase, &max_gas_phase,
	 sizeof (struct gas_phase));
  gas_phase_copy (&gas_phase[n], &gas_phase[count_gas_phase], n_user);
  count_gas_phase++;
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
gas_phase_delete (int n_user_old)
/* ---------------------------------------------------------------------- */
/*
 *   Frees space for user number n_user_old, removes structure from
 *   array gas_phase.
 */
{
  int i;
  int n_old;
  struct gas_phase *gas_phase_ptr_old;
/*
 *   Find n_user_old in structure array
 */
  gas_phase_ptr_old = gas_phase_bsearch (n_user_old, &n_old);
  if (gas_phase_ptr_old != NULL)
  {
    /*
     *   Delete gas_phase
     */
    gas_phase_free (&gas_phase[n_old]);

    for (i = n_old + 1; i < count_gas_phase; i++)
    {
      memcpy ((void *) &gas_phase[i - 1], (void *) &gas_phase[i],
	      (size_t) sizeof (struct gas_phase));
    }
    count_gas_phase--;
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
gas_phase_duplicate (int n_user_old, int n_user_new)
/* ---------------------------------------------------------------------- */
{
/*
 *   Checks if n_user_new exists, and frees space if it does
 *   Copies gas_phase[n_user_old] to old n_user_new space if
 *   found or to gas_phase[count_gas_phase] if not found.
 *   Gas_Phase array may not be in sort order after the copy.
 */
  int n_old, n_new, sort;
  struct gas_phase *gas_phase_ptr_old, *gas_phase_ptr_new;
/*
 *   Find n_user_old in structure array gas_phase
 */
  if (n_user_old == n_user_new)
    return (OK);
  gas_phase_ptr_old = gas_phase_bsearch (n_user_old, &n_old);
  if (gas_phase_ptr_old == NULL)
  {
    sprintf (error_string, "Gas_Phase %d not found.", n_user_old);
    error_msg (error_string, CONTINUE);
    input_error++;
    return (ERROR);
  }
/*
 *   Find n_user_new in structure array gas_phase or make new space
 */
  sort = FALSE;
  gas_phase_ptr_new = gas_phase_bsearch (n_user_new, &n_new);
  if (gas_phase_ptr_new != NULL)
  {
    gas_phase_free (&gas_phase[n_new]);
#ifdef SKIP
    gas_phase[n_new].comps = free_check_null (gas_phase[n_new].comps);
#endif
  }
  else
  {
    space ((void **) ((void *) &gas_phase), count_gas_phase, &max_gas_phase,
	   sizeof (struct gas_phase));
    if (n_user_new < gas_phase[count_gas_phase - 1].n_user)
      sort = TRUE;
    n_new = count_gas_phase++;
  }
/*
 *   Copy data
 */
  gas_phase_copy (&gas_phase[n_old], &gas_phase[n_new], n_user_new);
  if (sort == TRUE)
    gas_phase_sort ();
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
gas_phase_free (struct gas_phase *gas_phase_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Frees space associated with gas_phase_ptr.
 */
  if (gas_phase_ptr == NULL)
    return (ERROR);
/*
 *   Free space allocated for gas_phase structure
 */
  gas_phase_ptr->description =
    (char *) free_check_null (gas_phase_ptr->description);
  gas_phase_ptr->comps =
    (struct gas_comp *) free_check_null (gas_phase_ptr->comps);
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
gas_phase_init (struct gas_phase *gas_phase_ptr, int n_user, int n_user_end,
		char *description)
/* ---------------------------------------------------------------------- */
{
/*
 *   Frees space associated with gas_phase_ptr.
 */
  if (gas_phase_ptr == NULL)
    return (ERROR);
  gas_phase_ptr->n_user = n_user;
  gas_phase_ptr->n_user_end = n_user_end;
  gas_phase_ptr->description = string_duplicate (description);
  gas_phase_ptr->new_def = TRUE;
  gas_phase_ptr->solution_equilibria = FALSE;
  gas_phase_ptr->n_solution = 0;
  gas_phase_ptr->type = PRESSURE;
  gas_phase_ptr->total_p = 1.0;
  gas_phase_ptr->total_moles = 0.0;
  gas_phase_ptr->volume = 1.0;
  gas_phase_ptr->temperature = 298.15;
  gas_phase_ptr->count_comps = 0;
  gas_phase_ptr->comps =
    (struct gas_comp *) PHRQ_malloc ((size_t) sizeof (struct gas_comp));
  if (gas_phase_ptr->comps == NULL)
    malloc_error ();

  return (OK);
}

/* ---------------------------------------------------------------------- */
int
gas_phase_ptr_to_user (struct gas_phase *gas_phase_ptr_old, int n_user_new)
/* ---------------------------------------------------------------------- */
{
/*
 *   Checks if n_user_new exists, and frees space if it does
 *   Copies gas_phase_ptr_old to old n_user_new space if
 *   found or to gas_phase[count_gas_phase] if not found.
 *   Gas_Phase array may not be in sort order after the copy.
 */
  int n_new, sort;
  struct gas_phase *gas_phase_ptr_new;
/*
 *   Find n_user_old in structure array gas_phase
 */
  if (gas_phase_ptr_old == NULL)
  {
    sprintf (error_string, "Gas_Phase pointer is NULL.");
    error_msg (error_string, CONTINUE);
    input_error++;
    return (ERROR);
  }
/*
 *   Find n_user_new in structure array gas_phase or make new space
 */
  sort = FALSE;
  gas_phase_ptr_new = gas_phase_bsearch (n_user_new, &n_new);
  if (gas_phase_ptr_new == gas_phase_ptr_old)
    return (OK);
  if (gas_phase_ptr_new != NULL)
  {
    gas_phase_free (&gas_phase[n_new]);
#ifdef SKIP
    gas_phase[n_new].comps = free_check_null (gas_phase[n_new].comps);
#endif
  }
  else
  {
    space ((void **) ((void *) &gas_phase), count_gas_phase, &max_gas_phase,
	   sizeof (struct gas_phase));
    if (n_user_new < gas_phase[count_gas_phase - 1].n_user)
      sort = TRUE;
    n_new = count_gas_phase++;
  }
/*
 *   Copy data
 */
  gas_phase_copy (gas_phase_ptr_old, &gas_phase[n_new], n_user_new);
  if (sort == TRUE)
    gas_phase_sort ();
  return (OK);
}

/* ---------------------------------------------------------------------- */
struct gas_phase *
gas_phase_replicate (struct gas_phase *gas_phase_old_ptr, int n_user_new)
/* ---------------------------------------------------------------------- */
{
  struct gas_phase *gas_phase_ptr;
  gas_phase_ptr = gas_phase_alloc ();
  gas_phase_copy (gas_phase_old_ptr, gas_phase_ptr, n_user_new);
  return (gas_phase_ptr);
}

/* ---------------------------------------------------------------------- */
struct gas_phase *
gas_phase_search (int n_user, int *n)
/* ---------------------------------------------------------------------- */
{
/*   Linear search of the structure array "gas_phase" for user number n_user
 *
 *   Arguments:
 *      n_user  input, user number
 *      n       output, position in gas_phase
 *
 *   Returns:
 *      if found, the address of the gas_phase element
 *      if not found, NULL
 *
 */
  int i;
  for (i = 0; i < count_gas_phase; i++)
  {
    if (gas_phase[i].n_user == n_user)
    {
      *n = i;
      return (&(gas_phase[i]));
    }
  }
/*
 *   Gas_phase not found
 */
  return (NULL);
}

/* ---------------------------------------------------------------------- */
int
gas_phase_sort (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Sort array of gas_phase structures
 */
  if (count_gas_phase > 0)
  {
    qsort (gas_phase, (size_t) count_gas_phase,
	   (size_t) sizeof (struct gas_phase), gas_phase_compare);
  }
  return (OK);
}

/* **********************************************************************
 *
 *   Routines related to structure "inverse"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
struct inverse *
inverse_alloc (void)
/* ---------------------------------------------------------------------- */
/*
 *   Allocates space for a new inverse structure at position count_inverse.
 *   Initializes structure.
 *      arguments
 *      input:  none
 *      output: pointer to an inverse structure
 *      return: OK
 */
{
  struct inverse *inverse_ptr;

  count_inverse++;
  inverse =
    (struct inverse *) PHRQ_realloc (inverse,
				     (size_t) count_inverse *
				     sizeof (struct inverse));
  if (inverse == NULL)
    malloc_error ();
  inverse_ptr = &(inverse[count_inverse - 1]);
/*
 *   Initialize variables
 */
  inverse_ptr->description = NULL;
  inverse_ptr->count_uncertainties = 0;
  inverse_ptr->count_solns = 0;
  inverse_ptr->count_elts = 0;
  inverse_ptr->count_isotopes = 0;
  inverse_ptr->count_i_u = 0;
  inverse_ptr->count_phases = 0;
  inverse_ptr->count_force_solns = 0;
/*
 *   allocate space for pointers in structure to NULL
 */

  inverse_ptr->uncertainties =
    (double *) PHRQ_malloc ((size_t) sizeof (LDBLE));
  if (inverse_ptr->uncertainties == NULL)
    malloc_error ();

  inverse_ptr->ph_uncertainties =
    (double *) PHRQ_malloc ((size_t) sizeof (LDBLE));
  if (inverse_ptr->ph_uncertainties == NULL)
    malloc_error ();

#ifdef SKIP
  inverse_ptr->alk_uncertainties = PHRQ_malloc ((size_t) sizeof (LDBLE));
  if (inverse_ptr->alk_uncertainties == NULL)
    malloc_error ();
#endif

  inverse_ptr->force_solns = (int *) PHRQ_malloc ((size_t) sizeof (int));
  if (inverse_ptr->force_solns == NULL)
    malloc_error ();

  inverse_ptr->dalk_dph = NULL;
  inverse_ptr->dalk_dc = NULL;

  inverse_ptr->solns = NULL;

  inverse_ptr->elts =
    (struct inv_elts *) PHRQ_malloc ((size_t) sizeof (struct inv_elts));
  if (inverse_ptr->elts == NULL)
    malloc_error ();
  inverse_ptr->elts[0].name = NULL;
  inverse_ptr->elts[0].uncertainties = NULL;

  inverse_ptr->isotopes =
    (struct inv_isotope *) PHRQ_malloc ((size_t) sizeof (struct inv_isotope));
  if (inverse_ptr->isotopes == NULL)
    malloc_error ();
  inverse_ptr->isotopes[0].isotope_name = NULL;
  inverse_ptr->isotopes[0].isotope_number = 0;
  inverse_ptr->isotopes[0].elt_name = NULL;

  inverse_ptr->i_u =
    (struct inv_isotope *) PHRQ_malloc ((size_t) sizeof (struct inv_isotope));
  if (inverse_ptr->i_u == NULL)
    malloc_error ();
  inverse_ptr->i_u[0].isotope_name = NULL;
  inverse_ptr->i_u[0].isotope_number = 0;
  inverse_ptr->i_u[0].elt_name = NULL;

  inverse_ptr->phases =
    (struct inv_phases *) PHRQ_malloc ((size_t) sizeof (struct inv_phases));
  if (inverse_ptr->phases == NULL)
    malloc_error ();

  return (inverse_ptr);
}

/* ---------------------------------------------------------------------- */
int
inverse_compare (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
/*
 *   Compare inverse values for n_user
 */
  const struct inverse *nptr1;
  const struct inverse *nptr2;

  nptr1 = (const struct inverse *) ptr1;
  nptr2 = (const struct inverse *) ptr2;
  if (nptr1->n_user > nptr2->n_user)
    return (1);
  if (nptr1->n_user < nptr2->n_user)
    return (-1);
  return (0);
}

/* ---------------------------------------------------------------------- */
int
inverse_delete (int i)
/* ---------------------------------------------------------------------- */
{
/*
 *   Deletes inverse i from list (i is not user number),
 *   Frees memory allocated to inverse struct
 *   Input: i, number of inverse struct to delete
 *   Return: OK
 */
  int j;

  inverse_free (&(inverse[i]));
  for (j = i; j < (count_inverse - 1); j++)
  {
    memcpy ((void *) &(inverse[j]), (void *) &(inverse[j + 1]),
	    (size_t) sizeof (struct inverse));
  }
  count_inverse--;
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
inverse_free (struct inverse *inverse_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Free all memory for an inverse structure.
 */
  int i;

  inverse_ptr->description =
    (char *) free_check_null (inverse_ptr->description);
/*   Free solns */
  inverse_ptr->solns = (int *) free_check_null (inverse_ptr->solns);

/*   Free uncertainties */
  inverse_ptr->uncertainties =
    (double *) free_check_null (inverse_ptr->uncertainties);
  inverse_ptr->ph_uncertainties =
    (double *) free_check_null (inverse_ptr->ph_uncertainties);
#ifdef SKIP
  inverse_ptr->alk_uncertainties =
    (double *) free_check_null (inverse_ptr->alk_uncertainties);
#endif

/*   Free force_solns */
  inverse_ptr->force_solns =
    (int *) free_check_null (inverse_ptr->force_solns);

/*   Free elts */
  for (i = 0; i < inverse_ptr->count_elts; i++)
  {
    inverse_ptr->elts[i].uncertainties =
      (double *) free_check_null (inverse_ptr->elts[i].uncertainties);
  };
  inverse_ptr->elts = (struct inv_elts *) free_check_null (inverse_ptr->elts);

/*   Free isotopes */
  for (i = 0; i < inverse_ptr->count_isotopes; i++)
  {
    inverse_ptr->isotopes[i].uncertainties =
      (double *) free_check_null (inverse_ptr->isotopes[i].uncertainties);
  };
  inverse_ptr->isotopes =
    (struct inv_isotope *) free_check_null (inverse_ptr->isotopes);

  for (i = 0; i < inverse_ptr->count_i_u; i++)
  {
    inverse_ptr->i_u[i].uncertainties =
      (double *) free_check_null (inverse_ptr->i_u[i].uncertainties);
  };
  inverse_ptr->i_u =
    (struct inv_isotope *) free_check_null (inverse_ptr->i_u);

/*   Free phases */
  for (i = 0; i < inverse_ptr->count_phases; i++)
  {
    inverse_ptr->phases[i].isotopes =
      (struct isotope *) free_check_null (inverse_ptr->phases[i].isotopes);
  }
  inverse_ptr->phases =
    (struct inv_phases *) free_check_null (inverse_ptr->phases);

/*   Free carbon derivatives */
  inverse_ptr->dalk_dph = (LDBLE *) free_check_null (inverse_ptr->dalk_dph);
  inverse_ptr->dalk_dc = (LDBLE *) free_check_null (inverse_ptr->dalk_dc);

  return (OK);
}

/* ---------------------------------------------------------------------- */
int
inverse_isotope_compare (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
  int i;
  const struct inv_isotope *iso_ptr1, *iso_ptr2;

  iso_ptr1 = (const struct inv_isotope *) ptr1;
  iso_ptr2 = (const struct inv_isotope *) ptr2;
  i = strcmp_nocase (iso_ptr1->elt_name, iso_ptr2->elt_name);
  if (i != 0)
    return (i);
  if (iso_ptr1->isotope_number < iso_ptr2->isotope_number)
  {
    return (-1);
  }
  else if (iso_ptr1->isotope_number > iso_ptr2->isotope_number)
  {
    return (1);
  }
  return (0);
}

/* ---------------------------------------------------------------------- */
struct inverse *
inverse_search (int n_user, int *n)
/* ---------------------------------------------------------------------- */
{
/*   Linear search of the structure array "inverse" for user number n_user.
 *
 *   Arguments:
 *      n_user  input, user number
 *      n       output, position in inverse
 *
 *   Returns:
 *      if found, the address of the inverse element
 *      if not found, NULL
 *
 */
  int i;
  for (i = 0; i < count_inverse; i++)
  {
    if (inverse[i].n_user == n_user)
    {
      *n = i;
      return (&(inverse[i]));
    }
  }
/*
 *   An inverse structure with n_user was not found
 */
  return (NULL);
}

/* ---------------------------------------------------------------------- */
int
inverse_sort (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Sort array of inverse structures
 */
  if (count_inverse > 0)
  {
    qsort (inverse, (size_t) count_inverse, (size_t) sizeof (struct inverse),
	   inverse_compare);
  }
  return (OK);
}

/* **********************************************************************
 *
 *   Routines related to structure "irrev"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
struct irrev *
irrev_bsearch (int k, int *n)
/* ---------------------------------------------------------------------- */
{
/*
 *   Search irrev to find if user irrev number k exists.
 *   Uses binary search, assumes array irrev is in sort order.
 */
  void *void_ptr;
  if (count_irrev == 0)
  {
    *n = -999;
    return (NULL);
  }
  void_ptr = (void *)
    bsearch ((char *) &k,
	     (char *) irrev,
	     (size_t) count_irrev,
	     (size_t) sizeof (struct irrev), irrev_compare_int);
  if (void_ptr == NULL)
  {
    *n = -999;
    return (NULL);
  }
  *n = (int) ((struct irrev *) void_ptr - irrev);
  return ((struct irrev *) void_ptr);
}

/* ---------------------------------------------------------------------- */
int
irrev_compare (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
  const struct irrev *irrev_ptr1, *irrev_ptr2;
  irrev_ptr1 = (const struct irrev *) ptr1;
  irrev_ptr2 = (const struct irrev *) ptr2;
  if (irrev_ptr1->n_user > irrev_ptr2->n_user)
    return (1);
  if (irrev_ptr1->n_user < irrev_ptr2->n_user)
    return (-1);
  return (0);

}

/* ---------------------------------------------------------------------- */
static int
irrev_compare_int (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
/*
 *   Compare irrev user numbers
 */
  const int *nptr1;
  const struct irrev *nptr2;

  nptr1 = (const int *) ptr1;
  nptr2 = (const struct irrev *) ptr2;
  if (*nptr1 > nptr2->n_user)
    return (1);
  if (*nptr1 < nptr2->n_user)
    return (-1);
  return (0);
}

/* ---------------------------------------------------------------------- */
int
irrev_copy (struct irrev *irrev_old_ptr, struct irrev *irrev_new_ptr,
	    int n_user_new)
/* ---------------------------------------------------------------------- */
{
/*
 *   Copies irrev structure from irrev_old_ptr to new location irrev_new_ptr.
 *   Space for a irrev structure must already be malloced
 *   Space for irrev components is alloced here.
 */
  int count_list, count_steps;
  char token[MAX_LENGTH];
/*
 *   Store data for structure gas
 */
  memcpy (irrev_new_ptr, irrev_old_ptr, (size_t) sizeof (struct irrev));
  irrev_new_ptr->n_user = n_user_new;
  irrev_new_ptr->n_user_end = n_user_new;
  sprintf (token, "Reaction defined in simulation %d.", simulation);
  irrev_new_ptr->description = string_duplicate (token);
  count_list = irrev_old_ptr->count_list;
  count_steps = irrev_old_ptr->count_steps;
/*
 *   Allocate space
 */
  irrev_new_ptr->list =
    (struct name_coef *) PHRQ_malloc ((size_t) (count_list) *
				      sizeof (struct name_coef));
  if (irrev_new_ptr->list == NULL)
    malloc_error ();
  memcpy (irrev_new_ptr->list, irrev_old_ptr->list,
	  (size_t) (count_list) * sizeof (struct name_coef));
  if (count_steps < 0)
    count_steps = 1;
  irrev_new_ptr->steps =
    (LDBLE *) PHRQ_malloc ((size_t) (count_steps) * sizeof (LDBLE));
  if (irrev_new_ptr->steps == NULL)
    malloc_error ();
  memcpy (irrev_new_ptr->steps, irrev_old_ptr->steps,
	  (size_t) (count_steps) * sizeof (LDBLE));

  if (irrev_old_ptr->elts != NULL)
  {
    irrev_new_ptr->elts = elt_list_dup (irrev_old_ptr->elts);
  }
  else
  {
    irrev_new_ptr->elts = NULL;
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
irrev_duplicate (int n_user_old, int n_user_new)
/* ---------------------------------------------------------------------- */
{
/*
 *   Finds irrev structure with user number n_user_old,
 *   copies data to irrev structure with user number n_user_new.
 *   If n_user_new already exists, it is overwritten and sort order
 *      is maintained.
 *   If n_user_new does not exist, a new structure is allocated; sort
 *      order may not be preserved.
 */
  int n_old, n_new, sort;
  char token[MAX_LENGTH];
  int count_steps, count_list;
  struct irrev *irrev_ptr_old, *irrev_ptr_new;
/*
 *   Find n_user_old in structure array irrev_assemblage or make new space
 */
  if (n_user_old == n_user_new)
    return (OK);
  irrev_ptr_old = irrev_bsearch (n_user_old, &n_old);
  if (irrev_ptr_old == NULL)
  {
    sprintf (error_string, "Irreversible reaction %d not found.", n_user_old);
    error_msg (error_string, CONTINUE);
    input_error++;
    return (ERROR);
  }
/*
 *   Find n_user_old in structure array irrev or make new space
 */
  sort = FALSE;
  irrev_ptr_new = irrev_bsearch (n_user_new, &n_new);
  if (irrev_ptr_new != NULL)
  {
    irrev_free (irrev_ptr_new);
  }
  else
  {
    irrev =
      (struct irrev *) PHRQ_realloc (irrev,
				     (size_t) (count_irrev +
					       1) * sizeof (struct irrev));
    if (irrev == NULL)
      malloc_error ();
    if (n_user_new < irrev[count_irrev - 1].n_user)
      sort = TRUE;
    n_new = count_irrev++;
  }
/*
 *   Store data for structure irrev
 */
  count_steps = irrev[n_old].count_steps;
  irrev[n_new].n_user = n_user_new;
  irrev[n_new].n_user_end = n_user_new;
  sprintf (token, "Irreversible reaction defined in simulation %d.",
	   simulation);
  irrev[n_new].description = string_duplicate (token);
  irrev[n_new].units = irrev[n_old].units;
  irrev[n_new].count_steps = count_steps;
  count_list = irrev[n_old].count_list;
  irrev[n_new].count_list = irrev[n_old].count_list;
/*
 *   Allocate space
 */
  irrev[n_new].list =
    (struct name_coef *) PHRQ_malloc ((size_t) (count_list) *
				      sizeof (struct name_coef));
  if (irrev[n_new].list == NULL)
    malloc_error ();
  memcpy (irrev[n_new].list, irrev[n_old].list,
	  (size_t) (count_list) * sizeof (struct name_coef));

  if (count_steps < 0)
    count_steps = 1;
  irrev[n_new].steps =
    (LDBLE *) PHRQ_malloc ((size_t) (count_steps) * sizeof (LDBLE));
  if (irrev[n_new].steps == NULL)
    malloc_error ();
  memcpy (irrev[n_new].steps, irrev[n_old].steps,
	  (size_t) (count_steps) * sizeof (LDBLE));

  if (irrev[n_old].elts != NULL)
  {
    irrev[n_new].elts = elt_list_dup (irrev[n_old].elts);
  }
  else
  {
    irrev[n_new].elts = NULL;
  }
  if (sort == TRUE)
    irrev_sort ();
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
irrev_free (struct irrev *irrev_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Free memory pointed to by an irrev structure.
 */
  if (irrev_ptr == NULL)
    return (ERROR);
/*
 *   Free space allocated for irrev structure
 */
  irrev_ptr->description = (char *) free_check_null (irrev_ptr->description);
  irrev_ptr->list = (struct name_coef *) free_check_null (irrev_ptr->list);
  irrev_ptr->elts = (struct elt_list *) free_check_null (irrev_ptr->elts);
  irrev_ptr->steps = (LDBLE *) free_check_null (irrev_ptr->steps);
  return (OK);
}

/* ---------------------------------------------------------------------- */
struct irrev *
irrev_search (int n_user, int *n)
/* ---------------------------------------------------------------------- */
{
/*
 *   Linear search for irrev user number n_user.
 *   If found, n is position in irrev array and pointer to irrev structure
 *      is returned
 *   If not found, NULL is returned.
 */
  int i;
  for (i = 0; i < count_irrev; i++)
  {
    if (n_user == irrev[i].n_user)
    {
      break;
    }
  }
  if (i >= count_irrev)
  {
    *n = -999;
    return (NULL);
  }
  *n = i;
  return (&irrev[i]);
}

/* ---------------------------------------------------------------------- */
int
irrev_sort (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Sort array of irrev structures
 */
  if (count_irrev > 0)
  {
    qsort (irrev, (size_t) count_irrev, (size_t) sizeof (struct irrev),
	   irrev_compare);
  }
  return (OK);
}

/* **********************************************************************
 *
 *   Routines related to structure "kinetics"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
struct kinetics *
kinetics_alloc (void)
/* ---------------------------------------------------------------------- */
{
  struct kinetics *kinetics_ptr;
  kinetics_ptr = (struct kinetics *) PHRQ_malloc (sizeof (struct kinetics));
  if (kinetics_ptr == NULL)
    malloc_error ();
  kinetics_ptr->n_user = -1;
  kinetics_ptr->n_user_end = -1;
  kinetics_ptr->description = NULL;
  kinetics_ptr->count_comps = 0;
  kinetics_ptr->comps = NULL;
  kinetics_ptr->count_steps = 0;
  kinetics_ptr->steps = NULL;
  kinetics_ptr->totals = NULL;
  kinetics_ptr->rk = 0;
  kinetics_ptr->bad_step_max = 0;
  kinetics_ptr->use_cvode = FALSE;

  return (kinetics_ptr);
}

/* ---------------------------------------------------------------------- */
struct kinetics *
kinetics_bsearch (int k, int *n)
/* ---------------------------------------------------------------------- */
{
/*
 *   Binary search kinetics to find if user kinetics number k exists.
 *   Kinetics is assumed to be in sort order by user number.
 */
  void *void_ptr;
  if (count_kinetics == 0)
  {
    *n = -999;
    return (NULL);
  }
  void_ptr = (void *)
    bsearch ((char *) &k,
	     (char *) kinetics,
	     (size_t) count_kinetics,
	     (size_t) sizeof (struct kinetics), kinetics_compare_int);
  if (void_ptr == NULL)
  {
    *n = -999;
    return (NULL);
  }
  *n = (int) ((struct kinetics *) void_ptr - kinetics);
  return ((struct kinetics *) void_ptr);
}

/* ---------------------------------------------------------------------- */
int
kinetics_compare (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
  const struct kinetics *kinetics_ptr1, *kinetics_ptr2;
  kinetics_ptr1 = (const struct kinetics *) ptr1;
  kinetics_ptr2 = (const struct kinetics *) ptr2;
  if (kinetics_ptr1->n_user > kinetics_ptr2->n_user)
    return (1);
  if (kinetics_ptr1->n_user < kinetics_ptr2->n_user)
    return (-1);
  return (0);

}

/* ---------------------------------------------------------------------- */
static int
kinetics_compare_int (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
/*
 *   Compare kinetics user numbers
 */
  const int *nptr1;
  const struct kinetics *nptr2;

  nptr1 = (const int *) ptr1;
  nptr2 = (const struct kinetics *) ptr2;
  if (*nptr1 > nptr2->n_user)
    return (1);
  if (*nptr1 < nptr2->n_user)
    return (-1);
  return (0);
}

/* ---------------------------------------------------------------------- */
int
kinetics_copy (struct kinetics *kinetics_old_ptr,
	       struct kinetics *kinetics_new_ptr, int n_user_new)
/* ---------------------------------------------------------------------- */
{
/*
 *   Copies kinetics data from kinetics_old_ptr to new location, kinetics_new_ptr.
 *   Space for the kinetics_new_ptr structure must already be malloced.
 *   Space for kinetics components is malloced here.
 */
  char token[MAX_LENGTH];
  int count_steps, i;
/*
 *   copy old to new
 */
  memcpy (kinetics_new_ptr, kinetics_old_ptr, sizeof (struct kinetics));
/*
 *   Store data for structure kinetics
 */
  kinetics_new_ptr->n_user = n_user_new;
  kinetics_new_ptr->n_user_end = n_user_new;
  sprintf (token, "Kinetics defined in simulation %d.", simulation);
  kinetics_new_ptr->description = string_duplicate (token);
  kinetics_new_ptr->totals = NULL;
  kinetics_new_ptr->totals = elt_list_dup (kinetics_old_ptr->totals);
/*
 *   Copy time steps
 */
  if (kinetics_new_ptr->count_steps > 0)
  {
    count_steps = kinetics_new_ptr->count_steps;
  }
  else if (kinetics_new_ptr->count_steps < 0)
  {
    count_steps = 1;
  }
  else
  {
    count_steps = 0;
  }
  if (count_steps > 0)
  {
    kinetics_new_ptr->steps =
      (LDBLE *) PHRQ_malloc ((size_t) (count_steps * sizeof (LDBLE)));
    if (kinetics_new_ptr->steps == NULL)
      malloc_error ();
    memcpy (kinetics_new_ptr->steps, kinetics_old_ptr->steps,
	    (size_t) count_steps * sizeof (LDBLE));
  }
/*
 *   Copy kinetic components
 */
  if (kinetics_new_ptr->count_comps > 0)
  {
    kinetics_new_ptr->comps =
      (struct kinetics_comp *) PHRQ_malloc ((size_t) kinetics_old_ptr->
					    count_comps *
					    sizeof (struct kinetics_comp));
    if (kinetics_new_ptr->comps == NULL)
      malloc_error ();
  }
  else
  {
    kinetics_new_ptr->comps = NULL;
  }
  for (i = 0; i < kinetics_old_ptr->count_comps; i++)
  {
    kinetics_comp_duplicate (&kinetics_new_ptr->comps[i],
			     &kinetics_old_ptr->comps[i]);
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
kinetics_comp_duplicate (struct kinetics_comp *kinetics_comp_new_ptr,
			 struct kinetics_comp *kinetics_comp_old_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Copies kinetics comp data from kinetics_comp_old_ptr to new location
 *   Space for the new structure must already be allocated
 */
  memcpy (kinetics_comp_new_ptr, kinetics_comp_old_ptr,
	  sizeof (struct kinetics_comp));
/*
 *   Kinetics character parameters
 */
  if (kinetics_comp_new_ptr->count_c_params > 0)
  {
    kinetics_comp_new_ptr->c_params =
      (char **) PHRQ_malloc ((size_t) (kinetics_comp_old_ptr->count_c_params)
			     * sizeof (char *));
    if (kinetics_comp_new_ptr->c_params == NULL)
      malloc_error ();
    memcpy (kinetics_comp_new_ptr->c_params, kinetics_comp_old_ptr->c_params,
	    (size_t) (kinetics_comp_old_ptr->count_c_params) *
	    sizeof (char *));
  }
  else
  {
    kinetics_comp_new_ptr->c_params = NULL;
  }
/*
 *   Kinetics LDBLE parameters
 */
  if (kinetics_comp_new_ptr->count_d_params > 0)
  {
    kinetics_comp_new_ptr->d_params =
      (LDBLE *) PHRQ_malloc ((size_t) (kinetics_comp_old_ptr->count_d_params)
			     * sizeof (LDBLE));
    if (kinetics_comp_new_ptr->d_params == NULL)
      malloc_error ();
    memcpy (kinetics_comp_new_ptr->d_params, kinetics_comp_old_ptr->d_params,
	    (size_t) (kinetics_comp_old_ptr->count_d_params) *
	    sizeof (LDBLE));
  }
  else
  {
    kinetics_comp_new_ptr->d_params =
      (LDBLE *) PHRQ_malloc ((size_t) sizeof (LDBLE));
  }
/*
 *   Kinetics list of formulae
 */
  if (kinetics_comp_new_ptr->count_list > 0)
  {
    kinetics_comp_new_ptr->list =
      (struct name_coef *)
      PHRQ_malloc ((size_t) (kinetics_comp_old_ptr->count_list) *
		   sizeof (struct name_coef));
    if (kinetics_comp_new_ptr->list == NULL)
      malloc_error ();
    memcpy (kinetics_comp_new_ptr->list, kinetics_comp_old_ptr->list,
	    (size_t) (kinetics_comp_old_ptr->count_list) *
	    sizeof (struct name_coef));
  }
  else
  {
    kinetics_comp_new_ptr->d_params = NULL;
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
kinetics_copy_to_last (int n, int n_user)
/* ---------------------------------------------------------------------- */
{
/*
 *   Copies an kinetics definition from position n to position count_kinetics.
 *   New kinetics structure is given user number n_user.
 */
  space ((void **) ((void *) &kinetics), count_kinetics, &max_kinetics,
	 sizeof (struct kinetics));
  kinetics_copy (&kinetics[n], &kinetics[count_kinetics], n_user);
  count_kinetics++;
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
kinetics_delete (int n_user_old)
/* ---------------------------------------------------------------------- */
/*
 *   Frees space for user number n_user_old, removes structure from
 *   array kinetics.
 */
{
  int i;
  int n_old;
  struct kinetics *kinetics_ptr_old;
/*
 *   Find n_user_old in structure array
 */
  kinetics_ptr_old = kinetics_bsearch (n_user_old, &n_old);
  if (kinetics_ptr_old != NULL)
  {
    /*
     *   Delete kinetics
     */
    kinetics_free (&kinetics[n_old]);

    for (i = n_old + 1; i < count_kinetics; i++)
    {
      memcpy ((void *) &kinetics[i - 1], (void *) &kinetics[i],
	      (size_t) sizeof (struct kinetics));
    }
    count_kinetics--;
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
kinetics_duplicate (int n_user_old, int n_user_new)
/* ---------------------------------------------------------------------- */
{
/*
 *   Checks if n_user_new exists, and frees space if it does
 *   Copies kinetics[n_user_old] to old n_user_new space if
 *   found or to kinetics[count_kinetics] if not found.
 *   Kinetics array may not be in sort order after the copy.
 */
  int n_old, n_new, sort;
  struct kinetics *kinetics_ptr_old, *kinetics_ptr_new;
/*
 *   Find n_user_old in structure array kinetics
 */
  if (n_user_old == n_user_new)
    return (OK);
  kinetics_ptr_old = kinetics_bsearch (n_user_old, &n_old);
  if (kinetics_ptr_old == NULL)
  {
    sprintf (error_string, "Kinetics %d not found.", n_user_old);
    error_msg (error_string, CONTINUE);
    input_error++;
    return (ERROR);
  }
/*
 *   Find n_user_new in structure array kinetics or make new space
 */
  sort = FALSE;
  kinetics_ptr_new = kinetics_bsearch (n_user_new, &n_new);
  if (kinetics_ptr_new != NULL)
  {
    kinetics_free (kinetics_ptr_new);
  }
  else
  {
    space ((void **) ((void *) &kinetics), count_kinetics, &max_kinetics,
	   sizeof (struct kinetics));
    if (n_user_new < kinetics[count_kinetics - 1].n_user)
      sort = TRUE;
    n_new = count_kinetics++;
  }
/*
 *   Copy data
 */
  kinetics_copy (&kinetics[n_old], &kinetics[n_new], n_user_new);
  if (sort == TRUE)
    kinetics_sort ();
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
kinetics_free (struct kinetics *kinetics_ptr)
/* ---------------------------------------------------------------------- */
{
  int i;
/*
 *   Frees all data associated with kinetics structure.
 */
  if (kinetics_ptr == NULL)
    return (ERROR);
/*
 *   Free space allocated for kinetics structure
 */
  for (i = 0; i < kinetics_ptr->count_comps; i++)
  {
    kinetics_ptr->comps[i].c_params =
      (char **) free_check_null (kinetics_ptr->comps[i].c_params);
    kinetics_ptr->comps[i].d_params =
      (LDBLE *) free_check_null (kinetics_ptr->comps[i].d_params);
    kinetics_ptr->comps[i].list =
      (struct name_coef *) free_check_null (kinetics_ptr->comps[i].list);
  }
  kinetics_ptr->description =
    (char *) free_check_null (kinetics_ptr->description);
  kinetics_ptr->steps = (LDBLE *) free_check_null (kinetics_ptr->steps);
  kinetics_ptr->comps =
    (struct kinetics_comp *) free_check_null (kinetics_ptr->comps);
  kinetics_ptr->totals =
    (struct elt_list *) free_check_null (kinetics_ptr->totals);
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
kinetics_init (struct kinetics *kinetics_ptr, int n_user, int n_user_end,
	       char *description)
/* ---------------------------------------------------------------------- */
{
  if (kinetics_ptr == NULL)
    return (ERROR);
  kinetics_ptr->n_user = n_user;
  kinetics_ptr->n_user_end = n_user_end;
  kinetics_ptr->description = string_duplicate (description);
  kinetics_ptr->count_comps = 0;
  kinetics_ptr->comps =
    (struct kinetics_comp *) PHRQ_malloc ((size_t)
					  sizeof (struct kinetics_comp));
  if (kinetics_ptr->comps == NULL)
    malloc_error ();
  kinetics_ptr->count_steps = 0;
  kinetics_ptr->steps = (LDBLE *) PHRQ_malloc ((size_t) sizeof (LDBLE));
  if (kinetics_ptr->steps == NULL)
    malloc_error ();
  kinetics_ptr->step_divide = 1.0;
  /*kinetics_ptr->units = string_hsave("sec"); */
  kinetics_ptr->totals = NULL;
  kinetics_ptr->rk = 3;
  kinetics_ptr->bad_step_max = 500;
  kinetics_ptr->use_cvode = FALSE;
  kinetics_ptr->cvode_order = 5;
  kinetics_ptr->cvode_steps = 100;

  return (OK);
}

/* ---------------------------------------------------------------------- */
int
kinetics_ptr_to_user (struct kinetics *kinetics_ptr_old, int n_user_new)
/* ---------------------------------------------------------------------- */
{
/*
 *   Checks if n_user_new exists, and frees space if it does
 *   Copies kinetics_ptr_old to old n_user_new space if
 *   found or to kinetics[count_kinetics] if not found.
 *   Kinetics array may not be in sort order after the copy.
 */
  int n_new, sort;
  struct kinetics *kinetics_ptr_new;
/*
 *   Find n_user_old in structure array kinetics
 */
  if (kinetics_ptr_old == NULL)
  {
    sprintf (error_string, "Kinetics pointer is NULL.");
    error_msg (error_string, CONTINUE);
    input_error++;
    return (ERROR);
  }
/*
 *   Find n_user_new in structure array kinetics or make new space
 */
  sort = FALSE;
  kinetics_ptr_new = kinetics_bsearch (n_user_new, &n_new);
  if (kinetics_ptr_new == kinetics_ptr_old)
    return (OK);
  if (kinetics_ptr_new != NULL)
  {
    kinetics_free (kinetics_ptr_new);
  }
  else
  {
    space ((void **) ((void *) &kinetics), count_kinetics, &max_kinetics,
	   sizeof (struct kinetics));
    if (n_user_new < kinetics[count_kinetics - 1].n_user)
      sort = TRUE;
    n_new = count_kinetics++;
  }
/*
 *   Copy data
 */
  kinetics_copy (kinetics_ptr_old, &kinetics[n_new], n_user_new);
  if (sort == TRUE)
    kinetics_sort ();
  return (OK);
}

/* ---------------------------------------------------------------------- */
struct kinetics *
kinetics_replicate (struct kinetics *kinetics_old_ptr, int n_user_new)
/* ---------------------------------------------------------------------- */
{
  struct kinetics *kinetics_ptr;
  kinetics_ptr = kinetics_alloc ();
  kinetics_copy (kinetics_old_ptr, kinetics_ptr, n_user_new);
  return (kinetics_ptr);
}

/* ---------------------------------------------------------------------- */
struct kinetics *
kinetics_search (int n_user, int *n, int print)
/* ---------------------------------------------------------------------- */
{
/*
 *  Does linear search for user number n_user.
 *  Returns pointer to kinetics structure and position number, n,
 *     in kinetics array if found.
 *  Returns NULL if not found.
 */
  int i;
  for (i = 0; i < count_kinetics; i++)
  {
    if (n_user == kinetics[i].n_user)
    {
      break;
    }
  }
  if (i >= count_kinetics)
  {
    if (print == TRUE)
    {
      sprintf (error_string, "Kinetics %d not found.", n_user);
      error_msg (error_string, CONTINUE);
    }
    *n = -999;
    return (NULL);
  }
  *n = i;
  return (&kinetics[i]);
}

/* ---------------------------------------------------------------------- */
int
kinetics_sort (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Sort array of kinetics structures
 */
  if (count_kinetics > 0)
  {
    qsort (kinetics, (size_t) count_kinetics,
	   (size_t) sizeof (struct kinetics), kinetics_compare);
  }
  return (OK);
}

/* **********************************************************************
 *
 *   Routines related to structure "master"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
struct master *
master_alloc (void)
/* ---------------------------------------------------------------------- */
/*
 *   Allocates space to a master structure and initializes the space.
 *      arguments: void
 *      return: pointer to a master structure
 */
{
  struct master *ptr;
  ptr = (struct master *) PHRQ_malloc (sizeof (struct master));
  if (ptr == NULL)
    malloc_error ();
/*
 *   set pointers in structure to NULL
 */
  ptr->in = FALSE;
  ptr->number = -1;
  ptr->last_model = -1;
  ptr->type = 0;
  ptr->primary = FALSE;
  ptr->coef = 0.0;
  ptr->total = 0.0;
  ptr->isotope_ratio = 0;
  ptr->isotope_ratio_uncertainty = 0;
  ptr->isotope = 0;
  ptr->total_primary = 0;
  ptr->elt = NULL;
  ptr->alk = 0.0;
  ptr->gfw = 0.0;
  ptr->gfw_formula = NULL;
  ptr->unknown = NULL;
  ptr->s = NULL;
  ptr->rxn_primary = NULL;
  ptr->rxn_secondary = NULL;
  ptr->pe_rxn = NULL;
  ptr->minor_isotope = FALSE;
  return (ptr);
}

/* ---------------------------------------------------------------------- */
int
master_delete (char *ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Delete master species:  Free memory of master species structure, free
 *   the structure, and remove from array master.
 *
 *   Input
 *        ptr  character string with name of element or valence state
 *   Returns
 *        TRUE if master species was deleted.
 *        FALSE if master species was not found.
 */
  int j, n;

  if (master_search (ptr, &n) == NULL)
    return (FALSE);
  master_free (master[n]);
  for (j = n; j < (count_master - 1); j++)
  {
    master[j] = master[j + 1];
  }
  count_master--;
  return (TRUE);
}

/* ---------------------------------------------------------------------- */
int
master_free (struct master *master_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Free memory pointed to by master species pointer, master_ptr.
 *   Frees master_ptr itself.
 */
  if (master_ptr == NULL)
    return (ERROR);
  rxn_free (master_ptr->rxn_primary);
  rxn_free (master_ptr->rxn_secondary);
  master_ptr = (struct master *) free_check_null (master_ptr);
  return (OK);
}

/* ---------------------------------------------------------------------- */
struct master *
master_bsearch (const char *ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Uses binary search. Assumes master is in sort order.
 *   Find master species for string (*ptr) containing name of element or valence state.
 *
 *   Input: ptr    pointer to string containing element name
 *
 *   Return: pointer to master structure containing name ptr or NULL.
 */
  void *void_ptr;
  if (count_master == 0)
  {
    return (NULL);
  }
  void_ptr = bsearch ((const char *) ptr,
		      (char *) master,
		      (unsigned) count_master,
		      sizeof (struct master *), master_compare_string);
  if (void_ptr == NULL)
  {
    return (NULL);
  }
  else
  {
    return (*(struct master **) void_ptr);
  }
}

/* ---------------------------------------------------------------------- */
int
master_compare_string (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
  const char *string_ptr;
  const struct master *master_ptr;

  string_ptr = (const char *) ptr1;
  master_ptr = *(const struct master **) ptr2;
  return (strcmp_nocase (string_ptr, master_ptr->elt->name));
}

/* ---------------------------------------------------------------------- */
int
master_compare (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
  const struct master *master_ptr1, *master_ptr2;
  master_ptr1 = *(const struct master **) ptr1;
  master_ptr2 = *(const struct master **) ptr2;
  return (strcmp_nocase (master_ptr1->elt->name, master_ptr2->elt->name));
}

/* ---------------------------------------------------------------------- */
struct master *
master_bsearch_primary (char *ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Find primary master species for first element in the string, ptr.
 *   Uses binary search. Assumes master is in sort order.
 */
  int l;
  char *ptr1;
  char elt[MAX_LENGTH];
  struct master *master_ptr_primary;
/*
 *   Find element name
 */
  ptr1 = ptr;
  get_elt (&ptr1, elt, &l);
/*
 *   Search master species list
 */
  master_ptr_primary = master_bsearch (elt);
  if (master_ptr_primary == NULL)
  {
    input_error++;
    sprintf (error_string, "Could not find primary master species for %s.",
	     ptr);
    error_msg (error_string, CONTINUE);
  }
  return (master_ptr_primary);
}

/* ---------------------------------------------------------------------- */
struct master *
master_search (char *ptr, int *n)
/* ---------------------------------------------------------------------- */
{
/*
 *   Linear search of master to find master species in string, ptr.
 *   Returns pointer if found. n contains position in array master.
 *   Returns NULL if not found.
 */
  int i;
  struct master *master_ptr;
/*
 *   Search master species list
 */
  *n = -999;
  for (i = 0; i < count_master; i++)
  {
    if (strcmp (ptr, master[i]->elt->name) == 0)
    {
      *n = i;
      master_ptr = master[i];
      return (master_ptr);
    }
  }
  return (NULL);
}

/* **********************************************************************
 *
 *   Routines related to structure "mix"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
struct mix *
mix_bsearch (int k, int *n)
/* ---------------------------------------------------------------------- */
{
/*
 *   Search mix to find if user mix number k exists.
 *   Uses binary search, assumes array mix is in sort order.
 */
  void *void_ptr;
  if (count_mix == 0)
  {
    *n = -999;
    return (NULL);
  }
  void_ptr = (void *)
    bsearch ((char *) &k,
	     (char *) mix,
	     (size_t) count_mix,
	     (size_t) sizeof (struct mix), mix_compare_int);
  if (void_ptr == NULL)
  {
    *n = -999;
    return (NULL);
  }
  *n = (int) ((struct mix *) void_ptr - mix);
  return ((struct mix *) void_ptr);
}

/* ---------------------------------------------------------------------- */
int
mix_compare (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
  const struct mix *mix_ptr1, *mix_ptr2;
  mix_ptr1 = (const struct mix *) ptr1;
  mix_ptr2 = (const struct mix *) ptr2;
  if (mix_ptr1->n_user > mix_ptr2->n_user)
    return (1);
  if (mix_ptr1->n_user < mix_ptr2->n_user)
    return (-1);
  return (0);

}

/* ---------------------------------------------------------------------- */
static int
mix_compare_int (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
/*
 *   Compare mix user numbers
 */
  const int *nptr1;
  const struct mix *nptr2;

  nptr1 = (const int *) ptr1;
  nptr2 = (const struct mix *) ptr2;
  if (*nptr1 > nptr2->n_user)
    return (1);
  if (*nptr1 < nptr2->n_user)
    return (-1);
  return (0);
}

/* ---------------------------------------------------------------------- */
int
mix_copy (struct mix *mix_old_ptr, struct mix *mix_new_ptr, int n_user_new)
/* ---------------------------------------------------------------------- */
{
/*
 *   Copies mix data from mix_old_ptr to new location, mix_new_ptr.
 *   Space for the mix_new_ptr structure must already be malloced.
 *   Space for mix components is malloced here.
 */
  char token[MAX_LENGTH];
/*
 *   Copy old to new
 */
  memcpy (mix_new_ptr, mix_old_ptr, sizeof (struct mix));
/*
 *   Store data for structure mix
 */
  mix_new_ptr->n_user = n_user_new;
  mix_new_ptr->n_user_end = n_user_new;
  sprintf (token, "Mix defined in simulation %d.", simulation);
  mix_new_ptr->description = string_duplicate (token);
/*
 *   Count mix components and allocate space
 */
  mix_new_ptr->comps =
    (struct mix_comp *) PHRQ_malloc ((size_t) (mix_old_ptr->count_comps) *
				     sizeof (struct mix_comp));
  if (mix_new_ptr->comps == NULL)
    malloc_error ();
  memcpy (mix_new_ptr->comps, mix_old_ptr->comps,
	  (size_t) (mix_old_ptr->count_comps) * sizeof (struct mix_comp));

  return (OK);
}

/* ---------------------------------------------------------------------- */
int
mix_duplicate (int n_user_old, int n_user_new)
/* ---------------------------------------------------------------------- */
{
/*
 *   Checks if n_user_new exists, and frees space if it does
 *   Copies mix[n_user_old] to old n_user_new space if
 *   found or to mix[count_mix] if not found.
 *   Mix array may not be in sort order after the copy.
 */
  int n_old, n_new, sort;
  char token[MAX_LENGTH];
  int count_comps;
  struct mix *mix_ptr_old, *mix_ptr_new;
/*
 *   Find n_user_old in structure array mix or make new space
 */
  if (n_user_old == n_user_new)
    return (OK);
  mix_ptr_old = mix_bsearch (n_user_old, &n_old);
  if (mix_ptr_old == NULL)
  {
    sprintf (error_string, "Mix %d not found.", n_user_old);
    error_msg (error_string, CONTINUE);
    input_error++;
    return (ERROR);
  }
/*
 *   Find n_user_old in structure array mix or make new space
 */
  sort = FALSE;
  mix_ptr_new = mix_bsearch (n_user_new, &n_new);
  if (mix_ptr_new != NULL)
  {
    mix_free (&(mix[n_new]));
  }
  else
  {
    mix =
      (struct mix *) PHRQ_realloc (mix,
				   (size_t) (count_mix +
					     1) * sizeof (struct mix));
    if (mix == NULL)
      malloc_error ();
    if (n_user_new < mix[count_mix - 1].n_user)
      sort = TRUE;
    n_new = count_mix++;
  }
/*
 *   Store data for structure mix
 */
  memcpy (&mix[n_new], &mix[n_old], sizeof (struct mix));
  count_comps = mix[n_old].count_comps;
  mix[n_new].n_user = n_user_new;
  mix[n_new].n_user_end = n_user_new;
  sprintf (token, "Mix defined in simulation %d.", simulation);
  mix[n_new].description = string_duplicate (token);
/*
 *   Count mix components and allocate space
 */
  mix[n_new].comps =
    (struct mix_comp *) PHRQ_malloc ((size_t) (count_comps) *
				     sizeof (struct mix_comp));
  if (mix[n_new].comps == NULL)
    malloc_error ();
/*
 *   Write mix_comp structure for each mix component
 */
  memcpy (mix[n_new].comps, mix[n_old].comps,
	  (size_t) (count_comps) * sizeof (struct mix_comp));
  if (sort == TRUE)
    mix_sort ();
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
mix_free (struct mix *mix_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Free space allocated for mixture structure.
 */
  if (mix_ptr == NULL)
    return (ERROR);

  mix_ptr->description = (char *) free_check_null (mix_ptr->description);
  mix_ptr->comps = (struct mix_comp *) free_check_null (mix_ptr->comps);
  return (OK);
}

/* ---------------------------------------------------------------------- */
struct mix *
mix_search (int n_user, int *n, int print)
/* ---------------------------------------------------------------------- */
{
/*
 *   Linear search of mix array for user number n_user.
 *   Returns pointer to mix structure if found, n is set to position in mix array.
 *   Returns NULL if not found.
 */
  int i;
  for (i = 0; i < count_mix; i++)
  {
    if (n_user == mix[i].n_user)
    {
      break;
    }
  }
  if (i >= count_mix)
  {
    if (print == TRUE)
    {
      sprintf (error_string, "Mix %d not found.", n_user);
      error_msg (error_string, CONTINUE);
    }
    *n = -999;
    return (NULL);
  }
  *n = i;
  return (&mix[i]);
}

/* ---------------------------------------------------------------------- */
int
mix_sort (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Sort array of mix structures
 */
  if (count_mix > 0)
  {
    qsort (mix, (size_t) count_mix, (size_t) sizeof (struct mix),
	   mix_compare);
  }
  return (OK);
}

/* **********************************************************************
 *
 *   Routines related to structure "pe_data"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
struct pe_data *
pe_data_alloc (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Allocate and initialize an array of struct pe_data
 */
  struct pe_data *pe_data_ptr;
  char token[MAX_LENGTH];

  pe_data_ptr =
    (struct pe_data *) PHRQ_malloc ((size_t) (2 * sizeof (struct pe_data)));
  if (pe_data_ptr == NULL)
    malloc_error ();
  pe_data_ptr[0].name = pe_string;
  if (s_eminus != NULL && s_eminus->rxn != NULL)
  {
    pe_data_ptr[0].rxn = rxn_dup (s_eminus->rxn);
  }
  else
  {
    pe_data_ptr[0].rxn = rxn_alloc (3);
    if (pe_data_ptr[0].rxn == NULL)
      malloc_error ();
    strcpy (token, "e-");
    s_eminus = s_store (token, -1.0, FALSE);
    pe_data_ptr[0].rxn->token[0].s = s_eminus;
    pe_data_ptr[0].rxn->token[0].coef = -1.0;
    pe_data_ptr[0].rxn->token[1].s = s_eminus;
    pe_data_ptr[0].rxn->token[1].coef = -1.0;
    pe_data_ptr[0].rxn->token[2].s = NULL;
  }
  pe_data_ptr[1].name = NULL;
  pe_data_ptr[1].rxn = NULL;

  return (pe_data_ptr);
}

/* ---------------------------------------------------------------------- */
struct pe_data *
pe_data_dup (struct pe_data *pe_ptr_old)
/* ---------------------------------------------------------------------- */
{
/*
 *   Allocate space and copy data from pe_ptr_old into the new space.
 *   Returns pointer to new space.
 */
  int i, count_pe;
  struct pe_data *pe_ptr_new;
/*
 *   Count pe data and copy into new structure
 */
  if (pe_ptr_old == NULL)
    return (ERROR);
  for (i = 0; pe_ptr_old[i].name != NULL; i++);
  count_pe = i + 1;
  pe_ptr_new =
    (struct pe_data *) PHRQ_malloc ((size_t) count_pe *
				    sizeof (struct pe_data));
  if (pe_ptr_new == NULL)
    malloc_error ();
  memcpy ((void *) pe_ptr_new, (void *) pe_ptr_old,
	  (size_t) count_pe * sizeof (struct pe_data));
  for (i = 0; i < count_pe - 1; i++)
  {
    pe_ptr_new[i].rxn = rxn_dup (pe_ptr_old[i].rxn);
  }
  pe_ptr_new[count_pe - 1].rxn = NULL;
  return (pe_ptr_new);
}

/* ---------------------------------------------------------------------- */
struct pe_data *
pe_data_free (struct pe_data *pe_data_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Frees all data for pe_data_ptr, including pe_data_ptr.
 */
  int i;

  if (pe_data_ptr == NULL)
    return (ERROR);
  for (i = 0; pe_data_ptr[i].name != NULL; i++)
  {
    rxn_free (pe_data_ptr[i].rxn);
  }
  pe_data_ptr = (struct pe_data *) free_check_null (pe_data_ptr);
  return (NULL);
}

/* ---------------------------------------------------------------------- */
int
pe_data_store (struct pe_data **pe, const char *token)
/* ---------------------------------------------------------------------- */
{
/*
 *   Find pe name in current list of pe_data structures
 *   or allocate new pe_data structure
 */
  int i;
  struct pe_data *pe_data_ptr;
/*
 *   Search for pe in list
 */
  pe_data_ptr = *pe;
  for (i = 0; pe_data_ptr[i].name != NULL; i++)
  {
    if (strcmp (token, pe_data_ptr[i].name) == 0)
    {
      return (i);
    }
  }
/*
 *   Save new struct pe_data in array
 */
  *pe =
    (struct pe_data *) PHRQ_realloc ((void *) *pe,
				     (size_t) (i +
					       2) * sizeof (struct pe_data));
  if (*pe == NULL)
    malloc_error ();
  pe_data_ptr = *pe;
  pe_data_ptr[i].name = string_hsave (token);
  pe_data_ptr[i].rxn = NULL;
  pe_data_ptr[i + 1].name = NULL;
  return (i);
}

/* **********************************************************************
 *
 *   Routines related to structure "phases"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
struct phase *
phase_alloc (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Allocates space to a phase structure and initializes
 *      arguments: void
 *      return: pointer to new phase structure
 */
  struct phase *phase_ptr;
/*
 *   Allocate space
 */
  phase_ptr = (struct phase *) PHRQ_malloc (sizeof (struct phase));
  if (phase_ptr == NULL)
    malloc_error ();
/*
 *   Initialize space
 */
  phase_init (phase_ptr);
  return (phase_ptr);
}

/* ---------------------------------------------------------------------- */
int
phase_compare (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
/*
 *   Compares names of phases for sort
 */
  const struct phase *phase_ptr1, *phase_ptr2;
  phase_ptr1 = *(const struct phase **) ptr1;
  phase_ptr2 = *(const struct phase **) ptr2;
  return (strcmp_nocase (phase_ptr1->name, phase_ptr2->name));
}

/* ---------------------------------------------------------------------- */
int
phase_compare_string (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
  const char *char_ptr;
  const struct phase *phase_ptr;
  char_ptr = (const char *) ptr1;
  phase_ptr = *(const struct phase **) ptr2;
  return (strcmp_nocase (char_ptr, phase_ptr->name));
}

/* ---------------------------------------------------------------------- */
int
phase_delete (int i)
/* ---------------------------------------------------------------------- */
{
/*
 *   Deletes phase i from list, phases
 *   Frees memory allocated to phase[i] and renumbers phases to remove number i.
 *   Input: i, number of phase
 *   Return: OK
 */
  int j;

  phase_free (phases[i]);
  phases[i] = (struct phase *) free_check_null (phases[i]);
  for (j = i; j < (count_phases - 1); j++)
  {
    phases[j] = phases[j + 1];
  }
  count_phases--;
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
phase_free (struct phase *phase_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Frees memory allocated within phase[i], does not free phase structure
 *   Input: i, number of phase
 *   Return: OK
 */
  if (phase_ptr == NULL)
    return (ERROR);
  phase_ptr->next_elt =
    (struct elt_list *) free_check_null (phase_ptr->next_elt);
  phase_ptr->next_sys_total =
    (struct elt_list *) free_check_null (phase_ptr->next_sys_total);
  rxn_free (phase_ptr->rxn);
  rxn_free (phase_ptr->rxn_s);
  rxn_free (phase_ptr->rxn_x);
  phase_ptr->add_logk =
    (struct name_coef *) free_check_null (phase_ptr->add_logk);
  return (OK);
}

/* ---------------------------------------------------------------------- */
struct phase *
phase_bsearch (const char *ptr, int *j, int print)
/* ---------------------------------------------------------------------- */
{
/*   Binary search the structure array "phases" for a name that is equal to
 *   ptr. Assumes array phases is in sort order.
 *
 *   Arguments:
 *      name  input, a character string to be located in phases.
 *      j            index number in array phases.
 *
 *   Returns:
 *      if found, pointer to phase structure.
 *      if not found, NULL
 *
 */
  void *void_ptr;

  void_ptr = NULL;
  if (count_phases > 0)
  {
    void_ptr = (void *)
      bsearch ((char *) ptr,
	       (char *) phases,
	       (size_t) count_phases,
	       (size_t) sizeof (struct phase *), phase_compare_string);
  }
  if (void_ptr == NULL && print == TRUE)
  {
    sprintf (error_string, "Could not find phase in list, %s.", ptr);
    error_msg (error_string, CONTINUE);
  }

  if (void_ptr == NULL)
  {
    *j = -1;
    return (NULL);
  }

  *j = (int) ((struct phase **) void_ptr - phases);
  return (*(struct phase **) void_ptr);
}

/* ---------------------------------------------------------------------- */
static int
phase_init (struct phase *phase_ptr)
/* ---------------------------------------------------------------------- */
/*
 *   set pointers in phase structure to NULL
 */
{
  int i;

  phase_ptr->name = NULL;
  phase_ptr->formula = NULL;
  phase_ptr->in = FALSE;
  phase_ptr->lk = 0.0;
  for (i = 0; i < 8; i++)
    phase_ptr->logk[i] = 0.0;
  phase_ptr->original_units = kjoules;
  phase_ptr->type = SOLID;
  phase_ptr->check_equation = TRUE;
  phase_ptr->next_elt = NULL;
  phase_ptr->next_sys_total = NULL;
  phase_ptr->rxn = NULL;
  phase_ptr->rxn_s = NULL;
  phase_ptr->rxn_x = NULL;
  phase_ptr->add_logk = NULL;
  phase_ptr->count_add_logk = 0;
  return (OK);
}

/* ---------------------------------------------------------------------- */
struct phase *
phase_store (char *name)
/* ---------------------------------------------------------------------- */
{
/*
 *   Function locates the string "name" in the hash table for phases.
 *
 *   If found, pointer to the appropriate phase structure is returned.
 *
 *   If the string is not found, a new entry is made at the end of
 *   the phases array (position count_phases) and count_phases is
 *   incremented. A new entry is made in the hash table. Pointer to
 *   the new structure is returned.
 *
 *   Arguments:
 *      name    input, character string to be located or stored.
 *
 *   Returns:
 *      The address of a phase structure that contains the phase data.
 *      If phase existed, it is reinitialized. The structure returned
 *      contains only the name of the phase.
 */
  int n;
  struct phase *phase_ptr;
  ENTRY item, *found_item;
  char token[MAX_LENGTH];
  char *ptr;
/*
 *   Search list
 */

  strcpy (token, name);
  str_tolower (token);
  ptr = string_hsave (token);

  item.key = ptr;
  item.data = NULL;
  found_item = hsearch_multi (phases_hash_table, item, FIND);
  if (found_item != NULL)
  {
    phase_ptr = (struct phase *) (found_item->data);
    phase_free (phase_ptr);
    phase_init (phase_ptr);
    phase_ptr->name = string_hsave (name);
    return (phase_ptr);
  }
/*
 *   Make new phase structure and return pointer to it
 */
  /* make sure there is space in phases */
  n = count_phases++;
  if (count_phases >= max_phases)
  {
    space ((void **) ((void *) &phases), count_phases, &max_phases,
	   sizeof (struct phase *));
  }
  phases[n] = phase_alloc ();
  /* set name in phase structure */
  phases[n]->name = string_hsave (name);
/*
 *   Update hash table
 */
  item.key = ptr;
  item.data = (void *) phases[n];
  found_item = hsearch_multi (phases_hash_table, item, ENTER);
  if (found_item == NULL)
  {
    sprintf (error_string, "Hash table error in phase_store.");
    error_msg (error_string, CONTINUE);
  }

  return (phases[n]);
}

#ifdef SKIP
/* ---------------------------------------------------------------------- */
struct phase *
phase_search (char *ptr, int *n)
/* ---------------------------------------------------------------------- */
{
/*
 *   Linear search of phases to find phase with name in ptr.
 *   Returns pointer if found. n contains position in array master.
 *   Returns NULL if not found.
 */
  int i;
  struct phase *phase_ptr;
/*
 *   Search phases list
 */
  *n = -999;
  for (i = 0; i < count_phases; i++)
  {
    if (strcmp_nocase (ptr, phases[i]->name) == 0)
    {
      *n = i;
      phase_ptr = phases[i];
      return (phase_ptr);
    }
  }
  return (NULL);
}

/* ---------------------------------------------------------------------- */
struct phase *
phase_store (char *name)
/* ---------------------------------------------------------------------- */
{
/*
 *   Function locates the string "name" in the hash table for phases.
 *
 *   If found, pointer to the appropriate phase structure is returned.
 *
 *   If the string is not found, a new entry is made at the end of
 *   the phases array (position count_phases) and count_phases is
 *   incremented. A new entry is made in the hash table. Pointer to
 *   the new structure is returned.
 *
 *   Arguments:
 *      name    input, character string to be located or stored.
 *
 *   Returns:
 *      The address of a phase structure that contains the phase data.
 *      If phase existed, it is reinitialized. The structure returned
 *      contains only the name of the phase.
 */
  int n;
  struct phase *phase_ptr;
/*
 *   Search list
 */
  phase_ptr = phase_search (name, &n);
  if (phase_ptr != NULL)
  {
    phase_free (phase_ptr);
    phase_init (phase_ptr);
    phase_ptr->name = string_hsave (name);
    return (phase_ptr);
  }
/*
 *   Make new phase structure and return pointer to it
 */
  /* make sure there is space in phases */
  n = count_phases++;
  if (count_phases >= max_phases)
  {
    space ((void **) &phases, count_phases, &max_phases,
	   sizeof (struct phase *));
  }
  phases[n] = phase_alloc ();
  /* set name in phase structure */
  phases[n]->name = string_hsave (name);
  return (phases[n]);
}
#endif
/* **********************************************************************
 *
 *   Routines related to structure "pp_assemblage"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
struct pp_assemblage *
pp_assemblage_alloc (void)
/* ---------------------------------------------------------------------- */
{
  struct pp_assemblage *pp_assemblage_ptr;
  pp_assemblage_ptr =
    (struct pp_assemblage *) PHRQ_malloc (sizeof (struct pp_assemblage));
  if (pp_assemblage_ptr == NULL)
    malloc_error ();
  pp_assemblage_ptr->n_user = -1;
  pp_assemblage_ptr->n_user_end = -1;
  pp_assemblage_ptr->description = NULL;
  pp_assemblage_ptr->new_def = 0;
  pp_assemblage_ptr->next_elt = NULL;
  pp_assemblage_ptr->count_comps = 0;
  pp_assemblage_ptr->pure_phases = NULL;
  return (pp_assemblage_ptr);
}

/* ---------------------------------------------------------------------- */
struct pp_assemblage *
pp_assemblage_bsearch (int k, int *n)
/* ---------------------------------------------------------------------- */
{
/*
 *   Binary search of array pp_assemblage to find user number k.
 *   Assumes pp_assemblage is in sort order.
 *
 *   Input: k, user number to find.
 *
 *   Output: n, position in array pp_assemblage for user number k.
 *
 *   Return: pointer to pp_assemblage structure for user number k, if found.
 *           NULL if not found.
 */
  void *void_ptr;
  if (count_pp_assemblage == 0)
  {
    *n = -999;
    return (NULL);
  }
  void_ptr = (void *)
    bsearch ((char *) &k,
	     (char *) pp_assemblage,
	     (size_t) count_pp_assemblage,
	     (size_t) sizeof (struct pp_assemblage),
	     pp_assemblage_compare_int);
  if (void_ptr == NULL)
  {
    *n = -999;
    return (NULL);
  }
  *n = (int) ((struct pp_assemblage *) void_ptr - pp_assemblage);
  return ((struct pp_assemblage *) void_ptr);
}

/* ---------------------------------------------------------------------- */
int
pp_assemblage_compare (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
  const struct pp_assemblage *pp_assemblage_ptr1, *pp_assemblage_ptr2;
  pp_assemblage_ptr1 = (const struct pp_assemblage *) ptr1;
  pp_assemblage_ptr2 = (const struct pp_assemblage *) ptr2;
  if (pp_assemblage_ptr1->n_user > pp_assemblage_ptr2->n_user)
    return (1);
  if (pp_assemblage_ptr1->n_user < pp_assemblage_ptr2->n_user)
    return (-1);
  return (0);

}

/* ---------------------------------------------------------------------- */
static int
pp_assemblage_compare_int (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
/*
 *   Compare pp_assemblage user numbers
 */
  const int *nptr1;
  const struct pp_assemblage *nptr2;

  nptr1 = (const int *) ptr1;
  nptr2 = (const struct pp_assemblage *) ptr2;
  if (*nptr1 > nptr2->n_user)
    return (1);
  if (*nptr1 < nptr2->n_user)
    return (-1);
  return (0);
}

/* ---------------------------------------------------------------------- */
int
pp_assemblage_copy (struct pp_assemblage *pp_assemblage_old_ptr,
		    struct pp_assemblage *pp_assemblage_new_ptr,
		    int n_user_new)
/* ---------------------------------------------------------------------- */
{
/*
 *   Copies pp_assemblage data from pp_assemblage_old_ptr to pp_assemblage_new_ptr.
 *   Space for a pp_assemblage_new_ptr must already be malloced.
 *   Space for list of pure phases in assemblage is malloced here.
 */
  int count_comps;
  char token[MAX_LENGTH];
/*
 *   Store data for structure pp_assemblage
 */
  pp_assemblage_new_ptr->n_user = n_user_new;
  pp_assemblage_new_ptr->n_user_end = n_user_new;
  pp_assemblage_new_ptr->new_def = pp_assemblage_old_ptr->new_def;
  sprintf (token, "Pp_Assemblage defined in simulation %d.", simulation);
  pp_assemblage_new_ptr->description = string_duplicate (token);
  pp_assemblage_new_ptr->next_elt =
    elt_list_dup (pp_assemblage_old_ptr->next_elt);
/*
 *   Count pure phases
 */
  count_comps = pp_assemblage_old_ptr->count_comps;
  pp_assemblage_new_ptr->count_comps = count_comps;
/*
 *   Malloc space and copy
 */
  pp_assemblage_new_ptr->pure_phases =
    (struct pure_phase *) PHRQ_malloc ((size_t) count_comps *
				       sizeof (struct pure_phase));
  if (pp_assemblage_new_ptr->pure_phases == NULL)
    malloc_error ();

  memcpy ((void *) pp_assemblage_new_ptr->pure_phases,
	  (void *) pp_assemblage_old_ptr->pure_phases,
	  (size_t) count_comps * sizeof (struct pure_phase));
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
pp_assemblage_copy_to_last (int n, int n_user)
/* ---------------------------------------------------------------------- */
{
/*
 *   Copies an pp_assemblage definition from position n
 *   to position count_pp_assemblage.
 */
  space ((void **) ((void *) &pp_assemblage), count_pp_assemblage,
	 &max_pp_assemblage, sizeof (struct pp_assemblage));
  pp_assemblage_copy (&pp_assemblage[n], &pp_assemblage[count_pp_assemblage],
		      n_user);
  count_pp_assemblage++;
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
pp_assemblage_duplicate (int n_user_old, int n_user_new)
/* ---------------------------------------------------------------------- */
{
/*
 *   Checks if n_user_new exists, and frees space if it does
 *   Copies pp_assemblage[n_user_old] to old n_user_new space if
 *   found or to pp_assemblage[count_pp_assemblage] if not found.
 *   Pp_Assemblage array may not be in sort order after the copy.
 */
  int n_old, n_new, sort;
  struct pp_assemblage *pp_assemblage_ptr_old, *pp_assemblage_ptr_new;
/*
 *   Find n_user_old in structure array pp_assemblage
 */
  if (n_user_old == n_user_new)
    return (OK);
  pp_assemblage_ptr_old = pp_assemblage_bsearch (n_user_old, &n_old);
  if (pp_assemblage_ptr_old == NULL)
  {
    sprintf (error_string, "Pp_Assemblage %d not found.", n_user_old);
    error_msg (error_string, CONTINUE);
    input_error++;
    return (ERROR);
  }
/*
 *   Find n_user_new in structure array pp_assemblage or make new space
 */
  sort = FALSE;
  pp_assemblage_ptr_new = pp_assemblage_bsearch (n_user_new, &n_new);
  if (pp_assemblage_ptr_new != NULL)
  {
    pp_assemblage_free (pp_assemblage_ptr_new);
  }
  else
  {
    space ((void **) ((void *) &pp_assemblage), count_pp_assemblage,
	   &max_pp_assemblage, sizeof (struct pp_assemblage));
    if (n_user_new < pp_assemblage[count_pp_assemblage - 1].n_user)
      sort = TRUE;
    n_new = count_pp_assemblage++;
  }
/*
 *   Copy data
 */
  pp_assemblage_copy (&pp_assemblage[n_old], &pp_assemblage[n_new],
		      n_user_new);
  if (sort == TRUE)
    pp_assemblage_sort ();
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
pp_assemblage_delete (int n_user_old)
/* ---------------------------------------------------------------------- */
/*
 *   Frees space for user number n_user_old, removes structure from
 *   array pp_assemblage.
 */
{
  int i;
  int n_old;
  struct pp_assemblage *pp_assemblage_ptr_old;
/*
 *   Find n_user_old in structure array
 */
  pp_assemblage_ptr_old = pp_assemblage_bsearch (n_user_old, &n_old);
  if (pp_assemblage_ptr_old != NULL)
  {
    /*
     *   Delete pp_assemblage
     */
    pp_assemblage_free (&pp_assemblage[n_old]);

    for (i = n_old + 1; i < count_pp_assemblage; i++)
    {
      memcpy ((void *) &pp_assemblage[i - 1], (void *) &pp_assemblage[i],
	      (size_t) sizeof (struct pp_assemblage));
    }
    count_pp_assemblage--;
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
pp_assemblage_free (struct pp_assemblage *pp_assemblage_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Free space allocated for pp_assemblage structure
 */
  if (pp_assemblage_ptr == NULL)
    return (ERROR);

  pp_assemblage_ptr->description =
    (char *) free_check_null (pp_assemblage_ptr->description);
  pp_assemblage_ptr->next_elt =
    (struct elt_list *) free_check_null (pp_assemblage_ptr->next_elt);
  pp_assemblage_ptr->pure_phases =
    (struct pure_phase *) free_check_null (pp_assemblage_ptr->pure_phases);
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
pp_assemblage_init (struct pp_assemblage *pp_assemblage_ptr, int n_user,
		    int n_user_end, char *description)
/* ---------------------------------------------------------------------- */
{
/*
 *   Provides initial values for a new pp_assemblage
 *      arguments:
 *            n position in array pp_assemblage
 *            n_user  user number for pp_assemblage
 *      return: OK
 */
  pp_assemblage_ptr->n_user = n_user;
  pp_assemblage_ptr->n_user_end = n_user_end;
  pp_assemblage_ptr->new_def = TRUE;
  pp_assemblage_ptr->description = string_duplicate (description);

  pp_assemblage_ptr->next_elt = NULL;
  pp_assemblage_ptr->count_comps = 0;
  pp_assemblage_ptr->pure_phases =
    (struct pure_phase *) PHRQ_malloc ((size_t) sizeof (struct pure_phase));
  if (pp_assemblage_ptr->pure_phases == NULL)
    malloc_error ();

  return (OK);
}

/* ---------------------------------------------------------------------- */
int
pp_assemblage_ptr_to_user (struct pp_assemblage *pp_assemblage_ptr_old,
			   int n_user_new)
/* ---------------------------------------------------------------------- */
{
/*
 *   Checks if n_user_new exists, and frees space if it does
 *   Copies pp_assemblage_ptr_old to old n_user_new space if
 *   found or to pp_assemblage[count_pp_assemblage] if not found.
 *   Pp_Assemblage array may not be in sort order after the copy.
 */
  int n_new, sort;
  struct pp_assemblage *pp_assemblage_ptr_new;
/*
 *   Find n_user_old in structure array pp_assemblage
 */
  if (pp_assemblage_ptr_old == NULL)
  {
    sprintf (error_string, "Pp_Assemblage pointer is NULL.");
    error_msg (error_string, CONTINUE);
    input_error++;
    return (ERROR);
  }
/*
 *   Find n_user_new in structure array pp_assemblage or make new space
 */
  sort = FALSE;
  pp_assemblage_ptr_new = pp_assemblage_bsearch (n_user_new, &n_new);
  if (pp_assemblage_ptr_new == pp_assemblage_ptr_old)
    return (OK);
  if (pp_assemblage_ptr_new != NULL)
  {
    pp_assemblage_free (pp_assemblage_ptr_new);
  }
  else
  {
    space ((void **) ((void *) &pp_assemblage), count_pp_assemblage,
	   &max_pp_assemblage, sizeof (struct pp_assemblage));
    if (n_user_new < pp_assemblage[count_pp_assemblage - 1].n_user)
      sort = TRUE;
    n_new = count_pp_assemblage++;
  }
/*
 *   Copy data
 */
  pp_assemblage_copy (pp_assemblage_ptr_old, &pp_assemblage[n_new],
		      n_user_new);
  if (sort == TRUE)
    pp_assemblage_sort ();
  return (OK);
}

/* ---------------------------------------------------------------------- */
struct pp_assemblage *
pp_assemblage_replicate (struct pp_assemblage *pp_assemblage_old_ptr,
			 int n_user_new)
/* ---------------------------------------------------------------------- */
{
  struct pp_assemblage *pp_assemblage_ptr;
  pp_assemblage_ptr = pp_assemblage_alloc ();
  pp_assemblage_copy (pp_assemblage_old_ptr, pp_assemblage_ptr, n_user_new);
  return (pp_assemblage_ptr);
}

/* ---------------------------------------------------------------------- */
struct pp_assemblage *
pp_assemblage_search (int n_user, int *n)
/* ---------------------------------------------------------------------- */
{
/*   Linear search of the structure array "pp_assemblage" for user number n_user.
 *
 *   Arguments:
 *      n_user  input, user number
 *      n       output, position in pp_assemblage
 *
 *   Returns:
 *      if found, the address of the pp_assemblage element
 *      if not found, NULL
 */
  int i;
  for (i = 0; i < count_pp_assemblage; i++)
  {
    if (pp_assemblage[i].n_user == n_user)
    {
      *n = i;
      return (&(pp_assemblage[i]));
    }
  }
/*
 *   Pp_Assemblage not found
 */
  return (NULL);
}

/* ---------------------------------------------------------------------- */
int
pp_assemblage_sort (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Sort array of pp_assemblage structures
 */
  if (count_pp_assemblage > 0)
  {
    qsort (pp_assemblage, (size_t) count_pp_assemblage,
	   (size_t) sizeof (struct pp_assemblage), pp_assemblage_compare);
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
pure_phase_compare (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
  const struct pure_phase *pure_phase_ptr1, *pure_phase_ptr2;
  pure_phase_ptr1 = (const struct pure_phase *) ptr1;
  pure_phase_ptr2 = (const struct pure_phase *) ptr2;
  return (strcmp_nocase (pure_phase_ptr1->name, pure_phase_ptr2->name));

}

/* **********************************************************************
 *
 *   Routines related to structure "rates"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
struct rate *
rate_bsearch (char *ptr, int *j)
/* ---------------------------------------------------------------------- */
{
/*   Binary search the structure array "rates" for a name that is equal to
 *   ptr. Assumes array rates is in sort order.
 *
 *   Arguments:
 *      name  input, a character string to be located in rates.
 *      j            index number in array rates.
 *
 *   Returns:
 *      if found, pointer to rate structure.
 *      if not found, NULL
 *
 */
  void *void_ptr;

  if (count_rates == 0)
  {
    *j = -1;
    return (NULL);
  }
  void_ptr = (void *)
    bsearch ((char *) ptr,
	     (char *) rates,
	     (size_t) count_rates,
	     (size_t) sizeof (struct rate *), rate_compare_string);

  if (void_ptr == NULL)
  {
    *j = -1;
    return (NULL);
  }

  *j = (int) ((struct rate *) void_ptr - rates);
  return ((struct rate *) void_ptr);
}

/* ---------------------------------------------------------------------- */
int
rate_compare (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
/*
 *   Compares names of rates for sort
 */
  const struct rate *rate_ptr1, *rate_ptr2;
  rate_ptr1 = *(const struct rate **) ptr1;
  rate_ptr2 = *(const struct rate **) ptr2;
  return (strcmp_nocase (rate_ptr1->name, rate_ptr2->name));
}

/* ---------------------------------------------------------------------- */
int
rate_compare_string (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
  const char *char_ptr;
  const struct rate *rate_ptr;
  char_ptr = (const char *) ptr1;
  rate_ptr = *(const struct rate **) ptr2;
  return (strcmp_nocase (char_ptr, rate_ptr->name));
}

/* ---------------------------------------------------------------------- */
int
rate_free (struct rate *rate_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Frees memory allocated within rate[i], does not free rate structure
 *   Input: i, number of rate
 *   Return: OK
 */
  char cmd[] = "new; quit";

  if (rate_ptr == NULL)
    return (ERROR);
  rate_ptr->commands = (char *) free_check_null (rate_ptr->commands);
  basic_run (cmd, rate_ptr->linebase, rate_ptr->varbase, rate_ptr->loopbase);
  rate_ptr->linebase = NULL;
  rate_ptr->varbase = NULL;
  rate_ptr->loopbase = NULL;
  return (OK);
}

/* ---------------------------------------------------------------------- */
struct rate *
rate_search (char *name, int *n)
/* ---------------------------------------------------------------------- */
{
/*   Linear search of the structure array "rates" for name.
 *
 *   Arguments:
 *     name     input, name of rate
 *      n       output, position in rates
 *
 *   Returns:
 *      if found, the address of the pp_assemblage element
 *      if not found, NULL
 */
  int i;
  *n = -1;
  for (i = 0; i < count_rates; i++)
  {
    if (strcmp_nocase (rates[i].name, name) == 0)
    {
      *n = i;
      return (&(rates[i]));
    }
  }
/*
 *   rate name not found
 */
  return (NULL);
}

/* ---------------------------------------------------------------------- */
int
rate_sort (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Sort array of rate structures
 */
  if (count_rates > 0)
  {
    qsort (rates, (size_t) count_rates, (size_t) sizeof (struct rate),
	   rate_compare);
  }
  return (OK);
}

/* **********************************************************************
 *
 *   Routines related to structure "reaction", balanced chemical reactions
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
struct reaction *
rxn_alloc (int ntokens)
/* ---------------------------------------------------------------------- */
{
  int i;
/*
 *   Allocates space to a rxn structure
 *      input: ntokens, number of tokens in reaction
 *      return: pointer to a species structure
 */
  struct reaction *rxn_ptr;
/*
 *   Malloc reaction structure
 */
  rxn_ptr = (struct reaction *) PHRQ_malloc (sizeof (struct reaction));
  if (rxn_ptr == NULL)
    malloc_error ();
/*
 *   zero log k data
 */
  for (i = 0; i < 7; i++)
  {
    rxn_ptr->logk[i] = 0.0;
  }
/*
 *   Malloc rxn_token structure
 */
  rxn_ptr->token =
    (struct rxn_token *) PHRQ_malloc ((size_t) ntokens *
				      sizeof (struct rxn_token));
  for (i = 0; i < ntokens; i++)
  {
    rxn_ptr->token[i].s = NULL;
    rxn_ptr->token[i].name = NULL;
    rxn_ptr->token[i].coef = 0.0;
  }
  if (rxn_ptr->token == NULL)
    malloc_error ();
  return (rxn_ptr);
}

/* ---------------------------------------------------------------------- */
struct reaction *
rxn_dup (struct reaction *rxn_ptr_old)
/* ---------------------------------------------------------------------- */
{
/*
 *   mallocs space for a reaction and copies the reaction
 *   input: rxn_ptr_old, pointer to a reaction structure to copy
 *
 *   Return: rxn_ptr_new,  pointer to duplicated structure to copy
 */
  int i;
  struct reaction *rxn_ptr_new;

  if (rxn_ptr_old == NULL)
    return (NULL);
  for (i = 0; rxn_ptr_old->token[i].s != NULL; i++);

  rxn_ptr_new = rxn_alloc (i + 1);
/*
 *   Copy logk data
 */
  memcpy (rxn_ptr_new->logk, rxn_ptr_old->logk, (size_t) 7 * sizeof (LDBLE));
/*
 *   Copy tokens
 */
  memcpy (rxn_ptr_new->token, rxn_ptr_old->token,
	  (size_t) (i + 1) * sizeof (struct rxn_token));

  return (rxn_ptr_new);
}

/* ---------------------------------------------------------------------- */
LDBLE
rxn_find_coef (struct reaction * r_ptr, const char *str)
/* ---------------------------------------------------------------------- */
{
/*
 *   Finds coefficient of token in reaction.
 *   input: r_ptr, pointer to a reaction structure
 *          str, string to find as reaction token
 *
 *   Return: 0.0, if token not found
 *           coefficient of token, if found.
 */
  struct rxn_token *r_token;
  LDBLE coef;

  r_token = r_ptr->token + 1;
  coef = 0.0;
  while (r_token->s != NULL)
  {
    if (strcmp (r_token->s->name, str) == 0)
    {
      coef = r_token->coef;
      break;
    }
    r_token++;
  }
  return (coef);
}

/* ---------------------------------------------------------------------- */
int
rxn_free (struct reaction *rxn_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Frees space allocated for a reaction structure
 *      input: rxn_ptr, pointer to reaction structure
 *      return: ERROR, if pointer is NULL
 *              OK, otherwise.
 */
  if (rxn_ptr == NULL)
    return (ERROR);
  rxn_ptr->token = (struct rxn_token *) free_check_null (rxn_ptr->token);
  rxn_ptr = (struct reaction *) free_check_null (rxn_ptr);
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
rxn_print (struct reaction *rxn_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Frees space allocated for a reaction structure
 *      input: rxn_ptr, pointer to reaction structure
 *      return: ERROR, if pointer is NULL
 *              OK, otherwise.
 */
  struct rxn_token *next_token;
  int i;
  if (rxn_ptr == NULL)
    return (ERROR);
  next_token = rxn_ptr->token;
  output_msg (OUTPUT_MESSAGE, "log k data:\n");
  for (i = 0; i < 7; i++)
  {
    output_msg (OUTPUT_MESSAGE, "\t%f\n", (double) rxn_ptr->logk[i]);
  }
  output_msg (OUTPUT_MESSAGE, "Reaction definition\n");
  while (next_token->s != NULL || next_token->name != NULL)
  {
    output_msg (OUTPUT_MESSAGE, "\tcoef %f ", next_token->coef);
    if (next_token->s != NULL)
    {
      output_msg (OUTPUT_MESSAGE, "\tspecies token: %s ",
		  next_token->s->name);
    }
    if (next_token->name != NULL)
    {
      output_msg (OUTPUT_MESSAGE, "\tname token: %s", next_token->name);
    }
    output_msg (OUTPUT_MESSAGE, "\n");
    next_token++;
  }
  return (OK);
}

/* **********************************************************************
 *
 *   Routines related to structure "species"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
struct species *
s_alloc (void)
/* ---------------------------------------------------------------------- */
/*
 *   Allocates space to a species structure, initializes
 *      arguments: void
 *      return: pointer to a species structure
 */
{
  struct species *s_ptr;
  s_ptr = (struct species *) PHRQ_malloc (sizeof (struct species));
  if (s_ptr == NULL)
    malloc_error ();
/*
 *   set pointers in structure to NULL, variables to zero
 */
  s_init (s_ptr);

  return (s_ptr);
}

/* ---------------------------------------------------------------------- */
int
s_compare (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
  const struct species *s_ptr1, *s_ptr2;
  s_ptr1 = *(const struct species **) ptr1;
  s_ptr2 = *(const struct species **) ptr2;
  return (strcmp (s_ptr1->name, s_ptr2->name));

}

/* ---------------------------------------------------------------------- */
int
s_delete (int i)
/* ---------------------------------------------------------------------- */
{
/*
 *   Delete species i: free memory and renumber array of pointers, s.
 */
  int j;

  s_free (s[i]);
  s[i] = (struct species *) free_check_null (s[i]);
  for (j = i; j < (count_s - 1); j++)
  {
    s[j] = s[j + 1];
  }
  count_s--;
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
s_free (struct species *s_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Free space allocated for species structure, s_ptr. Does not free s_ptr.
 */
  if (s_ptr == NULL)
    return (ERROR);
  s_ptr->next_elt = (struct elt_list *) free_check_null (s_ptr->next_elt);
  s_ptr->next_secondary =
    (struct elt_list *) free_check_null (s_ptr->next_secondary);
  s_ptr->next_sys_total =
    (struct elt_list *) free_check_null (s_ptr->next_sys_total);
  s_ptr->add_logk = (struct name_coef *) free_check_null (s_ptr->add_logk);
  rxn_free (s_ptr->rxn);
  rxn_free (s_ptr->rxn_s);
  rxn_free (s_ptr->rxn_x);
  s_ptr->diff_layer =
    (struct species_diff_layer *) free_check_null (s_ptr->diff_layer);
  return (OK);
}

/* ---------------------------------------------------------------------- */
static int
s_init (struct species *s_ptr)
/* ---------------------------------------------------------------------- */
/*
 *      return: pointer to a species structure
 */
{
  int i;
/*
 *   set pointers in structure to NULL
 */
  s_ptr->name = NULL;
  s_ptr->mole_balance = NULL;
  s_ptr->primary = NULL;
  s_ptr->secondary = NULL;
  s_ptr->next_elt = NULL;
  s_ptr->next_secondary = NULL;
  s_ptr->next_sys_total = NULL;
  s_ptr->rxn = NULL;
  s_ptr->rxn_s = NULL;
  s_ptr->rxn_x = NULL;
  s_ptr->diff_layer = NULL;
/*
 *   set varibles = 0
 */
  s_ptr->in = FALSE;
  s_ptr->gfw = 0.0;
  s_ptr->z = 0.0;
  s_ptr->dw = 0.0;
  s_ptr->alk = 0.0;
  s_ptr->carbon = 0.0;
  s_ptr->co2 = 0.0;
  s_ptr->h = 0.0;
  s_ptr->o = 0.0;
  s_ptr->dha = 0.0;
  s_ptr->dhb = 0.0;
  s_ptr->lk = 0.0;
  for (i = 0; i < 8; i++)
  {
    s_ptr->logk[i] = 0.0;
  }
  s_ptr->original_units = kjoules;
  s_ptr->count_add_logk = 0;
  s_ptr->add_logk = NULL;
  s_ptr->lg = 0.0;
  s_ptr->lm = 0.0;
  s_ptr->la = 0.0;
  s_ptr->dg = 0.0;
  s_ptr->moles = 0.0;
  s_ptr->type = 0;
  s_ptr->gflag = 0;
  s_ptr->check_equation = TRUE;

  for (i = 0; i < 5; i++)
  {
    s_ptr->cd_music[i] = 0.0;
  }
  for (i = 0; i < 3; i++)
  {
    s_ptr->dz[i] = 0.0;
  }


  return (OK);
}

/* ---------------------------------------------------------------------- */
struct species *
s_search (char *name)
/* ---------------------------------------------------------------------- */
{
/*
 *   Function locates the string "name" in the hash table for species.
 *
 *   Arguments:
 *      name  input, a character string to be located in species.
 *      i    is obsolete.
 *
 *   Returns:
 *   If found, pointer to the appropriate species structure is returned.
 *       else, NULL pointer is returned.
 */
  struct species *s_ptr;
  ENTRY item, *found_item;

  item.key = name;
  item.data = NULL;
  found_item = hsearch_multi (species_hash_table, item, FIND);
  if (found_item != NULL)
  {
    s_ptr = (struct species *) (found_item->data);
    return (s_ptr);
  }
  return (NULL);
}

/* ---------------------------------------------------------------------- */
struct species *
s_store (char *name, LDBLE z, int replace_if_found)
/* ---------------------------------------------------------------------- */
{
/*
 *   Function locates the string "name" in the hash table for species.
 *
 *   Pointer to a species structure is always returned.
 *
 *   If the string is not found, a new entry is made at the end of
 *      the elements array (position count_elements) and count_elements is
 *      incremented. A new entry is made in the hash table. Pointer to
 *      the new structure is returned.
 *   If "name" is found and replace is true, pointers in old species structure
 *      are freed and replaced with additional input.
 *   If "name" is found and replace is false, the old species structure is not
 *      modified and a pointer to it is returned.
 *
 *   Arguments:
 *      name    input, character string to be found in "species".
 *      z       input, charge on "name"
 *      replace_if_found input, TRUE means reinitialize species if found
 *                     FALSE means just return pointer if found.
 *
 *   Returns:
 *      pointer to species structure "s" where "name" can be found.
 */
  int n;
  struct species *s_ptr;
  ENTRY item, *found_item;
/*
 *   Search list
 */
  item.key = name;
  item.data = NULL;
  found_item = hsearch_multi (species_hash_table, item, FIND);

  if (found_item != NULL && replace_if_found == FALSE)
  {
    s_ptr = (struct species *) (found_item->data);
    return (s_ptr);
  }
  else if (found_item != NULL && replace_if_found == TRUE)
  {
    s_ptr = (struct species *) (found_item->data);
    s_free (s_ptr);
    s_init (s_ptr);
  }
  else
  {
    n = count_s++;
    /* make sure there is space in s */
    if (count_s >= max_s)
    {
      space ((void **) ((void *) &s), count_s, &max_s,
	     sizeof (struct species *));
    }
    /* Make new species structure */
    s[n] = s_alloc ();
    s_ptr = s[n];
  }
  /* set name and z in pointer in species structure */
  s_ptr->name = string_hsave (name);
  s_ptr->z = z;
/*
 *   Update hash table
 */
  item.key = s_ptr->name;
  item.data = (void *) s_ptr;
  found_item = hsearch_multi (species_hash_table, item, ENTER);
  if (found_item == NULL)
  {
    sprintf (error_string, "Hash table error in species_store.");
    error_msg (error_string, CONTINUE);
  }

  return (s_ptr);
}

/* **********************************************************************
 *
 *   Routines related to structure "s_s_assemblage"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
struct s_s_assemblage *
s_s_assemblage_alloc (void)
/* ---------------------------------------------------------------------- */
{
  struct s_s_assemblage *s_s_assemblage_ptr;
  s_s_assemblage_ptr =
    (struct s_s_assemblage *) PHRQ_malloc (sizeof (struct s_s_assemblage));
  if (s_s_assemblage_ptr == NULL)
    malloc_error ();
  s_s_assemblage_ptr->n_user = -1;
  s_s_assemblage_ptr->n_user_end = -1;
  s_s_assemblage_ptr->description = NULL;
  s_s_assemblage_ptr->new_def = 0;
  s_s_assemblage_ptr->count_s_s = 0;
  s_s_assemblage_ptr->s_s = NULL;
  return (s_s_assemblage_ptr);
}

/* ---------------------------------------------------------------------- */
struct s_s_assemblage *
s_s_assemblage_bsearch (int k, int *n)
/* ---------------------------------------------------------------------- */
{
/*
 *   Binary search of array s_s_assemblage to find user number k.
 *   Assumes s_s_assemblage is in sort order.
 *
 *   Input: k, user number to find.
 *
 *   Output: n, position in array s_s_assemblage for user number k.
 *
 *   Return: pointer to s_s_assemblage structure for user number k, if found.
 *           NULL if not found.
 */
  void *void_ptr;
  if (count_s_s_assemblage == 0)
  {
    *n = -999;
    return (NULL);
  }
  void_ptr = (void *)
    bsearch ((char *) &k,
	     (char *) s_s_assemblage,
	     (size_t) count_s_s_assemblage,
	     (size_t) sizeof (struct s_s_assemblage),
	     s_s_assemblage_compare_int);
  if (void_ptr == NULL)
  {
    *n = -999;
    return (NULL);
  }
  *n = (int) ((struct s_s_assemblage *) void_ptr - s_s_assemblage);
  return ((struct s_s_assemblage *) void_ptr);
}

/* ---------------------------------------------------------------------- */
int
s_s_assemblage_compare (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
  const struct s_s_assemblage *s_s_assemblage_ptr1, *s_s_assemblage_ptr2;
  s_s_assemblage_ptr1 = (const struct s_s_assemblage *) ptr1;
  s_s_assemblage_ptr2 = (const struct s_s_assemblage *) ptr2;
  if (s_s_assemblage_ptr1->n_user > s_s_assemblage_ptr2->n_user)
    return (1);
  if (s_s_assemblage_ptr1->n_user < s_s_assemblage_ptr2->n_user)
    return (-1);
  return (0);

}

/* ---------------------------------------------------------------------- */
static int
s_s_assemblage_compare_int (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
/*
 *   Compare s_s_assemblage user numbers
 */
  const int *nptr1;
  const struct s_s_assemblage *nptr2;

  nptr1 = (const int *) ptr1;
  nptr2 = (const struct s_s_assemblage *) ptr2;
  if (*nptr1 > nptr2->n_user)
    return (1);
  if (*nptr1 < nptr2->n_user)
    return (-1);
  return (0);
}

/* ---------------------------------------------------------------------- */
int
s_s_assemblage_copy (struct s_s_assemblage *s_s_assemblage_old_ptr,
		     struct s_s_assemblage *s_s_assemblage_new_ptr,
		     int n_user_new)
/* ---------------------------------------------------------------------- */
{
/*
 *   Copies s_s_assemblage data from s_s_assemblage_old_ptr to s_s_assemblage_new_ptr.
 *   Space for a s_s_assemblage_new_ptr must already be malloced.
 *   Space for list of pure phases in assemblage is malloced here.
 */
  int i, count_comps, count_s_s;
  char token[MAX_LENGTH];
/*
 *   Store data for structure s_s_assemblage
 */
  memcpy (s_s_assemblage_new_ptr, s_s_assemblage_old_ptr,
	  sizeof (struct s_s_assemblage));
  s_s_assemblage_new_ptr->n_user = n_user_new;
  s_s_assemblage_new_ptr->n_user_end = n_user_new;
  sprintf (token, "S_S_Assemblage defined in simulation %d.", simulation);
  s_s_assemblage_new_ptr->description = string_duplicate (token);
  s_s_assemblage_new_ptr->new_def = FALSE;
/*
 *   Malloc space for s_s structures and fill
 */
  count_s_s = s_s_assemblage_old_ptr->count_s_s;
  s_s_assemblage_new_ptr->s_s =
    (struct s_s *) PHRQ_malloc ((size_t) count_s_s * sizeof (struct s_s));
  if (s_s_assemblage_new_ptr->s_s == NULL)
    malloc_error ();
  memcpy ((void *) s_s_assemblage_new_ptr->s_s,
	  (void *) s_s_assemblage_old_ptr->s_s,
	  (size_t) count_s_s * sizeof (struct s_s));
/*
 *   Malloc space for components
 */
  for (i = 0; i < count_s_s; i++)
  {
    count_comps = s_s_assemblage_old_ptr->s_s[i].count_comps;
    s_s_assemblage_new_ptr->s_s[i].comps =
      (struct s_s_comp *) PHRQ_malloc ((size_t) count_comps *
				       sizeof (struct s_s_comp));
    if (s_s_assemblage_new_ptr->s_s[i].comps == NULL)
      malloc_error ();
    memcpy ((void *) s_s_assemblage_new_ptr->s_s[i].comps,
	    (void *) s_s_assemblage_old_ptr->s_s[i].comps,
	    (size_t) count_comps * sizeof (struct s_s_comp));
  }

  return (OK);
}

/* ---------------------------------------------------------------------- */
int
s_s_assemblage_copy_to_last (int n, int n_user)
/* ---------------------------------------------------------------------- */
{
/*
 *   Copies an s_s_assemblage definition from position n
 *   to position count_s_s_assemblage.
 */
  space ((void **) ((void *) &s_s_assemblage), count_s_s_assemblage,
	 &max_s_s_assemblage, sizeof (struct s_s_assemblage));
  s_s_assemblage_copy (&s_s_assemblage[n],
		       &s_s_assemblage[count_s_s_assemblage], n_user);
  count_s_s_assemblage++;
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
s_s_assemblage_duplicate (int n_user_old, int n_user_new)
/* ---------------------------------------------------------------------- */
{
/*
 *   Checks if n_user_new exists, and frees space if it does
 *   Copies s_s_assemblage[n_user_old] to old n_user_new space if
 *   found or to s_s_assemblage[count_s_s_assemblage] if not found.
 *   S_S_Assemblage array may not be in sort order after the copy.
 */
  int n_old, n_new, sort;
  struct s_s_assemblage *s_s_assemblage_ptr_old, *s_s_assemblage_ptr_new;
/*
 *   Find n_user_old in structure array s_s_assemblage
 */
  if (n_user_old == n_user_new)
    return (OK);
  s_s_assemblage_ptr_old = s_s_assemblage_bsearch (n_user_old, &n_old);
  if (s_s_assemblage_ptr_old == NULL)
  {
    sprintf (error_string, "S_S_Assemblage %d not found.", n_user_old);
    error_msg (error_string, CONTINUE);
    input_error++;
    return (ERROR);
  }
/*
 *   Find n_user_new in structure array s_s_assemblage or make new space
 */
  sort = FALSE;
  s_s_assemblage_ptr_new = s_s_assemblage_bsearch (n_user_new, &n_new);
  if (s_s_assemblage_ptr_new != NULL)
  {
    s_s_assemblage_free (s_s_assemblage_ptr_new);
  }
  else
  {
    space ((void **) ((void *) &s_s_assemblage), count_s_s_assemblage,
	   &max_s_s_assemblage, sizeof (struct s_s_assemblage));
    if (n_user_new < s_s_assemblage[count_s_s_assemblage - 1].n_user)
      sort = TRUE;
    n_new = count_s_s_assemblage++;
  }
/*
 *   Copy data
 */
  s_s_assemblage_copy (&s_s_assemblage[n_old], &s_s_assemblage[n_new],
		       n_user_new);
  if (sort == TRUE)
    s_s_assemblage_sort ();
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
s_s_assemblage_delete (int n_user_old)
/* ---------------------------------------------------------------------- */
/*
 *   Frees space for user number n_user_old, removes structure from
 *   array s_s_assemblage.
 */
{
  int i;
  int n_old;
  struct s_s_assemblage *s_s_assemblage_ptr_old;
/*
 *   Find n_user_old in structure array
 */
  s_s_assemblage_ptr_old = s_s_assemblage_bsearch (n_user_old, &n_old);
  if (s_s_assemblage_ptr_old != NULL)
  {
    /*
     *   Delete s_s_assemblage
     */
    s_s_assemblage_free (&s_s_assemblage[n_old]);

    for (i = n_old + 1; i < count_s_s_assemblage; i++)
    {
      memcpy ((void *) &s_s_assemblage[i - 1], (void *) &s_s_assemblage[i],
	      (size_t) sizeof (struct s_s_assemblage));
    }
    count_s_s_assemblage--;
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
s_s_assemblage_free (struct s_s_assemblage *s_s_assemblage_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Free space allocated for s_s_assemblage structure
 */
  int i;

  if (s_s_assemblage_ptr == NULL)
    return (ERROR);
  s_s_assemblage_ptr->description =
    (char *) free_check_null (s_s_assemblage_ptr->description);
  for (i = 0; i < s_s_assemblage_ptr->count_s_s; i++)
  {
    s_s_assemblage_ptr->s_s[i].comps =
      (struct s_s_comp *) free_check_null (s_s_assemblage_ptr->s_s[i].comps);
  }
  s_s_assemblage_ptr->s_s =
    (struct s_s *) free_check_null (s_s_assemblage_ptr->s_s);
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
s_s_assemblage_init (struct s_s_assemblage *s_s_assemblage_ptr, int n_user,
		     int n_user_end, char *description)
/* ---------------------------------------------------------------------- */
{
  s_s_assemblage_ptr->n_user = n_user;
  s_s_assemblage_ptr->n_user_end = n_user_end;
  s_s_assemblage_ptr->description = string_duplicate (description);
  s_s_assemblage_ptr->new_def = TRUE;
  s_s_assemblage_ptr->count_s_s = 0;
  s_s_assemblage_ptr->s_s =
    (struct s_s *) PHRQ_malloc ((size_t) sizeof (struct s_s));
  if (s_s_assemblage_ptr->s_s == NULL)
    malloc_error ();

  return (OK);
}

/* ---------------------------------------------------------------------- */
int
s_s_assemblage_ptr_to_user (struct s_s_assemblage *s_s_assemblage_ptr_old,
			    int n_user_new)
/* ---------------------------------------------------------------------- */
{
/*
 *   Checks if n_user_new exists, and frees space if it does
 *   Copies s_s_assemblage_ptr_old to old n_user_new space if
 *   found or to s_s_assemblage[count_s_s_assemblage] if not found.
 *   S_S_Assemblage array may not be in sort order after the copy.
 */
  int n_new, sort;
  struct s_s_assemblage *s_s_assemblage_ptr_new;
/*
 *   Find n_user_old in structure array s_s_assemblage
 */
  if (s_s_assemblage_ptr_old == NULL)
  {
    sprintf (error_string, "S_S_Assemblage pointer is NULL.");
    error_msg (error_string, CONTINUE);
    input_error++;
    return (ERROR);
  }
/*
 *   Find n_user_new in structure array s_s_assemblage or make new space
 */
  sort = FALSE;
  s_s_assemblage_ptr_new = s_s_assemblage_bsearch (n_user_new, &n_new);
  if (s_s_assemblage_ptr_new == s_s_assemblage_ptr_old)
    return (OK);
  if (s_s_assemblage_ptr_new != NULL)
  {
    s_s_assemblage_free (s_s_assemblage_ptr_new);
  }
  else
  {
    space ((void **) ((void *) &s_s_assemblage), count_s_s_assemblage,
	   &max_s_s_assemblage, sizeof (struct s_s_assemblage));
    if (n_user_new < s_s_assemblage[count_s_s_assemblage - 1].n_user)
      sort = TRUE;
    n_new = count_s_s_assemblage++;
  }
/*
 *   Copy data
 */
  s_s_assemblage_copy (s_s_assemblage_ptr_old, &s_s_assemblage[n_new],
		       n_user_new);
  if (sort == TRUE)
    s_s_assemblage_sort ();
  return (OK);
}

/* ---------------------------------------------------------------------- */
struct s_s_assemblage *
s_s_assemblage_replicate (struct s_s_assemblage *s_s_assemblage_old_ptr,
			  int n_user_new)
/* ---------------------------------------------------------------------- */
{
  struct s_s_assemblage *s_s_assemblage_ptr;
  s_s_assemblage_ptr = s_s_assemblage_alloc ();
  s_s_assemblage_copy (s_s_assemblage_old_ptr, s_s_assemblage_ptr,
		       n_user_new);
  return (s_s_assemblage_ptr);
}

/* ---------------------------------------------------------------------- */
struct s_s_assemblage *
s_s_assemblage_search (int n_user, int *n)
/* ---------------------------------------------------------------------- */
{
/*   Linear search of the structure array "s_s_assemblage" for user number n_user.
 *
 *   Arguments:
 *      n_user  input, user number
 *      n       output, position in s_s_assemblage
 *
 *   Returns:
 *      if found, the address of the s_s_assemblage element
 *      if not found, NULL
 */
  int i;
  for (i = 0; i < count_s_s_assemblage; i++)
  {
    if (s_s_assemblage[i].n_user == n_user)
    {
      *n = i;
      return (&(s_s_assemblage[i]));
    }
  }
/*
 *   S_S_Assemblage not found
 */
  return (NULL);
}

/* ---------------------------------------------------------------------- */
int
s_s_assemblage_sort (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Sort array of s_s_assemblage structures
 */
  if (count_s_s_assemblage > 0)
  {
    qsort (s_s_assemblage, (size_t) count_s_s_assemblage,
	   (size_t) sizeof (struct s_s_assemblage), s_s_assemblage_compare);
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
s_s_compare (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
  const struct s_s *s_s_ptr1, *s_s_ptr2;
  s_s_ptr1 = (const struct s_s *) ptr1;
  s_s_ptr2 = (const struct s_s *) ptr2;
  return (strcmp_nocase (s_s_ptr1->name, s_s_ptr2->name));

}

/* **********************************************************************
 *
 *   Routines related to structure "save_values"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
struct save_values *
save_values_bsearch (struct save_values *k, int *n)
/* ---------------------------------------------------------------------- */
{
/*
 *   Binary search save_values to find if one exists with given coefficients
 *   Save_Values is assumed to be in sort order by count_subscripts and
 *   values of subscripts
 */
  void *void_ptr;
  if (count_save_values == 0)
  {
    *n = -999;
    return (NULL);
  }
  void_ptr = (void *)
    bsearch ((char *) k,
	     (char *) save_values,
	     (size_t) count_save_values,
	     (size_t) sizeof (struct save_values), save_values_compare);
  if (void_ptr == NULL)
  {
    *n = -999;
    return (NULL);
  }
  *n = (int) ((struct save_values *) void_ptr - save_values);
  return ((struct save_values *) void_ptr);
}

/* ---------------------------------------------------------------------- */
int
save_values_compare (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
  int i;
  const struct save_values *save_values_ptr1, *save_values_ptr2;
  save_values_ptr1 = (const struct save_values *) ptr1;
  save_values_ptr2 = (const struct save_values *) ptr2;
  if (save_values_ptr1->count_subscripts < save_values_ptr2->count_subscripts)
  {
    return (-1);
  }
  else if (save_values_ptr1->count_subscripts >
	   save_values_ptr2->count_subscripts)
  {
    return (1);
  }
  else
  {
    for (i = 0; i < save_values_ptr1->count_subscripts; i++)
    {
      if (save_values_ptr1->subscripts[i] < save_values_ptr2->subscripts[i])
      {
	return (-1);
      }
      else if (save_values_ptr1->subscripts[i] >
	       save_values_ptr2->subscripts[i])
      {
	return (1);
      }
    }
  }
  return (0);
}

/* ---------------------------------------------------------------------- */
int
save_values_sort (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Sort array of save_values structures
 */
  if (count_save_values > 0)
  {
    qsort (save_values, (size_t) count_save_values,
	   (size_t) sizeof (struct save_values), save_values_compare);
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
save_values_store (struct save_values *s_v)
/* ---------------------------------------------------------------------- */
{
/*
 *   Look for subscripts
 */
  int n, i;
  struct save_values *s_v_ptr;

  s_v_ptr = save_values_bsearch (s_v, &n);
  if (s_v_ptr != NULL)
  {
    s_v_ptr->value = s_v->value;
  }
  else
  {
    save_values =
      (struct save_values *) PHRQ_realloc (save_values,
					   (size_t) (count_save_values +
						     1) *
					   sizeof (struct save_values));
    if (save_values == NULL)
      malloc_error ();
    save_values[count_save_values].value = s_v->value;
    save_values[count_save_values].count_subscripts = s_v->count_subscripts;
    i = s_v->count_subscripts;
    if (i == 0)
      i = 1;
    save_values[count_save_values].subscripts =
      (int *) PHRQ_malloc ((size_t) i * sizeof (int));
    if (save_values[count_save_values].subscripts == NULL)
      malloc_error ();
    save_values[count_save_values].subscripts =
      (int *) memcpy (save_values[count_save_values].subscripts,
		      s_v->subscripts, (size_t) i * sizeof (int));
    count_save_values++;
    save_values_sort ();
  }

  if (count_save_values > 0)
  {
    qsort (save_values, (size_t) count_save_values,
	   (size_t) sizeof (struct save_values), save_values_compare);
  }
  return (OK);
}

/* **********************************************************************
 *
 *   Routines related to structure "solution"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
int
conc_compare (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
  const struct conc *conc_ptr1, *conc_ptr2;
  conc_ptr1 = (const struct conc *) ptr1;
  conc_ptr2 = (const struct conc *) ptr2;
  return (strcmp_nocase (conc_ptr1->description, conc_ptr2->description));

}

/* ---------------------------------------------------------------------- */
int
conc_init (struct conc *conc_ptr)
/* ---------------------------------------------------------------------- */
{
  if (conc_ptr == NULL)
    return (ERROR);
  conc_ptr->equation_name = NULL;
  conc_ptr->phase = NULL;
  conc_ptr->phase_si = 0.0;
  conc_ptr->units = NULL;
  conc_ptr->n_pe = -1;
  conc_ptr->as = NULL;
  conc_ptr->gfw = 0.0;
  return OK;
}

/* ---------------------------------------------------------------------- */
int
isotope_compare (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
  int i;
  const struct isotope *iso_ptr1, *iso_ptr2;

  iso_ptr1 = (const struct isotope *) ptr1;
  iso_ptr2 = (const struct isotope *) ptr2;
  i = strcmp_nocase (iso_ptr1->elt_name, iso_ptr2->elt_name);
  if (i != 0)
    return (i);
  if (iso_ptr1->isotope_number < iso_ptr2->isotope_number)
  {
    return (-1);
  }
  else if (iso_ptr1->isotope_number > iso_ptr2->isotope_number)
  {
    return (1);
  }
  return (0);
}

/* ---------------------------------------------------------------------- */
struct solution *
solution_alloc (void)
/* ---------------------------------------------------------------------- */
/*
 *   Allocates space to a solution structure
 *      arguments
 *      input:  none
 *      return: pointer to a solution structure
 */
{
  int max_mass_balance;
  struct solution *solution_ptr;
  solution_ptr = (struct solution *) PHRQ_malloc (sizeof (struct solution));
  if (solution_ptr == NULL)
    malloc_error ();
/*
 *   set pointers in structure to NULL
 */
  solution_ptr->new_def = TRUE;
  solution_ptr->n_user = -1;
  solution_ptr->n_user_end = -1;
  solution_ptr->description = NULL;
  solution_ptr->tc = 0;
  solution_ptr->ph = 0;
  solution_ptr->solution_pe = 0;
  solution_ptr->mu = 0;
  solution_ptr->ah2o = 0;
  solution_ptr->density = 0;
  solution_ptr->total_h = 0;
  solution_ptr->total_o = 0;
  solution_ptr->cb = 0;
  solution_ptr->mass_water = 0;
  solution_ptr->total_alkalinity = 0;
  solution_ptr->units = NULL;
/*
 *   Initial allocation of space for pe's
 */
  solution_ptr->pe = pe_data_alloc ();
  solution_ptr->default_pe = 0;
/*
 *   Initial allocation of space for totals
 */
  max_mass_balance = MAX_MASS_BALANCE;
  space ((void **) ((void *) &(solution_ptr->totals)), INIT, &max_mass_balance,
	 sizeof (struct conc));
  solution_ptr->totals[0].description = NULL;
/*
 *   Master activities
 */
  space ((void **) ((void *) &(solution_ptr->master_activity)), INIT, &max_mass_balance,
	 sizeof (struct master_activity));
  solution_ptr->master_activity[0].description = NULL;
/*
 *   Initial allocation of space for isotopes
 */
  solution_ptr->isotopes = (struct isotope *) PHRQ_malloc (sizeof (struct isotope));
  if (solution_ptr->isotopes == NULL) malloc_error ();
  solution_ptr->count_isotopes = 0;
/*
 *   Initial allocation of space for species_gamma
 */
  solution_ptr->species_gamma = NULL;
  solution_ptr->count_species_gamma = 0;

  return (solution_ptr);
}

/* ---------------------------------------------------------------------- */
struct solution *
solution_bsearch (int k, int *n, int print)
/* ---------------------------------------------------------------------- */
{
/*
 *   Binary search of solution solution array to find user number k.
 *   Assumes array solution is in sort order.
 *
 *   Input: k user number to find.
 *          print, TRUE print warning if solution not found.
 *   Output: n, position of solution ptr in solution array.
 *   Return: pointer to solution structure if found
 *           NULL if not found.
 */
  void *void_ptr;
  void_ptr = NULL;
  if (count_solution > 0)
  {
    void_ptr = (void *)
      bsearch ((char *) &k,
	       (char *) solution,
	       (size_t) count_solution,
	       (size_t) sizeof (struct solution *), solution_compare_int);
  }
  if (void_ptr == NULL && print == TRUE)
  {
    sprintf (error_string, "Solution %d not found.", k);
    error_msg (error_string, CONTINUE);
  }
  if (void_ptr == NULL)
  {
    *n = -999;
    return (NULL);
  }

  *n = (int) ((struct solution **) void_ptr - solution);
  return (*(struct solution **) void_ptr);
}

/* ---------------------------------------------------------------------- */
int
solution_compare (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
/*
 *   Compare solution user numbers
 */
  const struct solution *nptr1;
  const struct solution *nptr2;

  nptr1 = *(const struct solution **) ptr1;
  nptr2 = *(const struct solution **) ptr2;
  if (nptr1->n_user > nptr2->n_user)
    return (1);
  if (nptr1->n_user < nptr2->n_user)
    return (-1);
  return (0);
}

/* ---------------------------------------------------------------------- */
int
solution_compare_int (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
/*
 *   Compare solution user number to an integer
 */
  const int *nptr1;
  const struct solution *nptr2;

  nptr1 = (const int *) ptr1;
  nptr2 = *(const struct solution **) ptr2;
  if (*nptr1 > nptr2->n_user)
    return (1);
  if (*nptr1 < nptr2->n_user)
    return (-1);
  return (0);
}

/* ---------------------------------------------------------------------- */
struct solution *
solution_copy (struct solution *solution_old_ptr, int n_user_new)
/* ---------------------------------------------------------------------- */
{
/*
 *   Copies solution data to new structure.
 *   Space for new structure is malloced here.
 *
 *   Input: solution_old_ptr, pointer to structure to be copied.
 *          n_user_new, user number to give to the new solution.
 *
 *   Return: pointer to new solution structure.
 */
  int i;
  int count_totals, count_master_activity;
  struct solution *solution_new_ptr;
/*
 *   malloc space for new solution structure
 */
  solution_new_ptr =
    (struct solution *) PHRQ_malloc (sizeof (struct solution));
  if (solution_new_ptr == NULL)
    malloc_error ();
/*
 *   Copy solution data, but set pointers pe and totals to NULL
 */
  memcpy ((void *) solution_new_ptr, (void *) solution_old_ptr,
	  (size_t) sizeof (struct solution));
  solution_new_ptr->n_user = n_user_new;
  solution_new_ptr->n_user_end = n_user_new;
  solution_new_ptr->description =
    string_duplicate (solution_old_ptr->description);
/*
 *   Count totals data and malloc space
 */
  for (i = 0; solution_old_ptr->totals[i].description != NULL; i++);
  count_totals = i + 1;
  solution_new_ptr->totals =
    (struct conc *) PHRQ_malloc ((size_t) count_totals *
				 sizeof (struct conc));
  if (solution_new_ptr->totals == NULL)
    malloc_error ();
/*
 *   Copy totals data
 */
  memcpy ((void *) solution_new_ptr->totals,
	  (void *) solution_old_ptr->totals,
	  (size_t) count_totals * sizeof (struct conc));
/*
 *   Count master activity guesses and malloc space
 */
  count_master_activity = solution_old_ptr->count_master_activity;
  solution_new_ptr->master_activity =
    (struct master_activity *) PHRQ_malloc ((size_t) count_master_activity *
					    sizeof (struct master_activity));
  if (solution_new_ptr->master_activity == NULL)
    malloc_error ();
/*
 *   Copy master activity guesses
 */
  memcpy ((void *) solution_new_ptr->master_activity,
	  (void *) solution_old_ptr->master_activity,
	  (size_t) count_master_activity * sizeof (struct master_activity));
/*
 *   malloc space for species gammas
 */
  if (solution_old_ptr->count_species_gamma > 0)
  {
    solution_new_ptr->species_gamma =
      (struct master_activity *) PHRQ_malloc ((size_t) solution_old_ptr->
					      count_species_gamma *
					      sizeof (struct
						      master_activity));
    if (solution_new_ptr->species_gamma == NULL)
      malloc_error ();
    /*
     *   Copy species gammas
     */
    memcpy ((void *) solution_new_ptr->species_gamma,
	    (void *) solution_old_ptr->species_gamma,
	    (size_t) solution_old_ptr->count_species_gamma *
	    sizeof (struct master_activity));
    solution_new_ptr->count_species_gamma =
      solution_old_ptr->count_species_gamma;
  }
  else
  {
    solution_new_ptr->species_gamma = NULL;
    solution_new_ptr->count_species_gamma = 0;
  }
/*
 *   Copy pe data
 */
  solution_new_ptr->pe = pe_data_dup (solution_old_ptr->pe);
/*
 *  Malloc and copy isotope data
 */
  if (solution_old_ptr->count_isotopes > 0)
  {
    solution_new_ptr->count_isotopes = solution_old_ptr->count_isotopes;
    solution_new_ptr->isotopes =
      (struct isotope *) PHRQ_malloc ((size_t) solution_old_ptr->
				      count_isotopes *
				      sizeof (struct isotope));
    if (solution_new_ptr->isotopes == NULL)
      malloc_error ();
    memcpy (solution_new_ptr->isotopes, solution_old_ptr->isotopes,
	    (size_t) solution_old_ptr->count_isotopes *
	    sizeof (struct isotope));
  }
  else
  {
    solution_new_ptr->count_isotopes = 0;
    solution_new_ptr->isotopes = NULL;
  }
/*
 *   Return pointer to new structure
 */
  return (solution_new_ptr);
}

/* ---------------------------------------------------------------------- */
int
solution_copy_to_last (int n, int n_user_new)
/* ---------------------------------------------------------------------- */
{
/*
 *   Copies solution from pointer solution[n]. Makes new structure
 *   and saves pointer in solution[count_solution].
 */
  /* Make sure solution array is large enough */
  if (count_solution + 1 >= max_solution)
  {
    space ((void **) ((void *) &(solution)), count_solution, &max_solution,
	   sizeof (struct solution *));
  }
  /* malloc space and copy solution */
  solution[count_solution] = solution_copy (solution[n], n_user_new);
  count_solution++;
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
solution_duplicate (int n_user_old, int n_user_new)
/* ---------------------------------------------------------------------- */
{
/*
 *   Checks if n_user_new exists, and frees structure if it does
 *   Copies structure for solution[n_user_old] to new space
 *   Saves pointer in old position if it existed, no sort necessary
 *   Otherwise saves in count_solution position, sort may be necessary
 *   Solution array is in sort order on exit
 */
  int n, n_old, n_new, sort;
  struct solution *solution_old_ptr, *solution_new_ptr;
/*
 *   Find solution n_user_old
 */
  if (n_user_old == n_user_new)
    return (OK);
  solution_old_ptr = solution_bsearch (n_user_old, &n_old, TRUE);
  if (solution_old_ptr == NULL)
  {
    sprintf (error_string, "Solution %d not found.", n_user_old);
    error_msg (error_string, STOP);
  }
/*
 *   Check if solution n_user_new already exists
 */
  solution_new_ptr = solution_bsearch (n_user_new, &n_new, FALSE);
  if (solution_new_ptr != NULL)
  {
    n = n_new;
    solution_free (solution_new_ptr);
    sort = FALSE;
  }
  else
  {
    n = count_solution;
    count_solution++;
    if (n_user_new > solution[n - 1]->n_user)
    {
      sort = FALSE;
    }
    else
    {
      sort = TRUE;
    }
  }
  solution[n] = solution_copy (solution_old_ptr, n_user_new);
/*
 *   Make sure surface array is large enough
 */
  if (count_solution >= max_solution)
  {
    space ((void **) ((void *) &(solution)), count_solution, &max_solution,
	   sizeof (struct solution *));
  }
/*
 *   Sort solution if necessary
 */
  if (sort == TRUE)
    solution_sort ();
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
solution_delete (int n_user)
/* ---------------------------------------------------------------------- */
{
/*
 *   Delete a solution structure: free memory and renumber solution.
 */
  struct solution *solution_ptr;
  int j, n;
  solution_ptr = solution_bsearch (n_user, &n, FALSE);
  if (solution_ptr != NULL)
  {
    solution_free (solution[n]);
    for (j = n; j < (count_solution - 1); j++)
    {
      solution[j] = solution[j + 1];
    }
    count_solution--;
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
solution_free (struct solution *solution_ptr)
/* ---------------------------------------------------------------------- */
{
  if (solution_ptr == NULL)
    return (OK);
/*
 *   Free all memory unique to a solution structure, free solution_ptr too.
 */
  solution_ptr->description =
    (char *) free_check_null (solution_ptr->description);
/*   Free struct conc array */
  solution_ptr->totals =
    (struct conc *) free_check_null (solution_ptr->totals);
  solution_ptr->master_activity =
    (struct master_activity *) free_check_null (solution_ptr->
						master_activity);
  solution_ptr->species_gamma =
    (struct master_activity *) free_check_null (solution_ptr->species_gamma);
  solution_ptr->count_species_gamma = 0;
/*   Free struct pe_data array */
  pe_data_free (solution_ptr->pe);
/*   Free isotope data */
  solution_ptr->isotopes =
    (struct isotope *) free_check_null (solution_ptr->isotopes);
/*   Free struct solution */
  solution_ptr = (struct solution *) free_check_null (solution_ptr);

  return (OK);
}

/* ---------------------------------------------------------------------- */
int
solution_ptr_to_user (struct solution *solution_old_ptr, int n_user_new)
/* ---------------------------------------------------------------------- */
{
/*
 *   Checks if n_user_new exists, and frees structure if it does
 *   Copies structure for solution[n_user_old] to new space
 *   Saves pointer in old position if it existed, no sort necessary
 *   Otherwise saves in count_solution position, sort may be necessary
 *   Solution array is in sort order on exit
 */
  int n, n_new, sort;
  struct solution *solution_new_ptr;
/*
 *   Find solution n_user_old
 */
  if (solution_old_ptr == NULL)
  {
    sprintf (error_string, "Solution pointer is NULL.");
    error_msg (error_string, STOP);
  }
/*
 *   Check if solution n_user_new already exists
 */
  solution_new_ptr = solution_bsearch (n_user_new, &n_new, FALSE);
  if (solution_new_ptr == solution_old_ptr)
    return (OK);
  if (solution_new_ptr != NULL)
  {
    n = n_new;
    solution_free (solution_new_ptr);
    sort = FALSE;
  }
  else
  {
    n = count_solution;
    count_solution++;
    /*
     *   Make sure surface array is large enough
     */
    if (count_solution >= max_solution)
    {
      space ((void **) ((void *) &(solution)), count_solution, &max_solution,
	     sizeof (struct solution *));
    }
    if (n_user_new > solution[n - 1]->n_user)
    {
      sort = FALSE;
    }
    else
    {
      sort = TRUE;
    }
  }
  solution[n] = solution_copy (solution_old_ptr, n_user_new);
/*
 *   Sort solution if necessary
 */
  if (sort == TRUE)
    solution_sort ();
  return (OK);
}

/* ---------------------------------------------------------------------- */
struct solution *
solution_replicate (struct solution *solution_old_ptr, int n_user_new)
/* ---------------------------------------------------------------------- */
{
  return (solution_copy (solution_old_ptr, n_user_new));
}

/* ---------------------------------------------------------------------- */
int
solution_sort (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Sort array of pointers to solution structures, solution.
 */
  if (count_solution > 0)
  {
    qsort (solution,
	   (size_t) count_solution,
	   (size_t) sizeof (struct solution *), solution_compare);
  }
  return (OK);
}

/* **********************************************************************
 *
 *   Routines related to structure "species_list"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
static int
species_list_compare (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
  int j;
  char *name1, *name2;
  const struct species_list *nptr1, *nptr2;

  nptr1 = (const struct species_list *) ptr1;
  nptr2 = (const struct species_list *) ptr2;

/*
 *   Put H+ first
 */
  if (nptr1->master_s != nptr2->master_s)
  {
    if (nptr1->master_s == s_hplus)
      return (-1);
    if (nptr2->master_s == s_hplus)
      return (1);
  }
/*
 *   Other element valence states
 */
  if (nptr1->master_s->secondary != NULL)
  {
    name1 = nptr1->master_s->secondary->elt->name;
  }
  else
  {
    name1 = nptr1->master_s->primary->elt->name;
  }
  if (nptr2->master_s->secondary != NULL)
  {
    name2 = nptr2->master_s->secondary->elt->name;
  }
  else
  {
    name2 = nptr2->master_s->primary->elt->name;
  }
/*
 *   Compare name of primary or secondary master species; log molality
 */

  j = strcmp (name1, name2);

/*
 *   Different master species
 */
  if (j != 0)
    return (j);

/*
 *   Else, descending order by log molality
 */
  if (nptr1->s->lm > nptr2->s->lm)
  {
    return (-1);
  }
  else if (nptr1->s->lm < nptr2->s->lm)
  {
    return (1);
  }
  else
  {
    return (0);
  }
}

/* ---------------------------------------------------------------------- */
int
species_list_compare_alk (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
  const struct species_list *nptr1, *nptr2;
  double alk1, alk2;

  nptr1 = (const struct species_list *) ptr1;
  nptr2 = (const struct species_list *) ptr2;
/*
 *   Else, descending order by log molality
 */
  alk1 = fabs (under (nptr1->s->lm) * nptr1->s->alk);
  alk2 = fabs (under (nptr2->s->lm) * nptr2->s->alk);

  if (alk1 > alk2)
  {
    return (-1);
  }
  else if (alk1 < alk2)
  {
    return (1);
  }
  else
  {
    return (0);
  }
}

/* ---------------------------------------------------------------------- */
int
species_list_compare_master (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
  char *name1, *name2;
  const struct species_list *nptr1, *nptr2;

  nptr1 = (const struct species_list *) ptr1;
  nptr2 = (const struct species_list *) ptr2;

/*
 *   Put H+ first
 */
  if (nptr1->master_s != nptr2->master_s)
  {
    if (nptr1->master_s == s_hplus)
      return (-1);
    if (nptr2->master_s == s_hplus)
      return (1);
  }
/*
 *   Other element valence states
 */
  if (nptr1->master_s->secondary != NULL)
  {
    name1 = nptr1->master_s->secondary->elt->name;
  }
  else
  {
    name1 = nptr1->master_s->primary->elt->name;
  }
  if (nptr2->master_s->secondary != NULL)
  {
    name2 = nptr2->master_s->secondary->elt->name;
  }
  else
  {
    name2 = nptr2->master_s->primary->elt->name;
  }
/*
 *   Compare name of primary or secondary master species; log molality
 */

  return (strcmp (name1, name2));
}

/* ---------------------------------------------------------------------- */
int
species_list_sort (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Sort list using rules in species_list_compare
 */
  if (count_species_list > 0)
  {
    qsort (&species_list[0], (size_t) count_species_list,
	   (size_t) sizeof (struct species_list), species_list_compare);
  }
  return (OK);
}

/* **********************************************************************
 *
 *   Routines related to structure "surface"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
struct Change_Surf *
change_surf_alloc (int count)
/* ---------------------------------------------------------------------- */
{
  if (count == 1)
    return (change_surf);
  change_surf =
    (struct Change_Surf *) PHRQ_realloc (change_surf,
					 (size_t) count *
					 sizeof (struct Change_Surf));
  if (change_surf == NULL)
    malloc_error ();
  change_surf[count - 1].cell_no = -99;
  change_surf[count - 1].next = FALSE;
  change_surf[count - 2].next = TRUE;

  return (change_surf);
}

/* ---------------------------------------------------------------------- */
struct surface *
surface_alloc (void)
/* ---------------------------------------------------------------------- */
{
  struct surface *surface_ptr;
  surface_ptr = (struct surface *) PHRQ_malloc (sizeof (struct surface));
  if (surface_ptr == NULL)
    malloc_error ();
  surface_ptr->n_user = -1;
  surface_ptr->n_user_end = -1;
  surface_ptr->new_def = 0;
  surface_ptr->only_counter_ions = 0;
  surface_ptr->dl_type = BORKOVEK_DL;
  surface_ptr->type = UNKNOWN_DL;
  surface_ptr->sites_units = SITES_ABSOLUTE;
  surface_ptr->thickness = 0;
  surface_ptr->debye_lengths = 0;
  surface_ptr->DDL_viscosity = 0;		/* viscosity relative to pure water */
  surface_ptr->DDL_limit = 0;		/* limits DDL water to this fraction of bulk water */
  surface_ptr->description = NULL;
  surface_ptr->solution_equilibria = 0;
  surface_ptr->n_solution = 0;
  surface_ptr->count_comps = 0;
  surface_ptr->comps = NULL;
  surface_ptr->count_charge = 0;
  surface_ptr->charge = NULL;
  surface_ptr->related_phases = 0;
  surface_ptr->related_rate = 0;
  surface_ptr->transport = 0;		/* transports comp's and charges if true */
  return (surface_ptr);
}

/* ---------------------------------------------------------------------- */
struct surface *
surface_bsearch (int k, int *n)
/* ---------------------------------------------------------------------- */
{
/*
 *   Binary search of array surface to find user number k.
 *   Assumes array surface is in sort order.
 *
 *   Input: k, user number to find.
 *   Output: n, position in array surface of user number k.
 *   Return: pointer to surface structure if found,
 *           NULL if not found.
 */
  void *void_ptr;
  if (count_surface == 0)
  {
    *n = -999;
    return (NULL);
  }
  void_ptr = (void *)
    bsearch ((char *) &k,
	     (char *) surface,
	     (size_t) count_surface,
	     (size_t) sizeof (struct surface), surface_compare_int);
  if (void_ptr == NULL)
  {
    *n = -999;
    return (NULL);
  }
  *n = (int) ((struct surface *) void_ptr - surface);
  return ((struct surface *) void_ptr);
}

/* ---------------------------------------------------------------------- */
struct master *
surface_get_psi_master (char *name, int plane)
/* ---------------------------------------------------------------------- */
{
  struct master *master_ptr;
  char token[MAX_LENGTH];

  if (name == NULL)
    return (NULL);
  strcpy (token, name);
  strcat (token, "_psi");
  switch (plane)
  {
  case SURF_PSI:
    break;
  case SURF_PSI1:
    strcat (token, "b");
    break;
  case SURF_PSI2:
    strcat (token, "d");
    break;
  default:
    error_msg ("Unknown plane for surface_get_psi_master", STOP);
  }
  master_ptr = master_bsearch (token);
  return (master_ptr);
}

/* ---------------------------------------------------------------------- */
int
surface_comp_compare (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
  const struct surface_comp *surface_comp_ptr1, *surface_comp_ptr2;
  surface_comp_ptr1 = (const struct surface_comp *) ptr1;
  surface_comp_ptr2 = (const struct surface_comp *) ptr2;
  return (strcmp_nocase
	  (surface_comp_ptr1->master->elt->name,
	   surface_comp_ptr2->master->elt->name));
}

/* ---------------------------------------------------------------------- */
int
surface_compare (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
  const struct surface *surface_ptr1, *surface_ptr2;
  surface_ptr1 = (const struct surface *) ptr1;
  surface_ptr2 = (const struct surface *) ptr2;
  if (surface_ptr1->n_user > surface_ptr2->n_user)
    return (1);
  if (surface_ptr1->n_user < surface_ptr2->n_user)
    return (-1);
  return (0);

}

/* ---------------------------------------------------------------------- */
static int
surface_compare_int (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
/*
 *   Compare surface user numbers to integer
 */
  const int *nptr1;
  const struct surface *nptr2;

  nptr1 = (const int *) ptr1;
  nptr2 = (const struct surface *) ptr2;
  if (*nptr1 > nptr2->n_user)
    return (1);
  if (*nptr1 < nptr2->n_user)
    return (-1);
  return (0);
}

/* ---------------------------------------------------------------------- */
int
surface_copy (struct surface *surface_old_ptr,
	      struct surface *surface_new_ptr, int n_user_new)
/* ---------------------------------------------------------------------- */
{
/*
 *   Copies data from surface_old_ptr to surface_new_ptr.
 *   Space for a surface structure must already be malloced.
 *   User number of new surface structure is n_user.
 */
  int count_comps, count_charge, i;
  char token[MAX_LENGTH];
/*
 *   Store data for structure surface
 */
  memcpy (surface_new_ptr, surface_old_ptr, sizeof (struct surface));
  count_comps = surface_old_ptr->count_comps;
  count_charge = surface_old_ptr->count_charge;
  surface_new_ptr->n_user = n_user_new;
  surface_new_ptr->n_user_end = n_user_new;
  surface_new_ptr->new_def = surface_old_ptr->new_def;
  sprintf (token, "Surface defined in simulation %d.", simulation);
  surface_new_ptr->description = string_duplicate (token);
  surface_new_ptr->solution_equilibria = FALSE;
  surface_new_ptr->n_solution = -10;
/*
 *   Write surface_comp structure for each surface component
 */
  surface_new_ptr->comps =
    (struct surface_comp *) PHRQ_malloc ((size_t) (count_comps) *
					 sizeof (struct surface_comp));
  if (surface_new_ptr->comps == NULL)
    malloc_error ();
  memcpy (surface_new_ptr->comps, surface_old_ptr->comps,
	  (size_t) (count_comps) * sizeof (struct surface_comp));
  for (i = 0; i < count_comps; i++)
  {
    surface_new_ptr->comps[i].formula_totals =
      elt_list_dup (surface_old_ptr->comps[i].formula_totals);
    surface_new_ptr->comps[i].totals =
      elt_list_dup (surface_old_ptr->comps[i].totals);
  }
/*
 *   Write surface_charge structure for each surface
 */
  /*if (surface_old_ptr->edl == TRUE) { */
  if (surface_old_ptr->type == DDL || surface_old_ptr->type == CD_MUSIC)
  {
    surface_new_ptr->charge =
      (struct surface_charge *) PHRQ_malloc ((size_t) (count_charge) *
					     sizeof (struct surface_charge));
    if (surface_new_ptr->charge == NULL)
      malloc_error ();
    memcpy (surface_new_ptr->charge, surface_old_ptr->charge,
	    (size_t) (count_charge) * sizeof (struct surface_charge));
    for (i = 0; i < count_charge; i++)
    {
      surface_new_ptr->charge[i].diffuse_layer_totals =
	elt_list_dup (surface_old_ptr->charge[i].diffuse_layer_totals);
      surface_new_ptr->charge[i].g = NULL;
    }
  }
  else
  {
    surface_new_ptr->count_charge = 0;
    surface_new_ptr->charge = NULL;
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
struct surface_charge *
surface_charge_duplicate (struct surface_charge *charge_old_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Duplicates data from surface_old_ptr
 *   Space is malloced.
 */
  struct surface_charge *charge;
/*
 *   Write surface_charge structure for each surface
 */
  charge =
    (struct surface_charge *) PHRQ_malloc (sizeof (struct surface_charge));
  if (charge == NULL)
    malloc_error ();
  memcpy (charge, charge_old_ptr, sizeof (struct surface_charge));
  charge->diffuse_layer_totals =
    elt_list_dup (charge_old_ptr->diffuse_layer_totals);
  charge->count_g = 0;
  charge->g = NULL;
  return (charge);
}

/* ---------------------------------------------------------------------- */
int
surface_charge_free (struct surface_charge *charge)
/* ---------------------------------------------------------------------- */
{
/*
 *    Frees all space related to surface_charge
 */
  if (charge == NULL)
    return (ERROR);
  charge->diffuse_layer_totals =
    (struct elt_list *) free_check_null (charge->diffuse_layer_totals);
  charge->g = (struct surface_diff_layer *) free_check_null (charge->g);
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
surface_copy_to_last (int n, int n_user)
/* ---------------------------------------------------------------------- */
{
/*
 *   Copies an surface definition to position count_surface.
 */
  space ((void **) ((void *) &surface), count_surface, &max_surface,
	 sizeof (struct surface));
  surface_copy (&surface[n], &surface[count_surface], n_user);
  count_surface++;
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
surface_duplicate (int n_user_old, int n_user_new)
/* ---------------------------------------------------------------------- */
{
/*
 *   Checks if n_user_new exists, and frees space if it does
 *   Copies surface[n_user_old] to old n_user_new space if
 *   found or to surface[count_surface] if not found.
 *   Surface array may not be in sort order after the copy.
 */
  int n_old, n_new, sort;
  struct surface *surface_ptr_old, *surface_ptr_new;
/*
 *   Find n_user_old in structure array surface
 */
  if (n_user_old == n_user_new)
    return (OK);
  surface_ptr_old = surface_bsearch (n_user_old, &n_old);
  if (surface_ptr_old == NULL)
  {
    sprintf (error_string, "Surface %d not found.", n_user_old);
    error_msg (error_string, CONTINUE);
    input_error++;
    return (ERROR);
  }
/*
 *   Find n_user_new in structure array surface or make new space
 */
  sort = FALSE;
  surface_ptr_new = surface_bsearch (n_user_new, &n_new);
  if (surface_ptr_new != NULL)
  {
    surface_free (surface_ptr_new);
  }
  else
  {
    space ((void **) ((void *) &surface), count_surface, &max_surface,
	   sizeof (struct surface));
    if (n_user_new < surface[count_surface - 1].n_user)
      sort = TRUE;
    n_new = count_surface++;
  }
/*
 *   Copy data
 */
  surface_copy (&surface[n_old], &surface[n_new], n_user_new);
  if (sort == TRUE)
    surface_sort ();
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
surface_delete (int n_user_old)
/* ---------------------------------------------------------------------- */
/*
 *   Frees space for user number n_user_old, removes structure from
 *   array surface.
 */
{
  int i;
  int n_old;
  struct surface *surface_ptr_old;
/*
 *   Find n_user_old in structure array
 */
  surface_ptr_old = surface_bsearch (n_user_old, &n_old);
  if (surface_ptr_old != NULL)
  {
    /*
     *   Delete surface
     */
    surface_free (&surface[n_old]);

    for (i = n_old + 1; i < count_surface; i++)
    {
      memcpy ((void *) &surface[i - 1], (void *) &surface[i],
	      (size_t) sizeof (struct surface));
    }
    count_surface--;
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
surface_free (struct surface *surface_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *    Frees all space related to surface_ptr, but not surface_ptr.
 */
  int k;

  if (surface_ptr == NULL)
    return (ERROR);
  surface_ptr->description =
    (char *) free_check_null (surface_ptr->description);
/*
 *   totals, then comps
 */
  for (k = 0; k < surface_ptr->count_comps; k++)
  {
    surface_ptr->comps[k].formula_totals =
      (struct elt_list *) free_check_null (surface_ptr->comps[k].
					   formula_totals);
    surface_ptr->comps[k].totals =
      (struct elt_list *) free_check_null (surface_ptr->comps[k].totals);
  }
  surface_ptr->comps =
    (struct surface_comp *) free_check_null (surface_ptr->comps);
/*
 *   diffuse_layer_totals and g, then charge
 */
  /*if (surface_ptr->edl == TRUE) { */
  if (surface_ptr->type == DDL || surface_ptr->type == CD_MUSIC)
  {
    for (k = 0; k < surface_ptr->count_charge; k++)
    {
      surface_ptr->charge[k].diffuse_layer_totals =
	(struct elt_list *) free_check_null (surface_ptr->charge[k].
					     diffuse_layer_totals);
      surface_ptr->charge[k].g =
	(struct surface_diff_layer *) free_check_null (surface_ptr->charge[k].
						       g);
    }
  }
  surface_ptr->charge =
    (struct surface_charge *) free_check_null (surface_ptr->charge);
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
surface_init (struct surface *surface_ptr, int n_user, int n_user_end,
	      char *description)
/* ---------------------------------------------------------------------- */
{

  if (surface_ptr == NULL)
    return (ERROR);

  surface_ptr->n_user = n_user;
  surface_ptr->n_user_end = n_user_end;
  surface_ptr->new_def = TRUE;
  surface_ptr->dl_type = NO_DL;
  surface_ptr->sites_units = SITES_ABSOLUTE;
  /*surface_ptr->edl = TRUE; */
  surface_ptr->type = DDL;
  surface_ptr->only_counter_ions = FALSE;
  /*surface_ptr->donnan = FALSE; */
  surface_ptr->thickness = 1e-8;
  surface_ptr->debye_lengths = 0.0;
  surface_ptr->DDL_viscosity = 1.0;
  surface_ptr->DDL_limit = 0.8;
  surface_ptr->description = string_duplicate (description);
  surface_ptr->solution_equilibria = FALSE;
  surface_ptr->n_solution = -999;
  surface_ptr->transport = FALSE;
/*
 *   Malloc one surface_comp structure
 */
  surface_ptr->count_comps = 0;
  surface_ptr->comps =
    (struct surface_comp *) PHRQ_malloc ((size_t)
					 sizeof (struct surface_comp));
  if (surface_ptr->comps == NULL)
    malloc_error ();
  surface_ptr->comps[0].master = NULL;
  surface_ptr->comps[0].totals = NULL;
  surface_ptr->comps[0].phase_name = NULL;
  surface_ptr->comps[0].rate_name = NULL;
  surface_ptr->comps[0].charge = -1;
  surface_ptr->comps[0].moles = 0;
  surface_ptr->comps[0].cb = 0;
  surface_ptr->comps[0].Dw = 0;
/*
 *   Malloc one charge structure
 */
  surface_ptr->count_charge = 0;
  surface_ptr->charge =
    (struct surface_charge *) PHRQ_malloc ((size_t)
					   sizeof (struct surface_charge));
  if (surface_ptr->charge == NULL)
    malloc_error ();

  surface_ptr->related_phases = FALSE;
  surface_ptr->related_rate = FALSE;

  return (OK);
}

/* ---------------------------------------------------------------------- */
int
surface_ptr_to_user (struct surface *surface_ptr_old, int n_user_new)
/* ---------------------------------------------------------------------- */
{
/*
 *   Checks if n_user_new exists, and frees space if it does
 *   Copies surface_ptr_old to old n_user_new space if
 *   found or to surface[count_surface] if not found.
 *   Surface array may not be in sort order after the copy.
 */
  int n_new, sort;
  struct surface *surface_ptr_new;
/*
 *   Find n_user_old in structure array surface
 */
  if (surface_ptr_old == NULL)
  {
    sprintf (error_string, "Surface pointer is NULL.");
    error_msg (error_string, CONTINUE);
    input_error++;
    return (ERROR);
  }
/*
 *   Find n_user_new in structure array surface or make new space
 */
  sort = FALSE;
  surface_ptr_new = surface_bsearch (n_user_new, &n_new);
  if (surface_ptr_new == surface_ptr_old)
    return (OK);
  if (surface_ptr_new != NULL)
  {
    surface_free (surface_ptr_new);
  }
  else
  {
    space ((void **) ((void *) &surface), count_surface, &max_surface,
	   sizeof (struct surface));
    if (n_user_new < surface[count_surface - 1].n_user)
      sort = TRUE;
    n_new = count_surface++;
  }
/*
 *   Copy data
 */
  surface_copy (surface_ptr_old, &surface[n_new], n_user_new);
  if (sort == TRUE)
    surface_sort ();
  return (OK);
}

/* ---------------------------------------------------------------------- */
struct surface *
surface_replicate (struct surface *surface_old_ptr, int n_user_new)
/* ---------------------------------------------------------------------- */
{
  struct surface *surface_ptr;
  surface_ptr = surface_alloc ();
  surface_copy (surface_old_ptr, surface_ptr, n_user_new);
  return (surface_ptr);
}

/* ---------------------------------------------------------------------- */
struct surface *
surface_search (int n_user, int *n, int print)
/* ---------------------------------------------------------------------- */
{
/*   Linear search of the structure array "surface" for user number n_user.
 *
 *   Arguments:
 *      n_user  input, user number.
 *      n       output, position in surface.
 *
 *   Returns:
 *      if found, pointer to surface structure
 *      if not found, NULL
 */
  int i;
  for (i = 0; i < count_surface; i++)
  {
    if (n_user == surface[i].n_user)
    {
      break;
    }
  }
  if (i >= count_surface)
  {
    if (print == TRUE)
    {
      sprintf (error_string, "Surface %d not found.", n_user);
      error_msg (error_string, CONTINUE);
    }
    *n = -999;
    return (NULL);
  }
  *n = i;
  return (&surface[i]);
}

/* ---------------------------------------------------------------------- */
int
surface_sort (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Sort array of surface structures
 */
  if (count_surface > 0)
  {
    qsort (surface, (size_t) count_surface, (size_t) sizeof (struct surface),
	   surface_compare);
  }
  return (OK);
}

/* **********************************************************************
 *
 *   Routines related to structure "temperature"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
struct temperature *
temperature_bsearch (int k, int *n)
/* ---------------------------------------------------------------------- */
{
/*
 *   Binary search of array temperature to find user number k.
 *   Assumes array temperature is in sort order.
 *
 *   Input: k, user number to find.
 *   Output: n, position in array temperature of user number k.
 *   Return: pointer to temperature structure if found,
 *           NULL if not found.
 */
  void *void_ptr;
  if (count_temperature <= 0)
  {
    *n = -999;
    return (NULL);
  }
  void_ptr = (void *)
    bsearch ((char *) &k,
	     (char *) temperature,
	     (size_t) count_temperature,
	     (size_t) sizeof (struct temperature), temperature_compare_int);
  if (void_ptr == NULL)
  {
    *n = -999;
    return (NULL);
  }
  *n = (int) ((struct temperature *) void_ptr - temperature);
  return ((struct temperature *) void_ptr);
}

/* ---------------------------------------------------------------------- */
int
temperature_compare (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
  const struct temperature *temperature_ptr1, *temperature_ptr2;
  temperature_ptr1 = (const struct temperature *) ptr1;
  temperature_ptr2 = (const struct temperature *) ptr2;
  if (temperature_ptr1->n_user > temperature_ptr2->n_user)
    return (1);
  if (temperature_ptr1->n_user < temperature_ptr2->n_user)
    return (-1);
  return (0);

}

/* ---------------------------------------------------------------------- */
static int
temperature_compare_int (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
/*
 *   Compare temperature user numbers
 */
  const int *nptr1;
  const struct temperature *nptr2;

  nptr1 = (const int *) ptr1;
  nptr2 = (const struct temperature *) ptr2;
  if (*nptr1 > nptr2->n_user)
    return (1);
  if (*nptr1 < nptr2->n_user)
    return (-1);
  return (0);
}

/* ---------------------------------------------------------------------- */
int
temperature_copy (struct temperature *temperature_old_ptr,
		  struct temperature *temperature_new_ptr, int n_user_new)
/* ---------------------------------------------------------------------- */
{
/*
 *   Copies temperature data from temperature_old_ptr to new location, temperature_new_ptr.
 *   Space for the temperature_new_ptr structure must already be malloced.
 *   Space for temperature components is malloced here.
 */
  int count;
  char token[MAX_LENGTH];
/*
 *   Copy old to new
 */
  memcpy (temperature_new_ptr, temperature_old_ptr,
	  sizeof (struct temperature));
/*
 *   Store data for structure temperature
 */
  temperature_new_ptr->n_user = n_user_new;
  temperature_new_ptr->n_user_end = n_user_new;
  sprintf (token, "Temperature defined in simulation %d.", simulation);
  temperature_new_ptr->description = string_duplicate (token);
/*
 *   Count temperature components and allocate space
 */
  if (temperature_old_ptr->count_t < 0)
  {
    count = 2;
  }
  else
  {
    count = temperature_old_ptr->count_t;
  }
  temperature_new_ptr->t =
    (double *) PHRQ_malloc ((size_t) (count) * sizeof (double));
  if (temperature_new_ptr->t == NULL)
    malloc_error ();
  memcpy (temperature_new_ptr->t, temperature_old_ptr->t,
	  (size_t) (count * sizeof (double)));

  return (OK);
}

/* ---------------------------------------------------------------------- */
int
temperature_duplicate (int n_user_old, int n_user_new)
/* ---------------------------------------------------------------------- */
{
/*
 *   Checks if n_user_new exists, and frees space if it does
 *   Copies temperature[n_user_old] to old n_user_new space if
 *   found or to temperature[count_temperature] if not found.
 *   Temperature array may not be in sort order after the copy.
 */
  int n_old, n_new, sort;
  char token[MAX_LENGTH];
  int count_t;
  struct temperature *temperature_ptr_old, *temperature_ptr_new;
/*
 *   Find n_user_old in structure array temperature or make new space
 */
  if (n_user_old == n_user_new)
    return (OK);
  temperature_ptr_old = temperature_bsearch (n_user_old, &n_old);
  if (temperature_ptr_old == NULL)
  {
    sprintf (error_string, "Temperature %d not found.", n_user_old);
    error_msg (error_string, CONTINUE);
    input_error++;
    return (ERROR);
  }
/*
 *   Find n_user_old in structure array temperature or make new space
 */
  sort = FALSE;
  temperature_ptr_new = temperature_bsearch (n_user_new, &n_new);
  if (temperature_ptr_new != NULL)
  {
    temperature_free (&temperature[n_new]);
  }
  else
  {
    temperature = (struct temperature *) PHRQ_realloc (temperature,
						       (size_t)
						       (count_temperature +
							1) *
						       sizeof (struct
							       temperature));
    if (temperature == NULL)
      malloc_error ();
    if (n_user_new < temperature[count_temperature - 1].n_user)
      sort = TRUE;
    n_new = count_temperature++;
  }
/*
 *   Store data for structure temperature
 */
  count_t = temperature[n_old].count_t;
  temperature[n_new].n_user = n_user_new;
  temperature[n_new].n_user_end = n_user_new;
  sprintf (token, "Temperature defined in simulation %d.", simulation);
  temperature[n_new].description = string_duplicate (token);
  temperature[n_new].count_t = count_t;
  if (count_t < 0)
    count_t = 2;
/*
 *   Count temperature components and allocate space
 */
  temperature[n_new].t = (LDBLE *) PHRQ_malloc ((size_t) (count_t) *
						sizeof (LDBLE));
  if (temperature[n_new].t == NULL)
    malloc_error ();
/*
 *   Write temperature_comp structure for each temperature component
 */
  memcpy (temperature[n_new].t, temperature[n_old].t,
	  (size_t) (count_t) * sizeof (LDBLE));
  if (sort == TRUE)
    temperature_sort ();
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
temperature_free (struct temperature *temperature_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Free space allocated for temperature structure, does not free temperature_ptr.
 */
  if (temperature_ptr == NULL)
    return (ERROR);

  temperature_ptr->description =
    (char *) free_check_null (temperature_ptr->description);
  temperature_ptr->t = (LDBLE *) free_check_null (temperature_ptr->t);
  return (OK);
}

/* ---------------------------------------------------------------------- */
struct temperature *
temperature_search (int n_user, int *n)
/* ---------------------------------------------------------------------- */
{
/*   Linear search of the structure array "temperature" for user number n_user.
 *
 *   Arguments:
 *      n_user  input, user number.
 *      n       output, position in temperature.
 *
 *   Returns:
 *      if found, pointer to temperature structure
 *      if not found, NULL
 */
  int i;
  for (i = 0; i < count_temperature; i++)
  {
    if (n_user == temperature[i].n_user)
    {
      break;
    }
  }
  if (i >= count_temperature)
  {
    *n = -999;
    return (NULL);
  }
  *n = i;
  return (&temperature[i]);
}

/* ---------------------------------------------------------------------- */
int
temperature_sort (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Sort array of temperature structures
 */
  if (count_temperature > 0)
  {
    qsort (temperature, (size_t) count_temperature,
	   (size_t) sizeof (struct temperature), temperature_compare);
  }
  return (OK);
}

/* **********************************************************************
 *
 *   Routines related to structure "trxn"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
int
rxn_token_temp_compare (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
  const struct rxn_token_temp *rxn_token_temp_ptr1, *rxn_token_temp_ptr2;
  rxn_token_temp_ptr1 = (const struct rxn_token_temp *) ptr1;
  rxn_token_temp_ptr2 = (const struct rxn_token_temp *) ptr2;
  return (strcmp (rxn_token_temp_ptr1->name, rxn_token_temp_ptr2->name));
}

/* ---------------------------------------------------------------------- */
int
trxn_add (struct reaction *r_ptr, LDBLE coef, int combine)
/* ---------------------------------------------------------------------- */
{
/*
 *   Adds reactions together.
 *
 *   Global variable count_trxn determines which position in trxn is used.
 *      If count_trxn=0, then the equation effectively is copied into trxn.
 *      If count_trxn>0, then new equation is added to existing equation.
 *
 *   Arguments:
 *      *r_ptr         points to rxn structure to add.
 *
 *       coef          added equation is multiplied by coef.
 *       combine       if TRUE, reaction is reaction is sorted and
 *                     like terms combined.
 */
  int i;
  struct rxn_token *next_token;
/*
 *   Accumulate log k for reaction
 */
  if (count_trxn == 0)
  {
    memcpy ((void *) trxn.logk, (void *) r_ptr->logk,
	    (size_t) 7 * sizeof (LDBLE));
  }
  else
  {
    for (i = 0; i < 7; i++)
    {
      trxn.logk[i] += coef * (r_ptr->logk[i]);
    }
  }
/*
 *   Copy  equation into work space
 */
  next_token = r_ptr->token;
  while (next_token->s != NULL)
  {
    if (count_trxn + 1 >= max_trxn)
    {
      space ((void **) ((void *) &(trxn.token)), count_trxn + 1, &max_trxn,
	     sizeof (struct rxn_token_temp));
    }
    trxn.token[count_trxn].name = next_token->s->name;
    trxn.token[count_trxn].s = next_token->s;
    trxn.token[count_trxn].coef = coef * next_token->coef;
    count_trxn++;
    next_token++;
  }
  if (combine == TRUE)
    trxn_combine ();
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
trxn_add_phase (struct reaction *r_ptr, LDBLE coef, int combine)
/* ---------------------------------------------------------------------- */
{
/*
 *   Adds reactions together.
 *
 *   Global variable count_trxn determines which position in trxn is used.
 *      If count_trxn=0, then the equation effectively is copied into trxn.
 *      If count_trxn>0, then new equation is added to existing equation.
 *
 *   Arguments:
 *      *r_ptr         points to rxn structure to add.
 *
 *       coef          added equation is multiplied by coef.
 *       combine       if TRUE, reaction is reaction is sorted and
 *                     like terms combined.
 */
  int i;
  struct rxn_token *next_token;
/*
 *   Accumulate log k for reaction
 */
  if (count_trxn == 0)
  {
    memcpy ((void *) trxn.logk, (void *) r_ptr->logk,
	    (size_t) 7 * sizeof (LDBLE));
  }
  else
  {
    for (i = 0; i < 7; i++)
    {
      trxn.logk[i] += coef * (r_ptr->logk[i]);
    }
  }
/*
 *   Copy  equation into work space
 */
  next_token = r_ptr->token;
  while (next_token->s != NULL || next_token->name != NULL)
  {
    if (count_trxn + 1 >= max_trxn)
    {
      space ((void **) ((void *) &(trxn.token)), count_trxn + 1, &max_trxn,
	     sizeof (struct rxn_token_temp));
    }
    if (next_token->s != NULL)
    {
      trxn.token[count_trxn].name = next_token->s->name;
      trxn.token[count_trxn].s = next_token->s;
    }
    else
    {
      trxn.token[count_trxn].name = next_token->name;
      trxn.token[count_trxn].s = NULL;
    }
    trxn.token[count_trxn].coef = coef * next_token->coef;
    count_trxn++;
    next_token++;
  }
  if (combine == TRUE)
    trxn_combine ();
  return (OK);
}

#ifdef SKIP
/* ---------------------------------------------------------------------- */
int
trxn_combine (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Combines coefficients of tokens that are equal in temporary
 *   reaction structure, trxn.
 */
  int j, k;
/*
 *   Sort trxn species
 */
  trxn_sort ();
/*
 *   Combine trxn tokens
 */
  j = 1;
  for (k = 2; k < count_trxn; k++)
  {
    if ((j > 0) && (trxn.token[k].s == trxn.token[j].s))
    {
      trxn.token[j].coef += trxn.token[k].coef;
      if (equal (trxn.token[j].coef, 0.0, 1e-5))
	j--;
    }
    else
    {
      j++;
      if (k != j)
      {
	trxn.token[j].name = trxn.token[k].name;
	trxn.token[j].s = trxn.token[k].s;
	trxn.token[j].coef = trxn.token[k].coef;
      }
    }
  }
  count_trxn = j + 1;		/* number excluding final NULL */
  return (OK);
}
#endif
/* ---------------------------------------------------------------------- */
int
trxn_combine (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Combines coefficients of tokens that are equal in temporary
 *   reaction structure, trxn.
 */
  int j, k;
/*
 *   Sort trxn species
 */
  trxn_sort ();
/*
 *   Combine trxn tokens
 */
  j = 1;
  for (k = 2; k < count_trxn; k++)
  {
    if (trxn.token[k].s != NULL)
    {
      if ((j > 0) && (trxn.token[k].s == trxn.token[j].s))
      {
	trxn.token[j].coef += trxn.token[k].coef;
	if (equal (trxn.token[j].coef, 0.0, 1e-5))
	  j--;
      }
      else
      {
	j++;
	if (k != j)
	{
	  trxn.token[j].name = trxn.token[k].name;
	  trxn.token[j].s = trxn.token[k].s;
	  trxn.token[j].coef = trxn.token[k].coef;
	}
      }
    }
    else
    {
      if ((j > 0) && (trxn.token[k].s == trxn.token[j].s)
	  && (trxn.token[k].name == trxn.token[j].name))
      {
	trxn.token[j].coef += trxn.token[k].coef;
	if (equal (trxn.token[j].coef, 0.0, 1e-5))
	  j--;
      }
      else
      {
	j++;
	if (k != j)
	{
	  trxn.token[j].name = trxn.token[k].name;
	  trxn.token[j].s = trxn.token[k].s;
	  trxn.token[j].coef = trxn.token[k].coef;
	}
      }
    }
  }
  count_trxn = j + 1;		/* number excluding final NULL */
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
trxn_copy (struct reaction *rxn_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Copies trxn to a reaction structure.
 *
 *   Input: rxn_ptr, pointer to reaction structure to copy trxn to.
 *
 */
  int i;
/*
 *   Copy logk data
 */
  for (i = 0; i < 7; i++)
  {
    rxn_ptr->logk[i] = trxn.logk[i];
  }
/*
 *   Copy tokens
 */
  for (i = 0; i < count_trxn; i++)
  {
    rxn_ptr->token[i].s = trxn.token[i].s;
    rxn_ptr->token[i].name = trxn.token[i].name;
    rxn_ptr->token[i].coef = trxn.token[i].coef;
  }
  rxn_ptr->token[count_trxn].s = NULL;

  return (OK);
}

/* ---------------------------------------------------------------------- */
LDBLE
trxn_find_coef (const char *str, int start)
/* ---------------------------------------------------------------------- */
{
/*
 *   Finds coefficient of specified token in trxn.
 *   Input: str, token name in reaction.
 *
 *   Return: 0.0, if token not found.
 *           coefficient of token, if token found.
 */
  int i;
  LDBLE coef;

  coef = 0.0;
  for (i = start; i < count_trxn; i++)
  {
    if (strcmp (trxn.token[i].s->name, str) == 0)
    {
      coef = trxn.token[i].coef;
      break;
    }
  }
  return (coef);
}

/* ---------------------------------------------------------------------- */
int
trxn_multiply (LDBLE coef)
/* ---------------------------------------------------------------------- */
{
/*
 *   Multiplies temporary reaction, trxn,  by a constant
 *
 *   Arguments:
 *       input: coef          multiplier.
 */
  int i;
/*
 *   Multiply log k for reaction
 */
  for (i = 0; i < 7; i++)
  {
    trxn.logk[i] *= coef;
  }
/*
 *   Multiply coefficients of reaction
 */
  for (i = 0; i < count_trxn; i++)
  {
    trxn.token[i].coef *= coef;
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
trxn_print (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Prints trxn
 */
  int i;
/*
 *   Print log k for reaction
 */

  output_msg (OUTPUT_MESSAGE, "\tlog k data:\n");
  for (i = 0; i < 7; i++)
  {
    output_msg (OUTPUT_MESSAGE, "\t\t%f\n", (double) trxn.logk[i]);
  }

/*
 *   Print stoichiometry
 */
  output_msg (OUTPUT_MESSAGE, "\tReaction stoichiometry\n");
  for (i = 0; i < count_trxn; i++)
  {
    output_msg (OUTPUT_MESSAGE, "\t\t%-20s\t%10.2f\n", trxn.token[i].name,
		(double) trxn.token[i].coef);
  }
  output_msg (OUTPUT_MESSAGE, "\n");
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
trxn_reverse_k (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Changes K from dissociation to association and back
 */
  int i;
/*
 *   Accumulate log k for reaction
 */
  for (i = 0; i < 7; i++)
  {
    trxn.logk[i] = -trxn.logk[i];
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
trxn_sort (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Compare names in tokens in trxn array for sorting
 */
  if (count_trxn - 1 > 0)
  {
    qsort (&trxn.token[1],
	   (size_t) count_trxn - 1,
	   (size_t) sizeof (struct rxn_token_temp), rxn_token_temp_compare);
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
trxn_swap (const char *token)
/* ---------------------------------------------------------------------- */
{
/*
 *   Moves specified token to initial position in reaction.
 *   Input: token, token name to move to initial position.
 *
 *   Return: ERROR, if token not found.
 *           OK, if token moved to initial position.
 */
  int i, j;
  LDBLE coef;
/*
 *   Locate token
 */
  for (j = 0; j < count_trxn; j++)
  {
    if (strcmp (trxn.token[j].s->name, token) == 0)
      break;
  }
  if (j >= count_trxn)
  {
    input_error++;
    sprintf (error_string, "Could not find token in equation, %s.", token);
    error_msg (error_string, CONTINUE);
    for (i = 0; i < count_trxn; i++)
    {
      output_msg (OUTPUT_MESSAGE, "%f\t%s\t", (double) trxn.token[i].coef,
		  trxn.token[i].name);
    }
    output_msg (OUTPUT_MESSAGE, "\n");
    return (ERROR);
  }
/*
 *   Swap token to first position
 */
  trxn.token[count_trxn].name = trxn.token[0].name;
  trxn.token[count_trxn].s = trxn.token[0].s;
  trxn.token[count_trxn].coef = trxn.token[0].coef;

  trxn.token[0].name = trxn.token[j].name;
  trxn.token[0].s = trxn.token[j].s;
  trxn.token[0].coef = trxn.token[j].coef;

  trxn.token[j].name = trxn.token[count_trxn].name;
  trxn.token[j].s = trxn.token[count_trxn].s;
  trxn.token[j].coef = trxn.token[count_trxn].coef;
/*
 *   Make coefficient of token -1.0
 */
  coef = -1.0 / trxn.token[0].coef;
  trxn_multiply (coef);
  return (OK);
}

/* **********************************************************************
 *
 *   Routines related to structure "unknown"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
struct unknown *
unknown_alloc (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Allocates space to an "unknown" structure
 *      arguments: void
 *      return: pointer to an "unknown" structure
 */
  struct unknown *unknown_ptr;
/*
 *   Allocate space
 */
  unknown_ptr = (struct unknown *) PHRQ_malloc (sizeof (struct unknown));
  if (unknown_ptr == NULL)
    malloc_error ();
/*
 *   set pointers in structure to NULL
 */
  unknown_ptr->type = 0;
  unknown_ptr->moles = 0.0;
  unknown_ptr->ln_moles = 0.0;
  unknown_ptr->f = 0.0;
  unknown_ptr->sum = 0.0;
  unknown_ptr->delta = 0.0;
  unknown_ptr->la = 0.0;
  unknown_ptr->number = 0;
  unknown_ptr->description = NULL;
  unknown_ptr->master = NULL;
  unknown_ptr->phase = NULL;
  unknown_ptr->si = 0.0;
  unknown_ptr->gas_phase = NULL;
  unknown_ptr->total = NULL;
  unknown_ptr->s = NULL;
  unknown_ptr->exch_comp = NULL;
  unknown_ptr->pure_phase = NULL;
  unknown_ptr->s_s = NULL;
  unknown_ptr->s_s_comp = NULL;
  unknown_ptr->s_s_comp_number = 0;
  unknown_ptr->s_s_in = FALSE;
  unknown_ptr->surface_comp = NULL;
  unknown_ptr->related_moles = 0.0;
  unknown_ptr->potential_unknown = NULL;
  unknown_ptr->potential_unknown1 = NULL;
  unknown_ptr->potential_unknown2 = NULL;
  unknown_ptr->count_comp_unknowns = 0;
  unknown_ptr->comp_unknowns = NULL;
  unknown_ptr->phase_unknown = NULL;
  unknown_ptr->surface_charge = NULL;
  unknown_ptr->mass_water = 0.0;
  unknown_ptr->dissolve_only = FALSE;

  return (unknown_ptr);
}

/* ---------------------------------------------------------------------- */
int
unknown_delete (int i)
/* ---------------------------------------------------------------------- */
{
/*
 *   Delete unknow from list x
 */
  int j;

  unknown_free (x[i]);
  for (j = i; j < (count_unknowns); j++)
  {
    x[j] = x[j + 1];
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
unknown_free (struct unknown *unknown_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Frees space allocated to an unknown structure, frees unknown_ptr.
 */
  if (unknown_ptr == NULL)
    return (ERROR);
  unknown_ptr->master =
    (struct master **) free_check_null (unknown_ptr->master);
  if (unknown_ptr->type == SURFACE_CB)
  {
    /*
       surface_charge_free(unknown_ptr->surface_charge);
       unknown_ptr->surface_charge = (struct surface_charge *) free_check_null(unknown_ptr->surface_charge);
     */
  }
  unknown_ptr->comp_unknowns = (struct unknown **) free_check_null(unknown_ptr->comp_unknowns); 
  unknown_ptr = (struct unknown *) free_check_null (unknown_ptr);
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
system_duplicate (int i, int save_old)
/* ---------------------------------------------------------------------- */
{
  int n;

  if (solution_bsearch (i, &n, TRUE) != NULL)
    solution_duplicate (i, save_old);
  if (pp_assemblage_bsearch (i, &n) != NULL)
    pp_assemblage_duplicate (i, save_old);
  if (exchange_bsearch (i, &n) != NULL)
    exchange_duplicate (i, save_old);
  if (surface_bsearch (i, &n) != NULL)
    surface_duplicate (i, save_old);
  if (gas_phase_bsearch (i, &n) != NULL)
    gas_phase_duplicate (i, save_old);
  if (kinetics_bsearch (i, &n) != NULL)
    kinetics_duplicate (i, save_old);
  if (s_s_assemblage_bsearch (i, &n) != NULL)
    s_s_assemblage_duplicate (i, save_old);
  return (OK);
}

/* ---------------------------------------------------------------------- */
struct logk *
logk_store (char *name, int replace_if_found)
/* ---------------------------------------------------------------------- */
{
/*
 *   Function locates the string "name" in the hash table for logk.
 *
 *   Pointer to a logk structure is always returned.
 *
 *   If the string is not found, a new entry is made in the hash table. Pointer to
 *      the new structure is returned.
 *   If "name" is found and replace is true, pointers in old logk structure
 *      are freed and replaced with additional input.
 *   If "name" is found and replace is false, the old logk structure is not
 *      modified and a pointer to it is returned.
 *
 *   Arguments:
 *      name    input, character string to be found in "logk".
 *      replace_if_found input, TRUE means reinitialize logk structure if found
 *                     FALSE means just return pointer if found.
 *
 *   Returns:
 *      pointer to logk structure "logk" where "name" can be found.
 */
  int n;
  struct logk *logk_ptr;
  ENTRY item, *found_item;
/*
 *   Search list
 */
  str_tolower (name);
  item.key = name;
  item.data = NULL;
  found_item = hsearch_multi (logk_hash_table, item, FIND);

  if (found_item != NULL && replace_if_found == FALSE)
  {
    logk_ptr = (struct logk *) (found_item->data);
    return (logk_ptr);
  }
  else if (found_item != NULL && replace_if_found == TRUE)
  {
    logk_ptr = (struct logk *) (found_item->data);
    logk_init (logk_ptr);
  }
  else
  {
    n = count_logk++;
    /* make sure there is space in s */
    if (count_logk >= max_logk)
    {
      space ((void **) ((void *) &logk), count_logk, &max_logk,
	     sizeof (struct logk *));
    }
    /* Make new logk structure */
    logk[n] = logk_alloc ();
    logk_ptr = logk[n];
  }
  /* set name and z in pointer in logk structure */
  logk_ptr->name = string_hsave (name);
/*
 *   Update hash table
 */
  item.key = logk_ptr->name;
  item.data = (void *) logk_ptr;
  found_item = hsearch_multi (logk_hash_table, item, ENTER);
  if (found_item == NULL)
  {
    sprintf (error_string, "Hash table error in logk_store.");
    error_msg (error_string, CONTINUE);
  }

  return (logk_ptr);
}

/* ---------------------------------------------------------------------- */
struct logk *
logk_alloc (void)
/* ---------------------------------------------------------------------- */
/*
 *   Allocates space to a logk structure, initializes
 *      arguments: void
 *      return: pointer to a logk structure
 */
{
  struct logk *logk_ptr;
  logk_ptr = (struct logk *) PHRQ_malloc (sizeof (struct logk));
  if (logk_ptr == NULL)
    malloc_error ();
/*
 *   set pointers in structure to NULL, variables to zero
 */
  logk_init (logk_ptr);

  return (logk_ptr);
}

/* ---------------------------------------------------------------------- */
static int
logk_init (struct logk *logk_ptr)
/* ---------------------------------------------------------------------- */
/*
 *      return: pointer to a logk structure
 */
{
  int i;
/*
 *   set pointers in structure to NULL
 */
  logk_ptr->name = NULL;
/*
 *   set varibles = 0
 */
  logk_ptr->lk = 0.0;
  for (i = 0; i < 8; i++)
  {
    logk_ptr->log_k[i] = 0.0;
    logk_ptr->log_k_original[i] = 0.0;
  }
  logk_ptr->count_add_logk = 0;
  logk_ptr->add_logk = NULL;
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
logk_copy2orig (struct logk *logk_ptr)
/* ---------------------------------------------------------------------- */
/*
 *   Copies log k data to logk_original
 */
{
  int i;
  for (i = 0; i < 7; i++)
  {
    logk_ptr->log_k_original[i] = logk_ptr->log_k[i];
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
struct logk *
logk_search (char *name)
/* ---------------------------------------------------------------------- */
{
/*
 *   Function locates the string "name" in the hash table for logk.
 *
 *   Arguments:
 *      name    input, character string to be found in "logk".
 *
 *   Returns:
 *      pointer to logk structure "logk" where "name" can be found.
 *      or NULL if not found.
 */
  struct logk *logk_ptr;
  ENTRY item, *found_item;
/*
 *   Search list
 */
  str_tolower (name);
  item.key = name;
  item.data = NULL;
  found_item = hsearch_multi (logk_hash_table, item, FIND);

  if (found_item != NULL)
  {
    logk_ptr = (struct logk *) (found_item->data);
    return (logk_ptr);
  }
  return (NULL);
}

/* ---------------------------------------------------------------------- */
int
entity_exists (char *name, int n_user)
/* ---------------------------------------------------------------------- */
{
/*
 *   Reads solution,                   0 Solution
 *         reaction,                   1 Reaction
 *         exchange,                   2 Exchange
 *         surface,                    3 Surface
 *         gas_phase,                  4 Gas_phase
 *         equilibrium_phases,         5 Pure_phase
 *         solid_solution,             6 Ss_phase
 *         kinetics,                   7 Kinetics
 *         mix,                        8 Mix
 *         reaction_temperature        9 Temperature
 *         unknown                     10 UnKnown
 */
  int i, return_value;
  char *ptr;
  char token[MAX_LENGTH];
  enum entity_type type;
/*
 *   Read keyword
 */
  ptr = name;
  copy_token (token, &ptr, &i);
  type = get_entity_enum (token);
  return_value = TRUE;
  switch (type)
  {
  case UnKnown:
    warning_msg
      ("EXISTS expecting keyword solution, mix, kinetics, reaction, reaction_temperature, equilibrium_phases, exchange, surface, gas_phase, or solid_solutions.");
    return_value = 2;
    break;
  case Solution:		/* Solution */
    if (solution_bsearch (n_user, &i, FALSE) == NULL)
    {
      return_value = FALSE;
    }
    break;
  case Pure_phase:		/* Pure phases */
    if (pp_assemblage_bsearch (n_user, &i) == NULL)
    {
      return_value = FALSE;
    }
    break;
  case Reaction:		/* Reaction */
    if (irrev_bsearch (n_user, &i) == NULL)
    {
      return_value = FALSE;
    }
    break;
  case Mix:			/* Mix */
    if (mix_bsearch (n_user, &i) == NULL)
    {
      return_value = FALSE;
    }
    break;
  case Exchange:		/* Ex */
    if (exchange_bsearch (n_user, &i) == NULL)
    {
      return_value = FALSE;
    }
    break;
  case Surface:		/* Surface */
    if (surface_bsearch (n_user, &i) == NULL)
    {
      return_value = FALSE;
    }
    break;
  case Temperature:		/* Temperature */
    if (temperature_bsearch (n_user, &i) == NULL)
    {
      return_value = FALSE;
    }
    break;
  case Gas_phase:		/* Gas */
    if (gas_phase_bsearch (n_user, &i) == NULL)
    {
      return_value = FALSE;
    }
    break;
  case Kinetics:		/* Kinetics */
    if (kinetics_bsearch (n_user, &i) == NULL)
    {
      return_value = FALSE;
    }
    break;
  case Ss_phase:		/* solid_solutions */
    if (s_s_assemblage_bsearch (n_user, &i) == NULL)
    {
      return_value = FALSE;
    }
    break;
  }
  return (return_value);
}

/* ---------------------------------------------------------------------- */
enum entity_type
get_entity_enum (char *name)
/* ---------------------------------------------------------------------- */
{
/*
 *   Reads solution,                   0 Solution
 *         reaction,                   1 Reaction
 *         exchange,                   2 Exchange
 *         surface,                    3 Surface
 *         gas_phase,                  4 Gas_phase
 *         equilibrium_phases,         5 Pure_phase
 *         solid_solution,             6 Ss_phase
 *         kinetics,                   7 Kinetics
 *         mix,                        8 Mix
 *         reaction_temperature        9 Temperature
 *         unknown                     10 UnKnown
 *
 */
  int i;
  char *ptr;
  char token[MAX_LENGTH];
/*
 *   Read keyword
 */
  ptr = name;
  copy_token (token, &ptr, &i);
  check_key (token);
  switch (next_keyword)
  {
  case -1:			/* Have not read line with keyword */
  case 0:			/* End encountered */
  case 1:			/* EOF encountered */
  case 2:			/* Read aqueous model */
  case 3:			/* Read master species */
  case 5:
  case 9:
  case 10:
  case 11:
  case 12:
  case 14:
  case 15:
  case 18:
  case 20:
  case 21:
  case 22:
  case 23:
  case 24:
  case 25:
  case 30:
  case 31:
  case 32:
  case 34:
  case 35:
  case 38:
  case 39:
    warning_msg
      ("EXISTS expecting keyword solution, mix, kinetics, reaction, reaction_temperature, equilibrium_phases, exchange, surface, gas_phase, or solid_solutions.");
    return (UnKnown);
    break;
  case 4:			/* Solution */
    return (Solution);
    break;
  case 6:			/* Pure phases */
  case 26:
  case 27:
  case 28:
  case 29:
    return (Pure_phase);
    break;
  case 7:			/* Reaction */
    return (Reaction);
    break;
  case 8:			/* Mix */
    return (Mix);
    break;
  case 13:			/* Ex */
    return (Exchange);
    break;
  case 16:			/* Surface */
    return (Surface);
    break;
  case 17:			/* Temperature */
    return (Temperature);
    break;
  case 19:			/* Gas */
    return (Gas_phase);
    break;
  case 33:			/* Kinetics */
    return (Kinetics);
    break;
  case 40:			/* solid_solutions */
  case 41:			/* solid_solution */
    return (Ss_phase);
    break;
  }
  return (UnKnown);
}

/*
 * copier routines
 */
/* ---------------------------------------------------------------------- */
int
copier_add (struct copier *copier_ptr, int n_user, int start, int end)
/* ---------------------------------------------------------------------- */
/*
 *   add new set of copy instructions
 */
{

  if (copier_ptr->count >= copier_ptr->max)
  {
    copier_ptr->max = copier_ptr->count * 2;
    copier_ptr->n_user =
      (int *) PHRQ_realloc (copier_ptr->n_user,
			    (size_t) (copier_ptr->max * sizeof (int)));
    if (copier_ptr->n_user == NULL)
      malloc_error ();
    copier_ptr->start =
      (int *) PHRQ_realloc (copier_ptr->start,
			    (size_t) (copier_ptr->max * sizeof (int)));
    if (copier_ptr->start == NULL)
      malloc_error ();
    copier_ptr->end =
      (int *) PHRQ_realloc (copier_ptr->end,
			    (size_t) (copier_ptr->max * sizeof (int)));
    if (copier_ptr->end == NULL)
      malloc_error ();
  }
  copier_ptr->n_user[copier_ptr->count] = n_user;
  copier_ptr->start[copier_ptr->count] = start;
  copier_ptr->end[copier_ptr->count] = end;
  copier_ptr->count++;
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
copier_free (struct copier *copier_ptr)
/* ---------------------------------------------------------------------- */
/*
 *   initialize copier structure
 */
{

  copier_ptr->n_user = (int *) free_check_null (copier_ptr->n_user);
  copier_ptr->start = (int *) free_check_null (copier_ptr->start);
  copier_ptr->end = (int *) free_check_null (copier_ptr->end);
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
copier_init (struct copier *copier_ptr)
/* ---------------------------------------------------------------------- */
/*
 *   initialize copier structure
 */
{

  copier_ptr->count = 0;
  copier_ptr->max = 10;
  copier_ptr->n_user =
    (int *) PHRQ_malloc ((size_t) (copier_ptr->max * sizeof (int)));
  copier_ptr->start =
    (int *) PHRQ_malloc ((size_t) (copier_ptr->max * sizeof (int)));
  copier_ptr->end =
    (int *) PHRQ_malloc ((size_t) (copier_ptr->max * sizeof (int)));
  return (OK);
}
