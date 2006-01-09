#define EXTERNAL extern
#define MAINSUBS
#include "global.h"
#include "phqalloc.h"
#include "output.h"
#include "phrqproto.h"
#include "input.h"

extern void test_classes (void);
/*#define PHREEQC_XML*/
#ifdef PHREEQC_XML
#include "SAXPhreeqc.h"
extern void SAX_cleanup (void);
#endif

static char const svnid[] =
  "$Id: mainsubs.c 4 2009-04-21 17:29:29Z delucia $";

#if defined(WINDOWS) || defined(_WINDOWS)
#include <windows.h>
#endif

static int copy_use (int i);
static int set_use (void);
static int xexchange_save (int n_user);
static int xgas_save (int n_user);
static int xpp_assemblage_save (int n_user);
static int xs_s_assemblage_save (int n_user);
static int xsurface_save (int n_user);

#ifdef PHREEQ98
extern int phreeq98_debug;
extern int AddSeries, connect_simulations;
#endif
/* ---------------------------------------------------------------------- */
void
initialize (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Initialize global variables
 */
  ENTRY item, *found_item;
  int i;
  struct logk *logk_ptr;
  char token[MAX_LENGTH];
  if (svnid == NULL)
    fprintf (stderr, " ");

  moles_per_kilogram_string = string_duplicate ("Mol/kgw");
  pe_string = string_duplicate ("pe");

  debug_model = FALSE;
  debug_prep = FALSE;
  debug_set = FALSE;
  debug_diffuse_layer = FALSE;
  debug_inverse = FALSE;
  itmax = 100;
#ifdef USE_LONG_LDBLE
  /* from float.h, sets tolerance for cl1 routine */
  ineq_tol = pow ((long double) 10, (long double) -LDBL_DIG);
#else
  ineq_tol = pow ((double) 10, (double) -DBL_DIG);
#endif
  convergence_tolerance = 1e-8;
#ifdef USE_LONG_LDBLE
  /* from float.h, sets tolerance for cl1 routine */
  inv_tol_default = pow ((long double) 10, (long double) -LDBL_DIG + 5);
#else
  inv_tol_default = pow ((double) 10, (double) -DBL_DIG + 5);
#endif
  step_size = 100.;
  pe_step_size = 10.;
  pp_scale = 1.0;
  pp_column_scale = 1.0;
  diagonal_scale = FALSE;
  censor = 0.0;
  mass_water_switch = FALSE;
  /*      mass_water_switch = TRUE; */
  delay_mass_water = FALSE;
  incremental_reactions = FALSE;
  aqueous_only = 0;
  negative_concentrations = FALSE;

  LOG_10 = log (10.0);
  /* Use space for all memory allocation */
  max_solution = MAX_SOLUTION;
  max_pp_assemblage = MAX_PP_ASSEMBLAGE;
  max_exchange = MAX_PP_ASSEMBLAGE;
  max_surface = MAX_PP_ASSEMBLAGE;
  max_gas_phase = MAX_PP_ASSEMBLAGE;
  max_kinetics = MAX_PP_ASSEMBLAGE;
  max_s_s_assemblage = MAX_PP_ASSEMBLAGE;

  max_elements = MAX_ELEMENTS;
  max_elts = MAX_ELTS;
  max_line = MAX_LINE;
  max_master = MAX_MASTER;
  max_mb_unknowns = MAX_TRXN;
  max_phases = MAX_PHASES;
  max_s = MAX_S;
  max_strings = MAX_STRINGS;
  max_trxn = MAX_TRXN;
  max_logk = MAX_S;
  max_master_isotope = MAX_ELTS;

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
/*
 *   Initialize advection
 */
  count_ad_cells = 1;
  count_ad_shifts = 1;
  print_ad_modulus = 1;
  punch_ad_modulus = 1;
  advection_punch = (int *) PHRQ_malloc (sizeof (int));
  if (advection_punch == NULL)
    malloc_error ();
  advection_punch[0] = TRUE;
  advection_kin_time = 0.0;
  advection_kin_time_defined = FALSE;
  advection_print = (int *) PHRQ_malloc (sizeof (int));
  if (advection_print == NULL)
    malloc_error ();
  advection_print[0] = TRUE;
  advection_warnings = TRUE;
/*
 *   Initialize transport
 */
  count_cells = 1;
  count_shifts = 1;
  ishift = 1;
  bcon_first = bcon_last = 3;
  diffc = 0.3e-9;
  simul_tr = 0;
  tempr = 2.0;
  heat_diffc = -0.1;
  timest = 0.0;
  multi_Dflag = FALSE;
/* !!!!        count_stag = 0; */

  mcd_substeps = 1.0;
  print_modulus = 1;
  punch_modulus = 1;
  dump_modulus = 0;
  dump_in = FALSE;
  transport_warnings = TRUE;
/*
 *   Allocate space
 */
  space ((void **) ((void *) &pp_assemblage), INIT, &max_pp_assemblage,
	 sizeof (struct pp_assemblage));

  space ((void **) ((void *) &exchange), INIT, &max_exchange,
	 sizeof (struct exchange));

  space ((void **) ((void *) &surface), INIT, &max_surface,
	 sizeof (struct surface));

  space ((void **) ((void *) &gas_phase), INIT, &max_gas_phase,
	 sizeof (struct gas_phase));

  space ((void **) ((void *) &kinetics), INIT, &max_kinetics,
	 sizeof (struct kinetics));

  space ((void **) ((void *) &s_s_assemblage), INIT, &max_s_s_assemblage,
	 sizeof (struct s_s_assemblage));

  space ((void **) ((void *) &cell_data), INIT, &count_cells,
	 sizeof (struct cell_data));

  space ((void **) ((void *) &elements), INIT, &max_elements,
	 sizeof (struct element *));

  space ((void **) ((void *) &elt_list), INIT, &max_elts,
	 sizeof (struct elt_list));


  inverse = (struct inverse *) PHRQ_malloc ((size_t) sizeof (struct inverse));
  if (inverse == NULL)
    malloc_error ();
  count_inverse = 0;

  irrev = (struct irrev *) PHRQ_malloc ((size_t) sizeof (struct irrev));
  if (irrev == NULL)
    malloc_error ();

  space ((void **) ((void *) &line), INIT, &max_line, sizeof (char));

  space ((void **) ((void *) &line_save), INIT, &max_line, sizeof (char));

  space ((void **) ((void *) &master), INIT, &max_master,
	 sizeof (struct master *));

  space ((void **) ((void *) &mb_unknowns), INIT, &max_mb_unknowns,
	 sizeof (struct unknown_list));

  mix = (struct mix *) PHRQ_malloc ((size_t) sizeof (struct mix));
  if (mix == NULL)
    malloc_error ();
  count_mix = 0;
/* !!!! */
  stag_data = (struct stag_data *) PHRQ_calloc (1, sizeof (struct stag_data));
  if (stag_data == NULL)
    malloc_error ();
  stag_data->count_stag = 0;
  stag_data->exch_f = 0;
  stag_data->th_m = 0;
  stag_data->th_im = 0;

  space ((void **) ((void *) &phases), INIT, &max_phases,
	 sizeof (struct phase *));

  space ((void **) ((void *) &trxn.token), INIT, &max_trxn,
	 sizeof (struct rxn_token_temp));

  space ((void **) ((void *) &s), INIT, &max_s, sizeof (struct species *));

  space ((void **) ((void *) &logk), INIT, &max_logk, sizeof (struct logk *));

  space ((void **) ((void *) &master_isotope), INIT, &max_master_isotope,
	 sizeof (struct master_isotope *));

  solution =
    (struct solution **) PHRQ_malloc ((size_t) MAX_SOLUTION *
				      sizeof (struct solution *));
  if (solution == NULL)
    malloc_error ();

  temperature =
    (struct temperature *) PHRQ_malloc ((size_t) sizeof (struct temperature));
  if (temperature == NULL)
    malloc_error ();

  title_x = NULL;
  pe_x = NULL;
  description_x = NULL;
  units_x = NULL;
  s_x = NULL;
/* SRC ADDED */
  sum_mb1 = NULL;
  sum_mb2 = NULL;
  sum_jacob0 = NULL;
  sum_jacob1 = NULL;
  sum_jacob2 = NULL;
  sum_delta = NULL;
/* SRC ADDED */
  isotopes_x =
    (struct isotope *) PHRQ_malloc ((size_t) sizeof (struct isotope));
  if (isotopes_x == NULL)
    malloc_error ();
  x = NULL;
  max_unknowns = 0;

  array = NULL;
  delta = NULL;
  residual = NULL;
  s_h2o = NULL;
  s_hplus = NULL;
  s_h3oplus = NULL;
  s_eminus = NULL;
  s_co3 = NULL;
  s_h2 = NULL;
  s_o2 = NULL;

  hcreate_multi ((unsigned) max_logk, &logk_hash_table);
  hcreate_multi ((unsigned) max_master_isotope, &master_isotope_hash_table);
/*
 *   Create hash table for strings
 */
  hcreate_multi ((unsigned) max_strings, &strings_hash_table);
  hcreate_multi ((unsigned) max_elements, &elements_hash_table);
  hcreate_multi ((unsigned) max_s, &species_hash_table);
  hcreate_multi ((unsigned) max_phases, &phases_hash_table);
  hcreate_multi ((unsigned) 2 * NKEYS, &keyword_hash_table);
/*
 *  Initialize use pointers
 */
  use.solution_in = FALSE;
  use.pp_assemblage_in = FALSE;
  use.mix_in = FALSE;
  use.irrev_in = FALSE;
/*
 *   Initialize punch
 */
  punch.in = FALSE;
  punch.count_totals = 0;
  punch.totals =
    (struct name_master *) PHRQ_malloc (sizeof (struct name_master));
  if (punch.totals == NULL)
    malloc_error ();
  punch.count_molalities = 0;
  punch.molalities =
    (struct name_species *) PHRQ_malloc (sizeof (struct name_species));
  if (punch.molalities == NULL)
    malloc_error ();
  punch.count_activities = 0;
  punch.activities =
    (struct name_species *) PHRQ_malloc (sizeof (struct name_species));
  if (punch.activities == NULL)
    malloc_error ();
  punch.count_pure_phases = 0;
  punch.pure_phases =
    (struct name_phase *) PHRQ_malloc (sizeof (struct name_phase));
  if (punch.pure_phases == NULL)
    malloc_error ();
  punch.count_si = 0;
  punch.si = (struct name_phase *) PHRQ_malloc (sizeof (struct name_phase));
  if (punch.si == NULL)
    malloc_error ();
  punch.count_gases = 0;
  punch.gases =
    (struct name_phase *) PHRQ_malloc (sizeof (struct name_phase));
  if (punch.gases == NULL)
    malloc_error ();
  punch.count_s_s = 0;
  punch.s_s = (struct name_phase *) PHRQ_malloc (sizeof (struct name_phase));
  if (punch.s_s == NULL)
    malloc_error ();

  punch.count_kinetics = 0;
  punch.kinetics =
    (struct name_phase *) PHRQ_malloc (sizeof (struct name_phase));
  if (punch.kinetics == NULL)
    malloc_error ();

  punch.count_isotopes = 0;
  punch.isotopes =
    (struct name_master *) PHRQ_malloc (sizeof (struct name_master));
  if (punch.isotopes == NULL)
    malloc_error ();

  punch.count_calculate_values = 0;
  punch.calculate_values =
    (struct name_master *) PHRQ_malloc (sizeof (struct name_master));
  if (punch.calculate_values == NULL)
    malloc_error ();

  count_save_values = 0;
  save_values =
    (struct save_values *) PHRQ_malloc (sizeof (struct save_values));
  if (save_values == NULL)
    malloc_error ();

  punch.inverse = TRUE;

  punch.sim = TRUE;
  punch.state = TRUE;
  punch.soln = TRUE;
  punch.dist = TRUE;
  punch.time = TRUE;
  punch.step = TRUE;
  punch.rxn = FALSE;
  punch.temp = FALSE;
  punch.ph = TRUE;
  punch.pe = TRUE;
  punch.alk = FALSE;
  punch.mu = FALSE;
  punch.water = FALSE;
  punch.high_precision = FALSE;
  punch.user_punch = TRUE;
  punch.charge_balance = FALSE;
  punch.percent_error = FALSE;
/*
 *   last model
 */
  last_model.exchange = NULL;
  last_model.gas_phase = NULL;
  last_model.s_s_assemblage = NULL;
  last_model.kinetics = NULL;
  last_model.pp_assemblage = NULL;
  last_model.add_formula = NULL;
  last_model.si = NULL;
  last_model.surface_comp = NULL;
  last_model.surface_charge = NULL;
/*
 *   Update hash table
 */
  keyword_hash =
    (struct key *) PHRQ_malloc ((size_t) NKEYS * sizeof (struct key));
  if (keyword_hash == NULL)
    malloc_error ();
  for (i = 0; i < NKEYS; i++)
  {
    keyword_hash[i].name = string_hsave (keyword[i].name);
    keyword_hash[i].keycount = i;
    item.key = keyword_hash[i].name;
    item.data = (void *) &keyword_hash[i];
    found_item = hsearch_multi (keyword_hash_table, item, ENTER);
    if (found_item == NULL)
    {
      sprintf (error_string, "Hash table error in keyword initialization.");
      error_msg (error_string, STOP);
    }
  }
/*
 *   rates
 */
  rates = (struct rate *) PHRQ_malloc (sizeof (struct rate));
  if (rates == NULL)
    malloc_error ();
  count_rates = 0;
  initial_total_time = 0;
  rate_m = 0;
  rate_m0 = 0;
  rate_p = NULL;
  rate_time = 0;
  rate_sim_time_start = 0;
  rate_sim_time_end = 0;
  rate_sim_time = 0;
  rate_moles = 0;
  initial_total_time = 0;

/*
 *   user_print, user_punch
 */
  user_print = (struct rate *) PHRQ_malloc ((size_t) sizeof (struct rate));
  if (user_print == NULL)
    malloc_error ();
  user_print->commands = NULL;
  user_print->linebase = NULL;
  user_print->varbase = NULL;
  user_print->loopbase = NULL;
  user_punch = (struct rate *) PHRQ_malloc ((size_t) sizeof (struct rate));
  if (user_punch == NULL)
    malloc_error ();
  user_punch->commands = NULL;
  user_punch->linebase = NULL;
  user_punch->varbase = NULL;
  user_punch->loopbase = NULL;
  user_punch_headings = (char **) PHRQ_malloc (sizeof (char *));
  if (user_punch_headings == NULL)
    malloc_error ();
  user_punch_count_headings = 0;
#ifdef PHREEQ98
/*
 *   user_graph
 */
  user_graph = PHRQ_malloc ((size_t) sizeof (struct rate));
  if (user_graph == NULL)
    malloc_error ();
  user_graph->commands = NULL;
  user_graph->linebase = NULL;
  user_graph->varbase = NULL;
  user_graph->loopbase = NULL;
  user_graph_headings = (char **) PHRQ_malloc (sizeof (char *));
  if (user_graph_headings == NULL)
    malloc_error ();
  user_graph_count_headings = 0;
#endif
  /*
     Initialize llnl aqueous model parameters
   */
  llnl_temp = (LDBLE *) PHRQ_malloc (sizeof (LDBLE));
  if (llnl_temp == NULL)
    malloc_error ();
  llnl_count_temp = 0;
  llnl_adh = (LDBLE *) PHRQ_malloc (sizeof (LDBLE));
  if (llnl_adh == NULL)
    malloc_error ();
  llnl_count_adh = 0;
  llnl_bdh = (LDBLE *) PHRQ_malloc (sizeof (LDBLE));
  if (llnl_bdh == NULL)
    malloc_error ();
  llnl_count_bdh = 0;
  llnl_bdot = (LDBLE *) PHRQ_malloc (sizeof (LDBLE));
  if (llnl_bdot == NULL)
    malloc_error ();
  llnl_count_bdot = 0;
  llnl_co2_coefs = (LDBLE *) PHRQ_malloc (sizeof (LDBLE));
  if (llnl_co2_coefs == NULL)
    malloc_error ();
  llnl_count_co2_coefs = 0;
/*
 *
 */
  cmd_initialize ();
  change_surf =
    (struct Change_Surf *)
    PHRQ_malloc ((size_t) (2 * sizeof (struct Change_Surf)));
  if (change_surf == NULL)
    malloc_error ();
  change_surf[0].cell_no = -99;
  change_surf[0].next = TRUE;
  change_surf[1].cell_no = -99;
  change_surf[1].next = FALSE;
  change_surf_count = 0;


#if defined(WINDOWS) || defined(_WINDOWS)
  /* SRC pr.status = FALSE; */
#endif

  // MDL: libphree can choose wheter or not print results on file
  if (Pfileprint == FALSE) 
    {
      /* Initialize print here, not in global.h */
      pr.all = FALSE;
      pr.initial_solutions = FALSE;
      pr.initial_exchangers = FALSE;
      pr.reactions = FALSE;
      pr.gas_phase = FALSE;
      pr.s_s_assemblage = FALSE;
      pr.pp_assemblage = FALSE;
      pr.surface = FALSE;
      pr.exchange = FALSE;
      pr.kinetics = FALSE;
      pr.totals = FALSE;
      pr.eh = FALSE;
      pr.species = FALSE;
      pr.saturation_indices = FALSE;
      pr.irrev = FALSE;
      pr.mix = FALSE;
      pr.reaction = FALSE;
      pr.use = FALSE;
      pr.logfile = FALSE;
      pr.punch = FALSE;
      pr.status = FALSE;
      pr.inverse = FALSE;
      pr.dump = FALSE;
      pr.user_print = FALSE;
      pr.headings = FALSE;
      pr.user_graph = FALSE;
      pr.echo_input = FALSE;
      count_warnings = 0;
      pr.warnings = 100;
      pr.initial_isotopes = FALSE;
      pr.isotope_ratios = FALSE;
      pr.isotope_alphas = FALSE;
      pr.hdf = FALSE;
      pr.alkalinity = FALSE;
    }  
  else
    {
      // MDL: this is the original initialization
      /* Initialize print here, not in global.h */
      pr.all = TRUE;
      pr.initial_solutions = TRUE;
      pr.initial_exchangers = TRUE;
      pr.reactions = TRUE;
      pr.gas_phase = TRUE;
      pr.s_s_assemblage = TRUE;
      pr.pp_assemblage = TRUE;
      pr.surface = TRUE;
      pr.exchange = TRUE;
      pr.kinetics = TRUE;
      pr.totals = TRUE;
      pr.eh = TRUE;
      pr.species = TRUE;
      pr.saturation_indices = TRUE;
      pr.irrev = TRUE;
      pr.mix = TRUE;
      pr.reaction = TRUE;
      pr.use = TRUE;
      pr.logfile = FALSE;
      pr.punch = TRUE;
      if (phast == TRUE)
	{
	  pr.status = FALSE;
	}
      else
	{
	  pr.status = TRUE;
	}
      pr.inverse = TRUE;
      pr.dump = TRUE;
      pr.user_print = TRUE;
      pr.headings = TRUE;
      pr.user_graph = TRUE;
      pr.echo_input = TRUE;
      count_warnings = 0;
      pr.warnings = 100;
      pr.initial_isotopes = TRUE;
      pr.isotope_ratios = TRUE;
      pr.isotope_alphas = TRUE;
      pr.hdf = FALSE;
      pr.alkalinity = FALSE;
    }

  species_list = NULL;

  user_database = NULL;
  first_read_input = TRUE;
  have_punch_name = FALSE;
  selected_output_file_name = NULL;
  dump_file_name = NULL;

#ifdef PHREEQCI_GUI
  g_spread_sheet.heading = NULL;
  g_spread_sheet.units = NULL;
  g_spread_sheet.count_rows = 0;
  g_spread_sheet.rows = NULL;
  g_spread_sheet.defaults.units = NULL;
  g_spread_sheet.defaults.count_iso = 0;
  g_spread_sheet.defaults.iso = NULL;
#endif
  /* calculate_value */
  max_calculate_value = MAX_ELTS;
  count_calculate_value = 0;
  space ((void **) ((void *) &calculate_value), INIT, &max_calculate_value,
	 sizeof (struct calculate_value *));
  hcreate_multi ((unsigned) max_calculate_value, &calculate_value_hash_table);

  /* isotope_ratio */
  max_isotope_ratio = MAX_ELTS;
  count_isotope_ratio = 0;
  space ((void **) ((void *) &isotope_ratio), INIT, &max_isotope_ratio,
	 sizeof (struct isotope_ratio *));
  hcreate_multi ((unsigned) max_isotope_ratio, &isotope_ratio_hash_table);

  /* isotope_value */
  max_isotope_alpha = MAX_ELTS;
  count_isotope_alpha = 0;
  space ((void **) ((void *) &isotope_alpha), INIT, &max_isotope_alpha,
	 sizeof (struct isotope_alpha *));
  hcreate_multi ((unsigned) max_isotope_alpha, &isotope_alpha_hash_table);

  /* 
   * define constant named log_k
   */
  strcpy (token, "XconstantX");
  logk_ptr = logk_store (token, TRUE);
  read_log_k_only ("1.0", &logk_ptr->log_k[0]);

  phreeqc_mpi_myself = 0;

  copier_init (&copy_solution);
  copier_init (&copy_pp_assemblage);
  copier_init (&copy_exchange);
  copier_init (&copy_surface);
  copier_init (&copy_s_s_assemblage);
  copier_init (&copy_gas_phase);
  copier_init (&copy_kinetics);
  copier_init (&copy_mix);
  copier_init (&copy_irrev);
  copier_init (&copy_temperature);

  set_forward_output_to_log (FALSE);
  simulation = 0;
  /*
   *  cvode
   */
  cvode_init ();
  /*
   *  Pitzer
   */
  pitzer_init ();
  /*
   * to facilitate debuging
   */
  dbg_use = &use;
  dbg_solution = solution;
  dbg_exchange = exchange;
  dbg_surface = surface;
  dbg_pp_assemblage = pp_assemblage;
  dbg_kinetics = kinetics;
  dbg_irrev = irrev;
  dbg_mix = mix;
  dbg_master = master;
  calculating_deriv = FALSE;
  numerical_deriv = FALSE;

  zeros = (LDBLE *) PHRQ_malloc(sizeof(LDBLE));
  if (zeros == NULL) malloc_error();
  zeros[0] = 0.0;
  zeros_max = 1;

  pore_volume = 0;
  return;
}

/* ---------------------------------------------------------------------- */
static int
set_use (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Structure "use" has list of solution, ex, surf, pp_assemblage, 
 *   gas_phase and solid solution to use in current calculations, 
 *   also mix, irrev, and temp. 
 *   This routine searches for the user numbers in each list
 *   (solution, ex, ...) and sets a pointer in structure use
 */

/*
 *   Initial solution case
 */
  use.pp_assemblage_ptr = NULL;
  use.mix_ptr = NULL;
  use.irrev_ptr = NULL;
  use.exchange_ptr = NULL;
  use.kinetics_ptr = NULL;
  use.surface_ptr = NULL;
  use.temperature_ptr = NULL;
  use.gas_phase_ptr = NULL;
  use.s_s_assemblage_ptr = NULL;

  if (state < REACTION)
  {
    return (OK);
  }
/*
 *   Reaction case
 */
  if (use.pp_assemblage_in == FALSE &&
      use.irrev_in == FALSE &&
      use.mix_in == FALSE &&
      use.exchange_in == FALSE &&
      use.kinetics_in == FALSE &&
      use.surface_in == FALSE &&
      use.temperature_in == FALSE &&
      use.gas_phase_in == FALSE && use.s_s_assemblage_in == FALSE)
  {
    return (FALSE);
  }
  if (use.solution_in == FALSE && use.mix_in == FALSE)
    return (FALSE);
/*
 *   Find solution
 */
  if (use.solution_in == TRUE)
  {
    use.solution_ptr =
      solution_bsearch (use.n_solution_user, &use.n_solution, FALSE);
    if (use.solution_ptr == NULL)
    {
      sprintf (error_string, "Solution %d not found.", use.n_solution_user);
      error_msg (error_string, STOP);
    }
  }
/*
 *   Find mixture
 */
  if (use.mix_in == TRUE)
  {
    use.mix_ptr = mix_bsearch (use.n_mix_user, &use.n_mix);
    use.n_mix_user_orig = use.n_mix_user;
    if (use.mix_ptr == NULL)
    {
      sprintf (error_string, "Mixture %d not found.", use.n_mix_user);
      error_msg (error_string, STOP);
    }
  }
  else
  {
    use.mix_ptr = NULL;
  }
/*
 *   Find pure phase assemblage
 */
  if (use.pp_assemblage_in == TRUE)
  {
    use.pp_assemblage_ptr =
      pp_assemblage_bsearch (use.n_pp_assemblage_user, &use.n_pp_assemblage);
    if (use.pp_assemblage_ptr == NULL)
    {
      sprintf (error_string, "Pure phase assemblage %d not found.",
	       use.n_pp_assemblage_user);
      error_msg (error_string, STOP);
    }
  }
  else
  {
    use.pp_assemblage_ptr = NULL;
  }
/*
 *   Find irrev reaction
 */
  if (use.irrev_in == TRUE)
  {
    use.irrev_ptr = irrev_bsearch (use.n_irrev_user, &use.n_irrev);
    if (use.irrev_ptr == NULL)
    {
      sprintf (error_string, "Irreversible reaction %d not found.",
	       use.n_irrev_user);
      error_msg (error_string, STOP);
    }
  }
  else
  {
    use.irrev_ptr = NULL;
  }
/*
 *   Find exchange
 */
  if (use.exchange_in == TRUE)
  {
    use.exchange_ptr =
      exchange_bsearch (use.n_exchange_user, &use.n_exchange);
    if (use.exchange_ptr == NULL)
    {
      sprintf (error_string, "Exchange %d not found.", use.n_exchange_user);
      error_msg (error_string, STOP);
    }
  }
  else
  {
    use.exchange_ptr = NULL;
  }
/*
 *   Find kinetics
 */
  if (use.kinetics_in == TRUE)
  {
    use.kinetics_ptr =
      kinetics_bsearch (use.n_kinetics_user, &use.n_kinetics);
    if (use.kinetics_ptr == NULL)
    {
      sprintf (error_string, "Kinetics %d not found.", use.n_kinetics_user);
      error_msg (error_string, STOP);
    }
  }
  else
  {
    use.kinetics_ptr = NULL;
  }
/*
 *   Find surface
 */
  dl_type_x = NO_DL;
  if (use.surface_in == TRUE)
  {
    use.surface_ptr = surface_bsearch (use.n_surface_user, &use.n_surface);
    if (use.surface_ptr == NULL)
    {
      sprintf (error_string, "Surface %d not found.", use.n_surface_user);
      error_msg (error_string, STOP);
    }
  }
  else
  {
    use.surface_ptr = NULL;
  }
/*
 *   Find temperature
 */
  if (use.temperature_in == TRUE)
  {
    use.temperature_ptr =
      temperature_bsearch (use.n_temperature_user, &use.n_temperature);
    if (use.temperature_ptr == NULL)
    {
      sprintf (error_string, "Temperature %d not found.",
	       use.n_temperature_user);
      error_msg (error_string, STOP);
    }
  }
  else
  {
    use.temperature_ptr = NULL;
  }
/*
 *   Find gas
 */
  if (use.gas_phase_in == TRUE)
  {
    use.gas_phase_ptr =
      gas_phase_bsearch (use.n_gas_phase_user, &use.n_gas_phase);
    if (use.gas_phase_ptr == NULL)
    {
      sprintf (error_string, "Gas_phase %d not found.", use.n_gas_phase_user);
      error_msg (error_string, STOP);
    }
  }
  else
  {
    use.gas_phase_ptr = NULL;
  }
/*
 *   Find s_s_assemblage
 */
  if (use.s_s_assemblage_in == TRUE)
  {
    use.s_s_assemblage_ptr =
      s_s_assemblage_bsearch (use.n_s_s_assemblage_user,
			      &use.n_s_s_assemblage);
    if (use.s_s_assemblage_ptr == NULL)
    {
      sprintf (error_string, "s_s_assemblage %d not found.",
	       use.n_s_s_assemblage_user);
      error_msg (error_string, STOP);
    }
  }
  else
  {
    use.s_s_assemblage_ptr = NULL;
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
initial_solutions (int print)
/* ---------------------------------------------------------------------- */
{
/*
 *   Go through list of solutions, make initial solution calculations
 *   for any marked "new".
 */
  int i, converge, converge1;
  int n, last, n_user, print1;
  char token[2 * MAX_LENGTH];

  state = INITIAL_SOLUTION;
  set_use ();
  print1 = TRUE;
  dl_type_x = NO_DL;
  for (n = 0; n < count_solution; n++)
  {
    initial_solution_isotopes = FALSE;
    if (solution[n] != NULL && solution[n]->new_def == TRUE)
    {
      if (print1 == TRUE && print == TRUE)
      {
	dup_print ("Beginning of initial solution calculations.", TRUE);
	print1 = FALSE;
      }
      if (print == TRUE)
      {
	sprintf (token, "Initial solution %d.\t%.350s", solution[n]->n_user,
		 solution[n]->description);
	dup_print (token, FALSE);
      }
      use.solution_ptr = solution[n];
      prep ();
      k_temp (solution[n]->tc);
      set (TRUE);
      converge = model ();
      if (converge == ERROR && diagonal_scale == FALSE)
      {
	diagonal_scale = TRUE;
	set (TRUE);
	converge = model ();
	diagonal_scale = FALSE;
      }
      converge1 = check_residuals ();
      sum_species ();
      add_isotopes (solution[n]);
      punch_all ();
      print_all ();
      /* free_model_allocs(); */
      if (converge == ERROR || converge1 == ERROR)
      {
	error_msg ("Model failed to converge for initial solution.", STOP);
      }
      n_user = solution[n]->n_user;
      last = solution[n]->n_user_end;
      /* copy isotope data */
      if (solution[n]->count_isotopes > 0)
      {
	count_isotopes_x = solution[n]->count_isotopes;
	isotopes_x =
	  (struct isotope *) PHRQ_realloc (isotopes_x,
					   (size_t) count_isotopes_x *
					   sizeof (struct isotope));
	if (isotopes_x == NULL)
	  malloc_error ();
	memcpy (isotopes_x, solution[n]->isotopes,
		(size_t) count_isotopes_x * sizeof (struct isotope));
      }
      else
      {
	count_isotopes_x = 0;
      }
      xsolution_save (n_user);
      for (i = n_user + 1; i <= last; i++)
      {
	solution_duplicate (n_user, i);
      }
    }
  }
  initial_solution_isotopes = FALSE;
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
initial_exchangers (int print)
/* ---------------------------------------------------------------------- */
{
/*
 *   Go through list of exchangers, make initial calculations
 *   for any marked "new" that are defined to be in equilibrium with a 
 *   solution.
 */
  int i, converge, converge1;
  int n, last, n_user, print1;
  char token[2 * MAX_LENGTH];

  state = INITIAL_EXCHANGE;
  set_use ();
  print1 = TRUE;
  dl_type_x = NO_DL;
  for (n = 0; n < count_exchange; n++)
  {
    if (exchange[n].new_def != TRUE)
      continue;
    n_user = exchange[n].n_user;
    last = exchange[n].n_user_end;
    exchange[n].n_user_end = n_user;
    if (exchange[n].solution_equilibria == TRUE)
    {
      if (print1 == TRUE && print == TRUE)
      {
	dup_print ("Beginning of initial exchange"
		   "-composition calculations.", TRUE);
	print1 = FALSE;
      }
      if (print == TRUE)
      {
	sprintf (token, "Exchange %d.\t%.350s",
		 exchange[n].n_user, exchange[n].description);
	dup_print (token, FALSE);
      }
      use.exchange_ptr = &(exchange[n]);
      use.solution_ptr = solution_bsearch (exchange[n].n_solution, &i, TRUE);
      if (use.solution_ptr == NULL)
      {
	error_msg ("Solution not found for initial exchange calculation",
		   STOP);
      }
#ifdef PHREEQC_XML
      {
	/*
	   int n;
	   SAX_StartSystem();
	   SAX_AddSolution(use.solution_ptr);
	   SAX_EndSystem();
	   SAX_UnpackSolutions(SAX_GetXMLStr(), SAX_GetXMLLength());
	   SAX_cleanup();
	 */
      }
#endif
      prep ();
      k_temp (use.solution_ptr->tc);
      set (TRUE);
      converge = model ();
      converge1 = check_residuals ();
      sum_species ();
      species_list_sort ();
      print_exchange ();
      xexchange_save (n_user);
      punch_all ();
      /* free_model_allocs(); */
      if (converge == ERROR || converge1 == ERROR)
      {
	error_msg
	  ("Model failed to converge for initial exchange calculation.",
	   STOP);
      }
    }
    for (i = n_user + 1; i <= last; i++)
    {
      exchange_duplicate (n_user, i);
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
initial_gas_phases (int print)
/* ---------------------------------------------------------------------- */
{
/*
 *   Go through list of gas_phases, make initial calculations
 *   for any marked "new" that are defined to be in equilibrium with a 
 *   solution.
 */
  int i, converge, converge1;
  int n, last, n_user, print1;
  char token[2 * MAX_LENGTH];
  struct gas_comp *gas_comp_ptr;
  struct phase *phase_ptr;
  struct rxn_token *rxn_ptr;
  LDBLE lp;

  state = INITIAL_GAS_PHASE;
  set_use ();
  print1 = TRUE;
  dl_type_x = NO_DL;
  for (n = 0; n < count_gas_phase; n++)
  {
    if (gas_phase[n].new_def != TRUE)
      continue;
    n_user = gas_phase[n].n_user;
    last = gas_phase[n].n_user_end;
    gas_phase[n].n_user_end = n_user;
    if (gas_phase[n].solution_equilibria == TRUE)
    {
      if (print1 == TRUE && print == TRUE)
      {
	dup_print ("Beginning of initial gas_phase"
		   "-composition calculations.", TRUE);
	print1 = FALSE;
      }
      if (print == TRUE)
      {
	sprintf (token, "Gas_Phase %d.\t%.350s",
		 gas_phase[n].n_user, gas_phase[n].description);
	dup_print (token, FALSE);
      }
      use.solution_ptr = solution_bsearch (gas_phase[n].n_solution, &i, TRUE);
      prep ();
      k_temp (use.solution_ptr->tc);
      set (TRUE);
      converge = model ();
      converge1 = check_residuals ();
      if (converge == ERROR || converge1 == ERROR)
      {
	/* free_model_allocs(); */
	error_msg
	  ("Model failed to converge for initial gas phase calculation.",
	   STOP);
      }
      use.gas_phase_ptr = &(gas_phase[n]);
      use.gas_phase_ptr->total_p = 0;
      use.gas_phase_ptr->total_moles = 0;
      for (i = 0; i < use.gas_phase_ptr->count_comps; i++)
      {
	gas_comp_ptr = &(use.gas_phase_ptr->comps[i]);
	phase_ptr = gas_comp_ptr->phase;
	if (phase_ptr->in == TRUE)
	{
	  lp = -phase_ptr->lk;
	  for (rxn_ptr = phase_ptr->rxn_x->token + 1; rxn_ptr->s != NULL;
	       rxn_ptr++)
	  {
	    lp += rxn_ptr->s->la * rxn_ptr->coef;
	  }
	  gas_comp_ptr->phase->p_soln_x = exp (lp * LOG_10);
	  use.gas_phase_ptr->total_p += gas_comp_ptr->phase->p_soln_x;
	  gas_comp_ptr->phase->moles_x = gas_comp_ptr->phase->p_soln_x *
	    use.gas_phase_ptr->volume / (R_LITER_ATM * tk_x);
	  gas_comp_ptr->moles = gas_comp_ptr->phase->moles_x;
	  use.gas_phase_ptr->total_moles += gas_comp_ptr->moles;
	}
	else
	{
	  gas_comp_ptr->phase->moles_x = 0;
	}
      }
      print_gas_phase ();
      xgas_save (n_user);
      punch_all ();
      /* free_model_allocs(); */
    }
    for (i = n_user + 1; i <= last; i++)
    {
      gas_phase_duplicate (n_user, i);
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
initial_surfaces (int print)
/* ---------------------------------------------------------------------- */
{
/*
 *   Go through list of surfaces, make initial calculations
 *   for any marked "new" that are defined to be in equilibrium with a 
 *   solution.
 */
  int i;
  int n, last, n_user, print1;
  char token[2 * MAX_LENGTH];

  state = INITIAL_SURFACE;
  set_use ();
  print1 = TRUE;
  for (n = 0; n < count_surface; n++)
  {
    if (surface[n].new_def != TRUE)
      continue;
    n_user = surface[n].n_user;
    last = surface[n].n_user_end;
    surface[n].n_user_end = n_user;
    if (surface[n].solution_equilibria == TRUE)
    {
      if (print1 == TRUE && print == TRUE)
      {
	dup_print ("Beginning of initial surface-composition calculations.",
		   TRUE);
	print1 = FALSE;
      }
      if (print == TRUE)
      {
	sprintf (token, "Surface %d.\t%.350s",
		 surface[n].n_user, surface[n].description);
	dup_print (token, FALSE);
      }
      use.surface_ptr = &(surface[n]);
      dl_type_x = use.surface_ptr->dl_type;
      use.solution_ptr = solution_bsearch (surface[n].n_solution, &i, TRUE);
      if (use.solution_ptr == NULL)
      {
	error_msg ("Solution not found for initial surface calculation",
		   STOP);
      }
#ifdef SKIP
      if (diffuse_layer_x == TRUE)
      {
	converge = surface_model ();
      }
      else
      {
	prep ();
	k_temp (use.solution_ptr->tc);
	set (TRUE);
	converge = model ();
      }
      converge1 = check_residuals ();
      sum_species ();
#endif
      set_and_run_wrapper (-1, FALSE, FALSE, -1, 0.0);
      species_list_sort ();
      print_surface ();
      /*print_all(); */
      punch_all ();
#ifdef SKIP
      if (converge == ERROR || converge1 == ERROR)
      {
	/* free_model_allocs(); */
	error_msg
	  ("Model failed to converge for initial surface calculation.", STOP);
      }
#endif
      xsurface_save (n_user);
      /* free_model_allocs(); */
    }
    for (i = n_user + 1; i <= last; i++)
    {
      surface_duplicate (n_user, i);
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
reactions (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Make all reaction calculation which could include:
 *      equilibrium with a pure-phase assemblage,
 *      equilibrium with an exchanger,
 *      equilibrium with an surface,
 *      equilibrium with a gas phase,
 *      equilibrium with a solid solution assemblage,
 *      kinetics,
 *      change of temperature,
 *      mixture,
 *      or irreversible reaction.
 */
  int count_steps, use_mix, m;
  char token[2 * MAX_LENGTH];
  struct save save_data;
  LDBLE kin_time;
  struct kinetics *kinetics_ptr;

  state = REACTION;
  /* last_model.force_prep = TRUE; */
  if (set_use () == FALSE)
    return (OK);
/*
 *   Find maximum number of steps
 */
  dup_print ("Beginning of batch-reaction calculations.", TRUE);
  count_steps = 1;
  if (use.irrev_in == TRUE && use.irrev_ptr != NULL)
  {
    if (abs (use.irrev_ptr->count_steps) > count_steps)
      count_steps = abs (use.irrev_ptr->count_steps);
  }
  if (use.kinetics_in == TRUE && use.kinetics_ptr != NULL)
  {
    if (abs (use.kinetics_ptr->count_steps) > count_steps)
      count_steps = abs (use.kinetics_ptr->count_steps);
  }
  if (use.temperature_in == TRUE && use.temperature_ptr != NULL)
  {
    if (abs (use.temperature_ptr->count_t) > count_steps)
      count_steps = abs (use.temperature_ptr->count_t);
  }
  count_total_steps = count_steps;
/*
 *  save data for saving solutions
 */
  memcpy (&save_data, &save, sizeof (struct save));
  /* 
   *Copy everything to -2
   */
  copy_use (-2);
  rate_sim_time_start = 0;
  rate_sim_time = 0;
  for (reaction_step = 1; reaction_step <= count_steps; reaction_step++)
  {
    sprintf (token, "Reaction step %d.", reaction_step);
    if (reaction_step > 1 && incremental_reactions == FALSE)
    {
      copy_use (-2);
    }
    set_initial_moles (-2);
    dup_print (token, FALSE);
/*
 *  Determine time step for kinetics
 */
    kin_time = 0.0;
    if (use.kinetics_in == TRUE)
    {
      kinetics_ptr = kinetics_bsearch (-2, &m);
      if (incremental_reactions == FALSE)
      {
	if (kinetics_ptr->count_steps > 0)
	{
	  if (reaction_step > kinetics_ptr->count_steps)
	  {
	    kin_time = kinetics_ptr->steps[kinetics_ptr->count_steps - 1];
	  }
	  else
	  {
	    kin_time = kinetics_ptr->steps[reaction_step - 1];
	  }
	}
	else if (kinetics_ptr->count_steps < 0)
	{
	  if (reaction_step > -kinetics_ptr->count_steps)
	  {
	    kin_time = kinetics_ptr->steps[0];
	  }
	  else
	  {
	    kin_time =
	      reaction_step * kinetics_ptr->steps[0] /
	      ((LDBLE) (-kinetics_ptr->count_steps));
	  }
	}
      }
      else
      {
	/* incremental reactions */
	if (kinetics_ptr->count_steps > 0)
	{
	  if (reaction_step > kinetics_ptr->count_steps)
	  {
	    kin_time = kinetics_ptr->steps[kinetics_ptr->count_steps - 1];
	  }
	  else
	  {
	    kin_time = kinetics_ptr->steps[reaction_step - 1];
	  }
	}
	else if (kinetics_ptr->count_steps < 0)
	{
	  if (reaction_step > -kinetics_ptr->count_steps)
	  {
	    kin_time = 0;
	  }
	  else
	  {
	    kin_time =
	      kinetics_ptr->steps[0] / ((LDBLE) (-kinetics_ptr->count_steps));
	  }
	}
      }
    }
    if (incremental_reactions == FALSE ||
	(incremental_reactions == TRUE && reaction_step == 1))
    {
      use_mix = TRUE;
    }
    else
    {
      use_mix = FALSE;
    }
/*
 *   Run reaction step
 */
    run_reactions (-2, kin_time, use_mix, 1.0);
    if (incremental_reactions == TRUE)
    {
      rate_sim_time_start += kin_time;
      rate_sim_time = rate_sim_time_start;
    }
    else
    {
      rate_sim_time = kin_time;
    }
    if (state != ADVECTION)
    {
      punch_all ();
      print_all ();
    }
    /* saves back into -2 */
    if (reaction_step < count_steps)
    {
      saver ();
    }
  }
/*
 *   save end of reaction
 */
  memcpy (&save, &save_data, sizeof (struct save));
  if (use.kinetics_in == TRUE)
  {
    kinetics_duplicate (-2, use.n_kinetics_user);
  }
  saver ();

  /* free_model_allocs(); */
  /* last_model.force_prep = TRUE; */
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
saver (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Save results of calcuations (data in variables with _x,
 *   in unknown structure x, in master, or s) into structure 
 *   arrays.  Structure "save" has info on whether to save
 *   data for each entity (solution, ex, surf, pp, gas, or s_s).
 *   Initial calculation may be saved into multiple "n_user"
 *   slots.
 */
  int i, n;
  char token[MAX_LENGTH];

  if (save.solution == TRUE)
  {
    sprintf (token, "Solution after simulation %d.", simulation);
    description_x = (char *) free_check_null (description_x);
    description_x = string_duplicate (token);
    n = save.n_solution_user;
    xsolution_save (n);
    for (i = save.n_solution_user + 1; i <= save.n_solution_user_end; i++)
    {
      solution_duplicate (n, i);
    }
  }
  if (save.pp_assemblage == TRUE)
  {
    n = save.n_pp_assemblage_user;
    xpp_assemblage_save (n);
    for (i = save.n_pp_assemblage_user + 1;
	 i <= save.n_pp_assemblage_user_end; i++)
    {
      pp_assemblage_duplicate (n, i);
    }
  }
  if (save.exchange == TRUE)
  {
    n = save.n_exchange_user;
    xexchange_save (n);
    for (i = save.n_exchange_user + 1; i <= save.n_exchange_user_end; i++)
    {
      exchange_duplicate (n, i);
    }
  }
  if (save.surface == TRUE)
  {
    n = save.n_surface_user;
    xsurface_save (n);
    for (i = save.n_surface_user + 1; i <= save.n_surface_user_end; i++)
    {
      surface_duplicate (n, i);
    }
  }
  if (save.gas_phase == TRUE)
  {
    n = save.n_gas_phase_user;
    xgas_save (n);
    for (i = save.n_gas_phase_user + 1; i <= save.n_gas_phase_user_end; i++)
    {
      gas_phase_duplicate (n, i);
    }
  }
  if (save.s_s_assemblage == TRUE)
  {
    n = save.n_s_s_assemblage_user;
    xs_s_assemblage_save (n);
    for (i = save.n_s_s_assemblage_user + 1;
	 i <= save.n_s_s_assemblage_user_end; i++)
    {
      s_s_assemblage_duplicate (n, i);
    }
  }
  if (save.kinetics == TRUE && use.kinetics_in == TRUE
      && use.kinetics_ptr != NULL)
  {
    n = use.kinetics_ptr->n_user;
    for (i = save.n_kinetics_user; i <= save.n_kinetics_user_end; i++)
    {
      kinetics_duplicate (n, i);
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
static int
xexchange_save (int n_user)
/* ---------------------------------------------------------------------- */
{
/*
 *   Save exchanger assemblage into structure exchange with user 
 *   number n_user.
 */
  int i, j, n;
  char token[MAX_LENGTH];
  int count_comps;
  struct exchange temp_exchange, *exchange_ptr;
  LDBLE charge;
  if (use.exchange_ptr == NULL)
    return (OK);
/*
 *   Store data for structure exchange
 */
  memcpy (&temp_exchange, use.exchange_ptr, sizeof (struct exchange));
  temp_exchange.n_user = n_user;
  temp_exchange.n_user_end = n_user;
  temp_exchange.new_def = FALSE;
  sprintf (token, "Exchange assemblage after simulation %d.", simulation);
  temp_exchange.description = string_duplicate (token);
  temp_exchange.solution_equilibria = FALSE;
  temp_exchange.n_solution = -2;
  temp_exchange.count_comps = 0;
/*
 *   Count exchange components and allocate space
 */
  count_comps = 0;
  for (i = 0; i < count_unknowns; i++)
  {
    if (x[i]->type == EXCH)
      count_comps++;
  }
  temp_exchange.comps =
    (struct exch_comp *) PHRQ_malloc ((size_t) (count_comps) *
				      sizeof (struct exch_comp));
  if (temp_exchange.comps == NULL)
    malloc_error ();
  temp_exchange.count_comps = count_comps;
/*
 *   Write exch_comp structure for each exchange component
 */
  count_comps = 0;
  for (i = 0; i < count_unknowns; i++)
  {
    if (x[i]->type == EXCH)
    {
      memcpy (&temp_exchange.comps[count_comps], x[i]->exch_comp,
	      sizeof (struct exch_comp));
      temp_exchange.comps[count_comps].master = x[i]->master[0];
      temp_exchange.comps[count_comps].la = x[i]->master[0]->s->la;
      temp_exchange.comps[count_comps].formula_totals =
	elt_list_dup (x[i]->exch_comp->formula_totals);
      temp_exchange.comps[count_comps].moles = 0.;
/*
 *   Save element concentrations on exchanger
 */
      count_elts = 0;
      paren_count = 0;
      charge = 0.0;
      for (j = 0; j < count_species_list; j++)
      {
	if (species_list[j].master_s == x[i]->master[0]->s)
	{
	  add_elt_list (species_list[j].s->next_elt,
			species_list[j].s->moles);
	  charge += species_list[j].s->moles * species_list[j].s->z;
	}
      }
/*
 *   Keep exchanger related to phase even if none currently in solution
 */
      if (x[i]->exch_comp->phase_name != NULL && count_elts == 0)
      {
	add_elt_list (x[i]->master[0]->s->next_elt, 1e-20);
      }
/*
 *   Store list
 */
      temp_exchange.comps[count_comps].charge_balance = charge;
      temp_exchange.comps[count_comps].totals = elt_list_save ();
/* debug
                        output_msg(OUTPUT_MESSAGE, "Exchange charge_balance: %e\n", charge);
 */
      /* update unknown pointer */
      x[i]->exch_comp = &(temp_exchange.comps[count_comps]);
      count_comps++;
    }
  }
/*
 *   Finish up
 */
  exchange_ptr = exchange_bsearch (n_user, &n);
  if (exchange_ptr == NULL)
  {
    n = count_exchange++;
    space ((void **) ((void *) &exchange), count_exchange, &max_exchange,
	   sizeof (struct exchange));
  }
  else
  {
    exchange_free (&exchange[n]);
  }
  memcpy (&exchange[n], &temp_exchange, sizeof (struct exchange));
  /* sort only if necessary */
  if (n == count_exchange - 1 && count_exchange > 1)
  {
    if (exchange[n].n_user < exchange[n - 1].n_user)
    {
      qsort (exchange,
	     (size_t) count_exchange,
	     (size_t) sizeof (struct exchange), exchange_compare);
    }
  }
  use.exchange_ptr = NULL;
  return (OK);
}

/* ---------------------------------------------------------------------- */
static int
xgas_save (int n_user)
/* ---------------------------------------------------------------------- */
{
/*
 *   Save gas composition into structure gas_phase with user 
 *   number n_user.
 */
  int count_comps, n, i;
  struct gas_phase temp_gas_phase, *gas_phase_ptr;
  char token[MAX_LENGTH];

  if (use.gas_phase_ptr == NULL)
    return (OK);
/*
 *   Count gases
 */
  count_comps = use.gas_phase_ptr->count_comps;
/*
 *   Malloc space and copy
 */
  temp_gas_phase.comps =
    (struct gas_comp *) PHRQ_malloc ((size_t) count_comps *
				     sizeof (struct gas_comp));
  if (temp_gas_phase.comps == NULL)
    malloc_error ();
  memcpy ((void *) temp_gas_phase.comps, (void *) use.gas_phase_ptr->comps,
	  (size_t) count_comps * sizeof (struct gas_comp));
/*
 *   Store in gas_phase
 */

  temp_gas_phase.n_user = n_user;
  temp_gas_phase.n_user_end = n_user;
  sprintf (token, "Gas phase after simulation %d.", simulation);
  temp_gas_phase.description = string_duplicate (token);
  temp_gas_phase.new_def = FALSE;
  temp_gas_phase.solution_equilibria = FALSE;
  temp_gas_phase.n_solution = -99;
  temp_gas_phase.type = use.gas_phase_ptr->type;
  temp_gas_phase.total_p = use.gas_phase_ptr->total_p;
  temp_gas_phase.total_moles = use.gas_phase_ptr->total_moles;
  temp_gas_phase.volume = use.gas_phase_ptr->volume;
  temp_gas_phase.temperature = use.gas_phase_ptr->temperature;
  temp_gas_phase.count_comps = count_comps;
/*
 *   Update amounts
 */
  for (i = 0; i < count_comps; i++)
  {
    temp_gas_phase.comps[i].moles =
      use.gas_phase_ptr->comps[i].phase->moles_x;
  }
/*
 *   Finish up
 */
  gas_phase_ptr = gas_phase_bsearch (n_user, &n);
  if (gas_phase_ptr == NULL)
  {
    n = count_gas_phase++;
    space ((void **) ((void *) &gas_phase), count_gas_phase, &max_gas_phase,
	   sizeof (struct gas_phase));
  }
  else
  {
    gas_phase_free (&gas_phase[n]);
  }
  memcpy (&gas_phase[n], &temp_gas_phase, sizeof (struct gas_phase));

  /* update unknown pointer */
  if (gas_unknown != NULL)
  {
    gas_unknown->gas_phase = &(gas_phase[n]);
  }
  /* sort only if necessary */
  if (n == count_gas_phase - 1 && count_gas_phase > 1)
  {
    if (gas_phase[n].n_user < gas_phase[n - 1].n_user)
    {
      qsort (gas_phase,
	     (size_t) count_gas_phase,
	     (size_t) sizeof (struct gas_phase), gas_phase_compare);
    }
  }
  use.gas_phase_ptr = NULL;
  return (OK);
}

/* ---------------------------------------------------------------------- */
static int
xs_s_assemblage_save (int n_user)
/* ---------------------------------------------------------------------- */
{
/*
 *   Save s_s_assemblage composition into structure s_s_assemblage with user 
 *   number n_user.
 */
  int i, j, n;
  int count_comps, count_s_s;
  struct s_s_assemblage temp_s_s_assemblage, *s_s_assemblage_ptr;
  char token[MAX_LENGTH];

  if (use.s_s_assemblage_ptr == NULL)
    return (OK);
/*
 *   Set s_s_assemblage
 */
  temp_s_s_assemblage.n_user = n_user;
  temp_s_s_assemblage.n_user_end = n_user;
  sprintf (token, "Solid solution assemblage after simulation %d.",
	   simulation);
  temp_s_s_assemblage.description = string_duplicate (token);
  temp_s_s_assemblage.new_def = FALSE;
#ifdef SKIP
  temp_s_s_assemblage.type = use.s_s_assemblage_ptr->type;
  temp_s_s_assemblage.solution_equilibria = FALSE;
  temp_s_s_assemblage.n_solution = -99;
#endif
  count_s_s = use.s_s_assemblage_ptr->count_s_s;
  temp_s_s_assemblage.count_s_s = count_s_s;
/*
 *   Malloc space for solid solutions
 */
  /* s_s_assemblage->s_s */
  temp_s_s_assemblage.s_s =
    (struct s_s *) PHRQ_malloc ((size_t) count_s_s * sizeof (struct s_s));
  if (temp_s_s_assemblage.s_s == NULL)
    malloc_error ();
  for (i = 0; i < count_s_s; i++)
  {
    memcpy (&(temp_s_s_assemblage.s_s[i]), &(use.s_s_assemblage_ptr->s_s[i]),
	    sizeof (struct s_s));
    /* 
     * Malloc space for solid soution components
     */
    count_comps = use.s_s_assemblage_ptr->s_s[i].count_comps;
    temp_s_s_assemblage.s_s[i].comps =
      (struct s_s_comp *) PHRQ_malloc ((size_t) count_comps *
				       sizeof (struct s_s_comp));
    if (temp_s_s_assemblage.s_s[i].comps == NULL)
      malloc_error ();
    memcpy ((void *) temp_s_s_assemblage.s_s[i].comps,
	    (void *) use.s_s_assemblage_ptr->s_s[i].comps,
	    (size_t) count_comps * sizeof (struct s_s_comp));

    /* set initial moles for quick setup */
    for (j = 0; j < count_comps; j++)
    {
      temp_s_s_assemblage.s_s[i].comps[j].initial_moles =
	temp_s_s_assemblage.s_s[i].comps[j].moles;
    }
  }
/*
 *   Finish up
 */
  s_s_assemblage_ptr = s_s_assemblage_bsearch (n_user, &n);
  if (s_s_assemblage_ptr == NULL)
  {
    space ((void **) ((void *) &s_s_assemblage), count_s_s_assemblage,
	   &max_s_s_assemblage, sizeof (struct s_s_assemblage));
    n = count_s_s_assemblage++;
  }
  else
  {
    s_s_assemblage_free (&s_s_assemblage[n]);
  }
  memcpy (&s_s_assemblage[n], &temp_s_s_assemblage,
	  sizeof (struct s_s_assemblage));
  /* sort only if necessary */
  if (n == count_s_s_assemblage - 1 && count_s_s_assemblage > 1)
  {
    if (s_s_assemblage[n].n_user < s_s_assemblage[n - 1].n_user)
    {
      qsort (s_s_assemblage,
	     (size_t) count_s_s_assemblage,
	     (size_t) sizeof (struct s_s_assemblage), s_s_assemblage_compare);
    }
  }
  use.s_s_assemblage_ptr = NULL;
  return (OK);
}

/* ---------------------------------------------------------------------- */
static int
xpp_assemblage_save (int n_user)
/* ---------------------------------------------------------------------- */
{
/*
 *   Save pure_phase assemblage into structure pp_assemblage with user 
 *   number n_user.
 */
  int count_comps, n, j, i;
  struct pp_assemblage temp_pp_assemblage, *pp_assemblage_ptr;
  char token[MAX_LENGTH];

  if (use.pp_assemblage_ptr == NULL)
    return (OK);
/*
 *   Count pure phases
 */
  count_comps = use.pp_assemblage_ptr->count_comps;

  temp_pp_assemblage.n_user = n_user;
  temp_pp_assemblage.n_user_end = n_user;
  sprintf (token, "Pure-phase assemblage after simulation %d.", simulation);
  temp_pp_assemblage.description = string_duplicate (token);
  temp_pp_assemblage.new_def = FALSE;
  temp_pp_assemblage.count_comps = count_comps;
  temp_pp_assemblage.next_elt =
    elt_list_dup (use.pp_assemblage_ptr->next_elt);
/*
 *   Malloc space and copy pure phase data
 */
  temp_pp_assemblage.pure_phases =
    (struct pure_phase *) PHRQ_malloc ((size_t) count_comps *
				       sizeof (struct pure_phase));
  if (temp_pp_assemblage.pure_phases == NULL)
    malloc_error ();
  memcpy ((void *) temp_pp_assemblage.pure_phases,
	  (void *) use.pp_assemblage_ptr->pure_phases,
	  (size_t) count_comps * sizeof (struct pure_phase));
/*
 *   Update amounts
 */
  i = 0;
  for (j = 0; j < count_unknowns; j++)
  {
    if (x[j]->type != PP)
      continue;
    temp_pp_assemblage.pure_phases[i].moles = x[j]->moles;
    temp_pp_assemblage.pure_phases[i].delta = 0.0;

    /* update unknown ptr, old may be freed later */
    x[j]->pure_phase = &(temp_pp_assemblage.pure_phases[i]);
    i++;
  }
/*
 *   Finish up
 */
  pp_assemblage_ptr = pp_assemblage_bsearch (n_user, &n);
  if (pp_assemblage_ptr == NULL)
  {
    space ((void **) ((void *) &pp_assemblage), count_pp_assemblage,
	   &max_pp_assemblage, sizeof (struct pp_assemblage));
    n = count_pp_assemblage++;
  }
  else
  {
    pp_assemblage_free (&pp_assemblage[n]);
  }
  memcpy (&pp_assemblage[n], &temp_pp_assemblage,
	  sizeof (struct pp_assemblage));
  /* sort only if necessary */
  if (n == count_pp_assemblage - 1 && count_pp_assemblage > 1)
  {
    if (pp_assemblage[n].n_user < pp_assemblage[n - 1].n_user)
    {
      qsort (pp_assemblage,
	     (size_t) count_pp_assemblage,
	     (size_t) sizeof (struct pp_assemblage), pp_assemblage_compare);
    }
  }
  use.pp_assemblage_ptr = NULL;
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
xsolution_save (int n_user)
/* ---------------------------------------------------------------------- */
{
/*
 *   Save solution composition into structure solution with user number
 *   n_user.
 *
 *   input:  n_user is user solution number of target
 */
  int i, j, n;
  int count_mass_balance, count_master_activity;
  int max_mass_balance, max_master_activity;
  struct solution *solution_ptr;
  struct master *master_i_ptr, *master_ptr;
/*
 *   Malloc space for solution data
 */
  solution_ptr = solution_alloc ();

  max_mass_balance = MAX_MASS_BALANCE;
  max_master_activity = MAX_MASS_BALANCE;

  solution_ptr->n_user = n_user;
  solution_ptr->n_user_end = n_user;
  solution_ptr->new_def = FALSE;
  solution_ptr->description = string_duplicate (description_x);
  solution_ptr->tc = tc_x;
  solution_ptr->ph = ph_x;
  solution_ptr->solution_pe = solution_pe_x;
  solution_ptr->mu = mu_x;
  solution_ptr->ah2o = ah2o_x;
  solution_ptr->density = density_x;
  solution_ptr->total_h = total_h_x;
  solution_ptr->total_o = total_o_x;
  solution_ptr->cb = cb_x;	/* cb_x does not include surface charge sfter sum_species */
  /* does include surface charge after step */

  solution_ptr->mass_water = mass_water_aq_x;
  solution_ptr->total_alkalinity = total_alkalinity;
  /*solution_ptr->total_co2 = total_co2; */
  solution_ptr->units = moles_per_kilogram_string;
/*
 *   Copy pe data
 */

  pe_data_free (solution_ptr->pe);
  /*solution_ptr->pe = pe_data_dup(pe_x); */
  solution_ptr->pe = pe_data_alloc ();
  solution_ptr->default_pe = 0;
  /*
   * Add in minor isotopes if initial solution calculation
   */
  if (initial_solution_isotopes == TRUE)
  {
    for (i = 0; i < count_master_isotope; i++)
    {
      if (master_isotope[i]->moles > 0)
      {
	master_i_ptr = master_bsearch (master_isotope[i]->name);
	master_ptr = master_isotope[i]->elt->master;
	if (master_isotope[i]->minor_isotope == TRUE)
	{
	  master_i_ptr->total = master_isotope[i]->moles;
	  if (master_ptr->total > 0)
	  {
	    master_i_ptr->s->la =
	      master_ptr->s->la +
	      log10 (master_i_ptr->total / master_ptr->total);
	  }
	  else
	  {
	    master_i_ptr->s->la = master_ptr->s->la;
	  }
	}
	else if (master_isotope[i]->minor_isotope == FALSE
		 && master_ptr->s != s_hplus && master_ptr->s != s_h2o)
	{
	  if (master_ptr->s->secondary != NULL)
	  {
	    master_ptr->s->secondary->total = master_isotope[i]->moles;
	  }
	  else
	  {
	    master_ptr->s->primary->total = master_isotope[i]->moles;
	  }
	}
      }
    }
  }
/*
 *   Copy totals data
 */
  count_mass_balance = 0;
  count_master_activity = 0;
  for (i = 0; i < count_master; i++)
  {
    if (master[i]->in == FALSE) continue;
    if (master[i]->s->type == EX ||
	master[i]->s->type == SURF || master[i]->s->type == SURF_PSI)
      continue;
    if (master[i]->s == s_hplus) continue;
    if (master[i]->s == s_h2o) continue;
/*
 *   Save list of log activities
 */
    if (master[i]->in != FALSE)
    {
      count_master_activity++;
    }
  }
  solution_ptr->master_activity =
    (struct master_activity *) PHRQ_realloc (solution_ptr->master_activity,
					     (size_t) (count_master_activity +
						       1) *
					     sizeof (struct master_activity));
  solution_ptr->count_master_activity = count_master_activity;
  count_master_activity = 0;
  for (i = 0; i < count_master; i++)
  {
    if (master[i]->s->type == EX ||
	master[i]->s->type == SURF || master[i]->s->type == SURF_PSI)
      continue;
    if (master[i]->s == s_hplus) continue;
    if (master[i]->s == s_h2o) continue;
/*
 *   Save list of log activities
 */
    if (master[i]->in != FALSE)
    {
      solution_ptr->master_activity[count_master_activity].description =
	master[i]->elt->name;
      solution_ptr->master_activity[count_master_activity++].la =
	master[i]->s->la;
      /*
         if (count_master_activity + 2 >= max_master_activity) {
         space ((void *) &(solution_ptr->master_activity), count_master_activity + 2,
         &max_master_activity, sizeof (struct master_activity));
         }
       */
    }
    if (master[i]->total <= MIN_TOTAL)
    {
      master[i]->total = 0.0;
      master[i]->total_primary = 0.0;
      continue;
    }
/*              if (master[i]->total <= 0.0) continue;  */
/*
 *   Save list of concentrations
 */
    solution_ptr->totals[count_mass_balance].description =
      master[i]->elt->name;
    solution_ptr->totals[count_mass_balance].input_conc = master[i]->total;
    solution_ptr->totals[count_mass_balance].moles = master[i]->total;
    solution_ptr->totals[count_mass_balance].units = solution_ptr->units;
    solution_ptr->totals[count_mass_balance].equation_name = NULL;
    solution_ptr->totals[count_mass_balance].n_pe = 0;
    solution_ptr->totals[count_mass_balance].phase = NULL;
    solution_ptr->totals[count_mass_balance].phase_si = 0.0;
    solution_ptr->totals[count_mass_balance].as = NULL;
    solution_ptr->totals[count_mass_balance].gfw = 0.0;
    count_mass_balance++;
/*
 *   Make space
 */
    if (count_mass_balance + 2 >= max_mass_balance)
    {
      space ((void **) ((void *) &(solution_ptr->totals)), count_mass_balance + 2,
	     &max_mass_balance, sizeof (struct conc));
    }
  }
  if (pitzer_model == TRUE)
  {
    i = 0;
    for (j = 0; j < count_s; j++)
    {
      if (s[j]->lg != 0.0)
	i++;
    }
    solution_ptr->species_gamma =
      (struct master_activity *) PHRQ_realloc (solution_ptr->species_gamma,
					       (size_t) (i *
							 sizeof (struct
								 master_activity)));
    i = 0;
    for (j = 0; j < count_s; j++)
    {
      if (s[j]->lg != 0.0)
      {
	solution_ptr->species_gamma[i].la = s[j]->lg;
	solution_ptr->species_gamma[i].description = s[j]->name;
	i++;
      }
    }
    solution_ptr->count_species_gamma = i;
  }
  else
  {
    solution_ptr->species_gamma = NULL;
    solution_ptr->count_species_gamma = 0;
  }
/*
 *   Mark end of totals
 */
  solution_ptr->totals[count_mass_balance].description = NULL;
  count_mass_balance++;
  solution_ptr->master_activity[count_master_activity].description = NULL;
  count_master_activity++;
  solution_ptr->totals =
    (struct conc *) PHRQ_realloc (solution_ptr->totals,
				  (size_t) count_mass_balance *
				  sizeof (struct conc));
  if (solution_ptr->totals == NULL)
    malloc_error ();
  solution_ptr->master_activity =
    (struct master_activity *) PHRQ_realloc (solution_ptr->master_activity,
					     (size_t) count_master_activity *
					     sizeof (struct master_activity));
  if (solution_ptr->master_activity == NULL)
    malloc_error ();
  solution_ptr->count_master_activity = count_master_activity;
/*
 *   Save isotope data
 */
  if (count_isotopes_x > 0)
  {
    solution_ptr->count_isotopes = count_isotopes_x;
    solution_ptr->isotopes =
      (struct isotope *) PHRQ_realloc (solution_ptr->isotopes,
				       (size_t) count_isotopes_x *
				       sizeof (struct isotope));
    if (solution_ptr->isotopes == NULL)
      malloc_error ();
    memcpy (solution_ptr->isotopes, isotopes_x,
	    (size_t) count_isotopes_x * sizeof (struct isotope));
    for (i = 0; i < count_isotopes_x; i++)
    {
      solution_ptr->isotopes[i].total =
	solution_ptr->isotopes[i].master->total;
      if (solution_ptr->isotopes[i].master == s_hplus->secondary)
      {
	solution_ptr->isotopes[i].total = 2 * mass_water_aq_x / gfw_water;
      }
      if (solution_ptr->isotopes[i].master == s_h2o->secondary)
      {
	solution_ptr->isotopes[i].total = mass_water_aq_x / gfw_water;
      }
    }
  }
  else
  {
    solution_ptr->count_isotopes = 0;
    solution_ptr->isotopes =
      (struct isotope *) free_check_null (solution_ptr->isotopes);
    solution_ptr->isotopes = NULL;
  }
/*
 *   Save solution in solution
 */
  if (solution_bsearch (n_user, &n, FALSE) != NULL)
  {
    solution_free (solution[n]);
  }
  else
  {
    n = count_solution++;
    if (count_solution >= max_solution)
    {
      space ((void **) ((void *) &(solution)), count_solution, &max_solution,
	     sizeof (struct solution *));
    }
  }
  solution[n] = solution_ptr;
  /* sort only if necessary */
  if (count_solution > 1
      && (solution[count_solution - 1]->n_user <
	  solution[count_solution - 2]->n_user))
  {
    solution_sort ();
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
static int
xsurface_save (int n_user)
/* ---------------------------------------------------------------------- */
{
/*
 *   Save surface data into structure surface with user 
 *   number n_user.
 */
  int i, j, n, last_charge;
  int count_comps, count_charge;
  char token[MAX_LENGTH];
  struct surface temp_surface, *surface_ptr;
  LDBLE charge;
  if (use.surface_ptr == NULL)
    return (OK);
/*
 *   Store data for structure surface
 */
  memcpy (&temp_surface, use.surface_ptr, sizeof (struct surface));
  temp_surface.n_user = n_user;
  temp_surface.n_user_end = n_user;
  temp_surface.new_def = FALSE;
  temp_surface.dl_type = dl_type_x;
  sprintf (token, "Surface assemblage after simulation %d.", simulation);
  temp_surface.description = string_duplicate (token);
  temp_surface.solution_equilibria = FALSE;
  temp_surface.n_solution = -10;
/*
 *   Allocate space to pointer comps
 */

  count_comps = use.surface_ptr->count_comps;
  count_charge = use.surface_ptr->count_charge;
  temp_surface.count_comps = count_comps;
  temp_surface.comps =
    (struct surface_comp *) PHRQ_malloc ((size_t) (count_comps) *
					 sizeof (struct surface_comp));
  if (temp_surface.comps == NULL)
    malloc_error ();
  /*if (temp_surface.edl == FALSE) { */
  if (temp_surface.type == NO_EDL)
  {
    temp_surface.charge = NULL;
    temp_surface.count_charge = 0;
  }
  else
  {
    temp_surface.count_charge = count_charge;
    temp_surface.charge =
      (struct surface_charge *) PHRQ_malloc ((size_t) (count_charge) *
					     sizeof (struct surface_charge));
    if (temp_surface.charge == NULL)
      malloc_error ();
  }
/*
 *   Write surface_comp structure for each surf component into comps_ptr
 */
  count_comps = 0;
  count_charge = 0;
  /*
   *  Initial entry of surface sites is random
   *  Charge balance numbering follows the initial entry
   *  Surface sites are then sorted alphabetically
   *  Now when we save, the site order differs from the charge order
   *  last_charge sets up logic to renumber charge balance equations.
   */
  last_charge = -1;
  for (i = 0; i < count_unknowns; i++)
  {
    if (x[i]->type == SURFACE)
    {
      memcpy (&temp_surface.comps[count_comps], x[i]->surface_comp,
	      sizeof (struct surface_comp));

      temp_surface.comps[count_comps].master = x[i]->master[0];
      temp_surface.comps[count_comps].la = x[i]->master[0]->s->la;
      /* temp_surface.comps[count_comps].formula = NULL; */
      temp_surface.comps[count_comps].moles = 0.;
      if (x[i]->surface_comp->charge == last_charge)
      {
	temp_surface.comps[count_comps].charge = count_charge - 1;
      }
      else
      {
	temp_surface.comps[count_comps].charge = count_charge;
      }
      last_charge = x[i]->surface_comp->charge;
/*
 *   Save element concentrations on surface
 */
      count_elts = 0;
      paren_count = 0;
      charge = 0.0;
      for (j = 0; j < count_species_list; j++)
      {
	if (species_list[j].master_s == x[i]->master[0]->s)
	{
	  add_elt_list (species_list[j].s->next_elt,
			species_list[j].s->moles);
	  charge += species_list[j].s->moles * species_list[j].s->z;
	}
      }
      temp_surface.comps[count_comps].totals = elt_list_save ();
      temp_surface.comps[count_comps].formula_totals =
	elt_list_dup (x[i]->surface_comp->formula_totals);
      temp_surface.comps[count_comps].cb = charge;

      /* update unknown pointer */
      x[i]->surface_comp = &(temp_surface.comps[count_comps]);
      count_comps++;
    }
    else if (x[i]->type == SURFACE_CB && use.surface_ptr->type == DDL)
    {
      memcpy (&temp_surface.charge[count_charge], x[i]->surface_charge,
	      sizeof (struct surface_charge));
      temp_surface.charge[count_charge].charge_balance = x[i]->f;
      temp_surface.charge[count_charge].mass_water =
	x[i]->surface_charge->mass_water;
      temp_surface.charge[count_charge].diffuse_layer_totals = NULL;
      temp_surface.charge[count_charge].count_g = 0;
      temp_surface.charge[count_charge].g = NULL;
/*
 *   Added code to save g
 */

      if (x[i]->surface_charge->count_g >
	  0 /*&& use.surface_ptr->type != CD_MUSIC */ )
      {
	temp_surface.charge[count_charge].count_g =
	  x[i]->surface_charge->count_g;
	temp_surface.charge[count_charge].g =
	  (struct surface_diff_layer *) PHRQ_malloc ((size_t) x[i]->
						     surface_charge->count_g *
						     sizeof (struct
							     surface_diff_layer));
	if (temp_surface.charge[count_charge].g == NULL)
	  malloc_error ();
	memcpy (temp_surface.charge[count_charge].g, x[i]->surface_charge->g,
		(size_t) x[i]->surface_charge->count_g *
		sizeof (struct surface_diff_layer));
      }

      /*temp_surface.charge[count_charge].psi_master = x[i]->master[0]; */
      temp_surface.charge[count_charge].la_psi = x[i]->master[0]->s->la;
/*
 *   Store moles from diffuse_layer
 */
      if (dl_type_x != NO_DL)
      {
	sum_diffuse_layer (x[i]->surface_charge);
	temp_surface.charge[count_charge].diffuse_layer_totals =
	  elt_list_save ();
      }

      /* update unknown pointer */
      x[i]->surface_charge = &(temp_surface.charge[count_charge]);
      x[i]->surface_comp = x[i - 1]->surface_comp;

      count_charge++;
    }
    else if (x[i]->type == SURFACE_CB && use.surface_ptr->type == CD_MUSIC)
    {
      memcpy (&temp_surface.charge[count_charge], x[i]->surface_charge,
	      sizeof (struct surface_charge));
      if (dl_type_x != NO_DL)
      {
	temp_surface.charge[count_charge].charge_balance =
	  (x[i]->surface_charge->sigma0 +
	   x[i]->surface_charge->sigma1 +
	   x[i]->surface_charge->sigma2 +
	   x[i]->surface_charge->sigmaddl)
	  * (x[i]->surface_charge->specific_area *
	     x[i]->surface_charge->grams) / F_C_MOL;
      }
      else
      {
	temp_surface.charge[count_charge].charge_balance =
	  (x[i]->surface_charge->sigma0 +
	   x[i]->surface_charge->sigma1 +
	   x[i]->surface_charge->sigma2)
	  * (x[i]->surface_charge->specific_area *
	     x[i]->surface_charge->grams) / F_C_MOL;
      }
      temp_surface.charge[count_charge].mass_water =
	x[i]->surface_charge->mass_water;
      temp_surface.charge[count_charge].diffuse_layer_totals = NULL;
      temp_surface.charge[count_charge].count_g = 0;
      temp_surface.charge[count_charge].g = NULL;
/*
 *   Added code to save g
 */

      if (x[i]->surface_charge->count_g > 0)
      {
	temp_surface.charge[count_charge].count_g =
	  x[i]->surface_charge->count_g;
	temp_surface.charge[count_charge].g =
	  (struct surface_diff_layer *) PHRQ_malloc ((size_t) x[i]->
						     surface_charge->count_g *
						     sizeof (struct
							     surface_diff_layer));
	if (temp_surface.charge[count_charge].g == NULL)
	  malloc_error ();
	memcpy (temp_surface.charge[count_charge].g, x[i]->surface_charge->g,
		(size_t) x[i]->surface_charge->count_g *
		sizeof (struct surface_diff_layer));
      }

      /*temp_surface.charge[count_charge].psi_master = x[i]->master[0]; */
      temp_surface.charge[count_charge].la_psi = x[i]->master[0]->s->la;
/*
 *   Store moles from diffuse_layer
 */
      if (dl_type_x != NO_DL)
      {
	sum_diffuse_layer (x[i]->surface_charge);
	temp_surface.charge[count_charge].diffuse_layer_totals =
	  elt_list_save ();
      }

      /* update unknown pointer */
      x[i]->surface_charge = &(temp_surface.charge[count_charge]);
      x[i]->surface_comp = x[i - 1]->surface_comp;

      count_charge++;
    }
  }
/*
 *   Finish up
 */
  surface_ptr = surface_bsearch (n_user, &n);
  if (surface_ptr == NULL)
  {
    n = count_surface++;
    space ((void **) ((void *) &surface), count_surface, &max_surface,
	   sizeof (struct surface));
  }
  else
  {
    surface_free (&surface[n]);
  }
  memcpy (&surface[n], &temp_surface, sizeof (struct surface));
  if (n == count_surface - 1 && count_surface > 1)
  {
    if (surface[n].n_user < surface[n - 1].n_user)
    {
      qsort (surface,
	     (size_t) count_surface,
	     (size_t) sizeof (struct surface), surface_compare);
    }
  }
  use.surface_ptr = NULL;
  return (OK);
}

/* ---------------------------------------------------------------------- */
FILE *
file_open (char *query, char *default_name, const char *status, int batch)
/* ---------------------------------------------------------------------- */
{
  char name[MAX_LENGTH], answer[MAX_LENGTH];
  FILE *new_file;
  int l;

#ifdef PHREEQ98
    strcpy (name, default_name);
    if (status[0] == 'r') {
        new_file = fopen(name, "r");
    } else if (status[0] == 'w') {
        new_file = fopen(name, "w");
    }
    if (new_file == NULL) {
        	sprintf (error_string, "Can't open file, %s.", name);
            error_msg (error_string, STOP);
    }

    return (new_file);
#endif

  for (;;)
  {
/*
 *   Get file name
 */
    strcpy (name, default_name);
    if (batch == FALSE)
    {
      output_msg (OUTPUT_SCREEN, "%s\n", query);
      if (default_name[0] != '\0')
      {
	output_msg (OUTPUT_SCREEN, "Default: %s\n", default_name);
      }
      fgets (name, MAX_LENGTH, stdin);
      l = (int) strlen (name);
      name[l - 1] = '\0';
      if (name[0] == '\0')
      {
	strcpy (name, default_name);
      }
    }
/*
 *   Open existing file to read
 */
    if (status[0] == 'r')
    {
      if ((new_file = fopen (name, "r")) == NULL)
      {
	sprintf (error_string, "Can't open file, %s.", name);
	output_msg (OUTPUT_SCREEN, "\nERROR: %s\n", error_string);
	output_fflush (OUTPUT_SCREEN);
	batch = FALSE;
	continue;
      }
      break;
/*
 *   Open file to write, no check
 */
    }
    else if (status[0] == 'w')
    {
      if ((new_file = fopen (name, "w")) == NULL)
      {
	sprintf (error_string, "Error opening file, %s.", name);
	output_msg (OUTPUT_SCREEN, "\nERROR: %s\n", error_string);
	output_fflush (OUTPUT_SCREEN);
	batch = FALSE;
	continue;
      }
      break;
/*
 *   Open file to write, check for existence
 */
    }
    else if (status[0] == 'c')
    {
      for (;;)
      {
	if ((new_file = fopen (name, "r")) != NULL)
	{
	  fclose (new_file);
	  output_msg (OUTPUT_SCREEN, "Warning: File already exists, %s.\n"
		      "Enter new file name or <Enter> to overwrite:", name);
	  fgets (answer, MAX_LENGTH, stdin);
/*	  l = (int) strlen (answer);*/
	  replace("\n","\0", answer);
/*	  answer[l - 1] = '\0';*/
	  if (answer[0] != '\0')
	  {
	    strcpy (name, answer);
	    continue;
	  }
	}
	break;
      }
      if ((new_file = fopen (name, "w")) == NULL)
      {
	sprintf (error_string, "Error opening file, %s.", name);
	output_msg (OUTPUT_SCREEN, "\nERROR: %s\n", error_string);
	output_fflush (OUTPUT_SCREEN);
	batch = FALSE;
	continue;
      }
      break;
    }
  }
  strncpy (default_name, name, MAX_LENGTH);
  return (new_file);
}

/* ---------------------------------------------------------------------- */
int
copy_use (int i)
/* ---------------------------------------------------------------------- */
{
/*
 *   Find mixture
 */
  if (use.mix_in == TRUE)
  {
    mix_duplicate (use.n_mix_user, i);
  }
/*
 *   Find solution
 */
  if (use.solution_in == TRUE)
  {
    solution_duplicate (use.n_solution_user, i);
  }
/*
 *   Always save solution to i, mixing or not
 */
  save.solution = TRUE;
  save.n_solution_user = i;
  save.n_solution_user_end = i;
/*
 *   Find pure phase assemblage
 */
  if (use.pp_assemblage_in == TRUE)
  {
    pp_assemblage_duplicate (use.n_pp_assemblage_user, i);
    save.pp_assemblage = TRUE;
    save.n_pp_assemblage_user = i;
    save.n_pp_assemblage_user_end = i;
  }
  else
  {
    save.pp_assemblage = FALSE;
  }
/*
 *   Find irrev reaction
 */
  if (use.irrev_in == TRUE)
  {
    irrev_duplicate (use.n_irrev_user, i);
    save.irrev = TRUE;
    save.n_irrev_user = i;
    save.n_irrev_user_end = i;
  }
  else
  {
    save.irrev = FALSE;
  }
/*
 *   Find exchange
 */
  if (use.exchange_in == TRUE)
  {
    exchange_duplicate (use.n_exchange_user, i);
    save.exchange = TRUE;
    save.n_exchange_user = i;
    save.n_exchange_user_end = i;
  }
  else
  {
    save.exchange = FALSE;
  }
/*
 *   Find kinetics
 */
  if (use.kinetics_in == TRUE)
  {
    kinetics_duplicate (use.n_kinetics_user, i);
    save.kinetics = TRUE;
    save.n_kinetics_user = i;
    save.n_kinetics_user_end = i;
  }
  else
  {
    save.kinetics = FALSE;
  }
/*
 *   Find surface
 */
  dl_type_x = NO_DL;
  if (use.surface_in == TRUE)
  {
    surface_duplicate (use.n_surface_user, i);
    save.surface = TRUE;
    save.n_surface_user = i;
    save.n_surface_user_end = i;
  }
  else
  {
    save.surface = FALSE;
  }
/*
 *   Find temperature
 */
  if (use.temperature_in == TRUE)
  {
    temperature_duplicate (use.n_temperature_user, i);
  }
/*
 *   Find gas
 */
  if (use.gas_phase_in == TRUE)
  {
    gas_phase_duplicate (use.n_gas_phase_user, i);
    save.gas_phase = TRUE;
    save.n_gas_phase_user = i;
    save.n_gas_phase_user_end = i;
  }
  else
  {
    save.gas_phase = FALSE;
  }
/*
 *   Find solid solution
 */
  if (use.s_s_assemblage_in == TRUE)
  {
    s_s_assemblage_duplicate (use.n_s_s_assemblage_user, i);
    save.s_s_assemblage = TRUE;
    save.n_s_s_assemblage_user = i;
    save.n_s_s_assemblage_user_end = i;
  }
  else
  {
    save.s_s_assemblage = FALSE;
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
step_save_exch (int n_user)
/* ---------------------------------------------------------------------- */
{
/*
 *   Save solution composition into structure solution with user number
 *   n_user.
 *
 *   input:  n_user is user solution number of target
 */
  int i, j, k, n;
  int found;
  struct exchange *exchange_ptr;
/*
 *   Malloc space for solution data
 */
  if (use.exchange_ptr == NULL)
    return (OK);
  exchange_duplicate (use.exchange_ptr->n_user, n_user);
  exchange_ptr = exchange_bsearch (n_user, &n);
  for (i = 0; i < count_master; i++)
  {
    if (master[i]->s->type != EX)
      continue;
    found = FALSE;
    for (j = 0; j < exchange_ptr->count_comps; j++)
    {
      for (k = 0; exchange_ptr->comps[j].totals[k].elt != NULL; k++)
      {
	if (exchange_ptr->comps[j].totals[k].elt == master[i]->elt)
	{
	  if (found == FALSE)
	  {
	    found = TRUE;
	    if (master[i]->total <= MIN_TOTAL)
	    {
	      exchange_ptr->comps[j].totals[k].coef = MIN_TOTAL;
	    }
	    else
	    {
	      exchange_ptr->comps[j].totals[k].coef = master[i]->total;
	    }
	    break;
	  }
	  else
	  {
	    exchange_ptr->comps[j].totals[k].coef = 0;
	  }
	}
      }
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
step_save_surf (int n_user)
/* ---------------------------------------------------------------------- */
{
/*
 *   Save surface for intermediate calculation
 *   Amt of surface may have changed due to reaction or surface related
 *   to kinetic reactant.
 *
 *   input:  n_user is user solution number of target
 */
  int i, j, k, n, m;
  struct surface *surface_ptr;
/*
 *   Malloc space for solution data
 */
  if (use.surface_ptr == NULL)
    return (OK);
  surface_duplicate (use.surface_ptr->n_user, n_user);
  surface_ptr = surface_bsearch (n_user, &n);
  for (i = 0; i < count_master; i++)
  {
    if (master[i]->s->type != SURF)
      continue;
    for (j = 0; j < surface_ptr->count_comps; j++)
    {
      for (k = 0; surface_ptr->comps[j].totals[k].elt != NULL; k++)
      {
	if (surface_ptr->comps[j].totals[k].elt == master[i]->elt)
	{
	  if (master[i]->total <= MIN_TOTAL)
	  {
	    surface_ptr->comps[j].totals[k].coef = MIN_TOTAL;
	  }
	  else
	  {
	    surface_ptr->comps[j].totals[k].coef = master[i]->total;
	  }
	  break;
	}
      }
    }
  }
  /*
   *   Update grams
   */
  /*if (surface_ptr->edl == TRUE && surface_ptr->related_rate == TRUE && use.kinetics_ptr != NULL) { */
  if ((surface_ptr->type == DDL || surface_ptr->type == CD_MUSIC)
      && surface_ptr->related_rate == TRUE && use.kinetics_ptr != NULL)
  {
    for (j = 0; j < surface_ptr->count_comps; j++)
    {
      if (surface_ptr->comps[j].rate_name != NULL)
      {
	for (m = 0; m < use.kinetics_ptr->count_comps; m++)
	{
	  if (strcmp_nocase
	      (use.kinetics_ptr->comps[m].rate_name,
	       surface_ptr->comps[j].rate_name) != 0)
	    continue;
	  surface_ptr->charge[surface_ptr->comps[j].charge].grams =
	    use.kinetics_ptr->comps[m].m;
	  break;
	}
      }
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
copy_entities (void)
/* ---------------------------------------------------------------------- */
{
  int i, j, n, return_value;
  return_value = OK;
  if (copy_solution.count > 0)
  {
    for (j = 0; j < copy_solution.count; j++)
    {
      if (solution_bsearch (copy_solution.n_user[j], &n, FALSE) != NULL)
      {
	for (i = copy_solution.start[j]; i <= copy_solution.end[j]; i++)
	{
	  if (i == copy_solution.n_user[j])
	    continue;
	  solution_duplicate (copy_solution.n_user[j], i);
	}
      }
      else
      {
	warning_msg ("SOLUTION to copy not found.");
	return_value = ERROR;
      }
    }
  }
  if (copy_pp_assemblage.count > 0)
  {
    for (j = 0; j < copy_pp_assemblage.count; j++)
    {
      if (pp_assemblage_bsearch (copy_pp_assemblage.n_user[j], &n) != NULL)
      {
	for (i = copy_pp_assemblage.start[j]; i <= copy_pp_assemblage.end[j];
	     i++)
	{
	  if (i == copy_pp_assemblage.n_user[j])
	    continue;
	  pp_assemblage_duplicate (copy_pp_assemblage.n_user[j], i);
	}
      }
      else
      {
	warning_msg ("EQUILIBRIUM_PHASES to copy not found.");
	return_value = ERROR;
      }
    }
  }
  if (copy_irrev.count > 0)
  {
    for (j = 0; j < copy_irrev.count; j++)
    {
      if (irrev_bsearch (copy_irrev.n_user[j], &n) != NULL)
      {
	for (i = copy_irrev.start[j]; i <= copy_irrev.end[j]; i++)
	{
	  if (i == copy_irrev.n_user[j])
	    continue;
	  irrev_duplicate (copy_irrev.n_user[j], i);
	}
      }
      else
      {
	warning_msg ("REACTION to copy not found.");
	return_value = ERROR;
      }
    }
  }
  if (copy_mix.count > 0)
  {
    for (j = 0; j < copy_mix.count; j++)
    {
      if (mix_bsearch (copy_mix.n_user[j], &n) != NULL)
      {
	for (i = copy_mix.start[j]; i <= copy_mix.end[j]; i++)
	{
	  if (i == copy_mix.n_user[j])
	    continue;
	  mix_duplicate (copy_mix.n_user[j], i);
	}
      }
      else
      {
	warning_msg ("MIX to copy not found.");
	return_value = ERROR;
      }
    }
  }
  if (copy_exchange.count > 0)
  {
    for (j = 0; j < copy_exchange.count; j++)
    {
      if (exchange_bsearch (copy_exchange.n_user[j], &n) != NULL)
      {
	for (i = copy_exchange.start[j]; i <= copy_exchange.end[j]; i++)
	{
	  if (i == copy_exchange.n_user[j])
	    continue;
	  exchange_duplicate (copy_exchange.n_user[j], i);
	}
      }
      else
      {
	warning_msg ("EXCHANGE to copy not found.");
	return_value = ERROR;
      }
    }
  }
  if (copy_surface.count > 0)
  {
    for (j = 0; j < copy_surface.count; j++)
    {
      if (surface_bsearch (copy_surface.n_user[j], &n) != NULL)
      {
	for (i = copy_surface.start[j]; i <= copy_surface.end[j]; i++)
	{
	  if (i == copy_surface.n_user[j])
	    continue;
	  surface_duplicate (copy_surface.n_user[j], i);
	}
      }
      else
      {
	warning_msg ("SURFACE to copy not found.");
	return_value = ERROR;
      }
    }
  }
  if (copy_temperature.count > 0)
  {
    for (j = 0; j < copy_temperature.count; j++)
    {
      if (temperature_bsearch (copy_temperature.n_user[j], &n) != NULL)
      {
	for (i = copy_temperature.start[j]; i <= copy_temperature.end[j]; i++)
	{
	  if (i == copy_temperature.n_user[j])
	    continue;
	  temperature_duplicate (copy_temperature.n_user[j], i);
	}
      }
      else
      {
	warning_msg ("TEMPERATURE to copy not found.");
	return_value = ERROR;
      }
    }
  }
  if (copy_gas_phase.count > 0)
  {
    for (j = 0; j < copy_gas_phase.count; j++)
    {
      if (gas_phase_bsearch (copy_gas_phase.n_user[j], &n) != NULL)
      {
	for (i = copy_gas_phase.start[j]; i <= copy_gas_phase.end[j]; i++)
	{
	  if (i == copy_gas_phase.n_user[j])
	    continue;
	  gas_phase_duplicate (copy_gas_phase.n_user[j], i);
	}
      }
      else
      {
	warning_msg ("GAS_PHASE to copy not found.");
	return_value = ERROR;
      }
    }
  }
  if (copy_kinetics.count > 0)
  {
    for (j = 0; j < copy_kinetics.count; j++)
    {
      if (kinetics_bsearch (copy_kinetics.n_user[j], &n) != NULL)
      {
	for (i = copy_kinetics.start[j]; i <= copy_kinetics.end[j]; i++)
	{
	  if (i == copy_kinetics.n_user[j])
	    continue;
	  kinetics_duplicate (copy_kinetics.n_user[j], i);
	}
      }
      else
      {
	warning_msg ("KINETICS to copy not found.");
	return_value = ERROR;
      }
    }
  }
  if (copy_s_s_assemblage.count > 0)
  {
    for (j = 0; j < copy_s_s_assemblage.count; j++)
    {
      if (s_s_assemblage_bsearch (copy_s_s_assemblage.n_user[j], &n) != NULL)
      {
	for (i = copy_s_s_assemblage.start[j];
	     i <= copy_s_s_assemblage.end[j]; i++)
	{
	  if (i == copy_s_s_assemblage.n_user[j])
	    continue;
	  s_s_assemblage_duplicate (copy_s_s_assemblage.n_user[j], i);
	}
      }
      else
      {
	warning_msg ("SOLID_SOLUTIONS to copy not found.");
	return_value = ERROR;
      }
    }
  }
  copy_solution.count = 0;
  copy_pp_assemblage.count = 0;
  copy_exchange.count = 0;
  copy_surface.count = 0;
  copy_s_s_assemblage.count = 0;
  copy_gas_phase.count = 0;
  copy_kinetics.count = 0;
  copy_mix.count = 0;
  copy_irrev.count = 0;
  copy_temperature.count = 0;
  new_copy = FALSE;
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
read_database (PFN_READ_CALLBACK pfn, void *cookie)
/* ---------------------------------------------------------------------- */
{
  int errors;
  simulation = 0;

/*
 *   Prepare error handling
 */
  errors = setjmp (mark);
  if (errors != 0)
  {
    return errors;
  }

  set_read_callback (pfn, cookie, TRUE);
  dup_print ("Reading data base.", TRUE);
  read_input ();
  tidy_model ();
  status (0, NULL);
  return 0;
}

/* ---------------------------------------------------------------------- */
int
run_simulations (PFN_READ_CALLBACK pfn, void *cookie)
/* ---------------------------------------------------------------------- */
{
  int errors;
  char token[MAX_LENGTH];

/*
 *   Prepare error handling
 */
  errors = setjmp (mark);
  if (errors != 0)
  {
    // MDL: output if error occurs
    printf(":: run_simulation setjmp errors %i\n",errors);
    return errors;
  }

  set_read_callback (pfn, cookie, FALSE);

/*
 *   Read input data for simulation
 */
// MDL: added condition simulation <2 (just the first simulation!)
//      to avoid problems with punching
  for (simulation = 1; simulation <2; simulation++)
  {

#ifdef PHREEQ98
    AddSeries = !connect_simulations;
#endif
    sprintf (token, "Reading input data for simulation %d.", simulation);

    output_msg (OUTPUT_GUI_ERROR, "\nSimulation %d\n", simulation);

    dup_print (token, TRUE);
    if (read_input () == EOF)
      break;

    if (title_x != NULL)
    {
      sprintf (token, "TITLE");
      dup_print (token, TRUE);
      if (pr.headings == TRUE)
	output_msg (OUTPUT_MESSAGE, "%s\n\n", title_x);
    }
    tidy_model ();
#ifdef PHREEQC_CPP
    /*test_classes(); */
#endif
#ifdef PHREEQ98
    if (!phreeq98_debug)
    {
#endif

/*
 *   Calculate distribution of species for initial solutions
 */
      if (new_solution)
	initial_solutions (TRUE);
/*
 *   Calculate distribution for exchangers
 */
      if (new_exchange)
	initial_exchangers (TRUE);
/*
 *   Calculate distribution for surfaces
 */
      if (new_surface)
	initial_surfaces (TRUE);
/*
 *   Calculate initial gas composition
 */
      if (new_gas_phase)
	initial_gas_phases (TRUE);
/*
 *   Calculate reactions
 */
      reactions ();
/*
 *   Calculate inverse models
 */
      inverse_models ();
/*
 *   Calculate advection
 */
      if (use.advect_in == TRUE)
      {
	dup_print ("Beginning of advection calculations.", TRUE);
	advection ();
      }
/*
 *   Calculate transport
 */
      if (use.trans_in == TRUE)
      {
	dup_print ("Beginning of transport calculations.", TRUE);
	transport ();
      }
/*
 *   Copy
 */
      if (new_copy)
	copy_entities ();
/*
 *   End of simulation
 */
      dup_print ("End of simulation.", TRUE);
#ifdef PHREEQ98
    }				/* if (!phreeq98_debug) */
#endif
#ifdef SKIP

    {
      int n;
      SAX_StartSystem ();
      for (n = 0; n < count_solution; ++n)
      {
	SAX_AddSolution (solution[n]);
      }
      SAX_EndSystem ();
      SAX_UnpackSolutions (SAX_GetXMLStr (), SAX_GetXMLLength ());
    }
#endif

  }


  return 0;
}

/* ---------------------------------------------------------------------- */
int
do_initialize (void)
/* ---------------------------------------------------------------------- */
{
  int errors;
/*
 *   Prepare error handling
 */
  errors = setjmp (mark);
  if (errors != 0)
  {
    return errors;
  }

  state = INITIALIZE;
  initialize ();

  return 0;
}

/* ---------------------------------------------------------------------- */
int
do_status (void)
/* ---------------------------------------------------------------------- */
{
  int errors;
/*
 *   Prepare error handling
 */
  errors = setjmp (mark);
  if (errors != 0)
  {
    return errors;
  }

  if (pr.status == TRUE)
  {
#if defined(PHREEQCI_GUI)
    state = -1;
    status (0, "\r\nDone.");
#else
    status (0, "\nDone.");
#endif
    output_msg (OUTPUT_SCREEN, "\n");
    output_msg (OUTPUT_SEND_MESSAGE, "\r\n");
  }
  dup_print ("End of run.", TRUE);
  output_msg (OUTPUT_SCREEN, "\nEnd of Run.\n");
  output_msg (OUTPUT_SEND_MESSAGE, "\r\nEnd of Run.\r\n");

  return 0;
}



/* ---------------------------------------------------------------------- */
int
P_run_simulations (PFN_READ_CALLBACK pfn, void *cookie)
/* ---------------------------------------------------------------------- */
{
  int errors, i;
  char token[MAX_LENGTH];

/*
 *   Prepare error handling
 */
  errors = setjmp (mark);
  if (errors != 0)
  {
    printf(":: R_run_simulation setjmp errors = %i\n",errors);
    return (errors);
  }


#ifdef MDL_DEBUG   
  printf(": deb ::  in P_run_simulation, Pnsim=%i\n",Pnsim);
#endif

  // MDL: alloc the total list
  set_read_callback (pfn, cookie, FALSE);


#ifdef MDL_DEBUG   
    printf(": deb ::  after set_read_callback\n");
#endif


/*
 *   Read input data for simulation
 */
  // MDL: added condition 
  for (simulation = 1; simulation <= Pnsim; simulation++)
  { 

#ifdef MDL_DEBUG   
    printf("---> : deb ::  simulation      %i!\n",simulation);
#endif
#ifdef PHREEQ98
    AddSeries = !connect_simulations;
#endif
    sprintf (token, "Reading input data for simulation %d.", simulation);

    output_msg (OUTPUT_GUI_ERROR, "\nSimulation %d\n", simulation);

    dup_print (token, TRUE);

    // MDL: reads the input
    if (read_input () == EOF)
      break;

    // MDL: simulation gets the number of simulations, not more
    //    printf(":: Simulation %i\n",simulation);


    if (title_x != NULL)
    {
      sprintf (token, "TITLE");
      dup_print (token, TRUE);
      if (pr.headings == TRUE)
	output_msg (OUTPUT_MESSAGE, "%s\n\n", title_x);
    }
    tidy_model ();
#ifdef PHREEQC_CPP
    /*test_classes(); */
#endif
#ifdef PHREEQ98
    if (!phreeq98_debug)
    {
#endif

/*
 *   Calculate distribution of species for initial solutions
 */
      if (new_solution)
	initial_solutions (TRUE);
/*
 *   Calculate distribution for exchangers
 */
      if (new_exchange)
	initial_exchangers (TRUE);
/*
 *   Calculate distribution for surfaces
 */
      if (new_surface)
	initial_surfaces (TRUE);
/*
 *   Calculate initial gas composition
 */
      if (new_gas_phase)
	initial_gas_phases (TRUE);
/*
 *   Calculate reactions
 */
      reactions ();

#ifdef MDL_DEBUG   
    printf(": deb :: reactions() OK!\n");
#endif

/*
 *   Calculate inverse models
 */
      inverse_models ();
/*
 *   Calculate advection
 */
      if (use.advect_in == TRUE)
      {
	dup_print ("Beginning of advection calculations.", TRUE);
	advection ();
      }
/*
 *   Calculate transport
 */
      if (use.trans_in == TRUE)
      {
	dup_print ("Beginning of transport calculations.", TRUE);
	transport ();
      }
/*
 *   Copy
 */
      if (new_copy)
	copy_entities ();

#ifdef MDL_DEBUG   
    printf(": deb ::  after copy entities... \n");
#endif

/*
 *   End of simulation
 */
      dup_print ("End of simulation.", TRUE);
#ifdef PHREEQ98
    }				/* if (!phreeq98_debug) */
#endif
#ifdef SKIP
    
    {
      int n;
      SAX_StartSystem ();
      for (n = 0; n < count_solution; ++n)
      {
	SAX_AddSolution (solution[n]);
      }
      SAX_EndSystem ();
      SAX_UnpackSolutions (SAX_GetXMLStr (), SAX_GetXMLLength ());
    }
#endif

    // MDL: save this simulation
#ifdef MDL_DEBUG
    printf(": deb :: R_sol_punch called, sim=%i, Ppunch_dim %i\n",simulation,Ppunch_dim);
#endif
    // MDL: Pdo_punch = TRUE means we are calling this punch
    Pdo_punch = TRUE;
    punch.user_punch = TRUE;
    punch_user_punch ();
    punch.user_punch = FALSE;
    Pdo_punch = FALSE;

    // MDL: copy the vector from the global pointer
    for (i=0; i<Ppunch_dim; i++)
      {
	Pout[(simulation-1)*Ppunch_dim+i] = Ppunch[i];
      }
#ifdef MDL_DEBUG   
    printf(": deb :: mainsub: P_do_out [OK] \n");
#endif
  }

  return(OK);
}
