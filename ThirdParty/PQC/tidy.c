#define EXTERNAL extern
#include "global.h"
#include "phqalloc.h"
#include "output.h"
#include "phrqproto.h"

static char const svnid[] = "$Id: tidy.c 4 2009-04-21 17:29:29Z delucia $";

static int check_species_input (void);
static LDBLE coef_in_master (struct master *master_ptr);
static int phase_rxn_to_trxn (struct phase *phase_ptr,
			      struct reaction *rxn_ptr);
static int reset_last_model (void);
static int rewrite_eqn_to_primary (void);
static int rewrite_eqn_to_secondary (void);
static int species_rxn_to_trxn (struct species *s_ptr);
static int tidy_logk (void);
static int tidy_min_exchange (void);
static int tidy_kin_exchange (void);
static int tidy_gas_phase (void);
static int tidy_inverse (void);
static int tidy_isotopes (void);
static int tidy_isotope_ratios (void);
static int tidy_isotope_alphas (void);
static int tidy_kin_surface (void);
static int tidy_master_isotope (void);
static int tidy_min_surface (void);
static int tidy_phases (void);
static int tidy_pp_assemblage (void);
static int tidy_solutions (void);
static int tidy_s_s_assemblage (void);
static int tidy_species (void);
static int tidy_surface (void);

static LDBLE a0, a1, kc, kb;
static int scan (LDBLE f (LDBLE x), LDBLE * xx0, LDBLE * xx1);
static LDBLE halve (LDBLE f (LDBLE x), LDBLE x0, LDBLE x1, LDBLE tol);
static LDBLE f_spinodal (LDBLE x);
static int solve_misc (LDBLE * xxc1, LDBLE * xxc2, LDBLE tol);
static int s_s_calc_a0_a1 (struct s_s *s_s_ptr);
#define ZERO_TOL 1.0e-30

/* ---------------------------------------------------------------------- */
int
tidy_model (void)
/* ---------------------------------------------------------------------- */
{
  int i, j;
  int n_user, last;
  int new_named_logk;
  if (svnid == NULL)
    fprintf (stderr, " ");
  /*
   * Determine if any new elements, species, phases have been read
   */
  new_model = FALSE;
  new_pp_assemblage = FALSE;
  new_surface = FALSE;
  new_exchange = FALSE;
  new_reaction = FALSE;
  new_temperature = FALSE;
  new_mix = FALSE;
  new_solution = FALSE;
  new_gas_phase = FALSE;
  new_inverse = FALSE;
  new_punch = FALSE;
  new_surface = FALSE;
  new_s_s_assemblage = FALSE;
  new_kinetics = FALSE;
  new_pitzer = FALSE;
  new_named_logk = FALSE;
  if (keyword[2].keycount > 0 ||	/*"species" */
      keyword[3].keycount > 0 ||	/*"master" */
      keyword[5].keycount > 0 ||	/*"phases" */
      keyword[11].keycount > 0 ||	/*"exchange_species" */
      keyword[12].keycount > 0 ||	/*"master_exchange_species" */
      keyword[14].keycount > 0 ||	/*"surface_species" */
      keyword[15].keycount > 0 ||	/*"master_surface_species" */
      keyword[36].keycount > 0 ||	/*"rates" */
      keyword[47].keycount > 0 ||	/*"llnl_aqueous_model_parameters" */
      (keyword[49].keycount > 0 && simulation == 0) ||	/*"database" */
      keyword[51].keycount > 0 ||	/*"named_analytical_expressions" */
      keyword[54].keycount > 0 ||	/*"isotopes" */
      keyword[55].keycount > 0 ||	/*"calculate_values" */
      keyword[56].keycount > 0 ||	/*"isotopes_ratios", */
      keyword[57].keycount > 0 ||	/*"isotopes_alphas" */
      keyword[59].keycount > 0)
  {				/*"pitzer" */
    new_model = TRUE;
  }
  if (keyword[6].keycount > 0)
    new_pp_assemblage = TRUE;	/*"pure_phases" */
  if (keyword[16].keycount > 0)
    new_surface = TRUE;		/*"surface" */
  if (keyword[13].keycount > 0)
    new_exchange = TRUE;	/*"exchange" */
  if (keyword[7].keycount > 0)
    new_reaction = TRUE;	/*"reaction" */
  if (keyword[17].keycount > 0)
    new_temperature = TRUE;	/*"reacton_temperature" */
  if (keyword[8].keycount > 0)
    new_mix = TRUE;		/*"mix" */
  if (keyword[4].keycount > 0 ||	/*"solution" */
      keyword[43].keycount > 0)
  {				/*"spread_solution" */
    new_solution = TRUE;
  }
  if (keyword[19].keycount > 0)
    new_gas_phase = TRUE;	/*"gas_phase" */
  if (keyword[18].keycount > 0)
    new_inverse = TRUE;		/*"inverse_modeling" */
  if (keyword[22].keycount > 0 ||	/*"selected_output" */
      keyword[39].keycount > 0)
  {				/*"user_punch" */
    new_punch = TRUE;
  }
  if (keyword[40].keycount > 0)
    new_s_s_assemblage = TRUE;	/*"solid_solutions" */
  if (keyword[33].keycount > 0)
    new_kinetics = TRUE;	/*"kinetics" */
  if (keyword[58].keycount > 0)
    new_copy = TRUE;		/*"copy" */
  if (keyword[59].keycount > 0)
    new_pitzer = TRUE;		/*"pitzer" */
  if (keyword[50].keycount > 0 ||
      keyword[51].keycount > 0 ||
      keyword[52].keycount > 0 || keyword[53].keycount > 0)
    new_named_logk = TRUE;	/*"named_log_k" */

  /*
     0      "eof"
     1      "end"
     2      "species"
     3      "master"
     4      "solution"
     5       "phases"
     6       "pure_phases"
     7       "reaction"
     8       "mix"
     9       "use"
     10      "save"
     11      "exchange_species"
     12      "master_exchange_species"
     13      "exchange"
     14      "surface_species"
     15      "master_surface_species"
     16      "surface"
     17      "reacton_temperature"
     18      "inverse_modeling"
     19      "gas_phase"
     20      "transport"
     21      "debug"
     22      "selected_output"
     23      "select_output"
     24      "knobs"
     25      "print"
     26      "equilibrium_phases"  
     27      "equilibria"
     28      "equilibrium"         
     29      "pure"                
     30      "title"
     31      "comment"
     32      "advection"
     33      "kinetics"
     34      "incremental_reactions"
     35      "incremental"
     36      "rates"
     37      "solution_s"
     38      "user_print"
     39      "user_punch"
     40      "solid_solutions"
     41      "solid_solution"
     42      "solution_spread"
     43      "spread_solution"
     44      "selected_out"
     45      "select_out"
     46      "user_graph"
     47      "llnl_aqueous_model_parameters"
     48      "llnl_aqueous_model"
     49      "database"
     50      "named_analytical_expression"
     51      "named_analytical_expressions"
     52      "named_expressions"
     53      "named_log_k"
     54      "isotopes"
     55      "calculate_values"
     56      "isotopes_ratios",
     57      "isotopes_alphas"
     58      "copy"
     59      "pitzer"
   */

/*
 *   Sort arrays
 */

/* species */
  if (new_model == TRUE)
  {
    qsort (s,
	   (size_t) count_s, (size_t) sizeof (struct species *), s_compare);

/* master species */
    qsort (master,
	   (unsigned) count_master, sizeof (struct master *), master_compare);

/* elements */
    qsort (elements,
	   (size_t) count_elements,
	   (size_t) sizeof (struct element *), element_compare);
/* phases */
    qsort (phases,
	   (size_t) count_phases,
	   (size_t) sizeof (struct phase *), phase_compare);

  }
/* pure_phases */
  if (new_pp_assemblage)
    pp_assemblage_sort ();

/* solid solutions */
  if (new_s_s_assemblage)
    s_s_assemblage_sort ();

/* exchangers */
  if (new_exchange)
    exchange_sort ();
/* surfaces */
  if (new_surface)
    surface_sort ();
#ifdef SKIP
/* mixtures */
  qsort (mix, (size_t) count_mix, (size_t) sizeof (struct mix), mix_compare);
#endif

/* mixtures */
/* !!!!!
 * In transport mode, with stagnant cells, cell number is 'n_mix_user' for
 *  both mix with stagnant cells and for dispersive mix. 
 *  qsort(mix) then fails. Hence ...
 */

  if ((state != TRANSPORT) || (simul_tr < 2) || (stag_data->count_stag == 0))
  {
    mix_sort ();
  }

/* gas_phase */
  if (new_gas_phase)
    gas_phase_sort ();

/* kinetics */
  if (new_kinetics)
    kinetics_sort ();

  /* named_log_k */
  if (new_named_logk)
    tidy_logk ();
/*
 *   Check pointers, write reactions for species
 */
  if (new_model)
  {
    tidy_species ();

    tidy_phases ();

    tidy_master_isotope ();
/*
 *   calculate gfw of water, kg/mole
 */
    compute_gfw ("H2O", &gfw_water);
    gfw_water *= 0.001;
  }
/*
 *   tidy surface data
 */
  if (new_model || new_surface)
    tidy_surface ();
/*
 *   tidy inverse data
 */
  if (new_inverse)
    tidy_inverse ();
/*
 *   tidy gas phase data
 */
  if (new_gas_phase)
    tidy_gas_phase ();
/*
 *   tidy pp_assemblage data
 */
  if (new_model || new_pp_assemblage)
    tidy_pp_assemblage ();
/*
 *   tidy s_s_assemblage data
 */
  if (new_model || new_s_s_assemblage)
    tidy_s_s_assemblage ();
/*
 *   tidy exchange data, after pp_assemblages
 */
  if (new_exchange)
    tidy_min_exchange ();
  if (new_exchange)
    tidy_kin_exchange ();
/*
 *   tidy surface data
 */
  if (new_surface)
    tidy_min_surface ();
  if (new_surface)
    tidy_kin_surface ();
/*
 *   tidy solution isotope data
 */
  if (new_solution)
    tidy_isotopes ();
  if (new_model)
    tidy_isotope_ratios ();
  if (new_model)
    tidy_isotope_alphas ();
/*
 *   Duplicate reaction
 */
  if (new_reaction)
  {
    for (i = 0; i < count_irrev; i++)
    {
      if (irrev[i].n_user_end > irrev[i].n_user)
      {
	n_user = irrev[i].n_user;
	last = irrev[i].n_user_end;
	irrev[i].n_user_end = irrev[i].n_user;
	for (j = n_user + 1; j <= last; j++)
	{
	  irrev_duplicate (n_user, j);
	}
      }
    }
  }
/*
 *   Duplicate kinetics
 */
  if (new_kinetics)
  {
    for (i = 0; i < count_kinetics; i++)
    {
      if (kinetics[i].n_user_end > kinetics[i].n_user)
      {
	n_user = kinetics[i].n_user;
	last = kinetics[i].n_user_end;
	kinetics[i].n_user_end = kinetics[i].n_user;
	for (j = n_user + 1; j <= last; j++)
	{
	  kinetics_duplicate (n_user, j);
	}
      }
    }
  }
/*
 *   Duplicate temperature
 */
  if (new_temperature)
  {
    temperature_sort ();
    for (i = 0; i < count_temperature; i++)
    {
      if (temperature[i].n_user_end > temperature[i].n_user)
      {
	n_user = temperature[i].n_user;
	last = temperature[i].n_user_end;
	temperature[i].n_user_end = temperature[i].n_user;
	for (j = n_user + 1; j <= last; j++)
	{
	  temperature_duplicate (n_user, j);
	}
      }
    }
  }
/*
 *   Tidy pitzer information
 */
  if (pitzer_model && new_model)
    pitzer_tidy ();
/*
 *   Tidy punch information
 */
  if (input_error == 0 && (new_punch || new_model))
    tidy_punch ();
/*
 *   Tidy solution information
 */
  if (new_solution)
    tidy_solutions ();

  /*      if (new_model || new_exchange || new_pp_assemblage || new_surface || new_gas_phase || new_kinetics) reset_last_model(); */
  if (new_model)
    reset_last_model ();
/*
 *   make sure essential species are defined
 */
  if (s_h2o == NULL)
  {
    input_error++;
    error_msg ("H2O not defined.", STOP);
  }
  if (s_h2o->primary == NULL)
  {
    input_error++;
    error_msg ("H2O, primary master species for O, not defined.", CONTINUE);
  }
  if (s_h2o->secondary == NULL)
  {
    input_error++;
    error_msg ("H2O, secondary master species for O(-2), not defined.",
	       CONTINUE);
  }
  if (s_hplus == NULL && s_h3oplus == NULL)
  {
    input_error++;
    error_msg ("Neither H+ nor H3O+ are defined in solution_species.", STOP);
  }
  else if (s_hplus == NULL && s_h3oplus != NULL)
  {
    s_hplus = s_h3oplus;
    s_h3oplus = NULL;
  }
  else if (s_hplus != NULL && s_h3oplus == NULL)
  {
  }
  else if (s_hplus != NULL && s_h3oplus != NULL)
  {
    input_error++;
    error_msg ("Can not define both H+ and H3O+ in solution_species.", STOP);
  }
  if (s_hplus->primary == NULL)
  {
    input_error++;
    error_msg ("H3O+, primary master species for H, not defined.", CONTINUE);
  }
  if (s_hplus->secondary == NULL)
  {
    input_error++;
    error_msg ("H3O+, secondary master species for H(1), not defined.",
	       CONTINUE);
  }
  if (s_eminus == NULL)
  {
    input_error++;
    error_msg ("e- not defined in solution_species.", CONTINUE);
  }
  if (s_eminus->primary == NULL)
  {
    input_error++;
    error_msg ("e-, primary master species for E-, not defined.", CONTINUE);
  }
  if (pitzer_model == FALSE || pitzer_pe == TRUE)
  {
    if (s_h2 == NULL)
    {
      input_error++;
      error_msg ("H2(aq) not defined in solution_species.", CONTINUE);
    }
    if (s_o2 == NULL)
    {
      input_error++;
      error_msg ("O2(aq) not defined in solution_species.", CONTINUE);
    }
  }
  element_h_one = element_store ("H(1)");
  if (element_h_one == NULL)
  {
    input_error++;
    error_msg ("H(1) not defined in solution_master_species.", CONTINUE);
  }
/*
 *   Error check, program termination
 */
  if (input_error > 0 || parse_error > 0)
  {
    error_msg ("Calculations terminating due to input errors.", STOP);
  }

  return (OK);
}

/* ---------------------------------------------------------------------- */
int
check_species_input (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Check species data for completeness
 */
  int i;
  int return_value;

  return_value = OK;
  for (i = 0; i < count_s; i++)
  {
    if (s[i]->next_elt == NULL)
    {
      input_error++;
      return_value = ERROR;
      sprintf (error_string,
	       "Elements in species have not been tabulated, %s.",
	       s[i]->name);
      error_msg (error_string, CONTINUE);
    }
    if (s[i]->rxn == NULL)
    {
      input_error++;
      return_value = ERROR;
      sprintf (error_string, "Reaction for species has not been defined, %s.",
	       s[i]->name);
      error_msg (error_string, CONTINUE);
    }
    else
    {
      select_log_k_expression (s[i]->logk, s[i]->rxn->logk);
      add_other_logk (s[i]->rxn->logk, s[i]->count_add_logk, s[i]->add_logk);
    }
  }
  return (return_value);
}

/* ---------------------------------------------------------------------- */
int
select_log_k_expression (LDBLE * source_k, LDBLE * target_k)
/* ---------------------------------------------------------------------- */
{
  int j, analytic;

  analytic = FALSE;
  for (j = 2; j < 7; j++)
  {
    if (source_k[j] != 0.0)
    {
      analytic = TRUE;
      break;
    }
  }
  if (analytic == TRUE)
  {
    target_k[0] = 0.0;
    target_k[1] = 0.0;
    for (j = 2; j < 7; j++)
    {
      target_k[j] = source_k[j];
    }
  }
  else
  {
    target_k[0] = source_k[0];
    target_k[1] = source_k[1];
    for (j = 2; j < 7; j++)
    {
      target_k[j] = 0.0;
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
tidy_logk (void)
/* ---------------------------------------------------------------------- */
/*
 *  Picks log k expression
 */
{
  int i;
  for (i = 0; i < count_logk; i++)
  {
    select_log_k_expression (logk[i]->log_k_original, logk[i]->log_k);
    logk[i]->done = FALSE;
  }
  for (i = 0; i < count_logk; i++)
  {
    if (logk[i]->done == FALSE)
    {
      add_logks (logk[i], 0);
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
add_other_logk (LDBLE * source_k, int count_add_logk,
		struct name_coef *add_logk)
/* ---------------------------------------------------------------------- */
{
  int i, j, analytic;
  struct logk *logk_ptr;
  char token[MAX_LENGTH];
  LDBLE coef;
  ENTRY item, *found_item;

  if (count_add_logk == 0)
    return (OK);
  for (i = 0; i < count_add_logk; i++)
  {
    coef = add_logk[i].coef;
    strcpy (token, add_logk[i].name);
    str_tolower (token);
    item.key = token;
    item.data = NULL;
    found_item = hsearch_multi (logk_hash_table, item, FIND);
    if (found_item == NULL)
    {
      input_error++;
      sprintf (error_string,
	       "Could not find named temperature expression, %s\n", token);
      error_msg (error_string, CONTINUE);
      return (ERROR);
    }
    logk_ptr = (struct logk *) found_item->data;
    analytic = FALSE;
    for (j = 2; j < 7; j++)
    {
      if (logk_ptr->log_k[j] != 0.0)
      {
	analytic = TRUE;
	break;
      }
    }
    if (analytic == TRUE)
    {
      for (j = 2; j < 7; j++)
      {
	source_k[j] += logk_ptr->log_k[j] * coef;
      }
    }
    else
    {
      source_k[0] += logk_ptr->log_k[0] * coef;
      source_k[1] += logk_ptr->log_k[1] * coef;
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
add_logks (struct logk *logk_ptr, int repeats)
/* ---------------------------------------------------------------------- */
{
  int i, j;
  struct logk *next_logk_ptr;
  char token[MAX_LENGTH];
  LDBLE coef;
  ENTRY item, *found_item;
  /*
   *  Adds in other named_expressions to get complete log K
   *  Evaluates others recursively if necessary
   */
  if (repeats > 15)
  {
    input_error++;
    sprintf (error_string, "Circular definition of named_logk? %s\n",
	     logk_ptr->name);
    error_msg (error_string, CONTINUE);
    return (ERROR);
  }
  for (i = 0; i < logk_ptr->count_add_logk; i++)
  {
    coef = logk_ptr->add_logk[i].coef;
    strcpy (token, logk_ptr->add_logk[i].name);
    str_tolower (token);
    item.key = token;
    item.data = NULL;
    found_item = hsearch_multi (logk_hash_table, item, FIND);
    if (found_item == NULL)
    {
      input_error++;
      sprintf (error_string,
	       "Could not find named temperature expression, %s\n", token);
      error_msg (error_string, CONTINUE);
      return (ERROR);
    }
    next_logk_ptr = (struct logk *) found_item->data;
    if (next_logk_ptr->done == FALSE)
    {
      /*output_msg(OUTPUT_MESSAGE, "Done == FALSE\n", token); */
      if (add_logks (next_logk_ptr, repeats + 1) == ERROR)
      {
	return (ERROR);
      }
    }
    for (j = 0; j < 7; j++)
    {
      logk_ptr->log_k[j] += next_logk_ptr->log_k[j] * coef;
    }
  }
  logk_ptr->done = TRUE;
  return (OK);
}

/* ---------------------------------------------------------------------- */
LDBLE
coef_in_master (struct master * master_ptr)
/* ---------------------------------------------------------------------- */
{
  int l;
  LDBLE coef;
  char *ptr;
  char elt_name[MAX_LENGTH];
  struct elt_list *next_elt;

  coef = 0.0;
  ptr = master_ptr->elt->name;
  get_elt (&ptr, elt_name, &l);
  for (next_elt = master_ptr->s->next_elt; next_elt->elt != NULL; next_elt++)
  {
    if (strcmp (elt_name, next_elt->elt->name) == 0)
    {
      coef = next_elt->coef;
      break;
    }
  }
  return (coef);
}

/* ---------------------------------------------------------------------- */
int
rewrite_eqn_to_secondary (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Write equation for species in terms of secondary species
 *   Result is in trxn.
 */
  LDBLE coef;
  int repeat, i, add_count;
  struct rxn_token_temp *token_ptr;
/*
 *
 */
  add_count = 0;
  repeat = TRUE;
/*
 *   Reduce reaction equation to primary and secondary species
 */
  while (repeat == TRUE)
  {
    repeat = FALSE;
    /*   Check for too many iterations */
    if (++add_count > MAX_ADD_EQUATIONS)
    {
      parse_error++;
      sprintf (error_string,
	       "Could not reduce equation to secondary master species, %s.",
	       trxn.token[0].name);
      error_msg (error_string, CONTINUE);
      break;
    }

    for (i = 1; i < count_trxn; i++)
    {
      token_ptr = &(trxn.token[i]);
      if (token_ptr->s == NULL)
      {
	sprintf (error_string, "NULL species pointer for species, %s.",
		 token_ptr->name);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      if (token_ptr->s->secondary == NULL && token_ptr->s->primary == NULL)
      {
	coef = token_ptr->coef;
	trxn_add (token_ptr->s->rxn, coef, TRUE);
	repeat = TRUE;
	break;
      }
    }
  }
  trxn_combine ();
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
replace_solids_gases (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Write equation for species in terms of secondary species
 *   Result is in trxn.
 */
  LDBLE coef;
  int n;
  int repeat, i, add_count;
  struct rxn_token_temp *token_ptr;
  struct phase *phase_ptr;
  int replaced;
  char token[MAX_LENGTH];
/*
 *
 */
  add_count = 0;
  repeat = TRUE;
  replaced = FALSE;
/*
 *   Reduce reaction equation to primary and secondary species
 */
  while (repeat == TRUE)
  {
    repeat = FALSE;
    /*   Check for too many iterations */
    if (++add_count > MAX_ADD_EQUATIONS)
    {
      parse_error++;
      sprintf (error_string,
	       "Could not remove all solids and gases from equation, %s.",
	       trxn.token[0].name);
      error_msg (error_string, CONTINUE);
      break;
    }

    for (i = 1; i < count_trxn; i++)
    {
      token_ptr = &(trxn.token[i]);
      if (token_ptr->s == NULL)
      {
	phase_ptr = phase_bsearch (token_ptr->name, &n, FALSE);
	/* try phase name without (g) or  (s) */
	if (phase_ptr == NULL)
	{
	  strcpy (token, token_ptr->name);
	  replace ("(g)", "", token);
	  replace ("(s)", "", token);
	  replace ("(G)", "", token);
	  replace ("(S)", "", token);
	  phase_ptr = phase_bsearch (token, &n, FALSE);
	}
	if (phase_ptr == NULL)
	{
	  input_error++;
	  sprintf (error_string, "Phase not found, %s.", token_ptr->name);
	  error_msg (error_string, CONTINUE);
	  break;
	}
	coef = token_ptr->coef;
	/* add reaction for solid/gas */
	/* debug
	   output_msg(OUTPUT_MESSAGE, "Reaction to add.\n");
	   rxn_print(phase_ptr->rxn);
	 */
	trxn_add_phase (phase_ptr->rxn, coef, FALSE);

	/* remove solid/gas from trxn list */
	trxn.token[i].name = phase_ptr->rxn->token[0].name;
	trxn.token[i].s = phase_ptr->rxn->token[0].s;
	trxn.token[i].coef = -coef * phase_ptr->rxn->token[0].coef;
	repeat = TRUE;
	replaced = TRUE;
	/* debug
	   output_msg(OUTPUT_MESSAGE, "Before combined.\n");
	   trxn_print();
	 */
	/* combine */
	trxn_combine ();
	/* debug
	   output_msg(OUTPUT_MESSAGE, "Combined.\n");
	   trxn_print();
	 */
	break;
      }
    }
  }
  trxn_combine ();
  return (replaced);
}

/* ---------------------------------------------------------------------- */
int
rewrite_eqn_to_primary (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Write equation for secondary master species in terms of primary master species
 *   Store result in reaction structure for master species
 *   rewrite if necessary.
 *
 */
  int repeat, j, add_count;

/*
 *   Check secondary master species
 */
  repeat = TRUE;
  add_count = 0;
/*
 *   Check if reaction contains only primary master species
 */
  while (repeat == TRUE)
  {
    repeat = FALSE;
/*
 *   Check for too many iterations
 */
    if (++add_count > MAX_ADD_EQUATIONS)
    {
      parse_error++;
      sprintf (error_string,
	       "Could not reduce equation to primary master species, %s.",
	       trxn.token[0].s->name);
      error_msg (error_string, CONTINUE);
      break;
    }
/*
 *   Go through species in reaction for secondary master species, look for non-primary
 *   species as reactants, rewrite
 */
    for (j = 1; j < count_trxn; j++)
    {
      if (trxn.token[j].s->primary == NULL)
      {
	trxn_add (trxn.token[j].s->rxn, trxn.token[j].coef, TRUE);
	repeat = TRUE;
	break;
      }
    }
  }
  trxn_combine ();
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
tidy_gas_phase (void)
/* ---------------------------------------------------------------------- */
{
  int i, j, k, n_user, last;
  struct phase *phase_ptr;
/*
 *   Find all gases for each gas_phase in phase list
 */
  for (i = 0; i < count_gas_phase; i++)
  {
    if (gas_phase[i].new_def != TRUE)
      continue;
    gas_phase[i].new_def = FALSE;
    for (j = 0; j < gas_phase[i].count_comps; j++)
    {
      phase_ptr = phase_bsearch (gas_phase[i].comps[j].name, &k, FALSE);
      if (phase_ptr == NULL)
      {
	input_error++;
	sprintf (error_string, "Gas not found in PHASES data base, %s.",
		 gas_phase[i].comps[j].name);
	error_msg (error_string, CONTINUE);
	continue;
      }
      else
      {
	gas_phase[i].comps[j].phase = phase_ptr;
      }
/*
 *   Fixed pressure
 */
      if (gas_phase[i].type == PRESSURE)
      {
	if (gas_phase[i].solution_equilibria == TRUE)
	{
	  input_error++;
	  sprintf (error_string,
		   "Gas phase %d: can not use '-equilibrium' option with fixed pressure gas phase.",
		   gas_phase[i].n_user);
	  error_msg (error_string, CONTINUE);
	}
	/* calculate moles */
	if (gas_phase[i].comps[j].p_read != NAN)
	{
	  gas_phase[i].comps[j].moles = gas_phase[i].comps[j].p_read *
	    gas_phase[i].volume / R_LITER_ATM / gas_phase[i].temperature;
	}
	else
	{
	  input_error++;
	  sprintf (error_string,
		   "Gas phase %d: partial pressure of gas component %s not defined.",
		   gas_phase[i].n_user, gas_phase[i].comps[j].name);
	  error_msg (error_string, CONTINUE);
	}
      }
      else
      {
/*
 *   Fixed volume
 */
	if (gas_phase[i].solution_equilibria == FALSE)
	{
	  if (gas_phase[i].comps[j].p_read != NAN)
	  {
	    gas_phase[i].comps[j].moles =
	      gas_phase[i].comps[j].p_read * gas_phase[i].volume /
	      R_LITER_ATM / gas_phase[i].temperature;
	  }
	  else
	  {
	    input_error++;
	    sprintf (error_string,
		     "Gas phase %d: moles of gas component %s not defined.",
		     gas_phase[i].n_user, gas_phase[i].comps[j].name);
	    error_msg (error_string, CONTINUE);
	  }
	}
      }
/* 
 *   Duplicate gas phase, only if not solution equilibria
 */
    }
    if (gas_phase[i].solution_equilibria == FALSE)
    {
      n_user = gas_phase[i].n_user;
      last = gas_phase[i].n_user_end;
      gas_phase[i].n_user_end = gas_phase[i].n_user;
      for (j = n_user + 1; j <= last; j++)
      {
	gas_phase_duplicate (n_user, j);
      }
    }
    else
    {
      gas_phase[i].new_def = TRUE;
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
tidy_inverse (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   After all of data are read, fill in data for an inverse structure,
 *   including master species pointers, phase pointers, and uncertainties
 *   and a list of all elements from phases or -balance input.
 */
  int i, j, k, l;
  int count_in;
  LDBLE value;
  struct inv_elts *inv_elts;
  struct master *master_ptr;
  struct master *master_alk_ptr;
  struct elt_list *elt_list_ptr;
  master_alk_ptr = master_bsearch ("Alkalinity");
  for (i = 0; i < count_inverse; i++)
  {
    if (inverse[i].new_def != TRUE)
      continue;
/*
 *   Set default uncertainties for all solutions, if necessary
 */
    if (inverse[i].count_uncertainties < inverse[i].count_solns)
    {
      inverse[i].uncertainties =
	(LDBLE *) PHRQ_realloc (inverse[i].uncertainties,
				(size_t) inverse[i].count_solns *
				sizeof (LDBLE));
      if (inverse[i].uncertainties == NULL)
	malloc_error ();
      for (j = inverse[i].count_uncertainties; j < inverse[i].count_solns;
	   j++)
      {
	inverse[i].uncertainties[j] =
	  inverse[i].uncertainties[inverse[i].count_uncertainties - 1];
      }
    }
/*
 *   Set default ph uncertainties for all solutions, if necessary
 */
    if (inverse[i].count_ph_uncertainties < inverse[i].count_solns)
    {
      inverse[i].ph_uncertainties =
	(LDBLE *) PHRQ_realloc (inverse[i].ph_uncertainties,
				(size_t) inverse[i].count_solns *
				sizeof (LDBLE));
      if (inverse[i].ph_uncertainties == NULL)
	malloc_error ();
      for (j = inverse[i].count_ph_uncertainties; j < inverse[i].count_solns;
	   j++)
      {
	inverse[i].ph_uncertainties[j] =
	  inverse[i].ph_uncertainties[inverse[i].count_ph_uncertainties - 1];
      }
    }
/*
 *   Set default force for all solutions
 */
    if (inverse[i].count_force_solns < inverse[i].count_solns)
    {
      inverse[i].force_solns =
	(int *) PHRQ_realloc (inverse[i].force_solns,
			      (size_t) inverse[i].count_solns * sizeof (int));
      if (inverse[i].force_solns == NULL)
	malloc_error ();
      for (j = inverse[i].count_force_solns; j < inverse[i].count_solns; j++)
      {
	inverse[i].force_solns[j] = FALSE;
      }
    }
/*
 *   Find master species for element, set uncertainties
 */
    for (j = 0; j < inverse[i].count_elts; j++)
    {
      inverse[i].elts[j].master =
	master_bsearch_primary (inverse[i].elts[j].name);
      if (inverse[i].elts[j].master == NULL)
      {
	input_error++;
	sprintf (error_string, "No master species for element, %s.",
		 inverse[i].elts[j].name);
	error_msg (error_string, CONTINUE);
	continue;
      }
      inverse[i].elts[j].uncertainties =
	(double *) PHRQ_realloc (inverse[i].elts[j].uncertainties,
				 (size_t) inverse[i].count_solns *
				 sizeof (LDBLE));
      if (inverse[i].elts[j].uncertainties == NULL)
	malloc_error ();
      if (inverse[i].elts[j].count_uncertainties == 0)
      {
/* use default uncertainties for element */
	for (k = 0; k < inverse[i].count_solns; k++)
	{
	  inverse[i].elts[j].uncertainties[k] = inverse[i].uncertainties[k];
	}
      }
      else if (inverse[i].elts[j].count_uncertainties <
	       inverse[i].count_solns)
      {
/* use input uncertainties, fill in any missing at end */
	value =
	  inverse[i].elts[j].uncertainties[inverse[i].elts[j].
					   count_uncertainties - 1];
	for (k = inverse[i].elts[j].count_uncertainties;
	     k < inverse[i].count_solns; k++)
	{
	  inverse[i].elts[j].uncertainties[k] = value;
	}
      }
    }
/*
 *   Find phase
 */
    count_elts = 0;
    paren_count = 0;
    for (j = 0; j < inverse[i].count_phases; j++)
    {
      inverse[i].phases[j].phase =
	phase_bsearch (inverse[i].phases[j].name, &k, FALSE);
      if (inverse[i].phases[j].phase == NULL)
      {
	input_error++;
	sprintf (error_string, "Could not find phase, %s.",
		 inverse[i].phases[j].name);
	error_msg (error_string, CONTINUE);
	continue;
      }
/*
 *   Find isotope elements
 */
      if (inverse[i].phases[j].count_isotopes > 0)
      {
	for (k = 0; k < inverse[i].phases[j].count_isotopes; k++)
	{
	  inverse[i].phases[j].isotopes[k].primary = NULL;
	  inverse[i].phases[j].isotopes[k].master = NULL;
	  master_ptr =
	    master_bsearch (inverse[i].phases[j].isotopes[k].elt_name);
	  if (master_ptr == NULL)
	  {
	    input_error++;
	    sprintf (error_string,
		     "Element not found for isotope calculation: %s.",
		     inverse[i].phases[j].isotopes[k].elt_name);
	    error_msg (error_string, CONTINUE);
	    continue;
	  }
	  if (master_ptr->primary != TRUE)
	  {
	    input_error++;
	    sprintf (error_string, "Isotope ratio may only be used"
		     " for total element in phase.\n"
		     "Secondary species not allowed: %s.",
		     master_ptr->elt->name);
	    error_msg (error_string, CONTINUE);
	    continue;
	  }
	  inverse[i].phases[j].isotopes[k].primary = master_ptr;
	  inverse[i].phases[j].isotopes[k].master = master_ptr;
	  /* find coefficient for element */
	  for (elt_list_ptr = inverse[i].phases[j].phase->next_elt;
	       elt_list_ptr->elt != NULL; elt_list_ptr++)
	  {
	    if (elt_list_ptr->elt == master_ptr->elt)
	    {
	      inverse[i].phases[j].isotopes[k].coef = elt_list_ptr->coef;
	      break;
	    }
	  }
	  if (elt_list_ptr == NULL)
	  {
	    input_error++;
	    sprintf (error_string,
		     "Element, %s,for which isotope ratio was defined is not found in phase, %s",
		     master_ptr->elt->name, inverse[i].phases[j].phase->name);
	    error_msg (error_string, CONTINUE);
	    continue;
	  }
	}
	qsort (inverse[i].phases[j].isotopes,
	       (size_t) inverse[i].phases[j].count_isotopes,
	       (size_t) sizeof (struct isotope), isotope_compare);
      }
      add_elt_list (inverse[i].phases[j].phase->next_elt, 1.0);

    }
    if (input_error > 0)
      return (ERROR);
/*
 *   Sort elements in reaction and combine
 */
    if (count_elts > 0)
    {
      qsort (elt_list, (size_t) count_elts, (size_t) sizeof (struct elt_list),
	     elt_list_compare);
      elt_list_combine ();
    }
/*
 *   Mark master species list
 */
    for (j = 0; j < count_master; j++)
      master[j]->in = FALSE;
    for (j = 0; j < count_elts; j++)
    {
      elt_list[j].elt->master->in = TRUE;
    }
    /* Include all input elements */
    for (j = 0; j < inverse[i].count_elts; j++)
    {
      inverse[i].elts[j].master->in = TRUE;
    }
    s_eminus->primary->in = TRUE;	/* Include electrons */
    master_alk_ptr->in = TRUE;	/* Include alkalinity */
/*
 *   Unmark primary and mark secondary master species for redox elements
 */
    count_in = 0;
    inverse[i].count_redox_rxns = 0;
    for (j = 0; j < count_master; j++)
    {
      /*   skip all secondary master species in this loop */
      if (master[j]->primary == FALSE || master[j]->in == FALSE)
	continue;
      count_in++;
      if (j + 1 == count_master)
	continue;
      /*   if next master species is secondary, mark all 
         secondary master species until a primary is found */
      if (master[j + 1]->primary == FALSE)
      {
	master[j]->in = FALSE;
	count_in--;
	for (k = j + 1; k < count_master; k++)
	{
	  if (master[k]->primary == FALSE)
	  {
	    count_in++;
	    master[k]->in = TRUE;
	    if (master[k]->s->primary == NULL)
	    {
	      inverse[i].count_redox_rxns++;
	    }
	  }
	  else
	  {
	    break;
	  }
	}
      }
    }
/*
 *   Save list of master species in inv_elts structure
 */
    inv_elts =
      (struct inv_elts *) PHRQ_malloc ((size_t) (count_in) *
				       sizeof (struct inv_elts));
    if (inv_elts == NULL)
      malloc_error ();
    count_in = 0;
    for (j = 0; j < count_master; j++)
    {
      /* skip H(1) and O(-2) */
      if (master[j]->s == s_hplus || master[j]->s == s_h2o)
	continue;
      if (master[j]->in == TRUE)
      {
	/* set master */
	inv_elts[count_in].master = master[j];
	/* alloc uncertainties and set default */
	inv_elts[count_in].uncertainties =
	  (double *) PHRQ_malloc ((size_t) inverse[i].count_solns *
				  sizeof (LDBLE));
	if (inv_elts[count_in].uncertainties == NULL)
	  malloc_error ();
	for (k = 0; k < inverse[i].count_solns; k++)
	{
	  inv_elts[count_in].uncertainties[k] = inverse[i].uncertainties[k];
	}
	count_in++;
      }
    }
    if (s_co3->secondary->in == TRUE)
    {
      inverse[i].carbon = TRUE;
    }
    else
    {
      inverse[i].carbon = FALSE;
    }
/*
 *   copy in input uncertainties 
 */
    /* copy primary redox to all secondary redox */
    for (j = 0; j < inverse[i].count_elts; j++)
    {
      master_ptr = master_bsearch (inverse[i].elts[j].name);
      if (master_ptr == NULL)
      {
	input_error++;
	sprintf (error_string, "Element not found, %s.",
		 inverse[i].elts[j].name);
	error_msg (error_string, CONTINUE);
	continue;
      }
      if (master_ptr->primary == FALSE || master_ptr->s->secondary == NULL)
	continue;
      for (k = 0; k < count_in; k++)
      {
	if (master_ptr == inv_elts[k].master->elt->primary)
	{
	  for (l = 0; l < inverse[i].count_solns; l++)
	  {
	    inv_elts[k].uncertainties[l] =
	      inverse[i].elts[j].uncertainties[l];
	  }
	}
      }
      inverse[i].elts[j].uncertainties =
	(double *) free_check_null (inverse[i].elts[j].uncertainties);
    }
    /* copy masters that are not primary redox */
    for (j = 0; j < inverse[i].count_elts; j++)
    {
      master_ptr = master_bsearch (inverse[i].elts[j].name);
      if (master_ptr == NULL)
      {
	input_error++;
	sprintf (error_string, "Element not found, %s.",
		 inverse[i].elts[j].name);
	error_msg (error_string, CONTINUE);
	continue;
      }
      if (master_ptr->primary == TRUE && master_ptr->s->secondary != NULL)
	continue;
      for (k = 0; k < count_in; k++)
      {
	if (master_ptr == inv_elts[k].master)
	{
	  for (l = 0; l < inverse[i].count_solns; l++)
	  {
	    inv_elts[k].uncertainties[l] =
	      inverse[i].elts[j].uncertainties[l];
	  }
	  break;
	}
      }
      inverse[i].elts[j].uncertainties =
	(double *) free_check_null (inverse[i].elts[j].uncertainties);
    }
/*
 *   replace elts in inverse struct
 */
    inverse[i].elts = (struct inv_elts *) free_check_null (inverse[i].elts);
    inverse[i].elts = inv_elts;
    inverse[i].count_elts = count_in;
    for (j = 0; j < inverse[i].count_elts; j++)
    {
#ifdef SKIP
/*  make another pointer alkalinity uncertainties */
      if (inverse[i].elts[j].master == master_alk_ptr)
      {
	inverse[i].alk_uncertainties =
	  free_check_null (inverse[i].alk_uncertainties);
	inverse[i].alk_uncertainties = inverse[i].elts[j].uncertainties;
      }
#endif
/* debug
			output_msg(OUTPUT_MESSAGE, "\t%d\t%s", j, inverse[i].elts[j].master->elt->name);
			for (k = 0; k < inverse[i].count_solns; k++) {
				output_msg(OUTPUT_MESSAGE, "\t%f", inverse[i].elts[j].uncertainties[k]);
			}
			output_msg(OUTPUT_MESSAGE,"\n");
 */
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
tidy_phases (void)
/* ---------------------------------------------------------------------- */
{
  int i;
  int replaced;
  /*
   *  Fix log Ks first, so they can possibly be added to other phase equations
   */
  for (i = 0; i < count_phases; i++)
  {
    select_log_k_expression (phases[i]->logk, phases[i]->rxn->logk);
    add_other_logk (phases[i]->rxn->logk, phases[i]->count_add_logk,
		    phases[i]->add_logk);
    phases[i]->rxn->token[0].name = phases[i]->name;
    phases[i]->rxn->token[0].s = NULL;
  }
  /*
   *   Rewrite all phases to secondary species
   */
  for (i = 0; i < count_phases; i++)
  {
    /*
     *   Rewrite equation
     */
    count_trxn = 0;
    trxn_add_phase (phases[i]->rxn, 1.0, FALSE);
    trxn.token[0].name = phases[i]->name;
    /* debug 
       output_msg(OUTPUT_MESSAGE, "%s PHASE.\n", phases[i]->name);
       trxn_print();
     */
    replaced = replace_solids_gases ();
    /*  save rxn */
    rxn_free (phases[i]->rxn);
    phases[i]->rxn = rxn_alloc (count_trxn + 1);
    trxn_copy (phases[i]->rxn);
    /*  save rxn_s */
    trxn_reverse_k ();
    rewrite_eqn_to_secondary ();
    trxn_reverse_k ();
    rxn_free (phases[i]->rxn_s);
    phases[i]->rxn_s = rxn_alloc (count_trxn + 1);
    trxn_copy (phases[i]->rxn_s);
    /*
     *   Check equation
     */
    if (phases[i]->check_equation == TRUE)
    {
      if (replaced == FALSE)
      {
	phase_rxn_to_trxn (phases[i], phases[i]->rxn);
      }
      else
      {
	phase_rxn_to_trxn (phases[i], phases[i]->rxn_s);
      }
      if (check_eqn (FALSE) == ERROR)
      {
	input_error++;
	sprintf (error_string, "Equation for phase %s does not balance.",
		 phases[i]->name);
	error_msg (error_string, CONTINUE);
      }
    }
  }

  return (OK);
}

/* ---------------------------------------------------------------------- */
int
tidy_pp_assemblage (void)
/* ---------------------------------------------------------------------- */
{
  int i, j, k, l, n_user, first, last;
  struct phase *phase_ptr;
  LDBLE coef;
  char *ptr;
/*
 *   Find pointers for pure phases
 */
  for (i = 0; i < count_pp_assemblage; i++)
  {
    count_elts = 0;
    paren_count = 0;
    coef = 1.0;
    pp_assemblage[i].new_def = FALSE;
    for (j = 0; j < pp_assemblage[i].count_comps; j++)
    {
      phase_ptr =
	phase_bsearch (pp_assemblage[i].pure_phases[j].name, &k, FALSE);
      if (phase_ptr == NULL)
      {
	input_error++;
	sprintf (error_string, "Phase not found in data base, %s.",
		 pp_assemblage[i].pure_phases[j].name);
	error_msg (error_string, CONTINUE);
	continue;
      }
      else
      {
	pp_assemblage[i].pure_phases[j].phase = phase_ptr;
	add_elt_list (phase_ptr->next_elt, coef);

      }
      if (pp_assemblage[i].pure_phases[j].add_formula != NULL)
      {
	first = count_elts;
	phase_ptr =
	  phase_bsearch (pp_assemblage[i].pure_phases[j].add_formula, &k,
			 FALSE);
	if (phase_ptr != NULL)
	{
	  pp_assemblage[i].pure_phases[j].add_formula = phase_ptr->formula;
	}
	ptr = pp_assemblage[i].pure_phases[j].add_formula;
	get_elts_in_species (&ptr, coef);
	/* check that all elements are in the database */
	for (l = first; l < count_elts; l++)
	{
	  if (elt_list[l].elt->master == NULL)
	  {
	    input_error++;
	    sprintf (error_string,
		     "Element \"%s\" in alternative phase for \"%s\" in EQUILIBRIUM_PHASES not found in database.",
		     elt_list[l].elt->name,
		     pp_assemblage[i].pure_phases[j].name);
	    error_msg (error_string, CONTINUE);
	  }
	}
      }
    }
    if (count_elts > 0)
    {
      qsort (elt_list, (size_t) count_elts, (size_t) sizeof (struct elt_list),
	     elt_list_compare);
      elt_list_combine ();
    }
    pp_assemblage[i].next_elt =
      (struct elt_list *) free_check_null (pp_assemblage[i].next_elt);
    pp_assemblage[i].next_elt = elt_list_save ();

/*
 *   Store list with all elements in phases and add formulae
 */

/*
 *   Duplicate pure phases if necessary
 */
    if (pp_assemblage[i].n_user_end > pp_assemblage[i].n_user)
    {
      n_user = pp_assemblage[i].n_user;
      last = pp_assemblage[i].n_user_end;
      pp_assemblage[i].n_user_end = pp_assemblage[i].n_user;
      for (j = n_user + 1; j <= last; j++)
      {
	pp_assemblage_duplicate (n_user, j);
      }
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
tidy_s_s_assemblage (void)
/* ---------------------------------------------------------------------- */
{
  int i, j, k, k1, n_user, last;
  struct phase *phase_ptr;
  struct s_s *s_s_ptr;
  LDBLE nb, nc, n_tot, xb, xc, dnb, dnc, a0, a1;
  LDBLE xb2, xb3, xb4, xc2, xc3;
  LDBLE moles;
/*
 *   Find pointers for pure phases
 */
  for (i = 0; i < count_s_s_assemblage; i++)
  {
    count_elts = 0;
    paren_count = 0;
    for (j = 0; j < s_s_assemblage[i].count_s_s; j++)
    {
      for (k = 0; k < s_s_assemblage[i].s_s[j].count_comps; k++)
      {
	phase_ptr =
	  phase_bsearch (s_s_assemblage[i].s_s[j].comps[k].name, &k1, FALSE);
	if (phase_ptr == NULL)
	{
	  input_error++;
	  sprintf (error_string,
		   "Phase not found in data base, %s, assemblage %d.",
		   s_s_assemblage[i].s_s[j].comps[k].name,
		   s_s_assemblage[i].n_user);
	  error_msg (error_string, CONTINUE);
	  s_s_assemblage[i].s_s[j].comps[k].phase = NULL;
	  continue;
	}
	else
	{
	  s_s_assemblage[i].s_s[j].comps[k].phase = phase_ptr;
	  s_s_assemblage[i].s_s[j].comps[k].phase->moles_x = 0;
	  s_s_assemblage[i].s_s[j].comps[k].phase->fraction_x = 0;
	}
	if (s_s_assemblage[i].s_s[j].comps[k].moles == NAN)
	{
	  input_error++;
	  sprintf (error_string,
		   "Moles of solid solution component not defined, %s, assemblage %d.",
		   s_s_assemblage[i].s_s[j].comps[k].name,
		   s_s_assemblage[i].n_user);
	  error_msg (error_string, CONTINUE);
	  continue;
	}
      }

      if (s_s_assemblage[i].new_def == TRUE)
      {
	/*
	 *  Calculate a0 and a1 first
	 */
	s_s_calc_a0_a1 (&(s_s_assemblage[i].s_s[j]));
	s_s_ptr = &(s_s_assemblage[i].s_s[j]);

	n_tot = 0;
	for (k = 0; k < s_s_ptr->count_comps; k++)
	{
	  moles = s_s_assemblage[i].s_s[j].comps[k].moles;
	  if (s_s_assemblage[i].s_s[j].comps[k].moles <= 0.0)
	  {
	    moles = MIN_TOTAL_SS;
	    s_s_assemblage[i].s_s[j].comps[k].initial_moles = moles;
	  }
	  n_tot += moles;
	}

	for (k = 0; k < s_s_ptr->count_comps; k++)
	{
	  moles = s_s_assemblage[i].s_s[j].comps[k].moles;
	  if (s_s_assemblage[i].s_s[j].comps[k].moles <= 0.0)
	  {
	    moles = MIN_TOTAL_SS;
	  }
	  s_s_assemblage[i].s_s[j].comps[k].fraction_x = moles / n_tot;
	  s_s_assemblage[i].s_s[j].comps[k].log10_fraction_x =
	    log10 (moles / n_tot);
	}
	a0 = s_s_assemblage[i].s_s[j].a0;
	a1 = s_s_assemblage[i].s_s[j].a1;

/*
 *   Binary solid solution
 */
	if (a0 != 0.0 || a1 != 0)
	{
	  s_s_assemblage[i].s_s[j].dn = 1.0 / n_tot;
	  nc = s_s_assemblage[i].s_s[j].comps[0].moles;
	  if (nc == 0)
	    nc = MIN_TOTAL_SS;
	  nb = s_s_assemblage[i].s_s[j].comps[1].moles;
	  if (nb == 0)
	    nb = MIN_TOTAL_SS;
	  xc = nc / n_tot;
	  xb = nb / n_tot;

	  /* lambdas */
	  s_s_assemblage[i].s_s[j].comps[0].log10_lambda =
	    xb * xb * (a0 - a1 * (3 - 4 * xb)) / LOG_10;
	  s_s_assemblage[i].s_s[j].comps[1].log10_lambda =
	    xc * xc * (a0 + a1 * (4 * xb - 1)) / LOG_10;

	  /* derivatives wrt nc and nb */
	  xc2 = xc * xc;
	  xc3 = xc2 * xc;
	  xb2 = xb * xb;
	  xb3 = xb2 * xb;
	  xb4 = xb3 * xb;

	  /* component 1 */
	  dnb =
	    -2 * a0 * xb * xc2 - 8 * a1 * xb2 * xc2 + 6 * a1 * xb * xc2 -
	    4 * a1 * xc * xb4 - 8 * a1 * xb3 * xc2 - 4 * a1 * xb2 * xc3 -
	    2 * a0 * xc * xb2 - 8 * a1 * xc * xb3 + 6 * a1 * xc * xb2 + 1;
	  s_s_assemblage[i].s_s[j].comps[0].dnb = dnb / n_tot;
	  dnc =
	    2 * a0 * xb3 + 2 * a0 * xc * xb2 + 8 * a1 * xb4 +
	    8 * a1 * xc * xb3 - 2 * a1 * xb3 - 6 * a1 * xc * xb2;
	  s_s_assemblage[i].s_s[j].comps[0].dnc = -xb / nc + dnc / n_tot;
	  s_s_assemblage[i].s_s[j].comps[0].dn = 1.0 / n_tot;

	  /* component 2 */
	  dnb =
	    2 * a0 * xb * xc2 + 2 * a0 * xc3 + 8 * a1 * xb2 * xc2 +
	    8 * a1 * xb * xc3 - 2 * a1 * xb * xc2 - 6 * a1 * xc3;
	  s_s_assemblage[i].s_s[j].comps[1].dnb = -xc / nb + dnb / n_tot;
	  dnc =
	    -2 * a0 * xc * xb2 - 8 * a1 * xc * xb3 + 2 * a1 * xc * xb2 -
	    2 * a0 * xb * xc2 - 8 * a1 * xb2 * xc2 + 6 * a1 * xb * xc2 + 1;
	  s_s_assemblage[i].s_s[j].comps[1].dnc = dnc / n_tot;
	  s_s_prep (s_s_assemblage[i].s_s[j].tk, &(s_s_assemblage[i].s_s[j]),
		    TRUE);
	  s_s_assemblage[i].s_s[j].comps[1].dn = 1.0 / n_tot;
/*
 *   Ideal solid solution
 */
	}
	else
	{
	  s_s_assemblage[i].s_s[j].dn = 1.0 / n_tot;
	  for (k = 0; k < s_s_ptr->count_comps; k++)
	  {
	    s_s_assemblage[i].s_s[j].comps[k].log10_lambda = 0;
	    moles = s_s_assemblage[i].s_s[j].comps[k].moles;
	    if (moles <= 0.0)
	      moles = MIN_TOTAL_SS;
	    s_s_assemblage[i].s_s[j].comps[k].dnb =
	      (n_tot - moles) / (moles * n_tot);
	    s_s_assemblage[i].s_s[j].comps[k].dn = 1.0 / n_tot;
	  }
	}
      }
    }
    s_s_assemblage[i].new_def = FALSE;

/*
 *   Duplicate s_s_assemblage if necessary
 */
    n_user = s_s_assemblage[i].n_user;
    last = s_s_assemblage[i].n_user_end;
    s_s_assemblage[i].n_user_end = s_s_assemblage[i].n_user;
    for (j = n_user + 1; j <= last; j++)
    {
      s_s_assemblage_duplicate (n_user, j);
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
tidy_punch (void)
/* ---------------------------------------------------------------------- */
{
  int i, j, l;
  int punch_save;
  char token[MAX_LENGTH];
/*
 *   tidy punch information
 */
  if (punch.high_precision == FALSE)
  {
    l = 12;
  }
  else
  {
    l = 20;
  }
  if (punch.in == TRUE)
  {
    /* totals */

    for (i = 0; i < punch.count_totals; i++)
    {
      punch.totals[i].master = master_bsearch (punch.totals[i].name);
    }

    /* molalities */

    for (i = 0; i < punch.count_molalities; i++)
    {
      punch.molalities[i].s = s_search (punch.molalities[i].name);
    }

    /* log activities */

    for (i = 0; i < punch.count_activities; i++)
    {
      punch.activities[i].s = s_search (punch.activities[i].name);
    }

    /* equilibrium phases */

    for (i = 0; i < punch.count_pure_phases; i++)
    {
      punch.pure_phases[i].phase =
	phase_bsearch (punch.pure_phases[i].name, &j, FALSE);
    }

    /* saturation indices */

    for (i = 0; i < punch.count_si; i++)
    {
      punch.si[i].phase = phase_bsearch (punch.si[i].name, &j, FALSE);
    }

    /* gases */

    for (i = 0; i < punch.count_gases; i++)
    {
      punch.gases[i].phase = phase_bsearch (punch.gases[i].name, &j, FALSE);
    }
  }
  /*
   *  Always write new headings when SELECTED_OUTPUT is read
   */
  if (punch.new_def == TRUE && punch.in == TRUE)
  {
    punch_save = pr.punch;
    pr.punch = TRUE;

    /* constant stuff, sim, pH, etc. */

    if (punch.sim == TRUE)
    {
      output_msg (OUTPUT_PUNCH, "%*s\t", l, "sim");
    }
    if (punch.state == TRUE)
    {
      output_msg (OUTPUT_PUNCH, "%*s\t", l, "state");
    }
    if (punch.soln == TRUE)
    {
      output_msg (OUTPUT_PUNCH, "%*s\t", l, "soln");
    }
    if (punch.dist == TRUE)
    {
      output_msg (OUTPUT_PUNCH, "%*s\t", l, "dist_x");
    }
    if (punch.time == TRUE)
    {
      output_msg (OUTPUT_PUNCH, "%*s\t", l, "time");
    }
    if (punch.step == TRUE)
    {
      output_msg (OUTPUT_PUNCH, "%*s\t", l, "step");
    }
    if (punch.ph == TRUE)
    {
      output_msg (OUTPUT_PUNCH, "%*s\t", l, "pH");
    }
    if (punch.pe == TRUE)
    {
      output_msg (OUTPUT_PUNCH, "%*s\t", l, "pe");
    }
    if (punch.rxn == TRUE)
    {
      output_msg (OUTPUT_PUNCH, "%*s\t", l, "reaction");
    }
    if (punch.temp == TRUE)
    {
      output_msg (OUTPUT_PUNCH, "%*s\t", l, "temp");
    }
    if (punch.alk == TRUE)
    {
      output_msg (OUTPUT_PUNCH, "%*s\t", l, "Alk");
    }
    if (punch.mu == TRUE)
    {
      output_msg (OUTPUT_PUNCH, "%*s\t", l, "mu");
    }
    if (punch.water == TRUE)
    {
      output_msg (OUTPUT_PUNCH, "%*s\t", l, "mass_H2O");
    }
    if (punch.charge_balance == TRUE)
    {
      output_msg (OUTPUT_PUNCH, "%*s\t", l, "charge");
    }
    if (punch.percent_error == TRUE)
    {
      output_msg (OUTPUT_PUNCH, "%*s\t", l, "pct_err");
    }
    /* totals */

    for (i = 0; i < punch.count_totals; i++)
    {
      output_msg (OUTPUT_PUNCH, "%*s\t", l, punch.totals[i].name);
      if (punch.totals[i].master == NULL)
      {
	sprintf (error_string, "Did not find master species,"
		 " %s.", punch.totals[i].name);
	warning_msg (error_string);
      }
    }

    /* molalities */

    for (i = 0; i < punch.count_molalities; i++)
    {
      strcpy (token, "m_");
      strcat (token, punch.molalities[i].name);
      output_msg (OUTPUT_PUNCH, "%*s\t", l, token);
      if (punch.molalities[i].s == NULL)
      {
	sprintf (error_string, "Did not find species,"
		 " %s.", punch.molalities[i].name);
	warning_msg (error_string);
      }
    }

    /* log activities */

    for (i = 0; i < punch.count_activities; i++)
    {
      strcpy (token, "la_");
      strcat (token, punch.activities[i].name);
      output_msg (OUTPUT_PUNCH, "%*s\t", l, token);
      if (punch.activities[i].s == NULL)
      {
	sprintf (error_string, "Did not find species, "
		 "%s.", punch.activities[i].name);
	warning_msg (error_string);
      }
    }

    /* equilibrium phases */

    for (i = 0; i < punch.count_pure_phases; i++)
    {
      strcpy (token, "d_");
      strcat (token, punch.pure_phases[i].name);
      output_msg (OUTPUT_PUNCH, "%*s\t%*s\t", l, punch.pure_phases[i].name, l,
		  token);
      if (punch.pure_phases[i].phase == NULL)
      {
	sprintf (error_string, "Did not find phase, "
		 "%s.", punch.pure_phases[i].name);
	warning_msg (error_string);
      }
    }

    /* saturation indices */

    for (i = 0; i < punch.count_si; i++)
    {
      strcpy (token, "si_");
      strcat (token, punch.si[i].name);
      output_msg (OUTPUT_PUNCH, "%*s\t", l, token);
      if (punch.si[i].phase == NULL)
      {
	sprintf (error_string, "Did not find phase, "
		 "%s.", punch.si[i].name);
	warning_msg (error_string);
      }
    }

    /* gases */

    if (punch.count_gases > 0)
    {
      output_msg (OUTPUT_PUNCH, "%*s\t%*s\t%*s\t", l, "pressure", l,
		  "total mol", l, "volume");
    }
    for (i = 0; i < punch.count_gases; i++)
    {
      strcpy (token, "g_");
      strcat (token, punch.gases[i].name);
      output_msg (OUTPUT_PUNCH, "%*s\t", l, token);
      if (punch.gases[i].phase == NULL)
      {
	sprintf (error_string, "Did not find phase, "
		 "%s.", punch.gases[i].name);
	warning_msg (error_string);
      }
    }

    /* kinetics */

    for (i = 0; i < punch.count_kinetics; i++)
    {
      strcpy (token, "k_");
      strcat (token, punch.kinetics[i].name);
      output_msg (OUTPUT_PUNCH, "%*s\t", l, token);
      strcpy (token, "dk_");
      strcat (token, punch.kinetics[i].name);
      output_msg (OUTPUT_PUNCH, "%*s\t", l, token);
    }

    /* solid solutions */

    for (i = 0; i < punch.count_s_s; i++)
    {
      strcpy (token, "s_");
      strcat (token, punch.s_s[i].name);
      output_msg (OUTPUT_PUNCH, "%*s\t", l, token);
    }

    /* isotopes */

    for (i = 0; i < punch.count_isotopes; i++)
    {
      if (isotope_ratio_search (punch.isotopes[i].name) == NULL)
      {
	sprintf (error_string, "Did not find isotope_ratio definition for "
		 "%s in -isotopes of SELECTED_OUTPUT.\n%s must be defined in ISOTOPE_RATIO data block.",
		 punch.isotopes[i].name, punch.isotopes[i].name);
	warning_msg (error_string);
      }
      strcpy (token, "I_");
      strcat (token, punch.isotopes[i].name);
      output_msg (OUTPUT_PUNCH, "%*s\t", l, token);
    }

    /* calculate_values */

    for (i = 0; i < punch.count_calculate_values; i++)
    {
      if (calculate_value_search (punch.calculate_values[i].name) == NULL)
      {
	sprintf (error_string, "Did not find calculate_values definition for "
		 "%s in -calculate_values of SELECTED_OUTPUT.\n%s must be defined in CALCULATE_VALUES data block.",
		 punch.calculate_values[i].name,
		 punch.calculate_values[i].name);
	warning_msg (error_string);
      }
      strcpy (token, "V_");
      strcat (token, punch.calculate_values[i].name);
      output_msg (OUTPUT_PUNCH, "%*s\t", l, token);
    }

    /* user_punch */
    if (punch.user_punch == TRUE)
    {
      for (i = 0; i < user_punch_count_headings; i++)
      {
	output_msg (OUTPUT_PUNCH, "%*s\t", l, user_punch_headings[i]);
      }
    }
    output_msg (OUTPUT_PUNCH, "\n");
    punch.new_def = FALSE;
    pr.punch = punch_save;
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
tidy_species (void)
/* ---------------------------------------------------------------------- */
{
  int i, j;
  struct master *master_ptr;
  char c, *ptr;
/*
 *   Make sure species pointers are ok
 */
  if (check_species_input () == ERROR)
  {
    error_msg ("Calculations terminating due to input errors.", STOP);
  }
/*
 *   Set secondary and primary pointers in species structures
 */
  for (i = 0; i < count_s; i++)
  {
    s[i]->number = i;
    s[i]->primary = NULL;
    s[i]->secondary = NULL;
    if (s[i]->check_equation == TRUE)
    {
      species_rxn_to_trxn (s[i]);
      if (check_eqn (TRUE) == ERROR)
      {
	input_error++;
	sprintf (error_string, "Equation for species %s does not balance.",
		 s[i]->name);
	error_msg (error_string, CONTINUE);
      }
    }
  }
  for (i = 0; i < count_master; i++)
  {
    ptr = master[i]->elt->name;
    if (ptr[0] != '[')
    {
      while ((c = (int) *(++ptr)) != '\0')
      {
	if (isupper ((int) c))
	{
	  input_error++;
	  sprintf (error_string,
		   "Element or valence name in SOLUTION_MASTER_SPECIES should include only one element, %s.",
		   master[i]->elt->name);
	  error_msg (error_string, CONTINUE);
	  break;
	}
      }
    }
    /* store sequence number in master structure */
    master[i]->number = i;
    if (master[i]->primary == TRUE)
    {
      master[i]->s->primary = master[i];
    }
    else
    {
      master[i]->s->secondary = master[i];
    }
    if (strcmp (master[i]->elt->name, "C") == 0)
    {
      s_co3 = master[i]->s;
    }
    if (master[i]->gfw_formula != NULL)
    {
      if (compute_gfw (master[i]->gfw_formula, &master[i]->gfw) == ERROR)
      {
	input_error++;
	sprintf (error_string,
		 "Calculating gfw for master species, %s, formula %s.",
		 master[i]->elt->name, master[i]->gfw_formula);
	error_msg (error_string, CONTINUE);
      }
    }
  }
/*
 *   Write equations for all master species in terms of primary
 *   master species, set coefficient of element in master species
 */
  for (i = 0; i < count_master; i++)
  {
    count_trxn = 0;
    if (master[i]->s->primary != NULL)
    {
      trxn_add (master[i]->s->rxn, 1.0, FALSE);
      trxn_add (master[i]->s->rxn, -1.0, TRUE);
    }
    else
    {
      trxn_add (master[i]->s->rxn, 1.0, FALSE);
      rewrite_eqn_to_primary ();
    }
    rxn_free (master[i]->rxn_primary);
    master[i]->rxn_primary = rxn_alloc (count_trxn + 1);
    trxn_copy (master[i]->rxn_primary);
    master[i]->coef = coef_in_master (master[i]);
  }
/*
 *   Rewrite all species to secondary species
 */
  for (i = 0; i < count_s; i++)
  {
    count_trxn = 0;
    if (s[i]->primary != NULL || s[i]->secondary != NULL)
    {
      trxn_add (s[i]->rxn, 1.0, FALSE);
      trxn_add (s[i]->rxn, -1.0, TRUE);
    }
    else
    {
      trxn_add (s[i]->rxn, 1.0, FALSE);
      rewrite_eqn_to_secondary ();
    }
    rxn_free (s[i]->rxn_s);
    s[i]->rxn_s = rxn_alloc (count_trxn + 1);
    trxn_copy (s[i]->rxn_s);
    /* calculate alkalinity */
    s[i]->alk = calc_alk (s[i]->rxn_s);
    /* set co2 coefficient */
    s[i]->co2 = 0.0;
    for (j = 1; j < count_trxn; j++)
    {
      if (trxn.token[j].s == s_co3)
      {
	s[i]->co2 = trxn.token[j].coef;
	break;
      }
    }
  }
/*
 *   Set pointer in element to master species
 */
  for (i = 0; i < count_elements; i++)
  {
    elements[i]->master = master_bsearch (elements[i]->name);
    if (elements[i]->master == NULL)
    {
      input_error++;
      sprintf (error_string, "No master species for element %s.",
	       elements[i]->name);
      error_msg (error_string, CONTINUE);
    }
    elements[i]->primary = master_bsearch_primary (elements[i]->name);
    if (elements[i]->primary == NULL)
    {
      input_error++;
      sprintf (error_string, "No master species for element %s.",
	       elements[i]->name);
      error_msg (error_string, CONTINUE);
    }
  }
/*
 *   Make sure all primary master species for redox elements
 *   are also secondary master species
 */
  for (i = 0; i < count_master; i++)
  {
    if (master[i]->primary == FALSE)
    {
      master_ptr = master[i]->s->secondary->elt->primary;
      if (master_ptr == NULL)
      {
	input_error++;
	sprintf (error_string,
		 "Every primary master species for a redox element\n"
		 "\tmust also be a secondary master species.\n"
		 "\tError in definitions related to %s .\n",
		 master[i]->s->name);
	error_msg (error_string, CONTINUE);

      }
      else if (master_ptr->s->secondary == NULL)
      {
	input_error++;
	sprintf (error_string,
		 "Every primary master species for a redox element\n"
		 "\tmust also be a secondary master species.\n"
		 "\t%s is the primary master species for element %s.\n"
		 "\tAnother entry in SOLUTION_MASTER_SPECIES is needed.\n"
		 "\tDefine species %s as a secondary master species for a valence state.\n"
		 "\tFor example: \n" "\t%s(0)\t%s alk gfw",
		 master_ptr->s->name, master_ptr->elt->name,
		 master_ptr->s->name, master_ptr->elt->name,
		 master_ptr->s->name);
	error_msg (error_string, CONTINUE);
      }
    }
  }
/*
 *   Calculate H and O if alternate mass balance is given
 */
  for (i = 0; i < count_s; i++)
  {
    if (s[i]->next_secondary != NULL)
    {
      s[i]->h = 0.0;
      s[i]->o = 0.0;
      for (j = 0; s[i]->next_secondary[j].elt != NULL; j++)
      {
	if (s[i]->next_secondary[j].elt->primary == NULL)
	  continue;
	if (s[i]->next_secondary[j].elt->primary->s == s_hplus)
	{
	  s[i]->h += s[i]->next_secondary[j].coef;
	}
	else if (s[i]->next_secondary[j].elt->primary->s == s_h2o)
	{
	  s[i]->o += s[i]->next_secondary[j].coef;
	}
	else if (s[i]->mole_balance != NULL)
	{
	  master_ptr = s[i]->next_secondary[j].elt->master;
	  if (master_ptr->primary == TRUE)
	  {
	    if (master_ptr->s->secondary != NULL)
	    {
	      master_ptr = master_ptr->s->secondary;
	    }
	  }
	  if (master_ptr->coef != 1)
	  {
	    s[i]->next_secondary[j].coef /= master_ptr->coef;
	  }
	}
      }
      if (s[i]->type == EX)
      {
	for (j = 0; s[i]->next_secondary[j].elt != NULL; j++)
	{
	  if (s[i]->next_secondary[j].elt->primary->s->type == EX)
	  {
	    s[i]->equiv = s[i]->next_secondary[j].coef;
	    break;
	  }
	}
      }
    }
  }
#ifdef SKIP
  for (i = 0; i < count_s; i++)
  {
    if (match_elts_in_species (s[i]->name, "*{C,[13C]}{O,[18O]}3*") == TRUE)
    {
      output_msg (OUTPUT_MESSAGE, "Match: %s\n", s[i]->name);
    }
  }
#endif
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
tidy_surface (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   After all of data are read, fill in master species for surface comps
 *   Sort surface
 */
  int i, j, k;
  char *ptr1;
  struct surface *surface_ptr;
  struct master *master_ptr;

  for (k = 0; k < count_surface; k++)
  {
    surface_ptr = &surface[k];

    for (i = 0; i < surface_ptr->count_comps; i++)
    {
/*
 *   Find master species for each surface
 */
      for (j = 0; surface_ptr->comps[i].totals[j].elt != NULL; j++)
      {
	master_ptr = surface_ptr->comps[i].totals[j].elt->master;
	if (master_ptr == NULL)
	{
	  input_error++;
	  sprintf (error_string, "Master species not in data base for %s, "
		   "skipping element.",
		   surface_ptr->comps[i].totals[j].elt->name);
	  error_msg (error_string, CONTINUE);
	  continue;
	}
	if (master_ptr->type != SURF)
	  continue;
/*
 *   Set flags
 */
	surface_ptr->comps[i].master = master_ptr;
	/*
	 * Calculate moles of sites
	 */
	if (surface_ptr->new_def == TRUE
	    && surface_ptr->sites_units == SITES_DENSITY
	    && surface_ptr->comps[i].phase_name == NULL)
	{
	  surface_ptr->comps[i].moles =
	    surface_ptr->comps[i].moles * 1.0e18 *
	    surface_ptr->charge[surface_ptr->comps[i].charge].specific_area *
	    surface_ptr->charge[surface_ptr->comps[i].charge].grams /
	    AVOGADRO;
	  /*
	   *  Calculate totals
	   */
	  count_elts = 0;
	  paren_count = 0;
	  ptr1 = surface_ptr->comps[i].formula;
	  get_elts_in_species (&ptr1, surface_ptr->comps[i].moles);
	  surface_ptr->comps[i].formula_totals =
	    (struct elt_list *) free_check_null (surface_ptr->comps[i].
						 formula_totals);
	  surface_ptr->comps[i].formula_totals = elt_list_save ();
	  surface_ptr->comps[i].totals =
	    (struct elt_list *) free_check_null (surface_ptr->comps[i].
						 totals);
	  surface_ptr->comps[i].totals = elt_list_save ();
	}
	if (surface_ptr->type == CD_MUSIC)
	{
	  /*
	     surface_ptr->charge[surface_ptr->comps[i].charge].charge_balance += surface_ptr->comps[i].moles*surface_ptr->comps[i].master->s->z;
	   */
	  surface_ptr->charge[surface_ptr->comps[i].charge].charge_balance +=
	    surface_ptr->comps[i].moles * surface_ptr->comps[i].formula_z;

	}
	break;
      }
#ifdef SKIP_MUSIC
      /*
       *  If charge of formula is non zero
       */
      if (surface_ptr->type == CD_MUSIC)
      {
	surface_ptr->comps[i].cb =
	  surface_ptr->comps[i].formula_z * surface_ptr->comps[i].moles;
      }
#endif
    }
/*
 *   Sort components
 */
    if (input_error == 0)
    {
      qsort (surface[k].comps,
	     (size_t) surface_ptr->count_comps,
	     (size_t) sizeof (struct surface_comp), surface_comp_compare);
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
tidy_solutions (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Define n_user for any solutions read by solution_spread that
 *   don't have n_user defined
 */
  int i, l, n, last;
  struct conc *tot_ptr;
  char *ptr;
  struct master *master_ptr;
  char token[MAX_LENGTH];
  for (n = 0; n < count_solution; n++)
  {
    /*
     * Check that elements are in database
     */
    if (solution[n] != NULL && solution[n]->new_def == TRUE)
    {
      /*
       *   Sort totals by description
       */
      i = 0;
      while (solution[n]->totals[i].description != NULL)
	i++;
      qsort (solution[n]->totals,
	     (size_t) i, (size_t) sizeof (struct conc), conc_compare);
      /*
       *  sort isotopes
       */
      if (solution[n]->count_isotopes > 0)
      {
	qsort (solution[n]->isotopes,
	       (size_t) solution[n]->count_isotopes,
	       (size_t) sizeof (struct isotope), isotope_compare);
      }
      else
      {
	solution[n]->isotopes =
	  (struct isotope *) free_check_null (solution[n]->isotopes);
      }
      for (i = 0; solution[n]->totals[i].description != NULL; i++)
      {
	tot_ptr = &(solution[n]->totals[i]);
	if (strcmp (tot_ptr->description, "H(1)") == 0 ||
	    strcmp (tot_ptr->description, "E") == 0)
	{
	  tot_ptr->moles = 0.0;
	  continue;
	}
	ptr = tot_ptr->description;
	copy_token (token, &ptr, &l);
	master_ptr = master_bsearch (token);
	if (master_ptr == NULL)
	{
	  sprintf (error_string,
		   "Could not find element in database, %s.\n\tConcentration is set to zero.",
		   tot_ptr->description);
	  warning_msg (error_string);
	  tot_ptr->input_conc = 0.0;
	  continue;
	}
      }
    }
  }
  /*
   *  Calculate solution numbers
   */
  for (n = 0; n < count_solution; n++)
  {
    if (solution[n] != NULL && solution[n]->new_def == TRUE
	&& solution[n]->n_user < 0)
    {
      last = 0;
      for (i = 0; i < count_solution; i++)
      {
	if (solution[i]->n_user > last)
	  last = solution[i]->n_user;
	if (solution[i]->n_user_end > last)
	  last = solution[i]->n_user_end;
      }
      if (save.solution == TRUE)
      {
	if (save.n_solution_user > last)
	  last = save.n_solution_user;
	if (save.n_solution_user_end > last)
	  last = save.n_solution_user_end;
      }
      for (i = 0; i < count_solution; i++)
      {
	if (solution[i]->new_def == TRUE && solution[i]->n_user < 0)
	{
	  solution[i]->n_user = ++last;
	  solution[i]->n_user_end = last;
	  if (use.solution_in == TRUE && use.n_solution_user == -1)
	  {
	    use.n_solution_user = last;
	  }
	}
      }
      break;
    }
  }
  solution_sort ();
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
species_rxn_to_trxn (struct species *s_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Copy reaction from reaction structure to 
 *   temp reaction structure.
 */
  int i;

  for (i = 0; s_ptr->rxn->token[i].s != NULL; i++)
  {
    trxn.token[i].name = s_ptr->rxn->token[i].s->name;
    trxn.token[i].z = s_ptr->rxn->token[i].s->z;
    trxn.token[i].s = s_ptr->rxn->token[i].s;
    trxn.token[i].unknown = NULL;
    trxn.token[i].coef = s_ptr->rxn->token[i].coef;
    count_trxn = i + 1;
    if (count_trxn + 1 >= max_trxn)
    {
      space ((void **) ((void *) &(trxn.token)), count_trxn + 1, &max_trxn,
	     sizeof (struct rxn_token_temp));
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
phase_rxn_to_trxn (struct phase *phase_ptr, struct reaction *rxn_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Copy reaction from reaction structure to 
 *   temp reaction structure.
 */
  int i, l;
  char *ptr;
  char token[MAX_LENGTH];
  LDBLE z;

  trxn.token[0].name = phase_ptr->formula;
  /* charge */
  ptr = phase_ptr->formula;
  get_token (&ptr, token, &z, &l);
  trxn.token[0].z = z;
  trxn.token[0].s = NULL;
  trxn.token[0].unknown = NULL;
  /*trxn.token[0].coef = -1.0; */
  /* check for leading coefficient of 1.0 for phase did not work */
  trxn.token[0].coef = phase_ptr->rxn->token[0].coef;
  for (i = 1; rxn_ptr->token[i].s != NULL; i++)
  {
    trxn.token[i].name = rxn_ptr->token[i].s->name;
    trxn.token[i].z = rxn_ptr->token[i].s->z;
    trxn.token[i].s = NULL;
    trxn.token[i].unknown = NULL;
    trxn.token[i].coef = rxn_ptr->token[i].coef;
    count_trxn = i + 1;
    if (count_trxn + 1 >= max_trxn)
    {
      space ((void **) ((void *) &(trxn.token)), count_trxn + 1, &max_trxn,
	     sizeof (struct rxn_token_temp));
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
tidy_isotopes (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Isotope ratios for each element or element valence state
 */
  int i, j, k, l, n;
  int count_isotopes, primary_number, count_primary;
  LDBLE isotope_number;
  struct master *master_ptr, *primary_ptr;
  struct solution *solution_ptr;
  struct isotope *new_isotopes;
  struct isotope *primary_isotopes;
  char token[MAX_LENGTH];

  primary_number = 0;
  primary_ptr = NULL;
  for (n = 0; n < count_solution; n++)
  {
    if (solution[n]->new_def != TRUE)
      continue;
    if (solution[n]->count_isotopes == 0)
      continue;
    solution_ptr = solution[n];
    primary_isotopes =
      (struct isotope *) PHRQ_malloc ((size_t) (sizeof (struct isotope)));
    if (primary_isotopes == NULL)
      malloc_error ();
    count_primary = 0;
    new_isotopes =
      (struct isotope *) PHRQ_malloc ((size_t) (sizeof (struct isotope)));
    if (new_isotopes == NULL)
      malloc_error ();
    count_isotopes = 0;
/*
 *   Make list of primary master species for isotopes
 */
    for (i = 0; i < solution_ptr->count_isotopes; i++)
    {
      master_ptr =
	master_bsearch_primary (solution_ptr->isotopes[i].elt_name);
      isotope_number = solution_ptr->isotopes[i].isotope_number;
      if (master_ptr == NULL)
      {
	input_error++;
	sprintf (error_string,
		 "In isotope calculation: element not defined: %s.",
		 solution_ptr->isotopes[i].elt_name);
	error_msg (error_string, CONTINUE);
	continue;
      }
      for (j = 0; j < count_primary; j++)
      {
	if (primary_isotopes[j].master == master_ptr &&
	    primary_isotopes[j].isotope_number == isotope_number)
	  break;
      }
      if (j == count_primary)
      {
	primary_isotopes =
	  (struct isotope *) PHRQ_realloc (primary_isotopes,
					   (size_t) (count_primary +
						     1) *
					   sizeof (struct isotope));
	if (primary_isotopes == NULL)
	  malloc_error ();
	primary_isotopes[count_primary].master = master_ptr;
	primary_isotopes[count_primary].isotope_number = isotope_number;
	count_primary++;
      }
    }
    if (input_error > 0)
      return (ERROR);
/*
 *   Go through all redox states of the list of primary species and isotope number
 */
    for (j = 0; j < count_primary; j++)
    {

      /* find index number of master species, set flag to FALSE */
      master_ptr = primary_isotopes[j].master;
      isotope_number = primary_isotopes[j].isotope_number;
      for (k = 0; k < count_master; k++)
      {
	if (master[k] == master_ptr)
	{
	  primary_number = k;
	  primary_ptr = master[k];
	}
	master[k]->isotope = FALSE;
      }

      /* go through isotopes of solution and fill in master species */
      for (l = 0; l < solution_ptr->count_isotopes; l++)
      {
	master_ptr = master_bsearch (solution_ptr->isotopes[l].elt_name);
	if (master_ptr == NULL)
	{
	  input_error++;
	  sprintf (error_string,
		   "In isotope calculation: element not defined: %s.",
		   solution_ptr->isotopes[l].elt_name);
	  error_msg (error_string, CONTINUE);
	  continue;
	}

	/* only fill for pertinent isotope */
	if (master_ptr->elt->primary != primary_ptr)
	  continue;
	if (solution_ptr->isotopes[l].isotope_number != isotope_number)
	  continue;

	/* for primary, fill in ratio for all secondary species */
	if (master_ptr->primary == TRUE && master_ptr->s->secondary != NULL)
	{
	  for (k = primary_number + 1; k < count_master; k++)
	  {
	    if (master[k]->elt->primary != primary_ptr)
	      break;
	    master[k]->isotope_ratio = solution_ptr->isotopes[l].ratio;
	    master[k]->isotope_ratio_uncertainty =
	      solution_ptr->isotopes[l].ratio_uncertainty;
	    if (master[k]->isotope == TRUE)
	    {
	      sprintf (error_string,
		       "In isotope calculation: redefinition of isotope ratio for %s.",
		       solution_ptr->isotopes[l].elt_name);
	      error_msg (error_string, CONTINUE);
	    }
	    master[k]->isotope = TRUE;
	  }

	  /* for secondary and non redox, set ratio */
	}
	else
	{
	  master_ptr->isotope_ratio = solution_ptr->isotopes[l].ratio;
	  master_ptr->isotope_ratio_uncertainty =
	    solution_ptr->isotopes[l].ratio_uncertainty;
	  if (master_ptr->isotope == TRUE)
	  {
	    sprintf (error_string,
		     "In isotope calculation: redefinition of isotope ratio for %s.",
		     solution_ptr->isotopes[l].elt_name);
	    error_msg (error_string, CONTINUE);
	  }
	  master_ptr->isotope = TRUE;
	}
      }
/*
 *   Write new isotope structure
 */
      for (k = 0; k < count_master; k++)
      {
	/* skip primary master species of redox elements */
	if (master[k]->primary == TRUE && master[k]->s->secondary != NULL)
	  continue;
	if (master[k]->elt->primary == primary_ptr
	    && master[k]->isotope == FALSE)
	{
	  input_error++;
	  sprintf (error_string,
		   "Isotopic ratio not defined for element or valence state %g%s, using 0.",
		   (double) isotope_number, master[k]->elt->name);
	  warning_msg (error_string);
	  master[k]->isotope = TRUE;
	  master[k]->isotope_ratio = 0.0;
	  master[k]->isotope_ratio_uncertainty = 0.001;
	}
	if (master[k]->isotope == FALSE)
	  continue;
	new_isotopes =
	  (struct isotope *) PHRQ_realloc (new_isotopes,
					   (size_t) (count_isotopes +
						     1) *
					   sizeof (struct isotope));
	if (new_isotopes == NULL)
	  malloc_error ();
	new_isotopes[count_isotopes].master = master[k];
	new_isotopes[count_isotopes].primary = primary_ptr;
	new_isotopes[count_isotopes].isotope_number = isotope_number;
	new_isotopes[count_isotopes].elt_name = master[k]->elt->name;
	new_isotopes[count_isotopes].total = 0;
	new_isotopes[count_isotopes].ratio = master[k]->isotope_ratio;
	new_isotopes[count_isotopes].ratio_uncertainty =
	  master[k]->isotope_ratio_uncertainty;
	sprintf (token, "%d%s", (int) isotope_number, master[k]->elt->name);
	new_isotopes[count_isotopes].isotope_name = string_hsave (token);
	count_isotopes++;
      }
    }
    primary_isotopes = (struct isotope *) free_check_null (primary_isotopes);
    solution_ptr->isotopes =
      (struct isotope *) free_check_null (solution_ptr->isotopes);
    solution_ptr->isotopes = new_isotopes;
    solution_ptr->count_isotopes = count_isotopes;
    qsort (solution_ptr->isotopes,
	   (size_t) count_isotopes,
	   (size_t) sizeof (struct isotope), isotope_compare);
  }

  return (OK);
}

/* ---------------------------------------------------------------------- */
int
tidy_kin_exchange (void)
/* ---------------------------------------------------------------------- */
/*
 *  If exchanger is related to mineral, exchanger amount is 
 *  set in proportion
 */
{
  int i, j, k, n;
  struct kinetics *kinetics_ptr;
  struct exch_comp *comp_ptr;
  struct master *master_ptr;
  char *ptr;
  LDBLE conc;

  for (i = 0; i < count_exchange; i++)
  {
    if (exchange[i].new_def == FALSE)
      continue;
    if (exchange[i].n_user < 0)
      continue;
    for (j = 0; j < exchange[i].count_comps; j++)
    {
      if (exchange[i].comps[j].rate_name == NULL)
	continue;
      comp_ptr = &exchange[i].comps[j];
      comp_ptr->master = NULL;
      n = exchange[i].n_user;

      /* First find exchange master species */

      for (k = 0; comp_ptr->totals[k].elt != NULL; k++)
      {
	/* Find master species */
	master_ptr = comp_ptr->totals[k].elt->master;
	if (master_ptr == NULL)
	{
	  input_error++;
	  sprintf (error_string, "Master species not in data "
		   "base for %s, skipping element.",
		   comp_ptr->totals[k].elt->name);
	  error_msg (error_string, CONTINUE);
	  continue;
	}
	if (master_ptr->type != EX)
	  continue;
	comp_ptr->master = master_ptr;
	break;
      }
      if (comp_ptr->master == NULL)
      {
	input_error++;
	sprintf (error_string,
		 "Exchange formula does not contain an exchange master species, %s",
		 comp_ptr->formula);
	error_msg (error_string, CONTINUE);
	continue;
      }

      /* Now find associated kinetic reaction ...  */
      if ((kinetics_ptr = kinetics_bsearch (n, &k)) == NULL)
      {
	input_error++;
	sprintf (error_string,
		 "Kinetics %d must be defined to use exchange related to kinetic reaction, %s",
		 n, comp_ptr->formula);
	error_msg (error_string, CONTINUE);
	continue;
      }
      for (k = 0; k < kinetics_ptr->count_comps; k++)
      {
	if (strcmp_nocase
	    (comp_ptr->rate_name, kinetics_ptr->comps[k].rate_name) == 0)
	{
	  break;
	}
      }
      if (k == kinetics_ptr->count_comps)
      {
	input_error++;
	sprintf (error_string,
		 "Kinetic reaction, %s, related to exchanger, %s, not found in KINETICS %d",
		 comp_ptr->rate_name, comp_ptr->formula, n);
	error_msg (error_string, CONTINUE);
	continue;
      }

      /* use database name for phase */
      comp_ptr->rate_name = kinetics_ptr->comps[k].rate_name;

      /* make exchanger concentration proportional to mineral ... */
      conc = kinetics_ptr->comps[k].m * comp_ptr->phase_proportion;

      count_elts = 0;
      paren_count = 0;
      ptr = comp_ptr->formula;
      get_elts_in_species (&ptr, conc);
      comp_ptr->totals =
	(struct elt_list *) free_check_null (comp_ptr->totals);
      comp_ptr->totals = elt_list_save ();
/*
 *   No check on availability of exchange elements 
 */
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
tidy_min_exchange (void)
/* ---------------------------------------------------------------------- */
/*
 *  If exchanger is related to mineral, exchanger amount is 
 *  set in proportion
 */
{
  int i, j, k, n, jj;
  struct pp_assemblage *pp_a_ptr;
  struct exch_comp *comp_ptr;
  struct master *master_ptr;
  char *ptr;
  LDBLE conc;

  for (i = 0; i < count_exchange; i++)
  {
    if (exchange[i].new_def == FALSE)
      continue;
    if (exchange[i].n_user < 0)
      continue;
    for (j = 0; j < exchange[i].count_comps; j++)
    {
      if (exchange[i].comps[j].phase_name == NULL)
	continue;
      comp_ptr = &exchange[i].comps[j];
      comp_ptr->master = NULL;
      n = exchange[i].n_user;

      /* First find exchange master species */

      for (k = 0; comp_ptr->totals[k].elt != NULL; k++)
      {
	/* Find master species */
	master_ptr = comp_ptr->totals[k].elt->master;
	if (master_ptr == NULL)
	{
	  input_error++;
	  sprintf (error_string, "Master species not in data "
		   "base for %s, skipping element.",
		   comp_ptr->totals[k].elt->name);
	  error_msg (error_string, CONTINUE);
	  continue;
	}
	if (master_ptr->type != EX)
	  continue;
	comp_ptr->master = master_ptr;
	break;
      }
      if (comp_ptr->master == NULL)
      {
	input_error++;
	sprintf (error_string,
		 "Exchange formula does not contain an exchange master species, %s",
		 comp_ptr->formula);
	error_msg (error_string, CONTINUE);
	continue;
      }

      /* Now find the mineral on which exchanger depends...  */
      if ((pp_a_ptr = pp_assemblage_bsearch (n, &k)) == NULL)
      {
	input_error++;
	sprintf (error_string,
		 "Equilibrium_phases %d must be defined to use exchange related to mineral phase, %s",
		 n, comp_ptr->formula);
	error_msg (error_string, CONTINUE);
	continue;
      }
      for (k = 0; k < pp_a_ptr->count_comps; k++)
      {
	if (strcmp_nocase
	    (comp_ptr->phase_name, pp_a_ptr->pure_phases[k].name) == 0)
	{
	  break;
	}
      }
      if (k == pp_a_ptr->count_comps)
      {
	input_error++;
	sprintf (error_string,
		 "Mineral, %s, related to exchanger, %s, not found in Equilibrium_Phases %d",
		 comp_ptr->phase_name, comp_ptr->formula, n);
	error_msg (error_string, CONTINUE);
	continue;
      }
      /* use database name for phase */
      comp_ptr->phase_name = pp_a_ptr->pure_phases[k].phase->name;
      /* make exchanger concentration proportional to mineral ... */
      conc = pp_a_ptr->pure_phases[k].moles * comp_ptr->phase_proportion;
      count_elts = 0;
      paren_count = 0;
      ptr = comp_ptr->formula;
      get_elts_in_species (&ptr, conc);
      comp_ptr->totals =
	(struct elt_list *) free_check_null (comp_ptr->totals);
      comp_ptr->totals = elt_list_save ();
/*
 *   make sure exchange elements are in phase
 */
      count_elts = 0;
      paren_count = 0;
      ptr = comp_ptr->formula;
      get_elts_in_species (&ptr, -comp_ptr->phase_proportion);
      ptr = pp_a_ptr->pure_phases[k].phase->formula;
      get_elts_in_species (&ptr, 1.0);
      qsort (elt_list, (size_t) count_elts,
	     (size_t) sizeof (struct elt_list), elt_list_compare);
      elt_list_combine ();
      for (jj = 0; jj < count_elts; jj++)
      {
	if (elt_list[jj].elt->primary->s->type != EX && elt_list[jj].coef < 0)
	{
	  input_error++;
	  sprintf (error_string,
		   "Stoichiometry of exchanger, %s * %g mol sites/mol phase,\n\tmust be a subset of the related phase %s, %s.",
		   comp_ptr->formula, (double) comp_ptr->phase_proportion,
		   pp_a_ptr->pure_phases[k].phase->name,
		   pp_a_ptr->pure_phases[k].phase->formula);
	  error_msg (error_string, CONTINUE);
	  break;
	}
      }
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
tidy_min_surface (void)
/* ---------------------------------------------------------------------- */
/*
 *  If surface is related to mineral, surface amount is 
 *  set in proportion
 */
{
  int i, j, k, n, jj;
  struct pp_assemblage *pp_a_ptr;
  struct surface_comp *comp_ptr;
  struct master *master_ptr;
  char *ptr;
  LDBLE conc;

  for (i = 0; i < count_surface; i++)
  {
    if (surface[i].new_def == FALSE)
      continue;
    if (surface[i].n_user < 0)
      continue;
    for (j = 0; j < surface[i].count_comps; j++)
    {
      if (surface[i].comps[j].phase_name == NULL)
	continue;
      comp_ptr = &surface[i].comps[j];
      comp_ptr->master = NULL;
      n = surface[i].n_user;

      /* First find surface master species */

      for (k = 0; comp_ptr->totals[k].elt != NULL; k++)
      {
	/* Find master species */
	master_ptr = comp_ptr->totals[k].elt->master;
	if (master_ptr == NULL)
	{
	  input_error++;
	  sprintf (error_string, "Master species not in data "
		   "base for %s, skipping element.",
		   comp_ptr->totals[k].elt->name);
	  error_msg (error_string, CONTINUE);
	  continue;
	}
	if (master_ptr->type != SURF)
	  continue;
	comp_ptr->master = master_ptr;
	break;
      }
      if (comp_ptr->master == NULL)
      {
	input_error++;
	sprintf (error_string,
		 "Surface formula does not contain a surface master species, %s",
		 comp_ptr->formula);
	error_msg (error_string, CONTINUE);
	continue;
      }

      /* Now find the mineral on which surface depends...  */
      if ((pp_a_ptr = pp_assemblage_bsearch (n, &k)) == NULL)
      {
	input_error++;
	sprintf (error_string,
		 "Equilibrium_phases %d must be defined to use surface related to mineral phase, %s",
		 n, comp_ptr->formula);
	error_msg (error_string, CONTINUE);
	continue;
      }
      for (k = 0; k < pp_a_ptr->count_comps; k++)
      {
	if (strcmp_nocase
	    (comp_ptr->phase_name, pp_a_ptr->pure_phases[k].name) == 0)
	{
	  break;
	}
      }
      if (k == pp_a_ptr->count_comps)
      {
	input_error++;
	sprintf (error_string,
		 "Mineral, %s, related to surface, %s, not found in Equilibrium_Phases %d",
		 comp_ptr->phase_name, comp_ptr->formula, n);
	error_msg (error_string, CONTINUE);
	continue;
      }
      if (pp_a_ptr->pure_phases[k].phase == NULL)
      {
	input_error++;
	sprintf (error_string,
		 "Mineral, %s, related to surface, %s, not found in database.",
		 pp_a_ptr->pure_phases[k].name, comp_ptr->formula);
	error_msg (error_string, CONTINUE);
	continue;
      }
      /* use database name for phase */
      comp_ptr->phase_name = pp_a_ptr->pure_phases[k].phase->name;
      /* make surface concentration proportional to mineral ... */
      conc = pp_a_ptr->pure_phases[k].moles * comp_ptr->phase_proportion;
#ifdef SKIP_MUSIC
      comp_ptr->cb = conc * comp_ptr->formula_z;
#endif
/*			if (conc < MIN_RELATED_SURFACE) conc = 0.0; */
      ptr = comp_ptr->formula;
      count_elts = 0;
      paren_count = 0;
      get_elts_in_species (&ptr, conc);
      comp_ptr->totals =
	(struct elt_list *) free_check_null (comp_ptr->totals);
      comp_ptr->totals = elt_list_save ();

      /* area */
      surface[i].charge[comp_ptr->charge].grams =
	pp_a_ptr->pure_phases[k].moles;
/*
 *   make sure surface elements are in phase
 *   logically necessary for mass balance and to avoid negative concentrations when dissolving phase
 */
      count_elts = 0;
      paren_count = 0;
      ptr = pp_a_ptr->pure_phases[k].phase->formula;
      get_elts_in_species (&ptr, 1.0);
      for (jj = 0; jj < surface[i].count_comps; jj++)
      {
	if (comp_ptr->charge != surface[i].comps[jj].charge)
	  continue;
	if (surface[i].type == CD_MUSIC)
	{
	  ptr = surface[i].comps[jj].formula;
	  get_elts_in_species (&ptr, -surface[i].comps[jj].phase_proportion);
	}
	else
	{
	  if (surface[i].comps[jj].master->s->z != 0.0)
	  {
	    input_error++;
	    sprintf (error_string,
		     "Master species of surface, %s, must be uncharged if the number of sites is related to a phase.",
		     surface[i].comps[jj].master->s->name);
	    error_msg (error_string, CONTINUE);
	  }
	  ptr = surface[i].comps[jj].master->s->name;
	  get_elts_in_species (&ptr, -surface[i].comps[jj].phase_proportion);
	}
      }
      qsort (elt_list, (size_t) count_elts,
	     (size_t) sizeof (struct elt_list), elt_list_compare);
      elt_list_combine ();
      /* Makes no sense: sorbed species need not be in mineral structure... */
      /* But elements that can desorb into solution must be in mineral */
      /* If you precipitate Ca-Mont, and make SurfMg (assuming this is the 
         formula in SURFACE), where does the Mg come from? 
         Further, if you precipitate Ca-Mont, make SurfCa, desorb
         all the Ca, then dissolve the "Ca-Mont", you must remove SurfCa, or you
         will end up with Ca in solution. H and O are excluded */
      for (jj = 0; jj < count_elts; jj++)
      {
	if (elt_list[jj].elt->primary->s->type != SURF
	    && elt_list[jj].coef < 0
	    && elt_list[jj].elt->primary->s != s_hplus
	    && elt_list[jj].elt->primary->s != s_h2o)
	{
	  input_error++;
	  sprintf (error_string, "Element %s in sum of surface sites,\n"
		   "\tincluding %s * %g mol sites/mol phase,\n"
		   "\texceeds stoichiometry in the related phase %s, %s.",
		   elt_list[jj].elt->name,
		   comp_ptr->master->s->name,
		   (double) comp_ptr->phase_proportion,
		   pp_a_ptr->pure_phases[k].phase->name,
		   pp_a_ptr->pure_phases[k].phase->formula);
	  error_msg (error_string, CONTINUE);
	  break;
	}
      }
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
tidy_kin_surface (void)
/* ---------------------------------------------------------------------- */
/*
 *  If surface is related to mineral, surface amount is 
 *  set in proportion
 */
{
  int i, j, k, n, jj, l;
  struct kinetics *kinetics_ptr;
  struct surface_comp *comp_ptr;
  struct master *master_ptr;
  struct phase *phase_ptr;
  char token[MAX_LENGTH];
  char *ptr;
  LDBLE conc;
  struct elt_list *elt_list_kinetics;
  int count_elts_kinetics;

  n = -999;
  comp_ptr = NULL;
  for (i = 0; i < count_surface; i++)
  {
    if (surface[i].new_def == FALSE)
      continue;
    if (surface[i].n_user < 0)
      continue;
    for (j = 0; j < surface[i].count_comps; j++)
    {
      if (surface[i].comps[j].rate_name == NULL)
	continue;
      comp_ptr = &surface[i].comps[j];
      comp_ptr->master = NULL;
      n = surface[i].n_user;

      /* First find surface master species */

      for (k = 0; comp_ptr->totals[k].elt != NULL; k++)
      {
	/* Find master species */
	master_ptr = comp_ptr->totals[k].elt->master;
	if (master_ptr == NULL)
	{
	  input_error++;
	  sprintf (error_string, "Master species not in data "
		   "base for %s, skipping element.",
		   comp_ptr->totals[k].elt->name);
	  error_msg (error_string, CONTINUE);
	  continue;
	}
	if (master_ptr->type != SURF)
	  continue;
	comp_ptr->master = master_ptr;
	break;
      }
      if (comp_ptr->master == NULL)
      {
	input_error++;
	sprintf (error_string,
		 "Surface formula does not contain a surface master species, %s",
		 comp_ptr->formula);
	error_msg (error_string, CONTINUE);
	continue;
      }

      /* Now find the kinetic reaction on which surface depends...  */
      if ((kinetics_ptr = kinetics_bsearch (n, &k)) == NULL)
      {
	input_error++;
	sprintf (error_string,
		 "Kinetics %d must be defined to use surface related to kinetic reaction, %s",
		 n, comp_ptr->formula);
	error_msg (error_string, CONTINUE);
	continue;
      }
      for (k = 0; k < kinetics_ptr->count_comps; k++)
      {
	if (strcmp_nocase
	    (comp_ptr->rate_name, kinetics_ptr->comps[k].rate_name) == 0)
	{
	  break;
	}
      }
      if (k == kinetics_ptr->count_comps)
      {
	input_error++;
	sprintf (error_string,
		 "Kinetic reaction, %s, related to surface, %s, not found in Kinetics %d",
		 comp_ptr->rate_name, comp_ptr->formula, n);
	error_msg (error_string, CONTINUE);
	continue;
      }

      /* use database name for phase */
      comp_ptr->rate_name = kinetics_ptr->comps[k].rate_name;

      /* make surface concentration proportional to mineral ... */
      conc = kinetics_ptr->comps[k].m * comp_ptr->phase_proportion;

/*			if (conc < MIN_RELATED_SURFACE) conc = 0.0; */
      ptr = comp_ptr->formula;
      count_elts = 0;
      paren_count = 0;
      get_elts_in_species (&ptr, conc);
      comp_ptr->totals =
	(struct elt_list *) free_check_null (comp_ptr->totals);
      comp_ptr->totals = elt_list_save ();

      /* area */
      surface[i].charge[comp_ptr->charge].grams = kinetics_ptr->comps[k].m;
    }
/*
 *   check on elements
 */
    /* Go through each kinetic reaction, add all related surface compositions
     * check for negative values
     */
    if (surface[i].related_rate == FALSE)
      continue;
    kinetics_ptr = kinetics_bsearch (n, &k);
    for (k = 0; k < kinetics_ptr->count_comps; k++)
    {
      count_elts = 0;
      paren_count = 0;

      /* added in kinetics formula */
      for (j = 0; j < kinetics_ptr->comps[k].count_list; j++)
      {
	phase_ptr = NULL;
	strcpy (token, kinetics_ptr->comps[k].list[j].name);
	phase_ptr = phase_bsearch (token, &jj, FALSE);
	if (phase_ptr != NULL)
	{
	  add_elt_list (phase_ptr->next_elt, 1.0);
	}
	else
	{
	  ptr = kinetics_ptr->comps[k].list[j].name;
	  get_elts_in_species (&ptr, kinetics_ptr->comps[k].list[j].coef);
	}
      }
      /* save kinetics formula */
      if (count_elts > 0)
      {
	qsort (elt_list, (size_t) count_elts,
	       (size_t) sizeof (struct elt_list), elt_list_compare);
	elt_list_combine ();
      }
      elt_list_kinetics = elt_list_save ();
      count_elts_kinetics = count_elts;

      /* get surface formulas */
      count_elts = 0;
      paren_count = 0;
      for (j = 0; j < surface[i].count_comps; j++)
      {
	comp_ptr = &surface[i].comps[j];
	if (comp_ptr->rate_name == NULL)
	  continue;
	if (strcmp_nocase
	    (comp_ptr->rate_name, kinetics_ptr->comps[k].rate_name) == 0)
	{
	  ptr = comp_ptr->formula;
	  get_elts_in_species (&ptr, -1 * comp_ptr->phase_proportion);
	}
      }
      if (count_elts > 0)
      {
	qsort (elt_list, (size_t) count_elts,
	       (size_t) sizeof (struct elt_list), elt_list_compare);
	elt_list_combine ();
      }
      for (j = 0; j < count_elts; j++)
      {
	if (elt_list[j].elt->primary->s->type <= H2O)
	{
	  for (l = 0; l < count_elts_kinetics; l++)
	  {
	    if (elt_list[j].elt == elt_list_kinetics[l].elt)
	    {
	      break;
	    }
	  }
	  if (l == count_elts_kinetics)
	  {
	    input_error++;
	    sprintf (error_string,
		     "Stoichiometry of surface, %s * %g mol sites/mol reactant,\n\tmust be a subset of the formula defined for the related reactant %s.\n\tElement %s is not present in reactant formula.",
		     comp_ptr->formula, (double) comp_ptr->phase_proportion,
		     comp_ptr->rate_name, elt_list[j].elt->name);
	    error_msg (error_string, CONTINUE);
	  }
	  else if (fabs (elt_list[j].coef) > fabs (elt_list_kinetics[l].coef))
	  {
	    input_error++;
	    sprintf (error_string,
		     "Stoichiometry of surface, %s * %g mol sites/mol reactant,\n\tmust be a subset of the formula defined for the related reactant %s.\n\tCoefficient of element %s in surface exceeds amount present in reactant formula.",
		     comp_ptr->formula, (double) comp_ptr->phase_proportion,
		     comp_ptr->rate_name, elt_list[j].elt->name);
	    error_msg (error_string, CONTINUE);
	  }
	}
      }
      elt_list_kinetics =
	(struct elt_list *) free_check_null (elt_list_kinetics);
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
s_s_prep (LDBLE t, struct s_s *s_s_ptr, int print)
/* ---------------------------------------------------------------------- */
{
  int i, j, k, converged, divisions;
  LDBLE r, rt, ag0, ag1, crit_pt;
  LDBLE xc, tc;
  LDBLE x, x0, x1, xsm1, xsm2, xb1, xb2;
  LDBLE xc1, xc2;
  LDBLE facb1, faca1, spim1, xblm1, acrae, acrael, xliapt, xliapm;
  LDBLE xaly, xaly1, xaly2;
  LDBLE faca, facb, spialy, facal, facbl;
  LDBLE tol;

  if (pr.s_s_assemblage == FALSE)
    print = FALSE;
  tol = 1e-6;
  r = R_KJ_DEG_MOL;
  rt = r * t;
  a0 = s_s_ptr->ag0 / rt;
  a1 = s_s_ptr->ag1 / rt;
  s_s_ptr->a0 = a0;
  s_s_ptr->a1 = a1;
  ag0 = a0 * rt;
  ag1 = a1 * rt;
  kc = exp (k_calc (s_s_ptr->comps[0].phase->rxn->logk, t) * LOG_10);
  kb = exp (k_calc (s_s_ptr->comps[1].phase->rxn->logk, t) * LOG_10);
  crit_pt = fabs (a0) + fabs (a1);
/*
 *   Default, no miscibility or spinodal gaps
 */
  s_s_ptr->miscibility = FALSE;
  s_s_ptr->spinodal = FALSE;
  xsm1 = 0.5;
  xsm2 = 0.5;
  xb1 = 0.5;
  xb2 = 0.5;
  xc1 = 0;
  xc2 = 0;

  if (crit_pt >= tol)
  {
/*
 *   Miscibility gap information
 */
    if (fabs (a1) < tol)
    {
      xc = 0.5;
      tc = ag0 / (2 * r);
    }
    else
    {
      xc = 0.5 + (pow ((ag0 * ag0 + 27 * ag1 * ag1), 0.5) - ag0) / (18 * ag1);
      tc = (12 * ag1 * xc - 6 * ag1 + 2 * ag0) * (xc - xc * xc) / r;
    }
    if (print == TRUE)
    {
      sprintf (error_string, "Description of Solid Solution %s",
	       s_s_ptr->name);
      dup_print (error_string, TRUE);
    }
    if (print == TRUE)
    {
      output_msg (OUTPUT_MESSAGE,
		  "\t                              Temperature: %g kelvin\n",
		  (double) t);
      output_msg (OUTPUT_MESSAGE,
		  "\t                       A0 (dimensionless): %g\n",
		  (double) a0);
      output_msg (OUTPUT_MESSAGE,
		  "\t                       A1 (dimensionless): %g\n",
		  (double) a1);
      output_msg (OUTPUT_MESSAGE,
		  "\t                              A0 (kJ/mol): %g\n",
		  (double) ag0);
      output_msg (OUTPUT_MESSAGE,
		  "\t                              A1 (kJ/mol): %g\n\n",
		  (double) ag1);
    }
    if (xc < 0 || xc > 1)
    {
      if (print == TRUE)
	output_msg (OUTPUT_MESSAGE,
		    "No miscibility gap above 0 degrees kelvin.\n");
    }
    else
    {
      if (print == TRUE)
      {
	output_msg (OUTPUT_MESSAGE,
		    "\t    Critical mole-fraction of component 2: %g\n",
		    (double) xc);
	output_msg (OUTPUT_MESSAGE,
		    "\t                     Critical temperature: %g kelvin\n",
		    (double) tc);
	output_msg (OUTPUT_MESSAGE,
		    "\n(The critical temperature calculation assumes that the Guggenheim model\ndefined at %g kelvin is valid at the critical temperature.)\n\n\n",
		    (double) t);
      }
    }
/*
 *   Calculate miscibility and spinodal gaps
 */
    if (tc >= t)
    {

      /* search for sign changes */
      x0 = 0;
      x1 = 1;
      if (scan (f_spinodal, &x0, &x1) == TRUE)
      {

	/* find first spinodal pt */
	xsm1 = halve (f_spinodal, x0, x1, tol);
	s_s_ptr->spinodal = TRUE;

	/* find second spinodal pt */
	x0 = x1;
	x1 = 1;
	if (scan (f_spinodal, &x0, &x1) == TRUE)
	{
	  xsm2 = halve (f_spinodal, x0, x1, tol);
	}
	else
	{
	  error_msg ("Failed to find second spinodal point.", STOP);
	}
      }
    }
  }
/*
 *   Now find Miscibility gap
 */
  if (s_s_ptr->spinodal == TRUE)
  {
    if (print == TRUE)
      output_msg (OUTPUT_MESSAGE,
		  "\t Spinodal-gap mole fractions, component 2: %g\t%g\n",
		  (double) xsm1, (double) xsm2);
    converged = FALSE;
    if (converged == FALSE)
    {
      for (i = 1; i < 3; i++)
      {
	divisions = (int) pow (10., i);
	for (j = 0; j < divisions; j++)
	{
	  for (k = divisions; k > 0; k--)
	  {
	    xc1 = (LDBLE) j / divisions + 0.001;
	    xc2 = (LDBLE) k / divisions;
	    converged = solve_misc (&xc1, &xc2, tol);
	    if (converged == TRUE)
	      break;
	  }
	  if (converged == TRUE)
	    break;
	}
	if (converged == TRUE)
	  break;
      }
    }
    if (converged == FALSE)
    {
      error_msg ("Failed to find miscibility gap.", STOP);
    }
    s_s_ptr->miscibility = TRUE;
    if (xc1 < xc2)
    {
      xb1 = 1 - xc2;
      xb2 = 1 - xc1;
      xc1 = 1 - xb1;
      xc2 = 1 - xb2;
    }
    else
    {
      xb1 = 1 - xc1;
      xb2 = 1 - xc2;
    }
    facb1 = kb * xb1 * exp (xc1 * xc1 * (a0 + a1 * (4 * xb1 - 1)));
    faca1 = kc * xc1 * exp (xb1 * xb1 * (a0 - a1 * (3 - 4 * xb1)));
    spim1 = log10 (faca1 + facb1);
    xblm1 = 1. / (1. + faca1 / facb1);
    acrae = facb1 / faca1;
    acrael = log10 (acrae);
    xliapt = log10 (facb1);
    xliapm = log10 (faca1);

    if (print == TRUE)
    {
      output_msg (OUTPUT_MESSAGE,
		  "\t   Miscibility-gap fractions, component 2: %g\t%g\n",
		  (double) xb1, (double) xb2);
      output_msg (OUTPUT_MESSAGE, "\n\t\t\tEutectic Point Calculations\n\n");
      output_msg (OUTPUT_MESSAGE,
		  "\t     Aqueous activity ratio (comp2/comp1): %g\n",
		  (double) acrae);
      output_msg (OUTPUT_MESSAGE,
		  "\t Log aqueous activity ratio (comp2/comp1): %g\n",
		  (double) acrael);
      output_msg (OUTPUT_MESSAGE,
		  "\t Aqueous activity fraction of component 2: %g\n",
		  (double) xblm1);
      output_msg (OUTPUT_MESSAGE,
		  "\t                    Log IAP (component 2): %g\n",
		  (double) xliapt);
      output_msg (OUTPUT_MESSAGE,
		  "\t                    Log IAP (component 1): %g\n",
		  (double) xliapm);
      output_msg (OUTPUT_MESSAGE,
		  "\t                               Log Sum Pi: %g\n",
		  (double) spim1);
    }
    s_s_ptr->tk = t;
    s_s_ptr->xb1 = xb1;
    s_s_ptr->xb2 = xb2;
  }
/*
 *   Alyotropic point calculation
 */
  xaly = -1.0;
  x = a0 * a0 + 3 * a1 * a1 + 6 * a1 * log (kb / kc);
  if (x > 0)
  {
    if (fabs (x - a0 * a0) >= tol)
    {
      xaly1 = (-(a0 - 3 * a1) + pow (x, 0.5)) / (6 * a1);
      xaly2 = (-(a0 - 3 * a1) - pow (x, 0.5)) / (6 * a1);
      if (xaly1 >= 0 && xaly1 <= 1)
      {
	xaly = xaly1;
      }
      if (xaly2 >= 0 && xaly2 <= 1)
      {
	xaly = xaly2;
      }
    }
    else
    {
      xaly = 0.5 + log (kb / kc) / (2 * a0);
    }
    if (xaly > 0 && xaly < 1)
    {
      faca = kc * (1 - xaly) * exp (xaly * xaly * (a0 - a1 * (3 - 4 * xaly)));
      facb =
	kb * xaly * exp ((1 - xaly) * (1 - xaly) *
			 (a0 + a1 * (4 * xaly - 1.0)));
      spialy = log10 (faca + facb);
      facal = log10 (faca);
      facbl = log10 (facb);
      if (xaly > xb1 && xaly < xb2)
      {
	if (print == TRUE)
	  output_msg (OUTPUT_MESSAGE,
		      "\nLocal minimum in the solidus curve coresponding to a maximum\nin the minimum stoichiometric saturation curve.\n\n");
      }
      else
      {
	if (print == TRUE)
	  output_msg (OUTPUT_MESSAGE, "\n\t\t\tAlyotropic Point\n\n");
      }
      if (print == TRUE)
      {
	output_msg (OUTPUT_MESSAGE,
		    "\t       Solid mole fraction of component 2: %g\n",
		    (double) xaly);
	output_msg (OUTPUT_MESSAGE,
		    "\t                    Log IAP (component 2): %g\n",
		    (double) facbl);
	output_msg (OUTPUT_MESSAGE,
		    "\t                    Log IAP (component 1): %g\n",
		    (double) facal);
	output_msg (OUTPUT_MESSAGE,
		    "\t                               Log Sum Pi: %g\n",
		    (double) spialy);
      }
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
static LDBLE
halve (LDBLE f (LDBLE x), LDBLE x0, LDBLE x1, LDBLE tol)
/* ---------------------------------------------------------------------- */
{
  int i;
  LDBLE x, y, y0, dx;

  y0 = f (x0);
  dx = (x1 - x0);
/*
 *  Loop for interval halving
 */
  for (i = 0; i < 100; i++)
  {
    dx *= 0.5;
    x = x0 + dx;
    y = f (x);
    if (dx < tol || y == 0)
    {
      break;
    }
#ifdef SKIP
    if (y0 * y < 0)
    {
      x1 = x;
    }
    else
    {
      x0 = x;
      y0 = y;
    }
#endif
    if (y0 * y >= 0)
    {
      x0 = x;
      y0 = y;
    }
  }
  return (x0 + dx);
}

/* ---------------------------------------------------------------------- */
static int
scan (LDBLE f (LDBLE x), LDBLE * xx0, LDBLE * xx1)
/* ---------------------------------------------------------------------- */
{
  int i, j;
  LDBLE fx0, fx1, divisions;
  LDBLE x0, x1, diff;

  x0 = *xx0;
  x1 = *xx1;
  diff = x1 - x0;
  for (j = 0; j < 3; j++)
  {
    fx0 = f (x0);
    divisions = (int) pow ((double) 10, (double) j);
    for (i = 1; i < divisions; i++)
    {
      x1 = *xx0 + diff * (LDBLE) i / divisions;
      fx1 = f (x1);
      if (fx0 * fx1 <= 0)
      {
	*xx0 = x0;
	*xx1 = x1;
	return (TRUE);
      }
      x0 = x1;
      fx0 = fx1;
    }
  }
  return (FALSE);
}

/* ---------------------------------------------------------------------- */
static LDBLE
f_spinodal (LDBLE x)
/* ---------------------------------------------------------------------- */
{
  LDBLE fx;
  fx =
    -12 * a1 * x * x * x + (18 * a1 - 2 * a0) * x * x + (2 * a0 -
							 6 * a1) * x - 1.0;
  return (fx);
}

/* ---------------------------------------------------------------------- */
int
slnq (int n, LDBLE * a, LDBLE * delta, int ncols, int print)
/* ---------------------------------------------------------------------- */
{
  int i, j, k, m;
/* debug
 */
  int row;

  LDBLE b;
/* Debug
*/
  if (print == TRUE)
  {
    output_msg (OUTPUT_MESSAGE, "\nArray in slnq: \n\n");
    for (i = 0; i < ncols - 1; i++)
    {
      row = i * (n + 1);
      for (j = 0; j < ncols; j++)
      {
	output_msg (OUTPUT_MESSAGE, "%10.2e", (double) a[row + j]);
      }
      output_msg (OUTPUT_MESSAGE, "\n");
    }
    output_msg (OUTPUT_MESSAGE, "\n");
  }

  if (n == 0)
    return (OK);
/*   Trivial case */
  if (n == 1)
  {
    if (fabs (a[0]) < ZERO_TOL)
      goto slnq_error;
    delta[0] = a[1] / a[0];
    return (OK);
  }

/*   Reduction loop */
  for (i = 0; i < n - 1; i++)
  {
    b = fabs (a[i * ncols + i]);
    m = i;

/*   Find maximum value in column */
    for (j = i + 1; j < n; j++)
    {
      if (fabs (a[j * ncols + i]) > b)
      {
	b = fabs (a[j * ncols + i]);
	m = j;
      }
    }

/*   Check for singularity */
    if (b < ZERO_TOL)
      goto slnq_error;

/*   Exchange rows if necessary */
    if (m != i)
    {
      for (j = i; j <= n; j++)
      {
	b = a[i * ncols + j];
	a[i * ncols + j] = a[m * ncols + j];
	a[m * ncols + j] = b;
      }
    }

/*   Make a[i][i]=1.0 */
    for (j = n; j >= i; j--)
    {
      a[i * ncols + j] /= a[i * ncols + i];
    }

/*   Reduction step */
    for (j = i + 1; j < n; j++)
    {
      if (a[j * ncols + i] == 0.0)
	continue;
      b = -a[j * ncols + i];
      for (k = i + 1; k <= n; k++)
      {
	a[j * ncols + k] += b * a[i * ncols + k];
      }
    }
  }

/*   Calculation of delta[n] */
  if (fabs (a[(n - 1) * ncols + n - 1]) > ZERO_TOL)
  {
    delta[n - 1] = a[(n - 1) * ncols + n] / a[(n - 1) * ncols + n - 1];
  }
  else
  {
    output_msg (OUTPUT_MESSAGE, "Error: Divide by zero in slnq.\n");
    delta[n] = 0.0;
    goto slnq_error;
  }

/*   Back substitution for other delta values */
  for (i = n - 2; i >= 0; i--)
  {
    delta[i] = a[i * ncols + n];
    for (j = i + 1; j < n; j++)
    {
      delta[i] -= a[i * ncols + j] * delta[j];
    }
  }
  if (print == TRUE)
  {
    output_msg (OUTPUT_MESSAGE, "\nResults from slnq: \n\n");
    for (i = 0; i < n; i++)
    {
      output_msg (OUTPUT_MESSAGE, "%10.2e", (double) delta[i]);
    }
    output_msg (OUTPUT_MESSAGE, "\n");
  }
  return (OK);

slnq_error:{
    sprintf (error_string, "Error: Singular matrix in subroutine slnq. \n");
    warning_msg (error_string);
  }
  return (ERROR);
}

/* ---------------------------------------------------------------------- */
static int
solve_misc (LDBLE * xxc1, LDBLE * xxc2, LDBLE tol)
/* ---------------------------------------------------------------------- */
{
  int i, repeat, converged, max_iter;
  LDBLE x1, x2, xb1, xb2;
  LDBLE xc1, xc1_2, xc1_3, xc2, xc2_2, xc2_3;
  LDBLE lc1, lc2, lb1, lb2;
  LDBLE a[6], d[2];
  LDBLE t;

  xc1 = *xxc1;
  xc2 = *xxc2;
  x1 = 0;
  x2 = 0;
  converged = TRUE;
  max_iter = 25;
  for (i = 0; i < max_iter; i++)
  {
    /*
       output_msg(OUTPUT_MESSAGE, "%d  xc1: %g\txc2 %g\n", i, xc1, xc2);
     */
    xb1 = 1 - xc1;
    xb2 = 1 - xc2;
    xc1_2 = xc1 * xc1;
    xc1_3 = xc1_2 * xc1;
    xc2_2 = xc2 * xc2;
    xc2_3 = xc2_2 * xc2;

    lc1 = exp (xb1 * xb1 * (a0 - a1 * (3 - 4 * xb1)));
    lb1 = exp (xc1 * xc1 * (a0 + a1 * (4 * xb1 - 1)));
    lc2 = exp (xb2 * xb2 * (a0 - a1 * (3 - 4 * xb2)));
    lb2 = exp (xc2 * xc2 * (a0 + a1 * (4 * xb2 - 1)));

    /* -fb */
    a[2] = -(xb1 * lb1 - xb2 * lb2);

    /* -fc */
    a[5] = -(xc1 * lc1 - xc2 * lc2);

    if (fabs (a[2]) < tol && fabs (a[5]) < tol)
      break;

    /* dfb/dxc1 */
    t = exp (a0 * xc1_2 - 4 * a1 * xc1_3 + 3 * a1 * xc1_2);
    a[0] =
      (2 * a0 * xc1 + 6 * a1 * xc1 - 2 * a0 * xc1_2 + 12 * a1 * xc1_3 -
       18 * a1 * xc1_2 - 1) * t;

    /* dfb/dxc2 */
    t = exp (a0 * xc2_2 - 4 * a1 * xc2_3 + 3 * a1 * xc2_2);
    a[1] =
      (2 * a0 * xc2_2 - 12 * a1 * xc2_3 - 2 * a0 * xc2 + 18 * a1 * xc2_2 -
       6 * a1 * xc2 + 1) * t;


    /* dfc/dxc1 */
    t =
      exp (a0 * xc1_2 - 2 * a0 * xc1 + a0 - 4 * a1 * xc1_3 + 9 * a1 * xc1_2 -
	   6 * a1 * xc1 + a1);
    a[3] =
      (2 * a0 * xc1_2 - 2 * a0 * xc1 - 12 * a1 * xc1_3 + 18 * a1 * xc1_2 -
       6 * a1 * xc1 + 1) * t;

    /* dfc/dxc2 */
    t =
      exp (a0 * xc2_2 - 2 * a0 * xc2 + a0 - 4 * a1 * xc2_3 + 9 * a1 * xc2_2 -
	   6 * a1 * xc2 + a1);
    a[4] =
      (-2 * a0 * xc2_2 + 2 * a0 * xc2 + 12 * a1 * xc2_3 - 18 * a1 * xc2_2 +
       6 * a1 * xc2 - 1) * t;


    /* solve for dxc1 and dxc2 */
    slnq (2, a, d, 3, FALSE);

    repeat = TRUE;
    while (repeat == TRUE)
    {
      x1 = xc1 + d[0];
      x2 = xc2 + d[1];
      if (x1 > 1 || x1 < 0 || x2 > 1 || x2 < 0)
      {
	d[0] *= 0.5;
	d[1] *= 0.5;
      }
      else
      {
	repeat = FALSE;
      }
    };
    xc1 = x1;
    xc2 = x2;

    if (fabs (xc1 - xc2) < .01)
    {
      converged = FALSE;
      break;
    }
  }
  if (i == max_iter)
    converged = FALSE;
  *xxc1 = xc1;
  *xxc2 = xc2;
  return (converged);
}

/* ---------------------------------------------------------------------- */
static int
s_s_calc_a0_a1 (struct s_s *s_s_ptr)
/* ---------------------------------------------------------------------- */
{
  int i, done;
  LDBLE r, rt, *p;
  LDBLE q1, q2, xbq1, xbq2, xb1, xb2, xc1, xc2;
  LDBLE r1, r2, pa1, pb1, pa2, pb2, xsm1, xsm2;
  LDBLE pn9, pn10, c5, c6, pl9, pl10, pj9, pj10;
  LDBLE xc, tc;
  LDBLE spialy, azero, phi1, phi2, test;
  LDBLE dq1, dq2, denom, ratio, dr1, dr2, x21, x22, x61, x62;
  LDBLE a0, a1, ag0, ag1;
  LDBLE wg2, wg1, alpha2, alpha3;
  LDBLE kc, kb;
  LDBLE xaly, xcaly, alpha0, alpha1, fx, fx1;
  LDBLE tol;

  tol = 1e-6;
  rt = s_s_ptr->tk * R_KJ_DEG_MOL;
  if (s_s_ptr->comps[0].phase == NULL || s_s_ptr->comps[1].phase == NULL)
  {
    input_error++;
    sprintf (error_string,
	     "Two components were not defined for %s solid solution",
	     s_s_ptr->name);
    error_msg (error_string, CONTINUE);
    return (ERROR);
  }
  kc =
    exp (k_calc (s_s_ptr->comps[0].phase->rxn->logk, s_s_ptr->tk) * LOG_10);
  kb =
    exp (k_calc (s_s_ptr->comps[1].phase->rxn->logk, s_s_ptr->tk) * LOG_10);

  p = s_s_ptr->p;

  a0 = 0;
  a1 = 0;
  ag0 = 0;
  ag1 = 0;
  dq2 = 0;
  switch (s_s_ptr->input_case)
  {
    /*
     *  dimensionless a0 and a1
     */
  case 0:
    a0 = p[0];
    a1 = p[1];
    ag0 = a0 * rt;
    ag1 = a1 * rt;
    break;
    /*
     *  two activity coefficients
     *  q1, q2, xbq1, xbq2
     */
  case 1:
    q1 = p[0];
    q2 = p[1];
    xbq1 = p[2];
    xbq2 = p[3];
    done = FALSE;
    if (fabs (1 - xbq1) > 0 && q1 > 0)
    {
      dq1 = log (q1) / ((1 - xbq1) * (1 - xbq1));
      if (xbq2 <= 0 || xbq2 > 1)
      {
	a0 = dq1;
	a1 = 0;
	done = TRUE;
      }
    }
    if (done == FALSE)
    {
      if (fabs (xbq2) < 0 || q2 <= 0)
      {
	input_error++;
	sprintf (error_string,
		 "No solution possible for A0 and A1 calculation from two activity coefficients, %s.\n",
		 s_s_ptr->name);
	error_msg (error_string, CONTINUE);
	done = TRUE;
      }
    }
    if (done == FALSE)
    {
      dq2 = log (q2) / (xbq2 * xbq2);
      if (xbq1 < 0. || xbq2 > 1.)
      {
	a0 = dq2;
	a1 = 0;
	done = TRUE;
      }
    }
    if (done == FALSE)
    {
      denom = 4 * (xbq1 - xbq2) + 2;
      if (fabs (denom) >= tol)
      {
	if (fabs (1 - xbq1) > 0 && q1 > 0)
	{
	  dq1 = log (q1) / ((1 - xbq1) * (1 - xbq1));
	  a0 = (dq1 * (3 - 4 * xbq2) + dq2 * (4 * xbq1 - 1)) / denom;
	  a1 = (dq1 - dq2) / denom;
	  done = TRUE;
	}
      }
    }
    if (done == FALSE)
    {
      input_error++;
      sprintf (error_string,
	       "No solution possible for A0 and A1 calculation from two activity coefficients, %s.\n",
	       s_s_ptr->name);
      error_msg (error_string, CONTINUE);
    }
    /* io = 1 */
    ag0 = a0 * rt;
    ag1 = a1 * rt;
    break;
    /*
     *  two distribution coefficients
     *  q1, q2, xbq1, xbq2
     */
  case 2:
    q1 = p[0];
    q2 = p[1];
    xbq1 = p[2];
    xbq2 = p[3];
    ratio = kc / kb;
    dr1 = log (q1 / ratio);
    x21 = 2 * xbq1 - 1;
    if (fabs (xbq1 - xbq2) < tol || xbq2 < 0)
    {
      a0 = dr1 / x21;
      a1 = 0;
    }
    else
    {
      dr2 = log (q2 / ratio);
      x22 = 2 * xbq2 - 1;
      if (xbq1 < 0.)
      {
	a0 = dr2 / x22;
	a1 = 0;
      }
      else
      {
	x61 = 6 * xbq1 * xbq1 - 6 * xbq1 + 1;
	x62 = 6 * xbq2 * xbq2 - 6 * xbq2 + 1;
	if (fabs (x22 * x61 - x21 * x62) < tol)
	{
	  input_error++;
	  sprintf (error_string,
		   "No solution possible for A0 and A1 calculation from two distribution coefficients, %s.\n",
		   s_s_ptr->name);
	  error_msg (error_string, CONTINUE);
	}
	a0 = (x61 * dr2 - x62 * dr1) / (x22 * x61 - x21 * x62);
	a1 = (x21 * dr2 - x22 * dr1) / (x21 * x62 - x22 * x61);
      }
    }

    /* io = 1 */
    ag0 = a0 * rt;
    ag1 = a1 * rt;
    break;
    /*
     *  from miscibility gap fractions
     *  q1, q2
     */
  case 3:
    q1 = p[0];
    q2 = p[1];
    xb1 = q1;
    xb2 = q2;
    xc1 = 1 - xb1;
    xc2 = 1 - xb2;
    r1 = log (xb1 / xb2);
    r2 = log (xc1 / xc2);
    pa1 = xc2 * xc2 - xc1 * xc1;
    pb1 =
      3 * (xc2 * xc2 - xc1 * xc1) - 4 * (xc2 * xc2 * xc2 - xc1 * xc1 * xc1);
    pa2 = xb2 * xb2 - xb1 * xb1;
    pb2 =
      -(3 * (xb2 * xb2 - xb1 * xb1) -
	4 * (xb2 * xb2 * xb2 - xb1 * xb1 * xb1));
    a0 = (r1 - pb1 / pb2 * r2) / (pa1 - pa2 * pb1 / pb2);
    a1 = (r1 - pa1 / pa2 * r2) / (pb1 - pb2 * pa1 / pa2);

    /* io = 1 */
    ag0 = a0 * rt;
    ag1 = a1 * rt;
    break;
    /*
     *  from spinodal gap fractions
     *  q1, q2
     */
  case 4:
    q1 = p[0];
    q2 = p[1];
    xsm1 = q1;
    xsm2 = q2;
    pn9 = 1 / xsm1;
    pn10 = 1 / xsm2;
    c5 = 1 - xsm1;
    c6 = 1 - xsm2;
    pl9 = 6 * c5 - 12 * c5 * c5;
    pl10 = 6 * c6 - 12 * c6 * c6;
    pj9 = 2 * c5;
    pj10 = 2 * c6;
    a0 = (pn9 - pl9 / pl10 * pn10) / (pj9 - pl9 / pl10 * pj10);
    a1 = (pn9 - pj9 / pj10 * pn10) / (pl9 - pj9 / pj10 * pl10);

    /* io = 1 */
    ag0 = a0 * rt;
    ag1 = a1 * rt;
    break;
    /*
     *  from critical point
     *  q1, q2
     */
  case 5:
    xc = p[0];
    tc = p[1];
    r = R_KJ_DEG_MOL;
    ag1 = r * tc * (2 * xc - 1) / (12 * xc * xc * (1 - xc) * (1 - xc));
    ag0 = (r * tc / (xc * (1 - xc)) - (12 * xc - 6) * ag1) / 2;

    /* io = 0 */
    a0 = ag0 / rt;
    a1 = ag1 / rt;
    break;
    /*
     *  from alyotropic point
     *  q1, q2
     */
  case 6:
    q1 = p[0];
    q2 = p[1];
    xaly = q1;
    r = log (kb / kc);
    alpha0 = 2 * xaly - 1;
    alpha1 = 6 * xaly * (xaly - 1) + 1;
    spialy = pow (10., q2);
    a0 = -999.;
    a1 = -999.;
    if (fabs (alpha0) < tol)
    {
      input_error++;
      sprintf (error_string,
	       "No solution possible for A0 and A1 calculation from alyotropic point, %s.\n",
	       s_s_ptr->name);
      error_msg (error_string, CONTINUE);
    }
    else
    {
      azero = 1;
      if (fabs (alpha0) > tol)
	azero = r / alpha0;
      xcaly = 1 - xaly;
/*
 *  Solve for a0 by Newton's method
 */
      for (i = 0; i < 50; i++)
      {
	phi1 =
	  xcaly * xcaly * (azero +
			   (r - azero * alpha0) * (4 * xaly - 1) / alpha1);
	phi2 =
	  xaly * xaly * (azero +
			 (3 - 4 * xaly) * (azero * alpha0 - r) / alpha1);
	phi1 = xaly * kb * exp (phi1);
	phi2 = xcaly * kc * exp (phi2);
	fx = phi1 + phi2 - spialy;
	fx1 = xcaly * xcaly * (1 - alpha0 * (4 * xaly - 1) / alpha1) * phi1 +
	  xaly * xaly * (1 + alpha0 * (3 - 4 * xaly) / alpha1) * phi2;
	if (fabs (fx1) < 1e-10)
	{
	  input_error++;
	  sprintf (error_string,
		   "Could not find A0 and A1 calculation from alyotropic point, %s.\n",
		   s_s_ptr->name);
	  error_msg (error_string, CONTINUE);
	  break;
	}
	a0 = azero - fx / fx1;
	test = fabs (a0 - azero) + fabs (fx);
	azero = a0;
	if (test < tol)
	  break;
      }
      if (i == 50)
      {
	input_error++;
	sprintf (error_string,
		 "Too many iterations, could not find A0 and A1 calculation from alyotropic point, %s.\n",
		 s_s_ptr->name);
	error_msg (error_string, CONTINUE);
      }
      else
      {
	a1 = (r - a0 * alpha0) / alpha1;

	/* io = 0 */
	ag0 = a0 * rt;
	ag1 = a1 * rt;
      }

    }
    break;
    /*
     *  dimensional (kJ/mol) Guggenheim parameters
     *  ag0, ag1
     */
  case 7:
    ag0 = p[0];
    ag1 = p[1];
    a0 = ag0 / rt;
    a1 = ag1 / rt;
    break;
    /*
     *  Waldbaum-Thompson
     *  wg2, wg1
     */
  case 8:
    wg2 = p[0];
    wg1 = p[1];
    ag0 = (wg2 + wg1) / 2;
    ag1 = (wg2 - wg1) / 2;
    a0 = ag0 / rt;
    a1 = ag1 / rt;
    break;
    /*
     *  Margules
     *  alpha2, alpha3
     */
  case 9:
    alpha2 = p[0];
    alpha3 = p[1];
    a0 = alpha2 + 3 * alpha3 / 4;
    a1 = alpha3 / 4;
    ag0 = a0 * rt;
    ag1 = a1 * rt;
    break;
  }

  s_s_ptr->ag0 = ag0;
  s_s_ptr->ag1 = ag1;
  s_s_ptr->a0 = a0;
  s_s_ptr->a1 = a1;
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
tidy_master_isotope (void)
/* ---------------------------------------------------------------------- */
{
  int i;
  struct master *master_ptr;

  for (i = 0; i < count_master_isotope; i++)
  {
    /*
     * Mark master species list as minor isotope
     */
    if (master_isotope[i]->minor_isotope == TRUE)
    {
      master_ptr = master_bsearch (master_isotope[i]->name);
      if (master_ptr == NULL)
      {
	input_error++;
	sprintf (error_string, "Did not find master species for isotope, %s",
		 master_isotope[i]->name);
	error_msg (error_string, CONTINUE);
	master_isotope[i]->master = NULL;
	continue;
      }
      else
      {
	master_isotope[i]->master = master_ptr;
      }
      master_ptr->minor_isotope = TRUE;
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
tidy_isotope_ratios (void)
/* ---------------------------------------------------------------------- */
{
  int i;
  struct master *master_ptr;
  struct master_isotope *master_isotope_ptr;
  struct calculate_value *calculate_value_ptr;

  for (i = 0; i < count_isotope_ratio; i++)
  {
    /*
     * Mark master species list as minor isotope
     */
    master_isotope_ptr =
      master_isotope_search (isotope_ratio[i]->isotope_name);
    if (master_isotope_ptr == NULL)
    {
      input_error++;
      sprintf (error_string,
	       "For ISOTOPE_RATIO %s, did not find ISOTOPE definition for this isotope, %s",
	       isotope_ratio[i]->name, isotope_ratio[i]->isotope_name);
      error_msg (error_string, CONTINUE);
    }
    master_ptr = master_bsearch (isotope_ratio[i]->isotope_name);
    if (master_ptr == NULL)
    {
      input_error++;
      sprintf (error_string,
	       "For ISOTOPE_RATIO %s, did not find SOLUTION_MASTER_SPECIES for isotope, %s",
	       isotope_ratio[i]->name, isotope_ratio[i]->isotope_name);
      error_msg (error_string, CONTINUE);
    }
    calculate_value_ptr = calculate_value_search (isotope_ratio[i]->name);
    if (calculate_value_ptr == NULL)
    {
      input_error++;
      sprintf (error_string,
	       "For ISOTOPE_RATIOS %s, did not find corresponding CALCULATE_VALUE definition",
	       isotope_ratio[i]->name);
      error_msg (error_string, CONTINUE);
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
tidy_isotope_alphas (void)
/* ---------------------------------------------------------------------- */
{
  int i;
  struct calculate_value *calculate_value_ptr;
  struct logk *logk_ptr;

  for (i = 0; i < count_isotope_alpha; i++)
  {
    /*
     * Mark master species list as minor isotope
     */
    calculate_value_ptr = calculate_value_search (isotope_alpha[i]->name);
    if (calculate_value_ptr == NULL)
    {
      input_error++;
      sprintf (error_string,
	       "For ISOTOPE_ALPHAS %s, did not find corresponding CALCULATE_VALUE definition",
	       isotope_alpha[i]->name);
      error_msg (error_string, CONTINUE);
    }
    if (isotope_alpha[i]->named_logk != NULL)
    {
      logk_ptr = logk_search (isotope_alpha[i]->named_logk);
      if (logk_ptr == NULL)
      {
	input_error++;
	sprintf (error_string,
		 "For ISOTOPE_ALPHAS %s, did not find corresponding NAMED_EXPRESSION definition %s.",
		 isotope_alpha[i]->name, isotope_alpha[i]->named_logk);
	error_msg (error_string, CONTINUE);
      }
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
reset_last_model (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Initialize model
 */
  last_model.force_prep = TRUE;
  last_model.count_exchange = 0;
  last_model.exchange =
    (struct master **) free_check_null (last_model.exchange);
  last_model.count_gas_phase = 0;
  last_model.gas_phase =
    (struct phase **) free_check_null (last_model.gas_phase);
  last_model.count_s_s_assemblage = 0;
  last_model.s_s_assemblage =
    (char **) free_check_null (last_model.s_s_assemblage);
  last_model.count_pp_assemblage = 0;
  last_model.pp_assemblage =
    (struct phase **) free_check_null (last_model.pp_assemblage);
  last_model.add_formula = (char **) free_check_null (last_model.add_formula);
  last_model.si = (LDBLE *) free_check_null (last_model.si);
  last_model.dl_type = NO_DL;
  last_model.count_surface_comp = 0;
  last_model.surface_comp =
    (char **) free_check_null (last_model.surface_comp);
  last_model.count_surface_charge = 0;
  last_model.surface_charge =
    (char **) free_check_null (last_model.surface_charge);
  return (OK);
}
