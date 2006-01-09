#define EXTERNAL extern
#include "global.h"
#include "phqalloc.h"
#include "output.h"
#include "phrqproto.h"
#include "step.h"

static char const svnid[] = "$Id: step.c 4 2009-04-21 17:29:29Z delucia $";

static int check_pp_assemblage (struct pp_assemblage *pp_assemblage_ptr);
static int gas_phase_check (struct gas_phase *gas_phase_ptr);
static int pp_assemblage_check (struct pp_assemblage *pp_assemblage_ptr);
static int reaction_calc (struct irrev *irrev_ptr);
static int solution_check (void);
static int s_s_assemblage_check (struct s_s_assemblage *s_s_assemblage_ptr);

/* ---------------------------------------------------------------------- */
int
P_step (LDBLE step_fraction)
/* ---------------------------------------------------------------------- */
{
/*
 *   zero global solution, add solution or mixture, add exchange,
 *   add surface, add gas phase, add solid solutions,
 *   set temperature, and add reaction. 
 *   Ensure all elements
 *   included in any of these are present in small amounts.
 *   Save result as n_user -1.
 */
  LDBLE difftemp;
  int step_number;
  struct pp_assemblage *pp_assemblage_save = NULL;
  struct s_s_assemblage *s_s_assemblage_save = NULL;

  if (svnid == NULL)
    fprintf (stderr, " ");

/*
 *   Zero out global solution data
 */

  xsolution_zero ();
/*
 *   Set reaction to zero
 */
  step_x = 0.0;
  step_number = reaction_step;
/*
 *   Mixing or solution
 */
  if (use.mix_ptr != NULL)
  {
    add_mix (use.mix_ptr);
  }
  else if (use.solution_ptr != NULL)
  {
    add_solution (use.solution_ptr, 1.0, 1.0);
  }
  else
  {
    input_error++;
    error_msg ("Neither mixing nor an initial solution have "
	       "been defined in reaction step.", STOP);
  }
/*
 *   Reaction
 */
  if (use.irrev_ptr != NULL)
  {
    add_reaction (use.irrev_ptr, step_number, step_fraction);
  }
/*
 *   Kinetics
 */
  if (use.kinetics_ptr != NULL)
  {
    add_kinetics (use.kinetics_ptr);
    /*
       master_ptr =master_bsearch("S(6)");
       output_msg(OUTPUT_MESSAGE,"Added kinetics, S(6) %e\n", master_ptr->total);
       master_ptr =master_bsearch("S");
       output_msg(OUTPUT_MESSAGE,"Added kinetics, S %e\n", master_ptr->total);
     */
  }
/*
 *   Exchange
 */
  if (use.exchange_ptr != NULL)
  {
    add_exchange (use.exchange_ptr);
  }
/*
 *   Surface
 */
  if (use.surface_ptr != NULL)
  {
    add_surface (use.surface_ptr);
  }
/*
 *   Gases
 */
  if (use.gas_phase_ptr != NULL)
  {
    add_gas_phase (use.gas_phase_ptr);
  }
/*
 *   Temperature
 */
  if (use.temperature_ptr != NULL)
  {
    add_temperature (use.temperature_ptr, step_number);
  }
  if ((state == TRANSPORT) && (transport_step != 0) &&
      (cell > 0) && (cell != count_cells + 1))
  {
    difftemp = tc_x - cell_data[cell - 1].temp;
    cell_data[cell - 1].temp += difftemp / tempr;
    tc_x = cell_data[cell - 1].temp;
  }
/*
 *   Pure phases and solid solutions are added to avoid
 *   zero or negative concentrations
 */
/*
 *   Pure phases
 */
  if (use.pp_assemblage_ptr != NULL)
  {
    pp_assemblage_save =
      (struct pp_assemblage *) PHRQ_malloc (sizeof (struct pp_assemblage));
    if (pp_assemblage_save == NULL)
      malloc_error ();
    pp_assemblage_copy (use.pp_assemblage_ptr, pp_assemblage_save,
			use.pp_assemblage_ptr->n_user);
    add_pp_assemblage (use.pp_assemblage_ptr);
  }
/*
 *   Solid solutions
 */
  if (use.s_s_assemblage_ptr != NULL)
  {
    s_s_assemblage_save =
      (struct s_s_assemblage *) PHRQ_malloc (sizeof (struct s_s_assemblage));
    if (s_s_assemblage_save == NULL)
      malloc_error ();
    s_s_assemblage_copy (use.s_s_assemblage_ptr, s_s_assemblage_save,
			 use.s_s_assemblage_ptr->n_user);
    add_s_s_assemblage (use.s_s_assemblage_ptr);
  }
/*
 *   Check that elements are available for gas components,
 *   pure phases, and solid solutions
 */
  if (use.gas_phase_ptr != NULL)
  {
    gas_phase_check (use.gas_phase_ptr);
  }
  if (use.pp_assemblage_ptr != NULL)
  {
    pp_assemblage_check (use.pp_assemblage_ptr);
  }
  if (use.s_s_assemblage_ptr != NULL)
  {
    s_s_assemblage_check (use.s_s_assemblage_ptr);
  }
/*
 *   Check that element moles are >= zero
 */
  if (solution_check () == MASS_BALANCE)
  {
    /* reset moles and deltas */
    if (use.pp_assemblage_ptr != NULL)
    {
      pp_assemblage_free (use.pp_assemblage_ptr);
      pp_assemblage_copy (pp_assemblage_save, use.pp_assemblage_ptr,
			  use.pp_assemblage_ptr->n_user);
      pp_assemblage_free (pp_assemblage_save);
      pp_assemblage_save =
	(struct pp_assemblage *) free_check_null (pp_assemblage_save);
    }
    if (use.s_s_assemblage_ptr != NULL)
    {
      s_s_assemblage_free (use.s_s_assemblage_ptr);
      s_s_assemblage_copy (s_s_assemblage_save, use.s_s_assemblage_ptr,
			   use.s_s_assemblage_ptr->n_user);
      s_s_assemblage_free (s_s_assemblage_save);
      s_s_assemblage_save =
	(struct s_s_assemblage *) free_check_null (s_s_assemblage_save);
    }
    return (MASS_BALANCE);
  }
/*
 *   Copy global into solution n_user = -1
 */
  xsolution_save (-1);
  step_save_surf (-1);
  step_save_exch (-1);
/*
 *   Clean up temporary space
 */
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
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
xsolution_zero (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Zero out _x variables, master->totals, and species->la
 */
  int i;
/*
 *   Zero totals in master structures
 */
  new_x = FALSE;

  tc_x = 0.0;
  ph_x = 0.0;
  solution_pe_x = 0.0;
  mu_x = 0.0;
  ah2o_x = 0.0;
  density_x = 0.0;
  total_h_x = 0.0;
  total_o_x = 0.0;
  cb_x = 0.0;
  mass_water_aq_x = 0.0;
  units_x = moles_per_kilogram_string;

  for (i = 0; i < count_master; i++)
  {
    master[i]->total = 0.0;
    master[i]->total_primary = 0.0;
    master[i]->s->la = 0.0;
  }
  if (pitzer_model == TRUE)
  {
    for (i = 0; i < count_s; i++)
    {
      s[i]->lg = 0.0;
    }
  }
/*
 *   Copy pe data (not sure this will be used
 */
/*
	pe_data_free (pe_x);
	pe_x = pe_data_alloc();
 */
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
add_solution (struct solution *solution_ptr, LDBLE extensive, LDBLE intensive)
/* ---------------------------------------------------------------------- */
{
/*
 *   Accumulate solution data in master->totals and _x variables.
 *
 *   extensive is multiplication factor for solution
 *   intensive is fraction of all multiplication factors for all solutions
 */
  int i;
  struct master *master_ptr;
  struct species *species_ptr;
/*
 *   Add solution to global variables
 */
  tc_x += solution_ptr->tc * intensive;
  ph_x += solution_ptr->ph * intensive;
  solution_pe_x += solution_ptr->solution_pe * intensive;
  mu_x += solution_ptr->mu * intensive;
  ah2o_x += solution_ptr->ah2o * intensive;
  density_x += solution_ptr->density * intensive;

  total_h_x += solution_ptr->total_h * extensive;
  total_o_x += solution_ptr->total_o * extensive;
  cb_x += solution_ptr->cb * extensive;
  mass_water_aq_x += solution_ptr->mass_water * extensive;
/*
 *   Copy totals data into primary master species
 */
  for (i = 0; solution_ptr->totals[i].description != NULL; i++)
  {
    master_ptr = master_bsearch_primary (solution_ptr->totals[i].description);
    master_ptr->total += solution_ptr->totals[i].moles * extensive;
  }
/*
 *   Accumulate initial guesses for activities
 */
  /*for (i=0; solution_ptr->master_activity[i].description != NULL; i++) { */
  for (i = 0; i < solution_ptr->count_master_activity; i++)
  {
    if (solution_ptr->master_activity[i].description != NULL)
    {
      master_ptr =
	master_bsearch (solution_ptr->master_activity[i].description);
      if (master_ptr != NULL)
      {
	master_ptr->s->la += solution_ptr->master_activity[i].la * intensive;
      }
    }
  }
/*
 *   Accumulate initial guesses for log gamma
 */
  if (pitzer_model == TRUE)
  {
    for (i = 0; i < solution_ptr->count_species_gamma; i++)
    {
      species_ptr = s_search (solution_ptr->species_gamma[i].description);
      if (species_ptr != NULL)
      {
	species_ptr->lg += solution_ptr->species_gamma[i].la * intensive;
      }
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
add_exchange (struct exchange *exchange_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Accumulate exchange data in master->totals and _x variables.
 */
  int i, j;
  struct master *master_ptr;

  if (exchange_ptr == NULL)
    return (OK);
/*
 *   Add element concentrations on exchanger to master species totals
 */
  for (i = 0; i < exchange_ptr->count_comps; i++)
  {
    for (j = 0; exchange_ptr->comps[i].totals[j].elt != NULL; j++)
    {
      master_ptr = exchange_ptr->comps[i].totals[j].elt->primary;
      if (master_ptr->s == s_hplus)
      {
	total_h_x += exchange_ptr->comps[i].totals[j].coef;
      }
      else if (master_ptr->s == s_h2o)
      {
	total_o_x += exchange_ptr->comps[i].totals[j].coef;
      }
      else
      {
	master_ptr->total += exchange_ptr->comps[i].totals[j].coef;
      }
    }
  }
  if (exchange_ptr->new_def == TRUE)
  {
    for (i = 0; i < count_master; i++)
    {
      if (master[i]->type == EX && master[i]->total > 0)
      {
	master[i]->s->la = log10 (0.1 * master[i]->total);
      }
    }
  }
  else
  {
    for (i = 0; i < exchange_ptr->count_comps; i++)
    {
      exchange_ptr->comps[i].master->s->la = exchange_ptr->comps[i].la;
      cb_x += exchange_ptr->comps[i].charge_balance;
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
add_surface (struct surface *surface_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Accumulate surface data in master->totals and _x variables.
 */
  int i, j;
  struct master *master_ptr;

  if (surface_ptr == NULL)
    return (OK);
/*
 *   Add element concentrations on surface to master species totals
 */
  dl_type_x = surface_ptr->dl_type;
  for (i = 0; i < surface_ptr->count_comps; i++)
  {
    /*if(surface_ptr->edl == FALSE) { */
    if (surface_ptr->type == NO_EDL)
    {
      cb_x += surface_ptr->comps[i].cb;
    }
#ifdef SKIP_MUSIC
    if (surface_ptr->type == CD_MUSIC)
    {
      cb_x += surface_ptr->comps[i].cb;
    }
#endif
    if (surface_ptr->new_def == FALSE)
    {
      surface_ptr->comps[i].master->s->la = surface_ptr->comps[i].la;
    }
/*
 *   Add surface and specifically sorbed elements
 */
    for (j = 0; surface_ptr->comps[i].totals[j].elt != NULL; j++)
    {
      master_ptr = surface_ptr->comps[i].totals[j].elt->primary;
      if (master_ptr == NULL)
      {
	input_error++;
	sprintf (error_string, "Element not defined in database, %s.",
		 surface_ptr->comps[i].totals[j].elt->name);
	error_msg (error_string, STOP);
      }
      if (master_ptr->s == s_hplus)
      {
	total_h_x += surface_ptr->comps[i].totals[j].coef;
      }
      else if (master_ptr->s == s_h2o)
      {
	total_o_x += surface_ptr->comps[i].totals[j].coef;
      }
      else
      {
	master_ptr->total += surface_ptr->comps[i].totals[j].coef;
      }
    }
  }
  /*if (surface_ptr->edl == FALSE) return(OK); */
  if (surface_ptr->type != DDL && surface_ptr->type != CD_MUSIC)
    return (OK);
  for (i = 0; i < surface_ptr->count_charge; i++)
  {
    /*if (surface_ptr->edl == TRUE) { */
    /*cb_x += surface_ptr->charge[i].charge_balance; */
    if (surface_ptr->type == DDL || surface_ptr->type == CD_MUSIC)
    {
      cb_x += surface_ptr->charge[i].charge_balance;
    }
    if (surface_ptr->new_def == FALSE)
    {
      master_ptr =
	surface_get_psi_master (surface_ptr->charge[i].name, SURF_PSI);
      master_ptr->s->la = surface_ptr->charge[i].la_psi;
      /*surface_ptr->charge[i].psi_master->s->la = surface_ptr->charge[i].la_psi; */
    }
/*
 *   Add diffuse layer elements (including water in Debye layer)
 */
    if (surface_ptr->dl_type != NO_DL && surface_ptr->new_def == FALSE)
    {
      for (j = 0; surface_ptr->charge[i].diffuse_layer_totals[j].elt != NULL;
	   j++)
      {
	master_ptr =
	  surface_ptr->charge[i].diffuse_layer_totals[j].elt->primary;
	if (master_ptr->s == s_hplus)
	{
	  total_h_x += surface_ptr->charge[i].diffuse_layer_totals[j].coef;
	}
	else if (master_ptr->s == s_h2o)
	{
	  total_o_x += surface_ptr->charge[i].diffuse_layer_totals[j].coef;
	}
	else
	{
	  master_ptr->total +=
	    surface_ptr->charge[i].diffuse_layer_totals[j].coef;
	}
      }
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
add_mix (struct mix *mix_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   calls add_solution to accumulate all data in master->totals
 *   and other variables.
 */
  int i;
  int n;
  LDBLE sum_fractions, intensive, extensive;
  struct solution *solution_ptr;
  int count_positive;
  LDBLE sum_positive;

  if (mix_ptr == NULL)
    return (OK);
  if (mix_ptr->count_comps <= 0)
    return (OK);
  sum_fractions = 0.0;
  sum_positive = 0.0;
  count_positive = 0;
  for (i = 0; i < mix_ptr->count_comps; i++)
  {
    sum_fractions += mix_ptr->comps[i].fraction;
    if (mix_ptr->comps[i].fraction > 0)
    {
      sum_positive += mix_ptr->comps[i].fraction;
      count_positive++;
    }
  }
  for (i = 0; i < mix_ptr->count_comps; i++)
  {
    solution_ptr = solution_bsearch (mix_ptr->comps[i].n_solution, &n, TRUE);
    if (solution_ptr == NULL)
    {
      input_error++;
      continue;
    }
    extensive = mix_ptr->comps[i].fraction;
    intensive = extensive / sum_fractions;
    if (count_positive < mix_ptr->count_comps)
    {
      if (mix_ptr->comps[i].fraction > 0)
      {
	intensive = extensive / sum_positive;
      }
      else
      {
	intensive = 0;
      }
    }
    add_solution (solution_ptr, extensive, intensive);
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
add_pp_assemblage (struct pp_assemblage *pp_assemblage_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Add a small amount of each phase if necessary to insure
 *   all elements exist in solution.
 */
  int i, j;
  LDBLE amount_to_add, total;
  char token[MAX_LENGTH];
  char *ptr;
  struct pure_phase *pure_phase_ptr;
  struct master *master_ptr;

  if (check_pp_assemblage (pp_assemblage_ptr) == OK)
    return (OK);
/*
 *   Go through list and generate list of elements and
 *   coefficient of elements in reaction
 */
  count_elts = 0;
  paren_count = 0;
/*
 *   Check that all elements are in solution for phases with greater than zero mass
 */
  pure_phase_ptr = pp_assemblage_ptr->pure_phases;
  for (j = 0; j < pp_assemblage_ptr->count_comps; j++)
  {
    count_elts = 0;
    paren_count = 0;
    amount_to_add = 0.0;
    pure_phase_ptr[j].delta = 0.0;
    if (pure_phase_ptr[j].add_formula != NULL)
    {
      strcpy (token, pure_phase_ptr[j].add_formula);
      ptr = &(token[0]);
      get_elts_in_species (&ptr, 1.0);
    }
    else
    {
      strcpy (token, pure_phase_ptr[j].phase->formula);
      add_elt_list (pure_phase_ptr[j].phase->next_elt, 1.0);
    }
    if (pure_phase_ptr[j].moles > 0.0)
    {
      for (i = 0; i < count_elts; i++)
      {
	master_ptr = elt_list[i].elt->primary;
	if (master_ptr->s == s_hplus)
	{
	  continue;
	}
	else if (master_ptr->s == s_h2o)
	{
	  continue;
	}
	else if (master_ptr->total > MIN_TOTAL)
	{
	  continue;
	}
	else
	{
	  total = (-master_ptr->total + 1e-10) / elt_list[i].coef;
	  if (amount_to_add < total)
	  {
	    amount_to_add = total;
	  }
	}
      }
      if (pure_phase_ptr[j].moles < amount_to_add)
      {
	amount_to_add = pure_phase_ptr[j].moles;
      }
    }
    if (amount_to_add > 0.0)
    {
      pure_phase_ptr[j].moles -= amount_to_add;
      pure_phase_ptr[j].delta = amount_to_add;
/*
 *   Add reaction to totals
 */
      for (i = 0; i < count_elts; i++)
      {
	master_ptr = elt_list[i].elt->primary;
	if (master_ptr->s == s_hplus)
	{
	  total_h_x += elt_list[i].coef * amount_to_add;
	}
	else if (master_ptr->s == s_h2o)
	{
	  total_o_x += elt_list[i].coef * amount_to_add;
	}
	else
	{
	  master_ptr->total += elt_list[i].coef * amount_to_add;
	}
      }
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
static int
check_pp_assemblage (struct pp_assemblage *pp_assemblage_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Check list of all elements in pure_phase assemblage to see
 *   if all are in model. Return true if all are present,
 *   Return false if one or more is missing.
 */
  int j;
  struct master *master_ptr;
  for (j = 0; pp_assemblage_ptr->next_elt[j].elt != NULL; j++)
  {
    master_ptr = pp_assemblage_ptr->next_elt[j].elt->primary;
    if (master_ptr->s == s_h2o || master_ptr->s == s_hplus)
      continue;
    if (master_ptr->total > MIN_TOTAL)
      continue;
    return (FALSE);
  }
  return (TRUE);
}

/* ---------------------------------------------------------------------- */
int
add_reaction (struct irrev *irrev_ptr, int step_number, LDBLE step_fraction)
/* ---------------------------------------------------------------------- */
{
/*
 *   Add irreversible reaction
 */
  int i;
  char c;
  struct master *master_ptr;
/*
 *   Calculate and save reaction
 */
/* !!!!! with kinetics reaction, coeff's may change 
 *       and reaction_calc must be called ....
 */
  if (irrev_ptr->elts == NULL)
  {
    if (reaction_calc (irrev_ptr) == ERROR)
    {
      return (ERROR);
    }
  }
/*
 *   Step size
 */
  if (incremental_reactions == FALSE)
  {
    if (irrev_ptr->count_steps > 0)
    {
      if (step_number > irrev_ptr->count_steps)
      {
	step_x = irrev_ptr->steps[irrev_ptr->count_steps - 1];
      }
      else
      {
	step_x = irrev_ptr->steps[step_number - 1];
      }
    }
    else if (irrev_ptr->count_steps < 0)
    {
      if (step_number > -irrev_ptr->count_steps)
      {
	step_x = irrev_ptr->steps[0];
      }
      else
      {
	step_x = irrev_ptr->steps[0] *
	  ((LDBLE) step_number) / ((LDBLE) (-irrev_ptr->count_steps));
      }
    }
    else
    {
      step_x = 0.0;
    }
  }
  else
  {
    /* Incremental reactions */
    if (irrev_ptr->count_steps > 0)
    {
      if (step_number > irrev_ptr->count_steps)
      {
	step_x = irrev_ptr->steps[irrev_ptr->count_steps - 1];
      }
      else
      {
	step_x = irrev_ptr->steps[step_number - 1];
      }
    }
    else if (irrev_ptr->count_steps < 0)
    {
      if (step_number > -irrev_ptr->count_steps)
      {
	step_x = 0;
      }
      else
      {
	step_x = irrev_ptr->steps[0] / ((LDBLE) (-irrev_ptr->count_steps));
      }
    }
    else
    {
      step_x = 0.0;
    }
  }
/*
 *   Convert units
 */
  c = irrev_ptr->units[0];
  if (c == 'm')
  {
    step_x *= 1e-3;
  }
  else if (c == 'u')
  {
    step_x *= 1e-6;
  }
  else if (c == 'n')
  {
    step_x *= 1e-9;
  }
/*
 *   Add reaction to totals
 */
  for (i = 0; irrev_ptr->elts[i].elt != NULL; i++)
  {
    master_ptr = irrev_ptr->elts[i].elt->primary;
    if (master_ptr->s == s_hplus)
    {
      total_h_x += irrev_ptr->elts[i].coef * step_x * step_fraction;
    }
    else if (master_ptr->s == s_h2o)
    {
      total_o_x += irrev_ptr->elts[i].coef * step_x * step_fraction;
    }
    else
    {
      master_ptr->total += irrev_ptr->elts[i].coef * step_x * step_fraction;
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
static int
reaction_calc (struct irrev *irrev_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *    Go through irreversible reaction initially to
 *    determine a list of elements and amounts in 
 *    the reaction.
 */
  int i, j, return_value;
  LDBLE coef;
  char token[MAX_LENGTH];
  char *ptr;
  struct phase *phase_ptr;
/*
 *   Go through list and generate list of elements and
 *   coefficient of elements in reaction
 */
  return_value = OK;
  count_elts = 0;
  paren_count = 0;

  for (i = 0; i < irrev_ptr->count_list; i++)
  {
    coef = irrev_ptr->list[i].coef;
    strcpy (token, irrev_ptr->list[i].name);
    phase_ptr = phase_bsearch (token, &j, FALSE);
/*
 *   Reactant is a pure phase, copy formula into token
 */
    if (phase_ptr != NULL)
    {
      add_elt_list (phase_ptr->next_elt, coef);
    }
    else
    {
      ptr = &(token[0]);
      get_elts_in_species (&ptr, coef);
    }
  }
/*
 *   Check that all elements are in database
 */
  for (i = 0; i < count_elts; i++)
  {
    if (elt_list[i].elt->master == NULL)
    {
      sprintf (error_string, "Element or phase not defined in database, %s.",
	       elt_list[i].elt->name);
      error_msg (error_string, CONTINUE);
      input_error++;
      return_value = ERROR;
    }
  }
  irrev_ptr->elts = elt_list_save ();

  return (return_value);
}

/* ---------------------------------------------------------------------- */
int
add_temperature (struct temperature *temperature_ptr, int step_number)
/* ---------------------------------------------------------------------- */
{
/*
 *   Determine temperature of reaction step, if reaction_temperature
 *   information is present.
 */
  int denom;
  LDBLE tc_temp;
/*
 *   Find temperature
 */
  if (temperature_ptr == NULL)
    return (ERROR);
  if (temperature_ptr->count_t > 0)
  {
    if (step_number > temperature_ptr->count_t)
    {
      tc_temp = temperature_ptr->t[temperature_ptr->count_t - 1];
    }
    else
    {
      tc_temp = temperature_ptr->t[step_number - 1];
    }
  }
  else if (temperature_ptr->count_t < 0)
  {
    if (step_number > -temperature_ptr->count_t)
    {
      tc_temp = temperature_ptr->t[1];
    }
    else
    {
      if (-temperature_ptr->count_t <= 1)
      {
	denom = 1;
      }
      else
      {
	denom = -temperature_ptr->count_t - 1;
      }
      tc_temp =
	temperature_ptr->t[0] + (temperature_ptr->t[1] -
				 temperature_ptr->t[0]) *
	((LDBLE) (step_number - 1)) / ((LDBLE) denom);
    }
  }
  else
  {
    tc_temp = 25.0;
  }
  tc_x = tc_temp;
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
add_gas_phase (struct gas_phase *gas_phase_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Accumulate gas data in master->totals and _x variables.
 */
  int i;

  struct gas_comp *gas_comp_ptr;
  struct master *master_ptr;

  if (gas_phase_ptr == NULL)
    return (OK);
  gas_comp_ptr = gas_phase_ptr->comps;
/*
 *   calculate reaction
 */
  count_elts = 0;
  paren_count = 0;
  for (i = 0; i < gas_phase_ptr->count_comps; i++)
  {
    add_elt_list (gas_comp_ptr[i].phase->next_elt, gas_comp_ptr[i].moles);
  }
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
 *   Add gas elements to totals
 */
  for (i = 0; i < count_elts; i++)
  {
    master_ptr = elt_list[i].elt->primary;
    if (master_ptr->s == s_hplus)
    {
      total_h_x += elt_list[i].coef;
    }
    else if (master_ptr->s == s_h2o)
    {
      total_o_x += elt_list[i].coef;
    }
    else
    {
      master_ptr->total += elt_list[i].coef;
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
add_s_s_assemblage (struct s_s_assemblage *s_s_assemblage_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Accumulate solid_solution data in master->totals and _x variables.
 */
  int i, j, k;
  LDBLE amount_to_add, total;
  struct s_s *s_s_ptr;
  struct master *master_ptr;
  char token[MAX_LENGTH];
  char *ptr;

  if (s_s_assemblage_ptr == NULL)
    return (OK);
  count_elts = 0;
  paren_count = 0;
/*
 *   Check that all elements are in solution for phases with greater than zero mass
 */
  for (i = 0; i < s_s_assemblage_ptr->count_s_s; i++)
  {
    count_elts = 0;
    paren_count = 0;
    s_s_ptr = &(s_s_assemblage_ptr->s_s[i]);
    for (j = 0; j < s_s_ptr->count_comps; j++)
    {
      amount_to_add = 0.0;
      s_s_ptr->comps[j].delta = 0.0;
      if (s_s_ptr->comps[j].moles > 0.0)
      {
	strcpy (token, s_s_ptr->comps[j].phase->formula);
	ptr = &(token[0]);
	get_elts_in_species (&ptr, 1.0);
	for (k = 0; k < count_elts; k++)
	{
	  master_ptr = elt_list[k].elt->primary;
	  if (master_ptr->s == s_hplus)
	  {
	    continue;
	  }
	  else if (master_ptr->s == s_h2o)
	  {
	    continue;
	  }
	  else if (master_ptr->total > MIN_TOTAL_SS)
	  {
	    continue;
	  }
	  else
	  {
	    total = (-master_ptr->total + 1e-10) / elt_list[k].coef;
	    if (amount_to_add < total)
	    {
	      amount_to_add = total;
	    }
	  }
	}
      }
      if (s_s_ptr->comps[j].moles < amount_to_add)
      {
	amount_to_add = s_s_ptr->comps[j].moles;
      }
      if (amount_to_add > 0.0)
      {
	s_s_ptr->comps[j].moles -= amount_to_add;
	s_s_ptr->comps[j].delta = amount_to_add;
/*
 *   Add reaction to totals
 */
	for (k = 0; k < count_elts; k++)
	{
	  master_ptr = elt_list[k].elt->primary;
	  if (master_ptr->s == s_hplus)
	  {
	    total_h_x += elt_list[k].coef * amount_to_add;
	  }
	  else if (master_ptr->s == s_h2o)
	  {
	    total_o_x += elt_list[k].coef * amount_to_add;
	  }
	  else
	  {
	    master_ptr->total += elt_list[k].coef * amount_to_add;
	  }
	}
      }
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
add_kinetics (struct kinetics *kinetics_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Add kinetic reaction
 */
  int i;
  struct master *master_ptr;
/*
 *   Add reaction to totals
 */
  if (kinetics_ptr->totals == NULL)
    return (OK);
  for (i = 0; kinetics_ptr->totals[i].elt != NULL; i++)
  {
    master_ptr = kinetics_ptr->totals[i].elt->primary;
    if (master_ptr == NULL)
    {
      input_error++;
      sprintf (error_string,
	       "Element %s in kinetic reaction not found in database.",
	       kinetics_ptr->totals[i].elt->name);
      error_msg (error_string, STOP);
    }
    if (master_ptr->s == s_hplus)
    {
      total_h_x += kinetics_ptr->totals[i].coef;
    }
    else if (master_ptr->s == s_h2o)
    {
      total_o_x += kinetics_ptr->totals[i].coef;
    }
    else
    {
      master_ptr->total += kinetics_ptr->totals[i].coef;
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
gas_phase_check (struct gas_phase *gas_phase_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Check for missing elements
 */
  int i, j;

  struct gas_comp *gas_comp_ptr;
  struct master *master_ptr;

  if (gas_phase_ptr == NULL)
    return (OK);
  gas_comp_ptr = gas_phase_ptr->comps;
/*
 *   Check that all elements are in solution for phases with zero mass
 */
  for (i = 0; i < gas_phase_ptr->count_comps; i++)
  {
    count_elts = 0;
    paren_count = 0;
    if (gas_comp_ptr[i].moles <= 0.0)
    {
      add_elt_list (gas_comp_ptr[i].phase->next_elt, 1.0);
      for (j = 0; j < count_elts; j++)
      {
	master_ptr = elt_list[j].elt->primary;
	if (master_ptr->s == s_hplus)
	{
	  continue;
	}
	else if (master_ptr->s == s_h2o)
	{
	  continue;
	}
	else if (master_ptr->total > MIN_TOTAL)
	{
	  continue;
	}
	else
	{
	  if (state != ADVECTION && state != TRANSPORT && state != PHAST)
	  {
	    sprintf (error_string,
		     "Element %s is contained in gas %s (which has 0.0 mass),\nbut is not in solution or other phases.",
		     elt_list[j].elt->name, gas_comp_ptr[i].phase->name);
	    warning_msg (error_string);
	  }
	}
      }
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
pp_assemblage_check (struct pp_assemblage *pp_assemblage_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Check for missing elements
 */
  int i, j, k;
  char token[MAX_LENGTH];
  char *ptr;
  struct pure_phase *pure_phase_ptr;
  struct master *master_ptr;

  if (check_pp_assemblage (pp_assemblage_ptr) == OK)
    return (OK);
/*
 *   Check that all elements are in solution for phases with zero mass
 */
  pure_phase_ptr = pp_assemblage_ptr->pure_phases;
  for (j = 0; j < pp_assemblage_ptr->count_comps; j++)
  {
    count_elts = 0;
    paren_count = 0;
    if (pure_phase_ptr[j].moles <= 0.0)
    {
      pure_phase_ptr[j].delta = 0.0;
      if (pure_phase_ptr[j].add_formula != NULL)
      {
	strcpy (token, pure_phase_ptr[j].add_formula);
	ptr = &(token[0]);
	get_elts_in_species (&ptr, 1.0);
      }
      else
      {
	strcpy (token, pure_phase_ptr[j].phase->formula);
	add_elt_list (pure_phase_ptr[j].phase->next_elt, 1.0);
      }
      for (i = 0; i < count_elts; i++)
      {
	master_ptr = elt_list[i].elt->primary;
	if (master_ptr->s == s_hplus)
	{
	  continue;
	}
	else if (master_ptr->s == s_h2o)
	{
	  continue;
	}
	else if (master_ptr->total > MIN_TOTAL)
	{
	  continue;
	}
	else
	{
	  if (state != ADVECTION && state != TRANSPORT && state != PHAST)
	  {
	    sprintf (error_string,
		     "Element %s is contained in %s (which has 0.0 mass),"
		     "\t\nbut is not in solution or other phases.",
		     elt_list[i].elt->name, pure_phase_ptr[j].phase->name);
	    warning_msg (error_string);
	  }
/*
 *   Make la's of all master species for the element small, so SI will be small
 *   and no mass transfer will be calculated
 */
	  for (k = 0; k < count_master; k++)
	  {
	    if (master[k]->elt->primary == master_ptr)
	    {
	      master[k]->s->la = -9999.999;
	    }
	  }
	}
      }
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
s_s_assemblage_check (struct s_s_assemblage *s_s_assemblage_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Check for missing elements
 */
  int i, j, k, l;
  struct master *master_ptr;

  if (s_s_assemblage_ptr == NULL)
    return (OK);
/*
 *   Check that all elements are in solution for phases with zero mass
 */
  for (i = 0; i < s_s_assemblage_ptr->count_s_s; i++)
  {
    for (j = 0; j < s_s_assemblage_ptr->s_s[i].count_comps; j++)
    {
      count_elts = 0;
      paren_count = 0;
      if (s_s_assemblage_ptr->s_s[i].comps[j].moles <= 0.0)
      {
	add_elt_list (s_s_assemblage_ptr->s_s[i].comps[j].phase->next_elt,
		      1.0);
	for (l = 0; l < count_elts; l++)
	{
	  master_ptr = elt_list[l].elt->primary;
	  if (master_ptr->s == s_hplus)
	  {
	    continue;
	  }
	  else if (master_ptr->s == s_h2o)
	  {
	    continue;
	  }
	  else if (master_ptr->total > MIN_TOTAL_SS)
	  {
	    continue;
	  }
	  else
	  {
	    if (state != ADVECTION && state != TRANSPORT && state != PHAST)
	    {
	      sprintf (error_string,
		       "Element %s is contained in solid solution %s (which has 0.0 mass),\nbut is not in solution or other phases.",
		       elt_list[l].elt->name,
		       s_s_assemblage_ptr->s_s[i].comps[j].phase->name);
	      warning_msg (error_string);
	    }
	  }
	  /*
	   *   Make la's of all master species for the element small, 
	   *   so SI will be small
	   *   and no mass transfer will be calculated
	   */
	  for (k = 0; k < count_master; k++)
	  {
	    if (master[k]->elt->primary == master_ptr)
	    {
	      master[k]->s->la = -9999.999;
	    }
	  }
	}
      }
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
solution_check (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Check for missing elements
 */
  int i;
  struct master *master_ptr;

/*
 *   Check that all elements are in solution for phases with zero mass
 */
  for (i = 0; i < count_master; i++)
  {
    master_ptr = master[i];
    if (master_ptr->total >= 0.0)
      continue;
    if (master_ptr->total > -MIN_TOTAL)
    {
      master_ptr->total = 0;
      continue;
    }
    if (master_ptr->s == s_eminus || master_ptr->s == s_h2o
	|| master_ptr->s == s_hplus || master_ptr->s == s_h3oplus)
    {
      master_ptr->total = 0;
      continue;
    }
    /*
    sprintf (error_string,
	     "Element %s has negative moles in solution, %e. \n\tErroneous mole balance occurs as moles are added to produce zero moles.\n\tUsually caused by KINETICS, REACTION, or diffuse layer calculation.\n\tMay be due to large time steps in early part of KINETICS simulation or negative concentrations in the diffuse layer.",
	     master_ptr->elt->name, (double) master_ptr->total);
    */
    sprintf (error_string,
	     "Negative moles in solution for %s, %e. Recovering...", master_ptr->elt->name, (double) master_ptr->total);
    warning_msg (error_string);
    return (MASS_BALANCE);
  }
  return (OK);
}
