#define EXTERNAL extern
#include "global.h"
#include "phqalloc.h"
#include "output.h"
#include "phrqproto.h"
#define PITZER_EXTERNAL extern
#include "pitzer.h"
#include <assert.h>

static char const svnid[] = "$Id: print.c 4 2009-04-21 17:29:29Z delucia $";

static int print_alkalinity (void);
static int print_diffuse_layer (struct surface_charge *surface_charge_ptr);
static int print_eh (void);
static int print_irrev (void);
static int print_kinetics (void);
static int print_mix (void);
static int print_pp_assemblage (void);
static int print_s_s_assemblage (void);
static int print_saturation_indices (void);
static int print_surface_cd_music (void);
static int print_totals (void);
static int print_using (void);
/*static int print_user_print(void);*/
static int punch_gas_phase (void);
static int punch_identifiers (void);
static int punch_kinetics (void);
static int punch_molalities (void);
static int punch_activities (void);
static int punch_pp_assemblage (void);
static int punch_s_s_assemblage (void);
static int punch_saturation_indices (void);
static int punch_totals (void);
// MDL: prototype to phrqproto.h
// static int punch_user_punch (void);

#ifdef PHREEQ98
static int punch_user_graph (void);
extern int colnr, rownr;
extern int graph_initial_solutions;
extern int prev_advection_step, prev_transport_step;	/*, prev_reaction_step */
/* extern int shifts_as_points; */
extern int chart_type;
extern int AddSeries;
extern int FirstCallToUSER_GRAPH;
#endif

#if defined(SWIG_SHARED_OBJ)
extern int EndRow (void);
extern void AddSelectedOutput (const char *name, const char *format,
			       va_list argptr);
#endif

/* ---------------------------------------------------------------------- */
int
array_print (LDBLE * array_l, int row_count, int column_count,
	     int max_column_count)
/* ---------------------------------------------------------------------- */
{
  int i, j, k;
  if (svnid == NULL)
    fprintf (stderr, " ");


  for (i = 0; i < row_count; i++)
  {
    k = 0;
    output_msg (OUTPUT_MESSAGE, "%d\n", i);
    for (j = 0; j < column_count; j++)
    {
      if (k > 7)
      {
	output_msg (OUTPUT_MESSAGE, "\n");
	k = 0;
      }
      output_msg (OUTPUT_MESSAGE, "%11.2e",
		  (double) array_l[i * max_column_count + j]);
      k++;
    }
    if (k != 0)
    {
      output_msg (OUTPUT_MESSAGE, "\n");
    }
    output_msg (OUTPUT_MESSAGE, "\n");
  }
  output_msg (OUTPUT_MESSAGE, "\n");
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
print_all (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Print each block of results in sequence to file output
 *   also print selected output to punch_file
 *   Each routine is controlled by a variable in structure print.
 *   print.all == FALSE will turn off all prints.
 */
/*
 *   Makes sorted list including all species for each valence state
 */
  if (pr.all == FALSE)
    return (OK);
  if (pr.surface == TRUE || pr.exchange == TRUE || pr.species == TRUE)
  {
    species_list_sort ();
  }
/*
 *   Print results
 */
  s_h2o->lm = s_h2o->la;
  print_using ();
  print_mix ();
  print_irrev ();
  print_kinetics ();
  print_user_print ();
  print_gas_phase ();
  print_pp_assemblage ();
  print_s_s_assemblage ();
  print_surface ();
  print_exchange ();
  print_initial_solution_isotopes ();
  print_isotope_ratios ();
  print_isotope_alphas ();
  print_totals ();
  print_eh ();
  print_species ();
  print_alkalinity ();
  print_saturation_indices ();
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
punch_all (void)
/* ---------------------------------------------------------------------- */
{
  int i;

#ifndef PHREEQ98		/* if not PHREEQ98 use the standard declaration */
  if (pr.hdf == FALSE && (punch.in == FALSE || pr.punch == FALSE))
    return (OK);
/*
 *   Punch results
 */
  if (state == TRANSPORT || state == PHAST || state == ADVECTION)
  {
    use.kinetics_ptr = kinetics_bsearch (use.n_kinetics_user, &i);
  }
  else if (use.kinetics_in != FALSE)
  {
    use.kinetics_ptr = kinetics_bsearch (-2, &i);
  }
#else /* if PHREEQ98 execute punch_user_graph first, that is, before punch.in and pr.punch is checked */
  if (state == TRANSPORT || state == PHAST || state == ADVECTION)
  {
    use.kinetics_ptr = kinetics_bsearch (use.n_kinetics_user, &i);
  }
  else if (use.kinetics_in != FALSE)
  {
    use.kinetics_ptr = kinetics_bsearch (-2, &i);
  }
  if (pr.user_graph == TRUE)
    punch_user_graph ();

  if (pr.hdf == FALSE && (punch.in == FALSE || pr.punch == FALSE))
    return (OK);
#endif
  punch_identifiers ();
  punch_totals ();
  punch_molalities ();
  punch_activities ();
  punch_pp_assemblage ();
  punch_saturation_indices ();
  punch_gas_phase ();
  punch_kinetics ();
  punch_s_s_assemblage ();
  punch_isotopes ();
  punch_calculate_values ();
  punch_user_punch ();
/*
 *   new line for punch_file
 */
  fpunchf (NULL, "\n");

#if defined(SWIG_SHARED_OBJ)
  EndRow ();
#endif

  return (OK);
}

/* ---------------------------------------------------------------------- */
int
print_diffuse_layer (struct surface_charge *surface_charge_ptr1)
/* ---------------------------------------------------------------------- */
{
/*
 *   Prints total moles of each element in diffuse layer
 *   Remove comment to print total moles of each species
 */
  int i, j, count_g;
  LDBLE mass_water_surface, r, sum_surfs;
  LDBLE molality, moles_excess, moles_surface;

  if (use.surface_ptr == NULL)
    return (OK);
/*
 *   Find position of component in surface charge data
 */
  i = 0;
  for (j = 0; j < count_unknowns; j++)
  {
    if (x[j]->type != SURFACE_CB)
      continue;
    if (x[j]->surface_charge == surface_charge_ptr1)
    {
      break;
    }
    i++;
  }
  if (j >= count_unknowns)
  {
    sprintf (error_string, "In print_diffuse_layer: component not found, %s.",
	     surface_charge_ptr1->name);
    error_msg (error_string, STOP);
  }
/*
 *   Loop through all surface components, calculate each H2O surface (diffuse layer),
 *   H2O aq, and H2O bulk (diffuse layers plus aqueous).
 */

  output_msg (OUTPUT_MESSAGE,
	      "\tWater in diffuse layer: %8.3e kg, %4.1f%% of total DDL-water.\n",
	      (double) surface_charge_ptr1->mass_water,
	      (double) 100 * surface_charge_ptr1->mass_water / mass_water_surfaces_x);
  if (use.surface_ptr->debye_lengths > 0)
  {
    sum_surfs = 0.0;
    for (j = 0; j < count_unknowns; j++)
    {
      if (x[j]->type != SURFACE_CB)
        continue;
      sum_surfs += x[j]->surface_charge->specific_area * x[j]->surface_charge->grams;
    }
    r = 0.002 * mass_water_bulk_x / sum_surfs;
    output_msg (OUTPUT_MESSAGE,
		"\tRadius of total pore:   %8.3e m; of free pore: %8.3e m.\n",
		(double) r, (double) r - use.surface_ptr->thickness);
  }

  if (debug_diffuse_layer == TRUE)
  {
    output_msg (OUTPUT_MESSAGE,
		"\n\t\tDistribution of species in diffuse layer\n\n");
    output_msg (OUTPUT_MESSAGE,
		"\n\tSpecies     \t    Moles   \tMoles excess\t      g\n");
  }
  mass_water_surface = surface_charge_ptr1->mass_water;
  count_elts = 0;
  paren_count = 0;
  for (j = 0; j < count_s_x; j++)
  {
    if (s_x[j]->type > HPLUS)
      continue;
    molality = under (s_x[j]->lm);
    count_g = s_x[j]->diff_layer[i].count_g;
#ifdef SKIP
    moles_excess = mass_water_bulk_x *
/*                      s_x[j]->diff_layer[i].charge->g[count_g].g * molality; */
      surface_charge_ptr1->g[count_g].g * molality;
#endif
    moles_excess =
      mass_water_aq_x * molality * surface_charge_ptr1->g[count_g].g;
    moles_surface = mass_water_surface * molality + moles_excess;
    if (debug_diffuse_layer == TRUE)
    {
      output_msg (OUTPUT_MESSAGE, "\t%-12s\t%12.3e\t%12.3e\t%12.3e\n",
		  s_x[j]->name, moles_surface, moles_excess,
		  s_x[j]->diff_layer[i].charge->g[count_g].g);
    }
/*
 *   Accumulate elements in diffuse layer
 */
    add_elt_list (s_x[j]->next_elt, moles_surface);
  }
/*
        strcpy(token, s_h2o->name);
        ptr = &(token[0]);
        get_elts_in_species (&ptr, mass_water_surface / gfw_water);
 */
  if (count_elts > 0)
  {
    qsort (elt_list, (size_t) count_elts,
	   (size_t) sizeof (struct elt_list), elt_list_compare);
    elt_list_combine ();
  }
/*
 *   Print totals
 */
  if (use.surface_ptr->dl_type != DONNAN_DL)
  {
    output_msg (OUTPUT_MESSAGE,
		"\n\tTotal moles in diffuse layer (excluding water)\n\n");
  }
  else
  {
    output_msg (OUTPUT_MESSAGE,
		"\n\tTotal moles in diffuse layer (excluding water), Donnan calculation \n\n");
  }
  output_msg (OUTPUT_MESSAGE, "\tElement       \t     Moles\n");
  for (j = 0; j < count_elts; j++)
  {
    output_msg (OUTPUT_MESSAGE, "\t%-14s\t%12.4e\n",
		elt_list[j].elt->name, (double) elt_list[j].coef);
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
print_eh (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Prints eh calculated from redox couples
 *   Only calculates eh if two redox states of an element have mass balance
 *   equations.
 */
  int i, j, k, first;
  LDBLE pe, eh;
  struct master *master_ptr0, *master_ptr1;
  char token[MAX_LENGTH];

  if (pr.eh == FALSE || pr.all == FALSE)
    return (OK);

  tk_x = tc_x + 273.15;

  first = TRUE;
  for (i = 0; i < count_master; i++)
  {
    if (master[i]->in != TRUE)
      continue;
    if (master[i]->primary == TRUE)
      continue;
/*
 *   Secondary master species has mass balance equation
 */
    master_ptr0 = master[i]->elt->primary;
    for (k = i + 1; k < count_master; k++)
    {
      if (master[k]->in != TRUE)
	continue;
      master_ptr1 = master[k]->elt->primary;
      if (master_ptr1 != master_ptr0)
	break;
/*
 *  Another secondary master species of same element has mass balance equation
 *  Rewrite equations to calculate pe
 */
      rewrite_master_to_secondary (master[k], master[i]);
      trxn_swap ("e-");
/* debug
                        trxn_print();
 */
/*
 *   Calculate pe, eh
 */
      pe = -k_calc (trxn.logk, tk_x);
      for (j = 1; j < count_trxn; j++)
      {
	pe -= trxn.token[j].s->la * trxn.token[j].coef;
      }
      eh = ((LOG_10 * R_KJ_DEG_MOL * tk_x) / F_KJ_V_EQ) * pe;
/*
 *   Print heading
 */
      if (first == TRUE)
      {
	print_centered ("Redox couples");
	output_msg (OUTPUT_MESSAGE, "\t%-15s%12s%12s\n\n", "Redox couple",
		    "pe", "Eh (volts)");
	first = FALSE;
      }
/*
 *   Print result
 */
      strcpy (token, master[i]->elt->name);
      strcat (token, "/");
      strcat (token, master[k]->elt->name);
      output_msg (OUTPUT_MESSAGE, "\t%-15s%12.4f%12.4f\n", token, (double) pe,
		  (double) eh);
    }
  }
  if (first == FALSE)
    output_msg (OUTPUT_MESSAGE, "\n");
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
print_exchange (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Print moles of each exchange species
 */
  int i;
  struct exchange *exchange_ptr;
  char *name, *name1;
  struct master *master_ptr;
  LDBLE dum, dum2;
/*
 *  Print exchange data
 */
  exchange_ptr = use.exchange_ptr;
  if (exchange_ptr == NULL || pr.exchange == FALSE || pr.all == FALSE)
    return (OK);

  if (state >= REACTION)
  {
    print_centered ("Exchange composition");
  }
/*
 *   Print list of species
 */

  s_h2o->lm = s_h2o->la;
  name = s_hplus->secondary->elt->name;
  for (i = 0; i < count_species_list; i++)
  {
/*
 *   Get name of master species
 */
    if (species_list[i].s->type != EX)
      continue;
    if (species_list[i].master_s->secondary != NULL)
    {
      master_ptr = species_list[i].master_s->secondary;
      name1 = species_list[i].master_s->secondary->elt->name;
    }
    else
    {
      master_ptr = species_list[i].master_s->primary;
      name1 = species_list[i].master_s->primary->elt->name;
    }
/*
 *   Check if new master species, print total molality
 */
    if (name1 != name)
    {
      name = name1;
      output_msg (OUTPUT_MESSAGE, "%-14s%12.3e mol", name,
		  (double) master_ptr->unknown->moles);
      if (master_ptr->unknown->exch_comp->phase_name != NULL)
      {
	output_msg (OUTPUT_MESSAGE, "\t[%g (mol %s)/(mol %s)]",
		    (double) master_ptr->unknown->exch_comp->phase_proportion,
		    master_ptr->unknown->exch_comp->formula,
		    master_ptr->unknown->exch_comp->phase_name);
      }
      else if (master_ptr->unknown->exch_comp->rate_name != NULL)
      {
	output_msg (OUTPUT_MESSAGE,
		    "\t[%g (mol %s)/(mol kinetic reactant %s)]",
		    (double) master_ptr->unknown->exch_comp->phase_proportion,
		    master_ptr->unknown->exch_comp->formula,
		    master_ptr->unknown->exch_comp->rate_name);
      }
      output_msg (OUTPUT_MESSAGE, "\n\n");
      /* Heading for species */
      output_msg (OUTPUT_MESSAGE, "\t%-15s%12s%12s%12s%10s\n", " ", " ",
		  "Equiv-  ", "Equivalent", "Log ");
      output_msg (OUTPUT_MESSAGE, "\t%-15s%12s%12s%12s%10s\n\n", "Species",
		  "Moles  ", "alents  ", "Fraction", "Gamma");
    }
/*
 *   Print species data
 */
/* !!!!! */
    if (master_ptr->total > 1.0e-10)
    {
      if (species_list[i].s->equiv != 0.0)
      {
	dum = fabs (species_list[i].s->equiv) / master_ptr->total;
      }
      else
      {
	if (species_list[i].master_s->z == 0)
	{
	  dum = 1 / master_ptr->total;
	}
	else
	{
	  dum = 1;
	}
      }
      if (species_list[i].master_s->z != 0.0)
      {
	dum2 = fabs (species_list[i].master_s->z);
      }
      else
      {
	dum2 = 1;
      }
      output_msg (OUTPUT_MESSAGE, "\t%-15s%12.3e%12.3e%12.3e%10.3f\n",
		  species_list[i].s->name,
		  (double) species_list[i].s->moles,
		  (double) (species_list[i].s->moles * dum2 *
			    species_list[i].s->equiv),
		  (double) (species_list[i].s->moles * dum /* / dum2 */ ),
		  (double) (species_list[i].s->lg - log10 (dum)));
    }
  }
  output_msg (OUTPUT_MESSAGE, "\n");
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
print_gas_phase (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Prints gas phase composition if present
 */
  int j;
  LDBLE lp, moles, initial_moles, delta_moles;
  struct rxn_token *rxn_ptr;

  if (pr.gas_phase == FALSE || pr.all == FALSE)
    return (OK);
  if (use.gas_phase_ptr == NULL)
    return (OK);
  if (use.gas_phase_ptr->type == PRESSURE)
  {
    if (gas_unknown == NULL)
      return (OK);
    if (gas_unknown->moles < 1e-12)
      return (OK);
    use.gas_phase_ptr->total_moles = gas_unknown->moles;
    use.gas_phase_ptr->volume =
      use.gas_phase_ptr->total_moles * R_LITER_ATM * tk_x /
      use.gas_phase_ptr->total_p;
  }
/*
 *   Print heading
 */
  print_centered ("Gas phase");
  output_msg (OUTPUT_MESSAGE, "\n");
  output_msg (OUTPUT_MESSAGE, "Total pressure: %8.4f   atmospheres\n",
	      (double) use.gas_phase_ptr->total_p);
  output_msg (OUTPUT_MESSAGE, "    Gas volume: %10.2e liters\n",
	      (double) use.gas_phase_ptr->volume);

  output_msg (OUTPUT_MESSAGE, "\n%68s\n%78s\n", "Moles in gas",
	      "----------------------------------");
  output_msg (OUTPUT_MESSAGE, "%-18s%12s%12s%12s%12s%12s\n\n", "Component",
	      "log P", "P", "Initial", "Final", "Delta");

  for (j = 0; j < use.gas_phase_ptr->count_comps; j++)
  {
/*
 *   Calculate partial pressure
 */
    if (use.gas_phase_ptr->comps[j].phase->in == TRUE)
    {
      lp = -use.gas_phase_ptr->comps[j].phase->lk;
      for (rxn_ptr = use.gas_phase_ptr->comps[j].phase->rxn_x->token + 1;
	   rxn_ptr->s != NULL; rxn_ptr++)
      {
	lp += rxn_ptr->s->la * rxn_ptr->coef;
      }
      moles = use.gas_phase_ptr->comps[j].phase->moles_x;
    }
    else
    {
      lp = -99.99;
      moles = 0;
      use.gas_phase_ptr->comps[j].phase->p_soln_x = 0;
    }
/*
 *   Print gas composition 
 */
    if (state != TRANSPORT && state != PHAST)
    {
      initial_moles = use.gas_phase_ptr->comps[j].moles;
      delta_moles = use.gas_phase_ptr->comps[j].phase->moles_x -
	use.gas_phase_ptr->comps[j].moles;
    }
    else
    {
      initial_moles = use.gas_phase_ptr->comps[j].initial_moles;
      delta_moles = use.gas_phase_ptr->comps[j].initial_moles -
	use.gas_phase_ptr->comps[j].moles;
    }
    if (moles <= MIN_TOTAL)
      moles = 0.0;
    if (fabs (delta_moles) <= MIN_TOTAL)
      delta_moles = 0.0;
    output_msg (OUTPUT_MESSAGE, "%-18s%12.2f%12.3e%12.3e%12.3e%12.3e\n",
		use.gas_phase_ptr->comps[j].phase->name,
		(double) lp,
		(double) use.gas_phase_ptr->comps[j].phase->p_soln_x,
		(double) initial_moles, (double) moles, (double) delta_moles);
  }
  output_msg (OUTPUT_MESSAGE, "\n");
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
print_s_s_assemblage (void)
/* ---------------------------------------------------------------------- */
{
  /*
   *   Prints solid solution composition if present
   */
  int i, j;
  LDBLE delta_moles;
  LDBLE nb, nc, xb, xb1, xb2, xb1moles, xb2moles;
  struct s_s *s_s_ptr;

  if (pr.s_s_assemblage == FALSE || pr.all == FALSE)
    return (OK);
  if (use.s_s_assemblage_ptr == NULL)
    return (OK);
  /*
   *   Print heading
   */
  print_centered ("Solid solutions");
  output_msg (OUTPUT_MESSAGE, "\n");
  output_msg (OUTPUT_MESSAGE, "%-15s  %22s  %11s  %11s  %11s\n\n",
	      "Solid solution", "Component", "Moles", "Delta moles",
	      "Mole fract");
  /*
   *   Print solid solutions
   */
  for (j = 0; j < use.s_s_assemblage_ptr->count_s_s; j++)
  {
    if (use.s_s_assemblage_ptr->s_s[j].s_s_in == TRUE)
    {
      /* solid solution name, moles */
      output_msg (OUTPUT_MESSAGE, "%-15s  %22s  %11.2e\n",
		  use.s_s_assemblage_ptr->s_s[j].name, "  ",
		  (double) use.s_s_assemblage_ptr->s_s[j].total_moles);
      /* component name, moles, delta moles, mole fraction */
      for (i = 0; i < use.s_s_assemblage_ptr->s_s[j].count_comps; i++)
      {
	if (state != TRANSPORT && state != PHAST)
	{
	  delta_moles = use.s_s_assemblage_ptr->s_s[j].comps[i].moles -
	    use.s_s_assemblage_ptr->s_s[j].comps[i].initial_moles -
	    use.s_s_assemblage_ptr->s_s[j].comps[i].delta;
	}
	else
	{
	  delta_moles = use.s_s_assemblage_ptr->s_s[j].comps[i].moles -
	    use.s_s_assemblage_ptr->s_s[j].comps[i].init_moles;
	}
	output_msg (OUTPUT_MESSAGE, "%15s  %22s  %11.2e  %11.2e  %11.2e\n",
		    " ", use.s_s_assemblage_ptr->s_s[j].comps[i].name,
		    (double) use.s_s_assemblage_ptr->s_s[j].comps[i].moles,
		    (double) delta_moles,
		    (double) (use.s_s_assemblage_ptr->s_s[j].comps[i].moles /
			      use.s_s_assemblage_ptr->s_s[j].total_moles));
      }
      s_s_ptr = &(use.s_s_assemblage_ptr->s_s[j]);
      if (s_s_ptr->miscibility == TRUE)
      {
	nc = s_s_ptr->comps[0].moles;
	nb = s_s_ptr->comps[1].moles;
	xb = nb / (nb + nc);
	xb1 = s_s_ptr->xb1;
	xb2 = s_s_ptr->xb2;

	if (xb > xb1 && xb < xb2)
	{
	  xb2moles = (xb1 - 1) / xb1 * nb + nc;
	  xb2moles = xb2moles / ((xb1 - 1) / xb1 * xb2 + (1 - xb2));
	  xb1moles = (nb - xb2moles * xb2) / xb1;
	  output_msg (OUTPUT_MESSAGE,
		      "\n%14s  Solid solution is in miscibility gap\n", " ");
	  output_msg (OUTPUT_MESSAGE, "%14s  End members in pct of %s\n\n",
		      " ", s_s_ptr->comps[1].name);
	  output_msg (OUTPUT_MESSAGE, "%22s  %11g pct  %11.2e\n", " ",
		      (double) xb1, (double) xb1moles);
	  output_msg (OUTPUT_MESSAGE, "%22s  %11g pct  %11.2e\n", " ",
		      (double) xb2, (double) xb2moles);
	}
      }
    }
    else
    {
      /* solid solution name, moles */
      output_msg (OUTPUT_MESSAGE, "%-15s  %22s  %11.2e\n",
		  use.s_s_assemblage_ptr->s_s[j].name, "  ", (double) 0.0);
      /* component name, moles, delta moles, mole fraction */
      for (i = 0; i < use.s_s_assemblage_ptr->s_s[j].count_comps; i++)
      {
	if (state != TRANSPORT && state != PHAST)
	{
	  delta_moles = use.s_s_assemblage_ptr->s_s[j].comps[i].moles -
	    use.s_s_assemblage_ptr->s_s[j].comps[i].initial_moles -
	    use.s_s_assemblage_ptr->s_s[j].comps[i].delta;
	}
	else
	{
	  delta_moles = use.s_s_assemblage_ptr->s_s[j].comps[i].moles -
	    use.s_s_assemblage_ptr->s_s[j].comps[i].init_moles;
	}
	output_msg (OUTPUT_MESSAGE, "%15s  %22s  %11.2e  %11.2e  %11.2e\n",
		    " ", use.s_s_assemblage_ptr->s_s[j].comps[i].name,
		    (double) 0, (double) delta_moles, (double) 0);
      }
    }
  }
  output_msg (OUTPUT_MESSAGE, "\n");
  return (OK);
}

/* ---------------------------------------------------------------------- */
static int
print_irrev (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   prints irreversible reaction as defined and as
 *   relative moles of each element and total amount
 *   of reaction
 */
  int i;
  struct irrev *irrev_ptr;

  if (pr.use == FALSE || pr.all == FALSE)
    return (OK);
  if (state < REACTION || use.irrev_in == FALSE)
    return (OK);
  if (state == TRANSPORT && transport_step == 0)
    return (OK);
  irrev_ptr = use.irrev_ptr;
/*
 *  Print amount of reaction
 */
  output_msg (OUTPUT_MESSAGE, "Reaction %d.\t%s\n\n", use.n_irrev_user,
	      irrev_ptr->description);
  output_msg (OUTPUT_MESSAGE,
	      "\t%11.3e moles of the following reaction have been added:\n\n",
	      (double) step_x);
/*
 *  Print reaction
 */
  output_msg (OUTPUT_MESSAGE, "\t%-15s%10s\n", " ", "Relative");
  output_msg (OUTPUT_MESSAGE, "\t%-15s%10s\n\n", "Reactant", "moles");
  for (i = 0; i < irrev_ptr->count_list; i++)
  {
    output_msg (OUTPUT_MESSAGE, "\t%-15s%13.5f\n", irrev_ptr->list[i].name,
		(double) irrev_ptr->list[i].coef);
  }
  output_msg (OUTPUT_MESSAGE, "\n");
/*
 *   Debug
 */

  output_msg (OUTPUT_MESSAGE, "\t%-15s%10s\n", " ", "Relative");
  output_msg (OUTPUT_MESSAGE, "\t%-15s%10s\n", "Element", "moles");
  for (i = 0; irrev_ptr->elts[i].elt != NULL; i++)
  {
    output_msg (OUTPUT_MESSAGE, "\t%-15s%13.5f\n",
		irrev_ptr->elts[i].elt->name,
		(double) irrev_ptr->elts[i].coef);
  }
  output_msg (OUTPUT_MESSAGE, "\n");
  return (OK);
}

/* ---------------------------------------------------------------------- */
static int
print_kinetics (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   prints kinetic reaction,
 *   should be called only on final kinetic step
 */
  int i, j;
  LDBLE sim_time;
  struct kinetics *kinetics_ptr;
  if (pr.kinetics == FALSE || pr.all == FALSE)
    return (OK);
  if (state < REACTION)
    return (OK);
  kinetics_ptr = NULL;
  if (use.kinetics_in == TRUE)
  {
    if (state == TRANSPORT || state == PHAST || state == ADVECTION)
    {
      kinetics_ptr = kinetics_bsearch (use.n_kinetics_user, &i);
    }
    else
    {
      kinetics_ptr = kinetics_bsearch (-2, &i);
    }
  }
  if (kinetics_ptr == NULL)
    return (OK);
/*
 *   determine time step
 */
  if (state == TRANSPORT || state == PHAST)
  {
    kin_time_x = timest;
  }
  else if (state == ADVECTION)
  {
    kin_time_x = advection_kin_time;
  }
  sim_time = 0.;
  if (incremental_reactions == TRUE)
  {
    if (kinetics_ptr->count_steps > 0)
    {
      for (i = 0; i < reaction_step; i++)
      {
	if (i < kinetics_ptr->count_steps)
	{
	  sim_time += kinetics_ptr->steps[i];
	}
	else
	{
	  sim_time += kinetics_ptr->steps[kinetics_ptr->count_steps - 1];
	}
      }
    }
    else if (kinetics_ptr->count_steps < 0)
    {
      if (reaction_step > -kinetics_ptr->count_steps)
      {
	sim_time = kinetics_ptr->steps[0];
      }
      else
      {
	sim_time =
	  reaction_step * kinetics_ptr->steps[0] /
	  ((LDBLE) (-kinetics_ptr->count_steps));
      }
    }
  }
/*
 *  Print amount of reaction
 */
  if (phast == FALSE)
  {
    output_msg (OUTPUT_MESSAGE, "Kinetics %d.\t%s\n\n", use.n_kinetics_user,
		kinetics_ptr->description);
  }
  else
  {
    output_msg (OUTPUT_MESSAGE, "Kinetics.\n\n");
  }
/*
 *  Print reaction
 */
  if (state == TRANSPORT)
  {
    output_msg (OUTPUT_MESSAGE, "\tTime:      %g seconds\n",
		(double) (initial_total_time + transport_step * timest));
    output_msg (OUTPUT_MESSAGE, "\tTime step: %g seconds\n\n",
		(double) kin_time_x);
  }
  else if (state == ADVECTION)
  {
    output_msg (OUTPUT_MESSAGE, "\tTime:      %g seconds\n",
		(double) (initial_total_time +
			  advection_step * advection_kin_time));
    output_msg (OUTPUT_MESSAGE, "\tTime step: %g seconds\n\n",
		(double) kin_time_x);
  }
  else if (state == PHAST)
  {
    output_msg (OUTPUT_MESSAGE, "\tTime:      %g seconds\n",
		(double) rate_sim_time_end);
    output_msg (OUTPUT_MESSAGE, "\tTime step: %g seconds\n\n",
		(double) kin_time_x);
  }
  else if (state == REACTION)
  {
    if (incremental_reactions == FALSE)
    {
      output_msg (OUTPUT_MESSAGE, "\tTime step: %g seconds\n\n",
		  (double) kin_time_x);
    }
    else
    {
      output_msg (OUTPUT_MESSAGE,
		  "\tTime step: %g seconds  (Incremented time: %g seconds)\n\n",
		  (double) kin_time_x, (double) sim_time);
    }
  }
  output_msg (OUTPUT_MESSAGE, "\t%-15s%12s%12s   %-15s%12s\n\n", "Rate name",
	      "Delta Moles", "Total Moles", "Reactant", "Coefficient");
  for (i = 0; i < kinetics_ptr->count_comps; i++)
  {
    if (state != TRANSPORT && state != PHAST)
    {
      output_msg (OUTPUT_MESSAGE, "\t%-15s%12.3e%12.3e",
		  kinetics_ptr->comps[i].rate_name,
		  (double) -kinetics_ptr->comps[i].moles,
		  (double) kinetics_ptr->comps[i].m);
    }
    else
    {
      output_msg (OUTPUT_MESSAGE, "\t%-15s%12.3e%12.3e",
		  kinetics_ptr->comps[i].rate_name,
		  (double) (kinetics_ptr->comps[i].m -
			    kinetics_ptr->comps[i].initial_moles),
		  (double) kinetics_ptr->comps[i].m);
    }
    for (j = 0; j < kinetics_ptr->comps[i].count_list; j++)
    {
      if (j == 0)
      {
	output_msg (OUTPUT_MESSAGE, "   %-15s%12g\n",
		    kinetics_ptr->comps[i].list[j].name,
		    (double) kinetics_ptr->comps[i].list[j].coef);
      }
      else
      {
	output_msg (OUTPUT_MESSAGE, "\t%39s   %-15s%12g\n", " ",
		    kinetics_ptr->comps[i].list[j].name,
		    (double) kinetics_ptr->comps[i].list[j].coef);
      }
    }
  }
  output_msg (OUTPUT_MESSAGE, "\n");
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
print_master_reactions (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Debugging print routine to test primary and secondary reactions
 */
  int i;
  struct rxn_token *next_token;

  for (i = 0; i < count_master; i++)
  {
    output_msg (OUTPUT_MESSAGE, "%s\t%s\n\tPrimary reaction\n",
		master[i]->elt->name, master[i]->s->name);
    next_token = master[i]->rxn_primary->token;
    for (; next_token->s != NULL; next_token++)
    {
      output_msg (OUTPUT_MESSAGE, "\t\t%s\t%f\n", next_token->s->name,
		  (double) next_token->coef);
    }
    output_msg (OUTPUT_MESSAGE, "\n\tSecondary reaction:\n");
    if (master[i]->rxn_secondary != NULL)
    {
      next_token = master[i]->rxn_secondary->token;
      for (; next_token->s != NULL; next_token++)
      {
	output_msg (OUTPUT_MESSAGE, "\t\t%s\t%f\n", next_token->s->name,
		    (double) next_token->coef);
      }
    }
    output_msg (OUTPUT_MESSAGE, "\n\tRedox reaction:\n");
    if (*(master[i]->pe_rxn) != NULL)
    {
      next_token = (*(master[i]->pe_rxn))->token;
      for (; next_token->s != NULL; next_token++)
      {
	output_msg (OUTPUT_MESSAGE, "\t\t%s\t%f\n", next_token->s->name,
		    (double) next_token->coef);
      }
    }
    output_msg (OUTPUT_MESSAGE, "\n");
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
static int
print_mix (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   prints definition of mixing, solution number and multiplier
 */
  int i, n;
  struct mix *mix_ptr;
  struct solution *solution_ptr;

  if (pr.use == FALSE || pr.all == FALSE)
    return (OK);
  if (use.mix_in == FALSE || state < REACTION)
    return (OK);
  if (state == TRANSPORT)
  {
    mix_ptr = mix_bsearch (use.n_mix_user, &i);
  }
  else
  {
    mix_ptr = mix_bsearch (use.n_mix_user_orig, &i);
  }
  if (mix_ptr == NULL)
  {
    mix_ptr = use.mix_ptr;
  }
/*
 *  Print mixture data
 */
  if (mix_ptr == NULL)
  {
    return (OK);

  }
  if (state == TRANSPORT)
  {
    output_msg (OUTPUT_MESSAGE, "Mixture %d.\t%s\n\n", use.n_mix_user,
		mix_ptr->description);
  }
  else
  {
    output_msg (OUTPUT_MESSAGE, "Mixture %d.\t%s\n\n", mix_ptr->n_user,
		mix_ptr->description);
  }
  for (i = 0; i < mix_ptr->count_comps; i++)
  {
    solution_ptr = solution_bsearch (mix_ptr->comps[i].n_solution, &n, TRUE);
    if (solution_ptr == NULL)
    {
      input_error++;
      return (ERROR);
    }
    output_msg (OUTPUT_MESSAGE, "\t%11.3e Solution %d\t%-55s\n",
		(double) mix_ptr->comps[i].fraction,
		mix_ptr->comps[i].n_solution, solution[n]->description);
  }
  output_msg (OUTPUT_MESSAGE, "\n");
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
print_reaction (struct reaction *rxn_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Debugging print of individual chemical reactions for
 *   species or phases
 */
  int j;
  struct rxn_token *next_token;

  if (pr.use == FALSE || pr.all == FALSE)
    return (OK);

  output_msg (OUTPUT_MESSAGE, "%s\t\n", rxn_ptr->token[0].s->name);
  output_msg (OUTPUT_MESSAGE, "\n\tlog k:\n");
  for (j = 0; j < 7; j++)
  {
    output_msg (OUTPUT_MESSAGE, "\t%f", (double) rxn_ptr->logk[j]);
  }
  output_msg (OUTPUT_MESSAGE, "\n\nReaction:\n");
  for (next_token = rxn_ptr->token; next_token->s != NULL; next_token++)
  {
    output_msg (OUTPUT_MESSAGE, "\t\t%s\t%f\n", next_token->s->name,
		(double) next_token->coef);
  }
  output_msg (OUTPUT_MESSAGE, "\n");
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
print_saturation_indices (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Prints saturation indices of all applicable pure_phases
 */
  int i;
  LDBLE si, iap;
  LDBLE lk;
  LDBLE la_eminus;
  struct rxn_token *rxn_ptr;

  if (pr.saturation_indices == FALSE || pr.all == FALSE)
    return (OK);
  if (state == INITIAL_SOLUTION)
  {
    iap = 0;
    for (rxn_ptr = pe_x[default_pe_x].rxn->token + 1; rxn_ptr->s != NULL;
	 rxn_ptr++)
    {
      iap += rxn_ptr->coef * rxn_ptr->s->la;
      /* fprintf(output,"\t%s\t%f\t%f\n", rxn_ptr->s->name, rxn_ptr->coef, rxn_ptr->s->la ); */
    }
    lk = k_calc (pe_x[default_pe_x].rxn->logk, tk_x);
    la_eminus = lk + iap;
    /* fprintf(output,"\t%s\t%f\n", "pe", si ); */
  }
  else
  {
    la_eminus = s_eminus->la;
  }
/*
 *   Print heading
 */
  print_centered ("Saturation indices");
  output_msg (OUTPUT_MESSAGE, "\t%-15s%7s%8s%8s\n\n", "Phase", "SI",
	      "log IAP", "log KT");

  for (i = 0; i < count_phases; i++)
  {
    if (phases[i]->in == FALSE || phases[i]->type != SOLID)
      continue;
/*
 *   Print saturation index
 */
#ifdef SKIP
    iap = 0.0;
    for (rxn_ptr = phases[i]->rxn_x->token + 1; rxn_ptr->s != NULL; rxn_ptr++)
    {
      iap += rxn_ptr->s->la * rxn_ptr->coef;
      /* fprintf(output,"\t%s\t%f\t%e\n", rxn_ptr->s->name, rxn_ptr->s->la, rxn_ptr->coef); */
    }
    si = -phases[i]->lk + iap;
    output_msg (OUTPUT_MESSAGE, "\t%-15s%7.2f%8.2f%8.2f  %s\n",
		phases[i]->name, (double) si, (double) iap,
		(double) phases[i]->lk, phases[i]->formula);
    si1 = si;
#endif

    /* pe_value_x = pe_x[default_pe_x].rxn; */
    lk = k_calc (phases[i]->rxn->logk, tk_x);
    iap = 0.0;
    for (rxn_ptr = phases[i]->rxn->token + 1; rxn_ptr->s != NULL; rxn_ptr++)
    {
      /* iap += rxn_ptr->s->la * rxn_ptr->coef; */
      if (rxn_ptr->s != s_eminus)
      {
	iap += (rxn_ptr->s->lm + rxn_ptr->s->lg) * rxn_ptr->coef;
      }
      else
      {
	iap += la_eminus * rxn_ptr->coef;
      }
      /* output_msg(OUTPUT_MESSAGE,"\t%s\t%f\t%e\t%e\t%e\n", rxn_ptr->s->name, rxn_ptr->s->la, rxn_ptr->coef, rxn_ptr->s->lm, rxn_ptr->s->lg ); */
    }
    si = -lk + iap;
    output_msg (OUTPUT_MESSAGE, "\t%-15s%7.2f%8.2f%8.2f  %s\n",
		phases[i]->name, (double) si, (double) iap, (double) lk,
		phases[i]->formula);

#ifdef SKIP
    if (fabs (si1 - si) > convergence_tolerance)
    {
      sprintf (error_string, "SI1 != SI, %s SI1 %e, SI %e, diff %e\n",
	       phases[i]->name, si1, si, si1 - si);
      error_msg (error_string, CONTINUE);
    }
#endif
  }
  output_msg (OUTPUT_MESSAGE, "\n");
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
print_pp_assemblage (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Prints saturation indices and masses of pure_phases in pp_assemblage
 */
  int j, k;
  LDBLE si, iap, lk;
  char token[MAX_LENGTH];
  struct rxn_token *rxn_ptr;
  struct phase *phase_ptr;

  if (pr.pp_assemblage == FALSE || pr.all == FALSE)
    return (OK);
  if (pure_phase_unknown == NULL)
    return (OK);
/*
 *   Print heading
 */
  print_centered ("Phase assemblage");
  output_msg (OUTPUT_MESSAGE, "%73s\n", "Moles in assemblage");
  output_msg (OUTPUT_MESSAGE, "%-18s%7s%8s%8s", "Phase", "SI", "log IAP",
	      "log KT");
  output_msg (OUTPUT_MESSAGE, " %12s%12s%12s", " Initial", " Final",
	      " Delta");
  output_msg (OUTPUT_MESSAGE, "\n\n");

  for (j = 0; j < count_unknowns; j++)
  {
    if (x[j]->type != PP)
      continue;
/*
 *   Print saturation index
 */
    iap = 0.0;
    phase_ptr = x[j]->phase;
    if (x[j]->phase->rxn_x == NULL || phase_ptr->in == FALSE)
    {
      output_msg (OUTPUT_MESSAGE, "%-18s%23s", x[j]->phase->name,
		  "Element not present.");
    }
    else
    {
      phase_ptr = x[j]->phase;
      lk = k_calc (phase_ptr->rxn->logk, tk_x);
      for (rxn_ptr = phase_ptr->rxn->token + 1; rxn_ptr->s != NULL; rxn_ptr++)
      {
	if (rxn_ptr->s != s_eminus)
	{
	  iap += (rxn_ptr->s->lm + rxn_ptr->s->lg) * rxn_ptr->coef;
	}
	else
	{
	  iap += s_eminus->la * rxn_ptr->coef;
	}
      }
      si = -lk + iap;
      /*
         for (rxn_ptr = x[j]->phase->rxn_x->token + 1; rxn_ptr->s != NULL; rxn_ptr++) {
         iap += rxn_ptr->s->la * rxn_ptr->coef;
         }
         si = -x[j]->phase->lk + iap;
         output_msg(OUTPUT_MESSAGE,"\t%-15s%7.2f%8.2f%8.2f", x[j]->phase->name, (double) si, (double) iap, (double) x[j]->phase->lk);
       */
      output_msg (OUTPUT_MESSAGE, "%-18s%7.2f%8.2f%8.2f", x[j]->phase->name,
		  (double) si, (double) iap, (double) lk);
    }
/*
 *   Print pure phase assemblage data
 */
    if (x[j]->moles < 0.0)
      x[j]->moles = 0.0;
    if (state != TRANSPORT && state != PHAST)
    {
      sprintf (token, " %11.3e %11.3e %11.3e",
	       (double) (x[j]->pure_phase->moles + x[j]->pure_phase->delta),
	       (double) x[j]->moles,
	       (double) (x[j]->moles - x[j]->pure_phase->moles -
			 x[j]->pure_phase->delta));
    }
    else
    {
      sprintf (token, " %11.3e %11.3e %11.3e",
	       (double) x[j]->pure_phase->initial_moles,
	       (double) x[j]->moles,
	       (double) (x[j]->moles - x[j]->pure_phase->initial_moles));
    }
    if (x[j]->moles <= 0.0)
    {
      for (k = 0; k < 11; k++)
      {
	token[12 + k] = ' ';
      }
    }
    if (x[j]->pure_phase->add_formula == NULL)
    {
      output_msg (OUTPUT_MESSAGE, "%37s\n", token);
    }
    else
    {
      output_msg (OUTPUT_MESSAGE, "\n         %-18s%-15s%36s\n",
		  x[j]->pure_phase->add_formula, " is reactant", token);
    }
  }
  output_msg (OUTPUT_MESSAGE, "\n");
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
print_species (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Prints description of solution, uses array species_list for
 *   order of aqueous species.
 */
  int i;
  char *name, *name1;
  struct master *master_ptr;
  LDBLE min;
  LDBLE lm;

  if (pr.species == FALSE || pr.all == FALSE)
    return (OK);
  min = -1000;
  print_centered ("Distribution of species");
/*
 *   Heading for species
 */
  if (pitzer_model == TRUE)
  {
    if (ICON == TRUE)
    {
      output_msg (OUTPUT_MESSAGE, "%67s%11s\n", "MacInnes", "MacInnes");
      output_msg (OUTPUT_MESSAGE, "\t%-15s%12s%12s%10s%10s%10s\n", " ", " ",
		  "MacInnes", "Log   ", "Log   ", "Log ");
    }
    else
    {
      output_msg (OUTPUT_MESSAGE, "%67s%11s\n", "Unscaled", "Unscaled");
      output_msg (OUTPUT_MESSAGE, "\t%-15s%12s%12s%10s%10s%10s\n", " ", " ",
		  "Unscaled", "Log   ", "Log   ", "Log ");
    }
  }
  else
  {
    output_msg (OUTPUT_MESSAGE, "   %-20s%12s%12s%10s%10s%10s\n", " ", " ",
		" ", "Log   ", "Log   ", "Log ");
  }
  output_msg (OUTPUT_MESSAGE, "   %-20s%12s%12s%10s%10s%10s\n\n", "Species",
	      "Molality", "Activity", "Molality", "Activity", "Gamma");
/*
 *   Print list of species
 */
  s_h2o->lm = s_h2o->la;
  name = s_hplus->secondary->elt->name;
  for (i = 0; i < count_species_list; i++)
  {
/*
 *   Get name of master species
 */
    if (species_list[i].s->type == EX)
      continue;
    if (species_list[i].s->type == SURF)
      continue;
    if (species_list[i].master_s->secondary != NULL)
    {
      master_ptr = species_list[i].master_s->secondary;
      name1 = species_list[i].master_s->secondary->elt->name;
    }
    else
    {
      master_ptr = species_list[i].master_s->primary;
      name1 = species_list[i].master_s->primary->elt->name;
    }
/*
 *   Check if new master species, print total molality
 */
    if (name1 != name)
    {
      name = name1;
      output_msg (OUTPUT_MESSAGE, "%-14s%12.3e\n", name,
		  (double) (master_ptr->total / mass_water_aq_x));
      min = censor * master_ptr->total / mass_water_aq_x;
      if (min > 0)
      {
	min = log10 (min);
      }
      else
      {
	min = -1000.;
      }
    }
/*
 *   Print species data
 */
    if (species_list[i].s->lm > min)
    {
      if (species_list[i].s == s_h2o)
      {
	lm = log10 (s_h2o->moles / mass_water_aq_x);
      }
      else
      {
	lm = species_list[i].s->lm;
      }
      output_msg (OUTPUT_MESSAGE, "   %-20s%12.3e%12.3e%10.3f%10.3f%10.3f\n",
		  species_list[i].s->name,
		  (double) ((species_list[i].s->moles) / mass_water_aq_x),
		  (double) under (species_list[i].s->lm +
				  species_list[i].s->lg), (double) lm,
		  (double) (species_list[i].s->lm + species_list[i].s->lg),
		  (double) species_list[i].s->lg);
    }
  }
  output_msg (OUTPUT_MESSAGE, "\n");
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
print_surface (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Prints description of surface, including charge and potential,
 *   grams and specific area, moles of each species on surface sites,
 *   and description of diffuse layer if applicable.
 */
  int i, j, k;
  struct surface *surface_ptr;
  char name[MAX_LENGTH], token[MAX_LENGTH];
  struct master *master_ptr;
  LDBLE molfrac, charge;
  char *ptr;
/*
 *  Print surface speciation
 */
  surface_ptr = use.surface_ptr;
  if (surface_ptr == NULL || pr.surface == FALSE || pr.all == FALSE)
    return (OK);
  if (surface_ptr->type == CD_MUSIC)
    return (print_surface_cd_music ());

  if (state >= REACTION)
  {
    print_centered ("Surface composition");
  }
/*
 *   Print list of species
 */

  s_h2o->lm = s_h2o->la;
  for (j = 0; j < count_unknowns; j++)
  {
    /*if (use.surface_ptr->edl == TRUE) { */
    if (use.surface_ptr->type == DDL)
    {
      if (x[j]->type != SURFACE_CB)
	continue;
      strcpy (name, x[j]->master[0]->elt->name);
      replace ("_psi", "", name);
    }
    else
    {
      if (x[j]->type != SURFACE)
	continue;
      strcpy (token, x[j]->master[0]->elt->name);
      replace ("_", " ", token);
      ptr = token;
      copy_token (name, &ptr, &k);
    }
    output_msg (OUTPUT_MESSAGE, "%-14s\n", name);
/*
 *   Description of surface
 */
    if (dl_type_x != NO_DL)
    {
      output_msg (OUTPUT_MESSAGE,
		  "\t%11.3e  Surface + diffuse layer charge, eq\n",
		  (double) x[j]->f);
    }
    /*if (use.surface_ptr->edl == TRUE && diffuse_layer_x == FALSE) { */
    if (use.surface_ptr->type == DDL && dl_type_x == NO_DL)
    {
      charge = x[j]->f;
    }
    else
    {
      charge = calc_surface_charge (name);
    }
    output_msg (OUTPUT_MESSAGE, "\t%11.3e  Surface charge, eq\n",
		(double) charge);
    if (x[j]->type == SURFACE_CB)
    {
      if ((x[j]->surface_charge->specific_area *
	   x[j]->surface_charge->grams) > 0)
      {
	output_msg (OUTPUT_MESSAGE, "\t%11.3e  sigma, C/m**2\n",
		    (double) (charge * F_C_MOL /
			      (x[j]->surface_charge->specific_area *
			       x[j]->surface_charge->grams)));
      }
      else
      {
	output_msg (OUTPUT_MESSAGE, "\tundefined  sigma, C/m**2\n");
      }
      output_msg (OUTPUT_MESSAGE, "\t%11.3e  psi, V\n",
		  (double) (x[j]->master[0]->s->la * 2 * R_KJ_DEG_MOL *
			    tk_x * LOG_10 / F_KJ_V_EQ));
      output_msg (OUTPUT_MESSAGE, "\t%11.3e  -F*psi/RT\n",
		  (double) (x[j]->master[0]->s->la * (-2) * LOG_10));
      output_msg (OUTPUT_MESSAGE, "\t%11.3e  exp(-F*psi/RT)\n",
		  exp (x[j]->master[0]->s->la * (-2) * LOG_10));
      if (x[j]->surface_comp->phase_name != NULL)
      {
	output_msg (OUTPUT_MESSAGE, "\t%11.3e  specific area, m**2/mol %s\n",
		    (double) x[j]->surface_charge->specific_area,
		    x[j]->surface_comp->phase_name);
	output_msg (OUTPUT_MESSAGE,
		    "\t%11.3e  m**2 for %11.3e moles of %s\n\n",
		    (double) (x[j]->surface_charge->grams *
			      x[j]->surface_charge->specific_area),
		    (double) x[j]->surface_charge->grams,
		    x[j]->surface_comp->phase_name);
      }
      else if (x[j]->surface_comp->rate_name != NULL)
      {
	output_msg (OUTPUT_MESSAGE, "\t%11.3e  specific area, m**2/mol %s\n",
		    (double) x[j]->surface_charge->specific_area,
		    x[j]->surface_comp->rate_name);
	output_msg (OUTPUT_MESSAGE,
		    "\t%11.3e  m**2 for %11.3e moles of %s\n\n",
		    (double) (x[j]->surface_charge->grams *
			      x[j]->surface_charge->specific_area),
		    (double) x[j]->surface_charge->grams,
		    x[j]->surface_comp->rate_name);
      }
      else
      {
	output_msg (OUTPUT_MESSAGE, "\t%11.3e  specific area, m**2/g\n",
		    (double) x[j]->surface_charge->specific_area);
	output_msg (OUTPUT_MESSAGE, "\t%11.3e  m**2 for %11.3e g\n\n",
		    (double) (x[j]->surface_charge->specific_area *
			      x[j]->surface_charge->grams),
		    (double) x[j]->surface_charge->grams);
      }
      if (dl_type_x != NO_DL)
	print_diffuse_layer (x[j]->surface_charge);
      output_msg (OUTPUT_MESSAGE, "\n");
/*
 *   Heading for species
 */
      for (k = j - 1; k < count_unknowns; k++)
      {
	if (x[k]->type != SURFACE)
	  continue;
	if (x[j] != x[k]->potential_unknown)
	  continue;
	master_ptr = x[k]->master[0];
	output_msg (OUTPUT_MESSAGE, "%-14s\n", x[k]->master[0]->elt->name);
	output_msg (OUTPUT_MESSAGE, "\t%11.3e  moles", (double) x[k]->moles);
	if (x[k]->surface_comp->phase_name != NULL)
	{
	  output_msg (OUTPUT_MESSAGE, "\t[%g mol/(mol %s)]\n",
		      (double) x[k]->surface_comp->phase_proportion,
		      x[k]->surface_comp->phase_name);
	}
	else if (x[k]->surface_comp->rate_name != NULL)
	{
	  output_msg (OUTPUT_MESSAGE,
		      "\t[%g mol/(mol kinetic reactant %s)]\n",
		      (double) x[k]->surface_comp->phase_proportion,
		      x[k]->surface_comp->rate_name);
	}
	else
	{
	  output_msg (OUTPUT_MESSAGE, "\n");
	}
	output_msg (OUTPUT_MESSAGE, "\t%-15s%12s%12s%12s%12s\n", " ", " ",
		    "Mole", " ", "Log");
	output_msg (OUTPUT_MESSAGE, "\t%-15s%12s%12s%12s%12s\n\n", "Species",
		    "Moles", "Fraction", "Molality", "Molality");
	for (i = 0; i < count_species_list; i++)
	{
	  if (species_list[i].master_s != master_ptr->s)
	    continue;
/*
 *   Print species data
 */
	  if (x[k]->moles >= MIN_RELATED_SURFACE)
	  {
	    molfrac =
	      (LDBLE) (species_list[i].s->moles) / x[k]->moles *
	      species_list[i].s->equiv;
	  }
	  else
	  {
	    molfrac = 0.0;
	  }
	  output_msg (OUTPUT_MESSAGE, "\t%-15s%12.3e%12.3f%12.3e%12.3f\n",
		      species_list[i].s->name,
		      (double) species_list[i].s->moles,
		      (double) molfrac,
		      (double) (species_list[i].s->moles / mass_water_aq_x),
		      log10 (species_list[i].s->moles / mass_water_aq_x));
	}
	output_msg (OUTPUT_MESSAGE, "\n");
      }
    }
    else
    {
      k = j;
      master_ptr = x[k]->master[0];
      output_msg (OUTPUT_MESSAGE, "%-14s\n", x[k]->master[0]->elt->name);
      output_msg (OUTPUT_MESSAGE, "\t%11.3e  moles\n", (double) x[k]->moles);
      output_msg (OUTPUT_MESSAGE, "\t%-15s%12s%12s%12s%12s\n", " ", " ",
		  "Mole", " ", "Log");
      output_msg (OUTPUT_MESSAGE, "\t%-15s%12s%12s%12s%12s\n\n", "Species",
		  "Moles", "Fraction", "Molality", "Molality");
      for (i = 0; i < count_species_list; i++)
      {
	if (species_list[i].master_s != master_ptr->s)
	  continue;
/*
 *   Print species data
 */
	if (x[k]->moles >= MIN_RELATED_SURFACE)
	{
	  molfrac =
	    (double) (species_list[i].s->moles) / x[k]->moles *
	    species_list[i].s->equiv;
	}
	else
	{
	  molfrac = 0.0;
	}
	output_msg (OUTPUT_MESSAGE, "\t%-15s%12.3e%12.3f%12.3e%12.3f\n",
		    species_list[i].s->name,
		    (double) species_list[i].s->moles,
		    (double) molfrac,
		    (double) (species_list[i].s->moles / mass_water_aq_x),
		    log10 (species_list[i].s->moles / mass_water_aq_x));
      }
      output_msg (OUTPUT_MESSAGE, "\n");
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
print_surface_cd_music (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Prints description of cd music surfaces, including charge and potential,
 *   grams and specific area, moles of each species on surface sites,
 *   and description of diffuse layer if applicable.
 */
  int i, j, k;
  struct surface *surface_ptr;
  char name[MAX_LENGTH];
  struct master *master_ptr, *master_ptr0, *master_ptr1, *master_ptr2;
  struct unknown *unknown_ptr0, *unknown_ptr1, *unknown_ptr2;
  LDBLE molfrac, charge0, charge1, charge2, sum;


/*
 *  Print surface speciation
 */
  surface_ptr = use.surface_ptr;
  if (surface_ptr == NULL || pr.surface == FALSE || pr.all == FALSE)
    return (OK);

  if (state >= REACTION)
  {
    print_centered ("Surface composition");
  }
/*
 *   Print list of species
 */

  s_h2o->lm = s_h2o->la;
  for (j = 0; j < count_unknowns; j++)
  {
    if (x[j]->type != SURFACE_CB)
      continue;
    strcpy (name, x[j]->master[0]->elt->name);
    replace ("_psi", "", name);
    output_msg (OUTPUT_MESSAGE, "%-14s\n", name);
/*
 *   Description of surface
 */
    if (dl_type_x != NO_DL)
    {
      output_msg (OUTPUT_MESSAGE,
		  "\t%11.3e  Surface + diffuse layer charge, eq\n\n",
		  (double) x[j + 2]->f + (x[j]->surface_charge->sigma0 +
					  x[j]->surface_charge->sigma1) *
		  (x[j]->surface_charge->specific_area *
		   x[j]->surface_charge->grams) / F_C_MOL);
    }
    master_ptr0 =
      surface_get_psi_master (x[j]->surface_charge->name, SURF_PSI);
    master_ptr1 =
      surface_get_psi_master (x[j]->surface_charge->name, SURF_PSI1);
    master_ptr2 =
      surface_get_psi_master (x[j]->surface_charge->name, SURF_PSI2);
    unknown_ptr0 = x[master_ptr0->unknown->number];
    unknown_ptr1 = x[master_ptr1->unknown->number];
    unknown_ptr2 = x[master_ptr2->unknown->number];

    charge0 = unknown_ptr0->f;
    charge1 = unknown_ptr1->f;
    if (dl_type_x != NO_DL)
    {
      charge2 =
	x[j]->surface_charge->sigma2 * (x[j]->surface_charge->specific_area *
					x[j]->surface_charge->grams) /
	F_C_MOL;
    }
    else
    {
      charge2 = unknown_ptr2->f;
    }
    sum = 0;
    for (k = 0; k < x[j]->count_comp_unknowns; k++)
    {
      sum +=
	x[j]->comp_unknowns[k]->moles *
	x[j]->comp_unknowns[k]->master[0]->s->z;
    }
    output_msg (OUTPUT_MESSAGE, "\t%11.3e  Surface charge, plane 0, eq\n",
		(double) charge0 + sum);
    output_msg (OUTPUT_MESSAGE, "\t%11.3e  Surface charge, plane 1, eq\n",
		(double) charge1);
    output_msg (OUTPUT_MESSAGE, "\t%11.3e  Surface charge, plane 2, eq\n",
		(double) charge2);
    output_msg (OUTPUT_MESSAGE,
		"\t%11.3e  Sum of surface charge, all planes, eq\n\n",
		(double) (charge0 + sum + charge1 + charge2));
    if (x[j]->type == SURFACE_CB)
    {
      if ((x[j]->surface_charge->specific_area *
	   x[j]->surface_charge->grams) > 0)
      {
	output_msg (OUTPUT_MESSAGE, "\t%11.3e  sigma, plane 0, C/m**2\n",
		    x[j]->surface_charge->sigma0);
	output_msg (OUTPUT_MESSAGE, "\t%11.3e  sigma, plane 1, C/m**2\n",
		    x[j]->surface_charge->sigma1);
	output_msg (OUTPUT_MESSAGE, "\t%11.3e  sigma, plane 2, C/m**2\n",
		    x[j]->surface_charge->sigma2);
	output_msg (OUTPUT_MESSAGE,
		    "\t%11.3e  sigma, diffuse layer, C/m**2\n\n",
		    x[j]->surface_charge->sigmaddl);
      }
      else
      {
	output_msg (OUTPUT_MESSAGE, "\tundefined  sigma, C/m**2\n");
      }
      output_msg (OUTPUT_MESSAGE, "\t%11.3e  psi, plane 0, V\n",
		  -master_ptr0->s->la * LOG_10 * R_KJ_DEG_MOL * tk_x /
		  F_KJ_V_EQ);
      output_msg (OUTPUT_MESSAGE, "\t%11.3e  psi, plane 1, V\n",
		  -master_ptr1->s->la * LOG_10 * R_KJ_DEG_MOL * tk_x /
		  F_KJ_V_EQ);
      output_msg (OUTPUT_MESSAGE, "\t%11.3e  psi, plane 2, V\n\n",
		  -master_ptr2->s->la * LOG_10 * R_KJ_DEG_MOL * tk_x /
		  F_KJ_V_EQ);
      output_msg (OUTPUT_MESSAGE, "\t%11.3e  exp(-F*psi/RT), plane 0\n",
		  exp (master_ptr0->s->la * LOG_10));
      output_msg (OUTPUT_MESSAGE, "\t%11.3e  exp(-F*psi/RT), plane 1\n",
		  exp (master_ptr1->s->la * LOG_10));
      output_msg (OUTPUT_MESSAGE, "\t%11.3e  exp(-F*psi/RT), plane 2\n\n",
		  exp (master_ptr2->s->la * LOG_10));

      output_msg (OUTPUT_MESSAGE, "\t%11.3e  capacitance 0-1, F/m^2\n",
		  x[j]->surface_charge->capacitance[0]);
      output_msg (OUTPUT_MESSAGE, "\t%11.3e  capacitance 1-2, F/m^2\n",
		  x[j]->surface_charge->capacitance[1]);

      if (x[j]->surface_comp->phase_name != NULL)
      {
	output_msg (OUTPUT_MESSAGE, "\t%11.3e  specific area, m^2/mol %s\n",
		    (double) x[j]->surface_charge->specific_area,
		    x[j]->surface_comp->phase_name);
	output_msg (OUTPUT_MESSAGE,
		    "\t%11.3e  m^2 for %11.3e moles of %s\n\n",
		    (double) (x[j]->surface_charge->grams *
			      x[j]->surface_charge->specific_area),
		    (double) x[j]->surface_charge->grams,
		    x[j]->surface_comp->phase_name);
      }
      else if (x[j]->surface_comp->rate_name != NULL)
      {
	output_msg (OUTPUT_MESSAGE, "\t%11.3e  specific area, m^2/mol %s\n",
		    (double) x[j]->surface_charge->specific_area,
		    x[j]->surface_comp->rate_name);
	output_msg (OUTPUT_MESSAGE,
		    "\t%11.3e  m^2 for %11.3e moles of %s\n\n",
		    (double) (x[j]->surface_charge->grams *
			      x[j]->surface_charge->specific_area),
		    (double) x[j]->surface_charge->grams,
		    x[j]->surface_comp->rate_name);
      }
      else
      {
	output_msg (OUTPUT_MESSAGE, "\t%11.3e  specific area, m^2/g\n",
		    (double) x[j]->surface_charge->specific_area);
	output_msg (OUTPUT_MESSAGE, "\t%11.3e  m^2 for %11.3e g\n\n",
		    (double) (x[j]->surface_charge->specific_area *
			      x[j]->surface_charge->grams),
		    (double) x[j]->surface_charge->grams);
      }
      if (dl_type_x != NO_DL)
	print_diffuse_layer (x[j]->surface_charge);
      output_msg (OUTPUT_MESSAGE, "\n");
/*
 *   Heading for species
 */
      for (k = j - 1; k < count_unknowns; k++)
      {
	if (x[k]->type != SURFACE)
	  continue;
	if (x[j] != x[k]->potential_unknown)
	  continue;
	master_ptr = x[k]->master[0];
	output_msg (OUTPUT_MESSAGE, "%-14s\n", x[k]->master[0]->elt->name);
	output_msg (OUTPUT_MESSAGE, "\t%11.3e  moles", (double) x[k]->moles);
	if (x[k]->surface_comp->phase_name != NULL)
	{
	  output_msg (OUTPUT_MESSAGE, "\t[%g mol/(mol %s)]\n",
		      (double) x[k]->surface_comp->phase_proportion,
		      x[k]->surface_comp->phase_name);
	}
	else if (x[k]->surface_comp->rate_name != NULL)
	{
	  output_msg (OUTPUT_MESSAGE,
		      "\t[%g mol/(mol kinetic reactant %s)]\n",
		      (double) x[k]->surface_comp->phase_proportion,
		      x[k]->surface_comp->rate_name);
	}
	else
	{
	  output_msg (OUTPUT_MESSAGE, "\n");
	}
	output_msg (OUTPUT_MESSAGE, "\t%-20s%12s%12s%12s%12s\n", " ", " ",
		    "Mole", " ", "Log");
	output_msg (OUTPUT_MESSAGE, "\t%-20s%12s%12s%12s%12s\n\n", "Species",
		    "Moles", "Fraction", "Molality", "Molality");
	for (i = 0; i < count_species_list; i++)
	{
	  if (species_list[i].master_s != master_ptr->s)
	    continue;
/*
 *   Print species data
 */
	  if (x[k]->moles >= MIN_RELATED_SURFACE)
	  {
	    molfrac =
	      (LDBLE) (species_list[i].s->moles) / x[k]->moles *
	      species_list[i].s->equiv;
	  }
	  else
	  {
	    molfrac = 0.0;
	  }
	  output_msg (OUTPUT_MESSAGE, "\t%-20s%12.3e%12.3f%12.3e%12.3f\n",
		      species_list[i].s->name,
		      (double) species_list[i].s->moles,
		      (double) molfrac,
		      (double) (species_list[i].s->moles / mass_water_aq_x),
		      log10 (species_list[i].s->moles / mass_water_aq_x));
	}
	output_msg (OUTPUT_MESSAGE, "\n");
      }
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
print_totals (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Print total concentrations of elements, molality and moles.
 */
  int i, pure_water;
  LDBLE EC;

  if (pr.totals == FALSE || pr.all == FALSE)
    return (OK);
  print_centered ("Solution composition");
  pure_water = TRUE;
  output_msg (OUTPUT_MESSAGE, "\t%-15s%12s%12s\n\n", "Elements", "Molality",
	      "Moles");
  for (i = 0; i < count_unknowns; i++)
  {
    if (x[i] == alkalinity_unknown)
    {
      output_msg (OUTPUT_MESSAGE, "\t%-15s%12.3e%12.3e\n",
		  x[i]->total->description,
		  (double) (x[i]->f / mass_water_aq_x), (double) x[i]->f);
      pure_water = FALSE;
    }
    if (x[i] == ph_unknown)
      continue;
    if (x[i] == pe_unknown)
      continue;
    if (x[i] == charge_balance_unknown)
    {
      output_msg (OUTPUT_MESSAGE, "\t%-15s%12.3e%12.3e", x[i]->description,
		  (double) (x[i]->sum / mass_water_aq_x), (double) x[i]->sum);
      output_msg (OUTPUT_MESSAGE, "  Charge balance\n");
      pure_water = FALSE;
      continue;
    }
    if (x[i]->type == SOLUTION_PHASE_BOUNDARY)
    {
      output_msg (OUTPUT_MESSAGE, "\t%-15s%12.3e%12.3e", x[i]->description,
		  (double) (x[i]->sum / mass_water_aq_x), (double) x[i]->sum);
      output_msg (OUTPUT_MESSAGE, "  Equilibrium with %s\n",
		  x[i]->phase->name);
      pure_water = FALSE;
      continue;
    }
    if (x[i]->type == MB)
    {
      output_msg (OUTPUT_MESSAGE, "\t%-15s%12.3e%12.3e\n", x[i]->description,
		  (double) (x[i]->sum / mass_water_aq_x), (double) x[i]->sum);
      pure_water = FALSE;
    }
  }

  if (pure_water == TRUE)
  {
    output_msg (OUTPUT_MESSAGE, "\t%-15s\n", "Pure water");
  }
/*
 *   Description of solution
 */
  output_msg (OUTPUT_MESSAGE, "\n");
  print_centered ("Description of solution");
/*
 *   pH
 */
  output_msg (OUTPUT_MESSAGE, "%45s%7.3f    ", "pH  = ",
	      (double) (-(s_hplus->la)));
  if (ph_unknown == NULL)
  {
    output_msg (OUTPUT_MESSAGE, "\n");
  }
  else if (ph_unknown == charge_balance_unknown)
  {
    output_msg (OUTPUT_MESSAGE, "  Charge balance\n");
  }
  else if (ph_unknown->type == SOLUTION_PHASE_BOUNDARY)
  {
    output_msg (OUTPUT_MESSAGE, "  Equilibrium with %s\n",
		ph_unknown->phase->name);
  }
  else if (ph_unknown->type == ALK)
  {
    output_msg (OUTPUT_MESSAGE, "  Adjust alkalinity\n");
  }
/*
 *    pe
 */
  output_msg (OUTPUT_MESSAGE, "%45s%7.3f    ", "pe  = ",
	      (double) (-(s_eminus->la)));
  if (pe_unknown == NULL)
  {
    output_msg (OUTPUT_MESSAGE, "\n");
  }
  else if (pe_unknown == charge_balance_unknown)
  {
    output_msg (OUTPUT_MESSAGE, "  Charge balance\n");
  }
  else if (pe_unknown->type == SOLUTION_PHASE_BOUNDARY)
  {
    output_msg (OUTPUT_MESSAGE, "  Equilibrium with %s\n",
		pe_unknown->phase->name);
  }
  else if (pe_unknown->type == MH)
  {
    output_msg (OUTPUT_MESSAGE, "  Adjusted to redox equilibrium\n");
  }
/*
 *   Others
 */
  EC = calc_SC ();
  if (EC > 0) {
  output_msg (OUTPUT_MESSAGE, "%45s%i\n", "Specific Conductance (uS/cm, 25 oC) = ",
	      (int) EC);
  }
  output_msg (OUTPUT_MESSAGE, "%45s%7.3f\n", "Activity of water  = ",
	      exp (s_h2o->la * LOG_10));
  output_msg (OUTPUT_MESSAGE, "%45s%11.3e\n", "Ionic strength  = ",
	      (double) mu_x);
  output_msg (OUTPUT_MESSAGE, "%45s%11.3e\n", "Mass of water (kg)  = ",
	      (double) mass_water_aq_x);
  if (alkalinity_unknown == NULL)
  {
    output_msg (OUTPUT_MESSAGE, "%45s%11.3e\n",
		"Total alkalinity (eq/kg)  = ",
		(double) (total_alkalinity / mass_water_aq_x));
  }
  if (carbon_unknown == NULL)
  {
    output_msg (OUTPUT_MESSAGE, "%45s%11.3e\n", "Total carbon (mol/kg)  = ",
		(double) (total_carbon / mass_water_aq_x));
  }
  output_msg (OUTPUT_MESSAGE, "%45s%11.3e\n", "Total CO2 (mol/kg)  = ",
	      (double) (total_co2 / mass_water_aq_x));
  output_msg (OUTPUT_MESSAGE, "%45s%7.3f\n", "Temperature (deg C)  = ",
	      (double) tc_x);
  output_msg (OUTPUT_MESSAGE, "%45s%11.3e\n", "Electrical balance (eq)  = ",
	      (double) cb_x);
  output_msg (OUTPUT_MESSAGE, "%45s%6.2f\n",
	      "Percent error, 100*(Cat-|An|)/(Cat+|An|)  = ",
	      (double) (100 * cb_x / total_ions_x));
  output_msg (OUTPUT_MESSAGE, "%45s%3d\n", "Iterations  = ", iterations);
  if (pitzer_model == TRUE)
  {
    output_msg (OUTPUT_MESSAGE, "%45s%3d\n", "Gamma iterations  = ",
		gamma_iterations);
    output_msg (OUTPUT_MESSAGE, "%45s%9.5f\n", "Osmotic coefficient  = ",
		COSMOT);
    output_msg (OUTPUT_MESSAGE, "%45s%9.5f\n", "Density of water  = ", DW0);
  }
  output_msg (OUTPUT_MESSAGE, "%45s%e\n", "Total H  = ", (double) total_h_x);
  output_msg (OUTPUT_MESSAGE, "%45s%e\n", "Total O  = ", (double) total_o_x);
  output_msg (OUTPUT_MESSAGE, "\n");

  return (OK);
}

/* ---------------------------------------------------------------------- */
int
print_user_print (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Print with user defined BASIC print routine
 */
  int i;
  struct kinetics *kinetics_ptr;

  char command[] = "run";

  if (pr.user_print == FALSE || pr.all == FALSE)
    return (OK);
  if (user_print->commands == NULL)
    return (OK);
  kinetics_ptr = NULL;
  if (use.kinetics_in == TRUE)
  {
    kinetics_ptr = use.kinetics_ptr;
    if (state == TRANSPORT || state == PHAST || state == ADVECTION)
    {
      use.kinetics_ptr = kinetics_bsearch (use.n_kinetics_user, &i);
    }
    else
    {
      use.kinetics_ptr = kinetics_bsearch (-2, &i);
    }
  }
  print_centered ("User print");
  if (user_print->new_def == TRUE)
  {
    /*            basic_renumber(user_print->commands, &user_print->linebase, &user_print->varbase, &user_print->loopbase); */
    if (basic_compile
	(user_print->commands, &user_print->linebase, &user_print->varbase,
	 &user_print->loopbase) != 0)
    {
      error_msg ("Fatal Basic error in USER_PRINT.", STOP);
    }
    user_print->new_def = FALSE;
  }
  if (basic_run
      (command, user_print->linebase, user_print->varbase,
       user_print->loopbase) != 0)
  {
    error_msg ("Fatal Basic error in USER_PRINT.", STOP);
  }
  output_msg (OUTPUT_MESSAGE, "\n");
  if (use.kinetics_in == TRUE)
  {
    use.kinetics_ptr = kinetics_ptr;
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
print_using (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Print entities used in calculation
 */
  struct mix *mix_ptr;
  struct solution *solution_ptr;
  struct exchange *exchange_ptr;
  struct surface *surface_ptr;
  struct pp_assemblage *pp_assemblage_ptr;
  struct s_s_assemblage *s_s_assemblage_ptr;
  struct gas_phase *gas_phase_ptr;
  struct temperature *temperature_ptr;
  struct irrev *irrev_ptr;
  struct kinetics *kinetics_ptr;
  int n;

  if (pr.use == FALSE || pr.all == FALSE)
    return (OK);
  if (state < REACTION || phast == TRUE)
    return (OK);
/*
 *   Mixture or Solution
 */
  if (use.mix_in == TRUE)
  {
    if (state == TRANSPORT)
    {
      mix_ptr = mix_bsearch (use.n_mix_user, &n);
    }
    else
    {
      mix_ptr = mix_bsearch (use.n_mix_user_orig, &n);
    }
    if (mix_ptr == NULL)
    {
      mix_ptr = use.mix_ptr;
    }
    if (mix_ptr != NULL)
    {
      if (state == TRANSPORT)
      {
	output_msg (OUTPUT_MESSAGE, "Using mix %d.\t%s\n", use.n_mix_user,
		    mix_ptr->description);
      }
      else
      {
	output_msg (OUTPUT_MESSAGE, "Using mix %d.\t%s\n",
		    use.n_mix_user_orig, mix_ptr->description);
      }

    }
  }
  else
  {
    solution_ptr = solution_bsearch (use.n_solution_user, &n, TRUE);
    output_msg (OUTPUT_MESSAGE, "Using solution %d.\t%s\n",
		use.n_solution_user, solution_ptr->description);
  }
/*
 *   Exchange and surface
 */
  if (use.exchange_in == TRUE)
  {
    exchange_ptr = exchange_bsearch (use.n_exchange_user, &n);
    output_msg (OUTPUT_MESSAGE, "Using exchange %d.\t%s\n",
		use.n_exchange_user, exchange_ptr->description);
  }
  if (use.surface_in == TRUE)
  {
    surface_ptr = surface_bsearch (use.n_surface_user, &n);
    output_msg (OUTPUT_MESSAGE, "Using surface %d.\t%s\n", use.n_surface_user,
		surface_ptr->description);
  }
  if (use.pp_assemblage_in == TRUE)
  {
    pp_assemblage_ptr = pp_assemblage_bsearch (use.n_pp_assemblage_user, &n);
    output_msg (OUTPUT_MESSAGE, "Using pure phase assemblage %d.\t%s\n",
		use.n_pp_assemblage_user, pp_assemblage_ptr->description);
  }
  if (use.s_s_assemblage_in == TRUE)
  {
    s_s_assemblage_ptr =
      s_s_assemblage_bsearch (use.n_s_s_assemblage_user, &n);
    output_msg (OUTPUT_MESSAGE, "Using solid solution assemblage %d.\t%s\n",
		use.n_s_s_assemblage_user, s_s_assemblage_ptr->description);
  }
  if (use.gas_phase_in == TRUE)
  {
    gas_phase_ptr = gas_phase_bsearch (use.n_gas_phase_user, &n);
    output_msg (OUTPUT_MESSAGE, "Using gas phase %d.\t%s\n",
		use.n_gas_phase_user, gas_phase_ptr->description);
  }
  if (use.temperature_in == TRUE)
  {
    temperature_ptr = temperature_bsearch (use.n_temperature_user, &n);
    output_msg (OUTPUT_MESSAGE, "Using temperature %d.\t%s\n",
		use.n_temperature_user, temperature_ptr->description);
  }
  if (use.irrev_in == TRUE)
  {
    if (state != TRANSPORT || transport_step > 0)
    {
      irrev_ptr = irrev_bsearch (use.n_irrev_user, &n);
      output_msg (OUTPUT_MESSAGE, "Using reaction %d.\t%s\n",
		  use.n_irrev_user, irrev_ptr->description);
    }
  }
  if (use.kinetics_in == TRUE)
  {
    if (state == TRANSPORT || state == PHAST || state == ADVECTION)
    {
      kinetics_ptr = kinetics_bsearch (use.n_kinetics_user, &n);
    }
    else
    {
      kinetics_ptr = kinetics_bsearch (-2, &n);
    }
    output_msg (OUTPUT_MESSAGE, "Using kinetics %d.\t%s\n",
		use.n_kinetics_user, kinetics_ptr->description);
  }
  output_msg (OUTPUT_MESSAGE, "\n");
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
punch_gas_phase (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Prints selected gas phase data
 */
  int i, j;
  LDBLE p, total_moles, volume;
  LDBLE moles;

  if (punch.count_gases <= 0)
    return (OK);
  p = 0.0;
  total_moles = 0.0;
  volume = 0.0;
  if (gas_unknown != NULL && use.gas_phase_ptr != NULL)
  {
    if (use.gas_phase_ptr->type == PRESSURE)
    {
      use.gas_phase_ptr->total_moles = gas_unknown->moles;
      use.gas_phase_ptr->volume =
	use.gas_phase_ptr->total_moles * R_LITER_ATM * tk_x /
	use.gas_phase_ptr->total_p;
    }
    p = use.gas_phase_ptr->total_p;
    total_moles = use.gas_phase_ptr->total_moles;
    volume = total_moles * R_LITER_ATM * tk_x / use.gas_phase_ptr->total_p;
  }
  if (punch.high_precision == FALSE)
  {
    fpunchf ("pressure", "%12.4e\t", p);
    fpunchf ("total mol", "%12.4e\t", total_moles);
    fpunchf ("volume", "%12.4e\t", volume);
  }
  else
  {
    fpunchf ("pressure", "%20.12e\t", p);
    fpunchf ("total mol", "%20.12e\t", total_moles);
    fpunchf ("volume", "%20.12e\t", volume);
  }
  for (i = 0; i < punch.count_gases; i++)
  {
    moles = 0.0;
    if (use.gas_phase_ptr != NULL && punch.gases[i].phase != NULL)
    {
      for (j = 0; j < use.gas_phase_ptr->count_comps; j++)
      {

	if (use.gas_phase_ptr->comps[j].phase != punch.gases[i].phase)
	  continue;
	moles = use.gas_phase_ptr->comps[j].phase->moles_x;
	if (moles <= MIN_TOTAL)
	  moles = 0.0;
	break;
      }
    }
    if (punch.high_precision == FALSE)
    {
      fpunchf (sformatf ("g_%s", punch.gases[i].name), "%12.4e\t", moles);
    }
    else
    {
      fpunchf (sformatf ("g_%s", punch.gases[i].name), "%20.12e\t", moles);
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
punch_s_s_assemblage (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Prints solid solution composition if present
 */
  int i, j, k;
  int found;
  LDBLE moles;

/*
 *   Print solid solutions
 */
  for (k = 0; k < punch.count_s_s; k++)
  {
    found = FALSE;
    if (use.s_s_assemblage_ptr != NULL)
    {
      for (j = 0; j < use.s_s_assemblage_ptr->count_s_s; j++)
      {
	for (i = 0; i < use.s_s_assemblage_ptr->s_s[j].count_comps; i++)
	{
	  if (strcmp_nocase
	      (punch.s_s[k].name,
	       use.s_s_assemblage_ptr->s_s[j].comps[i].name) == 0)
	  {
	    if (use.s_s_assemblage_ptr->s_s[j].s_s_in == TRUE)
	    {
	      moles = use.s_s_assemblage_ptr->s_s[j].comps[i].moles;
	    }
	    else
	    {
	      moles = 0;
	    }
	    if (punch.high_precision == FALSE)
	    {
	      fpunchf (sformatf ("s_%s", punch.s_s[k].name), "%12.4e\t",
		       (double) moles);
	    }
	    else
	    {
	      fpunchf (sformatf ("s_%s", punch.s_s[k].name), "%20.12e\t",
		       (double) moles);
	    }
	    found = TRUE;
	    break;
	  }
	}
	if (found == TRUE)
	  break;
      }
    }
    if (found == FALSE)
    {
      if (punch.high_precision == FALSE)
      {
	fpunchf (sformatf ("s_%s", punch.s_s[k].name), "%12.4e\t", 0.0);
      }
      else
      {
	fpunchf (sformatf ("s_%s", punch.s_s[k].name), "%20.12e\t", 0.0);
      }
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
punch_totals (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Print total concentrations of elements, molality and moles.
 */
  int j;
  LDBLE molality;

  for (j = 0; j < punch.count_totals; j++)
  {
    if (punch.totals[j].master == NULL)
    {
      molality = 0.0;
    }
    else if (punch.totals[j].master->primary == TRUE)
    {
      molality = punch.totals[j].master->total_primary / mass_water_aq_x;
    }
    else
    {
      molality = punch.totals[j].master->total / mass_water_aq_x;
    }
    if (punch.high_precision == FALSE)
    {
      fpunchf (sformatf ("%s(mol/kgw)", punch.totals[j].name), "%12.4e\t",
	       molality);
    }
    else
    {
      fpunchf (sformatf ("%s(mol/kgw)", punch.totals[j].name), "%20.12e\t",
	       molality);
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
punch_molalities (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Print concentrations of species (aq, ex, surf)
 */
  int j;
  LDBLE molality;

  for (j = 0; j < punch.count_molalities; j++)
  {
    molality = 0.0;
    if (punch.molalities[j].s != NULL && punch.molalities[j].s->in == TRUE)
    {
      molality = punch.molalities[j].s->moles / mass_water_aq_x;
    }
    if (punch.high_precision == FALSE)
    {
      fpunchf (sformatf ("m_%s(mol/kgw)", punch.molalities[j].name),
	       "%12.4e\t", molality);
    }
    else
    {
      fpunchf (sformatf ("m_%s(mol/kgw)", punch.molalities[j].name),
	       "%20.12e\t", molality);
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
punch_activities (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Print concentrations of species (aq, ex, surf)
 */
  int j;
  LDBLE la;

  for (j = 0; j < punch.count_activities; j++)
  {
    la = -999.999;
    if (punch.activities[j].s != NULL && punch.activities[j].s->in == TRUE)
    {
      /*la = punch.activities[j].s->lm + punch.activities[j].s->lg;*/
      la = log_activity(punch.activities[j].s->name);
    }
    if (punch.high_precision == FALSE)
    {
      fpunchf (sformatf ("la_%s", punch.activities[j].name), "%12.4e\t", la);
    }
    else
    {
      fpunchf (sformatf ("la_%s", punch.activities[j].name), "%20.12e\t", la);
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
punch_pp_assemblage (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Prints masses of selected pure_phases in pp_assemblage
 */
  int i, j;
  LDBLE moles, delta_moles;

  for (i = 0; i < punch.count_pure_phases; i++)
  {
    delta_moles = 0;
    moles = 0.0;
    if (punch.pure_phases[i].phase != NULL)
    {
      for (j = 0; j < count_unknowns; j++)
      {
	if (x == NULL || x[j]->type != PP)
	  continue;
/*
 *   Print pure phase assemblage data
 */
	if (punch.pure_phases[i].phase != x[j]->pure_phase->phase)
	  continue;
	if (state != TRANSPORT && state != PHAST)
	{
	  moles = x[j]->moles;
	  delta_moles =
	    x[j]->moles - x[j]->pure_phase->moles - x[j]->pure_phase->delta;
	}
	else
	{
	  moles = x[j]->moles;
	  delta_moles = x[j]->moles - x[j]->pure_phase->initial_moles;
	}
	break;
      }
    }
    if (punch.high_precision == FALSE)
    {
      fpunchf (punch.pure_phases[i].name, "%12.4e\t", moles);
      fpunchf (sformatf ("d_%s", punch.pure_phases[i].name), "%12.4e\t",
	       delta_moles);
    }
    else
    {
      fpunchf (punch.pure_phases[i].name, "%20.12e\t", moles);
      fpunchf (sformatf ("d_%s", punch.pure_phases[i].name), "%20.12e\t",
	       delta_moles);
    }
  }
  return (OK);
}

#define PHAST_NULL(x)          (phast ? NULL : x)
/* ---------------------------------------------------------------------- */
int
punch_identifiers (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   prints series of integers to identify simulation number,
 *   state of calculations, reaction or transport step number,
 *   and temp, ph, pe, and mass of water for each line
 *   of selected output.
 */
  const char *sformat;
  const char *dformat;
  const char *gformat;
  int i, l;
  char token[MAX_LENGTH];
  if (punch.in == FALSE)
    return (OK);
  if (punch.high_precision == FALSE)
  {
    l = 12;
    sformat = "%12s\t";
    dformat = "%12d\t";
    gformat = "%12g\t";
  }
  else
  {
    l = 20;
    sformat = "%20s\t";
    dformat = "%20d\t";
    gformat = "%20g\t";
  }

/* 
 *   simulation or simul_tr
 */

  if (punch.sim == TRUE)
  {
    if (state != TRANSPORT && state != PHAST)
    {
      fpunchf (PHAST_NULL ("sim"), dformat, simulation);
    }
    else
    {
      fpunchf (PHAST_NULL ("sim"), dformat, simul_tr);
    }
  }
  if (punch.state == TRUE)
  {
    switch (state)
    {
    case 0:
      strcpy (token, "init");
      break;
    case 1:
      strcpy (token, "i_soln");
      break;
    case 2:
      strcpy (token, "i_exch");
      break;
    case 3:
      strcpy (token, "i_surf");
      break;
    case 4:
      strcpy (token, "i_gas");
      break;
    case 5:
      strcpy (token, "react");
      break;
    case 6:
      strcpy (token, "inverse");
      break;
    case 7:
      strcpy (token, "advect");
      break;
    case 8:
      strcpy (token, "transp");
      break;
    }
    fpunchf (PHAST_NULL ("state"), sformat, token);

  }
/* 
 *   solution number or cell number and time
 */
  if (punch.soln == TRUE)
  {
    if (state == TRANSPORT || state == PHAST)
    {
      fpunchf (PHAST_NULL ("soln"), dformat, cell);
    }
    else if (state == ADVECTION)
    {
      fpunchf (PHAST_NULL ("soln"), dformat, use.n_solution_user);
    }
    else if (state < REACTION)
    {
      fpunchf (PHAST_NULL ("soln"), dformat, use.solution_ptr->n_user);
    }
    else
    {
      if (use.mix_in == TRUE)
      {
	if (state != TRANSPORT)
	{
	  fpunchf (PHAST_NULL ("soln"), dformat, use.n_mix_user_orig);
	}
	else
	{
	  fpunchf (PHAST_NULL ("soln"), dformat, use.n_mix_user);
	}

      }
      else
      {
	fpunchf (PHAST_NULL ("soln"), dformat, use.n_solution_user);
      }
    }
  }
  if (punch.dist == TRUE)
  {
    if (state == ADVECTION)
    {
      fpunchf (PHAST_NULL ("dist_x"), dformat, use.n_solution_user);
    }
    else if (state == TRANSPORT)
    {
      fpunchf (PHAST_NULL ("dist_x"), gformat,
	       cell_data[cell - 1].mid_cell_x);
    }
    else
    {
      fpunchf (PHAST_NULL ("dist_x"), dformat, -99);
    }
  }
  if (punch.time == TRUE)
  {
    if (state == REACTION && incremental_reactions == TRUE
	&& use.kinetics_ptr != NULL)
    {
      if (use.kinetics_ptr->count_steps > 0)
      {
	kin_time_x = 0.;
	for (i = 0; i < reaction_step; i++)
	{
	  if (i < use.kinetics_ptr->count_steps)
	  {
	    kin_time_x += use.kinetics_ptr->steps[i];
	  }
	  else
	  {
	    kin_time_x +=
	      use.kinetics_ptr->steps[use.kinetics_ptr->count_steps - 1];
	  }
	}
      }
      else if (use.kinetics_ptr->count_steps < 0)
      {
	if (reaction_step > -use.kinetics_ptr->count_steps)
	{
	  kin_time_x = use.kinetics_ptr->steps[0];
	}
	else
	{
	  kin_time_x =
	    reaction_step * use.kinetics_ptr->steps[0] /
	    ((LDBLE) (-use.kinetics_ptr->count_steps));
	}
      }
    }
    if (state == REACTION)
    {
      fpunchf (PHAST_NULL ("time"), gformat, kin_time_x);
    }
    else if (state == TRANSPORT || state == PHAST)
    {
      fpunchf (PHAST_NULL ("time"), gformat,
	       initial_total_time + rate_sim_time);
    }
    else if (state == ADVECTION)
    {
      if (advection_kin_time_defined == TRUE)
      {
	fpunchf (PHAST_NULL ("time"), gformat,
		 initial_total_time + rate_sim_time);
      }
      else
      {
	fpunchf (PHAST_NULL ("time"), dformat, advection_step);
      }
    }
    else
    {
      fpunchf (PHAST_NULL ("time"), dformat, -99);
    }
  }

/*
 *   reaction or transport step
 */
  if (punch.step == TRUE)
  {
    if (state == REACTION)
    {
      fpunchf (PHAST_NULL ("step"), dformat, reaction_step);
    }
    else if (state == ADVECTION)
    {
      fpunchf (PHAST_NULL ("step"), dformat, advection_step);
    }
    else if (state == TRANSPORT)
    {
      fpunchf (PHAST_NULL ("step"), dformat, transport_step);
    }
    else
    {
      fpunchf (PHAST_NULL ("step"), dformat, -99);
    }
  }
  if (punch.ph == TRUE)
  {
    if (punch.high_precision == FALSE)
    {
      fpunchf ("pH", "%12g\t", (double) (-s_hplus->la));
    }
    else
    {
      fpunchf ("pH", "%20.12e\t", (double) (-s_hplus->la));
    }
  }
  if (punch.pe == TRUE)
  {
    if (punch.high_precision == FALSE)
    {
      fpunchf ("pe", "%12g\t", (double) (-s_eminus->la));
    }
    else
    {
      fpunchf ("pe", "%20.12e\t", (double) (-s_eminus->la));
    }
  }
  if (punch.rxn == TRUE)
  {
    if (state >= REACTION && use.irrev_in == TRUE)
    {
      if (punch.high_precision == FALSE)
      {
	fpunchf ("reaction", "%12.4e\t", step_x);
      }
      else
      {
	fpunchf ("reaction", "%20.12e\t", step_x);
      }
    }
    else
    {
      if (punch.high_precision == FALSE)
      {
	fpunchf ("reaction", "%12d\t", -99);
      }
      else
      {
	fpunchf ("reaction", "%20d\t", -99);
      }
    }
  }
  if (punch.temp == TRUE)
  {
    if (punch.high_precision == FALSE)
    {
      fpunchf ("temp(C)", "%12.3f\t", tc_x);
    }
    else
    {
      fpunchf ("temp(C)", "%20.12e\t", tc_x);
    }
  }
  if (punch.alk == TRUE)
  {
    if (punch.high_precision == FALSE)
    {
      fpunchf ("Alk(eq/kgw)", "%12g\t",
	       (double) (total_alkalinity / mass_water_aq_x));
    }
    else
    {
      fpunchf ("Alk(eq/kgw)", "%20.12e\t",
	       (double) (total_alkalinity / mass_water_aq_x));
    }
  }
  if (punch.mu == TRUE)
  {
    if (punch.high_precision == FALSE)
    {
      fpunchf ("mu", "%12g\t", (double) mu_x);
    }
    else
    {
      fpunchf ("mu", "%20.12e\t", (double) mu_x);
    }
  }
  if (punch.water == TRUE)
  {
    if (punch.high_precision == FALSE)
    {
      fpunchf ("mass_H2O", "%12g\t", (double) mass_water_aq_x);
    }
    else
    {
      fpunchf ("mass_H2O", "%20.12e\t", (double) mass_water_aq_x);
    }
  }
  if (punch.charge_balance == TRUE)
  {
    if (punch.high_precision == FALSE)
    {
      fpunchf ("charge(eq)", "%12g\t", (double) cb_x);
    }
    else
    {
      fpunchf ("charge(eq)", "%20.12e\t", (double) cb_x);
    }
  }
  if (punch.percent_error == TRUE)
  {
    if (punch.high_precision == FALSE)
    {
      fpunchf ("pct_err", "%12g\t", (double) (100 * cb_x / total_ions_x));
    }
    else
    {
      fpunchf ("pct_err", "%20.12e\t", (double) (100 * cb_x / total_ions_x));
    }
  }

  return (OK);
}

#undef PHAST_NULL
/* ---------------------------------------------------------------------- */
int
punch_saturation_indices (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Prints saturation indices of selected phases
 */
  int i;
  LDBLE si, iap;
  struct rxn_token *rxn_ptr;

  for (i = 0; i < punch.count_si; i++)
  {
    if (punch.si[i].phase == NULL || punch.si[i].phase->in == FALSE)
    {
      si = -999.999;
    }
    else
    {
/*
 *   Print saturation index
 */
      iap = 0.0;
      for (rxn_ptr = punch.si[i].phase->rxn_x->token + 1; rxn_ptr->s != NULL;
	   rxn_ptr++)
      {
	iap += rxn_ptr->s->la * rxn_ptr->coef;
      }
      si = -punch.si[i].phase->lk + iap;
    }
    if (punch.high_precision == FALSE)
    {
      fpunchf (sformatf ("si_%s", punch.si[i].name), "%12.4f\t", si);
    }
    else
    {
      fpunchf (sformatf ("si_%s", punch.si[i].name), "%20.12e\t", si);
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
punch_kinetics (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   prints kinetic reaction,
 *   should be called only on final kinetic step
 */
  int i, j;
  struct kinetics *kinetics_ptr;
  LDBLE moles, delta_moles;

  kinetics_ptr = NULL;
  if (use.kinetics_in == TRUE)
  {
    if (state == TRANSPORT || state == PHAST || state == ADVECTION)
    {
      kinetics_ptr = kinetics_bsearch (use.n_kinetics_user, &i);
    }
    else
    {
      kinetics_ptr = kinetics_bsearch (-2, &i);
    }
  }
  for (i = 0; i < punch.count_kinetics; i++)
  {
    moles = 0.0;
    delta_moles = 0.0;
    if (kinetics_ptr != NULL)
    {
      for (j = 0; j < kinetics_ptr->count_comps; j++)
      {
	if (strcmp_nocase
	    (punch.kinetics[i].name, kinetics_ptr->comps[j].rate_name) == 0)
	{
	  if (state != TRANSPORT && state != PHAST)
	  {
	    moles = kinetics_ptr->comps[j].m;
	    delta_moles = -kinetics_ptr->comps[j].moles;
	  }
	  else
	  {
	    moles = kinetics_ptr->comps[j].m;
	    delta_moles =
	      kinetics_ptr->comps[j].m - kinetics_ptr->comps[j].initial_moles;
	  }
	  break;
	}
      }
    }
    if (punch.high_precision == FALSE)
    {
      fpunchf (sformatf ("k_%s", punch.kinetics[i].name), "%12.4e\t", moles);
      fpunchf (sformatf ("dk_%s", punch.kinetics[i].name), "%12.4e\t",
	       -delta_moles);
    }
    else
    {
      fpunchf (sformatf ("k_%s", punch.kinetics[i].name), "%20.12e\t", moles);
      fpunchf (sformatf ("dk_%s", punch.kinetics[i].name), "%20.12e\t",
	       -delta_moles);
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
punch_user_punch (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Punch with user defined BASIC print routine
 */

  char command[] = "run";
  extern int n_user_punch_index;

  // MDL: did we call this punch?
  if (Pdo_punch == FALSE)
    return(OK);

  n_user_punch_index = 0;
  if (punch.user_punch == FALSE)
    return (OK);
  if (user_punch->commands == NULL)
    return (OK);
  if (user_punch->new_def == TRUE)
  {
    if (basic_compile
	(user_punch->commands, &user_punch->linebase, &user_punch->varbase,
	 &user_punch->loopbase) != 0)
    {
      error_msg ("Fatal Basic error in USER_PUNCH.", STOP);
    }
    user_punch->new_def = FALSE;
  }
  if (basic_run
      (command, user_punch->linebase, user_punch->varbase,
       user_punch->loopbase) != 0)
  {
    error_msg ("Fatal Basic error in USER_PUNCH.", STOP);
  }
  // MDL: set this flag back to FALSE
  Pdo_punch = FALSE;
  return (OK);
}

#ifdef PHREEQ98
/* ---------------------------------------------------------------------- */
int
punch_user_graph (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Graph with user defined BASIC print routine
 */
  char command[] = "run";

  colnr = 0;
/*      //if (pr.user_graph == FALSE || pr.all == FALSE) return(OK); */
/*      //if (punch.user_punch == FALSE) return(OK); */
/*          //if (punch.in == FALSE) return(OK); */
  if (user_graph->commands == NULL)
    return (OK);
  if (((state == INITIAL_SOLUTION) || (state == INITIAL_EXCHANGE)
       || (state == INITIAL_SURFACE) || (state == INITIAL_GAS_PHASE))
      && (graph_initial_solutions == FALSE))
    return (OK);
  if (FirstCallToUSER_GRAPH)
    AddSeries = TRUE;
  if (state == REACTION)
  {
    /*if (reaction_step == 1) AddSeries = TRUE;
       else AddSeries = FALSE; */
    if (reaction_step > 1)
      AddSeries = FALSE;
  }
  if (state == ADVECTION)
  {
    if (advection_step == 0 && graph_initial_solutions == FALSE)
      return (OK);
    if (((chart_type == 1) && (advection_step == punch_ad_modulus)) ||
	((chart_type == 0) && (advection_step != prev_advection_step)))
      AddSeries = TRUE;
    else
      AddSeries = FALSE;
  }
  if (state == TRANSPORT)
  {
    if (transport_step == 0 && graph_initial_solutions == FALSE)
      return (OK);
    if (((chart_type == 1) && (transport_step == punch_modulus)) ||
	((chart_type == 0) && (transport_step != prev_transport_step)))
      AddSeries = TRUE;
    else
      AddSeries = FALSE;
  }
  if (user_graph->new_def == TRUE)
  {
    if (basic_compile
	(user_graph->commands, &user_graph->linebase, &user_graph->varbase,
	 &user_graph->loopbase) != 0)
    {
      error_msg ("Fatal Basic error in USER_GRAPH.", STOP);
    }
    user_graph->new_def = FALSE;
  }
  if (basic_run
      (command, user_graph->linebase, user_graph->varbase,
       user_graph->loopbase) != 0)
  {
    error_msg ("Fatal Basic error in USER_GRAPH.", STOP);
  }
  if (state == ADVECTION)
    prev_advection_step = advection_step;
  if (state == TRANSPORT)
    prev_transport_step = transport_step;
  /*if (state == REACTION) prev_reaction_step = reaction_step; */

  FirstCallToUSER_GRAPH = FALSE;

  return (OK);
}
#endif

#if defined(HDF5_CREATE)
extern void HDFWriteHyperSlabV (const char *name, const char *format,
				va_list argptr);
#endif

#if defined(USE_MPI) && defined(HDF5_CREATE) && defined(MERGE_FILES)
extern int Merge_fpunchf (const int length, const char *format,
			  va_list argptr);
#endif

int output_message (const int type, const char *err_str, const int stop,
		    const char *format, va_list args);

int
fpunchf (const char *name, const char *format, ...)
{
  int retval = 0;
  va_list args;

  va_start (args, format);
  retval = output_message (OUTPUT_PUNCH, name, CONTINUE, format, args);
  va_end (args);

  return retval;
}

int
fpunchf_user (int user_index, const char *format, ...)
{
  static int s_warning = 0;
  int retval = 0;
  va_list args;
  static char buffer[80];
  char *name;

  if (user_index < user_punch_count_headings)
  {
    name = user_punch_headings[user_index];
  }
  else
  {
    if (s_warning == 0)
    {
      sprintf (error_string,
	       "USER_PUNCH: Headings count doesn't match number of calls to PUNCH.\n");
      warning_msg (error_string);
      s_warning = 1;
    }
    sprintf (buffer, "no_heading_%d",
	     (user_index - user_punch_count_headings) + 1);
    name = buffer;
  }

  va_start (args, format);
  retval = output_message (OUTPUT_PUNCH, name, CONTINUE, format, args);
  va_end (args);

  return retval;
}

#ifdef SAVE
int
fpunchf (const char *name, const char *format, ...)
{
  int retval = 0;
  va_list args;
  static char big_buffer[200];

#if defined(HDF5_CREATE)
  if (pr.hdf == TRUE)
  {
    va_start (args, format);
    HDFWriteHyperSlabV (name, format, args);
    va_end (args);
  }
#endif
  if (pr.punch == TRUE && punch.in == TRUE)
  {
    va_start (args, format);
    retval = output_message (OUTPUT_PUNCH, NULL, CONTINUE, format, args);
    va_end (args);
#if defined(USE_MPI) && defined(HDF5_CREATE) && defined(MERGE_FILES)
    va_start (args, format);
    retval = vsprintf (big_buffer, format, args);
    Merge_fpunchf (retval, format, args);
    va_end (args);
#endif
#if defined(SWIG_SHARED_OBJ)
    va_start (args, format);
    AddSelectedOutput (name, format, args);
    va_end (args);
#endif
  }
  return retval;
}

int
fpunchf_user (int user_index, const char *format, ...)
{
  extern int warning_msg (const char *err_str);
  static int s_warning = 0;
  int retval = 0;
  va_list args;
  static char buffer[80];
  char *name;

  if (user_index < user_punch_count_headings)
  {
    name = user_punch_headings[user_index];
  }
  else
  {
    if (s_warning == 0)
    {
      sprintf (error_string,
	       "USER_PUNCH: Headings count doesn't match number of calls to PUNCH.\n");
      warning_msg (error_string);
      s_warning = 1;
    }
    sprintf (buffer, "no_heading_%d",
	     (user_index - user_punch_count_headings) + 1);
    name = buffer;
  }

#if defined(HDF5_CREATE)
  if (pr.hdf == TRUE)
  {
    va_start (args, format);
    HDFWriteHyperSlabV (name, format, args);
    va_end (args);
  }
#endif

  if (punch.in == TRUE && pr.punch == TRUE)
  {
    va_start (args, format);
    retval = output_message (OUTPUT_PUNCH, NULL, CONTINUE, format, args);
    va_end (args);
#if defined(USE_MPI) && defined(HDF5_CREATE) && defined(MERGE_FILES)
    va_start (args, format);
    Merge_fpunchf (retval, format, args);
    va_end (args);
#endif
#if defined(SWIG_SHARED_OBJ)
    va_start (args, format);
    AddSelectedOutput (name, format, args);
    va_end (args);
#endif
  }
  return retval;
}
#endif

char *
sformatf (const char *format, ...)
{
#if defined(HDF5_CREATE) || defined SWIG_SHARED_OBJ
  static char scratch[240];
  va_list args;

  va_start (args, format);
  if (vsprintf (scratch, format, args) > 239)
  {
    error_msg ("buffer overwrite in sformatf", STOP);
  }
  va_end (args);
  return scratch;
#else
  return NULL;
#endif
}

/* ---------------------------------------------------------------------- */
int
print_alkalinity (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Prints description of solution, uses array species_list for
 *   order of aqueous species.
 */
  int i, j;
  struct species_list *alk_list;
  int count_alk_list;
  double min;

  if (pr.alkalinity == FALSE || pr.all == FALSE)
    return (OK);
  print_centered ("Distribution of alkalinity");
  alk_list =
    (struct species_list *)
    PHRQ_malloc ((size_t) (count_s_x * sizeof (struct species_list)));
  if (alk_list == NULL)
    malloc_error ();
  j = 0;
  for (i = 0; i < count_s_x; i++)
  {
    if (s_x[i]->alk == 0.0)
      continue;
    alk_list[j].master_s = s_hplus;
    alk_list[j].s = s_x[i];
    alk_list[j].coef = s_x[i]->alk;
    j++;
  }
  count_alk_list = j;
  min = fabs (censor * total_alkalinity / mass_water_aq_x);
  if (count_alk_list > 0)
  {
    output_msg (OUTPUT_MESSAGE, "\t%26s%11.3e\n\n",
		"Total alkalinity (eq/kgw)  = ",
		(double) (total_alkalinity / mass_water_aq_x));
    output_msg (OUTPUT_MESSAGE, "\t%-15s%12s%12s%10s\n\n", "Species",
		"Alkalinity", "Molality", "Alk/Mol");
    qsort (&alk_list[0], (size_t) count_alk_list,
	   (size_t) sizeof (struct species_list), species_list_compare_alk);
    for (i = 0; i < count_alk_list; i++)
    {
      if (fabs (alk_list[i].s->alk * (alk_list[i].s->moles) / mass_water_aq_x)
	  < min)
	continue;
      output_msg (OUTPUT_MESSAGE, "\t%-15s%12.3e%12.3e%10.2f\n",
		  alk_list[i].s->name,
		  (double) (alk_list[i].s->alk * (alk_list[i].s->moles) /
			    mass_water_aq_x),
		  (double) ((alk_list[i].s->moles) / mass_water_aq_x),
		  (double) (alk_list[i].s->alk));
    }
  }

  output_msg (OUTPUT_MESSAGE, "\n");
  alk_list = (struct species_list *) free_check_null (alk_list);
  return (OK);
}
