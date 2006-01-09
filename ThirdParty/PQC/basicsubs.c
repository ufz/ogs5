#define  EXTERNAL extern
#include "global.h"
#include "phqalloc.h"
#include "output.h"
#include "phrqproto.h"

static char const svnid[] =
  "$Id: basicsubs.c 4 2009-04-21 17:29:29Z delucia $";

struct system_species
{
  char *name;
  char *type;
  LDBLE moles;
};
struct system_species *sys;
int count_sys, max_sys;
LDBLE sys_tot;

/* ---------------------------------------------------------------------- */
LDBLE
activity (char *species_name)
/* ---------------------------------------------------------------------- */
{
  struct species *s_ptr;
  LDBLE a;
  if (svnid == NULL)
    fprintf (stderr, " ");

  s_ptr = s_search (species_name);
  if (s_ptr == s_h2o)
  {
    a = pow (10., s_h2o->la);
  }
  else if (s_ptr == s_eminus)
  {
    a = pow (10., s_eminus->la);
  }
  else if (s_ptr == NULL || s_ptr->in == FALSE)
  {
    a = 1e-99;
  }
  else
  {
    a = pow (10., s_ptr->lm + s_ptr->lg);
  }
  return (a);
}

/* ---------------------------------------------------------------------- */
LDBLE
calc_SC (void)
/* ---------------------------------------------------------------------- */
{
  int i;
  LDBLE lm, a, z, Dw, SC, ff;

  SC = 0;
  for (i = 0; i < count_species_list; i++)
  {
    if (species_list[i].s->type == EX)
      continue;
    if (species_list[i].s->type == SURF)
      continue;
    if (i > 0
	&& strcmp (species_list[i].s->name, species_list[i - 1].s->name) == 0)
      continue;
    if (species_list[i].s == s_h2o)
      continue;
    if ((Dw = species_list[i].s->dw) == 0)
      continue;
    if ((z = fabs(species_list[i].s->z)) == 0)
      continue;

    lm = species_list[i].s->lm;
    if (lm > -9)
    {
/*
      if (z < 1.5) {
        ff = (mu_x < 0.36 ? 0.6 :
        sqrt(mu_x));
      }
      else {
        ff = (mu_x < pow(0.4*z, 2.0) ? 0.4 :
        sqrt(mu_x) / z);
      }
*/
      ff = (mu_x < .36 * z ? 0.6 / sqrt(z) :
	    sqrt(mu_x) / z);
	  
      a = under (lm + ff * species_list[i].s->lg);
      SC += a * z * z * Dw;
    }
  }
  SC *= 1e7 * F_C_MOL * F_C_MOL / (R_KJ_DEG_MOL * 298160.0);
  return (SC);
}

/* ---------------------------------------------------------------------- */
LDBLE
calc_logk_n (char *name)
/* ---------------------------------------------------------------------- */
{
  char token[MAX_LENGTH];
  int i;
  LDBLE lk;
  struct logk *logk_ptr;
  LDBLE logk[8];
  struct name_coef add_logk;

  for (i = 0; i < 8; i++)
  {
    logk[i] = 0.0;
  }
  strcpy (token, name);
  logk_ptr = logk_search (token);
  if (logk_ptr != NULL)
  {
    add_logk.name = token;
    add_logk.coef = 1.0;
    add_other_logk (logk, 1, &add_logk);
    lk = k_calc (logk, tk_x);
    return (lk);
  }
  return (-999.99);
}

/* ---------------------------------------------------------------------- */
LDBLE
calc_logk_p (char *name)
/* ---------------------------------------------------------------------- */
{
  int i, j;
  char token[MAX_LENGTH];
  struct phase *phase_ptr;
  LDBLE lk;
  LDBLE logk[8];

  strcpy (token, name);
  phase_ptr = phase_bsearch (token, &j, FALSE);
  if (phase_ptr != NULL)
  {
    for (i = 0; i < 8; i++)
    {
      logk[i] = 0.0;
    }
    select_log_k_expression (phase_ptr->logk, logk);
    add_other_logk (logk, phase_ptr->count_add_logk, phase_ptr->add_logk);
    lk = k_calc (logk, tk_x);
    return (lk);
  }
  return (-999.99);
}

/* ---------------------------------------------------------------------- */
LDBLE
calc_logk_s (char *name)
/* ---------------------------------------------------------------------- */
{
  int i;
  char token[MAX_LENGTH];
  struct species *s_ptr;
  LDBLE lk, logk[8];

  strcpy (token, name);
  s_ptr = s_search (token);
  if (s_ptr != NULL)
  {
    for (i = 0; i < 8; i++)
    {
      logk[i] = 0.0;
    }
    select_log_k_expression (s_ptr->logk, logk);
    add_other_logk (logk, s_ptr->count_add_logk, s_ptr->add_logk);
    lk = k_calc (logk, tk_x);
    return (lk);
  }
  return (-999.99);
}

/* ---------------------------------------------------------------------- */
LDBLE
calc_surface_charge (char *surface_name)
/* ---------------------------------------------------------------------- */
{
  char token[MAX_LENGTH], token1[MAX_LENGTH];
  char *ptr;
  int i, j, k;
  LDBLE charge;
  struct rxn_token_temp *token_ptr;
  struct master *master_ptr;

  /*
   *   Go through species, sum charge
   */
  charge = 0;
  for (k = 0; k < count_s_x; k++)
  {
    if (s_x[k]->type != SURF)
      continue;
    /*
     *   Match surface_name
     */
    count_trxn = 0;
    trxn_add (s_x[k]->rxn_s, 1.0, FALSE);	/* rxn_s is set in tidy_model */
    for (i = 1; i < count_trxn; i++)
    {
      token_ptr = &(trxn.token[i]);
      if (token_ptr->s->type != SURF)
	continue;
      master_ptr = trxn.token[i].s->primary;
      strcpy (token, master_ptr->elt->name);
      replace ("_", " ", token);
      ptr = token;
      copy_token (token1, &ptr, &j);
      if (strcmp (surface_name, token1) == 0)
      {
	charge += s_x[k]->moles * s_x[k]->z;
      }
    }
  }
  return (charge);
}

/* ---------------------------------------------------------------------- */
LDBLE
diff_layer_total (char *total_name, char *surface_name)
/* ---------------------------------------------------------------------- */
{
/*
 *   Provides total moles in DDL layer
 */
  int i, j, k, count_g;
  struct surface_charge *surface_charge_ptr1;
  char name[MAX_LENGTH], token[MAX_LENGTH];
  char surface_name_local[MAX_LENGTH];
  char *ptr;
  struct master *master_ptr;

  LDBLE mass_water_surface;
  LDBLE molality, moles_excess, moles_surface, charge;

  if (use.surface_ptr == NULL || (dl_type_x == NO_DL &&
				  strcmp_nocase ("psi", total_name) != 0 &&
				  strcmp_nocase ("psi1", total_name) != 0 &&
				  strcmp_nocase ("psi2", total_name) != 0 &&
				  strcmp_nocase ("charge", total_name) != 0 &&
				  strcmp_nocase ("charge1", total_name) != 0
				  && strcmp_nocase ("charge2",
						    total_name) != 0
				  && strcmp_nocase ("sigma", total_name) != 0
				  && strcmp_nocase ("sigma1", total_name) != 0
				  && strcmp_nocase ("sigma2",
						    total_name) != 0))
    return (0);

/*
 *   Find surface...
 */
  i = 0;
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
    else if (use.surface_ptr->type == CD_MUSIC)
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
    if (surface_name != NULL)
    {
      if (strcmp (name, surface_name) == 0)
	break;
    }
    else
    {
      break;
    }
  }
  if (j >= count_unknowns)
    return (0);
  strcpy (surface_name_local, name);
  /*
   *   Psi, charge, sigma
   */
  if (strcmp_nocase ("psi", total_name) == 0)
  {
    /*if (use.surface_ptr->edl == TRUE) { */
    if (use.surface_ptr->type == DDL)
    {
      return ((LDBLE) (x[j]->master[0]->s->la * 2 * R_KJ_DEG_MOL *
		       tk_x * LOG_10 / F_KJ_V_EQ));
    }
    else if (use.surface_ptr->type == CD_MUSIC)
    {
      master_ptr = surface_get_psi_master (surface_name, SURF_PSI);
      if (master_ptr != NULL)
      {
	/*return((LDBLE) (-x[j]->surface_charge->psi_master->s->la * R_KJ_DEG_MOL * tk_x * LOG_10 / F_KJ_V_EQ)); */
	return ((LDBLE)
		(-master_ptr->s->la * R_KJ_DEG_MOL * tk_x * LOG_10 /
		 F_KJ_V_EQ));
      }
      else
      {
	return (0.0);
      }
    }
    else
    {
      return (0);
    }
  }
  else if (strcmp_nocase ("psi1", total_name) == 0)
  {
    master_ptr = surface_get_psi_master (surface_name, SURF_PSI1);
    if (master_ptr != NULL)
    {
      return ((LDBLE)
	      (-master_ptr->s->la * R_KJ_DEG_MOL * tk_x * LOG_10 /
	       F_KJ_V_EQ));
    }
    else
    {
      return (0.0);
    }
  }
  else if (strcmp_nocase ("psi2", total_name) == 0)
  {
    master_ptr = surface_get_psi_master (surface_name, SURF_PSI2);
    if (master_ptr != NULL)
    {
      return ((LDBLE)
	      (-master_ptr->s->la * R_KJ_DEG_MOL * tk_x * LOG_10 /
	       F_KJ_V_EQ));
    }
    else
    {
      return (0.0);
    }
  }
  else if (strcmp_nocase ("charge", total_name) == 0)
  {
    if (use.surface_ptr->type == DDL && dl_type_x == NO_DL)
    {
      return ((LDBLE) (x[j]->f));
    }
    else if (use.surface_ptr->type == CD_MUSIC)
    {
      return ((LDBLE)
	      (x[j]->surface_charge->sigma0 *
	       (x[j]->surface_charge->specific_area *
		x[j]->surface_charge->grams) / F_C_MOL));
    }
    else
    {
      return (calc_surface_charge (surface_name_local));
    }
  }
  else if (strcmp_nocase ("charge1", total_name) == 0)
  {
    if (use.surface_ptr->type == CD_MUSIC)
    {
      return ((LDBLE)
	      (x[j]->surface_charge->sigma1 *
	       (x[j]->surface_charge->specific_area *
		x[j]->surface_charge->grams) / F_C_MOL));
    }
    else
    {
      return (0);
    }
  }
  else if (strcmp_nocase ("charge2", total_name) == 0)
  {
    if (use.surface_ptr->type == CD_MUSIC)
    {
      return ((LDBLE)
	      (x[j]->surface_charge->sigma2 *
	       (x[j]->surface_charge->specific_area *
		x[j]->surface_charge->grams) / F_C_MOL));
    }
    else
    {
      return (0);
    }
  }
  else if (strcmp_nocase ("sigma", total_name) == 0)
  {
    if (use.surface_ptr->type == DDL)
    {
      if (dl_type_x != NO_DL)
      {
	charge = calc_surface_charge (surface_name_local);
      }
      else
      {
	charge = x[j]->f;
      }
      if ((x[j]->surface_charge->specific_area *
	   x[j]->surface_charge->grams) > 0)
      {
	return ((LDBLE)
		(charge * F_C_MOL /
		 (x[j]->surface_charge->specific_area *
		  x[j]->surface_charge->grams)));
      }
      else
      {
	return (0);
      }
    }
    else if (use.surface_ptr->type == CD_MUSIC)
    {
      return ((LDBLE) (x[j]->surface_charge->sigma0));
    }
    else
    {
      return (0);
    }
  }
  else if (strcmp_nocase ("sigma1", total_name) == 0)
  {
    if (use.surface_ptr->type == CD_MUSIC)
    {
      return ((LDBLE) (x[j]->surface_charge->sigma1));
    }
    else
    {
      return (0);
    }
  }
  else if (strcmp_nocase ("sigma2", total_name) == 0)
  {
    if (use.surface_ptr->type == CD_MUSIC)
    {
      return ((LDBLE) (x[j]->surface_charge->sigma2));
    }
    else
    {
      return (0);
    }
  }
  else if (strcmp_nocase ("water", total_name) == 0)
  {
    if (dl_type_x != NO_DL)
    {
      return (x[j]->surface_charge->mass_water);
    }
    else
    {
      return (0);
    }
  }
/*
 *   find total moles of each element in diffuse layer...
 */
  surface_charge_ptr1 = x[j]->surface_charge;
  mass_water_surface = surface_charge_ptr1->mass_water;
  count_elts = 0;
  paren_count = 0;
  for (j = 0; j < count_s_x; j++)
  {
    if (s_x[j]->type > HPLUS)
      continue;
    molality = under (s_x[j]->lm);
    count_g = s_x[j]->diff_layer[i].count_g;

    moles_excess =
      mass_water_aq_x * molality * surface_charge_ptr1->g[count_g].g;
    moles_surface = mass_water_surface * molality + moles_excess;
/*
 *   Accumulate elements in diffuse layer
 */
    add_elt_list (s_x[j]->next_elt, moles_surface);
  }
  if (count_elts > 0)
  {
    qsort (elt_list, (size_t) count_elts,
	   (size_t) sizeof (struct elt_list), elt_list_compare);
    elt_list_combine ();
  }
/*
 *   Return totals
 */
  for (j = 0; j < count_elts; j++)
  {
    if (strcmp (elt_list[j].elt->name, total_name) == 0)
    {
      return ((LDBLE) elt_list[j].coef);
    }
  }
  return (0);
}

/* ---------------------------------------------------------------------- */
LDBLE
equi_phase (char *phase_name)
/* ---------------------------------------------------------------------- */
{
  int i, j;

  if (use.pp_assemblage_in == FALSE || use.pp_assemblage_ptr == NULL)
    return (0);
  for (j = 0; j < count_unknowns; j++)
  {
    if (x[j]->type != PP)
      continue;
    if (strcmp_nocase (x[j]->pure_phase->name, phase_name) == 0)
    {
      break;
    }
  }
/*
 *   Print pure phase assemblage data
 */
  if (j == count_unknowns)
  {
    /* if not an unknown */
    for (i = 0; i < use.pp_assemblage_ptr->count_comps; i++)
    {
      if (strcmp_nocase
	  (use.pp_assemblage_ptr->pure_phases[i].name, phase_name) == 0)
      {
	return (use.pp_assemblage_ptr->pure_phases[i].moles);
      }
    }
  }
  else
  {
    /* if an unknown */
    if (x[j]->moles < 0.0)
      x[j]->moles = 0.0;
    return (x[j]->moles);
  }
  return (0);
}

/* ---------------------------------------------------------------------- */
LDBLE
find_gas_comp (char *gas_comp_name)
/* ---------------------------------------------------------------------- */
{
  int i;

  if (use.gas_phase_in == FALSE || use.gas_phase_ptr == NULL)
    return (0);
  for (i = 0; i < use.gas_phase_ptr->count_comps; i++)
  {
    if (strcmp_nocase (use.gas_phase_ptr->comps[i].name, gas_comp_name) == 0)
    {
      return (use.gas_phase_ptr->comps[i].phase->moles_x);
    }
  }
  return (0);
}

/* ---------------------------------------------------------------------- */
LDBLE
find_misc1 (char *s_s_name)
/* ---------------------------------------------------------------------- */
{
  int j;

  if (use.s_s_assemblage_in == FALSE || use.s_s_assemblage_ptr == NULL)
    return (0.0);
  for (j = 0; j < use.s_s_assemblage_ptr->count_s_s; j++)
  {
    if (strcmp_nocase (use.s_s_assemblage_ptr->s_s[j].name, s_s_name) == 0)
    {
      if (use.s_s_assemblage_ptr->s_s[j].miscibility == TRUE)
      {
	return (use.s_s_assemblage_ptr->s_s[j].xb1);
      }
      else
      {
	return (1.0);
      }
    }
  }
  return (0);
}

/* ---------------------------------------------------------------------- */
LDBLE
find_misc2 (char *s_s_name)
/* ---------------------------------------------------------------------- */
{
  int j;

  if (use.s_s_assemblage_in == FALSE || use.s_s_assemblage_ptr == NULL)
    return (0.0);
  for (j = 0; j < use.s_s_assemblage_ptr->count_s_s; j++)
  {
    if (strcmp_nocase (use.s_s_assemblage_ptr->s_s[j].name, s_s_name) == 0)
    {
      if (use.s_s_assemblage_ptr->s_s[j].miscibility == TRUE)
      {
	return (use.s_s_assemblage_ptr->s_s[j].xb2);
      }
      else
      {
	return (1.0);
      }
    }
  }
  return (0);
}

/* ---------------------------------------------------------------------- */
LDBLE
find_s_s_comp (char *s_s_comp_name)
/* ---------------------------------------------------------------------- */
{
  int i, j;

  if (use.s_s_assemblage_in == FALSE || use.s_s_assemblage_ptr == NULL)
    return (0);
  for (j = 0; j < use.s_s_assemblage_ptr->count_s_s; j++)
  {
    for (i = 0; i < use.s_s_assemblage_ptr->s_s[j].count_comps; i++)
    {
      if (strcmp_nocase
	  (use.s_s_assemblage_ptr->s_s[j].comps[i].name, s_s_comp_name) == 0)
      {
	if (use.s_s_assemblage_ptr->s_s[j].s_s_in == TRUE)
	{
	  return (use.s_s_assemblage_ptr->s_s[j].comps[i].moles);
	}
	else
	{
	  return (0.0);
	}
      }
    }
  }
  return (0);
}

/* ---------------------------------------------------------------------- */
LDBLE
get_calculate_value (char *name)
/* ---------------------------------------------------------------------- */
/*
 *   Gets value from a calclate_value structure
 *	  arguments: name
 *	  return: LDBLE of value
 */
{
  struct calculate_value *calculate_value_ptr;
  calculate_value_ptr = calculate_value_search (name);
  if (calculate_value_ptr == NULL)
  {
    sprintf (error_string, "CALC_VALUE Basic function, %s not found.", name);
    error_msg (error_string, CONTINUE);
    input_error++;
    return (MISSING);
  }
  if (name == NULL)
  {
    sprintf (error_string, "Definition for calculated value not found, %s",
	     name);
    input_error++;
    error_msg (error_string, CONTINUE);
    return (MISSING);
  }
  if (calculate_value_ptr->calculated == FALSE)
  {
    sprintf (error_string,
	     "Calculated value used before it has been calculated, %s.\nCalculated values are evalutated in the order in which they were defined.",
	     name);
    input_error++;
    error_msg (error_string, CONTINUE);
    return (MISSING);
  }
  return (calculate_value_ptr->value);
}

/* ---------------------------------------------------------------------- */
LDBLE
kinetics_moles (char *kinetics_name)
/* ---------------------------------------------------------------------- */
{
  int i;

  if (use.kinetics_in == FALSE || use.kinetics_ptr == NULL)
    return (0);
  for (i = 0; i < use.kinetics_ptr->count_comps; i++)
  {
    if (strcmp_nocase (use.kinetics_ptr->comps[i].rate_name, kinetics_name) ==
	0)
    {
      return (use.kinetics_ptr->comps[i].m);
    }
  }

  sprintf (error_string, "No data for rate %s in KINETICS keyword.",
	   kinetics_name);
  warning_msg (error_string);
  return (0);
}

/* ---------------------------------------------------------------------- */
LDBLE
log_activity (char *species_name)
/* ---------------------------------------------------------------------- */
{
  struct species *s_ptr;
  LDBLE la;

  s_ptr = s_search (species_name);

  if (s_ptr == s_eminus)
  {
    la = s_eminus->la;
  }
  else if (s_ptr == NULL || s_ptr->in == FALSE)
  {
    la = -99.99;
  }
  else if (s_ptr == s_h2o)
  {
    la = s_h2o->la;
  }
  else
  {
    la = s_ptr->lm + s_ptr->lg;
  }
  return (la);
}

/* ---------------------------------------------------------------------- */
LDBLE
log_molality (char *species_name)
/* ---------------------------------------------------------------------- */
{
  struct species *s_ptr;
  LDBLE lm;

  s_ptr = s_search (species_name);

  if (s_ptr == s_eminus)
  {
    lm = -99.99;
  }
  else if (s_ptr == NULL || s_ptr->in == FALSE)
  {
    lm = -99.99;
  }
  else if (s_ptr == s_h2o)
  {
    lm = log10 (s_ptr->moles / mass_water_aq_x);
  }
  else
  {
    lm = s_ptr->lm;
  }
  return (lm);
}

/* ---------------------------------------------------------------------- */
LDBLE
molality (char *species_name)
/* ---------------------------------------------------------------------- */
{
  struct species *s_ptr;
  LDBLE m;

  s_ptr = s_search (species_name);
  if (s_ptr == NULL || s_ptr == s_eminus || s_ptr->in == FALSE)
  {
    m = 1e-99;
  }
  else
  {
    /* m = pow(10., s_ptr->lm); */
    m = s_ptr->moles / mass_water_aq_x;
  }
  return (m);
}

/* ---------------------------------------------------------------------- */
LDBLE
saturation_ratio (char *phase_name)
/* ---------------------------------------------------------------------- */
{
  struct rxn_token *rxn_ptr;
  struct phase *phase_ptr;
  int l;
  LDBLE si, iap;

  iap = 0.0;
  phase_ptr = phase_bsearch (phase_name, &l, FALSE);
  if (phase_ptr == NULL)
  {
    sprintf (error_string, "Mineral %s, not found.", phase_name);
    warning_msg (error_string);
    return (1e-99);
  }
  else if (phase_ptr->in != FALSE)
  {
    for (rxn_ptr = phase_ptr->rxn_x->token + 1; rxn_ptr->s != NULL; rxn_ptr++)
    {
      iap += rxn_ptr->s->la * rxn_ptr->coef;
    }
    si = iap - phase_ptr->lk;
    return (pow (10.0, si));
  }
  return (0.0);

}

/* ---------------------------------------------------------------------- */
int
saturation_index (char *phase_name, LDBLE * iap, LDBLE * si)
/* ---------------------------------------------------------------------- */
{
  struct rxn_token *rxn_ptr;
  struct phase *phase_ptr;
  int l;

  *si = -99.99;
  *iap = 0.0;
  phase_ptr = phase_bsearch (phase_name, &l, FALSE);
  if (phase_ptr == NULL)
  {
    sprintf (error_string, "Mineral %s, not found.", phase_name);
    warning_msg (error_string);
    *si = -99;
  }
  else if (phase_ptr->in != FALSE)
  {
    for (rxn_ptr = phase_ptr->rxn_x->token + 1; rxn_ptr->s != NULL; rxn_ptr++)
    {
      *iap += rxn_ptr->s->la * rxn_ptr->coef;
    }
    *si = *iap - phase_ptr->lk;
  }
  else
  {
    return (ERROR);
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
LDBLE
sum_match_gases (char *mytemplate, char *name)
/* ---------------------------------------------------------------------- */
{
  int i;
  LDBLE tot;
  struct elt_list *next_elt;

  if (use.gas_phase_in == FALSE || use.gas_phase_ptr == NULL)
    return (0);
  tot = 0;
  for (i = 0; i < use.gas_phase_ptr->count_comps; i++)
  {
    if (match_elts_in_species
	(use.gas_phase_ptr->comps[i].phase->formula, mytemplate) == TRUE)
    {
      if (name == NULL)
      {
	tot += use.gas_phase_ptr->comps[i].phase->moles_x;
      }
      else
      {
	for (next_elt = use.gas_phase_ptr->comps[i].phase->next_elt;
	     next_elt->elt != NULL; next_elt++)
	{
	  if (strcmp (next_elt->elt->name, name) == 0)
	  {
	    tot +=
	      next_elt->coef * use.gas_phase_ptr->comps[i].phase->moles_x;
	    break;
	  }
	}
      }
    }
  }
  return (tot);
}

/* ---------------------------------------------------------------------- */
LDBLE
sum_match_species (char *mytemplate, char *name)
/* ---------------------------------------------------------------------- */
{
  int i;
  LDBLE tot;
  struct elt_list *next_elt;

  count_elts = 0;
  paren_count = 0;
  tot = 0;
  for (i = 0; i < count_s_x; i++)
  {
    if (match_elts_in_species (s_x[i]->name, mytemplate) == TRUE)
    {
      if (name == NULL)
      {
	tot += s_x[i]->moles;
      }
      else
      {
	for (next_elt = s_x[i]->next_elt; next_elt->elt != NULL; next_elt++)
	{
	  if (strcmp (next_elt->elt->name, name) == 0)
	  {
	    tot += next_elt->coef * s_x[i]->moles;
	    break;
	  }
	}
      }
    }
  }
  return (tot);
}

/* ---------------------------------------------------------------------- */
LDBLE
sum_match_s_s (char *mytemplate, char *name)
/* ---------------------------------------------------------------------- */
{
  int i, j;
  LDBLE tot;
  struct elt_list *next_elt;

  if (use.s_s_assemblage_in == FALSE || use.s_s_assemblage_ptr == NULL)
    return (0);
  tot = 0;
  for (j = 0; j < use.s_s_assemblage_ptr->count_s_s; j++)
  {
    if (strcmp_nocase (use.s_s_assemblage_ptr->s_s[j].name, mytemplate) == 0)
    {
      if (use.s_s_assemblage_ptr->s_s[j].s_s_in == FALSE)
      {
	tot = 0;
	break;
      }
      for (i = 0; i < use.s_s_assemblage_ptr->s_s[j].count_comps; i++)
      {
	if (name == NULL)
	{
	  tot += use.s_s_assemblage_ptr->s_s[j].comps[i].moles;
	}
	else
	{
	  for (next_elt =
	       use.s_s_assemblage_ptr->s_s[j].comps[i].phase->next_elt;
	       next_elt->elt != NULL; next_elt++)
	  {
	    if (strcmp (next_elt->elt->name, name) == 0)
	    {
	      tot +=
		next_elt->coef *
		use.s_s_assemblage_ptr->s_s[j].comps[i].moles;
	      break;
	    }
	  }
	}
      }
      break;
    }
  }
  return (tot);
}

#ifdef SKIP
/* ---------------------------------------------------------------------- */
LDBLE
sum_match_s_s (char *template, char *name)
/* ---------------------------------------------------------------------- */
{
  int i, j;
  LDBLE tot;
  struct elt_list *next_elt;

  if (use.s_s_assemblage_in == FALSE || use.s_s_assemblage_ptr == NULL)
    return (0);
  tot = 0;
  for (j = 0; j < use.s_s_assemblage_ptr->count_s_s; j++)
  {
    for (i = 0; i < use.s_s_assemblage_ptr->s_s[j].count_comps; i++)
    {
      if (match_elts_in_species
	  (use.s_s_assemblage_ptr->s_s[j].comps[i].phase->formula,
	   template) == TRUE)
      {
	if (name == NULL)
	{
	  tot += use.s_s_assemblage_ptr->s_s[j].comps[i].moles;
	}
	else
	{
	  for (next_elt =
	       use.s_s_assemblage_ptr->s_s[j].comps[i].phase->next_elt;
	       next_elt->elt != NULL; next_elt++)
	  {
	    if (strcmp (next_elt->elt->name, name) == 0)
	    {
	      tot +=
		next_elt->coef *
		use.s_s_assemblage_ptr->s_s[j].comps[i].moles;
	      break;
	    }
	  }
	}
      }
    }
  }
  return (tot);
}
#endif
/* ---------------------------------------------------------------------- */
int
match_elts_in_species (char *name, char *mytemplate)
/* ---------------------------------------------------------------------- */
{
/*
 *	Makes a list of elements with their coefficients, stores elements
 *	in elt_list at position count_elts.  Global variable count_elts is
 *	updated with each stored element.  Also uses static global variable
 *	paren_count.
 *
 *	Arguments:
 *	   **t_ptr	input, point in token string to start looking
 *				  output, is next position to start looking
 *		 coef	 input, coefficient to multiply subscripts by
 */
  int i, i1, l, case_no, match;
  char c, c1;
  char *ptr, *ptr1, *replace_name, *name_ptr;
  LDBLE d;
  char element[MAX_LENGTH];
  char token[MAX_LENGTH], equal_list[MAX_LENGTH], elt_name[MAX_LENGTH];
  char token1[MAX_LENGTH], template1[MAX_LENGTH], equal_list1[MAX_LENGTH];
  char str[2];

  strcpy (token, name);
  squeeze_white (token);
  while (replace ("(+", "(", token));
  while (replace ("++++++", "+6", token));
  while (replace ("+++++", "+5", token));
  while (replace ("++++", "+4", token));
  while (replace ("+++", "+3", token));
  while (replace ("++", "+2", token));
  while (replace ("------", "-6", token));
  while (replace ("-----", "-5", token));
  while (replace ("----", "-4", token));
  while (replace ("---", "-3", token));
  while (replace ("--", "-2", token));

  ptr = token;
  count_match_tokens = 0;
  while ((c = *ptr) != '\0')
  {
    c1 = *(ptr + 1);
    str[0] = c;
    str[1] = '\0';
/*
 * New element
 */
    if (isupper ((int) c) || (c == 'e' && c1 == '-') || (c == '['))
    {
      /*
       *   Get new element and subscript
       */
      if (get_elt (&ptr, element, &l) == ERROR)
      {
	return (ERROR);
      }
      match_tokens[count_match_tokens].name = string_hsave (element);
      if (get_num (&ptr, &d) == ERROR)
      {
	return (ERROR);
      }
      match_tokens[count_match_tokens++].coef = d;
    }
    else
    {
      match_tokens[count_match_tokens].name = string_hsave (str);
      match_tokens[count_match_tokens++].coef = 1.0;
      ptr += 1;
    }
  }
  /*
   *  Replace elements with first of equivalent elements
   */
  strcpy (template1, mytemplate);
  squeeze_white (template1);
  ptr = template1;
  while (extract_bracket (&ptr, equal_list) == TRUE)
  {
    replace ("{", "", equal_list);
    replace ("}", "", equal_list);
    while (replace (",", " ", equal_list) == TRUE);
    ptr1 = equal_list;
    /*
     *   Get first name in a list from template
     */
    if (copy_token (elt_name, &ptr1, &l) == EMPTY)
    {
      sprintf (error_string,
	       "Expecting a nonempty list of element names in isotope sum. %s",
	       mytemplate);
      error_msg (error_string, CONTINUE);
      return (ERROR);
    }
    replace_name = string_hsave (elt_name);
    /*
     *   Replace in species all equivalent names from template
     */
    while (copy_token (elt_name, &ptr1, &l) != EMPTY)
    {
      name_ptr = string_hsave (elt_name);
      for (i = 0; i < count_match_tokens; i++)
      {
	if (name_ptr == match_tokens[i].name)
	{
	  match_tokens[i].name = replace_name;
	}
      }
    }
  }
  /*
   *  Combine contiguous elements
   */
  i1 = 0;
  for (i = 1; i < count_match_tokens; i++)
  {
    if ((isupper ((int) (match_tokens[i].name[0])) != FALSE)
	&& (match_tokens[i].name == match_tokens[i1].name))
    {
      match_tokens[i1].coef += match_tokens[i].coef;
    }
    else
    {
      i1++;
      match_tokens[i1].name = match_tokens[i].name;
      match_tokens[i1].coef = match_tokens[i].coef;
    }
  }
  count_match_tokens = i1 + 1;
  /*
   * write out string
   */
  token[0] = '\0';
  for (i = 0; i < count_match_tokens; i++)
  {
    strcat (token, match_tokens[i].name);
    if (match_tokens[i].coef != 1.0)
    {
      sprintf (token1, "%g", (double) match_tokens[i].coef);
      strcat (token, token1);
    }
  }
  /*
   *  Write a template name using first of equivalent elements
   */
  strcpy (template1, mytemplate);
  squeeze_white (template1);
  ptr = template1;
  while (extract_bracket (&ptr, equal_list) == TRUE)
  {
    strcpy (equal_list1, equal_list);
    replace ("{", "", equal_list);
    replace ("}", "", equal_list);
    while (replace (",", " ", equal_list) == TRUE);
    ptr1 = equal_list;
    /*
     *   Get first name in a list
     */
    if (copy_token (elt_name, &ptr1, &l) == EMPTY)
    {
      sprintf (error_string,
	       "Expecting a nonempty list of element names in isotope sum. %s",
	       mytemplate);
      error_msg (error_string, CONTINUE);
      return (ERROR);
    }
    replace_name = string_hsave (elt_name);
    replace (equal_list1, replace_name, template1);
    squeeze_white (template1);
    ptr = template1;
  }
  /*
   *   Compare string
   */
  /* Cases: 0 exact match
   *              1 leading wild card
   *              2 trailing wild card
   *              3 leading and trailing wild card
   */
  case_no = 0;
  if (template1[0] == '*')
    case_no = 1;
  l = (int) strlen (template1);
  if (template1[l - 1] == '*')
  {
    if (case_no != 1)
    {
      case_no = 2;
    }
    else
    {
      case_no = 3;
    }
  }
  while (replace ("*", "", template1));
  match = FALSE;
  switch (case_no)
  {
  case 0:
    /* exact match */
    if (strcmp (token, template1) == 0)
      match = TRUE;
    break;
  case 1:
    /* leading wild card */
    if ((ptr = strstr (token, template1)) == NULL)
    {
      match = FALSE;
    }
    else
    {
      if (strcmp (ptr, template1) == 0)
	match = TRUE;
    }
    break;
  case 2:
    /* trailing wild card */
    if (strstr (token, template1) == token)
      match = TRUE;
    break;
  case 3:
    /* trailing wild card */
    if (strstr (token, template1) != NULL)
      match = TRUE;
    break;
  }
  return (match);
}

/* ---------------------------------------------------------------------- */
int
extract_bracket (char **string, char *bracket_string)
/* ---------------------------------------------------------------------- */
{
  char *ptr, *ptr1;

  if ((ptr = strstr (*string, "{")) == NULL)
    return (FALSE);
  strcpy (bracket_string, ptr);
  if ((ptr1 = strstr (bracket_string, "}")) == NULL)
  {
    sprintf (error_string,
	     "No matching bracket (}) in isotope template string %s",
	     *string);
    error_msg (error_string, CONTINUE);
    input_error++;
    return (FALSE);
  }
  ptr1[1] = '\0';
  *string = strstr (*string, "}");
  *string += 1;
  return (TRUE);
}

/* ---------------------------------------------------------------------- */
LDBLE
surf_total (char *total_name, char *surface_name)
/* ---------------------------------------------------------------------- */
{
/*
 *   Provides total moles in LDBLE layer
 */
  int i, j, k;
  char name[MAX_LENGTH], token[MAX_LENGTH];
  char surface_name_local[MAX_LENGTH];
  char *ptr;

  if (use.surface_ptr == NULL)
    return (0);

/*
 *   Find surface...
 */
  for (j = 0; j < count_unknowns; j++)
  {
    if (x[j]->type != SURFACE)
      continue;
    strcpy (token, x[j]->master[0]->elt->name);
    replace ("_", " ", token);
    ptr = token;
    copy_token (name, &ptr, &k);
    if (surface_name != NULL)
    {
      if (strcmp (name, surface_name) == 0)
	break;
    }
    else
    {
      break;
    }
  }
  if (j >= count_unknowns)
    return (0);
  strcpy (surface_name_local, name);
/*
 *   find total moles of each element in diffuse layer...
 */
  count_elts = 0;
  paren_count = 0;
  for (j = 0; j < count_s_x; j++)
  {
    if (s_x[j]->type != SURF)
      continue;
    for (i = 0; s_x[j]->next_elt[i].elt != NULL; i++)
    {
      strcpy (token, s_x[j]->next_elt[i].elt->name);
      replace ("_", " ", token);
      ptr = token;
      copy_token (name, &ptr, &k);
      if (strcmp (name, surface_name_local) == 0)
      {
/*
 *   Accumulate elements in diffuse layer
 */
	add_elt_list (s_x[j]->next_elt, s_x[j]->moles);
	break;
      }
    }
  }
  if (count_elts > 0)
  {
    qsort (elt_list, (size_t) count_elts,
	   (size_t) sizeof (struct elt_list), elt_list_compare);
    elt_list_combine ();
  }
/*
 *   Return totals
 */
  for (j = 0; j < count_elts; j++)
  {
    if (strcmp (elt_list[j].elt->name, total_name) == 0)
    {
      return ((LDBLE) elt_list[j].coef);
    }
  }
  return (0);
}

/* ---------------------------------------------------------------------- */
LDBLE
total (char *total_name)
/* ---------------------------------------------------------------------- */
{
  struct master *master_ptr;
  LDBLE t;
  int i;

  if (strcmp (total_name, "H") == 0)
  {
    return (total_h_x / mass_water_aq_x);
  }
  if (strcmp (total_name, "O") == 0)
  {
    return (total_o_x / mass_water_aq_x);
  }
  master_ptr = master_bsearch (total_name);
  t = 0.0;
  if (master_ptr == NULL)
  {
    if (strcmp_nocase (total_name, "water") == 0)
    {
      return (mass_water_aq_x);
    }
/*
    sprintf (error_string, "Can not find definition for master species, %s.",
	     total_name);
    warning_msg (error_string);
*/
/*
 *  Primary master species
 */
  }
  else if (master_ptr->primary == TRUE)
  {
    /*
     *  Not a redox element
     */
    if (master_ptr->s->secondary == NULL)
    {
      t = master_ptr->total / mass_water_aq_x;
      /*
       * Redox element, need to sum totals of all redox states
       */
    }
    else
    {
      t = 0;
      for (i = master_ptr->number + 1;
	   (i < count_master && master[i]->elt->primary == master_ptr); i++)
      {
	t += master[i]->total / mass_water_aq_x;
      }
    }
/*
 *  Secondary master species
 */
  }
  else
  {
    t = master_ptr->total / mass_water_aq_x;
  }
  return (t);
}

/* ---------------------------------------------------------------------- */
LDBLE
system_total (char *total_name, LDBLE * count, char ***names, char ***types,
	      LDBLE ** moles)
/* ---------------------------------------------------------------------- */
{
/*
 *   Provides total moles in system and lists of species/phases in sort order
 */
  int i;

  sys_tot = 0;
  count_sys = 0;
  max_sys = 100;
  space ((void **) ((void *) &sys), INIT, &max_sys,
	 sizeof (struct system_species));
  if (strcmp_nocase (total_name, "elements") == 0)
  {
    system_total_elements ();
  }
  else if (strcmp_nocase (total_name, "phases") == 0)
  {
    system_total_si ();
  }
  else if (strcmp_nocase (total_name, "aq") == 0)
  {
    system_total_aq ();
  }
  else if (strcmp_nocase (total_name, "ex") == 0)
  {
    system_total_ex ();
  }
  else if (strcmp_nocase (total_name, "surf") == 0)
  {
    system_total_surf ();
  }
  else if (strcmp_nocase (total_name, "s_s") == 0)
  {
    system_total_s_s ();
  }
  else if (strcmp_nocase (total_name, "gas") == 0)
  {
    system_total_gas ();
  }
  else
  {
    if (strstr (total_name, "(") == NULL)
    {
      system_total_elt (total_name);
    }
    else
    {
      system_total_elt_secondary (total_name);
    }
  }
  /*
   *   Sort system species
   */
  if (count_sys > 1)
  {
    qsort (sys, (size_t) count_sys,
	   (size_t) sizeof (struct system_species), system_species_compare);
  }
  /*
   * malloc space
   */
  *names = (char **) PHRQ_malloc ((size_t) (count_sys + 1) * sizeof (char *));
  if (names == NULL)
    malloc_error ();
  *types = (char **) PHRQ_malloc ((size_t) (count_sys + 1) * sizeof (char *));
  if (types == NULL)
    malloc_error ();
  *moles = (LDBLE *) PHRQ_malloc ((size_t) (count_sys + 1) * sizeof (LDBLE));
  if (moles == NULL)
    malloc_error ();

  (*names)[0] = NULL;
  (*types)[0] = NULL;
  (*moles)[0] = 0;
  for (i = 0; i < count_sys; i++)
  {
    /*              output_msg(OUTPUT_MESSAGE, "%20s\t%10s\t%e\n", sys[i].name, sys[i].type, sys[i].moles); */
    (*names)[i + 1] = sys[i].name;
    (*types)[i + 1] = sys[i].type;
    (*moles)[i + 1] = sys[i].moles;
  }
  *count = (LDBLE) count_sys;
  PHRQ_free (sys);
  if (strcmp_nocase (total_name, "elements") == 0)
  {
    sys_tot = 0;;
    for (i = 0; i < count_sys; i++)
    {
      if (strcmp (sys[i].type, "dis") == 0 &&
	  strstr (sys[i].name, "(") == NULL &&
	  strcmp (sys[i].name, "H") != 0 && strcmp (sys[i].name, "O") != 0)
      {
	sys_tot += sys[i].moles;
      }
    }
  }
  return (sys_tot);
}

/* ---------------------------------------------------------------------- */
int
system_total_elements (void)
/* ---------------------------------------------------------------------- */
{
  int i, j;
  LDBLE t;
  char name[MAX_LENGTH];
  struct master *master_ptr;

  /*
   * Include H and O
   */
  sys[count_sys].name = string_duplicate ("H");
  sys[count_sys].moles = total_h_x;
  sys_tot += sys[count_sys].moles;
  sys[count_sys].type = string_duplicate ("dis");
  count_sys++;
  space ((void **) ((void *) &sys), count_sys, &max_sys,
	 sizeof (struct system_species));
  sys[count_sys].name = string_duplicate ("O");
  sys[count_sys].moles = total_o_x;
  sys_tot += sys[count_sys].moles;
  sys[count_sys].type = string_duplicate ("dis");
  count_sys++;
  space ((void **) ((void *) &sys), count_sys, &max_sys,
	 sizeof (struct system_species));
  /*
   * Include H(1) and O(-2)
   */
  sys[count_sys].name = string_duplicate ("H(1)");
  sys[count_sys].moles = solution_sum_secondary ("H(1)");
  sys_tot += sys[count_sys].moles;
  sys[count_sys].type = string_duplicate ("dis");
  count_sys++;
  space ((void **) ((void *) &sys), count_sys, &max_sys,
	 sizeof (struct system_species));
  sys[count_sys].name = string_duplicate ("O(-2)");
  sys[count_sys].moles = solution_sum_secondary ("O(-2)");
  sys_tot += sys[count_sys].moles;
  sys[count_sys].type = string_duplicate ("dis");
  count_sys++;
  space ((void **) ((void *) &sys), count_sys, &max_sys,
	 sizeof (struct system_species));

  for (i = 0; i < count_master; i++)
  {
    master_ptr = master[i];
    if (master_ptr->primary == TRUE && master_ptr->total_primary <= 0)
      continue;
    if (master_ptr->in == FALSE
	&& (master_ptr->primary == FALSE || master_ptr->total_primary == 0))
      continue;
    /*
     *  H and O
     */
    if (master_ptr->s == s_hplus)
    {
      continue;
    }
    else if (master_ptr->s == s_h2o)
    {
      continue;
    }
    if (master_ptr->primary == TRUE)
    {
      if (master_ptr->total_primary > 0)
      {
	t = master_ptr->total_primary;
	/*
	 *  Not a redox element
	 */
      }
      else if (master_ptr->s->secondary == NULL)
      {
	t = master_ptr->total;
	/*
	 * Redox element, need to sum totals of all redox states
	 */
      }
      else
      {
	t = 0;
	for (j = master_ptr->number + 1;
	     master[j]->elt->primary == master_ptr; j++)
	{
	  t += master[j]->total;
	}
      }
      /*
       *  Secondary master species
       */
    }
    else
    {
      t = master_ptr->total;
    }
    strcpy (name, master[i]->elt->name);
    sys[count_sys].name = string_duplicate (name);
    sys[count_sys].moles = t;
    sys_tot += sys[count_sys].moles;
    if (master[i]->s->type <= SOLID)
    {
      sys[count_sys].type = string_duplicate ("dis");
    }
    else if (master[i]->s->type == EX)
    {
      sys[count_sys].type = string_duplicate ("ex");
    }
    else if (master[i]->s->type == SURF || master[i]->s->type == SURF_PSI)
    {
      sys[count_sys].type = string_duplicate ("surf");
    }
    count_sys++;
    space ((void **) ((void *) &sys), count_sys, &max_sys,
	   sizeof (struct system_species));

  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
system_total_si (void)
/* ---------------------------------------------------------------------- */
{
  int i;
  LDBLE si, iap;
  struct rxn_token *rxn_ptr;
  char name[MAX_LENGTH];

  for (i = 0; i < count_phases; i++)
  {
    if (phases[i]->in == FALSE || phases[i]->type != SOLID)
      continue;
/*
 *   Print saturation index
 */
    iap = 0.0;
    for (rxn_ptr = phases[i]->rxn_x->token + 1; rxn_ptr->s != NULL; rxn_ptr++)
    {
      iap += rxn_ptr->s->la * rxn_ptr->coef;
    }
    si = -phases[i]->lk + iap;
    /* output_msg(OUTPUT_MESSAGE,"\t%-15s%7.2f%8.2f%8.2f  %s\n",phases[i]->name, si, iap, (LDBLE) phases[i]->lk, phases[i]->formula); */
    strcpy (name, phases[i]->name);
    sys[count_sys].name = string_duplicate (name);
    sys[count_sys].moles = si;
    if (si > sys_tot)
      sys_tot = si;
    sys[count_sys].type = string_duplicate ("phase");
    count_sys++;
    space ((void **) ((void *) &sys), count_sys, &max_sys,
	   sizeof (struct system_species));
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
system_total_aq (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Provides total moles in system and lists of species/phases in sort order
 */
  int i;
/*
 *   find total moles in aq, surface, and exchange
 */
  for (i = 0; i < count_s_x; i++)
  {
    if (s_x[i]->type != AQ)
      continue;
    sys[count_sys].name = string_duplicate (s_x[i]->name);
    sys[count_sys].moles = s_x[i]->moles;
    sys_tot += sys[count_sys].moles;
    sys[count_sys].type = string_duplicate ("aq");
    count_sys++;
    space ((void **) ((void *) &sys), count_sys, &max_sys,
	   sizeof (struct system_species));
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
system_total_ex (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Provides total moles in system and lists of species/phases in sort order
 */
  int i;
/*
 *   find total moles in aq, surface, and exchange
 */
  for (i = 0; i < count_s_x; i++)
  {
    if (s_x[i]->type != EX)
      continue;
    if (s_x[i]->primary != NULL)
      continue;
    sys[count_sys].name = string_duplicate (s_x[i]->name);
    sys[count_sys].moles = s_x[i]->moles;
    sys_tot += sys[count_sys].moles;
    sys[count_sys].type = string_duplicate ("ex");
    count_sys++;
    space ((void **) ((void *) &sys), count_sys, &max_sys,
	   sizeof (struct system_species));
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
system_total_surf (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Provides total moles in system and lists of species/phases in sort order
 */
  int i;
/*
 *   find total moles in aq, surface, and exchange
 */
  for (i = 0; i < count_s_x; i++)
  {
    if (s_x[i]->type != SURF)
      continue;
    sys[count_sys].name = string_duplicate (s_x[i]->name);
    sys[count_sys].moles = s_x[i]->moles;
    sys_tot += sys[count_sys].moles;
    sys[count_sys].type = string_duplicate ("surf");
    count_sys++;
    space ((void **) ((void *) &sys), count_sys, &max_sys,
	   sizeof (struct system_species));
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
system_total_gas (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Provides total moles in system and lists of species/phases in sort order
 */
  int i;

/*
 *   find total in gas phase
 */
  if (use.gas_phase_ptr == NULL)
    return (OK);
  for (i = 0; i < use.gas_phase_ptr->count_comps; i++)
  {
    sys[count_sys].name =
      string_duplicate (use.gas_phase_ptr->comps[i].phase->name);
    sys[count_sys].moles = use.gas_phase_ptr->comps[i].phase->moles_x;
    sys_tot += sys[count_sys].moles;
    sys[count_sys].type = string_duplicate ("gas");
    count_sys++;
    space ((void **) ((void *) &sys), count_sys, &max_sys,
	   sizeof (struct system_species));
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
system_total_s_s (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Provides total moles in system and lists of species/phases in sort order
 */
  int i, k;

/*
 *  Solid solutions
 */
  if (use.s_s_assemblage_ptr == NULL)
    return (OK);
  for (k = 0; k < use.s_s_assemblage_ptr->count_s_s; k++)
  {
    for (i = 0; i < use.s_s_assemblage_ptr->s_s[k].count_comps; i++)
    {
      sys[count_sys].name =
	string_duplicate (use.s_s_assemblage_ptr->s_s[k].comps[i].phase->
			  name);
      sys[count_sys].moles = use.s_s_assemblage_ptr->s_s[k].comps[i].moles;
      sys_tot += sys[count_sys].moles;
      sys[count_sys].type = string_duplicate ("s_s");
      count_sys++;
      space ((void **) ((void *) &sys), count_sys, &max_sys,
	     sizeof (struct system_species));
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
system_total_elt (char *total_name)
/* ---------------------------------------------------------------------- */
{
/*
 *   Provides total moles in system and lists of species/phases in sort order
 */
  int i, j, k, count_g;
  LDBLE molality, moles_excess, moles_surface, mass_water_surface;
  char name[MAX_LENGTH];

/*
 *   find total moles in aq, surface, and exchange
 */
  for (i = 0; i < count_s_x; i++)
  {
    count_elts = 0;
    paren_count = 0;
    add_elt_list (s_x[i]->next_elt, s_x[i]->moles);

    if (count_elts > 0)
    {
      qsort (elt_list, (size_t) count_elts,
	     (size_t) sizeof (struct elt_list), elt_list_compare);
      elt_list_combine ();
    }
    /*
     *   Look for element
     */
    for (j = 0; j < count_elts; j++)
    {
      if (strcmp (elt_list[j].elt->name, total_name) == 0)
      {
	sys[count_sys].name = string_duplicate (s_x[i]->name);
	sys[count_sys].moles = elt_list[j].coef;
	sys_tot += sys[count_sys].moles;
	if (s_x[i]->type == AQ)
	{
	  sys[count_sys].type = string_duplicate ("aq");
	}
	else if (s_x[i]->type == EX)
	{
	  sys[count_sys].type = string_duplicate ("ex");
	}
	else if (s_x[i]->type == SURF)
	{
	  sys[count_sys].type = string_duplicate ("surf");
	}
	else if (s_x[i]->type == HPLUS)
	{
	  sys[count_sys].type = string_duplicate ("aq");
	  /* sys[count_sys].moles = total_h_x; */
	}
	else if (s_x[i]->type == H2O)
	{
	  sys[count_sys].type = string_duplicate ("aq");
	  /* sys[count_sys].moles = total_o_x; */
	}
	else
	{
	  error_msg ("System_total", STOP);
	}
	count_sys++;
	space ((void **) ((void *) &sys), count_sys, &max_sys,
	       sizeof (struct system_species));
	break;
      }
    }
  }
  if (use.surface_ptr != NULL && dl_type_x != NO_DL)
  {
    /*
     *   Find position of component in surface charge data
     */
    i = -1;
    for (k = 0; k < count_unknowns; k++)
    {
      if (x[k]->type != SURFACE_CB)
	continue;
      i++;
      /*
       *   Loop through all surface components, calculate each H2O surface (diffuse layer),
       *   H2O aq, and H2O bulk (diffuse layers plus aqueous).
       */
      mass_water_surface = x[k]->surface_charge->mass_water;
      count_elts = 0;
      paren_count = 0;
      for (j = 0; j < count_s_x; j++)
      {
	if (s_x[j]->type > HPLUS)
	  continue;
	molality = under (s_x[j]->lm);
	count_g = s_x[j]->diff_layer[i].count_g;
	moles_excess =
	  mass_water_aq_x * molality * x[k]->surface_charge->g[count_g].g;
	moles_surface = mass_water_surface * molality + moles_excess;
	/*
	 *   Accumulate elements in diffuse layer
	 */
	add_elt_list (s_x[j]->next_elt, moles_surface);
      }
      if (count_elts > 0)
      {
	qsort (elt_list, (size_t) count_elts,
	       (size_t) sizeof (struct elt_list), elt_list_compare);
	elt_list_combine ();
      }
      /*
       *   Print totals
       */
      for (j = 0; j < count_elts; j++)
      {
	if (strcmp (elt_list[j].elt->name, total_name) == 0)
	{
	  strcpy (name, x[k]->master[0]->elt->name);
	  replace ("_psi", "", name);
	  sys[count_sys].name = string_duplicate (name);
	  sys[count_sys].moles = elt_list[j].coef;
	  sys_tot += sys[count_sys].moles;
	  sys[count_sys].type = string_duplicate ("diff");
	  count_sys++;
	  space ((void **) ((void *) &sys), count_sys, &max_sys,
		 sizeof (struct system_species));
	  break;
	}
      }
    }
  }
/*
 *   find total moles in mineral phases
 */
  if (use.pp_assemblage_in == TRUE && use.pp_assemblage_ptr != NULL)
  {
    for (i = 0; i < count_unknowns; i++)
    {
      if (x[i]->type != PP)
	continue;
      if (x[i]->pure_phase->add_formula != NULL)
	continue;
      count_elts = 0;
      paren_count = 0;
      add_elt_list (x[i]->pure_phase->phase->next_elt, x[i]->moles);
      if (count_elts > 0)
      {
	qsort (elt_list, (size_t) count_elts,
	       (size_t) sizeof (struct elt_list), elt_list_compare);
	elt_list_combine ();
      }
      for (j = 0; j < count_elts; j++)
      {
	if (strcmp (elt_list[j].elt->name, total_name) == 0)
	{
	  sys[count_sys].name = string_duplicate (x[i]->pure_phase->name);
	  sys[count_sys].moles = elt_list[j].coef;
	  sys_tot += sys[count_sys].moles;
	  sys[count_sys].type = string_duplicate ("equi");
	  count_sys++;
	  space ((void **) ((void *) &sys), count_sys, &max_sys,
		 sizeof (struct system_species));
	  break;
	}
      }
    }
  }
/*
 *  Solid solutions
 */
  if (use.s_s_assemblage_ptr != NULL)
  {
    for (k = 0; k < use.s_s_assemblage_ptr->count_s_s; k++)
    {
      if (use.s_s_assemblage_ptr->s_s[k].s_s_in == TRUE)
      {
	for (i = 0; i < use.s_s_assemblage_ptr->s_s[k].count_comps; i++)
	{
	  count_elts = 0;
	  paren_count = 0;
	  add_elt_list (use.s_s_assemblage_ptr->s_s[k].comps[i].phase->
			next_elt,
			use.s_s_assemblage_ptr->s_s[k].comps[i].moles);
	  if (count_elts > 0)
	  {
	    qsort (elt_list, (size_t) count_elts,
		   (size_t) sizeof (struct elt_list), elt_list_compare);
	    elt_list_combine ();
	  }
	  for (j = 0; j < count_elts; j++)
	  {
	    if (strcmp (elt_list[j].elt->name, total_name) == 0)
	    {
	      sys[count_sys].name =
		string_duplicate (use.s_s_assemblage_ptr->s_s[k].comps[i].
				  phase->name);
	      sys[count_sys].moles = elt_list[j].coef;
	      sys_tot += sys[count_sys].moles;
	      sys[count_sys].type = string_duplicate ("s_s");
	      count_sys++;
	      space ((void **) ((void *) &sys), count_sys, &max_sys,
		     sizeof (struct system_species));
	      break;
	    }
	  }
	}
      }
    }
  }
/*
 *   find total in gas phase
 */
  if (use.gas_phase_ptr != NULL)
  {
    for (i = 0; i < use.gas_phase_ptr->count_comps; i++)
    {
      if (use.gas_phase_ptr->comps[i].phase->in == TRUE)
      {
	count_elts = 0;
	paren_count = 0;
	add_elt_list (use.gas_phase_ptr->comps[i].phase->next_elt,
		      use.gas_phase_ptr->comps[i].phase->moles_x);

	if (count_elts > 0)
	{
	  qsort (elt_list, (size_t) count_elts,
		 (size_t) sizeof (struct elt_list), elt_list_compare);
	  elt_list_combine ();
	}
	/*
	 *   Look for element
	 */
	for (j = 0; j < count_elts; j++)
	{
	  if (strcmp (elt_list[j].elt->name, total_name) == 0)
	  {
	    sys[count_sys].name =
	      string_duplicate (use.gas_phase_ptr->comps[i].phase->name);
	    sys[count_sys].moles = elt_list[j].coef;
	    sys_tot += sys[count_sys].moles;
	    sys[count_sys].type = string_duplicate ("gas");
	    count_sys++;
	    space ((void **) ((void *) &sys), count_sys, &max_sys,
		   sizeof (struct system_species));
	    break;
	  }
	}
      }
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
system_total_elt_secondary (char *total_name)
/* ---------------------------------------------------------------------- */
{
/*
 *   Provides total moles in system and lists of species/phases in sort order
 */
  int i, j, k, l, count_g;
  LDBLE molality, moles_excess, moles_surface, mass_water_surface, sum, coef;
  char name[MAX_LENGTH];
/*
 *   find total moles in aq, surface, and exchange
 */
  for (i = 0; i < count_s_x; i++)
  {
    count_elts = 0;
    paren_count = 0;
    if (s_x[i]->next_secondary != NULL)
    {
      add_elt_list (s_x[i]->next_secondary, s_x[i]->moles);
    }
    else
    {
      add_elt_list (s_x[i]->next_sys_total, s_x[i]->moles);
    }
    if (count_elts > 0)
    {
      qsort (elt_list, (size_t) count_elts,
	     (size_t) sizeof (struct elt_list), elt_list_compare);
      elt_list_combine ();
    }
    /* debug
       output_msg(OUTPUT_MESSAGE, "%s\n", s_x[i]->name);
       for ( j=0; j < count_elts; j++ ) {
       output_msg(OUTPUT_MESSAGE, "\t%10s\t%g\n", elt_list[j].elt->name, elt_list[j].coef);
       }
     */
    /*
     *   Look for element
     */
    for (j = 0; j < count_elts; j++)
    {
      if (strcmp (elt_list[j].elt->name, total_name) == 0)
      {
	sys[count_sys].name = string_duplicate (s_x[i]->name);
	sys[count_sys].moles = elt_list[j].coef;
	sys_tot += sys[count_sys].moles;
	if (s_x[i]->type == AQ)
	{
	  sys[count_sys].type = string_duplicate ("aq");
	}
	else if (s_x[i]->type == EX)
	{
	  sys[count_sys].type = string_duplicate ("ex");
	}
	else if (s_x[i]->type == SURF)
	{
	  sys[count_sys].type = string_duplicate ("surf");
	}
	else if (s_x[i]->type == HPLUS)
	{
	  sys[count_sys].type = string_duplicate ("aq");
	  /* sys[count_sys].moles = total_h_x; */
	}
	else if (s_x[i]->type == H2O)
	{
	  sys[count_sys].type = string_duplicate ("aq");
	  /* sys[count_sys].moles = total_o_x; */
	}
	else
	{
	  error_msg ("System_total", STOP);
	}
	count_sys++;
	space ((void **) ((void *) &sys), count_sys, &max_sys,
	       sizeof (struct system_species));
	break;
      }
    }
  }
  if (use.surface_ptr != NULL && dl_type_x != NO_DL)
  {
    /*
     *   Find position of component in surface charge data
     */
    i = -1;
    for (k = 0; k < count_unknowns; k++)
    {
      if (x[k]->type != SURFACE_CB)
	continue;
      i++;
      /*
       *   Loop through all surface components, calculate each H2O surface (diffuse layer),
       *   H2O aq, and H2O bulk (diffuse layers plus aqueous).
       */
      mass_water_surface = x[k]->surface_charge->mass_water;
      sum = 0;
      for (j = 0; j < count_s_x; j++)
      {
	count_elts = 0;
	paren_count = 0;
	if (s_x[i]->next_secondary != NULL)
	{
	  add_elt_list (s_x[i]->next_secondary, 1);
	}
	else
	{
	  add_elt_list (s_x[i]->next_sys_total, 1);
	}
	for (l = 0; l < count_elts; l++)
	{
	  if (strcmp (elt_list[l].elt->name, total_name) == 0)
	  {
	    coef = elt_list[l].coef;
	    if (s_x[j]->type > H2O)
	      continue;
	    molality = under (s_x[j]->lm);
	    count_g = s_x[j]->diff_layer[i].count_g;
	    moles_excess =
	      mass_water_aq_x * molality * x[k]->surface_charge->g[count_g].g;
	    moles_surface = mass_water_surface * molality + moles_excess;
	    sum += moles_surface * coef;
	    break;
	  }
	}
	if (l >= count_elts)
	  continue;
	strcpy (name, x[k]->master[0]->elt->name);
	replace ("_psi", "", name);
	sys[count_sys].name = string_duplicate (name);
	sys[count_sys].moles = sum;
	sys_tot += sys[count_sys].moles;
	sys[count_sys].type = string_duplicate ("diff");
	count_sys++;
	space ((void **) ((void *) &sys), count_sys, &max_sys,
	       sizeof (struct system_species));
	break;
      }
    }
  }
/*
 *   find total moles in mineral phases
 */
  if (use.pp_assemblage_in == TRUE && use.pp_assemblage_ptr != NULL)
  {
    for (i = 0; i < count_unknowns; i++)
    {
      if (x[i]->type != PP)
	continue;
      if (x[i]->pure_phase->add_formula != NULL)
	continue;
      count_elts = 0;
      paren_count = 0;
      add_elt_list (x[i]->pure_phase->phase->next_sys_total, x[i]->moles);
      if (count_elts > 0)
      {
	qsort (elt_list, (size_t) count_elts,
	       (size_t) sizeof (struct elt_list), elt_list_compare);
	elt_list_combine ();
      }
      /* debug
         output_msg(OUTPUT_MESSAGE, "%s\n", x[i]->pure_phase->phase->name);
         for ( j=0; j < count_elts; j++ ) {
         output_msg(OUTPUT_MESSAGE, "\t%10s\t%g\n", elt_list[j].elt->name, elt_list[j].coef);
         }
       */
      for (j = 0; j < count_elts; j++)
      {
	if (strcmp (elt_list[j].elt->name, total_name) == 0)
	{
	  sys[count_sys].name = string_duplicate (x[i]->pure_phase->name);
	  sys[count_sys].moles = elt_list[j].coef;
	  sys_tot += sys[count_sys].moles;
	  sys[count_sys].type = string_duplicate ("equi");
	  count_sys++;
	  space ((void **) ((void *) &sys), count_sys, &max_sys,
		 sizeof (struct system_species));
	  break;
	}
      }
    }
  }
/*
 *  Solid solutions
 */
  if (use.s_s_assemblage_ptr != NULL)
  {
    for (k = 0; k < use.s_s_assemblage_ptr->count_s_s; k++)
    {
      if (use.s_s_assemblage_ptr->s_s[k].s_s_in == TRUE)
      {
	for (i = 0; i < use.s_s_assemblage_ptr->s_s[k].count_comps; i++)
	{
	  count_elts = 0;
	  paren_count = 0;
	  add_elt_list (use.s_s_assemblage_ptr->s_s[k].comps[i].phase->
			next_sys_total,
			use.s_s_assemblage_ptr->s_s[k].comps[i].moles);
	  if (count_elts > 0)
	  {
	    qsort (elt_list, (size_t) count_elts,
		   (size_t) sizeof (struct elt_list), elt_list_compare);
	    elt_list_combine ();
	  }
	  for (j = 0; j < count_elts; j++)
	  {
	    if (strcmp (elt_list[j].elt->name, total_name) == 0)
	    {
	      sys[count_sys].name =
		string_duplicate (use.s_s_assemblage_ptr->s_s[k].comps[i].
				  phase->name);
	      sys[count_sys].moles = elt_list[j].coef;
	      sys_tot += sys[count_sys].moles;
	      sys[count_sys].type = string_duplicate ("s_s");
	      count_sys++;
	      space ((void **) ((void *) &sys), count_sys, &max_sys,
		     sizeof (struct system_species));
	      break;
	    }
	  }
	}
      }
    }
  }
/*
 *   find total in gas phase
 */
  if (use.gas_phase_ptr != NULL)
  {
    for (i = 0; i < use.gas_phase_ptr->count_comps; i++)
    {
      if (use.gas_phase_ptr->comps[i].phase->in == TRUE)
      {
	count_elts = 0;
	paren_count = 0;
	add_elt_list (use.gas_phase_ptr->comps[i].phase->next_sys_total,
		      use.gas_phase_ptr->comps[i].phase->moles_x);

	if (count_elts > 0)
	{
	  qsort (elt_list, (size_t) count_elts,
		 (size_t) sizeof (struct elt_list), elt_list_compare);
	  elt_list_combine ();
	}
	/*
	 *   Look for element
	 */
	for (j = 0; j < count_elts; j++)
	{
	  if (strcmp (elt_list[j].elt->name, total_name) == 0)
	  {
	    sys[count_sys].name =
	      string_duplicate (use.gas_phase_ptr->comps[i].phase->name);
	    sys[count_sys].moles = elt_list[j].coef;
	    sys_tot += sys[count_sys].moles;
	    sys[count_sys].type = string_duplicate ("gas");
	    count_sys++;
	    space ((void **) ((void *) &sys), count_sys, &max_sys,
		   sizeof (struct system_species));
	    break;
	  }
	}
      }
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
LDBLE
solution_sum_secondary (char *total_name)
/* ---------------------------------------------------------------------- */
{
/*
 *   Provides total moles in system and lists of species/phases in sort order
 */
  int i, j;
  LDBLE sum;
/*
 *   find total moles in aq, surface, and exchange
 */
  sum = 0;
  for (i = 0; i < count_s_x; i++)
  {
    if (s_x[i]->type > H2O)
      continue;
    count_elts = 0;
    paren_count = 0;
    if (s_x[i]->next_secondary != NULL)
    {
      add_elt_list (s_x[i]->next_secondary, s_x[i]->moles);
    }
    else
    {
      add_elt_list (s_x[i]->next_sys_total, s_x[i]->moles);
    }
    if (count_elts > 0)
    {
      qsort (elt_list, (size_t) count_elts,
	     (size_t) sizeof (struct elt_list), elt_list_compare);
      elt_list_combine ();
    }
    /*
     *   Look for element
     */
    for (j = 0; j < count_elts; j++)
    {
      if (strcmp (elt_list[j].elt->name, total_name) == 0)
      {
	sum += elt_list[j].coef;
	break;
      }
    }
  }
  return (sum);
}

/* ---------------------------------------------------------------------- */
int
system_species_compare (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
  const struct system_species *a, *b;

  a = (const struct system_species *) ptr1;
  b = (const struct system_species *) ptr2;
  if (a->moles < b->moles)
    return (1);
  if (a->moles > b->moles)
    return (-1);
  return (0);
}

/* ---------------------------------------------------------------------- */
int
system_total_solids (struct exchange *exchange_ptr,
		     struct pp_assemblage *pp_assemblage_ptr,
		     struct gas_phase *gas_phase_ptr,
		     struct s_s_assemblage *s_s_assemblage_ptr,
		     struct surface *surface_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Provides total moles in solid phases
 */
  int i, j;

  count_elts = 0;
  paren_count = 0;
/*
 *   find total moles in exchanger
 */
  if (exchange_ptr != NULL)
  {
    for (i = 0; i < exchange_ptr->count_comps; i++)
    {
      add_elt_list (exchange_ptr->comps[i].totals, 1.0);
    }
  }
  if (surface_ptr != NULL)
  {
    for (i = 0; i < surface_ptr->count_comps; i++)
    {
      add_elt_list (surface_ptr->comps[i].totals, 1.0);
    }
  }
  if (s_s_assemblage_ptr != NULL)
  {
    for (i = 0; i < s_s_assemblage_ptr->count_s_s; i++)
    {
      for (j = 0; j < s_s_assemblage_ptr->s_s[i].count_comps; j++)
      {
	add_elt_list (s_s_assemblage_ptr->s_s[i].comps[j].phase->next_elt,
		      s_s_assemblage_ptr->s_s[i].comps[j].moles);
      }
    }
  }
  if (gas_phase_ptr != NULL)
  {
    for (i = 0; i < gas_phase_ptr->count_comps; i++)
    {
      add_elt_list (gas_phase_ptr->comps[i].phase->next_elt,
		    gas_phase_ptr->comps[i].moles);
    }
  }
  if (pp_assemblage_ptr != NULL)
  {
    for (i = 0; i < pp_assemblage_ptr->count_comps; i++)
    {
      add_elt_list (pp_assemblage_ptr->pure_phases[i].phase->next_elt,
		    pp_assemblage_ptr->pure_phases[i].moles);
    }
  }

  if (count_elts > 0)
  {
    qsort (elt_list, (size_t) count_elts,
	   (size_t) sizeof (struct elt_list), elt_list_compare);
    elt_list_combine ();
  }
  return (OK);
}
