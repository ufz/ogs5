#define EXTERNAL extern
#include "global.h"
#include "phqalloc.h"
#include "output.h"
#include "phrqproto.h"
#ifdef PHREEQC_CPP
extern int read_solution_raw (void);
extern int read_exchange_raw (void);
extern int read_surface_raw (void);
extern int read_equilibrium_phases_raw (void);
extern int read_kinetics_raw (void);
extern int read_solid_solutions_raw (void);
extern int read_gas_phase_raw (void);
extern int read_reaction_raw (void);
extern int read_mix_raw (void);
extern int read_temperature_raw (void);
#endif
static char const svnid[] = "$Id: read.c 4 2009-04-21 17:29:29Z delucia $";

#if defined(SWIG_SHARED_OBJ)
#define STATIC
#else
#define STATIC static
#endif
STATIC int add_psi_master_species (char *token);
STATIC int read_advection (void);
STATIC int read_analytical_expression_only (char *ptr, LDBLE * log_k);
STATIC int read_copy (void);
STATIC int read_debug (void);
STATIC int read_delta_h_only (char *ptr, LDBLE * delta_h,
			      DELTA_H_UNIT * units);
STATIC int read_llnl_aqueous_model_parameters (void);
STATIC int read_exchange (void);
STATIC int read_exchange_master_species (void);
STATIC int read_exchange_species (void);
STATIC int read_gas_phase (void);
STATIC int read_incremental_reactions (void);
STATIC int read_inverse (void);
STATIC int read_inv_balances (struct inverse *inverse_ptr, char *next_char);
STATIC int read_inv_isotopes (struct inverse *inverse_ptr, char *ptr);
STATIC int read_inv_phases (struct inverse *inverse_ptr, char *next_char);
STATIC int read_kinetics (void);
STATIC int read_line_doubles (char *next_char, LDBLE ** d, int *count_d,
			      int *count_alloc);
STATIC int read_lines_doubles (char *next_char, LDBLE ** d, int *count_d,
			       int *count_alloc, const char **opt_list,
			       int count_opt_list, int *opt);
STATIC LDBLE *read_list_doubles (char **ptr, int *count_doubles);
STATIC int *read_list_ints (char **ptr, int *count_ints, int positive);
STATIC int *read_list_t_f (char **ptr, int *count_ints);
STATIC int read_master_species (void);
STATIC int read_mix (void);
STATIC int read_named_logk (void);
STATIC int read_phases (void);
STATIC int read_print (void);
STATIC int read_pure_phases (void);
STATIC int read_rates (void);
STATIC int read_reaction (void);
STATIC int read_reaction_reactants (struct irrev *irrev_ptr);
STATIC int read_reaction_steps (struct irrev *irrev_ptr);
STATIC int read_solid_solutions (void);
STATIC int read_temperature (void);
STATIC int read_reaction_temps (struct temperature *temperature_ptr);
STATIC int read_save (void);
STATIC int read_selected_output (void);
STATIC int read_solution (void);
STATIC int read_species (void);
STATIC int read_surf (void);
STATIC int read_surface_master_species (void);
STATIC int read_surface_species (void);
STATIC int read_use (void);
STATIC int read_title (void);
STATIC int read_user_print (void);
STATIC int read_user_punch (void);
#ifdef PHREEQ98
STATIC int read_user_graph (void);
extern int connect_simulations, graph_initial_solutions;
/*extern*/ int shifts_as_points;
extern int chart_type;
extern int ShowChart;
extern int RowOffset, ColumnOffset;
#endif

extern int reading_database (void);
extern int check_line (const char *string, int allow_empty, int allow_eof,
		       int allow_keyword, int print);

#ifdef PHREEQ98
extern int copy_title (char *token_ptr, char **ptr, int *length);
extern int OpenCSVFile (char file_name[MAX_LENGTH]);
void GridHeadings (char *s, int i);
void SetAxisTitles (char *s, int i);
void SetAxisScale (char *a, int c, char *v, int l);
void SetChartTitle (char *s);
#endif

static LDBLE dummy;

/* ---------------------------------------------------------------------- */
int
read_input (void)
/* ---------------------------------------------------------------------- */
{
  int i, j, l;
  char *ptr;
  char token[2 * MAX_LENGTH];
  if (svnid == NULL)
    fprintf (stderr, " ");

  parse_error = 0;
  input_error = 0;
  next_keyword = 0;
  count_warnings = 0;
/*
 *  Initialize keyword flags
 */
  for (i = 0; i < NKEYS; i++)
  {
    keyword[i].keycount = 0;
  }
/*
 *  Initialize use and save pointers
 */
  use.solution_in = FALSE;
  use.solution_ptr = NULL;
  use.pp_assemblage_in = FALSE;
  use.pp_assemblage_ptr = NULL;
  use.mix_in = FALSE;
  use.mix_ptr = NULL;
  use.irrev_in = FALSE;
  use.irrev_ptr = NULL;
  use.kinetics_in = FALSE;
  use.kinetics_ptr = NULL;
  use.exchange_in = FALSE;
  use.exchange_ptr = NULL;
  use.surface_in = FALSE;
  use.surface_ptr = NULL;
  use.temperature_in = FALSE;
  use.temperature_ptr = NULL;
  use.gas_phase_in = FALSE;
  use.gas_phase_ptr = NULL;
  use.s_s_assemblage_in = FALSE;
  use.s_s_assemblage_ptr = NULL;
  use.trans_in = FALSE;
  use.advect_in = FALSE;

  save.solution = FALSE;
  save.mix = FALSE;
  save.irrev = FALSE;
  save.pp_assemblage = FALSE;
  save.exchange = FALSE;
  save.surface = FALSE;
  save.gas_phase = FALSE;
  save.s_s_assemblage = FALSE;
  title_x = (char *) free_check_null (title_x);

  while ((i =
	  check_line ("Subroutine Read", FALSE, TRUE, TRUE, TRUE)) != KEYWORD)
  {
    /* empty, eof, keyword, print */

    if (i == EOF)
      return (EOF);
    sprintf (error_string, "Unknown input, no keyword has been specified.");
    warning_msg (error_string);
  }
/*
	  0	  "eof"
	  1	  "end"
	  2	  "species"
	  3	  "master"
	  4	  "solution"
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
	  60      "solution_raw"
	  61      "exchange_raw"
	  62      "surface_raw"
	  63      "equilibrium_phases_raw"
	  64      "kinetics_raw"
	  65      "solid_solutions_raw"
	  66      "gas_phase_raw"
	  67      "reaction_raw"
	  68      "mix_raw"
	  69      "reaction_temperature_raw"
  */
  for (;;)
  {
    if (next_keyword >= 0)
    {
      /* keyword[next_keyword].keycount++; */
      if (next_keyword != 49 && !reading_database ())
	first_read_input = FALSE;
    }
    switch (next_keyword)
    {
    case -1:			/* Have not read line with keyword */
      sprintf (error_string, "Unknown input, no keyword has been specified.");
      warning_msg (error_string);
      while ((j =
	      check_line ("No keyword", FALSE, TRUE, TRUE, TRUE)) != KEYWORD
	     && j != EOF)
      {
	warning_msg (error_string);
      }
      break;
    case 0:			/* End encountered */
    case 1:			/* EOF encountered */
      goto END_OF_SIMULATION_INPUT;
    case 2:			/* Read aqueous model */
      keyword[2].keycount++;
      read_species ();
      break;
    case 3:			/* Read master species */
      keyword[3].keycount++;
      read_master_species ();
      break;
    case 4:			/* Read solution data */
      keyword[4].keycount++;
      read_solution ();
      solution_sort ();
      break;
    case 5:
      keyword[5].keycount++;
      read_phases ();
      break;
    case 6:
    case 26:
    case 27:
    case 28:
    case 29:
      keyword[6].keycount++;
      keyword[26].keycount++;
      keyword[27].keycount++;
      keyword[28].keycount++;
      keyword[29].keycount++;
      read_pure_phases ();
      break;
    case 7:
      keyword[7].keycount++;
      read_reaction ();
      break;
    case 8:
      keyword[8].keycount++;
      read_mix ();
      break;
    case 9:
      keyword[9].keycount++;
      read_use ();
      break;
    case 10:
      keyword[10].keycount++;
      read_save ();
      break;
    case 11:
      keyword[11].keycount++;
      read_exchange_species ();
      break;
    case 12:
      keyword[12].keycount++;
      read_exchange_master_species ();
      break;
    case 13:
      keyword[13].keycount++;
      read_exchange ();
      break;
    case 14:
      keyword[14].keycount++;
      read_surface_species ();
      break;
    case 15:
      keyword[15].keycount++;
      read_surface_master_species ();
      break;
    case 16:
      keyword[16].keycount++;
      read_surf ();
      break;
    case 17:
      keyword[17].keycount++;
      read_temperature ();
      break;
    case 18:
      keyword[18].keycount++;
      read_inverse ();
      break;
    case 19:
      keyword[19].keycount++;
      read_gas_phase ();
      break;
    case 20:
      keyword[20].keycount++;
      read_transport ();
      break;
    case 21:
    case 24:
      keyword[21].keycount++;
      keyword[24].keycount++;
      read_debug ();
      break;
    case 22:
    case 23:
    case 44:
    case 45:
      keyword[22].keycount++;
      keyword[23].keycount++;
      keyword[44].keycount++;
      keyword[45].keycount++;
      read_selected_output ();
      break;
    case 25:
      keyword[25].keycount++;
      read_print ();
      break;
    case 30:
    case 31:
      keyword[30].keycount++;
      keyword[31].keycount++;
      read_title ();
      break;
    case 32:
      keyword[32].keycount++;
      read_advection ();
      break;
    case 33:
      keyword[33].keycount++;
      read_kinetics ();
      break;
    case 34:
    case 35:
      keyword[34].keycount++;
      keyword[35].keycount++;
      read_incremental_reactions ();
      break;
    case 36:
      keyword[36].keycount++;
      read_rates ();
      break;
    case 37:
    case 42:
    case 43:
      keyword[37].keycount++;
      keyword[42].keycount++;
      keyword[43].keycount++;
      read_solution_spread ();
#ifdef SKIP
      solution_sort ();
#endif
      break;
    case 38:
      keyword[38].keycount++;
      read_user_print ();
      break;
    case 39:
      keyword[39].keycount++;
      read_user_punch ();
      break;
    case 40:
    case 41:
      keyword[40].keycount++;
      keyword[41].keycount++;
      read_solid_solutions ();
      break;
    case 46:
      keyword[46].keycount++;
#ifdef PHREEQ98
      read_user_graph ();
# else
      for (;;)
      {
	j = check_line ("Reading user_graph", FALSE, TRUE, TRUE, TRUE);
	if (j == EOF || j == KEYWORD)
	{
	  break;
	}
      }
#endif
      break;
    case 47:
    case 48:
      keyword[47].keycount++;
      keyword[48].keycount++;
      read_llnl_aqueous_model_parameters ();
      break;
    case 49:
      keyword[49].keycount++;
      if (reading_database ())
      {
	/* warning_msg("DATABASE is ignored in the database file."); */
      }
      else if (first_read_input == FALSE)
      {
	error_msg ("DATABASE must be the first keyword in the input file.",
		   CONTINUE);
	input_error++;
      }
      else
      {
	ptr = line;
	copy_token (token, &ptr, &l);
	user_database = string_duplicate (ptr);
	if (string_trim (user_database) == EMPTY)
	{
	  error_msg ("DATABASE file name is missing.", CONTINUE);
	  input_error++;
	  user_database = (char *) free_check_null (user_database);
	}
	first_read_input = FALSE;
      }
      j = check_line ("Reading after DATABASE", FALSE, TRUE, TRUE, TRUE);
      break;
    case 50:
    case 51:
    case 52:
    case 53:
      keyword[50].keycount++;
      keyword[51].keycount++;
      keyword[52].keycount++;
      keyword[53].keycount++;
      read_named_logk ();
      break;
    case 54:
      keyword[54].keycount++;
      read_isotopes ();
      break;
    case 55:
      keyword[55].keycount++;
      read_calculate_values ();
      break;
    case 56:
      keyword[56].keycount++;
      read_isotope_ratios ();
      break;
    case 57:
      keyword[57].keycount++;
      read_isotope_alphas ();
      break;
    case 58:
      keyword[58].keycount++;
      read_copy ();
      break;
    case 59:
      keyword[59].keycount++;
      read_pitzer ();
      break;
#ifdef PHREEQC_CPP
    case 60:
      keyword[60].keycount++;
      read_solution_raw ();
      break;
    case 61:
      keyword[61].keycount++;
      read_exchange_raw ();
      break;
    case 62:
      keyword[62].keycount++;
      read_surface_raw ();
      break;
    case 63:
      keyword[63].keycount++;
      read_equilibrium_phases_raw ();
      break;
    case 64:
      keyword[64].keycount++;
      read_kinetics_raw ();
      break;
    case 65:
      keyword[65].keycount++;
      read_solid_solutions_raw ();
      break;
    case 66:
      keyword[66].keycount++;
      read_gas_phase_raw ();
      break;
    case 67:
      keyword[67].keycount++;
      read_reaction_raw ();
      break;
    case 68:
      keyword[68].keycount++;
      read_mix_raw ();
      break;
    case 69:
      keyword[69].keycount++;
      read_temperature_raw ();
      break;
#endif
    }
  }
END_OF_SIMULATION_INPUT:
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
read_conc (int n, int count_mass_balance, char *str)
/* ---------------------------------------------------------------------- */
{
  int j, l;
  int alk;
  int count_redox_states;

  char *ptr;
  char token[MAX_LENGTH], token1[MAX_LENGTH];
/*
 *   Set defaults
 */
  /*
     solution[n]->totals[count_mass_balance].equation_name = NULL;
     solution[n]->totals[count_mass_balance].phase = NULL;
     solution[n]->totals[count_mass_balance].phase_si = 0.0;
     solution[n]->totals[count_mass_balance].units=NULL;
     solution[n]->totals[count_mass_balance].n_pe=-1;
     solution[n]->totals[count_mass_balance].as=NULL;
     solution[n]->totals[count_mass_balance].gfw= 0.0;
   */
  conc_init (&(solution[n]->totals[count_mass_balance]));

/*
 *   Remove space between "kg" and "solution" or "water" in units
 */
  replace ("Kg", "kg", str);
  replace ("KG", "kg", str);
  while (replace ("kg ", "kg", str) == TRUE);
  ptr = str;
/*
 *   Read master species list for mass balance equation
 */
  token1[0] = '\0';
  count_redox_states = 0;
  while (((j = copy_token (token, &ptr, &l)) == UPPER) ||
	 (token[0] == '[') ||
	 (strcmp_nocase_arg1 (token, "ph") == 0) ||
	 (strcmp_nocase_arg1 (token, "pe") == 0))
  {
    count_redox_states++;
    replace ("(+", "(", token);
    if (count_redox_states > 1)
      strcat (token1, " ");
    strcat (token1, token);
  }
  if (count_redox_states == 0)
  {
    input_error++;
    error_msg ("No element or master species given for concentration input.",
	       CONTINUE);
    return (ERROR);
  }
  solution[n]->totals[count_mass_balance].description = string_hsave (token1);
/*
 *   Determine if reading alkalinity, allow equivalents for units
 */
  str_tolower (token1);
  if (strstr (token1, "alk") == token1)
  {
    alk = TRUE;
  }
  else
  {
    alk = FALSE;
  }
/*
 *   Read concentration
 */

  j =
    sscanf (token, SCANFORMAT,
	    &(solution[n]->totals[count_mass_balance].input_conc));
  if (j == 0)
  {
    sprintf (error_string,
	     "Concentration data error for %s in solution input.", token1);
    error_msg (error_string, CONTINUE);
    return (ERROR);
  }
  if ((j = copy_token (token, &ptr, &l)) == EMPTY)
    return (OK);
/*
 *   Read optional data
 */
  strcpy (token1, token);
/*
 *   Check for units info
 */
  if (check_units (token1, alk, FALSE, solution[n]->units, FALSE) == OK)
  {
    if (check_units (token1, alk, FALSE, solution[n]->units, TRUE) == OK)
    {
      solution[n]->totals[count_mass_balance].units = string_hsave (token1);
      if ((j = copy_token (token, &ptr, &l)) == EMPTY)
	return (OK);
    }
    else
    {
      return (ERROR);
    }
  }
/*
 *   Check for "as" followed by formula to be used for gfw
 */
  strcpy (token1, token);
  str_tolower (token1);
  if (strcmp (token1, "as") == 0)
  {
    copy_token (token, &ptr, &l);
    solution[n]->totals[count_mass_balance].as = string_hsave (token);
    if ((j = copy_token (token, &ptr, &l)) == EMPTY)
      return (OK);
/*
 *   Check for "gfw" followed by gram formula weight
 */
  }
  else if (strcmp (token1, "gfw") == 0)
  {
    if (copy_token (token, &ptr, &l) != DIGIT)
    {
      error_msg ("Expecting gram formula weight.", CONTINUE);
      return (ERROR);
    }
    else
    {
      sscanf (token, SCANFORMAT,
	      &(solution[n]->totals[count_mass_balance].gfw));
      if ((j = copy_token (token, &ptr, &l)) == EMPTY)
	return (OK);
    }
  }
/*
 *   Check for redox couple for pe
 */
  if (strcmp_nocase_arg1 (token, "pe") == 0)
  {
    solution[n]->totals[count_mass_balance].n_pe =
      pe_data_store (&(solution[n]->pe), token);
    if ((j = copy_token (token, &ptr, &l)) == EMPTY)
      return (OK);
  }
  else if (strstr (token, "/") != NULL)
  {
    if (parse_couple (token) == OK)
    {
      solution[n]->totals[count_mass_balance].n_pe =
	pe_data_store (&(solution[n]->pe), token);
      if ((j = copy_token (token, &ptr, &l)) == EMPTY)
	return (OK);
    }
    else
    {
      return (ERROR);
    }
  }
/*
 *   Must have phase
 */
  solution[n]->totals[count_mass_balance].equation_name =
    string_hsave (token);
  if ((j = copy_token (token, &ptr, &l)) == EMPTY)
    return (OK);
/*
 *   Check for saturation index
 */
  j =
    sscanf (token, SCANFORMAT,
	    &(solution[n]->totals[count_mass_balance].phase_si));
  if (j != 1)
  {
    error_msg ("Expected saturation index.", CONTINUE);
    return (ERROR);
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
read_exchange_species (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Read data for exchange species, parse equations
 */
  int i;
  int association;
  char token[MAX_LENGTH];
  char *ptr;
  struct phase *phase_ptr;

  struct species *s_ptr;
  struct elt_list *next_elt;
  struct rxn_token *token_ptr;
  LDBLE exchange_coef;
  LDBLE offset;

  int return_value, opt, opt_save;
  char *next_char;
  const char *opt_list[] = {
    "no_check",			/* 0 */
    "check",			/* 1 */
    "mb",			/* 2 */
    "mass_balance",		/* 3 */
    "log_k",			/* 4 */
    "logk",			/* 5 */
    "delta_h",			/* 6 */
    "deltah",			/* 7 */
    "analytical_expression",	/* 8 */
    "a_e",			/* 9 */
    "ae",			/* 10 */
    "mole_balance",		/* 11 */
    "gamma",			/* 12 */
    "davies",			/* 13 */
    "offset",			/* 14 */
    "llnl_gamma",		/* 15 */
    "add_logk",			/* 16 */
    "add_log_k",		/* 17 */
    "add_constant"		/* 18 */
  };
  int count_opt_list = 19;

  association = TRUE;
  s_ptr = NULL;
/*
 *   Read eqn from file and call parser
 */
  opt_save = OPTION_DEFAULT;
  return_value = UNKNOWN;
  for (;;)
  {
    opt = get_option (opt_list, count_opt_list, &next_char);
    if (opt == OPTION_DEFAULT)
    {
      opt = opt_save;
    }
    switch (opt)
    {
    case OPTION_EOF:		/* end of file */
      return_value = EOF;
      break;
    case OPTION_KEYWORD:	/* keyword */
      return_value = KEYWORD;
      break;
    case OPTION_ERROR:
      input_error++;
      error_msg ("Unknown input in EXCHANGE_SPECIES keyword.", CONTINUE);
      error_msg (line_save, CONTINUE);
      break;
    case 0:			/* no_check */
      if (s_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      s_ptr->check_equation = FALSE;
      break;
    case 1:			/* check */
      if (s_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      s_ptr->check_equation = TRUE;
      break;
    case 2:			/* mb */
    case 3:			/* mass_balance */
    case 11:			/* mole_balance */
      if (s_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      count_elts = 0;
      paren_count = 0;
      copy_token (token, &next_char, &i);
      s_ptr->mole_balance = string_hsave (token);
      ptr = token;
      get_secondary_in_species (&ptr, 1.0);
      s_ptr->next_secondary =
	(struct elt_list *) free_check_null (s_ptr->next_secondary);
      s_ptr->next_secondary = elt_list_save ();
/* debug
			for (i = 0; i < count_elts; i++) {
				output_msg(OUTPUT_MESSAGE,"%s\t%f\n", elt_list[i].elt->name,
					elt_list[i].coef);
			}
 */
      opt_save = OPTION_DEFAULT;
      break;
    case 4:			/* log_k */
    case 5:			/* logk */
      if (s_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      read_log_k_only (next_char, &s_ptr->logk[0]);
      opt_save = OPTION_DEFAULT;
      break;
    case 6:			/* delta_h */
    case 7:			/* deltah */
      if (s_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      read_delta_h_only (next_char, &s_ptr->logk[1], &s_ptr->original_units);
      opt_save = OPTION_DEFAULT;
      break;
    case 8:			/* analytical_expression */
    case 9:			/* a_e */
    case 10:			/* ae */
      if (s_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      read_analytical_expression_only (next_char, &(s_ptr->logk[2]));
      opt_save = OPTION_DEFAULT;
      break;
    case 12:			/* gamma data */
      if (s_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      s_ptr->exch_gflag = 2;
      i = sscanf (next_char, SCANFORMAT SCANFORMAT, &s_ptr->dha, &s_ptr->dhb);
      if (i < 2)
      {
	sprintf (error_string, "Expecting 2 activity-"
		 "coefficient parameters, a and b.");
	warning_msg (error_string);
      }
      if (s_ptr->dha == 0 && s_ptr->dhb == 0)
      {
	s_ptr->dhb = 99.9;
	s_ptr->exch_gflag = 1;
      }
      opt_save = OPTION_DEFAULT;
      break;
    case 13:			/* davies eqn */
      if (s_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      s_ptr->exch_gflag = 1;
      s_ptr->dha = 0;
      s_ptr->dhb = 99.9;
      opt_save = OPTION_DEFAULT;
      break;
    case 14:			/* offset */
      if (s_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      if (sscanf (next_char, SCANFORMAT, &offset) != 1)
      {
	error_msg ("No offset for log K given", STOP);
      }
      s_ptr->logk[0] += offset;
      opt_save = OPTION_DEFAULT;
      break;
    case 15:			/* llnl_gamma */
      if (s_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      s_ptr->exch_gflag = 7;	/* llnl D-H */
      i = sscanf (next_char, SCANFORMAT, &s_ptr->dha);
      if (i < 1)
      {
	sprintf (error_string,
		 "Expecting activity-coefficient parameter, a.");
	warning_msg (error_string);
      }
      opt_save = OPTION_DEFAULT;
      break;
    case 16:			/* add_logk */
    case 17:			/* add_log_k */
      if (s_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      if (s_ptr->count_add_logk == 0)
      {
	s_ptr->add_logk =
	  (struct name_coef *) PHRQ_malloc (sizeof (struct name_coef));
	if (s_ptr->add_logk == NULL)
	  malloc_error ();
      }
      else
      {
	s_ptr->add_logk =
	  (struct name_coef *) PHRQ_realloc (s_ptr->add_logk,
					     (size_t) ((s_ptr->
							count_add_logk +
							1) *
						       sizeof (struct
							       name_coef)));
	if (s_ptr->add_logk == NULL)
	  malloc_error ();
      }
      /* read name */
      if (copy_token (token, &next_char, &i) == EMPTY)
      {
	input_error++;
	sprintf (error_string, "Expected the name of a NAMED_EXPRESSION.");
	error_msg (error_string, CONTINUE);
	break;
      }
      s_ptr->add_logk[s_ptr->count_add_logk].name = string_hsave (token);
      /* read coef */
      i =
	sscanf (next_char, SCANFORMAT,
		&s_ptr->add_logk[s_ptr->count_add_logk].coef);
      if (i <= 0)
      {
	s_ptr->add_logk[s_ptr->count_add_logk].coef = 1;
      }
      s_ptr->count_add_logk++;
      opt_save = OPTION_DEFAULT;
      break;
    case 18:			/* add_constant */
      if (s_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      if (s_ptr->count_add_logk == 0)
      {
	s_ptr->add_logk =
	  (struct name_coef *) PHRQ_malloc (sizeof (struct name_coef));
	if (s_ptr->add_logk == NULL)
	  malloc_error ();
      }
      else
      {
	s_ptr->add_logk =
	  (struct name_coef *) PHRQ_realloc (s_ptr->add_logk,
					     (size_t) ((s_ptr->
							count_add_logk +
							1) *
						       sizeof (struct
							       name_coef)));
	if (s_ptr->add_logk == NULL)
	  malloc_error ();
      }
      i =
	sscanf (next_char, SCANFORMAT,
		&s_ptr->add_logk[s_ptr->count_add_logk].coef);
      if (i <= 0)
      {
	input_error++;
	sprintf (error_string,
		 "Expected the constant to add for log_K definition.");
	error_msg (error_string, CONTINUE);
	break;
      }
      /* set name */
      s_ptr->add_logk[s_ptr->count_add_logk].name =
	string_hsave ("XconstantX");
      /* read coef */
      s_ptr->count_add_logk++;
      opt_save = OPTION_DEFAULT;
      break;
    case OPTION_DEFAULT:
/*
 *   Get exchange species information and parse equation
 */
      s_ptr = NULL;
      if (parse_eq (line, &next_elt, association) == ERROR)
      {
	parse_error++;
	error_msg ("Parsing equation.", CONTINUE);
	error_msg (line_save, CONTINUE);
	break;
      }
/*
 *   Get pointer to each species in the reaction, store new species if necessary
 */
      trxn.token[0].s = s_store (trxn.token[0].name, trxn.token[0].z, TRUE);
      for (i = 1; i < count_trxn; i++)
      {
	trxn.token[i].s =
	  s_store (trxn.token[i].name, trxn.token[i].z, FALSE);
      }
/*
 *   Save element list and carbon, hydrogen, and oxygen in species
 */
      trxn.token[0].s->next_elt = next_elt;
      for (; next_elt->elt != NULL; next_elt++)
      {
	if (strcmp (next_elt->elt->name, "C") == 0)
	{
	  trxn.token[0].s->carbon = next_elt->coef;
	}
	if (strcmp (next_elt->elt->name, "H") == 0)
	{
	  trxn.token[0].s->h = next_elt->coef;
	}
	if (strcmp (next_elt->elt->name, "O") == 0)
	{
	  trxn.token[0].s->o = next_elt->coef;
	}
      }
/*
 *   Find valence of cation from coefficients of reaction components
 *   Changed to be coefficient of exchanger
 */
      exchange_coef = 0.0;
      for (i = 1; i < count_trxn; i++)
      {
	if (trxn.token[i].s->type == EX)
	{
	  exchange_coef = trxn.token[i].coef;
	}
      }
      trxn.token[0].s->equiv = exchange_coef;
/*
 *   Malloc space for species reaction
 */
      trxn.token[0].s->rxn = rxn_alloc (count_trxn + 1);
/*
 *   Copy reaction to reaction for species
 */
      token_ptr = trxn.token[0].s->rxn->token;
      for (i = 0; i < count_trxn; i++)
      {
	token_ptr[i].s = trxn.token[i].s;
	token_ptr[i].coef = trxn.token[i].coef;
      }
      token_ptr[i].s = NULL;
/*
 *   Set type for species
 */
      trxn.token[0].s->type = EX;
      s_ptr = trxn.token[0].s;
/*
 *   Set gamma data
 */
      s_ptr->gflag = 4;
      s_ptr->exch_gflag = 3;
      s_ptr->dha = 0.0;
      s_ptr->dhb = 0.0;
      opt_save = OPTION_DEFAULT;
/*
 *  Save as a phase for inverse modeling only
 */
      phase_ptr = phase_store (s_ptr->name);
      if (phase_ptr == NULL)
      {
	input_error++;
	sprintf (error_string, "Copying exchange to phases.");
	error_msg (error_string, CONTINUE);
      }
      phase_ptr->formula = s_ptr->name;
      phase_ptr->check_equation = FALSE;
      phase_ptr->type = EX;
      phase_ptr->next_elt = elt_list_dup (s_ptr->next_elt);
      phase_ptr->rxn = rxn_dup (s_ptr->rxn);
      break;
    }
    if (return_value == EOF || return_value == KEYWORD)
      break;
  }
  return (return_value);
}

/* ---------------------------------------------------------------------- */
int
read_exchange (void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads exchange data
 *
 *      Arguments:
 *         none
 *
 *      Returns:
 *         KEYWORD if keyword encountered, input_error may be incremented if
 *                    a keyword is encountered in an unexpected position
 *         EOF     if eof encountered while reading mass balance concentrations
 *         ERROR   if error occurred reading data
 *
 */
  int i, j, l, n, count_comps;
  int n_user, n_user_end;
  LDBLE conc;
  char *ptr;
  char *description;
  char token[MAX_LENGTH], token1[MAX_LENGTH];
  struct exchange *exchange_ptr;

  int return_value, opt;
  char *next_char;
  const char *opt_list[] = {
    "equilibrate",		/* 0 */
    "equil",			/* 1 */
    "pitzer_exchange_gammas",	/* 2 */
    "exchange_gammas",		/* 3 */
    "gammas"			/* 4 */
  };
  int count_opt_list = 5;
/*
 * kin_exch is for exchangers, related to kinetically reacting minerals
 *    they are defined if "sites" is followed by mineral name:
 *    Z         Manganite                ('equi' or 'kine')      0.25
 *    ^Name     ^equi or kinetic mineral ^switch                  ^prop.factor
 */
/*
 *   Read exchange number and description
 */
  ptr = line;
  read_number_description (ptr, &n_user, &n_user_end, &description);
/*
 *   Find space for exchange data
 */
  exchange_ptr = exchange_search (n_user, &n, FALSE);
  if (exchange_ptr != NULL)
  {
    exchange_free (exchange_ptr);
  }
  else
  {
    n = count_exchange++;
    space ((void **) ((void *) &exchange), count_exchange, &max_exchange,
	   sizeof (struct exchange));
  }
/*
 *   Default values
 */
  exchange_init (&(exchange[n]), n_user, n_user_end, description);
  free_check_null (description);
  if (use.exchange_in == FALSE)
  {
    use.exchange_in = TRUE;
    use.n_exchange_user = n_user;
  }
/*
 *   Read exchange data
 */
  count_comps = 0;
  return_value = UNKNOWN;
  for (;;)
  {
    opt = get_option (opt_list, count_opt_list, &next_char);
    switch (opt)
    {
    case OPTION_EOF:		/* end of file */
      return_value = EOF;
      break;
    case OPTION_KEYWORD:	/* keyword */
      return_value = KEYWORD;
      break;
    case OPTION_ERROR:
      input_error++;
      error_msg ("Unknown input in EXCHANGE keyword.", CONTINUE);
      error_msg (line_save, CONTINUE);
      break;
    case 0:			/* equilibrate */
    case 1:
      /*
       *   Read solution to equilibrate with
       */
      for (;;)
      {
	i = copy_token (token, &next_char, &l);
	if (i == DIGIT)
	{
	  sscanf (token, "%d", &exchange[n].n_solution);
	  exchange[n].new_def = TRUE;
	  exchange[n].solution_equilibria = TRUE;
	  break;
	}
	if (i == EMPTY)
	{
	  error_msg
	    ("Expected a solution number with which to equilibrate exchanger.",
	     CONTINUE);
	  error_msg (line_save, CONTINUE);
	  input_error++;
	  break;
	}
      }
      break;
    case 2:			/* pitzer_exchange_gammas */
    case 3:			/* exchange_gammas */
    case 4:			/* gammas */
      exchange[n].pitzer_exchange_gammas = get_true_false (next_char, TRUE);
      break;
    case OPTION_DEFAULT:
      exchange[n].comps =
	(struct exch_comp *) PHRQ_realloc (exchange[n].comps,
					   (size_t) (count_comps +
						     1) *
					   sizeof (struct exch_comp));
      if (exchange[n].comps == NULL)
	malloc_error ();
      exchange[n].comps[count_comps].formula = NULL;
      exchange[n].comps[count_comps].formula_z = 0.0;
      exchange[n].comps[count_comps].formula_totals = NULL;
      exchange[n].comps[count_comps].moles = 0.0;
      exchange[n].comps[count_comps].la = 0.0;
      exchange[n].comps[count_comps].charge_balance = 0.0;
      exchange[n].comps[count_comps].phase_name = NULL;
      exchange[n].comps[count_comps].phase_proportion = 0.0;
      exchange[n].comps[count_comps].rate_name = NULL;
      ptr = line;
      i = copy_token (token, &ptr, &l);
      /*
       *   Species formula is stored in token
       */
      if (i != UPPER && token[0] != '[')
      {
	error_msg ("Expected exchanger name to begin with a capital letter.",
		   CONTINUE);
	error_msg (line_save, CONTINUE);
	input_error++;
	break;
      }
      exchange[n].comps[count_comps].formula = string_hsave (token);
      i = copy_token (token1, &ptr, &l);
      if (i == DIGIT)
      {
	/*
	 *   Read exchange concentration
	 */

	/* exchanger conc. is read directly .. */
	if (sscanf (token1, SCANFORMAT, &conc) < 1)
	{
	  error_msg ("Expected concentration of exchanger.", CONTINUE);
	  error_msg (line_save, CONTINUE);
	  input_error++;
	  break;
	}
	j = copy_token (token1, &ptr, &l);
	if (j == UPPER || j == LOWER)
	{
	  exchange[n].comps[count_comps].rate_name = string_hsave (token1);
	  if (copy_token (token1, &ptr, &l) != DIGIT)
	  {
	    error_msg
	      ("Expected a coefficient to relate exchange to kinetic reaction.\n",
	       CONTINUE);
	    input_error++;
	    break;
	  }
	  sscanf (token1, SCANFORMAT,
		  &exchange[n].comps[count_comps].phase_proportion);
	  exchange[n].related_rate = TRUE;
	}
	/*
	 *   Read equilibrium phase name or kinetics rate name
	 */
      }
      else if (i != EMPTY)
      {

	/* exchanger conc. is related to mineral or kinetics */
	exchange[n].comps[count_comps].phase_name = string_hsave (token1);
	j = copy_token (token1, &ptr, &l);
	if (j == DIGIT)
	{
	  exchange[n].related_phases = TRUE;
	}
	else
	{
	  if (token1[0] == 'K' || token1[0] == 'k')
	  {
	    exchange[n].comps[count_comps].rate_name =
	      exchange[n].comps[count_comps].phase_name;
	    exchange[n].comps[count_comps].phase_name = NULL;
	    exchange[n].related_rate = TRUE;
	  }
	  else if (token1[0] == 'E' || token1[0] == 'e')
	  {
	    exchange[n].related_phases = TRUE;
	  }
	  else
	  {
	    error_msg
	      ("Character string expected to be 'equilibrium_phase' or 'kinetics' to relate exchange to mineral or kinetic reaction.\n",
	       CONTINUE);
	    input_error++;
	    break;
	  }
	  j = copy_token (token1, &ptr, &l);
	}


	if (j != DIGIT)
	{
	  error_msg
	    ("Expected a coefficient to relate exchanger to mineral or kinetic reaction.\n",
	     CONTINUE);
	  input_error++;
	  break;
	}
	sscanf (token1, SCANFORMAT,
		&exchange[n].comps[count_comps].phase_proportion);
	/* real conc must be defined in tidy_model */
	conc = 1.0;
      }
      else
      {
	error_msg
	  ("Expected concentration of exchanger, mineral name, or kinetic reaction name.",
	   CONTINUE);
	error_msg (line_save, CONTINUE);
	input_error++;
	break;
      }
      /*
       *   Accumulate elements in elt_list
       */
      count_elts = 0;
      paren_count = 0;
      ptr = token;
      get_elts_in_species (&ptr, conc);
      /*
       *   save formula for adjusting number of exchange sites
       */
      ptr = token;
      get_token (&ptr, token1, &exchange[n].comps[count_comps].formula_z, &l);
      exchange[n].comps[count_comps].formula_totals = elt_list_save ();
      /* 
       *   Save elt_list 
       */
      exchange[n].comps[count_comps].moles = conc;
      exchange[n].comps[count_comps].totals = elt_list_save ();
      exchange[n].comps[count_comps].charge_balance = 0.0;
      count_comps++;
    }
    if (return_value == EOF || return_value == KEYWORD)
      break;
  }
  exchange[n].count_comps = count_comps;
  return (return_value);
}

/* ---------------------------------------------------------------------- */
int
read_exchange_master_species (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Reads master species data from data file or input file
 */
  int j, l;
  char *ptr, *ptr1;
  LDBLE z;
  struct element *elts_ptr;
  struct species *s_ptr;
  char token[MAX_LENGTH], token1[MAX_LENGTH];
  for (;;)
  {
    j = check_line ("Exchange species equation", FALSE, TRUE, TRUE, TRUE);
    if (j == EOF || j == KEYWORD)
    {
      break;
    }
/*
 *   Get element name with valence, allocate space, store
 */
    ptr = line;
/*
 *   Get element name and save pointer to character string
 */
    if (copy_token (token, &ptr, &l) != UPPER && token[0] != '[')
    {
      parse_error++;
      error_msg ("Reading element for master species.", CONTINUE);
      error_msg (line_save, CONTINUE);
      continue;
    }
    /*
       if (token[0] == '[') {
       ptr1 = token;
       get_elt(&ptr, element, &l);
       strcpy(token, element);
       }
     */
    replace ("(+", "(", token);
/*
 *   Delete master if it exists
 */
    master_delete (token);
/*
 *   Increase pointer array, if necessary,  and malloc space
 */
    if (count_master >= max_master)
    {
      space ((void **) ((void *) &master), count_master + 1, &max_master,
	     sizeof (struct master *));
    }
    master[count_master] = master_alloc ();
/*
 *   Set type to EX
 */
    master[count_master]->type = EX;
/*
 *   Save element name
 */
    master[count_master]->elt = element_store (token);
/*
 *   Save pointer to species data for master species
 */
    if ((copy_token (token, &ptr, &l) != UPPER) &&
	token[0] != '[' && (strcmp_nocase_arg1 (token, "e-") != 0))
    {
      parse_error++;
      error_msg ("Reading master species name.", CONTINUE);
      error_msg (line_save, CONTINUE);
      continue;
    }
    s_ptr = s_search (token);
    if (s_ptr != NULL)
    {
      master[count_master]->s = s_ptr;
    }
    else
    {
      ptr1 = token;
      get_token (&ptr1, token1, &z, &l);
      master[count_master]->s = s_store (token1, z, FALSE);
    }
/*
 *   MAKE LISTS OF PRIMARY AND SECONDARY MASTER SPECIES
 */
    master[count_master]->primary = TRUE;
    if (strcmp (master[count_master]->elt->name, "E") != 0)
    {
      elts_ptr = element_store (master[count_master]->elt->name);
      elts_ptr->gfw = 0.0;
    }

    count_master++;
    if (count_master >= max_master)
    {
      space ((void **) ((void *) &master), count_master, &max_master,
	     sizeof (struct master *));
    }
  }
  return (j);
}

/* ---------------------------------------------------------------------- */
int
read_gas_phase (void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads gas phase data
 *
 *      Arguments:
 *         none
 *
 *      Returns:
 *         KEYWORD if keyword encountered, input_error may be incremented if
 *                    a keyword is encountered in an unexpected position
 *         EOF     if eof encountered while reading mass balance concentrations
 *         ERROR   if error occurred reading data
 *
 */
  int i, j, n, l;
  int count_comps;
  int n_user, n_user_end;
  char *ptr;
  char *description;
  char token[MAX_LENGTH];
  struct gas_phase *gas_phase_ptr;
  int return_value, opt;
  char *next_char;
  const char *opt_list[] = {
    "pressure",			/* 0 */
    "volume",			/* 1 */
    "temp",			/* 2 */
    "temperature",		/* 3 */
    "fixed_pressure",		/* 4 */
    "fixed_volume",		/* 5 */
    "equilibrium",		/* 6 */
    "equilibrate",		/* 7 */
    "equil"			/* 8 */
  };
  int count_opt_list = 9;
/*
 *   Read gas_phase number
 */
  ptr = line;
  read_number_description (ptr, &n_user, &n_user_end, &description);
/*
 *   Find old gas_phase or alloc space for new gas_phase
 */
  gas_phase_ptr = gas_phase_search (n_user, &n);
  if (gas_phase_ptr != NULL)
  {
    gas_phase_free (gas_phase_ptr);
  }
  else
  {
    n = count_gas_phase;
    count_gas_phase++;
    space ((void **) ((void *) &gas_phase), count_gas_phase, &max_gas_phase,
	   sizeof (struct gas_phase));
  }
/*
 *   Initialize
 */
  gas_phase_init (&(gas_phase[n]), n_user, n_user_end, description);
  free_check_null (description);
/*
 *   Set use data to first read
 */
  if (use.gas_phase_in == FALSE)
  {
    use.gas_phase_in = TRUE;
    use.n_gas_phase_user = n_user;
  }
/*
 *   Read phases
 */
  count_comps = 0;
  return_value = UNKNOWN;
  for (;;)
  {
    opt = get_option (opt_list, count_opt_list, &next_char);
    switch (opt)
    {
    case OPTION_EOF:		/* end of file */
      return_value = EOF;
      break;
    case OPTION_KEYWORD:	/* keyword */
      return_value = KEYWORD;
      break;
    case OPTION_ERROR:
      input_error++;
      error_msg ("Unknown input in GAS_PHASE keyword.", CONTINUE);
      error_msg (line_save, CONTINUE);
      break;
    case 0:			/* pressure */
      sscanf (next_char, SCANFORMAT, &dummy);
      gas_phase[n].total_p = (LDBLE) dummy;
      break;
    case 1:			/* Volume */
      sscanf (next_char, SCANFORMAT, &(gas_phase[n].volume));
      break;
    case 2:			/* Temperature */
    case 3:
      j = sscanf (next_char, SCANFORMAT, &(gas_phase[n].temperature));
      if (j == 1)
      {
	gas_phase[n].temperature += 273.15;
      }
      break;
    case 4:			/* fixed_pressure */
      gas_phase[n].type = PRESSURE;
      break;
    case 5:			/* fixed_volume */
      gas_phase[n].type = VOLUME;
      break;
    case 6:			/* equilibrate */
    case 7:			/* equilibrium */
    case 8:			/* equil */
/*
 *   Read solution to equilibrate with
 */
      for (;;)
      {
	i = copy_token (token, &next_char, &l);
	if (i == DIGIT)
	{
	  sscanf (token, "%d", &gas_phase[n].n_solution);
	  gas_phase[n].new_def = TRUE;
	  gas_phase[n].solution_equilibria = TRUE;
	  break;
	}
	if (i == EMPTY)
	{
	  error_msg
	    ("Expected a solution number with which to equilibrate gas phase.",
	     CONTINUE);
	  error_msg (line_save, CONTINUE);
	  input_error++;
	  break;
	}
      }
      break;
    case OPTION_DEFAULT:
/*
 *   Make space, set default
 */
      gas_phase[n].comps =
	(struct gas_comp *) PHRQ_realloc (gas_phase[n].comps,
					  (size_t) (count_comps +
						    1) *
					  sizeof (struct gas_comp));
      if (gas_phase[n].comps == NULL)
	malloc_error ();
      gas_phase[n].comps[count_comps].p_read = 0.0;
      gas_phase[n].comps[count_comps].moles = 0.0;
      count_comps++;
/*
 *   Read name
 */
      ptr = line;
      copy_token (token, &ptr, &l);
      gas_phase[n].comps[count_comps - 1].name = string_hsave (token);
      if ((j = copy_token (token, &ptr, &l)) == EMPTY)
      {
	gas_phase[n].comps[count_comps - 1].p_read = NAN;
	break;
      }
/*
 *   Read initial partial pressure of gas
 */

      j =
	sscanf (token, SCANFORMAT,
		&(gas_phase[n].comps[count_comps - 1].p_read));
      if (j != 1)
      {
	error_msg ("Expected partial pressure of gas in gas phase.",
		   CONTINUE);
	error_msg (line_save, CONTINUE);
	input_error++;
      }
      break;
    }
    if (return_value == EOF || return_value == KEYWORD)
      break;
  }
/*
 *   Sort components by name (lowercase)
 */
  gas_phase[n].count_comps = count_comps;
  qsort (gas_phase[n].comps,
	 (size_t) count_comps,
	 (size_t) sizeof (struct gas_comp), gas_comp_compare);
  return (return_value);
}

/* ---------------------------------------------------------------------- */
int
read_inverse (void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads data for mass_balance calculations
 *
 *      Arguments:
 *         none
 *
 *      Returns:
 *         KEYWORD if keyword encountered, input_error may be incremented if
 *                    a keyword is encountered in an unexpected position
 *         EOF     if eof encountered while reading mass balance concentrations
 *         ERROR   if error occurred reading data
 *
 */
  int n, j;
  int n_user, n_user_end;
  char *ptr;
  char *description;
  LDBLE range_max, inv_tol, water_uncertainty;

  int return_value, opt, opt_save;
  char *next_char;
  const char *opt_list[] = {
    "solutions",		/* 0 */
    "uncertainty",		/* 1 */
    "uncertainties",		/* 2 */
    "balances",			/* 3 */
    "phase_data",		/* 4 */
    "range",			/* 5 */
    "minimal",			/* 6 */
    "minimum",			/* 7 */
    "balance",			/* 8 */
    "bal",			/* 9 */
    "sol",			/* 10 */
    "phases",			/* 11 */
    "ranges",			/* 12 */
    "tolerance",		/* 13 */
    "u_water",			/* 14 */
    "uncertainty_water",	/* 15 */
    "force",			/* 16 */
    "force_solution",		/* 17 */
    "force_solutions",		/* 18 */
    "isotopes",			/* 19 */
    "mineral_water",		/* 20 */
    "phase",			/* 21 */
    "multiple_precision",	/* 22 */
    "mp_tolerance",		/* 23 */
    "censor_mp",		/* 24 */
    "lon_netpath",	        /* 25 */
    "pat_netpath"	        /* 26 */
  };
  int count_opt_list = 27;

  ptr = line;
/*
 *   Read solution number and description
 */
  read_number_description (ptr, &n_user, &n_user_end, &description);
/*
 *   Malloc space for solution data
 */
  if (inverse_search (n_user, &n) != NULL)
  {
    inverse_delete (n);
  }
  inverse_alloc ();
  n = count_inverse - 1;
/*
 *   Initialize structure and use
 */
  inverse[n].new_def = TRUE;
  inverse[n].n_user = n_user;
  inverse[n].range = FALSE;
  inverse[n].range_max = 1000.;
  inverse[n].tolerance = 1e-10;
  inverse[n].minimal = FALSE;
  inverse[n].description = description;
  inverse[n].count_uncertainties = 1;
  inverse[n].uncertainties[0] = 0.05;
  inverse[n].count_ph_uncertainties = 1;
  inverse[n].ph_uncertainties[0] = 0.05;
  inverse[n].water_uncertainty = 0.0;
  inverse[n].mineral_water = TRUE;
  inverse[n].mp = FALSE;
  inverse[n].mp_tolerance = 1e-12;
  inverse[n].mp_censor = 1e-20;
  inverse[n].netpath = NULL;
  inverse[n].pat = NULL;
/*
 *   Read data for inverse modeling
 */
  opt_save = OPTION_ERROR;
  return_value = UNKNOWN;
  for (;;)
  {
    opt = get_option (opt_list, count_opt_list, &next_char);
    if (opt == OPTION_DEFAULT)
    {
      opt = opt_save;
    }
    opt_save = opt;
    switch (opt)
    {
    case OPTION_EOF:		/* end of file */
      return_value = EOF;
      break;
    case OPTION_KEYWORD:	/* keyword */
      return_value = KEYWORD;
      break;
    case OPTION_DEFAULT:
    case OPTION_ERROR:
      input_error++;
      error_msg ("Unknown input in INVERSE_MODELING keyword.", CONTINUE);
      error_msg (line_save, CONTINUE);
      break;
    case 0:			/* solutions */
    case 10:			/* solution */
      inverse[n].solns =
	read_list_ints (&next_char, &inverse[n].count_solns, TRUE);
      opt_save = OPTION_ERROR;
      break;
    case 1:			/* uncertainty */
    case 2:			/* uncertainties */
      inverse[n].uncertainties =
	(double *) free_check_null (inverse[n].uncertainties);
      inverse[n].uncertainties =
	read_list_doubles (&next_char, &inverse[n].count_uncertainties);
      opt_save = OPTION_ERROR;
      break;
    case 3:			/* balances */
    case 8:			/* balance */
    case 9:			/* bal */
      read_inv_balances (&(inverse[n]), next_char);
      break;
    case 4:			/* phase_data */
    case 11:			/* phases */
    case 21:			/* phase */
      read_inv_phases (&(inverse[n]), next_char);
      break;
    case 5:			/* range */
    case 12:			/* ranges */
      inverse[n].range = TRUE;
      j = sscanf (next_char, SCANFORMAT, &range_max);
      if (j == 1)
      {
	inverse[n].range_max = range_max;
      }
      opt_save = OPTION_ERROR;
      break;
    case 6:			/* minimal */
    case 7:			/* minimum */
      inverse[n].minimal = TRUE;
      opt_save = OPTION_ERROR;
      break;
    case 13:			/* tolerance */
      j = sscanf (next_char, SCANFORMAT, &inv_tol);
      if (j == 1)
      {
	inverse[n].tolerance = inv_tol;
      }
      opt_save = OPTION_ERROR;
      break;
    case 14:			/* u_water */
    case 15:			/* uncertainty_water */
      j = sscanf (next_char, SCANFORMAT, &water_uncertainty);
      if (j == 1)
      {
	inverse[n].water_uncertainty = water_uncertainty;
      }
      opt_save = OPTION_ERROR;
      break;
    case 16:			/* force */
    case 17:			/* force_solution */
    case 18:			/* force_solutions */
      inverse[n].force_solns =
	(int *) free_check_null (inverse[n].force_solns);
      inverse[n].force_solns =
	read_list_t_f (&next_char, &inverse[n].count_force_solns);
      opt_save = OPTION_ERROR;
      break;
    case 19:			/* isotope values */
      read_inv_isotopes (&(inverse[n]), next_char);
      break;
    case 20:			/* mineral_water */
      inverse[n].mineral_water = get_true_false (next_char, TRUE);
      break;
    case 22:			/* multiple_precision */
      inverse[n].mp = get_true_false (next_char, TRUE);
      break;
    case 23:			/* mp_tolerance */
      j = sscanf (next_char, SCANFORMAT, &inv_tol);
      if (j == 1)
      {
	inverse[n].mp_tolerance = fabs (inv_tol);
      }
      break;
    case 24:			/* censor_mp */
      j = sscanf (next_char, SCANFORMAT, &inv_tol);
      if (j == 1)
      {
	inverse[n].mp_censor = fabs (inv_tol);
      }
      break;
    case 25:			/* lon_netpath */
      /*copy_token(file_name, &next_char, &l); */
      if (string_trim (next_char) != EMPTY) 
      {
	inverse[n].netpath = string_hsave(next_char);
      } else {
	inverse[n].netpath = string_hsave("netpath");
      }
      opt_save = OPTION_ERROR;
      break;
    case 26:			/* pat_netpath */
      /*copy_token(file_name, &next_char, &l); */
      if (string_trim (next_char) != EMPTY) 
      {
	inverse[n].pat = string_hsave(next_char);
      } else {
	inverse[n].pat = string_hsave("netpath");
      }
      opt_save = OPTION_ERROR;
      break;
    }
    if (return_value == EOF || return_value == KEYWORD)
      break;
  }
/*
 *  Default: soln 1 -> soln 2
 */
  if (inverse[n].count_solns == 0)
  {
    inverse[n].solns = (int *) PHRQ_malloc (2 * sizeof (int));
    if (inverse[n].solns == NULL)
      malloc_error ();
    inverse[n].solns[0] = 1;
    inverse[n].solns[1] = 2;
    inverse[n].count_solns = 2;
  }
/*
 *   Sort isotopes
 */
  if (inverse[n].count_isotopes > 0)
  {
    qsort (inverse[n].isotopes,
	   (size_t) inverse[n].count_isotopes,
	   (size_t) sizeof (struct inv_isotope), inverse_isotope_compare);
  }

  if (inverse[n].count_i_u > 0)
  {
    qsort (inverse[n].i_u,
	   (size_t) inverse[n].count_i_u,
	   (size_t) sizeof (struct inv_isotope), inverse_isotope_compare);
  }

  return (return_value);
}

/* ---------------------------------------------------------------------- */
int
read_inv_balances (struct inverse *inverse_ptr, char *ptr)
/* ---------------------------------------------------------------------- */
{
  int j, l, count;
  char token[MAX_LENGTH];
/*
 *   Read element name
 */
  j = copy_token (token, &ptr, &l);
  if (j == EMPTY)
  {
    return (OK);
  }
  else if (j == LOWER && strcmp_nocase_arg1 (token, "ph") != 0)
  {
    error_msg ("Expecting element name.", CONTINUE);
    error_msg (line_save, CONTINUE);
    input_error++;
  }
  else if (strcmp_nocase_arg1 (token, "ph") != 0)
  {
    inverse_ptr->elts =
      (struct inv_elts *) PHRQ_realloc (inverse_ptr->elts,
					(size_t) (inverse_ptr->count_elts +
						  1) *
					sizeof (struct inv_elts));
    if (inverse_ptr->elts == NULL)
      malloc_error ();
    replace ("(+", "(", token);
    inverse_ptr->elts[inverse_ptr->count_elts].name = string_hsave (token);
/*
 *   Read element uncertainties
 */
    inverse_ptr->elts[inverse_ptr->count_elts].uncertainties =
      read_list_doubles (&ptr, &count);
    inverse_ptr->elts[inverse_ptr->count_elts].count_uncertainties = count;
    inverse_ptr->count_elts++;
  }
  else if (strcmp_nocase_arg1 (token, "ph") == 0)
  {
    inverse_ptr->ph_uncertainties =
      (double *) free_check_null (inverse_ptr->ph_uncertainties);
    inverse_ptr->ph_uncertainties = read_list_doubles (&ptr, &count);
    inverse_ptr->count_ph_uncertainties = count;
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
read_inv_isotopes (struct inverse *inverse_ptr, char *ptr)
/* ---------------------------------------------------------------------- */
{
  int i, j, l, l1, l2, count;
  LDBLE isotope_number;
  char token[MAX_LENGTH], token1[MAX_LENGTH];
  char *ptr1, *ptr2, *redox_name, *element_name;
/*
 *   Read element name
 */
  ptr1 = ptr;
  j = copy_token (token, &ptr1, &l);
/*
 *   ptr1 is start of uncertainties
 */
  if (j == EMPTY)
  {
    return (OK);
  }
  else if (j != DIGIT)
  {
    error_msg ("Expecting isotope to begin with isotope number.", CONTINUE);
    error_msg (line_save, CONTINUE);
    input_error++;
    return (ERROR);
  }
/*
 *   Read isotope name
 */
  ptr2 = token;
  get_num (&ptr2, &isotope_number);
  if (ptr2[0] == '\0' || isupper ((int) ptr2[0]) == FALSE)
  {
    error_msg ("Expecting element name.", CONTINUE);
    error_msg (line_save, CONTINUE);
    input_error++;
    return (ERROR);
  }

  /* redox state name with parentheses */
  redox_name = string_hsave (ptr2);

  copy_token (token, &ptr2, &l1);
  replace ("(", " ", token);
  ptr2 = token;

  /* element name, without parentheses */
  copy_token (token1, &ptr2, &l2);
  element_name = string_hsave (token1);

/*
 *  add element name to inv_ptr->isotopes
 */
  for (i = 0; i < inverse_ptr->count_isotopes; i++)
  {
    if (element_name == inverse_ptr->isotopes[i].elt_name)
      break;
  }
  if (i == inverse_ptr->count_isotopes)
  {
    inverse_ptr->isotopes =
      (struct inv_isotope *) PHRQ_realloc (inverse_ptr->isotopes,
					   (size_t) (inverse_ptr->
						     count_isotopes +
						     1) *
					   sizeof (struct inv_isotope));
    if (inverse_ptr->isotopes == NULL)
      malloc_error ();
    inverse_ptr->isotopes[inverse_ptr->count_isotopes].isotope_number =
      isotope_number;
    inverse_ptr->isotopes[inverse_ptr->count_isotopes].elt_name =
      element_name;
    inverse_ptr->isotopes[inverse_ptr->count_isotopes].uncertainties =
      (double *) PHRQ_malloc ((size_t) sizeof (LDBLE));
    if (inverse_ptr->isotopes[inverse_ptr->count_isotopes].uncertainties ==
	NULL)
      malloc_error ();
    inverse_ptr->count_isotopes++;
  }
/*
 *  add redox state name to inv_ptr->i_u
 */
  inverse_ptr->i_u =
    (struct inv_isotope *) PHRQ_realloc (inverse_ptr->i_u,
					 (size_t) (inverse_ptr->count_i_u +
						   1) *
					 sizeof (struct inv_isotope));
  if (inverse_ptr->i_u == NULL)
    malloc_error ();
  inverse_ptr->i_u[inverse_ptr->count_i_u].elt_name = redox_name;
  inverse_ptr->i_u[inverse_ptr->count_i_u].isotope_number = isotope_number;
/*
 *   Read isotope uncertainties
 */
  inverse_ptr->i_u[inverse_ptr->count_i_u].uncertainties =
    read_list_doubles (&ptr1, &count);
  inverse_ptr->i_u[inverse_ptr->count_i_u].count_uncertainties = count;
  inverse_ptr->count_i_u++;
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
read_inv_phases (struct inverse *inverse_ptr, char *ptr)
/* ---------------------------------------------------------------------- */
{
  int j, l;
  int count_isotopes;
  char token[MAX_LENGTH], token1[MAX_LENGTH];
  char *ptr1;
  struct isotope *isotopes;
/*
 *   Read phase name
 */
  j = copy_token (token, &ptr, &l);
  if (j == EMPTY)
    return (OK);
  inverse_ptr->phases =
    (struct inv_phases *) PHRQ_realloc (inverse_ptr->phases,
					(size_t) (inverse_ptr->count_phases +
						  1) *
					sizeof (struct inv_phases));
  if (inverse_ptr->phases == NULL)
    malloc_error ();
  inverse_ptr->phases[inverse_ptr->count_phases].name = string_hsave (token);
/*
 *   Read constraint, force, and isotopes
 */
  inverse_ptr->phases[inverse_ptr->count_phases].constraint = EITHER;
  inverse_ptr->phases[inverse_ptr->count_phases].force = FALSE;
  count_isotopes = 0;
  isotopes = (struct isotope *) PHRQ_malloc (sizeof (struct isotope));
  if (isotopes == NULL)
    malloc_error ();

  for (;;)
  {
    j = copy_token (token, &ptr, &l);
    if (j == EMPTY)
      break;
    strcpy (token1, token);
    str_tolower (token1);
    if (token1[0] == 'p')
    {
      inverse_ptr->phases[inverse_ptr->count_phases].constraint = PRECIPITATE;
    }
    else if (token1[0] == 'd')
    {
      inverse_ptr->phases[inverse_ptr->count_phases].constraint = DISSOLVE;
    }
    else if (token[0] == 'f')
    {
      inverse_ptr->phases[inverse_ptr->count_phases].force = TRUE;
    }
    else if (j == DIGIT)
    {
/* 
 *   read isotope data
 */
      isotopes =
	(struct isotope *) PHRQ_realloc (isotopes,
					 (size_t) (count_isotopes +
						   1) *
					 sizeof (struct isotope));
      if (isotopes == NULL)
	malloc_error ();
      ptr1 = token;

      /* isotope number */
      get_num (&ptr1, &(isotopes[count_isotopes].isotope_number));
      if (ptr1[0] == '\0' || isupper ((int) ptr1[0]) == FALSE)
      {
	sprintf (error_string, "Expecting element name: %s.", ptr1);
	error_msg (error_string, CONTINUE);
	error_msg (line_save, CONTINUE);
	input_error++;
	break;
      }

      /* element name */
      isotopes[count_isotopes].elt_name = string_hsave (ptr1);

      /* ratio */
      j = copy_token (token, &ptr, &l);
      if (j != DIGIT)
      {
	error_msg ("Expecting isotope ratio for phase.", CONTINUE);
	error_msg (line_save, CONTINUE);
	input_error++;
	break;
      }
      sscanf (token, SCANFORMAT, &(isotopes[count_isotopes].ratio));

      /* read and store isotope ratio uncertainty */
      if (copy_token (token, &ptr, &l) != DIGIT)
      {
	input_error++;
	sprintf (error_string,
		 "Expected numeric value for uncertainty in isotope ratio.");
	error_msg (error_string, CONTINUE);
	continue;
      }
      sscanf (token, SCANFORMAT,
	      &(isotopes[count_isotopes].ratio_uncertainty));

      count_isotopes++;
    }
    else
    {
      sprintf (error_string, "Unknown option for inverse modeling phase.");
      warning_msg (error_string);
    }
  }
  if (count_isotopes > 0)
  {
    inverse_ptr->phases[inverse_ptr->count_phases].isotopes = isotopes;
    inverse_ptr->phases[inverse_ptr->count_phases].count_isotopes =
      count_isotopes;
  }
  else
  {
    inverse_ptr->phases[inverse_ptr->count_phases].isotopes = NULL;
    inverse_ptr->phases[inverse_ptr->count_phases].count_isotopes = 0;
    isotopes = (struct isotope *) free_check_null (isotopes);
  }
  inverse_ptr->count_phases++;
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
read_kinetics (void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads kinetics data
 *
 *      Arguments:
 *         none
 *
 *      Returns:
 *         KEYWORD if keyword encountered, input_error may be incremented if
 *                    a keyword is encountered in an unexpected position
 *         EOF     if eof encountered while reading mass balance concentrations
 *         ERROR   if error occurred reading data
 *
 */
/*
 *   Read kinetics
 */
  int i, j, k, l, count_comps, count_steps, count_list;
  char *ptr;
  char *description;
  char token[MAX_LENGTH];
  int n;
  int n_user, n_user_end;
  struct kinetics *kinetics_ptr;
  struct kinetics_comp *kinetics_comp_ptr;
  LDBLE step, coef;

  int return_value, opt;
  char *next_char;
  const char *opt_list[] = {
    "tol",			/* 0 */
    "m",			/* 1 */
    "m0",			/* 2 */
    "parms",			/* 3 */
    "formula",			/* 4 */
    "steps",			/* 5 */
    "step_divide",		/* 6 */
    "parameters",		/* 7 */
    "runge-kutta",		/* 8 */
    "runge_kutta",		/* 9 */
    "rk",			/* 10 */
    "bad_step_max",		/* 11 */
    "cvode",			/* 12 */
    "cvode_steps",		/* 13 */
    "cvode_order"		/* 14 */
  };
  int count_opt_list = 15;

/*
 *   Read kinetics number
 */
  ptr = line;
  read_number_description (ptr, &n_user, &n_user_end, &description);

/*
 *   Find space for kinetics data
 */
  kinetics_ptr = kinetics_search (n_user, &n, FALSE);
  if (kinetics_ptr != NULL)
  {
    kinetics_free (kinetics_ptr);
  }
  else
  {
    space ((void **) ((void *) &kinetics), count_kinetics, &max_kinetics,
	   sizeof (struct kinetics));
    n = count_kinetics++;
  }
/*
 *   Set use data to first read
 */
  if (use.kinetics_in == FALSE)
  {
    use.kinetics_in = TRUE;
    use.n_kinetics_user = n_user;
  }
/*
 *   Initialize
 */
  kinetics_init (&(kinetics[n]), n_user, n_user_end, description);
  free_check_null (description);
  kinetics_ptr = &kinetics[n];

  count_steps = 0;
/*
 *   Read kinetics data
 */
  return_value = UNKNOWN;
  kinetics_comp_ptr = NULL;
  for (;;)
  {
    opt = get_option (opt_list, count_opt_list, &next_char);
    switch (opt)
    {
    case OPTION_EOF:		/* end of file */
      return_value = EOF;
      break;
    case OPTION_KEYWORD:	/* keyword */
      return_value = KEYWORD;
      break;
    case OPTION_DEFAULT:	/* allocate space, read new name */
      count_comps = kinetics_ptr->count_comps++;
      kinetics_ptr->comps =
	(struct kinetics_comp *) PHRQ_realloc (kinetics_ptr->comps,
					       (size_t) (count_comps +
							 1) *
					       sizeof (struct kinetics_comp));
      if (kinetics_ptr->comps == NULL)
	malloc_error ();
      kinetics_ptr->comps[count_comps].moles = 0;
      ptr = line;
      copy_token (token, &ptr, &l);
      kinetics_ptr->comps[count_comps].rate_name = string_hsave (token);
#ifdef SKIP
      kinetics_ptr->comps[count_comps].formula =
	kinetics_ptr->comps[count_comps].rate_name;
#endif
      kinetics_ptr->comps[count_comps].tol = 1e-8;
      kinetics_ptr->comps[count_comps].m0 = -1.0;
      kinetics_ptr->comps[count_comps].m = -1.0;
      kinetics_ptr->comps[count_comps].count_c_params = 0;

      kinetics_comp_ptr = &kinetics_ptr->comps[count_comps];
      kinetics_comp_ptr->d_params = (double *) PHRQ_malloc (sizeof (LDBLE));
      if (kinetics_comp_ptr->d_params == NULL)
	malloc_error ();
      kinetics_comp_ptr->count_d_params = 0;

      kinetics_comp_ptr->c_params = (char **) PHRQ_malloc (sizeof (char *));
      if (kinetics_comp_ptr->c_params == NULL)
	malloc_error ();
      kinetics_comp_ptr->count_c_params = 0;

      kinetics_comp_ptr->count_list = 0;
      kinetics_comp_ptr->list =
	(struct name_coef *) PHRQ_malloc (sizeof (struct name_coef));
      if (kinetics_comp_ptr->list == NULL)
	malloc_error ();
      break;
    case OPTION_ERROR:
      input_error++;
      error_msg ("Unknown input in KINETICS keyword.", CONTINUE);
      error_msg (line_save, CONTINUE);
      break;
    case 0:			/* tolerance */
      if (kinetics_comp_ptr == NULL)
      {
	sprintf (error_string, "No rate name has been defined.");
	error_msg (error_string, CONTINUE);
	input_error++;
      }
      else
      {
	if (copy_token (token, &next_char, &l) == DIGIT)
	{
	  kinetics_comp_ptr->tol = strtod (token, &ptr);
	}
	else
	{
	  sprintf (error_string, "Expecting numerical value for tolerance.");
	  error_msg (error_string, CONTINUE);
	  input_error++;
	}
      }
      break;
    case 1:			/* m */
      if (kinetics_comp_ptr == NULL)
      {
	sprintf (error_string, "No rate name has been defined.");
	error_msg (error_string, CONTINUE);
	input_error++;
      }
      else
      {
	if (copy_token (token, &next_char, &l) == DIGIT)
	{
	  kinetics_comp_ptr->m = strtod (token, &ptr);
	}
	else
	{
	  sprintf (error_string,
		   "Expecting numerical value for moles of reactant.");
	  error_msg (error_string, CONTINUE);
	  input_error++;
	}
      }
      break;
    case 2:			/* m0 */
      if (kinetics_comp_ptr == NULL)
      {
	sprintf (error_string, "No rate name has been defined.");
	error_msg (error_string, CONTINUE);
	input_error++;
      }
      else
      {
	if (copy_token (token, &next_char, &l) == DIGIT)
	{
	  kinetics_comp_ptr->m0 = strtod (token, &ptr);
	}
	else
	{
	  sprintf (error_string,
		   "Expecting numerical value for initial moles of reactant.");
	  error_msg (error_string, CONTINUE);
	  input_error++;
	}
      }
      break;
    case 3:			/* parms */
    case 7:			/* parameters */
      if (kinetics_comp_ptr == NULL)
      {
	sprintf (error_string, "No rate name has been defined.");
	error_msg (error_string, CONTINUE);
	input_error++;
      }
      else
      {
	while ((j = copy_token (token, &next_char, &l)) != EMPTY)
	{
	  /*
	   *   Store a LDBLE parameter
	   */
	  if (j == DIGIT)
	  {
	    kinetics_comp_ptr->d_params =
	      (double *) PHRQ_realloc (kinetics_comp_ptr->d_params,
				       (size_t) (kinetics_comp_ptr->
						 count_d_params +
						 1) * sizeof (LDBLE));
	    if (kinetics_comp_ptr->d_params == NULL)
	      malloc_error ();
	    kinetics_comp_ptr->d_params[kinetics_comp_ptr->count_d_params] =
	      strtod (token, &ptr);
	    kinetics_comp_ptr->count_d_params++;
	  }
	  else
	  {
	    /*
	     *   Store a character parameter
	     */
	    kinetics_comp_ptr->c_params =
	      (char **) PHRQ_realloc (kinetics_comp_ptr->c_params,
				      (size_t) (kinetics_comp_ptr->
						count_c_params +
						1) * sizeof (char *));
	    if (kinetics_comp_ptr->c_params == NULL)
	      malloc_error ();
	    kinetics_comp_ptr->c_params[kinetics_comp_ptr->count_c_params] =
	      string_hsave (token);
	    kinetics_comp_ptr->count_c_params++;
	  }
	}
      }
      break;
    case 4:			/* formula */
      if (kinetics_comp_ptr == NULL)
      {
	sprintf (error_string, "No rate name has been defined.");
	error_msg (error_string, CONTINUE);
	input_error++;
      }
      else
      {
	/*
	 *   Store reactant name, default coefficient
	 */
	ptr = next_char;
	while (copy_token (token, &ptr, &l) != EMPTY)
	{
	  if (isalpha ((int) token[0]) || (token[0] == '(')
	      || (token[0] == '['))
	  {
	    count_list = kinetics_comp_ptr->count_list++;
	    kinetics_comp_ptr->list =
	      (struct name_coef *) PHRQ_realloc (kinetics_comp_ptr->list,
						 (size_t) (count_list +
							   1) *
						 sizeof (struct name_coef));
	    if (kinetics_comp_ptr->list == NULL)
	      malloc_error ();
	    kinetics_comp_ptr->list[count_list].name = string_hsave (token);
	    kinetics_comp_ptr->list[count_list].coef = 1.0;
	  }
	  else
	  {
	    /*
	     *   Store relative coefficient
	     */
	    j = sscanf (token, SCANFORMAT, &coef);
	    if (j == 1)
	    {
	      count_list = kinetics_comp_ptr->count_list - 1;
	      kinetics_comp_ptr->list[count_list].coef = coef;
	    }
	    else
	    {
	      error_msg ("Reading relative coefficient of reactant.",
			 CONTINUE);
	      error_msg (line_save, CONTINUE);
	      input_error++;
	    }
	  }
	}
      }
      break;
    case 5:			/* steps */
      /*
       *   Read one or more kinetics time increments
       */
      while ((j = copy_token (token, &next_char, &l)) == DIGIT)
      {
	/*  Read next step increment(s) */
/* multiple, equal timesteps 15 aug. 2005 */
	if (replace ("*", " ", token) == TRUE)
	{
	  if (sscanf (token, "%d" SCANFORMAT, &k, &step) == 2)
	  {
	    for (i = 0; i < k; i++)
	    {
	      count_steps++;
	      kinetics_ptr->steps =
		(double *) PHRQ_realloc (kinetics_ptr->steps,
					 (size_t) count_steps *
					 sizeof (LDBLE));
	      if (kinetics_ptr->steps == NULL)
		malloc_error ();
	      kinetics_ptr->steps[kinetics_ptr->count_steps] = step;
	      kinetics_ptr->count_steps = count_steps;
	    }
	  }
	  else
	  {
	    input_error++;
	    error_msg
	      ("Format error in multiple, equal KINETICS timesteps.\nCorrect is (for example): 20 4*10 2*5 3\n",
	       CONTINUE);
	  }
	}
	else
	{
	  step = strtod (token, &ptr);
	  count_steps++;
	  kinetics_ptr->steps =
	    (double *) PHRQ_realloc (kinetics_ptr->steps,
				     (size_t) count_steps * sizeof (LDBLE));
	  if (kinetics_ptr->steps == NULL)
	    malloc_error ();
	  kinetics_ptr->steps[kinetics_ptr->count_steps] = step;
	  kinetics_ptr->count_steps = count_steps;
	}
      }
      if (j == EMPTY)
	break;
      /*
       *   Read number of increments
       */
      if (kinetics_ptr->count_steps != 1)
      {
	error_msg
	  ("To define equal time increments, only one total time should be defined.",
	   CONTINUE);
	input_error++;
	break;
      }
      do
      {
	j = sscanf (token, "%d", &i);
	if (j == 1 && i > 0)
	{
	  kinetics_ptr->count_steps = -i;
	  break;
	}
	else if (j == 1 && i <= 0)
	{
	  error_msg ("Expecting positive number for number of equal "
		     "time increments for kinetics.", CONTINUE);
	  error_msg (line_save, CONTINUE);
	  input_error++;
	  break;
	}
      }
      while (copy_token (token, &next_char, &l) != EMPTY);
      break;
    case 6:			/* step_divide */
      if (copy_token (token, &next_char, &l) == DIGIT)
      {
	sscanf (token, SCANFORMAT, &kinetics_ptr->step_divide);
      }
      else
      {
	sprintf (error_string, "Expecting numerical value for step_divide.");
	error_msg (error_string, CONTINUE);
	input_error++;
      }
      break;
    case 8:			/* runge-kutta */
    case 9:			/* runge_kutta */
    case 10:			/* rk */
      j = copy_token (token, &next_char, &l);
      if (j == DIGIT)
      {
	kinetics_ptr->rk = (int) strtod (token, &ptr);
      }
      else if (j == EMPTY)
      {
      }
      else
      {
	sprintf (error_string, "Expecting order for Runge-Kutta method.");
	error_msg (error_string, CONTINUE);
	input_error++;
      }
      break;
    case 11:			/* bad_step_max */
      j = copy_token (token, &next_char, &l);
      if (j == DIGIT)
      {
	kinetics_ptr->bad_step_max = (int) strtod (token, &ptr);
      }
      else if (j == EMPTY)
      {
      }
      else
      {
	sprintf (error_string, "Expecting maximal bad steps number.");
	error_msg (error_string, CONTINUE);
	input_error++;
      }
      break;
    case 12:			/* cvode */
      kinetics[n].use_cvode = get_true_false (next_char, TRUE);
      break;
    case 13:			/* cvode_steps */
      j = copy_token (token, &next_char, &l);
      if (j == DIGIT)
      {
	kinetics_ptr->cvode_steps = (int) strtod (token, &ptr);
      }
      else if (j == EMPTY)
      {
      }
      else
      {
	sprintf (error_string, "Expecting maximum number of steps for one call to cvode.");
	error_msg (error_string, CONTINUE);
	input_error++;
      }
      break;
    case 14:			/* cvode_order */
      j = copy_token (token, &next_char, &l);
      if (j == DIGIT)
      {
	kinetics_ptr->cvode_order = (int) strtod (token, &ptr);
      }
      else if (j == EMPTY)
      {
      }
      else
      {
	sprintf (error_string, "Expecting number of terms (order) used in cvode (1-5).");
	error_msg (error_string, CONTINUE);
	input_error++;
      }
      break;
    }
    if (return_value == EOF || return_value == KEYWORD)
      break;
  }

/*
 *   Default reactant
 */
  for (i = 0; i < kinetics[n].count_comps; i++)
  {
    if (kinetics[n].comps[i].count_list == 0)
    {
      kinetics[n].comps[i].list[0].name = kinetics_ptr->comps[i].rate_name;
      kinetics[n].comps[i].list[0].coef = 1;
      kinetics[n].comps[i].count_list = 1;
    }
  }
/*
 *   Default 1 sec
 */
  if (kinetics[n].count_steps == 0)
  {
    kinetics[n].count_steps = 1;
    kinetics[n].steps[0] = 1.0;
  }
/*
 *   set defaults for moles
 */
  for (i = 0; i < kinetics[n].count_comps; i++)
  {
    if (kinetics[n].comps[i].m0 < 0)
    {
      if (kinetics[n].comps[i].m < 0)
      {
	kinetics[n].comps[i].m0 = 1;
      }
      else
      {
	kinetics[n].comps[i].m0 = kinetics[n].comps[i].m;
      }
    }
    if (kinetics[n].comps[i].m < 0)
    {
      kinetics[n].comps[i].m = kinetics[n].comps[i].m0;
    }
  }
  return (return_value);
}

/* ---------------------------------------------------------------------- */
LDBLE *
read_list_doubles (char **ptr, int *count_doubles)
/* ---------------------------------------------------------------------- */
{
/*
 *   Reads a list of LDBLE numbers until end of line is reached or 
 *   a LDBLE can not be read from a token.
 *
 *      Arguments:
 *         ptr    entry: points to line to read from
 *                exit:  points to next non-LDBLE token or end of line
 *
 *         count_doubles exit: number of LDBLEs read
 *
 *      Returns:
 *         pointer to a list of count_doubles LDBLEs.
 */

  LDBLE *LDBLE_list;
  char token[MAX_LENGTH];
  LDBLE value;
  char *ptr_save;
  int l;

  LDBLE_list = (LDBLE *) PHRQ_malloc (sizeof (LDBLE));
  if (LDBLE_list == NULL)
    malloc_error ();
  *count_doubles = 0;

  ptr_save = *ptr;
  while (copy_token (token, ptr, &l) != EMPTY)
  {
    if (sscanf (token, SCANFORMAT, &value) == 1)
    {
      *count_doubles = *count_doubles + 1;
      LDBLE_list =
	(LDBLE *) PHRQ_realloc (LDBLE_list,
				(size_t) (*count_doubles) * sizeof (LDBLE));
      if (LDBLE_list == NULL)
	malloc_error ();
      LDBLE_list[(*count_doubles) - 1] = value;
      ptr_save = *ptr;
    }
    else
    {
      *ptr = ptr_save;
      break;
    }
  }
  return (LDBLE_list);
}

/* ---------------------------------------------------------------------- */
int *
read_list_ints (char **ptr, int *count_ints, int positive)
/* ---------------------------------------------------------------------- */
{
/*
 *   Reads a list of int numbers until end of line is reached or 
 *   an int can not be read from a token.
 *
 *      Arguments:
 *         ptr    entry: points to line to read from
 *                exit:  points to next non-int token or end of line
 *
 *         count_ints exit: number of LDBLEs read
 *
 *         positive  entry: if TRUE, expects to read only positive integers
 *
 *      Returns:
 *         pointer to a list of count_ints ints.
 */
  int *int_list;
  char token[MAX_LENGTH];
  int value;
  int l;
  char *ptr_save;

  int_list = (int *) PHRQ_malloc (sizeof (int));
  if (int_list == NULL)
    malloc_error ();
  *count_ints = 0;

  ptr_save = *ptr;
  while (copy_token (token, ptr, &l) != EMPTY)
  {
    if (sscanf (token, "%d", &value) == 1)
    {
      (*count_ints)++;
      int_list =
	(int *) PHRQ_realloc (int_list,
			      (size_t) (*count_ints) * sizeof (int));
      if (int_list == NULL)
	malloc_error ();
      int_list[(*count_ints) - 1] = value;
      if (value <= 0 && positive == TRUE)
      {
	error_msg ("Expected an integer greater than zero.", CONTINUE);
	error_msg (line_save, CONTINUE);
	input_error++;
      }
      ptr_save = *ptr;
    }
    else
    {
      *ptr = ptr_save;
      break;
    }
  }
  return (int_list);
}

/* ---------------------------------------------------------------------- */
int *
read_list_ints_range (char **ptr, int *count_ints, int positive,
		      int *int_list)
/* ---------------------------------------------------------------------- */
{
/*
 *   Reads a list of int numbers until end of line is reached or 
 *   an int can not be read from a token.
 *
 *      Arguments:
 *         ptr    entry: points to line to read from
 *                exit:  points to next non-int token or end of line
 *
 *         count_ints entry: number of ints already in list
 *
 *         positive  entry: if TRUE, expects to read only positive integers
 *
 *      Returns:
 *         pointer to a list of count_ints ints 
 */
  char token[MAX_LENGTH];
  int value, value1, value2;
  int i, l;
  char *ptr_save;

  if (int_list == NULL)
  {
    int_list = (int *) PHRQ_malloc (sizeof (int));
    if (int_list == NULL)
      malloc_error ();
    *count_ints = 0;
  }
  ptr_save = *ptr;
  while (copy_token (token, ptr, &l) != EMPTY)
  {
    if (sscanf (token, "%d", &value) == 1)
    {
      /* Read an integer */
      (*count_ints)++;
      int_list =
	(int *) PHRQ_realloc (int_list,
			      (size_t) (*count_ints) * sizeof (int));
      if (int_list == NULL)
	malloc_error ();
      int_list[(*count_ints) - 1] = value;
      if (value <= 0 && positive == TRUE)
      {
	error_msg ("Expected an integer greater than zero.", CONTINUE);
	error_msg (line_save, CONTINUE);
	input_error++;
      }
      /* Read range of integers */
      if (replace ("-", " ", token) == TRUE)
      {
	if (sscanf (token, "%d %d", &value1, &value2) != 2)
	{
	  error_msg ("Expected an integer range n-m.", CONTINUE);
	  error_msg (line_save, CONTINUE);
	  input_error++;
	}
	else if (value2 < value1)
	{
	  error_msg ("Expected an integer range n-m, with n <= m.", CONTINUE);
	  error_msg (line_save, CONTINUE);
	  input_error++;
	}
	else if (value2 <= 0 && positive == TRUE)
	{
	  error_msg ("Expected an integer greater than zero.", CONTINUE);
	  error_msg (line_save, CONTINUE);
	  input_error++;
	}
	else
	{
	  for (i = value1 + 1; i <= value2; i++)
	  {
	    (*count_ints)++;
	    int_list =
	      (int *) PHRQ_realloc (int_list,
				    (size_t) (*count_ints) * sizeof (int));
	    if (int_list == NULL)
	      malloc_error ();
	    int_list[(*count_ints) - 1] = i;
	  }
	}
      }
      ptr_save = *ptr;
    }
    else
    {
      *ptr = ptr_save;
      break;
    }
  }
  return (int_list);
}

/* ---------------------------------------------------------------------- */
int *
read_list_t_f (char **ptr, int *count_ints)
/* ---------------------------------------------------------------------- */
{
/*
 *   Reads a list of true and false until end of line is reached or 
 *   until non- t or f is found
 *
 *      Arguments:
 *         ptr    entry: points to line to read from
 *                exit:  points to next non-int token or end of line
 *
 *         count_ints exit: number of LDBLEs read
 *
 *         positive  entry: if TRUE, expects to read only positive integers
 *
 *      Returns:
 *         pointer to a list of count_ints ints.
 */
  int *int_list;
  char token[MAX_LENGTH];
  int value;
  int l;

  int_list = (int *) PHRQ_malloc (sizeof (int));
  if (int_list == NULL)
    malloc_error ();
  *count_ints = 0;

  while (copy_token (token, ptr, &l) != EMPTY)
  {
    str_tolower (token);
    if (token[0] == 't')
    {
      value = TRUE;
    }
    else if (token[0] == 'f')
    {
      value = FALSE;
    }
    else
    {
      error_msg ("Expected TRUE or FALSE.", CONTINUE);
      error_msg (line_save, CONTINUE);
      input_error++;
      break;
    }
    (*count_ints)++;
    int_list =
      (int *) PHRQ_realloc (int_list, (size_t) (*count_ints) * sizeof (int));
    if (int_list == NULL)
      malloc_error ();
    int_list[(*count_ints) - 1] = value;
  }
  return (int_list);
}

/* ---------------------------------------------------------------------- */
int
read_log_k_only (char *ptr, LDBLE * log_k)
/* ---------------------------------------------------------------------- */
{
/*
 *   Read log k
 */
  *log_k = 0.0;
  replace ("=", " ", ptr);
  if (sscanf (ptr, SCANFORMAT, log_k) < 1)
  {
    input_error++;
    error_msg ("Expecting log k.", CONTINUE);
    return (ERROR);
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
read_delta_h_only (char *ptr, LDBLE * delta_h, DELTA_H_UNIT * units)
/* ---------------------------------------------------------------------- */
{
  int j, l, kilo, joul;
  char token[MAX_LENGTH];
/*
 *   Read delta H
 */
  *delta_h = 0.0;
  replace ("=", " ", ptr);
  j = copy_token (token, &ptr, &l);
  if (j == EMPTY)
  {
    input_error++;
    error_msg ("Expecting numeric value for delta H.", CONTINUE);
    return (ERROR);
  }
  if (sscanf (token, SCANFORMAT, delta_h) < 1)
  {
    input_error++;
    error_msg ("Expecting numeric value for delta H.", CONTINUE);
    return (ERROR);
  }
/*
 *   Read delta H units
 */
  j = copy_token (token, &ptr, &l);
  *units = kjoules;
  kilo = TRUE;
  joul = TRUE;
  if (j == EMPTY)
  {
    return (OK);
  }
  if (j == UPPER || j == LOWER)
  {
    str_tolower (token);
    if (strstr (token, "k") != token)
    {
      /* convert to kilo */
      kilo = FALSE;
      *delta_h /= 1000.;
    }
    if (strstr (token, "c") != NULL)
    {
      /* convert to joules */
      *delta_h *= JOULES_PER_CALORIE;
      joul = FALSE;
    }
  }
  if (kilo == TRUE && joul == TRUE)
  {
    *units = kjoules;
  }
  else if (kilo == FALSE && joul == TRUE)
  {
    *units = joules;
  }
  else if (kilo == TRUE && joul == FALSE)
  {
    *units = kcal;
  }
  else if (kilo == FALSE && joul == FALSE)
  {
    *units = cal;
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
read_analytical_expression_only (char *ptr, LDBLE * log_k)
/* ---------------------------------------------------------------------- */
{
  int j;
/*
 *   Read analytical expression
 */
  for (j = 0; j < 5; j++)
  {
    log_k[j] = 0.0;
  }
  j =
    sscanf (ptr, SCANFORMAT SCANFORMAT SCANFORMAT SCANFORMAT SCANFORMAT,
	    &(log_k[0]), &(log_k[1]), &(log_k[2]), &(log_k[3]), &(log_k[4]));
  if (j < 1)
  {
    input_error++;
    error_msg ("Expecting numeric values for analytical expression.",
	       CONTINUE);
    return (ERROR);
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
read_incremental_reactions (void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Define flow only
 *
 *      Arguments:
 *         none
 *
 *      Returns:
 *         KEYWORD if keyword encountered, input_error may be incremented if
 *                    a keyword is encountered in an unexpected position
 *         EOF     if eof encountered while reading mass balance concentrations
 *         ERROR   if error occurred reading data
 *
 */
  int j, l;
  char *ptr;
  char token[MAX_LENGTH];

  ptr = line;
  /* read keyword */
  copy_token (token, &ptr, &l);

  /* read true or false */
  incremental_reactions = get_true_false (ptr, TRUE);
/*
 *   find next keyword
 */
  while ((j =
	  check_line ("Subroutine Read", FALSE, TRUE, TRUE,
		      FALSE)) != KEYWORD)
  {
    /* empty, eof, keyword, print */
    if (j == EOF)
      return (EOF);
    sprintf (error_string, "Unknown input: %s", line);
    error_msg (error_string, CONTINUE);
    input_error++;
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
read_master_species (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Reads master species data from data file or input file
 */
  int j, i, l;
  char *ptr, *ptr1;
  LDBLE z;
  struct element *elts_ptr;
  struct species *s_ptr;
  char token[MAX_LENGTH], token1[MAX_LENGTH];

  elts_ptr = NULL;
  for (;;)
  {
    j = check_line ("Master species", FALSE, TRUE, TRUE, TRUE);
    if (j == EOF || j == KEYWORD)
    {
      break;
    }
/*
 *   Get element name with valence, allocate space, store
 */
    ptr = line;
/*
 *   Get element name and save pointer to character string
 */
    if (copy_token (token, &ptr, &l) != UPPER && token[0] != '[')
    {
      parse_error++;
      error_msg ("Reading element for master species.", CONTINUE);
      error_msg (line_save, CONTINUE);
      continue;
    }
    /*
       if (token[0] == '[') {
       ptr1 = token;
       get_elt(&ptr, element, &l);
       strcpy(token, element);
       }
     */
    replace ("(+", "(", token);
/*
 *   Delete master if it exists
 */
    master_delete (token);
/*
 *   Increase pointer array, if necessary,  and malloc space
 */
    if (count_master >= max_master)
    {
      space ((void **) ((void *) &master), count_master + 1, &max_master,
	     sizeof (struct master *));
    }
    master[count_master] = master_alloc ();
/*
 *   Set type to AQ
 */
    master[count_master]->type = AQ;
/*
 *   Save element name
 */
    master[count_master]->elt = element_store (token);
/*
 *   Save pointer to species data for master species
 */
    if ((copy_token (token, &ptr, &l) != UPPER) &&
	token[0] != '[' && (strcmp_nocase_arg1 (token, "e-") != 0))
    {
      parse_error++;
      error_msg ("Reading master species name.", CONTINUE);
      error_msg (line_save, CONTINUE);
      continue;
    }

    s_ptr = s_search (token);
    if (s_ptr != NULL)
    {
      master[count_master]->s = s_ptr;
    }
    else
    {
      ptr1 = token;
      get_token (&ptr1, token1, &z, &l);
      master[count_master]->s = s_store (token1, z, FALSE);
    }

/*
 *   Read alkalinity for species
 */
    copy_token (token, &ptr, &l);
    i = sscanf (token, SCANFORMAT, &master[count_master]->alk);
    if (i != 1)
    {
      input_error++;
      if (elts_ptr != NULL)
      {
	sprintf (error_string,
		 "Expected alkalinity for master species, %s, in master species input.",
		 elts_ptr->name);
      }
      else
      {
	sprintf (error_string,
		 "Expected alkalinity for master species in master species input.");
      }
      error_msg (error_string, CONTINUE);
      continue;
    }
/*
 *   Read default gfw for species
 */
    i = copy_token (token, &ptr, &l);
    if (i == DIGIT)
    {
      sscanf (token, SCANFORMAT, &master[count_master]->gfw);
    }
    else if (i == UPPER)
    {
      master[count_master]->gfw_formula = string_hsave (token);
    }
    else
    {
      input_error++;
      if (elts_ptr != NULL)
      {
	sprintf (error_string,
		 "Expected gram formula weight for master species, %s, in master species input.",
		 elts_ptr->name);
      }
      else
      {
	sprintf (error_string,
		 "Expected gram formula weight for master species in master species input.");
      }
      error_msg (error_string, CONTINUE);
      continue;
    }
/*
 *   MAKE LISTS OF PRIMARY AND SECONDARY MASTER SPECIES
 */
    if (strchr (master[count_master]->elt->name, '(') == NULL)
    {
      master[count_master]->primary = TRUE;
      /* Read gram formula weight for primary */
      if (strcmp (master[count_master]->elt->name, "E") != 0)
      {
	elts_ptr = master[count_master]->elt;
	i = copy_token (token, &ptr, &l);
	if (i == DIGIT)
	{
	  sscanf (token, SCANFORMAT, &elts_ptr->gfw);
	}
	else
	{
	  input_error++;
	  if (elts_ptr != NULL)
	  {
	    sprintf (error_string,
		     "Expected gram formula weight for element, %s.",
		     elts_ptr->name);
	  }
	  else
	  {
	    sprintf (error_string,
		     "Expected gram formula weight for element.");
	  }

	  error_msg (error_string, CONTINUE);
	  continue;
	}
      }
    }
    else
    {
      master[count_master]->primary = FALSE;
    }
    count_master++;
    if (count_master >= max_master)
    {
      space ((void **) ((void *) &master), count_master, &max_master,
	     sizeof (struct master *));
    }

  }
  return (j);
}

/* ---------------------------------------------------------------------- */
int
read_mix (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Reads mixing fractions
 */
  int n, n_user, n_user_end;
  int return_value;
  int count_comps;
  int n_solution;
  LDBLE fraction;
  int j, i, l;
  char *ptr;
  char token[MAX_LENGTH];
  char *description;
  struct mix *mix_ptr;
/*
 *   Read mix number
 */
  ptr = line;
  read_number_description (ptr, &n_user, &n_user_end, &description);
/*
 *   Update list of mixtures
 */
  mix_ptr = mix_search (n_user, &n, FALSE);
  if (mix_ptr != NULL)
  {
    mix_free (&mix[n]);
  }
  else
  {
    n = count_mix++;
    mix =
      (struct mix *) PHRQ_realloc (mix,
				   (size_t) count_mix * sizeof (struct mix));
    if (mix == NULL)
      malloc_error ();
  }
/*
 *   Set use data to first read
 */
  if (use.mix_in == FALSE)
  {
    use.mix_in = TRUE;
    use.n_mix_user = n_user;
  }
/*
 *   Defaults
 */
  mix[n].description = description;
  mix[n].n_user = n_user;
  mix[n].n_user_end = n_user_end;
  mix[n].comps = NULL;
  count_comps = 0;
/*
 *   Read mixture data
 */
  for (;;)
  {
    return_value = check_line ("Mixture data", FALSE, TRUE, TRUE, TRUE);
    /* empty, eof, keyword, print */
    if (return_value == EOF || return_value == KEYWORD)
    {
      break;
    }
    ptr = line;
/*
 *   Read n_user
 */
    i = copy_token (token, &ptr, &l);
    if (i == DIGIT)
    {
      sscanf (token, "%d ", &n_solution);
    }
    else
    {
      input_error++;
      error_msg ("Expected a solution number in mix input.", CONTINUE);
      error_msg (line_save, CONTINUE);
      continue;
    }
/*
 *   Read fraction for solution
 */
    copy_token (token, &ptr, &l);
    j = sscanf (token, SCANFORMAT, &fraction);
    if (j != 1)
    {
      input_error++;
      error_msg ("Expected a mixing fraction.", CONTINUE);
      error_msg (line_save, CONTINUE);
      continue;
    }
/*
 *   Malloc space
 */
    count_comps++;
    if (mix[n].comps == NULL)
    {
      mix[n].comps =
	(struct mix_comp *) PHRQ_malloc ((size_t) sizeof (struct mix_comp));
    }
    else
    {
      mix[n].comps = (struct mix_comp *) PHRQ_realloc (mix[n].comps,
						       (size_t) count_comps *
						       sizeof (struct
							       mix_comp));
    }
    if (mix[n].comps == NULL)
      malloc_error ();
/*
 *   Save data
 */
    mix[n].comps[count_comps - 1].n_solution = n_solution;
    mix[n].comps[count_comps - 1].fraction = fraction;
  }
  if (count_comps <= 0)
  {
    input_error++;
    error_msg
      ("Must define at least one solution number and mixing fraction for MIX input.",
       CONTINUE);
  }
  mix[n].count_comps = count_comps;
  return (return_value);
}

/* ---------------------------------------------------------------------- */
int
read_number_description (char *ptr, int *n_user,
			 int *n_user_end, char **description)
/* ---------------------------------------------------------------------- */
{
  int l, n;
  char token[MAX_LENGTH];
  char *ptr1;
/*
 *   Read user number
 */
  copy_token (token, &ptr, &l);
  ptr1 = ptr;
  if (copy_token (token, &ptr, &l) != DIGIT)
  {
/*
		sprintf(error_string, "No number given with %s keyword, "
			"%s 1 assumed.", str, str);
			warning_msg(error_string);
 */
    *n_user = 1;
    *n_user_end = 1;
  }
  else if (strstr (token, "-") != NULL)
  {
    replace ("-", " ", token);
    n = sscanf (token, "%d%d", n_user, n_user_end);
    if (n != 2)
    {
      if (next_keyword >= 0)
      {
	sprintf (error_string, "Reading number range for %s.",
		 keyword[next_keyword].name);
      }
      else
      {
	sprintf (error_string, "Reading number range for keyword.");
      }
      error_msg (error_string, CONTINUE);
      input_error++;
    }
    ptr1 = ptr;
  }
  else
  {
    sscanf (token, "%d", n_user);
    *n_user_end = *n_user;
    ptr1 = ptr;
  }
/*
 *   Read description
 */
  for (; isspace ((int) ptr1[0]); ptr1++);
  *description = string_duplicate (ptr1);
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
read_phases (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Read data for phases, parse equations
 */
  int j, i, l;
  int association;
  char *ptr;
  char token[MAX_LENGTH];
  char token1[MAX_LENGTH];
  struct phase *phase_ptr;
  struct elt_list *next_elt;
  struct rxn_token *token_ptr;

  int return_value, opt, opt_save;
  char *next_char;
  const char *opt_list[] = {
    "no_check",			/* 0 */
    "check",			/* 1 */
    "log_k",			/* 2 */
    "logk",			/* 3 */
    "delta_h",			/* 4 */
    "deltah",			/* 5 */
    "analytical_expression",	/* 6 */
    "a_e",			/* 7 */
    "ae",			/* 8 */
    "add_logk",			/* 9 */
    "add_log_k",		/* 10 */
    "add_constant"		/* 11 */
  };
  int count_opt_list = 12;

  association = FALSE;
/*
 *   Read eqn from file and call parser
 */
  opt_save = OPTION_DEFAULT;
  return_value = UNKNOWN;
  phase_ptr = NULL;
  for (;;)
  {
    opt = get_option (opt_list, count_opt_list, &next_char);
    if (opt == OPTION_DEFAULT)
    {
      opt = opt_save;
    }
    switch (opt)
    {
    case OPTION_EOF:		/* end of file */
      return_value = EOF;
      break;
    case OPTION_KEYWORD:	/* keyword */
      return_value = KEYWORD;
      break;
    case OPTION_ERROR:
      input_error++;
      error_msg ("Unknown input in PHASES keyword.", CONTINUE);
      error_msg (line_save, CONTINUE);
      break;
    case 0:			/* no_check */
      if (phase_ptr == NULL)
	break;
      phase_ptr->check_equation = FALSE;
      break;
    case 1:			/* check */
      if (phase_ptr == NULL)
	break;
      phase_ptr->check_equation = TRUE;
      break;
    case 2:			/* log_k */
    case 3:			/* logk */
      if (phase_ptr == NULL)
	break;
      read_log_k_only (next_char, &phase_ptr->logk[0]);
      opt_save = OPTION_DEFAULT;
      break;
    case 4:			/* delta_h */
    case 5:			/* deltah */
      if (phase_ptr == NULL)
	break;
      read_delta_h_only (next_char, &phase_ptr->logk[1],
			 &phase_ptr->original_units);
      opt_save = OPTION_DEFAULT;
      break;
    case 6:			/* analytical_expression */
    case 7:			/* a_e */
    case 8:			/* ae */
      if (phase_ptr == NULL)
	break;
      read_analytical_expression_only (next_char, &(phase_ptr->logk[2]));
      opt_save = OPTION_DEFAULT;
      break;
    case 9:			/* add_logk */
    case 10:			/* add_log_k */
      if (phase_ptr == NULL)
	break;
      if (phase_ptr->count_add_logk == 0)
      {
	phase_ptr->add_logk =
	  (struct name_coef *) PHRQ_malloc (sizeof (struct name_coef));
	if (phase_ptr->add_logk == NULL)
	  malloc_error ();
      }
      else
      {
	phase_ptr->add_logk =
	  (struct name_coef *) PHRQ_realloc (phase_ptr->add_logk,
					     (size_t) ((phase_ptr->
							count_add_logk +
							1) *
						       sizeof (struct
							       name_coef)));
	if (phase_ptr->add_logk == NULL)
	  malloc_error ();
      }
      /* read name */
      if (copy_token (token, &next_char, &i) == EMPTY)
      {
	input_error++;
	sprintf (error_string, "Expected the name of a NAMED_EXPRESSION.");
	error_msg (error_string, CONTINUE);
	break;
      }
      phase_ptr->add_logk[phase_ptr->count_add_logk].name =
	string_hsave (token);
      /* read coef */
      i =
	sscanf (next_char, SCANFORMAT,
		&phase_ptr->add_logk[phase_ptr->count_add_logk].coef);
      if (i <= 0)
      {
	phase_ptr->add_logk[phase_ptr->count_add_logk].coef = 1;
      }
      phase_ptr->count_add_logk++;
      opt_save = OPTION_DEFAULT;
      break;
    case 11:			/* add_constant */
      if (phase_ptr == NULL)
	break;
      if (phase_ptr->count_add_logk == 0)
      {
	phase_ptr->add_logk =
	  (struct name_coef *) PHRQ_malloc (sizeof (struct name_coef));
	if (phase_ptr->add_logk == NULL)
	  malloc_error ();
      }
      else
      {
	phase_ptr->add_logk =
	  (struct name_coef *) PHRQ_realloc (phase_ptr->add_logk,
					     (size_t) ((phase_ptr->
							count_add_logk +
							1) *
						       sizeof (struct
							       name_coef)));
	if (phase_ptr->add_logk == NULL)
	  malloc_error ();
      }
      i =
	sscanf (next_char, SCANFORMAT,
		&phase_ptr->add_logk[phase_ptr->count_add_logk].coef);
      if (i <= 0)
      {
	input_error++;
	sprintf (error_string,
		 "Expected the constant to add for log_K definition.");
	error_msg (error_string, CONTINUE);
	break;
      }
      /* set name */
      phase_ptr->add_logk[phase_ptr->count_add_logk].name =
	string_hsave ("XconstantX");
      /* read coef */
      phase_ptr->count_add_logk++;
      opt_save = OPTION_DEFAULT;
      break;
    case OPTION_DEFAULT:
/*
 *   Get element name and save pointer to character string
 */
      phase_ptr = NULL;
      ptr = line;
      copy_token (token, &ptr, &l);
/*
 *   Get and parse equation
 */
      j = check_line ("Phase equation", FALSE, TRUE, TRUE, TRUE);
      if (j == EOF || j == KEYWORD)
      {
	return_value = j;
	break;
      }
      if (parse_eq (line, &next_elt, association) == ERROR)
      {
	parse_error++;
	error_msg ("Parsing equation.", CONTINUE);
	error_msg (line_save, CONTINUE);
	break;
      }
      phase_ptr = phase_store (token);
/*
 *   Get pointer to each species in the reaction, store new species if necessary
 */
      strcpy (token1, trxn.token[0].name);
      replace ("(g)", "", token1);
      replace ("(s)", "", token1);
      replace ("(G)", "", token1);
      replace ("(S)", "", token1);
      phase_ptr->formula = string_hsave (token1);
      for (i = 1; i < count_trxn; i++)
      {
	if ((strstr (trxn.token[i].name, "(s)") == NULL) &&
	    (strstr (trxn.token[i].name, "(g)") == NULL) &&
	    (strstr (trxn.token[i].name, "(S)") == NULL) &&
	    (strstr (trxn.token[i].name, "(G)") == NULL))
	{
	  strcpy (token1, trxn.token[i].name);
	  replace ("(aq)", "", token1);
	  replace ("(AQ)", "", token1);
	  replace ("H2O(l)", "H2O", token1);
	  replace ("(H2O(L)", "H2O", token1);
	  trxn.token[i].s = s_store (token1, trxn.token[i].z, FALSE);
	}
	else
	{
	  trxn.token[i].s = NULL;
	}
      }
/*
 *   Save element list
 */
      phase_ptr->next_elt = next_elt;
/*
 *   Malloc space for phase reaction
 */
      phase_ptr->rxn = rxn_alloc (count_trxn + 1);
/*
 *   Copy reaction to reaction for phase, first token (token[0]) is not used
 *   except to check that coef of phase formula = 1.0
 */
      token_ptr = phase_ptr->rxn->token;
      /* token_ptr[0].coef=0; */
      token_ptr[0].coef = trxn.token[0].coef;
      token_ptr[0].s = trxn.token[1].s;
      for (i = 1; i < count_trxn; i++)
      {
	token_ptr[i].name = NULL;
	token_ptr[i].s = trxn.token[i].s;
	token_ptr[i].coef = trxn.token[i].coef;
	if (token_ptr[i].s == NULL)
	{
	  token_ptr[i].name = trxn.token[i].name;
	}
      }
      token_ptr[0].name = trxn.token[1].name;
      /*
         token_ptr[0].name=phase_ptr->name;
         token_ptr[0].s=NULL;
       */
      token_ptr[i].s = NULL;
      token_ptr[i].name = NULL;
/*
 *   Set type for phase
 */
      phase_ptr->type = SOLID;
      opt_save = OPTION_DEFAULT;
      break;
    }
    if (return_value == EOF || return_value == KEYWORD)
      break;
  }
  return (return_value);
}

#ifdef SKIP
/* ---------------------------------------------------------------------- */
int
read_pure_phases (void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads pure phase data
 *
 *      Arguments:
 *         none
 *
 *      Returns:
 *         KEYWORD if keyword encountered, input_error may be incremented if
 *                    a keyword is encountered in an unexpected position
 *         EOF     if eof encountered while reading mass balance concentrations
 *         ERROR   if error occurred reading data
 *
 */
  int j, n, l, return_value;
  int count_pure_phases;
  int n_user, n_user_end;
  char *ptr;
  char *description;
  char token[MAX_LENGTH];

  ptr = line;
/*
 *   Read pp_assemblage number
 */
  read_number_description (ptr, &n_user, &n_user_end, &description);
/*
 *   Find pp_assemblage or realloc space for pp_assemblage
 */

  if (pp_assemblage_search (n_user, &n) != NULL)
  {
    pp_assemblage_free (&pp_assemblage[n]);
  }
  else
  {
    n = count_pp_assemblage++;
    space ((void **) ((void *) &pp_assemblage), count_pp_assemblage,
	   &max_pp_assemblage, sizeof (struct pp_assemblage));
  }
/*
 *   Set use data to first read
 */
  if (use.pp_assemblage_in == FALSE)
  {
    use.pp_assemblage_in = TRUE;
    use.n_pp_assemblage_user = n_user;
  }

  pp_assemblage_init (&(pp_assemblage[n]), n_user, n_user_end, description);
  free_check_null (description);
/*
 *   Read phases
 */
  count_pure_phases = 0;
  for (;;)
  {
    return_value = check_line ("Pure phase data", FALSE, TRUE, TRUE, TRUE);
    /* empty, eof, keyword, print */
    if (return_value == EOF || return_value == KEYWORD)
    {
      break;
    }
/*
 *   Make space, set default
 */
    count_pure_phases++;
    pp_assemblage[n].pure_phases =
      (struct pure_phase *) PHRQ_realloc (pp_assemblage[n].pure_phases,
					  (size_t) (count_pure_phases) *
					  sizeof (struct pure_phase));
    if (pp_assemblage[n].pure_phases == NULL)
      malloc_error ();
    pp_assemblage[n].pure_phases[count_pure_phases - 1].si = 0.0;
    pp_assemblage[n].pure_phases[count_pure_phases - 1].add_formula = NULL;
    pp_assemblage[n].pure_phases[count_pure_phases - 1].moles = 10.0;
    pp_assemblage[n].pure_phases[count_pure_phases - 1].delta = 0.0;
    pp_assemblage[n].pure_phases[count_pure_phases - 1].initial_moles = 0.0;
    pp_assemblage[n].pure_phases[count_pure_phases - 1].dissolve_only = FALSE;
/*
 *   Read name
 */
    ptr = line;
    copy_token (token, &ptr, &l);
    pp_assemblage[n].pure_phases[count_pure_phases - 1].name =
      string_hsave (token);
    if ((j = copy_token (token, &ptr, &l)) == EMPTY)
      continue;
/*
 *   Read saturation index
 */
    j = sscanf (token, SCANFORMAT, &dummy);
    pp_assemblage[n].pure_phases[count_pure_phases - 1].si = (LDBLE) dummy;
    if (j != 1)
    {
      error_msg ("Expected saturation index.", CONTINUE);
      error_msg (line_save, CONTINUE);
      input_error++;
      continue;
    }
/*
 *   Adding a reaction to the phase boundary
 */
    if ((j = copy_token (token, &ptr, &l)) == EMPTY)
      continue;
    if (j == UPPER || j == LOWER)
    {
      pp_assemblage[n].pure_phases[count_pure_phases - 1].add_formula =
	string_hsave (token);
      j = copy_token (token, &ptr, &l);
    }
/*
 *   Read amount
 */
    if (j == EMPTY)
      continue;
    j = sscanf (token, SCANFORMAT, &dummy);
    pp_assemblage[n].pure_phases[count_pure_phases - 1].moles = (LDBLE) dummy;
    if (j != 1)
    {
      error_msg ("Expected amount of mineral.", CONTINUE);
      error_msg (line_save, CONTINUE);
      input_error++;
      continue;
    }
    if ((j = copy_token (token, &ptr, &l)) == EMPTY)
      continue;
    str_tolower (token);
    if (strstr (token, "d") == token)
    {
      pp_assemblage[n].pure_phases[count_pure_phases - 1].dissolve_only =
	TRUE;
    }
    else
    {
      error_msg ("Unexpected data at end of equilibrium-phase definition.",
		 CONTINUE);
      input_error++;
      continue;
    }
  }
  pp_assemblage[n].count_comps = count_pure_phases;
/*
 *   Sort phases by name (lowercase)
 */
  qsort (pp_assemblage[n].pure_phases,
	 (size_t) count_pure_phases,
	 (size_t) sizeof (struct pure_phase), pure_phase_compare);

  return (return_value);
}
#endif
/* ---------------------------------------------------------------------- */
int
read_pure_phases (void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads pure phase data
 *
 *      Arguments:
 *         none
 *
 *      Returns:
 *         KEYWORD if keyword encountered, input_error may be incremented if
 *                    a keyword is encountered in an unexpected position
 *         EOF     if eof encountered while reading mass balance concentrations
 *         ERROR   if error occurred reading data
 *
 */
  int j, n, l, return_value;
  int count_pure_phases;
  int n_user, n_user_end;
  char *ptr;
  char *description;
  char token[MAX_LENGTH];
  int opt, opt_save;
  char *next_char;
  const char *opt_list[] = {
    "force_equality"		/* 0 */
  };
  int count_opt_list = 1;

  ptr = line;
  /*
   *   Read pp_assemblage number
   */
  read_number_description (ptr, &n_user, &n_user_end, &description);
  /*
   *   Find pp_assemblage or realloc space for pp_assemblage
   */

  if (pp_assemblage_search (n_user, &n) != NULL)
  {
    pp_assemblage_free (&pp_assemblage[n]);
  }
  else
  {
    n = count_pp_assemblage++;
    space ((void **) ((void *) &pp_assemblage), count_pp_assemblage,
	   &max_pp_assemblage, sizeof (struct pp_assemblage));
  }
  /*
   *   Set use data to first read
   */
  if (use.pp_assemblage_in == FALSE)
  {
    use.pp_assemblage_in = TRUE;
    use.n_pp_assemblage_user = n_user;
  }

  pp_assemblage_init (&(pp_assemblage[n]), n_user, n_user_end, description);
  free_check_null (description);
  count_pure_phases = 0;
  /*
   *  Read equilibrium phase data
   */
  opt_save = OPTION_DEFAULT;
  return_value = UNKNOWN;
  for (;;)
  {
    opt = get_option (opt_list, count_opt_list, &next_char);
    if (opt == OPTION_DEFAULT)
    {
      opt = opt_save;
    }
    switch (opt)
    {
    case OPTION_EOF:		/* end of file */
      return_value = EOF;
      break;
    case OPTION_KEYWORD:	/* keyword */
      return_value = KEYWORD;
      break;
    case OPTION_ERROR:
      input_error++;
      error_msg ("Unknown input in PHASES keyword.", CONTINUE);
      error_msg (line_save, CONTINUE);
      break;
    case 0:			/* force_equality */
      if (count_pure_phases == 0)
      {
	error_msg
	  ("Force_equality defined before equilibrium phase has been defined.",
	   CONTINUE);
	error_msg (line_save, CONTINUE);
	input_error++;
      }
      else
      {
	pp_assemblage[n].pure_phases[count_pure_phases - 1].force_equality =
	  get_true_false (next_char, TRUE);
      }
      break;
    case OPTION_DEFAULT:
      /*
       *   Make space, set default
       */
      count_pure_phases++;
      pp_assemblage[n].pure_phases =
	(struct pure_phase *) PHRQ_realloc (pp_assemblage[n].pure_phases,
					    (size_t) (count_pure_phases) *
					    sizeof (struct pure_phase));
      if (pp_assemblage[n].pure_phases == NULL)
	malloc_error ();
      pp_assemblage[n].pure_phases[count_pure_phases - 1].si = 0.0;
      pp_assemblage[n].pure_phases[count_pure_phases - 1].add_formula = NULL;
      pp_assemblage[n].pure_phases[count_pure_phases - 1].moles = 10.0;
      pp_assemblage[n].pure_phases[count_pure_phases - 1].delta = 0.0;
      pp_assemblage[n].pure_phases[count_pure_phases - 1].initial_moles = 0.0;
      pp_assemblage[n].pure_phases[count_pure_phases - 1].force_equality =
	FALSE;
      pp_assemblage[n].pure_phases[count_pure_phases - 1].dissolve_only =
	FALSE;
      /*
       *   Read name
       */
      ptr = line;
      copy_token (token, &ptr, &l);
      pp_assemblage[n].pure_phases[count_pure_phases - 1].name =
	string_hsave (token);
      if ((j = copy_token (token, &ptr, &l)) == EMPTY)
	continue;
      /*
       *   Read saturation index
       */
      j = sscanf (token, SCANFORMAT, &dummy);
      pp_assemblage[n].pure_phases[count_pure_phases - 1].si = (LDBLE) dummy;
      if (j != 1)
      {
	error_msg ("Expected saturation index.", CONTINUE);
	error_msg (line_save, CONTINUE);
	input_error++;
	continue;
      }
      /*
       *   Adding a reaction to the phase boundary
       */
      if ((j = copy_token (token, &ptr, &l)) == EMPTY)
	continue;
      if (j == UPPER || j == LOWER)
      {
	pp_assemblage[n].pure_phases[count_pure_phases - 1].add_formula =
	  string_hsave (token);
	j = copy_token (token, &ptr, &l);
      }
      /*
       *   Read amount
       */
      if (j == EMPTY)
	continue;
      j = sscanf (token, SCANFORMAT, &dummy);
      pp_assemblage[n].pure_phases[count_pure_phases - 1].moles =
	(LDBLE) dummy;
      if (j != 1)
      {
	error_msg ("Expected amount of mineral.", CONTINUE);
	error_msg (line_save, CONTINUE);
	input_error++;
	continue;
      }
      if ((j = copy_token (token, &ptr, &l)) == EMPTY)
	continue;
      str_tolower (token);
      if (strstr (token, "d") == token)
      {
	pp_assemblage[n].pure_phases[count_pure_phases - 1].dissolve_only =
	  TRUE;
      }
      else
      {
	error_msg ("Unexpected data at end of equilibrium-phase definition.",
		   CONTINUE);
	input_error++;
	continue;
      }
      break;
    }
    if (return_value == EOF || return_value == KEYWORD)
      break;
  }
  pp_assemblage[n].count_comps = count_pure_phases;
  /*
   *   Sort phases by name (lowercase)
   */
  qsort (pp_assemblage[n].pure_phases,
	 (size_t) count_pure_phases,
	 (size_t) sizeof (struct pure_phase), pure_phase_compare);

  return (return_value);
}

/* ---------------------------------------------------------------------- */
int
read_reaction (void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads reaction data
 *
 *      Arguments:
 *         none
 *
 *      Returns:
 *         KEYWORD if keyword encountered, input_error may be incremented if
 *                    a keyword is encountered in an unexpected position
 *         EOF     if eof encountered while reading mass balance concentrations
 *         ERROR   if error occurred reading data
 *
 */
/*
 *   Read reaction
 */
  int l;
  char *ptr;
  char *description;
  char token[MAX_LENGTH];
  int n, return_value;
  int n_user, n_user_end;
  struct irrev *irrev_ptr;
/*
 *   Read reaction number
 */
  ptr = line;
  read_number_description (ptr, &n_user, &n_user_end, &description);

/*
 *   Read reaction or realloc space for irrev reaction
 */
  irrev_ptr = irrev_search (n_user, &n);
  if (irrev_ptr != NULL)
  {
    irrev_free (&irrev[n]);
  }
  else
  {
    n = count_irrev++;
    irrev =
      (struct irrev *) PHRQ_realloc (irrev,
				     (size_t) count_irrev *
				     sizeof (struct irrev));
    if (irrev == NULL)
      malloc_error ();
  }
/*
 *   Set use data to first read
 */
  if (use.irrev_in == FALSE)
  {
    use.irrev_in = TRUE;
    use.n_irrev_user = n_user;
  }
/*
 *   Defaults
 */
  irrev[n].n_user = n_user;
  irrev[n].n_user_end = n_user_end;
  irrev[n].description = description;
  irrev[n].steps = (double *) PHRQ_malloc ((size_t) sizeof (LDBLE));
  if (irrev[n].steps == NULL)
    malloc_error ();
  irrev[n].list =
    (struct name_coef *) PHRQ_malloc ((size_t) sizeof (struct name_coef));
  if (irrev[n].list == NULL)
    malloc_error ();
  irrev[n].count_steps = 0;
  irrev[n].count_list = 0;
  irrev[n].units = string_hsave ("Mol");
  irrev[n].elts = NULL;
/*
 *   Read reaction data
 */
  for (;;)
  {
/*
 *   Read line
 */
    return_value = check_line ("Pure phase data", FALSE, TRUE, TRUE, TRUE);
    /* empty, eof, keyword, print */
    if (return_value == EOF || return_value == KEYWORD)
    {
      break;
    }
    ptr = line;
    copy_token (token, &ptr, &l);
    if (isalpha ((int) token[0]) || (token[0] == '(') || (token[0] == '['))
    {
/*
 *   Read reactant information
 */
      read_reaction_reactants (&(irrev[n]));
    }
    else
    {
/*
 *   Read steps information
 */
      read_reaction_steps (&(irrev[n]));
    }
  }
/*
 *   Default 1 mol of reaction
 */
  if (irrev[n].count_steps == 0)
  {
    irrev[n].count_steps = 1;
    irrev[n].steps[0] = 1.0;
  }
  return (return_value);
}

/* ---------------------------------------------------------------------- */
int
read_reaction_reactants (struct irrev *irrev_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Read reactants, may be a chemical formula or a pure_phase name
 *   followed by relative reaction coefficient, default 1.0.
 *
 */
  int j, l;
  int count_list;
  char token[MAX_LENGTH];
  LDBLE coef;
  char *ptr;
/*
 *   Read one or more reactants
 */
  ptr = line;
  while (copy_token (token, &ptr, &l) != EMPTY)
  {
/*
 *   Store reactant name, default coefficient
 */
    if (isalpha ((int) token[0]) || (token[0] == '(') || (token[0] == '['))
    {
      irrev_ptr->count_list++;
      count_list = irrev_ptr->count_list - 1;
      irrev_ptr->list =
	(struct name_coef *) PHRQ_realloc (irrev_ptr->list,
					   (size_t) (count_list +
						     1) *
					   sizeof (struct name_coef));
      if (irrev_ptr->list == NULL)
	malloc_error ();
      irrev_ptr->list[count_list].name = string_hsave (token);
      irrev_ptr->list[count_list].coef = 1.0;
/*
 *   Store relative coefficient
 */
    }
    else
    {
      j = sscanf (token, SCANFORMAT, &coef);
      if (j == 1)
      {
	count_list = irrev_ptr->count_list - 1;
	irrev_ptr->list[count_list].coef = coef;
      }
      else
      {
	error_msg ("Reading relative coefficient of reactant.", CONTINUE);
	error_msg (line_save, CONTINUE);
	input_error++;
      }
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
read_reaction_steps (struct irrev *irrev_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Read amount(s) of irrev reactions in one of three forms:
 *
 *   6 millimoles in 6 steps   or
 *
 *   1 2 3 4 5 6 millimoles    or
 *
 *   6*1 millimoles
 *   INCREMENTAL_REACTIONS
 */
  int i, j, l, n;
  int count_steps;
  char *ptr;
  char token[MAX_LENGTH], token1[MAX_LENGTH];
  LDBLE step, value;

  ptr = line;
  count_steps = irrev_ptr->count_steps;
/*
 *   Read one or more reaction increments
 */
  for (;;)
  {
    if (copy_token (token, &ptr, &l) == EMPTY)
    {
      return (OK);
    }
/*
 *   Read next step increment
 */
/* begin modif 29 july 2005... */
    if (replace ("*", " ", token) == TRUE)
    {
      if (sscanf (token, "%d" SCANFORMAT, &n, &value) == 2)
      {
	for (i = 0; i < n; i++)
	{
	  count_steps++;
	  irrev_ptr->steps =
	    (double *) PHRQ_realloc (irrev_ptr->steps,
				     (size_t) count_steps * sizeof (LDBLE));
	  if (irrev_ptr->steps == NULL)
	    malloc_error ();
	  irrev_ptr->steps[irrev_ptr->count_steps] = value;
	  irrev_ptr->count_steps = count_steps;
	}
      }
      else
      {
	input_error++;
	error_msg
	  ("Format error in multiple, equal REACTION steps.\nCorrect is (for example): 0.2 4*0.1 2*0.5 0.3\n",
	   CONTINUE);
      }
    }
    else
    {
      j = sscanf (token, SCANFORMAT, &step);
      if (j == 1)
      {
	count_steps++;
	irrev_ptr->steps =
	  (double *) PHRQ_realloc (irrev_ptr->steps,
				   (size_t) count_steps * sizeof (LDBLE));
	if (irrev_ptr->steps == NULL)
	  malloc_error ();
	irrev_ptr->steps[irrev_ptr->count_steps] = step;
	irrev_ptr->count_steps = count_steps;
      }
      else
      {
	break;
      }
    }
/* ...end modif 29 july 2005 */
  }
/*
 *   Read units
 */
  strcpy (token1, token);
  strcat (token1, "/l");
  if (check_units (token1, FALSE, FALSE, NULL, FALSE) == OK)
  {
    replace ("/l", "", token1);
    if (strstr (token1, "Mol") == NULL)
    {
      sprintf (error_string, "Units of steps not in moles, %s.", token);
      error_msg (error_string, CONTINUE);
      input_error++;
      return (ERROR);
    }
    else
    {
      irrev_ptr->units = string_hsave (token1);
    }
    if (copy_token (token, &ptr, &l) == EMPTY)
    {
      return (OK);
    }
  }
/*
 *  Read number of equal increments, store as negative integer
 */
  if (count_steps != 1)
  {
    error_msg
      ("To define equal increments, only one reaction increment should be defined.",
       CONTINUE);
    input_error++;
    return (ERROR);
  }
  do
  {
    j = sscanf (token, "%d", &i);
    if (j == 1 && i > 0)
    {
      irrev_ptr->count_steps = -i;
      return (OK);
    }
    else if (j == 1 && i <= 0)
    {
      break;
    }
  }
  while (copy_token (token, &ptr, &l) != EMPTY);

  error_msg ("Expecting positive number for number of equal "
	     "increments to add.", CONTINUE);
  error_msg (line_save, CONTINUE);
  input_error++;
  return (ERROR);
}

/* ---------------------------------------------------------------------- */
int
read_save (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Reads solution, mix, irreversible reaction, and pure phases to use
 *   in reaction calculation
 */
  int i, l, n, n_user, n_user_end;
  char *ptr;
  char token[MAX_LENGTH];
/*
 *   Read "save"
 */
  ptr = line;
  copy_token (token, &ptr, &l);
/*
 *   Read keyword
 */
  copy_token (token, &ptr, &l);
  check_key (token);
  switch (next_keyword)
  {
  case -1:			/* Have not read line with keyword */
  case 0:			/* End encountered */
  case 1:			/* EOF encountered */
  case 2:			/* Read aqueous model */
  case 3:			/* Read master species */
  case 5:			/* Phases */
  case 7:			/* Reaction */
  case 8:			/* Mix */
  case 9:			/* Use */
  case 10:			/* Save */
  case 11:
  case 12:
  case 14:
  case 15:
  case 17:
  case 18:
  case 21:
  case 22:
  case 23:
  case 24:
  case 25:
  case 30:
  case 31:
  case 32:
  case 33:
  case 34:
  case 35:
  case 36:
  case 37:
  case 38:
  case 39:
    input_error++;
    error_msg
      ("Expecting keyword solution, equilibrium_phases, exchange, surface, gas_phase, or solid_solutions.",
       CONTINUE);
    error_msg (line_save, CONTINUE);
    check_line ("End of save", FALSE, TRUE, TRUE, TRUE);
    /* empty, eof, keyword, print */
    return (ERROR);
  }
/*
 *   Read number
 */
  for (;;)
  {
    i = copy_token (token, &ptr, &l);
    if (i == DIGIT)
    {
      replace ("-", " ", token);
      n = sscanf (token, "%d%d", &n_user, &n_user_end);
      if (n == 1)
      {
	n_user_end = n_user;
      }
      if (n_user < 0)
      {
	error_msg ("Number must be a positive integer.", CONTINUE);
	error_msg (line_save, CONTINUE);
	input_error++;
      }
      break;
    }
    else if (i == EMPTY)
    {
      sprintf (error_string, "No number given, 1 assumed.");
      warning_msg (error_string);
      n_user = 1;
      n_user_end = 1;
      break;
    }
  }
  switch (next_keyword)
  {
  case 4:			/* Solution */
    save.solution = TRUE;
    save.n_solution_user = n_user;
    save.n_solution_user_end = n_user_end;
    break;
  case 6:			/* Pure phases */
  case 26:
  case 27:
  case 28:
  case 29:
    save.pp_assemblage = TRUE;
    save.n_pp_assemblage_user = n_user;
    save.n_pp_assemblage_user_end = n_user_end;
    break;
  case 13:			/* exchange */
    save.exchange = TRUE;
    save.n_exchange_user = n_user;
    save.n_exchange_user_end = n_user_end;
    break;
  case 16:			/* surface */
    save.surface = TRUE;
    save.n_surface_user = n_user;
    save.n_surface_user_end = n_user_end;
    break;
  case 19:			/* gas_phase */
    save.gas_phase = TRUE;
    save.n_gas_phase_user = n_user;
    save.n_gas_phase_user_end = n_user_end;
    break;
  case 40:			/* solid_solutions */
  case 41:			/* solid_solution */
    save.s_s_assemblage = TRUE;
    save.n_s_s_assemblage_user = n_user;
    save.n_s_s_assemblage_user_end = n_user_end;
    break;
  }
  check_line ("End of save", FALSE, TRUE, TRUE, TRUE);
  /* empty, eof, keyword, print */
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
read_selected_output (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Read data for to output to flat file
 */
  int value;
  int return_value, opt, opt_save;
  char *next_char;
  const char *opt_list[] = {
    "file",			/* 0 */
    "totals",			/* 1 */
    "molalities",		/* 2 */
    "activities",		/* 3 */
    "pure_phases",		/* 4 */
    "si",			/* 5 */
    "saturation_indices",	/* 6 */
    "gases",			/* 7 */
    "equilibrium_phases",	/* 8 */
    "equilibria",		/* 9 */
    "equilibrium",		/* 10 */
    "pure",			/* 11 */
    "inverse",			/* 12 */
    "kinetic_reactants",	/* 13 */
    "kinetics",			/* 14 */
    "solid_solutions",		/* 15 */
    "inverse_modeling",		/* 16 */
    "reset",			/* 17 */
    "simulation",		/* 18 */
    "sim",			/* 19 */
    "state",			/* 20 */
    "solution",			/* 21 */
    "soln",			/* 22 */
    "distance",			/* 23 */
    "dist",			/* 24 */
    "time",			/* 25 */
    "step",			/* 26 */
    "reaction",			/* 27 */
    "rxn",			/* 28 */
    "temperature",		/* 29 */
    "temp",			/* 30 */
    "ph",			/* 31 */
    "pe",			/* 32 */
    "alkalinity",		/* 33 */
    "alk",			/* 34 */
    "ionic_strength",		/* 35 */
    "mu",			/* 36 */
    "water",			/* 37 */
    "high_precision",		/* 38 */
    "user_punch",		/* 39 */
    "mol",			/* 40 */
    "kin",			/* 41 */
    "charge_balance",		/* 42 */
    "percent_error",		/* 43 */
    "selected_out",		/* 44 */
    "selected_output",		/* 45 */
    "isotopes",			/* 46 */
    "calculate_values"		/* 47 */
  };
  int count_opt_list = 48;

  int i, l;
  char file_name[MAX_LENGTH], token[MAX_LENGTH];

  punch.in = TRUE;
  punch.new_def = TRUE;
  punch.count_totals = 0;
  punch.count_molalities = 0;
  punch.count_activities = 0;
  punch.count_pure_phases = 0;
  punch.count_si = 0;
  punch.count_gases = 0;
  punch.count_kinetics = 0;
  punch.count_s_s = 0;
/*
 *   Read eqn from file and call parser
 */
  opt_save = OPTION_ERROR;
  return_value = UNKNOWN;
  for (;;)
  {
    opt = get_option (opt_list, count_opt_list, &next_char);
    if (opt == OPTION_DEFAULT)
    {
      opt = opt_save;
    }
    opt_save = opt;
    switch (opt)
    {
    case OPTION_EOF:		/* end of file */
      return_value = EOF;
      break;
    case OPTION_KEYWORD:	/* keyword */
      return_value = KEYWORD;
      break;
    case OPTION_DEFAULT:
    case OPTION_ERROR:
      input_error++;
      error_msg ("Unknown input in SELECTED_OUTPUT keyword.", CONTINUE);
      error_msg (line_save, CONTINUE);
      break;
    case 0:			/* file name */
      /* copy_token(file_name, &next_char, &l); */
      if (string_trim (next_char) != EMPTY) 
	  {
		strcpy (file_name, next_char);
		have_punch_name = TRUE;
		if (output_open (OUTPUT_PUNCH, file_name) != OK)
		{
			sprintf (error_string, "Can't open file, %s.", file_name);
			input_error++;
			error_msg (error_string, CONTINUE);
		}
	  } 
      opt_save = OPTION_ERROR;
      break;
    case 1:			/* totals */
      while ((i = copy_token (token, &next_char, &l)) != EMPTY)
      {
	if (i != UPPER && token[0] != '[')
	{
	  sprintf (error_string, "Expected element name to"
		   " begin with upper case letter.");
	  warning_msg (error_string);
	}
	else
	{
	  punch.count_totals++;
	  punch.totals =
	    (struct name_master *) PHRQ_realloc (punch.totals,
						 (size_t) punch.count_totals *
						 sizeof (struct name_master));
	  if (punch.totals == NULL)
	    malloc_error ();
	  punch.totals[punch.count_totals - 1].name = string_hsave (token);
	}
      }
      break;
    case 2:			/* molalities */
    case 40:			/* mol */
      while ((i = copy_token (token, &next_char, &l)) != EMPTY)
      {
	if (i != UPPER && token[0] != '(' && (token[0] != '['))
	{
	  sprintf (error_string, "Expected species name to"
		   " begin with upper case letter.");
	  warning_msg (error_string);
	}
	else
	{
	  punch.count_molalities++;
	  punch.molalities =
	    (struct name_species *) PHRQ_realloc (punch.molalities,
						  (size_t) punch.
						  count_molalities *
						  sizeof (struct
							  name_species));
	  if (punch.molalities == NULL)
	    malloc_error ();
	  punch.molalities[punch.count_molalities - 1].name =
	    string_hsave (token);
	}
      }
      break;
    case 3:			/* activities */
      while ((i = copy_token (token, &next_char, &l)) != EMPTY)
      {
	if (i != UPPER && token[0] != '(' && (token[0] != '['))
	{
	  sprintf (error_string, "Expected species name to"
		   " begin with upper case letter.");
	  warning_msg (error_string);
	}
	else
	{
	  punch.count_activities++;
	  punch.activities =
	    (struct name_species *) PHRQ_realloc (punch.activities,
						  (size_t) punch.
						  count_activities *
						  sizeof (struct
							  name_species));
	  if (punch.activities == NULL)
	    malloc_error ();
	  punch.activities[punch.count_activities - 1].name =
	    string_hsave (token);
	}
      }
      break;
    case 4:			/* pure_phases */
    case 8:			/* equilibrium_phases */
    case 9:			/* equilibria */
    case 10:			/* equilibrium */
    case 11:			/* pure */
      while ((i = copy_token (token, &next_char, &l)) != EMPTY)
      {
	punch.count_pure_phases++;
	punch.pure_phases =
	  (struct name_phase *) PHRQ_realloc (punch.pure_phases,
					      (size_t) punch.
					      count_pure_phases *
					      sizeof (struct name_phase));
	if (punch.pure_phases == NULL)
	  malloc_error ();
	punch.pure_phases[punch.count_pure_phases - 1].name =
	  string_hsave (token);
      }
      break;
    case 5:			/* si */
    case 6:			/* saturation_index */
      while ((i = copy_token (token, &next_char, &l)) != EMPTY)
      {
	punch.count_si++;
	punch.si =
	  (struct name_phase *) PHRQ_realloc (punch.si,
					      (size_t) punch.count_si *
					      sizeof (struct name_phase));
	if (punch.si == NULL)
	  malloc_error ();
	punch.si[punch.count_si - 1].name = string_hsave (token);
      }
      break;
    case 7:			/* gases */
      while ((i = copy_token (token, &next_char, &l)) != EMPTY)
      {
	punch.count_gases++;
	punch.gases =
	  (struct name_phase *) PHRQ_realloc (punch.gases,
					      (size_t) punch.count_gases *
					      sizeof (struct name_phase));
	if (punch.gases == NULL)
	  malloc_error ();
	punch.gases[punch.count_gases - 1].name = string_hsave (token);
      }
      break;
    case 12:			/* inverse */
    case 16:			/* inverse_modeling */
      punch.inverse = get_true_false (next_char, TRUE);
      break;
    case 13:			/* kinetic_reactants */
    case 14:			/* kinetics */
    case 41:			/* kin */
      while ((i = copy_token (token, &next_char, &l)) != EMPTY)
      {
	punch.count_kinetics++;
	punch.kinetics =
	  (struct name_phase *) PHRQ_realloc (punch.kinetics,
					      (size_t) punch.count_kinetics *
					      sizeof (struct name_phase));
	if (punch.kinetics == NULL)
	  malloc_error ();
	punch.kinetics[punch.count_kinetics - 1].name = string_hsave (token);
      }
      break;
    case 15:			/* solid_solutions */
      while ((i = copy_token (token, &next_char, &l)) != EMPTY)
      {
	punch.count_s_s++;
	punch.s_s =
	  (struct name_phase *) PHRQ_realloc (punch.s_s,
					      (size_t) punch.count_s_s *
					      sizeof (struct name_phase));
	if (punch.s_s == NULL)
	  malloc_error ();
	punch.s_s[punch.count_s_s - 1].name = string_hsave (token);
      }
      break;
    case 46:			/* isotopes */
      while ((i = copy_token (token, &next_char, &l)) != EMPTY)
      {
	if (i != UPPER && token[0] != '[')
	{
	  sprintf (error_string, "Expected element name to"
		   " begin with upper case letter.");
	  warning_msg (error_string);
	}
	else
	{
	  punch.count_isotopes++;
	  punch.isotopes =
	    (struct name_master *) PHRQ_realloc (punch.isotopes,
						 (size_t) punch.
						 count_isotopes *
						 sizeof (struct name_master));
	  if (punch.isotopes == NULL)
	    malloc_error ();
	  punch.isotopes[punch.count_isotopes - 1].name =
	    string_hsave (token);
	}
      }
      break;
    case 47:			/* calculate_values */
      while ((i = copy_token (token, &next_char, &l)) != EMPTY)
      {
	punch.count_calculate_values++;
	punch.calculate_values =
	  (struct name_master *) PHRQ_realloc (punch.calculate_values,
					       (size_t) punch.
					       count_calculate_values *
					       sizeof (struct name_master));
	if (punch.calculate_values == NULL)
	  malloc_error ();
	punch.calculate_values[punch.count_calculate_values - 1].name =
	  string_hsave (token);
      }
      break;
    case 17:			/* reset */
      value = get_true_false (next_char, TRUE);
      punch.sim = value;
      punch.state = value;
      punch.soln = value;
      punch.dist = value;
      punch.time = value;
      punch.step = value;
      punch.rxn = value;
      punch.temp = value;
      punch.ph = value;
      punch.pe = value;
      punch.alk = value;
      punch.mu = value;
      punch.water = value;
      punch.charge_balance = value;
      punch.percent_error = value;
      break;
    case 18:			/* simulation */
    case 19:			/* sim */
      punch.sim = get_true_false (next_char, TRUE);
      break;
    case 20:			/* state */
      punch.state = get_true_false (next_char, TRUE);
      break;
    case 21:			/* solution */
    case 22:			/* soln */
      punch.soln = get_true_false (next_char, TRUE);
      break;
    case 23:			/* distance */
    case 24:			/* dist */
      punch.dist = get_true_false (next_char, TRUE);
      break;
    case 25:			/* time */
      punch.time = get_true_false (next_char, TRUE);
      break;
    case 26:			/* step */
      punch.step = get_true_false (next_char, TRUE);
      break;
    case 27:			/* reaction */
    case 28:			/* rxn */
      punch.rxn = get_true_false (next_char, TRUE);
      break;
    case 29:			/* temperature */
    case 30:			/* temp */
      punch.temp = get_true_false (next_char, TRUE);
      break;
    case 31:			/* ph */
      punch.ph = get_true_false (next_char, TRUE);
      break;
    case 32:			/* pe */
      punch.pe = get_true_false (next_char, TRUE);
      break;
    case 33:			/* alkalinity */
    case 34:			/* alk */
      punch.alk = get_true_false (next_char, TRUE);
      break;
    case 35:			/* ionic strength */
    case 36:			/* mu */
      punch.mu = get_true_false (next_char, TRUE);
      break;
    case 37:			/* water */
      punch.water = get_true_false (next_char, TRUE);
      break;
    case 38:			/* high_precision */
      punch.high_precision = get_true_false (next_char, TRUE);
      if (punch.high_precision == TRUE && convergence_tolerance > 1e-12)
      {
	convergence_tolerance = 1e-12;
      }
      break;
    case 39:			/* user_punch */
      punch.user_punch = get_true_false (next_char, TRUE);
      break;
    case 42:			/* charge_balance */
      punch.charge_balance = get_true_false (next_char, TRUE);
      break;
    case 43:			/* percent_error */
      punch.percent_error = get_true_false (next_char, TRUE);
      break;
    case 44:			/* selected_out */
    case 45:			/* selected_output */
      pr.punch = get_true_false (next_char, TRUE);
      break;
    }
    if (return_value == EOF || return_value == KEYWORD)
      break;
  }
  if (!have_punch_name)
  {
    if (output_open (OUTPUT_PUNCH, "selected.out") != OK)
    {
      sprintf (error_string, "Can't open file, %s.", "selected.out");
      input_error++;
      error_msg (error_string, CONTINUE);
    }
  }

  return (return_value);
}

/* ---------------------------------------------------------------------- */
int
read_solution (void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads solution data
 *
 *      Arguments:
 *         none
 *
 *      Returns:
 *         KEYWORD if keyword encountered, input_error may be incremented if
 *                    a keyword is encountered in an unexpected position
 *         EOF     if eof encountered while reading mass balance concentrations
 *         ERROR   if error occurred reading data
 *
 */
  int i, j, n, l;
  int n_user, n_user_end;
  int default_pe, alk;
  int count_isotopes;
  int max_mass_balance, count_mass_balance;
  char *ptr, *ptr1;
  char *description;
  char token[MAX_LENGTH], token1[MAX_LENGTH];

  int return_value, opt;
  char *next_char;
  const char *opt_list[] = {
    "temp",			/* 0 */
    "temperature",		/* 1 */
    "dens",			/* 2 */
    "density",			/* 3 */
    "units",			/* 4 */
    "redox",			/* 5 */
    "ph",			/* 6 */
    "pe",			/* 7 */
    "unit",			/* 8 */
    "isotope",			/* 9 */
    "water"			/* 10 */
  };
  int count_opt_list = 11;
/*
 *   Read solution number and description
 */
  ptr = line;
  read_number_description (ptr, &n_user, &n_user_end, &description);
/*
 *   Malloc space for solution data
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
  solution[n] = solution_alloc ();

  solution[n]->n_user = n_user;
  solution[n]->n_user_end = n_user_end;
  if (use.solution_in == FALSE)
  {
    use.solution_in = TRUE;
    use.n_solution_user = n_user;
  }
  max_mass_balance = MAX_MASS_BALANCE;
/*
 *   Set default ph, temp, density, pe, units
 */
  solution[n]->description = description;
  solution[n]->tc = 25.0;
  solution[n]->ph = 7.0;
  solution[n]->density = 1.0;
  solution[n]->solution_pe = 4.0;
  solution[n]->mass_water = 1.0;
  solution[n]->ah2o = 1.0;
  solution[n]->mu = 1e-7;
  solution[n]->cb = 0.0;
  default_pe = 0;
  solution[n]->units = string_hsave ("mMol/kgw");
  solution[n]->totals[0].description = NULL;
  count_mass_balance = 0;
  count_isotopes = 0;
/*
 *   Read concentration data
 */
  return_value = UNKNOWN;
  for (;;)
  {
    opt = get_option (opt_list, count_opt_list, &next_char);
    if (opt == OPTION_DEFAULT)
    {
      ptr = next_char;
      if (copy_token (token, &ptr, &l) == DIGIT)
      {
	opt = 9;
      }
    }
    switch (opt)
    {
    case OPTION_EOF:		/* end of file */
      return_value = EOF;
      break;
    case OPTION_KEYWORD:	/* keyword */
      return_value = KEYWORD;
      break;
    case OPTION_ERROR:
      input_error++;
      error_msg ("Unknown input in SOLUTION keyword.", CONTINUE);
      error_msg (line_save, CONTINUE);
      break;
    case 0:			/* temperature */
    case 1:
      if (sscanf (next_char, SCANFORMAT, &(solution[n]->tc)) != 1)
      {
	solution[n]->tc = 25;
      }
      break;
    case 2:			/* density */
    case 3:
      sscanf (next_char, SCANFORMAT, &(solution[n]->density));
      break;
    case 4:			/* units */
    case 8:			/* unit */
      if (copy_token (token, &next_char, &l) == EMPTY)
	break;
      if (check_units (token, FALSE, FALSE, solution[n]->units, TRUE) == OK)
      {
	solution[n]->units = string_hsave (token);
      }
      else
      {
	input_error++;
      }
      break;
    case 5:			/* redox */
      if (copy_token (token, &next_char, &l) == EMPTY)
	break;
      if (parse_couple (token) == OK)
      {
	default_pe = pe_data_store (&(solution[n]->pe), token);
      }
      else
      {
	input_error++;
      }
      break;
    case 6:			/* ph */
      if (read_conc (n, count_mass_balance, line) == ERROR)
      {
	input_error++;
	break;
      }
      solution[n]->ph = solution[n]->totals[count_mass_balance].input_conc;
      if (solution[n]->totals[count_mass_balance].equation_name == NULL)
      {
	break;
      }
      solution[n]->totals[count_mass_balance].description =
	string_hsave ("H(1)");
      count_mass_balance++;
      break;
    case 7:			/* pe */
      if (read_conc (n, count_mass_balance, line) == ERROR)
      {
	input_error++;
	break;
      }
      solution[n]->solution_pe =
	solution[n]->totals[count_mass_balance].input_conc;
      if (solution[n]->totals[count_mass_balance].equation_name == NULL)
      {
	break;
      }
      solution[n]->totals[count_mass_balance].description =
	string_hsave ("E");
      count_mass_balance++;
      break;
    case 9:			/* isotope */
      if (copy_token (token, &next_char, &l) != DIGIT)
      {
	input_error++;
	sprintf (error_string, "Expected isotope name to"
		 " begin with an isotopic number.");
	error_msg (error_string, CONTINUE);
	continue;
      }
      solution[n]->isotopes =
	(struct isotope *) PHRQ_realloc (solution[n]->isotopes,
					 (size_t) (count_isotopes +
						   1) *
					 sizeof (struct isotope));
      if (solution[n]->isotopes == NULL)
	malloc_error ();

      /* read and save element name */
      ptr1 = token;
      get_num (&ptr1,
	       &(solution[n]->isotopes[count_isotopes].isotope_number));
      if (ptr1[0] == '\0' || isupper ((int) ptr1[0]) == FALSE)
      {
	error_msg ("Expecting element name.", CONTINUE);
	error_msg (line_save, CONTINUE);
	input_error++;
	return (ERROR);
      }
      solution[n]->isotopes[count_isotopes].elt_name = string_hsave (ptr1);

      /* read and store isotope ratio */
      if (copy_token (token, &next_char, &l) != DIGIT)
      {
	input_error++;
	sprintf (error_string, "Expected numeric value for isotope ratio.");
	error_msg (error_string, CONTINUE);
	continue;
      }
      sscanf (token, SCANFORMAT,
	      &(solution[n]->isotopes[count_isotopes].ratio));

      solution[n]->isotopes[count_isotopes].ratio_uncertainty = NAN;

      /* read and store isotope ratio uncertainty */
      if ((j = copy_token (token, &next_char, &l)) != EMPTY)
      {
	if (j != DIGIT)
	{
	  input_error++;
	  sprintf (error_string,
		   "Expected numeric value for uncertainty in isotope ratio.");
	  error_msg (error_string, CONTINUE);
	  continue;
	}
	sscanf (token, SCANFORMAT,
		&(solution[n]->isotopes[count_isotopes].ratio_uncertainty));
      }
      count_isotopes++;
      break;
    case 10:			/* water */
      j = copy_token (token, &next_char, &l);
      if (j == EMPTY)
      {
	solution[n]->mass_water = 1.0;
      }
      else if (j != DIGIT)
      {
	input_error++;
	sprintf (error_string,
		 "Expected numeric value for mass of water in solution.");
	error_msg (error_string, CONTINUE);
      }
      else
      {
	sscanf (token, SCANFORMAT, &dummy);
	solution[n]->mass_water = (LDBLE) dummy;
      }
      break;
    case OPTION_DEFAULT:
/*
 *   Read concentration
 */
      if (read_conc (n, count_mass_balance, line) == ERROR)
      {
	input_error++;
	break;
      }
      count_mass_balance++;
      break;
    }
    if (count_mass_balance + 1 >= max_mass_balance)
    {
      space ((void **) ((void *) &(solution[n]->totals)), count_mass_balance + 1,
	     &max_mass_balance, sizeof (struct conc));
    }
    if (return_value == EOF || return_value == KEYWORD)
      break;
  }
#ifdef SKIP
/*
 *   Sort totals by description
 */
  qsort (solution[n]->totals,
	 (size_t) count_mass_balance,
	 (size_t) sizeof (struct conc), conc_compare);
#endif
/*
 *   fix up default units and default pe
 */
  for (i = 0; i < count_mass_balance; i++)
  {
    strcpy (token, solution[n]->totals[i].description);
    str_tolower (token);
    if (solution[n]->totals[i].units == NULL)
    {
      solution[n]->totals[i].units = solution[n]->units;
    }
    else
    {
      alk = FALSE;
      if (strstr (token, "alk") == token)
	alk = TRUE;
      strcpy (token1, solution[n]->totals[i].units);
      if (check_units (token1, alk, TRUE, solution[n]->units, TRUE) == ERROR)
      {
	input_error++;
      }
      else
      {
	solution[n]->totals[i].units = string_hsave (token1);
      }
    }
    if (solution[n]->totals[i].n_pe < 0)
    {
      solution[n]->totals[i].n_pe = default_pe;
    }
  }
  solution[n]->default_pe = default_pe;
/*
 *   Mark end of solution
 */
  solution[n]->totals[count_mass_balance].description = NULL;
  solution[n]->count_isotopes = count_isotopes;
#ifdef SKIP
  if (count_isotopes > 0)
  {
    qsort (solution[n]->isotopes,
	   (size_t) count_isotopes,
	   (size_t) sizeof (struct isotope), isotope_compare);
  }
  else
  {
    solution[n]->isotopes = free_check_null (solution[n]->isotopes);
  }
#endif
  return (return_value);
}

/* ---------------------------------------------------------------------- */
int
read_species (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Read data for aqueous species, parse equations
 */
  int i;
  int association;
  struct species *s_ptr;
  struct elt_list *next_elt;
  char *ptr, token[MAX_LENGTH];

  int return_value, opt, opt_save;
  char *next_char;
  const char *opt_list[] = {
    "no_check",			/* 0 */
    "check",			/* 1 */
    "gamma",			/* 2 */
    "mb",			/* 3 */
    "mass_balance",		/* 4 */
    "log_k",			/* 5 */
    "logk",			/* 6 */
    "delta_h",			/* 7 */
    "deltah",			/* 8 */
    "analytical_expression",	/* 9 */
    "a_e",			/* 10 */
    "ae",			/* 11 */
    "mole_balance",		/* 12 */
    "llnl_gamma",		/* 13 */
    "co2_llnl_gamma",		/* 14 */
    "activity_water",		/* 15 */
    "add_logk",			/* 16 */
    "add_log_k",		/* 17 */
    "add_constant",		/* 18 */
    "dw",			/* 19 */
  };
  int count_opt_list = 20;

  association = TRUE;
  s_ptr = NULL;
/*
 *   Read eqn from file and call parser
 */
  opt_save = OPTION_DEFAULT;
  return_value = UNKNOWN;
  for (;;)
  {
    opt = get_option (opt_list, count_opt_list, &next_char);
    if (opt == OPTION_DEFAULT)
    {
      opt = opt_save;
    }
    switch (opt)
    {
    case OPTION_EOF:		/* end of file */
      return_value = EOF;
      break;
    case OPTION_KEYWORD:	/* keyword */
      return_value = KEYWORD;
      break;
    case OPTION_ERROR:
      input_error++;
      error_msg ("Unknown input in SPECIES keyword.", CONTINUE);
      error_msg (line_save, CONTINUE);
      break;
    case 0:			/* no_check */
      if (s_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      s_ptr->check_equation = FALSE;
      break;
    case 1:			/* check */
      if (s_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      s_ptr->check_equation = TRUE;
      break;
    case 2:			/* gamma data */
      if (s_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      s_ptr->gflag = 2;		/* Wateq D-H */
      i = sscanf (next_char, SCANFORMAT SCANFORMAT, &s_ptr->dha, &s_ptr->dhb);
      if (i < 2)
      {
	sprintf (error_string, "Expecting 2 activity-"
		 "coefficient parameters, a and b.");
	warning_msg (error_string);
      }
      opt_save = OPTION_DEFAULT;
      break;
    case 3:			/* mb */
    case 4:			/* mass_balance */
    case 12:			/* mole_balance */
      if (s_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      count_elts = 0;
      paren_count = 0;
      copy_token (token, &next_char, &i);
      s_ptr->mole_balance = string_hsave (token);
      ptr = token;
      s_ptr->next_secondary =
	(struct elt_list *) free_check_null (s_ptr->next_secondary);
      get_secondary_in_species (&ptr, 1.0);
      s_ptr->next_secondary = elt_list_save ();
/* debug
			for (i = 0; i < count_elts; i++) {
				output_msg(OUTPUT_MESSAGE,"%s\t%f\n", elt_list[i].elt->name,
					elt_list[i].coef);
			}
 */
      opt_save = OPTION_DEFAULT;
      break;
    case 5:			/* log_k */
    case 6:			/* logk */
      if (s_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      read_log_k_only (next_char, &s_ptr->logk[0]);
      opt_save = OPTION_DEFAULT;
      break;
    case 7:			/* delta_h */
    case 8:			/* deltah */
      if (s_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      read_delta_h_only (next_char, &s_ptr->logk[1], &s_ptr->original_units);
      opt_save = OPTION_DEFAULT;
      break;
    case 9:			/* analytical_expression */
    case 10:			/* a_e */
    case 11:			/* ae */
      if (s_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      read_analytical_expression_only (next_char, &(s_ptr->logk[2]));
      opt_save = OPTION_DEFAULT;
      break;
    case 13:			/* llnl_gamma */
      if (s_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      s_ptr->gflag = 7;		/* llnl D-H */
      i = sscanf (next_char, SCANFORMAT, &s_ptr->dha);
      if (i < 1)
      {
	sprintf (error_string,
		 "Expecting activity-coefficient parameter, a.");
	warning_msg (error_string);
      }
      opt_save = OPTION_DEFAULT;
      break;
    case 14:			/* co2_llnl_gamma */
      if (s_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      s_ptr->gflag = 8;		/* llnl CO2 D-H */
      opt_save = OPTION_DEFAULT;
      break;
    case 15:			/* activity water */
      if (s_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      s_ptr->gflag = 9;		/* activity_water/55.5 for HDO, D2O, H2[O18], etc */
      opt_save = OPTION_DEFAULT;
      break;
    case 16:			/* add_logk */
    case 17:			/* add_log_k */
      if (s_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      if (s_ptr->count_add_logk == 0)
      {
	s_ptr->add_logk =
	  (struct name_coef *) PHRQ_malloc (sizeof (struct name_coef));
	if (s_ptr->add_logk == NULL)
	  malloc_error ();
      }
      else
      {
	s_ptr->add_logk =
	  (struct name_coef *) PHRQ_realloc (s_ptr->add_logk,
					     (size_t) ((s_ptr->
							count_add_logk +
							1) *
						       sizeof (struct
							       name_coef)));
	if (s_ptr->add_logk == NULL)
	  malloc_error ();
      }
      /* read name */
      if (copy_token (token, &next_char, &i) == EMPTY)
      {
	input_error++;
	sprintf (error_string, "Expected the name of a NAMED_EXPRESSION.");
	error_msg (error_string, CONTINUE);
	break;
      }
      s_ptr->add_logk[s_ptr->count_add_logk].name = string_hsave (token);
      /* read coef */
      i =
	sscanf (next_char, SCANFORMAT,
		&s_ptr->add_logk[s_ptr->count_add_logk].coef);
      if (i <= 0)
      {
	s_ptr->add_logk[s_ptr->count_add_logk].coef = 1;
      }
      s_ptr->count_add_logk++;
      opt_save = OPTION_DEFAULT;
      break;
    case 18:			/* add_constant */
      if (s_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      if (s_ptr->count_add_logk == 0)
      {
	s_ptr->add_logk =
	  (struct name_coef *) PHRQ_malloc (sizeof (struct name_coef));
	if (s_ptr->add_logk == NULL)
	  malloc_error ();
      }
      else
      {
	s_ptr->add_logk =
	  (struct name_coef *) PHRQ_realloc (s_ptr->add_logk,
					     (size_t) ((s_ptr->
							count_add_logk +
							1) *
						       sizeof (struct
							       name_coef)));
	if (s_ptr->add_logk == NULL)
	  malloc_error ();
      }
      i =
	sscanf (next_char, SCANFORMAT,
		&s_ptr->add_logk[s_ptr->count_add_logk].coef);
      if (i <= 0)
      {
	input_error++;
	sprintf (error_string,
		 "Expected the constant to add for log_K definition.");
	error_msg (error_string, CONTINUE);
	break;
      }
      /* set name */
      s_ptr->add_logk[s_ptr->count_add_logk].name =
	string_hsave ("XconstantX");
      /* read coef */
      s_ptr->count_add_logk++;
      opt_save = OPTION_DEFAULT;
      break;
    case 19:			/* tracer diffusion coefficient */
      if (s_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      i = sscanf (next_char, SCANFORMAT, &s_ptr->dw);
      opt_save = OPTION_DEFAULT;
      break;

    case OPTION_DEFAULT:
/*
 *   Get space for species information and parse equation
 */
      s_ptr = NULL;
      if (parse_eq (line, &next_elt, association) == ERROR)
      {
	parse_error++;
	error_msg ("Parsing equation.", CONTINUE);
	error_msg (line_save, CONTINUE);
	break;
      }
/*
 *   Get pointer to each species in the reaction, store new species if necessary
 */
      trxn.token[0].s = s_store (trxn.token[0].name, trxn.token[0].z, TRUE);
      for (i = 1; i < count_trxn; i++)
      {
	trxn.token[i].s =
	  s_store (trxn.token[i].name, trxn.token[i].z, FALSE);
      }
/*
 *   Save element list and carbon, hydrogen, and oxygen in species
 */
      trxn.token[0].s->next_elt = next_elt;
      trxn.token[0].s->next_secondary = NULL;
      for (; next_elt->elt != NULL; next_elt++)
      {
	if (strcmp (next_elt->elt->name, "C") == 0)
	{
	  trxn.token[0].s->carbon = next_elt->coef;
	}
	if (strcmp (next_elt->elt->name, "H") == 0)
	{
	  trxn.token[0].s->h = next_elt->coef;
	}
	if (strcmp (next_elt->elt->name, "O") == 0)
	{
	  trxn.token[0].s->o = next_elt->coef;
	}
      }
/*
 *   Malloc space for species reaction
 */
      trxn.token[0].s->rxn = rxn_alloc (count_trxn + 1);
/*
 *   Copy reaction to reaction for species
 */
      trxn_copy (trxn.token[0].s->rxn);
      s_ptr = trxn.token[0].s;
/*
 *   Default gamma data
 */
      s_ptr->dha = 0.0;
      s_ptr->dhb = 0.0;
      if (equal (s_ptr->z, 0.0, TOL) == TRUE)
      {
	s_ptr->gflag = 0;	/* Uncharged */
	s_ptr->dhb = 0.1;
      }
      else
      {
	s_ptr->gflag = 1;	/* Davies */
      }
/*
 *   Set type for species
 */
      if (strcmp (trxn.token[0].s->name, "H+") == 0)
      {
	s_hplus = trxn.token[0].s;
	s_hplus->type = HPLUS;
      }
      else if (strcmp (trxn.token[0].s->name, "H3O+") == 0)
      {
	s_h3oplus = trxn.token[0].s;
	s_h3oplus->type = HPLUS;
      }
      else if (strcmp (trxn.token[0].s->name, "e-") == 0)
      {
	s_eminus = trxn.token[0].s;
	s_eminus->type = EMINUS;
	s_eminus->gflag = 3;	/* Always 1 */
      }
      else if (strcmp (trxn.token[0].s->name, "H2O") == 0)
      {
	s_h2o = trxn.token[0].s;
	s_h2o->type = H2O;
	s_h2o->gflag = 3;	/* Always 1 */
      }
      else if (strstr (trxn.token[0].s->name, "(s)") != NULL)
      {
	trxn.token[0].s->type = SOLID;
      }
      else if (strcmp (trxn.token[0].s->name, "H2") == 0)
      {
	trxn.token[0].s->type = AQ;
	s_h2 = trxn.token[0].s;
      }
      else if (strcmp (trxn.token[0].s->name, "O2") == 0)
      {
	trxn.token[0].s->type = AQ;
	s_o2 = trxn.token[0].s;
      }
      else
      {
	trxn.token[0].s->type = AQ;
      }
      opt_save = OPTION_DEFAULT;
      break;
    }
    if (return_value == EOF || return_value == KEYWORD)
      break;
  }
  return (return_value);
}

/* ---------------------------------------------------------------------- */
int
read_use (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Reads solution, mix, irreversible reaction, and pure phases to use
 *   in reaction calculation
 */
  int i, l, n_user, return_value;
  char *ptr;
  char token[MAX_LENGTH], token1[MAX_LENGTH];;
/*
 *   Read "use"
 */
  ptr = line;
  copy_token (token, &ptr, &l);
/*
 *   Read keyword
 */
  copy_token (token, &ptr, &l);
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
    input_error++;
    error_msg
      ("Expecting keyword solution, mix, kinetics, reaction, reaction_temperature, equilibrium_phases, exchange, surface, gas_phase, or solid_solutions.",
       CONTINUE);
    error_msg (line_save, CONTINUE);
    check_line ("End of use", FALSE, TRUE, TRUE, TRUE);
    /* empty, eof, keyword, print */
    return (ERROR);
  }
/*
 *   Read number
 */
  strcpy (token1, token);
  for (;;)
  {
    i = copy_token (token, &ptr, &l);
    if (i == DIGIT)
    {
      sscanf (token, "%d", &n_user);
      if (n_user < 0)
      {
	error_msg ("Number must be a positive integer.", CONTINUE);
	error_msg (line_save, CONTINUE);
	input_error++;
      }
      if (strstr (token, "-") != NULL)
      {
	sprintf (error_string, "USE does not accept a range of numbers, %s.",
		 token);
	warning_msg (error_string);
	sprintf (error_string,
		 "Only %s %d will be used in the batch-reaction calculation.",
		 token1, n_user);
	warning_msg (error_string);
	sprintf (error_string,
		 "NOTE--USE is not needed for ADVECTION and TRANSPORT calculations.");
	warning_msg (error_string);

      }
      break;
    }
    else if (i == EMPTY)
    {
      sprintf (error_string, "No number given, 1 assumed.");
      warning_msg (error_string);
      n_user = 1;
      break;
    }
    else if (token[0] == 'N' || token[0] == 'n')
    {
      n_user = -2;
      break;
    }
  }
  switch (next_keyword)
  {
  case 4:			/* Solution */
    use.n_solution_user = n_user;
    if (n_user >= 0)
    {
      use.solution_in = TRUE;
    }
    else
    {
      use.solution_in = FALSE;
    }
    break;
  case 6:			/* Pure phases */
  case 26:
  case 27:
  case 28:
  case 29:
    use.n_pp_assemblage_user = n_user;
    if (n_user >= 0)
    {
      use.pp_assemblage_in = TRUE;
    }
    else
    {
      use.pp_assemblage_in = FALSE;
    }
    break;
  case 7:			/* Reaction */
    use.n_irrev_user = n_user;
    if (n_user >= 0)
    {
      use.irrev_in = TRUE;
    }
    else
    {
      use.irrev_in = FALSE;
    }
    break;
  case 8:			/* Mix */
    use.n_mix_user = n_user;
    if (n_user >= 0)
    {
      use.mix_in = TRUE;
    }
    else
    {
      use.mix_in = FALSE;
    }
    break;
  case 13:			/* Ex */
    use.n_exchange_user = n_user;
    if (n_user >= 0)
    {
      use.exchange_in = TRUE;
    }
    else
    {
      use.exchange_in = FALSE;
    }
    break;
  case 16:			/* Surface */
    use.n_surface_user = n_user;
    if (n_user >= 0)
    {
      use.surface_in = TRUE;
    }
    else
    {
      use.surface_in = FALSE;
    }
    break;
  case 17:			/* Temperature */
    use.n_temperature_user = n_user;
    if (n_user >= 0)
    {
      use.temperature_in = TRUE;
    }
    else
    {
      use.temperature_in = FALSE;
    }
    break;
  case 19:			/* Gas */
    use.n_gas_phase_user = n_user;
    if (n_user >= 0)
    {
      use.gas_phase_in = TRUE;
    }
    else
    {
      use.gas_phase_in = FALSE;
    }
    break;
  case 33:			/* Kinetics */
    use.n_kinetics_user = n_user;
    if (n_user >= 0)
    {
      use.kinetics_in = TRUE;
    }
    else
    {
      use.kinetics_in = FALSE;
    }
    break;
  case 40:			/* solid_solutions */
  case 41:			/* solid_solution */
    use.n_s_s_assemblage_user = n_user;
    if (n_user >= 0)
    {
      use.s_s_assemblage_in = TRUE;
    }
    else
    {
      use.s_s_assemblage_in = FALSE;
    }
    break;
  }
  return_value = check_line ("End of use", FALSE, TRUE, TRUE, TRUE);
  /* empty, eof, keyword, print */
  return (return_value);
}

/* ---------------------------------------------------------------------- */
int
read_surface_species (void)
/* ---------------------------------------------------------------------- */
{
  /*
   *   Read data for surface species, parse equations
   */
  int i, j;
  int association;
  char token[MAX_LENGTH];
  char *ptr;
  LDBLE offset;

  struct species *s_ptr;
  struct elt_list *next_elt;
  struct rxn_token *token_ptr;

  int return_value, opt, opt_save;
  char *next_char;
  const char *opt_list[] = {
    "no_check",			/* 0 */
    "check",			/* 1 */
    "mb",			/* 2 */
    "mass_balance",		/* 3 */
    "log_k",			/* 4 */
    "logk",			/* 5 */
    "delta_h",			/* 6 */
    "deltah",			/* 7 */
    "analytical_expression",	/* 8 */
    "a_e",			/* 9 */
    "ae",			/* 10 */
    "mole_balance",		/* 11 */
    "offset",			/* 12 */
    "add_logk",			/* 13 */
    "add_log_k",		/* 14 */
    "add_constant",		/* 15 */
    "cd_music",			/* 16 */
    "music"			/* 17 */
  };
  int count_opt_list = 18;

  association = TRUE;
  /*
   *   Read eqn from file and call parser
   */
  opt_save = OPTION_DEFAULT;
  return_value = UNKNOWN;
  s_ptr = NULL;
  for (;;)
  {
    opt = get_option (opt_list, count_opt_list, &next_char);
    if (opt == OPTION_DEFAULT)
    {
      opt = opt_save;
    }
    switch (opt)
    {
    case OPTION_EOF:		/* end of file */
      return_value = EOF;
      break;
    case OPTION_KEYWORD:	/* keyword */
      return_value = KEYWORD;
      break;
    case OPTION_ERROR:
      input_error++;
      error_msg ("Unknown input in SURFACE_SPECIES keyword.", CONTINUE);
      error_msg (line_save, CONTINUE);
      break;
    case 0:			/* no_check */
      if (s_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      s_ptr->check_equation = FALSE;
      break;
    case 1:			/* check */
      if (s_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      s_ptr->check_equation = TRUE;
      break;
    case 2:			/* mb */
    case 3:			/* mass_balance */
    case 11:			/* mole_balance */
      if (s_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      count_elts = 0;
      paren_count = 0;
      copy_token (token, &next_char, &i);
      s_ptr->mole_balance = string_hsave (token);
      ptr = token;
      s_ptr->next_secondary =
	(struct elt_list *) free_check_null (s_ptr->next_secondary);
      get_secondary_in_species (&ptr, 1.0);
      s_ptr->next_secondary = elt_list_save ();
      /* debug
         for (i = 0; i < count_elts; i++) {
         output_msg(OUTPUT_MESSAGE,"%s\t%f\n", elt_list[i].elt->name,
         elt_list[i].coef);
         }
       */
      opt_save = OPTION_DEFAULT;
      break;
    case 4:			/* log_k */
    case 5:			/* logk */
      if (s_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      read_log_k_only (next_char, &s_ptr->logk[0]);
      opt_save = OPTION_DEFAULT;
      break;
    case 6:			/* delta_h */
    case 7:			/* deltah */
      if (s_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      read_delta_h_only (next_char, &s_ptr->logk[1], &s_ptr->original_units);
      opt_save = OPTION_DEFAULT;
      break;
    case 8:			/* analytical_expression */
    case 9:			/* a_e */
    case 10:			/* ae */
      if (s_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      read_analytical_expression_only (next_char, &(s_ptr->logk[2]));
      opt_save = OPTION_DEFAULT;
      break;
    case 12:			/* offset */
      if (s_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      if (sscanf (next_char, SCANFORMAT, &offset) != 1)
      {
	error_msg ("No offset for log K given", STOP);
      }
      s_ptr->logk[0] += offset;
      opt_save = OPTION_DEFAULT;
      break;
    case 13:			/* add_logk */
    case 14:			/* add_log_k */
      if (s_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      if (s_ptr->count_add_logk == 0)
      {
	s_ptr->add_logk =
	  (struct name_coef *) PHRQ_malloc (sizeof (struct name_coef));
	if (s_ptr->add_logk == NULL)
	  malloc_error ();
      }
      else
      {
	s_ptr->add_logk =
	  (struct name_coef *) PHRQ_realloc (s_ptr->add_logk,
					     (size_t) ((s_ptr->
							count_add_logk +
							1) *
						       sizeof (struct
							       name_coef)));
	if (s_ptr->add_logk == NULL)
	  malloc_error ();
      }
      /* read name */
      if (copy_token (token, &next_char, &i) == EMPTY)
      {
	input_error++;
	sprintf (error_string, "Expected the name of a NAMED_EXPRESSION.");
	error_msg (error_string, CONTINUE);
	break;
      }
      s_ptr->add_logk[s_ptr->count_add_logk].name = string_hsave (token);
      /* read coef */
      i =
	sscanf (next_char, SCANFORMAT,
		&s_ptr->add_logk[s_ptr->count_add_logk].coef);
      if (i <= 0)
      {
	s_ptr->add_logk[s_ptr->count_add_logk].coef = 1;
      }
      s_ptr->count_add_logk++;
      opt_save = OPTION_DEFAULT;
      break;
    case 15:			/* add_constant */
      if (s_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      if (s_ptr->count_add_logk == 0)
      {
	s_ptr->add_logk =
	  (struct name_coef *) PHRQ_malloc (sizeof (struct name_coef));
	if (s_ptr->add_logk == NULL)
	  malloc_error ();
      }
      else
      {
	s_ptr->add_logk =
	  (struct name_coef *) PHRQ_realloc (s_ptr->add_logk,
					     (size_t) ((s_ptr->
							count_add_logk +
							1) *
						       sizeof (struct
							       name_coef)));
	if (s_ptr->add_logk == NULL)
	  malloc_error ();
      }
      i =
	sscanf (next_char, SCANFORMAT,
		&s_ptr->add_logk[s_ptr->count_add_logk].coef);
      if (i <= 0)
      {
	input_error++;
	sprintf (error_string,
		 "Expected the constant to add for log_K definition.");
	error_msg (error_string, CONTINUE);
	break;
      }
      /* set name */
      s_ptr->add_logk[s_ptr->count_add_logk].name =
	string_hsave ("XconstantX");
      /* read coef */
      s_ptr->count_add_logk++;
      opt_save = OPTION_DEFAULT;
      break;
    case 16:			/* cd_music */
    case 17:			/* music */
      if (s_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      for (j = 0; j < 5; j++)
      {
	if (copy_token (token, &next_char, &i) == EMPTY)
	  break;
	if (sscanf (token, SCANFORMAT, &s_ptr->cd_music[j]) != 1)
	  break;
      }
      s_ptr->dz[0] =
	s_ptr->cd_music[0] + s_ptr->cd_music[3] * s_ptr->cd_music[4];
      s_ptr->dz[1] =
	s_ptr->cd_music[1] + (1 - s_ptr->cd_music[3]) * s_ptr->cd_music[4];
      s_ptr->dz[2] = s_ptr->cd_music[2];
      opt_save = OPTION_DEFAULT;
      break;
    case OPTION_DEFAULT:
      /*
       *   Get surface species information and parse equation
       */
      s_ptr = NULL;
      if (parse_eq (line, &next_elt, association) == ERROR)
      {
	parse_error++;
	error_msg ("Parsing equation.", CONTINUE);
	error_msg (line_save, CONTINUE);
	break;
      }
      /*
       *   Get pointer to each species in the reaction, store new species if necessary
       */
      trxn.token[0].s = s_store (trxn.token[0].name, trxn.token[0].z, TRUE);
      for (i = 1; i < count_trxn; i++)
      {
	trxn.token[i].s =
	  s_store (trxn.token[i].name, trxn.token[i].z, FALSE);
      }
      /*
       *   Save element list and carbon, hydrogen, and oxygen in species
       */
      trxn.token[0].s->next_elt = next_elt;
      for (; next_elt->elt != NULL; next_elt++)
      {
	if (strcmp (next_elt->elt->name, "C") == 0)
	{
	  trxn.token[0].s->carbon = next_elt->coef;
	}
	if (strcmp (next_elt->elt->name, "H") == 0)
	{
	  trxn.token[0].s->h = next_elt->coef;
	}
	if (strcmp (next_elt->elt->name, "O") == 0)
	{
	  trxn.token[0].s->o = next_elt->coef;
	}
      }
      /*
       *   Find coefficient of surface in rxn, store in equiv
       */
      trxn.token[0].s->equiv = 0.0;
      for (i = 1; i < count_trxn; i++)
      {
	if (trxn.token[i].s->type == SURF)
	{
	  trxn.token[0].s->equiv = trxn.token[i].coef;
	}
      }
      if (trxn.token[0].s->equiv == 0.0)
	trxn.token[0].s->equiv = 1.0;
      /*
       *   Malloc space for species reaction
       */
      trxn.token[0].s->rxn = rxn_alloc (count_trxn + 1);
      /*
       *   Copy reaction to reaction for species
       */
      token_ptr = trxn.token[0].s->rxn->token;
      for (i = 0; i < count_trxn; i++)
      {
	token_ptr[i].s = trxn.token[i].s;
	token_ptr[i].coef = trxn.token[i].coef;
      }
      token_ptr[i].s = NULL;
      /*
       *   Set type for species
       */
      trxn.token[0].s->type = SURF;
      s_ptr = trxn.token[0].s;
      /*
       *   Read gamma data
       */
      s_ptr->gflag = 6;
      s_ptr->dha = 0.0;
      s_ptr->dhb = 0.0;
      opt_save = OPTION_DEFAULT;
      break;
    }
    if (return_value == EOF || return_value == KEYWORD)
      break;
  }
  return (return_value);
}

/* ---------------------------------------------------------------------- */
int
read_surf (void)
/* ---------------------------------------------------------------------- */
{
  /*
   *      Reads surface data
   *
   *      Arguments:
   *         none
   *
   *      Returns:
   *         KEYWORD if keyword encountered, input_error may be incremented if
   *                    a keyword is encountered in an unexpected position
   *         EOF     if eof encountered while reading mass balance concentrations
   *         ERROR   if error occurred reading data
   *
   */
  int i, j, n, l;
  int count_comps, count_charge;
  int n_user, n_user_end;
  LDBLE conc, area, grams, thickness;
  char *ptr, *ptr1;
  char *description;
  char token[MAX_LENGTH], token1[MAX_LENGTH], name[MAX_LENGTH];
  struct surface *surface_ptr;
  struct surface_charge *charge_ptr;
  struct surface_comp *comp_ptr;

  int return_value, opt;
  char *next_char;
  const char *opt_list[] = {
    "equilibrate",		/* 0 */
    "equil",			/* 1 */
    "diff",			/* 2 */
    "diffuse_layer",		/* 3 */
    "no_edl",			/* 4 */
    "no_electrostatic",		/* 5 */
    "only_counter_ions",	/* 6 */
    "donnan",			/* 7 */
    "cd_music",			/* 8 */
    "capacitances",		/* 9 */
    "sites",			/* 10 */
    "sites_units"		/* 11 */
  };
  int count_opt_list = 12;
  /*
   * kin_surf is for Surfaces, related to kinetically reacting minerals
   *    they are defined if "sites" is followed by mineral name:
   *    Surf_wOH  Manganite  [equilibrium_phases or kinetics]      0.25         4000
   *    ^Name     mineral    ^switch                               ^prop.factor ^m2/mol 
   */
  /*
   *   Read surface number and description
   */
  ptr = line;
  read_number_description (ptr, &n_user, &n_user_end, &description);
  /*
   *   Find space for surface data
   */
  surface_ptr = surface_search (n_user, &n, FALSE);
  if (surface_ptr != NULL)
  {
    surface_free (&surface[n]);
  }
  else
  {
    n = count_surface++;
    space ((void **) ((void *) &surface), count_surface, &max_surface,
	   sizeof (struct surface));
  }
  /*
   *   Initialize
   */
  surface_init (&(surface[n]), n_user, n_user_end, description);
  free_check_null (description);

  if (use.surface_in == FALSE)
  {
    use.surface_in = TRUE;
    use.n_surface_user = n_user;
  }
  /*
   *   Read surface data
   */
  count_charge = 0;
  count_comps = 0;
  return_value = UNKNOWN;
  comp_ptr = NULL;
  charge_ptr = NULL;
  for (;;)
  {
    opt = get_option (opt_list, count_opt_list, &next_char);
    switch (opt)
    {
    case OPTION_EOF:		/* end of file */
      return_value = EOF;
      break;
    case OPTION_KEYWORD:	/* keyword */
      return_value = KEYWORD;
      break;
    case OPTION_ERROR:
      input_error++;
      error_msg ("Unknown input in SURFACE keyword.", CONTINUE);
      error_msg (line_save, CONTINUE);
      break;
    case 0:			/* equilibrate */
    case 1:
      for (;;)
      {
	i = copy_token (token, &next_char, &l);
	if (i == DIGIT)
	{
	  sscanf (token, "%d", &surface[n].n_solution);
	  surface[n].solution_equilibria = TRUE;
	  break;
	}
	if (i == EMPTY)
	{
	  error_msg
	    ("Expected a solution number with which to equilibrate surface.",
	     CONTINUE);
	  error_msg (line_save, CONTINUE);
	  input_error++;
	  break;
	}
      }
      break;
    case 2:			/* diffuse_layer */
    case 3:
      surface[n].thickness = 1e-8;
      surface[n].dl_type = BORKOVEK_DL;
      sscanf (next_char, SCANFORMAT, &surface[n].thickness);
      /*                              surface[n].thickness = thickness;
         }
       */ break;
    case 4:			/* no electrostatic */
    case 5:
      /*surface[n].edl = FALSE; */
      surface[n].type = NO_EDL;
      break;
    case 6:
      surface[n].only_counter_ions = TRUE;
      break;
    case 7:			/* donnan for DL conc's */
      surface[n].dl_type = DONNAN_DL;
      /*surface[n].diffuse_layer = TRUE; */
      thickness = 0.0;
      for (;;)
      {
	i = copy_token (token, &next_char, &l);
	if (i == DIGIT)
	{
	  sscanf (token, SCANFORMAT, &surface[n].thickness);
	  thickness = 1;
	  continue;
	}
	else if (i != EMPTY)
	{
	  if (token[0] == 'D' || token[0] == 'd')
	  {
	    if (thickness != 0)
	    {
	      error_msg
		("You must enter EITHER thickness OR Debye lengths (1/k),\n	   and relative DDL viscosity, DDL limit.\nCorrect is (for example): -donnan 1e-8 viscosity 0.5\n or (default values):     -donnan debye_lengths 1 viscosity 1 limit 0.8",
		 CONTINUE);
	      error_msg (line_save, CONTINUE);
	      input_error++;
	      break;
	    }
	    j = copy_token (token1, &next_char, &l);
	    if (j == DIGIT)
	    {
	      sscanf (token1, SCANFORMAT, &surface[n].debye_lengths);
	      continue;
	    }
	    else if (j != EMPTY)
	    {
	      error_msg ("Expected number of Debye lengths (1/k).", CONTINUE);
	      error_msg (line_save, CONTINUE);
	      input_error++;
	      break;
	    }
	  }
	  else if (token[0] == 'V' || token[0] == 'v')
	  {
	    j = copy_token (token1, &next_char, &l);
	    if (j == DIGIT)
	    {
	      sscanf (token1, SCANFORMAT, &surface[n].DDL_viscosity);
	      continue;
	    }
	    else if (j != EMPTY)
	    {
	      error_msg ("Expected number for relative DDL viscosity.",
			 CONTINUE);
	      error_msg (line_save, CONTINUE);
	      input_error++;
	      break;
	    }
	  }
	  else if (token[0] == 'L' || token[0] == 'l')
	  {
	    j = copy_token (token1, &next_char, &l);
	    if (j == DIGIT)
	    {
	      sscanf (token1, SCANFORMAT, &surface[n].DDL_limit);
	      continue;
	    }
	    else if (j != EMPTY)
	    {
	      error_msg
		("Expected number for maximum of DDL water as fraction of bulk water.",
		 CONTINUE);
	      error_msg (line_save, CONTINUE);
	      input_error++;
	      break;
	    }
	  }
	  else
	  {
	    error_msg
	      ("Expected diffuse layer thickness (m) or Debye lengths (1/k) for \n\tcalculating the thickness, and relative DDL viscosity and DDL limit.\nCorrect is (for example): -donnan 1e-8 visc 0.5\n or (default values):     -donnan debye_lengths 1 visc 1 lim 0.8",
	       CONTINUE);
	    error_msg (line_save, CONTINUE);
	    input_error++;
	    break;
	  }
	}
	else
	  break;
      }
      break;
    case 8:			/* cd_music */
      /*surface[n].edl = FALSE; */
      surface[n].type = CD_MUSIC;
      break;
    case 9:			/* capacitances */
      /*
       *   Read capacitance for CD_MUSIC
       */
      if (charge_ptr == NULL)
      {
	error_msg ("Surface component has not been defined for capacitances.",
		   CONTINUE);
	error_msg (line_save, CONTINUE);
	input_error++;
	break;
      }
      copy_token (token1, &next_char, &l);
      if (sscanf (token1, SCANFORMAT, &grams) == 1)
      {
	charge_ptr->capacitance[0] = grams;
      }
      copy_token (token1, &next_char, &l);
      if (sscanf (token1, SCANFORMAT, &grams) == 1)
      {
	charge_ptr->capacitance[1] = grams;
      }
      break;
    case 10:			/* sites */
    case 11:			/* sites_units */
      j = copy_token (token1, &next_char, &l);
      if (token1[0] == 'A' || token1[0] == 'a')
      {
	surface[n].sites_units = SITES_ABSOLUTE;
      }
      else if (token1[0] == 'D' || token1[0] == 'd')
      {
	surface[n].sites_units = SITES_DENSITY;
      }
      else
      {
	error_msg
	  ("Character string expected to be 'Absolute' or 'Density' to define the units of the first item in the definition of a surface component.\n",
	   CONTINUE);
	input_error++;
	break;
      }
      break;
    case OPTION_DEFAULT:
      /*
       *   Read surface component
       */
      ptr = line;
      i = copy_token (token, &ptr, &l);
      if (i != UPPER && token[0] != '[')
      {
	error_msg ("Expected surface name to begin with a capital letter.",
		   CONTINUE);
	error_msg (line_save, CONTINUE);
	input_error++;
	break;
      }
      /*
       *   Realloc space for new surface component
       */
      surface[n].comps =
	(struct surface_comp *) PHRQ_realloc (surface[n].comps,
					      (size_t) (count_comps +
							1) *
					      sizeof (struct surface_comp));
      if (surface[n].comps == NULL)
	malloc_error ();
      surface[n].comps[count_comps].formula = string_hsave (token);
      surface[n].comps[count_comps].formula_totals = NULL;
      surface[n].comps[count_comps].formula_z = 0.0;
      surface[n].comps[count_comps].moles = 0;
      surface[n].comps[count_comps].la = 0;
      surface[n].comps[count_comps].charge = 0;
      surface[n].comps[count_comps].cb = 0;
      surface[n].comps[count_comps].phase_name = NULL;
      surface[n].comps[count_comps].phase_proportion = 0;
      surface[n].comps[count_comps].rate_name = NULL;
      comp_ptr = &(surface[n].comps[count_comps]);
      surface[n].comps[count_comps].Dw = 0;
      i = copy_token (token1, &ptr, &l);
      if (i == DIGIT)
      {
	/*
	 *   Read surface concentration
	 */
	/* surface concentration is read directly */
	if (sscanf (token1, SCANFORMAT, &conc) < 1)
	{
	  error_msg ("Expected number of surface sites in moles.", CONTINUE);
	  error_msg (line_save, CONTINUE);
	  input_error++;
	  break;
	}
	surface[n].comps[count_comps].moles = conc;
	/*
	 *   Read equilibrium phase name or kinetics rate name
	 */
      }
      else if (i != EMPTY)
      {

	/* surface conc. is related to mineral or kinetics */
	surface[n].comps[count_comps].phase_name = string_hsave (token1);
	j = copy_token (token1, &ptr, &l);

	/* read optional 'equilibrium_phases' or 'kinetics' */
	if (j == DIGIT)
	{
	  surface[n].related_phases = TRUE;
	}
	else
	{
	  if (token1[0] == 'K' || token1[0] == 'k')
	  {
	    surface[n].comps[count_comps].rate_name =
	      surface[n].comps[count_comps].phase_name;
	    surface[n].comps[count_comps].phase_name = NULL;
	    surface[n].related_rate = TRUE;
	  }
	  else if (token1[0] == 'E' || token1[0] == 'e')
	  {
	    surface[n].related_phases = TRUE;
	    surface[n].comps[count_comps].rate_name = NULL;
	  }
	  else
	  {
	    error_msg
	      ("Character string expected to be 'equilibrium_phase' or 'kinetics' to relate surface to mineral or kinetic reaction.\n",
	       CONTINUE);
	    input_error++;
	    break;
	  }
	  j = copy_token (token1, &ptr, &l);
	}

	/* read proportion */
	if (j != DIGIT)
	{
	  error_msg
	    ("Expected a coefficient to relate surface to mineral or kinetic reaction.\n",
	     CONTINUE);
	  input_error++;
	  break;
	}
	sscanf (token1, SCANFORMAT,
		&surface[n].comps[count_comps].phase_proportion);
	/* real conc must be defined in tidy_model */
	conc = 1.0;
      }
      else
      {
	error_msg
	  ("Expected concentration of surface, mineral name, or kinetic reaction name.",
	   CONTINUE);
	error_msg (line_save, CONTINUE);
	input_error++;
	break;
      }
      /*
       *   Accumulate elements in elt_list
       */
      count_elts = 0;
      paren_count = 0;
      ptr1 = token;
      get_elts_in_species (&ptr1, conc);

      /*
       *   save formula for adjusting number of exchange sites
       */

      ptr1 = token;
      get_token (&ptr1, token1, &surface[n].comps[count_comps].formula_z, &l);
      surface[n].comps[count_comps].formula_totals = elt_list_save ();
      surface[n].comps[count_comps].totals = elt_list_save ();
      /*
       *   Search for charge structure
       */
      ptr1 = token;
      get_elt (&ptr1, name, &l);
      ptr1 = strchr (name, '_');
      if (ptr1 != NULL)
	ptr1[0] = '\0';
      for (i = 0; i < count_charge; i++)
      {
	if (strcmp (surface[n].charge[i].name, name) == 0)
	  break;
      }
      if (i >= count_charge)
      {
	surface[n].charge =
	  (struct surface_charge *) PHRQ_realloc (surface[n].charge,
						  (size_t) (count_charge +
							    1) *
						  sizeof (struct
							  surface_charge));
	if (surface[n].charge == NULL)
	  malloc_error ();
	i = count_charge;
	surface[n].charge[i].name = string_hsave (name);
	if (surface[n].comps[count_comps].phase_name == NULL
	    && surface[n].comps[count_comps].rate_name == NULL)
	{
	  surface[n].charge[i].specific_area = 600.0;
	  surface[n].charge[i].grams = 0.0;
	}
	else
	{
	  surface[n].charge[i].specific_area = 0.0;
	  surface[n].charge[i].grams = 1.0;
	}
	surface[n].charge[i].charge_balance = 0.0;
	surface[n].charge[i].mass_water = 0.0;
	surface[n].charge[i].diffuse_layer_totals = NULL;
	surface[n].charge[i].count_g = 0;
	surface[n].charge[i].g = NULL;
	surface[n].charge[i].la_psi = 0;
	surface[n].charge[i].la_psi1 = 0;
	surface[n].charge[i].la_psi2 = 0;
	surface[n].charge[i].psi = 0;
	surface[n].charge[i].psi1 = 0;
	surface[n].charge[i].psi2 = 0;
	surface[n].charge[i].capacitance[0] = 1.0;
	surface[n].charge[i].capacitance[1] = 5.0;
	surface[n].charge[i].sigma0 = 0;
	surface[n].charge[i].sigma1 = 0;
	surface[n].charge[i].sigma2 = 0;
	count_charge++;
      }
      charge_ptr = &(surface[n].charge[i]);
      surface[n].comps[count_comps].charge = i;
      count_comps++;
      /*
       *   Read surface area (m2/g)
       */
      copy_token (token1, &ptr, &l);
      if (sscanf (token1, SCANFORMAT, &area) == 1)
      {
	surface[n].charge[i].specific_area = area;
      }
      else
      {
	break;
      }
      /*
       *   Read grams of solid (g)
       */
      copy_token (token1, &ptr, &l);
      if (sscanf (token1, SCANFORMAT, &grams) == 1)
      {
	surface[n].charge[i].grams = grams;
      }
      /* read Dw */
      copy_token (token1, &ptr, &l);
      str_tolower (token1);
      if (strcmp (token1, "dw") == 0)
      {
	j = copy_token (token1, &ptr, &l);
	if (j != DIGIT)
	{
	  error_msg ("Expected surface diffusion coefficient (m2/s).\n",
		     CONTINUE);
	  input_error++;
	  break;
	}
	else
	{
	  sscanf (token1, SCANFORMAT, &surface[n].comps[count_comps - 1].Dw);
	  if (surface[n].comps[count_comps - 1].Dw > 0)
	  {
	    surface[n].transport = TRUE;
	    if (surface[n].related_rate == TRUE
		|| surface[n].related_phases == TRUE)
	    {
	      error_msg
		("Can't transport surfaces related to phases or rates (yet).",
		 CONTINUE);
	      input_error++;
	    }
	  }
	  break;
	}
      }
      break;

    }
    if (return_value == EOF || return_value == KEYWORD)
      break;
  }
  surface[n].count_comps = count_comps;
  surface[n].count_charge = count_charge;
  /*
   * copy Dw > 0 to all comps with the same charge structure, check edl & dif for surface transport
   */
  if (surface[n].transport == TRUE)
  {
    /*
       if (surface[n].edl == FALSE || surface[n].diffuse_layer == FALSE) {
       surface[n].edl = TRUE;
       surface[n].diffuse_layer = TRUE;
       warning_msg("Can't transport without edl and diffuse layer, set to true");
       }
     */
    if (surface[n].type <= NO_EDL)
    {
      input_error++;
      error_msg
	("Must use default Dzombak and Morel or -cd_music for surface transport.",
	 CONTINUE);
    }
    if (surface[n].dl_type <= NO_DL)
    {
      input_error++;
      error_msg ("Must use -donnan or -diffuse_layer for surface transport.",
		 CONTINUE);
    }
    for (i = 0; i < count_comps; i++)
    {
      if (surface[n].comps[i].Dw > 0)
      {
	for (j = 0; j < count_comps; j++)
	{
	  if (surface[n].comps[j].charge == surface[n].comps[i].charge)
	    surface[n].comps[j].Dw = surface[n].comps[i].Dw;
	}
      }
    }
  }
  /*
   *   Make sure surface area is defined
   */
  /*if (surface[n].edl == TRUE) { */
  if (surface[n].type == DDL || surface[n].type == CD_MUSIC)
  {
    for (i = 0; i < count_charge; i++)
    {
      if (surface[n].charge[i].grams * surface[n].charge[i].specific_area <=
	  0)
      {
	sprintf (error_string, "Surface area not defined for %s.\n",
		 surface[n].charge[i].name);
	error_msg (error_string, CONTINUE);
	input_error++;
      }
    }
  }
  /*
   *  Logical checks
   */
  if (surface[n].type == UNKNOWN_DL)
  {
    sprintf (error_string, "Unknown surface type.\n");
    error_msg (error_string, CONTINUE);
    input_error++;
  }
  else if (surface[n].type == NO_EDL)
  {
    if (surface[n].dl_type != NO_DL)
    {
      sprintf (error_string,
	       "No electrostatic term calculations do not allow calculation of the diffuse layer composition.\n");
      warning_msg (error_string);
      surface[n].dl_type = NO_DL;
    }
  }
  else if (surface[n].type == DDL)
  {
    /* all values of dl_type are valid */
  }
  else if (surface[n].type == CD_MUSIC)
  {
    if (surface[n].dl_type == BORKOVEK_DL)
    {
      sprintf (error_string,
	       "Borkovec and Westall diffuse layer calculation is not allowed with a CD_MUSIC surface.\n\tUsing Donnan diffuse layer calculation.");
      warning_msg (error_string);
      surface[n].dl_type = DONNAN_DL;
    }
    if (surface[n].debye_lengths > 0)
    {
      error_msg
	("Variable DDL thickness is not permitted in CD_MUSIC. Fix DDL thickness\n\tfor example (default value): -donnan 1e-8",
	 CONTINUE);
      input_error++;
    }
  }
  return (return_value);
}

/* ---------------------------------------------------------------------- */
int
read_surface_master_species (void)
/* ---------------------------------------------------------------------- */
{
  /*
   *   Reads master species data from data file or input file
   */
  int l, return_value;
  char *ptr, *ptr1;
  LDBLE z;
  struct master *m_ptr;
  struct species *s_ptr;
  char token[MAX_LENGTH], token1[MAX_LENGTH];
  int opt, opt_save;
  char *next_char;
  const char *opt_list[] = {
    "capacitance",		/* 0 */
    "cd_music_capacitance"	/* 1 */
  };
  int count_opt_list = 0;
  opt_save = OPTION_DEFAULT;
  return_value = UNKNOWN;
  m_ptr = NULL;
  for (;;)
  {
    opt = get_option (opt_list, count_opt_list, &next_char);
    if (opt == OPTION_DEFAULT)
    {
      opt = opt_save;
    }
    switch (opt)
    {
    case OPTION_EOF:		/* end of file */
      return_value = EOF;
      break;
    case OPTION_KEYWORD:	/* keyword */
      return_value = KEYWORD;
      break;
    case OPTION_ERROR:
      input_error++;
      error_msg ("Unknown input in SURFACE_SPECIES keyword.", CONTINUE);
      error_msg (line_save, CONTINUE);
      break;
#ifdef SKIP
    case 0:			/* capacitance */
    case 1:			/* cd_music_capacitance */
      if (m_ptr == NULL)
      {
	sprintf (error_string, "No master_species defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      for (j = 0; j < 2; j++)
      {
	if (copy_token (token, &next_char, &l) == EMPTY)
	  break;
	if (sscanf (token, SCANFORMAT, &m_ptr->capacitance[j]) != 1)
	  break;
      }
      m_ptr->capacitance_defined = TRUE;
      opt_save = OPTION_DEFAULT;

      break;
#endif
    case OPTION_DEFAULT:
      /*
       *   Get "element" name with valence, allocate space, store
       */
      ptr = line;
      /*
       *   Get element name and save pointer to character string
       */
      if (copy_token (token, &ptr, &l) != UPPER && token[0] != '[')
      {
	parse_error++;
	error_msg ("Reading element for master species.", CONTINUE);
	error_msg (line_save, CONTINUE);
	continue;
      }
      replace ("(+", "(", token);
      /*
       *   Delete master if it exists
       */
      master_delete (token);
      /*
       *   Increase pointer array, if necessary,  and malloc space
       */
      if (count_master + 4 >= max_master)
      {
	space ((void **) ((void *) &master), count_master + 4, &max_master,
	       sizeof (struct master *));
      }
      /*
       *   Save values in master and species structure for surface sites
       */
      master[count_master] = master_alloc ();
      master[count_master]->type = SURF;
      master[count_master]->elt = element_store (token);
      m_ptr = master[count_master];
      if (copy_token (token, &ptr, &l) != UPPER && token[0] != '[')
      {
	parse_error++;
	error_msg ("Reading surface master species name.", CONTINUE);
	error_msg (line_save, CONTINUE);
	continue;
      }
      s_ptr = s_search (token);
      if (s_ptr != NULL)
      {
	master[count_master]->s = s_ptr;
      }
      else
      {
	ptr1 = token;
	get_token (&ptr1, token1, &z, &l);
	master[count_master]->s = s_store (token1, z, FALSE);
      }
      master[count_master]->primary = TRUE;
      strcpy (token, master[count_master]->elt->name);
      count_master++;
      /*
       *   Save values in master and species structure for surface psi
       */
      strcpy (token1, token);
      replace ("_", " ", token1);
      ptr1 = token1;
      copy_token (token, &ptr1, &l);
      strcat (token, "_psi");
      add_psi_master_species (token);
#ifdef SKIP
      if (master_ptr == NULL)
      {
	master[count_master] = master_alloc ();
	master[count_master]->type = SURF_PSI;
	master[count_master]->elt = element_store (token);
	s_ptr = s_search (token);
	if (s_ptr != NULL)
	{
	  master[count_master]->s = s_ptr;
	}
	else
	{
	  master[count_master]->s = s_store (token, 0.0, FALSE);
	}
	count_elts = 0;
	paren_count = 0;
	ptr = token;
	get_elts_in_species (&ptr, 1.0);
	master[count_master]->s->next_elt = elt_list_save ();
	master[count_master]->s->type = SURF_PSI;
	master[count_master]->primary = TRUE;
	master[count_master]->s->rxn = rxn_alloc (3);
	/*
	 *   Define reaction for psi
	 */
	for (i = 0; i < 8; i++)
	{
	  master[count_master]->s->rxn->logk[i] = 0.0;
	}
	master[count_master]->s->rxn->token[0].s = master[count_master]->s;
	master[count_master]->s->rxn->token[0].coef = -1.0;
	master[count_master]->s->rxn->token[1].s = master[count_master]->s;
	master[count_master]->s->rxn->token[1].coef = 1.0;
	master[count_master]->s->rxn->token[2].s = NULL;
	count_master++;
      }
#endif
      opt_save = OPTION_DEFAULT;
      break;
    }
    if (return_value == EOF || return_value == KEYWORD)
      break;
  }
  return (return_value);
}

/* ---------------------------------------------------------------------- */
int
add_psi_master_species (char *token)
/* ---------------------------------------------------------------------- */
{
  struct species *s_ptr;
  struct master *master_ptr;
  char *ptr;
  char token1[MAX_LENGTH];
  int i, n, plane;

  strcpy (token1, token);
  for (plane = SURF_PSI; plane <= SURF_PSI2; plane++)
  {
    strcpy (token, token1);
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
    }
    master_ptr = master_search (token, &n);
    if (master_ptr == NULL)
    {
      master[count_master] = master_alloc ();
      master[count_master]->type = plane;
      master[count_master]->elt = element_store (token);
      s_ptr = s_search (token);
      if (s_ptr != NULL)
      {
	master[count_master]->s = s_ptr;
      }
      else
      {
	master[count_master]->s = s_store (token, 0.0, FALSE);
      }
      count_elts = 0;
      paren_count = 0;
      ptr = token;
      get_elts_in_species (&ptr, 1.0);
      master[count_master]->s->next_elt = elt_list_save ();
      master[count_master]->s->type = plane;
      master[count_master]->primary = TRUE;
      master[count_master]->s->rxn = rxn_alloc (3);
      /*
       *   Define reaction for psi
       */
      for (i = 0; i < 8; i++)
      {
	master[count_master]->s->rxn->logk[i] = 0.0;
      }
      master[count_master]->s->rxn->token[0].s = master[count_master]->s;
      master[count_master]->s->rxn->token[0].coef = -1.0;
      master[count_master]->s->rxn->token[1].s = master[count_master]->s;
      master[count_master]->s->rxn->token[1].coef = 1.0;
      master[count_master]->s->rxn->token[2].s = NULL;
      count_master++;
    }
  }
  return (OK);
}

#ifdef SKIP
/* ---------------------------------------------------------------------- */
int
read_surface_master_species (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Reads master species data from data file or input file
 */
  int i, l, n, return_value;
  char *ptr, *ptr1;
  LDBLE z;
  struct master *master_ptr;
  struct species *s_ptr;
  char token[MAX_LENGTH], token1[MAX_LENGTH];

  for (;;)
  {
    return_value =
      check_line ("Surface species equation", FALSE, TRUE, TRUE, TRUE);
    if (return_value == EOF || return_value == KEYWORD)
      break;
/*
 *   Get "element" name with valence, allocate space, store
 */
    ptr = line;
/*
 *   Get element name and save pointer to character string
 */
    if (copy_token (token, &ptr, &l) != UPPER && token[0] != '[')
    {
      parse_error++;
      error_msg ("Reading element for master species.", CONTINUE);
      error_msg (line_save, CONTINUE);
      continue;
    }
    /*
       if (token[0] == '[') {
       ptr1 = token;
       get_elt(&ptr, element, &l);
       strcpy(token, element);
       }
     */
    replace ("(+", "(", token);
/*
 *   Delete master if it exists
 */
    master_delete (token);
/*
 *   Increase pointer array, if necessary,  and malloc space
 */
    if (count_master + 2 >= max_master)
    {
      space ((void **) ((void *) &master), count_master + 2, &max_master,
	     sizeof (struct master *));
    }
/*
 *   Save values in master and species structure for surface sites
 */
    master[count_master] = master_alloc ();
    master[count_master]->type = SURF;
    master[count_master]->elt = element_store (token);

    if (copy_token (token, &ptr, &l) != UPPER && token[0] != '[')
    {
      parse_error++;
      error_msg ("Reading surface master species name.", CONTINUE);
      error_msg (line_save, CONTINUE);
      continue;
    }
    s_ptr = s_search (token);
    if (s_ptr != NULL)
    {
      master[count_master]->s = s_ptr;
    }
    else
    {
      ptr1 = token;
      get_token (&ptr1, token1, &z, &l);
      master[count_master]->s = s_store (token1, z, FALSE);
    }
    master[count_master]->primary = TRUE;
    strcpy (token, master[count_master]->elt->name);
    count_master++;
/*
 *   Save values in master and species structure for surface psi
 */
    strcpy (token1, token);
    replace ("_", " ", token1);
    ptr1 = token1;
    copy_token (token, &ptr1, &l);
    strcat (token, "_psi");
    master_ptr = master_search (token, &n);
    if (master_ptr == NULL)
    {
      master[count_master] = master_alloc ();
      master[count_master]->type = SURF_PSI;
      master[count_master]->elt = element_store (token);
      s_ptr = s_search (token);
      if (s_ptr != NULL)
      {
	master[count_master]->s = s_ptr;
      }
      else
      {
	master[count_master]->s = s_store (token, 0.0, FALSE);
      }
      count_elts = 0;
      paren_count = 0;
      ptr = token;
      get_elts_in_species (&ptr, 1.0);
      master[count_master]->s->next_elt = elt_list_save ();

      master[count_master]->s->type = SURF_PSI;
      master[count_master]->primary = TRUE;
      master[count_master]->s->rxn = rxn_alloc (3);
/*
 *   Define reaction for psi
 */
      for (i = 0; i < 8; i++)
      {
	master[count_master]->s->rxn->logk[i] = 0.0;
      }
      master[count_master]->s->rxn->token[0].s = master[count_master]->s;
      master[count_master]->s->rxn->token[0].coef = -1.0;
      master[count_master]->s->rxn->token[1].s = master[count_master]->s;
      master[count_master]->s->rxn->token[1].coef = 1.0;
      master[count_master]->s->rxn->token[2].s = NULL;
      count_master++;
    }
  }
  return (return_value);
}
#endif
/* ---------------------------------------------------------------------- */
int
read_temperature (void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads temperature data for reaction steps
 *
 *      Arguments:
 *         none
 *
 *      Returns:
 *         KEYWORD if keyword encountered, input_error may be incremented if
 *                    a keyword is encountered in an unexpected position
 *         EOF     if eof encountered while reading mass balance concentrations
 *         ERROR   if error occurred reading data
 *
 */
  char *ptr;
  char *description;
  int n, return_value;
  int n_user, n_user_end;
  struct temperature *temperature_ptr;
/*
 *   Read reaction number
 */
  ptr = line;
  read_number_description (ptr, &n_user, &n_user_end, &description);
/*
 *   Read allocate space for temperature
 */
  temperature_ptr = temperature_search (n_user, &n);
  if (temperature_ptr != NULL)
  {
    temperature_free (&temperature[n]);
  }
  else
  {
    n = count_temperature++;
    temperature =
      (struct temperature *) PHRQ_realloc (temperature,
					   (size_t) count_temperature *
					   sizeof (struct temperature));
    if (temperature == NULL)
      malloc_error ();
  }
/*
 *   Set use data to first read
 */
  if (use.temperature_in == FALSE)
  {
    use.temperature_in = TRUE;
    use.n_temperature_user = n_user;
  }
/*
 *   Defaults
 */
  temperature[n].n_user = n_user;
  temperature[n].n_user_end = n_user_end;
  temperature[n].description = description;
  temperature[n].t = (LDBLE *) PHRQ_malloc ((size_t) sizeof (LDBLE));
  if (temperature[n].t == NULL)
    malloc_error ();
  temperature[n].t[0] = 25.0;
  temperature[n].count_t = 0;
/*
 *   Read temperature data
 */
  for (;;)
  {
    return_value =
      check_line ("reaction_temperature", FALSE, TRUE, TRUE, TRUE);
    /* empty, eof, keyword, print */
    if (return_value == EOF || return_value == KEYWORD)
      break;
/*
 *   Read steps information
 */
    read_reaction_temps (&(temperature[n]));
  }
  return (return_value);
}

/* ---------------------------------------------------------------------- */
int
read_reaction_temps (struct temperature *temperature_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Read temperature one of two forms:
 *
 *   (1) 25.0  30. 40. 50. 60.
 *
 *   (2) 25.0 100.0 in 4 steps # (25 50 75 100)
 *
 */
  int i, j, l;
  int count_t;
  char *ptr;
  char token[MAX_LENGTH];
  LDBLE step;

  ptr = line;
  count_t = temperature_ptr->count_t;
/*
 *   Read one or more reaction increments
 */
  for (;;)
  {
    if (copy_token (token, &ptr, &l) == EMPTY)
    {
      return (OK);
    }
/*
 *   Read next step increment
 */
    j = sscanf (token, SCANFORMAT, &step);
    if (j == 1)
    {
      count_t++;
      temperature_ptr->t =
	(LDBLE *) PHRQ_realloc (temperature_ptr->t,
				(size_t) count_t * sizeof (LDBLE));
      if (temperature_ptr->t == NULL)
	malloc_error ();
      temperature_ptr->t[temperature_ptr->count_t] = step;
      temperature_ptr->count_t = count_t;
    }
    else
    {
      break;
    }
  }
/*
 *  Read number of equal increments, store as negative integer
 */
  if (count_t != 2)
  {
    error_msg
      ("To define equal increments, exactly two temperatures should be defined.",
       CONTINUE);
    error_msg (line_save, CONTINUE);
    input_error++;
    return (ERROR);
  }
  do
  {
    j = sscanf (token, "%d", &i);
    if (j == 1 && i > 0)
    {
      temperature_ptr->count_t = -i;
      return (OK);
    }
    else if (j == 1 && i <= 0)
    {
      break;
    }
  }
  while (copy_token (token, &ptr, &l) != EMPTY);

  error_msg ("Expecting positive number for calculation of "
	     "temperature increments.", CONTINUE);
  error_msg (line_save, CONTINUE);
  input_error++;
  return (ERROR);
}

/* ---------------------------------------------------------------------- */
int
read_title (void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads title for simulation
 *
 *      Arguments:
 *         none
 *
 *      Returns:
 *         KEYWORD if keyword encountered, input_error may be incremented if
 *                    a keyword is encountered in an unexpected position
 *         EOF     if eof encountered while reading mass balance concentrations
 *         ERROR   if error occurred reading data
 *
 */
  char *ptr, *ptr1;
  int l, title_x_length, line_length;
  int return_value;
  char token[MAX_LENGTH];
/*
 *   Read anything after keyword
 */
  ptr = line;
  copy_token (token, &ptr, &l);
  ptr1 = ptr;
  title_x = (char *) free_check_null (title_x);
  if (copy_token (token, &ptr, &l) != EMPTY)
  {
    title_x = string_duplicate (ptr1);
  }
  else
  {
    title_x = (char *) PHRQ_malloc (sizeof (char));
    if (title_x == NULL)
      malloc_error ();
    title_x[0] = '\0';
  }

/*
 *   Read additonal lines
 */
  for (;;)
  {
    return_value = check_line ("title", TRUE, TRUE, TRUE, TRUE);
    /* empty, eof, keyword, print */
    if (return_value == EOF || return_value == KEYWORD)
      break;
/*
 *   append line to title_x
 */
    title_x_length = (int) strlen (title_x);
    line_length = (int) strlen (line);
    title_x =
      (char *) PHRQ_realloc (title_x,
			     (size_t) (title_x_length + line_length +
				       2) * sizeof (char));
    if (title_x == NULL)
      malloc_error ();
    if (title_x_length > 0)
    {
      title_x[title_x_length] = '\n';
      title_x[title_x_length + 1] = '\0';
    }
    strcat (title_x, line);
  }
  return (return_value);
}

/* ---------------------------------------------------------------------- */
int
read_advection (void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads advection information
 *
 *      Arguments:
 *         none
 *
 *      Returns:
 *         KEYWORD if keyword encountered, input_error may be incremented if
 *                    a keyword is encountered in an unexpected position
 *         EOF     if eof encountered while reading mass balance concentrations
 *         ERROR   if error occurred reading data
 *
 */
/*
 *   Read advection parameters: 
 *        number of cells;
 *        number of shifts;
 */
  char *ptr;
  char *description;
  int n_user, n_user_end, i;

  int count_punch, count_print;
  int *punch_temp, *print_temp;
  int return_value, opt, opt_save;
  char *next_char;
  const char *opt_list[] = {
    "cells",			/* 0 */
    "shifts",			/* 1 */
    "print",			/* 2 */
    "selected_output",		/* 3 */
    "punch",			/* 4 */
    "print_cells",		/* 5 */
    "selected_cells",		/* 6 */
    "time_step",		/* 7 */
    "timest",			/* 8 */
    "output",			/* 9 */
    "output_frequency",		/* 10 */
    "selected_output_frequency",	/* 11 */
    "punch_frequency",		/* 12 */
    "print_frequency",		/* 13 */
    "punch_cells",		/* 14 */
    "initial_time",		/* 15 */
    "warning",			/* 16 */
    "warnings"			/* 17 */
  };
  int count_opt_list = 18;
/*
 *   Read advection number (not currently used)
 */
  ptr = line;
  read_number_description (ptr, &n_user, &n_user_end, &description);
  description = (char *) free_check_null (description);
/*
 *   Set use data
 */
  use.advect_in = TRUE;
  count_ad_cells = 0;
  count_ad_shifts = 0;
  print_ad_modulus = 1;
  punch_ad_modulus = 1;
  count_punch = 0;
  count_print = 0;
  punch_temp = (int *) PHRQ_malloc (sizeof (int));
  if (punch_temp == NULL)
    malloc_error ();
  print_temp = (int *) PHRQ_malloc (sizeof (int));
  if (print_temp == NULL)
    malloc_error ();
/*
 *   Read lines
 */
  opt_save = OPTION_DEFAULT;
  return_value = UNKNOWN;
  for (;;)
  {
    opt = get_option (opt_list, count_opt_list, &next_char);
    if (opt == OPTION_DEFAULT)
    {
      opt = opt_save;
    }
    switch (opt)
    {
    case OPTION_EOF:		/* end of file */
      return_value = EOF;
      break;
    case OPTION_KEYWORD:	/* keyword */
      return_value = KEYWORD;
      break;
    case OPTION_DEFAULT:
    case OPTION_ERROR:
      input_error++;
      error_msg ("Unknown input in ADVECTION keyword.", CONTINUE);
      error_msg (line_save, CONTINUE);
      break;
    case 0:			/* cells */
      sscanf (next_char, "%d", &count_ad_cells);
      opt_save = OPTION_DEFAULT;
      break;
    case 1:			/* shifts */
      sscanf (next_char, "%d", &count_ad_shifts);
      opt_save = OPTION_DEFAULT;
      break;
    case 2:			/* print */
    case 5:			/* print_cells */
      print_temp =
	read_list_ints_range (&next_char, &count_print, TRUE, print_temp);
#ifdef SKIP
      while (copy_token (token, &next_char, &l) == DIGIT)
      {
	sscanf (token, "%d", &l);
	print_temp =
	  PHRQ_realloc (print_temp,
			(size_t) (count_print + 1) * sizeof (int));
	if (print_temp == NULL)
	  malloc_error ();
	print_temp[count_print] = l;
	count_print++;
      }
#endif
      opt_save = 2;
      break;
    case 3:			/* selected_output */
    case 11:			/* selected_output_frequency */
    case 12:			/* punch_frequency */
      sscanf (next_char, "%d", &punch_ad_modulus);
      opt_save = OPTION_DEFAULT;
      if (punch_ad_modulus <= 0) 
      {
	sprintf (error_string,
		 "Punch frequency must be greater than 0. Frequency set to 1000.");
	warning_msg (error_string);
	punch_ad_modulus = 1000;
      }
      break;
    case 4:			/* punch */
    case 14:			/* punch_cells */
    case 6:			/* selected_cells */
      punch_temp =
	read_list_ints_range (&next_char, &count_punch, TRUE, punch_temp);
#ifdef SKIP
      while (copy_token (token, &next_char, &l) == DIGIT)
      {
	sscanf (token, "%d", &l);
	punch_temp =
	  PHRQ_realloc (punch_temp,
			(size_t) (count_punch + 1) * sizeof (int));
	if (punch_temp == NULL)
	  malloc_error ();
	punch_temp[count_punch] = l;
	count_punch++;
      }
#endif
      opt_save = 4;
      break;
    case 7:			/* time_step */
    case 8:			/* timest */
      sscanf (next_char, SCANFORMAT, &advection_kin_time);
      advection_kin_time_defined = TRUE;
      opt_save = OPTION_DEFAULT;
      break;
    case 9:			/* output */
    case 10:			/* output_frequency */
    case 13:			/* print_frequency */
      sscanf (next_char, "%d", &print_ad_modulus);
      opt_save = OPTION_DEFAULT;
      if (print_ad_modulus <= 0) 
      {
	sprintf (error_string,
		 "Print frequency must be greater than 0. Frequency set to 1000.");
	warning_msg (error_string);
	print_ad_modulus = 1000;
      }
      break;
    case 15:			/* initial_time */
      sscanf (next_char, SCANFORMAT, &initial_total_time);
      opt_save = OPTION_DEFAULT;
      break;
    case 16:			/* warning */
    case 17:			/* warnings */
      advection_warnings = get_true_false (next_char, TRUE);
      break;
    }
    if (return_value == EOF || return_value == KEYWORD)
      break;
  }
/*
 *   Fill in data for punch
 */
  advection_punch =
    (int *) PHRQ_realloc (advection_punch,
			  (size_t) (count_ad_cells + 1) * sizeof (int));
  if (advection_punch == NULL)
    malloc_error ();
  if (count_punch != 0)
  {
    for (i = 0; i < count_ad_cells; i++)
      advection_punch[i] = FALSE;
    for (i = 0; i < count_punch; i++)
    {
      if (punch_temp[i] > count_ad_cells || punch_temp[i] < 1)
      {
	sprintf (error_string,
		 "Cell number for punch is out of range, %d. Request ignored.",
		 punch_temp[i]);
	warning_msg (error_string);
      }
      else
      {
	advection_punch[punch_temp[i] - 1] = TRUE;
      }
    }
  }
  else
  {
    for (i = 0; i < count_ad_cells; i++)
      advection_punch[i] = TRUE;
  }
  punch_temp = (int *) free_check_null (punch_temp);
/*
 *   Fill in data for print
 */
  advection_print =
    (int *) PHRQ_realloc (advection_print,
			  (size_t) (count_ad_cells + 1) * sizeof (int));
  if (advection_print == NULL)
    malloc_error ();
  if (count_print != 0)
  {
    for (i = 0; i < count_ad_cells; i++)
      advection_print[i] = FALSE;
    for (i = 0; i < count_print; i++)
    {
      if (print_temp[i] > count_ad_cells || print_temp[i] < 1)
      {
	sprintf (error_string,
		 "Cell number for print is out of range, %d. Request ignored.",
		 print_temp[i]);
	warning_msg (error_string);
      }
      else
      {
	advection_print[print_temp[i] - 1] = TRUE;
      }
    }
  }
  else
  {
    for (i = 0; i < count_ad_cells; i++)
      advection_print[i] = TRUE;
  }
  print_temp = (int *) free_check_null (print_temp);
  return (return_value);
}

/* ---------------------------------------------------------------------- */
int
read_debug (void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads knobs and debugging info
 *
 *      Arguments:
 *         none
 *
 *      Returns:
 *         KEYWORD if keyword encountered, input_error may be incremented if
 *                    a keyword is encountered in an unexpected position
 *         EOF     if eof encountered while reading mass balance concentrations
 *         ERROR   if error occurred reading data
 *
 */
  int return_value, opt;
  char *next_char;
  const char *opt_list[] = {
    "iterations",		/* 0 */
    "tolerance",		/* 1 */
    "step_size",		/* 2 */
    "pe_step_size",		/* 3 */
    "scale_pure_phases",	/* 4 */
    "diagonal_scale",		/* 5 */
    "debug_model",		/* 6 */
    "debug_prep",		/* 7 */
    "debug_set",		/* 8 */
    "debug_inverse",		/* 9 */
    "logfile",			/* 10 */
    "log_file",			/* 11 */
    "debug_diffuse_layer",	/* 12 */
    "delay_mass_water",		/* 13 */
    "convergence_tolerance",	/* 14 */
    "numerical_derivatives"	/* 15 */
  };
  int count_opt_list = 16;
/*
 *   Read parameters: 
 *        ineq_tol;
 *        step_size;
 *        pe_step_size;
 *        pp_scale;
 *        diagonal_scale;
 */
  return_value = UNKNOWN;
  for (;;)
  {
    opt = get_option (opt_list, count_opt_list, &next_char);
    switch (opt)
    {
    case OPTION_EOF:		/* end of file */
      return_value = EOF;
      break;
    case OPTION_KEYWORD:	/* keyword */
      return_value = KEYWORD;
      break;
    case OPTION_DEFAULT:
    case OPTION_ERROR:
      input_error++;
      error_msg ("Unknown input in KNOBS keyword.", CONTINUE);
      error_msg (line_save, CONTINUE);
      break;
    case 0:			/* iterations */
      sscanf (next_char, "%d", &itmax);
      break;
    case 1:			/* tolerance */
      sscanf (next_char, SCANFORMAT, &ineq_tol);
      break;
    case 2:			/* step_size */
      sscanf (next_char, SCANFORMAT, &step_size);
      break;
    case 3:			/* pe_step_size */
      sscanf (next_char, SCANFORMAT, &pe_step_size);
      break;
    case 4:			/* pp_scale */
      sscanf (next_char, SCANFORMAT, &pp_scale);
      break;
    case 5:			/* diagonal_scale */
      diagonal_scale = get_true_false (next_char, TRUE);
      break;
    case 6:			/* debug_model */
      debug_model = get_true_false (next_char, TRUE);
      break;
    case 7:			/* debug_prep */
      debug_prep = get_true_false (next_char, TRUE);
      break;
    case 8:			/* debug_set */
      debug_set = get_true_false (next_char, TRUE);
      break;
    case 9:			/* debug_inverse */
      debug_inverse = get_true_false (next_char, TRUE);
      break;
    case 10:			/* logfile */
    case 11:			/* log_file */
      pr.logfile = get_true_false (next_char, TRUE);
      if (phast == TRUE)
      {
	pr.logfile = FALSE;
	warning_msg ("PHREEQC log file is disabled in PHAST");
      }
      break;
    case 12:			/* debug_diffuse_layer */
      debug_diffuse_layer = get_true_false (next_char, TRUE);
      break;
    case 13:			/* delay_mass_water */
      delay_mass_water = get_true_false (next_char, TRUE);
      break;
    case 14:			/* convergence_tolerance */
      sscanf (next_char, SCANFORMAT, &convergence_tolerance);
      break;
    case 15:			/* numerical_derivatives */
      numerical_deriv = get_true_false (next_char, TRUE);
      break;
    }
    if (return_value == EOF || return_value == KEYWORD)
      break;
  }
  return (return_value);
}

/* ---------------------------------------------------------------------- */
int
read_print (void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads printing info
 *
 *      Arguments:
 *         none
 *
 *      Returns:
 *         KEYWORD if keyword encountered, input_error may be incremented if
 *                    a keyword is encountered in an unexpected position
 *         EOF     if eof encountered while reading mass balance concentrations
 *         ERROR   if error occurred reading data
 *
 */
  int return_value, opt, l;
  char *next_char;
  char token[MAX_LENGTH];
  const char *opt_list[] = {
    "reset",			/* 0 */
    "gas_phase",		/* 1 */
    "pure_phase",		/* 2 */
    "surface",			/* 3 */
    "exchange",			/* 4 */
    "totals",			/* 5 */
    "eh",			/* 6 */
    "species",			/* 7 */
    "saturation_indices",	/* 8 */
    "si",			/* 9 same a 8 */
    "reaction",			/* 10 irrev, not used, pr.use */
    "mix",			/* 11 */
    "use",			/* 12 */
    "selected_output",		/* 13 */
    "equilibrium_phases",	/* 14 same as 2 */
    "equilibria",		/* 15 same as 2 */
    "equilibrium",		/* 16 same as 2 */
    "pure",			/* 17 same as 2 */
    "other",			/* 18 same as 12 */
    "status",			/* 19 */
    "inverse",			/* 20 */
    "kinetics",			/* 21 */
    "dump",			/* 22 */
    "user_print",		/* 23 */
    "user_pr",			/* 24 */
    "solid_solution",		/* 25 */
    "solid_solutions",		/* 26 */
    "inverse_modeling",		/* 27 */
    "headings",			/* 28 */
    "heading",			/* 29 */
    "user_graph",		/* 30 */
    "echo_input",		/* 31 */
    "warning",			/* 32 */
    "warnings",			/* 33 */
    "initial_isotopes",		/* 34 */
    "isotope_ratios",		/* 35 */
    "isotope_alphas",		/* 36 */
    "censor_species",		/* 37 */
    "alkalinity"		/* 38 */
  };

  int count_opt_list = 39;
  int value;
/*
 *   Read flags:
 */
  return_value = UNKNOWN;
  for (;;)
  {
    opt = get_option (opt_list, count_opt_list, &next_char);
    switch (opt)
    {
    case OPTION_EOF:		/* end of file */
      return_value = EOF;
      break;
    case OPTION_KEYWORD:	/* keyword */
      return_value = KEYWORD;
      break;
    case OPTION_DEFAULT:
    case OPTION_ERROR:
      input_error++;
      error_msg ("Unknown input in PRINT keyword.", CONTINUE);
      error_msg (line_save, CONTINUE);
      break;
    case 0:			/* reset */
      value = get_true_false (next_char, TRUE);
      pr.kinetics = value;
      pr.gas_phase = value;
      pr.pp_assemblage = value;
      pr.surface = value;
      pr.exchange = value;
      pr.totals = value;
      pr.eh = value;
      pr.species = value;
      pr.saturation_indices = value;
      pr.irrev = value;
      pr.mix = value;
      pr.reaction = value;
      pr.use = value;
      pr.inverse = value;
      pr.user_print = value;
      pr.s_s_assemblage = value;
      pr.headings = value;
      pr.initial_isotopes = value;
      pr.isotope_ratios = value;
      pr.isotope_alphas = value;
      pr.echo_input = value;
      break;
    case 1:			/* gas_phase */
      pr.gas_phase = get_true_false (next_char, TRUE);
      break;
    case 2:			/* pure_phase */
    case 14:			/* equilibrium_phases */
    case 15:			/* equilibria */
    case 16:			/* equilibrium */
    case 17:			/* pure */
      pr.pp_assemblage = get_true_false (next_char, TRUE);
      break;
    case 3:			/* surface */
      pr.surface = get_true_false (next_char, TRUE);
      break;
    case 4:			/* exchange */
      pr.exchange = get_true_false (next_char, TRUE);
      break;
    case 5:			/* totals */
      pr.totals = get_true_false (next_char, TRUE);
      break;
    case 6:			/* eh */
      pr.eh = get_true_false (next_char, TRUE);
      break;
    case 7:			/* species */
      pr.species = get_true_false (next_char, TRUE);
      break;
    case 8:
    case 9:			/* saturation_indices */
      pr.saturation_indices = get_true_false (next_char, TRUE);
      break;
    case 10:			/* reaction, not used, pr.use controls */
      pr.irrev = get_true_false (next_char, TRUE);
      break;
    case 11:			/* mix, not used, pr.use controls */
      pr.mix = get_true_false (next_char, TRUE);
      break;
    case 12:			/* use */
    case 18:			/* other */
      pr.use = get_true_false (next_char, TRUE);
      break;
    case 13:			/* selected_output */
      pr.punch = get_true_false (next_char, TRUE);
      break;
    case 19:			/* status */
      pr.status = get_true_false (next_char, TRUE);
      break;
    case 20:			/* inverse */
    case 27:			/* inverse_modeling */
      pr.inverse = get_true_false (next_char, TRUE);
      break;
    case 21:			/* kinetics */
      pr.kinetics = get_true_false (next_char, TRUE);
      break;
    case 22:			/* dump */
      pr.dump = get_true_false (next_char, TRUE);
      break;
    case 23:			/* user_print */
    case 24:			/* user_pr */
      pr.user_print = get_true_false (next_char, TRUE);
      break;
    case 25:			/* solid_solution */
    case 26:			/* solid_solutions */
      pr.s_s_assemblage = get_true_false (next_char, TRUE);
      break;
    case 28:			/* headings */
    case 29:			/* heading */
      pr.headings = get_true_false (next_char, TRUE);
      break;
    case 30:			/* user_graph */
      pr.user_graph = get_true_false (next_char, TRUE);
      break;
    case 31:			/* echo_input */
      pr.echo_input = get_true_false (next_char, TRUE);
      break;
    case 32:			/* warning */
    case 33:			/* warnings */
      sscanf (next_char, "%d", &pr.warnings);
      break;
    case 34:			/* initial_isotopes */
      pr.initial_isotopes = get_true_false (next_char, TRUE);
      break;
    case 35:			/* isotope_ratios */
      pr.isotope_ratios = get_true_false (next_char, TRUE);
      break;
    case 36:			/* isotope_alphas */
      pr.isotope_alphas = get_true_false (next_char, TRUE);
      break;
    case 37:			/* censor_species */
      if (copy_token (token, &next_char, &l) != EMPTY)
      {
	sscanf (token, SCANFORMAT, &censor);
      }
      else
      {
	censor = 0;
      }
      break;
    case 38:			/* alkalinity */
      pr.alkalinity = get_true_false (next_char, TRUE);
      break;
    }
    if (return_value == EOF || return_value == KEYWORD)
      break;
  }
  return (return_value);
}

/* ---------------------------------------------------------------------- */
/*
 *   Utilitity routines for read
 */
/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */
int
check_key (char *str)
/* ---------------------------------------------------------------------- */
{
/*
 *   Check if string begins with a keyword, returns TRUE or FALSE
 *
 *   Arguments:
 *      *str is pointer to character string to be checked for keywords
 *   Returns:
 *      TRUE,
 *      FALSE.
 */
  int l;
  char *ptr;
  char token[MAX_LENGTH];
  ENTRY item, *found_item;

  ptr = str;
  copy_token (token, &ptr, &l);
  str_tolower (token);
  next_keyword = -1;
/*
 *   Search hash list
 */
  item.key = token;
  item.data = NULL;
  found_item = hsearch_multi (keyword_hash_table, item, FIND);
  if (found_item == NULL)
    return (FALSE);
  next_keyword = ((struct key *) (found_item->data))->keycount;
  return (TRUE);
}

/* ---------------------------------------------------------------------- */
int
check_units (char *tot_units, int alkalinity, int check_compatibility,
	     const char *default_units, int print)
/* ---------------------------------------------------------------------- */
{
#define NUNITS (sizeof(units) / sizeof(char *))
/*
 *   Check if legitimate units
 *   Input:
 *           tot_units           character string to check,
 *           alkalinity          TRUE if alkalinity, FALSE if any other total,
 *           check_compatibility TRUE check alk and default units, FALSE otherwise
 *           default_units       character string of default units (check /L, /kg, etc)
 *           print               TRUE print warning messages
 *   Output:
 *           tot_units           standard form for unit
 */
  int i, found;
  char *end;
  char string[MAX_LENGTH];
  const char *units[] = {
    "Mol/l",			/* 0 */
    "mMol/l",			/* 1 */
    "uMol/l",			/* 2 */
    "g/l",			/* 3 */
    "mg/l",			/* 4 */
    "ug/l",			/* 5 */
    "Mol/kgs",			/* 6 */
    "mMol/kgs",			/* 7 */
    "uMol/kgs",			/* 8 */
    "g/kgs",			/* 9 = ppt */
    "mg/kgs",			/* 10 = ppm */
    "ug/kgs",			/* 11 = ppb */
    "Mol/kgw",			/* 12 = mol/kg H2O */
    "mMol/kgw",			/* 13 = mmol/kg H2O */
    "uMol/kgw",			/* 14 = umol/kg H2O */
    "g/kgw",			/* 15 = mol/kg H2O */
    "mg/kgw",			/* 16 = mmol/kg H2O */
    "ug/kgw",			/* 17 = umol/kg H2O */
    "eq/l",			/* 18 */
    "meq/l",			/* 19 */
    "ueq/l",			/* 20 */
    "eq/kgs",			/* 21 */
    "meq/kgs",			/* 22 */
    "ueq/kgs",			/* 23 */
    "eq/kgw",			/* 24 */
    "meq/kgw",			/* 25 */
    "ueq/kgw",			/* 26 */
  };

  squeeze_white (tot_units);
  str_tolower (tot_units);
  replace ("milli", "m", tot_units);
  replace ("micro", "u", tot_units);
  replace ("grams", "g", tot_units);
  replace ("gram", "g", tot_units);
  replace ("moles", "Mol", tot_units);
  replace ("mole", "Mol", tot_units);
  replace ("mol", "Mol", tot_units);
  replace ("liter", "l", tot_units);
  replace ("kgh", "kgw", tot_units);
  replace ("ppt", "g/kgs", tot_units);
  replace ("ppm", "mg/kgs", tot_units);
  replace ("ppb", "ug/kgs", tot_units);
  replace ("equivalents", "eq", tot_units);
  replace ("equivalent", "eq", tot_units);
  replace ("equiv", "eq", tot_units);

  if ((end = strstr (tot_units, "/l")) != NULL)
  {
    *(end + 2) = '\0';
  }
  if ((end = strstr (tot_units, "/kgs")) != NULL)
  {
    *(end + 4) = '\0';
  }
  if ((end = strstr (tot_units, "/kgw")) != NULL)
  {
    *(end + 4) = '\0';
  }
/*
 *   Check if unit in list
 */
  found = FALSE;
  for (i = 0; i < (int) NUNITS; i++)
  {
    if (strcmp (tot_units, units[i]) == 0)
    {
      found = TRUE;
      break;
    }
  }
  if (found == FALSE)
  {
    if (print == TRUE)
    {
      sprintf (error_string, "Unknown unit, %s.", tot_units);
      error_msg (error_string, CONTINUE);
    }
    return (ERROR);
  }

/*
 *   Check if units are compatible with default_units
 */
  if (check_compatibility == FALSE)
    return (OK);
/*
 *   Special cases for alkalinity
 */
  if (alkalinity == TRUE && strstr (tot_units, "Mol") != NULL)
  {
    if (print == TRUE)
    {
      sprintf (error_string,
	       "Alkalinity given in moles, assumed to be equivalents.");
      warning_msg (error_string);
    }
    replace ("Mol", "eq", tot_units);
  }
  if (alkalinity == FALSE && strstr (tot_units, "eq") != NULL)
  {
    if (print == TRUE)
    {
      error_msg ("Only alkalinity can be entered in equivalents.", CONTINUE);
    }
    return (ERROR);
  }
/*
 *   See if default_units are compatible with tot_units
 */
  if (strstr (default_units, "/l") && strstr (tot_units, "/l"))
    return (OK);
  if (strstr (default_units, "/kgs") && strstr (tot_units, "/kgs"))
    return (OK);
  if (strstr (default_units, "/kgw") && strstr (tot_units, "/kgw"))
    return (OK);

  strcpy (string, default_units);
  replace ("kgs", "kg solution", string);
  replace ("kgs", "kg solution", tot_units);
  replace ("kgw", "kg water", string);
  replace ("kgw", "kg water", tot_units);
  replace ("/l", "/L", string);
  replace ("Mol", "mol", string);
  replace ("/l", "/L", tot_units);
  replace ("Mol", "mol", tot_units);
  if (print == TRUE)
  {
    sprintf (error_string,
	     "Units for master species, %s, are not compatible with default units, %s.",
	     tot_units, string);
    error_msg (error_string, CONTINUE);
  }
  return (ERROR);
}

/* ---------------------------------------------------------------------- */
int
find_option (char *item, int *n, const char **list, int count_list, int exact)
/* ---------------------------------------------------------------------- */
{
/*
 *   Compares a string value to match beginning letters of a list of options
 *
 *      Arguments:
 *         item    entry: pointer to string to compare
 *         n       exit:  item in list that was matched
 *         list    entry: pointer to list of character values, assumed to
 *                 be lower case
 *         count_list entry: number of character values in list
 *
 *      Returns:
 *         OK      item matched
 *         ERROR   item not matched
 *         n       -1      item not matched
 *                 i       position of match in list
 */
  int i;
  char token[MAX_LENGTH];

  strcpy (token, item);
  str_tolower (token);
  for (i = 0; i < count_list; i++)
  {
    if (exact == TRUE)
    {
      if (strcmp (list[i], token) == 0)
      {
	*n = i;
	return (OK);
      }
    }
    else
    {
      if (strstr (list[i], token) == list[i])
      {
	*n = i;
	return (OK);
      }
    }
  }
  *n = -1;
  return (ERROR);
}

/* ---------------------------------------------------------------------- */
int
get_true_false (char *string, int default_value)
/* ---------------------------------------------------------------------- */
{
/*
 *   Returns true unless string starts with "F" or "f"
 */
  int l;
  char token[MAX_LENGTH];
  char *ptr;

  ptr = string;

  if (copy_token (token, &ptr, &l) == EMPTY)
  {
    return (default_value);
  }
  else
  {
    if (token[0] == 'F' || token[0] == 'f')
    {
      return (FALSE);
    }
  }
  return (TRUE);
}

/* ---------------------------------------------------------------------- */
int
get_option (const char **opt_list, int count_opt_list, char **next_char)
/* ---------------------------------------------------------------------- */
{
/*
 *   Read a line and check for options
 */
  int j;
  int opt_l, opt;
  char *opt_ptr;
  char option[MAX_LENGTH];
/*
 *   Read line
 */
  j = check_line ("get_option", FALSE, TRUE, TRUE, FALSE);
  if (j == EOF)
  {
    j = OPTION_EOF;
  }
  else if (j == KEYWORD)
  {
    j = OPTION_KEYWORD;
  }
  else if (j == OPTION)
  {
    opt_ptr = line;
    copy_token (option, &opt_ptr, &opt_l);
    if (find_option (&(option[1]), &opt, opt_list, count_opt_list, FALSE) ==
	OK)
    {
      j = opt;
      replace (option, opt_list[j], line_save);
      replace (option, opt_list[j], line);
      opt_ptr = line;
      copy_token (option, &opt_ptr, &opt_l);
      *next_char = opt_ptr;
      if (pr.echo_input == TRUE)
      {
	if (!reading_database ())
	  output_msg (OUTPUT_MESSAGE, "\t%s\n", line_save);
      }
    }
    else
    {
      if (!reading_database ())
	output_msg (OUTPUT_MESSAGE, "\t%s\n", line_save);
      error_msg ("Unknown option.", CONTINUE);
      error_msg (line_save, CONTINUE);
      input_error++;
      j = OPTION_ERROR;
      *next_char = line;
    }
  }
  else
  {
    opt_ptr = line;
    copy_token (option, &opt_ptr, &opt_l);
    if (find_option (&(option[0]), &opt, opt_list, count_opt_list, TRUE) ==
	OK)
    {
      j = opt;
      *next_char = opt_ptr;
    }
    else
    {
      j = OPTION_DEFAULT;
      *next_char = line;
    }
    if (pr.echo_input == TRUE)
    {
      if (!reading_database ())
	output_msg (OUTPUT_MESSAGE, "\t%s\n", line_save);
    }
  }
  return (j);
}

/* ---------------------------------------------------------------------- */
int
read_rates (void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads basic code with which to calculate rates
 *
 *      Arguments:
 *         none
 *
 *      Returns:
 *         KEYWORD if keyword encountered, input_error may be incremented if
 *                    a keyword is encountered in an unexpected position
 *         EOF     if eof encountered while reading mass balance concentrations
 *         ERROR   if error occurred reading data
 *
 */
  char *ptr;
  int l, length, line_length, n;
  int return_value, opt, opt_save;
  char token[MAX_LENGTH];
  struct rate *rate_ptr;
  char *description;
  int n_user, n_user_end;
  char *next_char;
  const char *opt_list[] = {
    "start",			/* 0 */
    "end"			/* 1 */
  };
  int count_opt_list = 2;
/*
 *   Read advection number (not currently used)
 */
  n = -1;
  ptr = line;
  read_number_description (ptr, &n_user, &n_user_end, &description);
  description = (char *) free_check_null (description);
  opt_save = OPTION_DEFAULT;
/*
 *   Read lines
 */
  return_value = UNKNOWN;
  rate_ptr = NULL;
  for (;;)
  {
    opt = get_option (opt_list, count_opt_list, &next_char);
    if (opt == OPTION_DEFAULT)
    {
      opt = opt_save;
    }
    switch (opt)
    {
    case OPTION_EOF:		/* end of file */
      return_value = EOF;
      break;
    case OPTION_KEYWORD:	/* keyword */
      return_value = KEYWORD;
      break;
    case OPTION_ERROR:
      input_error++;
      error_msg ("Unknown input in RATES keyword.", CONTINUE);
      error_msg (line_save, CONTINUE);
      break;
    case 0:			/* start */
      opt_save = OPT_1;
      break;
    case 1:			/* end */
      opt_save = OPTION_DEFAULT;
      break;
    case OPTION_DEFAULT:	/* read rate name */
      ptr = line;
      copy_token (token, &ptr, &l);
      rate_ptr = rate_search (token, &n);
      if (rate_ptr == NULL)
      {
	rates =
	  (struct rate *) PHRQ_realloc (rates,
					(size_t) (count_rates +
						  1) * sizeof (struct rate));
	if (rates == NULL)
	  malloc_error ();
	rate_ptr = &rates[count_rates];
	count_rates++;
      }
      else
      {
	rate_free (rate_ptr);
      }
      rate_ptr->new_def = TRUE;
      rate_ptr->commands = (char *) PHRQ_malloc (sizeof (char));
      if (rate_ptr->commands == NULL)
	malloc_error ();
      rate_ptr->commands[0] = '\0';
      rate_ptr->name = string_hsave (token);
      rate_ptr->linebase = NULL;
      rate_ptr->varbase = NULL;
      rate_ptr->loopbase = NULL;
      opt_save = OPT_1;
      break;
    case OPT_1:		/* read command */
      if (rate_ptr == NULL)
      {
	input_error++;
	sprintf (error_string, "No rate name has been defined.");
	error_msg (error_string, CONTINUE);
	opt_save = OPT_1;
	break;
      }
      length = (int) strlen (rate_ptr->commands);
      line_length = (int) strlen (line);
      rate_ptr->commands =
	(char *) PHRQ_realloc (rate_ptr->commands,
			       (size_t) (length + line_length +
					 2) * sizeof (char));
      if (rate_ptr->commands == NULL)
	malloc_error ();
      rate_ptr->commands[length] = ';';
      rate_ptr->commands[length + 1] = '\0';
      strcat ((rate_ptr->commands), line);
      opt_save = OPT_1;
      break;
    }
    if (return_value == EOF || return_value == KEYWORD)
      break;
  }
/*	output_msg(OUTPUT_MESSAGE, "%s", rates[0].commands);
 */ return (return_value);
}

/* ---------------------------------------------------------------------- */
int
read_user_print (void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads basic code with which to calculate rates
 *
 *      Arguments:
 *         none
 *
 *      Returns:
 *         KEYWORD if keyword encountered, input_error may be incremented if
 *                    a keyword is encountered in an unexpected position
 *         EOF     if eof encountered while reading mass balance concentrations
 *         ERROR   if error occurred reading data
 *
 */
  int length, line_length;
  int return_value, opt, opt_save;
  char *next_char;
  const char *opt_list[] = {
    "start",			/* 0 */
    "end"			/* 1 */
  };
  int count_opt_list = 2;

  opt_save = OPTION_DEFAULT;
/*
 *   Read lines
 */
  return_value = UNKNOWN;
  for (;;)
  {
    opt = get_option (opt_list, count_opt_list, &next_char);
    if (opt == OPTION_DEFAULT)
    {
      opt = opt_save;
    }
    switch (opt)
    {
    case OPTION_EOF:		/* end of file */
      return_value = EOF;
      break;
    case OPTION_KEYWORD:	/* keyword */
      return_value = KEYWORD;
      break;
    case OPTION_ERROR:
      input_error++;
      error_msg ("Unknown input in RATES keyword.", CONTINUE);
      error_msg (line_save, CONTINUE);
      break;
    case 0:			/* start */
      opt_save = OPTION_DEFAULT;
      break;
    case 1:			/* end */
      opt_save = OPTION_DEFAULT;
      break;
    case OPTION_DEFAULT:	/* read first command */
      rate_free (user_print);
      user_print->new_def = TRUE;
      user_print->commands = (char *) PHRQ_malloc (sizeof (char));
      if (user_print->commands == NULL)
	malloc_error ();
      user_print->commands[0] = '\0';
      user_print->linebase = NULL;
      user_print->varbase = NULL;
      user_print->loopbase = NULL;
      user_print->name = string_hsave ("user defined Basic print routine");
    case OPT_1:		/* read command */
      length = (int) strlen (user_print->commands);
      line_length = (int) strlen (line);
      user_print->commands =
	(char *) PHRQ_realloc (user_print->commands,
			       (size_t) (length + line_length +
					 2) * sizeof (char));
      if (user_print->commands == NULL)
	malloc_error ();
      user_print->commands[length] = ';';
      user_print->commands[length + 1] = '\0';
      strcat ((user_print->commands), line);
      opt_save = OPT_1;
      break;
    }
    if (return_value == EOF || return_value == KEYWORD)
      break;
  }
/*	output_msg(OUTPUT_MESSAGE, "%s", rates[0].commands);
 */ return (return_value);
}

/* ---------------------------------------------------------------------- */
int
read_user_punch (void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads basic code with which to calculate rates
 *
 *      Arguments:
 *         none
 *
 *      Returns:
 *         KEYWORD if keyword encountered, input_error may be incremented if
 *                    a keyword is encountered in an unexpected position
 *         EOF     if eof encountered while reading mass balance concentrations
 *         ERROR   if error occurred reading data
 *
 */
  int l, length, line_length;
  int return_value, opt, opt_save;
  char token[MAX_LENGTH];
  char *next_char;
  const char *opt_list[] = {
    "start",			/* 0 */
    "end",			/* 1 */
    "heading",			/* 2 */
    "headings"			/* 3 */
  };
  int count_opt_list = 4;

  opt_save = OPTION_DEFAULT;
/*
 *   Read lines
 */
  user_punch_count_headings = 0;
  return_value = UNKNOWN;
  for (;;)
  {
    opt = get_option (opt_list, count_opt_list, &next_char);
    if (opt == OPTION_DEFAULT)
    {
      opt = opt_save;
    }
    switch (opt)
    {
    case OPTION_EOF:		/* end of file */
      return_value = EOF;
      break;
    case OPTION_KEYWORD:	/* keyword */
      return_value = KEYWORD;
      break;
    case OPTION_ERROR:
      input_error++;
      error_msg ("Unknown input in RATES keyword.", CONTINUE);
      error_msg (line_save, CONTINUE);
      break;
    case 0:			/* start */
      opt_save = OPTION_DEFAULT;
      break;
    case 1:			/* end */
      opt_save = OPTION_DEFAULT;
      break;
    case 2:			/* headings */
    case 3:			/* heading */
      while (copy_token (token, &next_char, &l) != EMPTY)
      {
	user_punch_headings =
	  (char **) PHRQ_realloc (user_punch_headings,
				  (size_t) (user_punch_count_headings +
					    1) * sizeof (char *));
	if (user_punch_headings == NULL)
	  malloc_error ();
	user_punch_headings[user_punch_count_headings] = string_hsave (token);
	user_punch_count_headings++;
      }
      break;
    case OPTION_DEFAULT:	/* read first command */
      rate_free (user_punch);
      user_punch->new_def = TRUE;
      user_punch->commands = (char *) PHRQ_malloc (sizeof (char));
      if (user_punch->commands == NULL)
	malloc_error ();
      user_punch->commands[0] = '\0';
      user_punch->linebase = NULL;
      user_punch->varbase = NULL;
      user_punch->loopbase = NULL;
      user_punch->name = string_hsave ("user defined Basic punch routine");
    case OPT_1:		/* read command */
      length = (int) strlen (user_punch->commands);
      line_length = (int) strlen (line);
      user_punch->commands =
	(char *) PHRQ_realloc (user_punch->commands,
			       (size_t) (length + line_length +
					 2) * sizeof (char));
      if (user_punch->commands == NULL)
	malloc_error ();
      user_punch->commands[length] = ';';
      user_punch->commands[length + 1] = '\0';
      strcat ((user_punch->commands), line);
      opt_save = OPT_1;
      break;
    }
    if (return_value == EOF || return_value == KEYWORD)
      break;
  }
  return (return_value);
}

#ifdef PHREEQ98
/* ---------------------------------------------------------------------- */
int
read_user_graph (void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads basic code with which to calculate rates
 *
 *      Arguments:
 *         none
 *
 *      Returns:
 *         KEYWORD if keyword encountered, input_error may be incremented if
 *                    a keyword is encountered in an unexpected position
 *         EOF     if eof encountered while reading mass balance concentrations
 *         ERROR   if error occurred reading data
 *
 */
  int l, length, line_length;
  int return_value, opt, opt_save;
  char file_name[MAX_LENGTH];
  char token[MAX_LENGTH];
  char *next_char;
  const char *opt_list[] = {
    "start",			/* 0 */
    "end",			/* 1 */
    "heading",			/* 2 */
    "headings",			/* 3 */
    "chart_title",		/* 4 */
    "axis_titles",		/* 5 */
    "axis_scale",		/* 6 */
    "initial_solutions",	/* 7 */
    "plot_concentration_vs",	/* 8 */
    "shifts_as_points",		/* 9 */
    "grid_offset",		/* 10 */
    "connect_simulations",	/* 11 */
    "plot_csv_file"		/* 12 */
  };
  int count_opt_list = 13;
  int i;

  opt_save = OPTION_DEFAULT;
/*
 *   Read lines
 */
  user_graph_count_headings = 0;
  return_value = UNKNOWN;
  for (;;)
  {
    opt = get_option (opt_list, count_opt_list, &next_char);
    if (opt == OPTION_DEFAULT)
    {
      opt = opt_save;
    }
    switch (opt)
    {
    case OPTION_EOF:		/* end of file */
      return_value = EOF;
      break;
    case OPTION_KEYWORD:	/* keyword */
      return_value = KEYWORD;
      break;
    case OPTION_ERROR:
      input_error++;
      error_msg ("Unknown input in USER_GRAPH keyword.", CONTINUE);
      error_msg (line_save, CONTINUE);
      break;
    case 0:			/* start */
      opt_save = OPTION_DEFAULT;
      break;
    case 1:			/* end */
      opt_save = OPTION_DEFAULT;
      break;
    case 2:			/* headings */
    case 3:			/* heading */
      while (copy_token (token, &next_char, &l) != EMPTY)
      {
	user_graph_headings =
	  PHRQ_realloc (user_graph_headings,
			(size_t) (user_graph_count_headings +
				  1) * sizeof (char *));
	if (user_graph_headings == NULL)
	  malloc_error ();
	user_graph_headings[user_graph_count_headings] = string_hsave (token);
	user_graph_count_headings++;
      }
      break;
/*Modifications of read_user_punch to change the chart's appearance */
    case 4:			/* chart title */
      copy_title (token, &next_char, &l);
      SetChartTitle (token);
      break;
    case 5:
      {				/* axis titles */
	i = 0;
	while (copy_title (token, &next_char, &l) != EMPTY)
	{
	  SetAxisTitles (token, i);
	  i++;
	}
      }
      break;
    case 6:
      {				/* axis scales */
	char *axis = "";
	int j = 0;
	copy_token (token, &next_char, &l);
	str_tolower (token);
	if (strstr (token, "x") == token)
	{
	  axis = "x";
	}
	else if (strstr (token, "y") == token)
	{
	  axis = "y";
	}
	else if (strstr (token, "s") == token)
	{
	  axis = "s";
	}
	else
	{
	  input_error++;
	  error_msg ("Expected axis type.", CONTINUE);
	}
	while ((j < 4) && (i = copy_token (token, &next_char, &l)) != EMPTY)
	{
	  str_tolower (token);
	  if ((i == DIGIT) || (strstr (token, "auto") == token))
	  {
	    SetAxisScale (axis, j, token, FALSE);
	  }
	  else
	  {
	    input_error++;
	    error_msg ("Expected numerical value or 'auto'.", CONTINUE);
	  }
	  j++;			/* counter for categories */
	}
	if (j == 4)
	  SetAxisScale (axis, j, 0, get_true_false (next_char, FALSE));
      }
      break;
    case 7:
      graph_initial_solutions = get_true_false (next_char, TRUE);
      break;
    case 8:
      copy_token (token, &next_char, &l);
      str_tolower (token);
      if (strstr (token, "x") == token)
      {
	chart_type = 0;
      }
      else if (strstr (token, "t") == token)
      {
	chart_type = 1;
      }
      else
      {
	input_error++;
	error_msg ("Expected simulation type.", CONTINUE);
      }
      break;
    case 9:
      shifts_as_points = get_true_false (next_char, TRUE);
      if (shifts_as_points == TRUE)
	chart_type = 0;
      else
	chart_type = 1;
      break;
    case 10:
      {
	i = copy_token (token, &next_char, &l);
	str_tolower (token);
	if (i == DIGIT)
	  sscanf (token, "%d", &RowOffset);
	i = copy_token (token, &next_char, &l);
	str_tolower (token);
	if (i == DIGIT)
	  sscanf (token, "%d", &ColumnOffset);
      }
      break;
    case 11:
      connect_simulations = get_true_false (next_char, TRUE);
      break;
    case 12:
      string_trim (next_char);
      strcpy (file_name, next_char);
      if (!OpenCSVFile (file_name))
      {
	sprintf (error_string, "Can't open file, %s.", file_name);
	input_error++;
	error_msg (error_string, CONTINUE);
      }
      break;
      /* End of modifications */
    case OPTION_DEFAULT:	/* read first command */
      rate_free (user_graph);
      user_graph->new_def = TRUE;
      user_graph->commands = PHRQ_malloc (sizeof (char));
      if (user_graph->commands == NULL)
	malloc_error ();
      user_graph->commands[0] = '\0';
      user_graph->linebase = NULL;
      user_graph->varbase = NULL;
      user_graph->loopbase = NULL;
      user_graph->name = string_hsave ("user defined Basic punch routine");
    case OPT_1:		/* read command */
      length = strlen (user_graph->commands);
      line_length = strlen (line);
      user_graph->commands =
	PHRQ_realloc (user_graph->commands,
		      (size_t) (length + line_length + 2) * sizeof (char));
      if (user_graph->commands == NULL)
	malloc_error ();
      user_graph->commands[length] = ';';
      user_graph->commands[length + 1] = '\0';
      strcat ((user_graph->commands), line);
      opt_save = OPT_1;
      break;
    }
    if (return_value == EOF || return_value == KEYWORD)
      break;
  }
  for (i = 0; i < user_graph_count_headings; i++)
    GridHeadings (user_graph_headings[i], i);
  return (return_value);
}
#endif
/* ---------------------------------------------------------------------- */
int
read_solid_solutions (void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads solid solution data
 *
 *      Arguments:
 *         none
 *
 *      Returns:
 *         KEYWORD if keyword encountered, input_error may be incremented if
 *                    a keyword is encountered in an unexpected position
 *         EOF     if eof encountered while reading mass balance concentrations
 *         ERROR   if error occurred reading data
 *
 */
  int i, j, n, l;
  int count_s_s, number_s_s, count_comps;
  int n_user, n_user_end;
  char *ptr;
  char *description;
  char token[MAX_LENGTH];

  int return_value, opt;
  char *next_char;
  const char *opt_list[] = {
    "component",		/* 0 */
    "comp",			/* 1 */
    "parms",			/* 2 */
    "gugg_nondimensional",	/* 3 */
    "gugg_kj",			/* 4 */
    "activity_coefficients",	/* 5 */
    "distribution_coefficients",	/* 6 */
    "miscibility_gap",		/* 7 */
    "spinodal_gap",		/* 8 */
    "critical_point",		/* 9 */
    "alyotropic_point",		/* 10 */
    "temp",			/* 11 */
    "tempk",			/* 12 */
    "tempc",			/* 13 */
    "thompson",			/* 14 */
    "margules",			/* 15 */
    "comp1",			/* 16 */
    "comp2"			/* 17 */
  };
  int count_opt_list = 18;
/*
 *   Read s_s_assemblage number
 */
  number_s_s = 0;
  ptr = line;
  read_number_description (ptr, &n_user, &n_user_end, &description);
/*
 *   Find old s_s_assemblage or alloc space for new s_s_assemblage
 */
  if (s_s_assemblage_search (n_user, &n) != NULL)
  {
    s_s_assemblage_free (&s_s_assemblage[n]);
  }
  else
  {
    n = count_s_s_assemblage;
    space ((void **) ((void *) &s_s_assemblage), count_s_s_assemblage,
	   &max_s_s_assemblage, sizeof (struct s_s_assemblage));
    count_s_s_assemblage++;
  }
/*
 *   Initialize
 */
  s_s_assemblage_init (&(s_s_assemblage[n]), n_user, n_user_end, description);
  free_check_null (description);
/*
 *   Set use data to first read
 */
  if (use.s_s_assemblage_in == FALSE)
  {
    use.s_s_assemblage_in = TRUE;
    use.n_s_s_assemblage_user = n_user;
  }
/*
 *   Read solid solutions
 */
  count_s_s = 0;
  count_comps = 0;
  return_value = UNKNOWN;
  for (;;)
  {
    opt = get_option (opt_list, count_opt_list, &next_char);
    switch (opt)
    {
    case OPTION_EOF:		/* end of file */
      return_value = EOF;
      break;
    case OPTION_KEYWORD:	/* keyword */
      return_value = KEYWORD;
      break;
    case OPTION_ERROR:
      input_error++;
      error_msg ("Unknown input in SOLID_SOLUTIONS keyword.", CONTINUE);
      error_msg (line_save, CONTINUE);
      break;

/* 
 * New component
 */
    case 0:			/* component */
    case 1:			/* comp */
      count_comps = s_s_assemblage[n].s_s[number_s_s].count_comps++;
      s_s_assemblage[n].s_s[number_s_s].comps =
	(struct s_s_comp *) PHRQ_realloc (s_s_assemblage[n].s_s[number_s_s].
					  comps,
					  (size_t) (count_comps +
						    1) *
					  sizeof (struct s_s_comp));
      if (s_s_assemblage[n].s_s[number_s_s].comps == NULL)
	malloc_error ();
      s_s_assemblage[n].s_s[number_s_s].comps[count_comps].initial_moles = 0;
      s_s_assemblage[n].s_s[number_s_s].comps[count_comps].delta = 0;
      /*
       *   Read phase name of component
       */
      ptr = next_char;
      copy_token (token, &ptr, &l);
      s_s_assemblage[n].s_s[number_s_s].comps[count_comps].name =
	string_hsave (token);
      /*
       *   Read moles of component
       */
      if ((j = copy_token (token, &ptr, &l)) == EMPTY)
      {
	s_s_assemblage[n].s_s[number_s_s].comps[count_comps].moles = NAN;
      }
      else
      {
	j = sscanf (token, SCANFORMAT, &dummy);
	s_s_assemblage[n].s_s[number_s_s].comps[count_comps].moles =
	  (LDBLE) dummy;
	if (j != 1)
	{
	  error_msg ("Expected moles of solid solution.", CONTINUE);
	  error_msg (line_save, CONTINUE);
	  input_error++;
	}
      }
      break;
    case 2:			/* parms */
    case 3:			/* gugg_nondimensional */
      /*
       *   Read parameters
       */

      ptr = next_char;
      if (copy_token (token, &ptr, &l) != EMPTY)
      {
	sscanf (token, SCANFORMAT, &(s_s_assemblage[n].s_s[number_s_s].p[0]));
      }
      if (copy_token (token, &ptr, &l) != EMPTY)
      {
	sscanf (token, SCANFORMAT, &(s_s_assemblage[n].s_s[number_s_s].p[1]));
      }
      s_s_assemblage[n].s_s[number_s_s].input_case = 0;
      break;
    case 4:			/* gugg_kj */
      ptr = next_char;
      if (copy_token (token, &ptr, &l) != EMPTY)
      {
	sscanf (token, SCANFORMAT, &(s_s_assemblage[n].s_s[number_s_s].p[0]));
      }
      if (copy_token (token, &ptr, &l) != EMPTY)
      {
	sscanf (token, SCANFORMAT, &(s_s_assemblage[n].s_s[number_s_s].p[1]));
      }
      s_s_assemblage[n].s_s[number_s_s].input_case = 7;
      break;
    case 5:			/* activity coefficients */
      ptr = next_char;
      j = 0;
      for (i = 0; i < 4; i++)
      {
	if (copy_token (token, &ptr, &l) != EMPTY)
	{
	  j +=
	    sscanf (token, SCANFORMAT,
		    &(s_s_assemblage[n].s_s[number_s_s].p[i]));
	}
      }
      if (j != 4)
      {
	sprintf (error_string,
		 "Expected 4 parameters to calculate a0 and a1 from two activity coefficients, assemblage %d, solid solution %s",
		 s_s_assemblage[n].n_user,
		 s_s_assemblage[n].s_s[number_s_s].name);
	error_msg (error_string, CONTINUE);
	input_error++;
      }
      s_s_assemblage[n].s_s[number_s_s].input_case = 1;
      break;
    case 6:			/* distribution coefficients */
      ptr = next_char;
      j = 0;
      for (i = 0; i < 4; i++)
      {
	if (copy_token (token, &ptr, &l) != EMPTY)
	{
	  j +=
	    sscanf (token, SCANFORMAT,
		    &(s_s_assemblage[n].s_s[number_s_s].p[i]));
	}
      }
      if (j != 4)
      {
	sprintf (error_string,
		 "Expected 4 parameters to calculate a0 and a1 from two distribution coefficients, assemblage %d, solid solution %s",
		 s_s_assemblage[n].n_user,
		 s_s_assemblage[n].s_s[number_s_s].name);
	error_msg (error_string, CONTINUE);
	input_error++;
      }
      s_s_assemblage[n].s_s[number_s_s].input_case = 2;
      break;
    case 7:			/* miscibility_gap */
      ptr = next_char;
      j = 0;
      for (i = 0; i < 2; i++)
      {
	if (copy_token (token, &ptr, &l) != EMPTY)
	{
	  j +=
	    sscanf (token, SCANFORMAT,
		    &(s_s_assemblage[n].s_s[number_s_s].p[i]));
	}
      }
      if (j != 2)
      {
	sprintf (error_string,
		 "Expected 2 miscibility gap fractions of component 2 to calculate a0 and a1, assemblage %d, solid solution %s",
		 s_s_assemblage[n].n_user,
		 s_s_assemblage[n].s_s[number_s_s].name);
	error_msg (error_string, CONTINUE);
	input_error++;
      }
      s_s_assemblage[n].s_s[number_s_s].input_case = 3;
      break;
    case 8:			/* spinodal_gap */
      ptr = next_char;
      j = 0;
      for (i = 0; i < 2; i++)
      {
	if (copy_token (token, &ptr, &l) != EMPTY)
	{
	  j +=
	    sscanf (token, SCANFORMAT,
		    &(s_s_assemblage[n].s_s[number_s_s].p[i]));
	}
      }
      if (j != 2)
      {
	sprintf (error_string,
		 "Expected 2 spinodal gap fractions of component 2 to calculate a0 and a1, assemblage %d, solid solution %s",
		 s_s_assemblage[n].n_user,
		 s_s_assemblage[n].s_s[number_s_s].name);
	error_msg (error_string, CONTINUE);
	input_error++;
      }
      s_s_assemblage[n].s_s[number_s_s].input_case = 4;
      break;
    case 9:			/* critical point */
      ptr = next_char;
      j = 0;
      for (i = 0; i < 2; i++)
      {
	if (copy_token (token, &ptr, &l) != EMPTY)
	{
	  j +=
	    sscanf (token, SCANFORMAT,
		    &(s_s_assemblage[n].s_s[number_s_s].p[i]));
	}
      }
      if (j != 2)
      {
	sprintf (error_string,
		 "Expected fraction of component 2 and critical temperature to calculate a0 and a1, assemblage %d, solid solution %s",
		 s_s_assemblage[n].n_user,
		 s_s_assemblage[n].s_s[number_s_s].name);
	error_msg (error_string, CONTINUE);
	input_error++;
      }
      s_s_assemblage[n].s_s[number_s_s].input_case = 5;
      break;
    case 10:			/* alyotropic point */
      ptr = next_char;
      j = 0;
      for (i = 0; i < 2; i++)
      {
	if (copy_token (token, &ptr, &l) != EMPTY)
	{
	  j +=
	    sscanf (token, SCANFORMAT,
		    &(s_s_assemblage[n].s_s[number_s_s].p[i]));
	}
      }
      if (j != 2)
      {
	sprintf (error_string,
		 "Expected fraction of component 2 and sigma pi at alyotropic point to calculate a0 and a1, assemblage %d, solid solution %s",
		 s_s_assemblage[n].n_user,
		 s_s_assemblage[n].s_s[number_s_s].name);
	error_msg (error_string, CONTINUE);
	input_error++;
      }
      s_s_assemblage[n].s_s[number_s_s].input_case = 6;
      break;
    case 12:			/* tempk */
      ptr = next_char;
      j = 0;
      if (copy_token (token, &ptr, &l) != EMPTY)
      {
	j =
	  sscanf (token, SCANFORMAT, &(s_s_assemblage[n].s_s[number_s_s].tk));
      }
      if (j != 1)
      {
	sprintf (error_string,
		 "Expected temperature (Kelvin) for parameters, assemblage %d, solid solution %s, using 298.15 K",
		 s_s_assemblage[n].n_user,
		 s_s_assemblage[n].s_s[number_s_s].name);
	warning_msg (error_string);
	s_s_assemblage[n].s_s[number_s_s].tk = 298.15;
      }
      break;
    case 11:			/* temp */
    case 13:			/* tempc */
      ptr = next_char;
      j = 0;
      if (copy_token (token, &ptr, &l) != EMPTY)
      {
	j =
	  sscanf (token, SCANFORMAT, &(s_s_assemblage[n].s_s[number_s_s].tk));
      }
      if (j != 1)
      {
	sprintf (error_string,
		 "Expected temperature (Celcius) for parameters, assemblage %d, solid solution %s, using 25 C",
		 s_s_assemblage[n].n_user,
		 s_s_assemblage[n].s_s[number_s_s].name);
	warning_msg (error_string);
	s_s_assemblage[n].s_s[number_s_s].tk = 25.;
      }
      s_s_assemblage[n].s_s[number_s_s].tk += 273.15;
      break;
    case 14:			/* Thompson and Waldbaum */
      ptr = next_char;
      j = 0;
      for (i = 0; i < 2; i++)
      {
	if (copy_token (token, &ptr, &l) != EMPTY)
	{
	  j +=
	    sscanf (token, SCANFORMAT,
		    &(s_s_assemblage[n].s_s[number_s_s].p[i]));
	}
      }
      if (j != 2)
      {
	sprintf (error_string,
		 "Expected Wg2 and Wg1 Thompson-Waldbaum parameters to calculate a0 and a1, assemblage %d, solid solution %s",
		 s_s_assemblage[n].n_user,
		 s_s_assemblage[n].s_s[number_s_s].name);
	error_msg (error_string, CONTINUE);
	input_error++;
      }
      s_s_assemblage[n].s_s[number_s_s].input_case = 8;
      break;
    case 15:			/* Margules */
      ptr = next_char;
      j = 0;
      for (i = 0; i < 2; i++)
      {
	if (copy_token (token, &ptr, &l) != EMPTY)
	{
	  j +=
	    sscanf (token, SCANFORMAT,
		    &(s_s_assemblage[n].s_s[number_s_s].p[i]));
	}
      }
      if (j != 2)
      {
	sprintf (error_string,
		 "Expected alpha2 and alpha3 Margules parameters to calculate a0 and a1, assemblage %d, solid solution %s",
		 s_s_assemblage[n].n_user,
		 s_s_assemblage[n].s_s[number_s_s].name);
	error_msg (error_string, CONTINUE);
	input_error++;
      }
      s_s_assemblage[n].s_s[number_s_s].input_case = 9;
      break;
    case 16:			/* comp1 */
      if (count_comps < 2)
      {
	s_s_assemblage[n].s_s[number_s_s].count_comps = 2;
	count_comps = 2;
	s_s_assemblage[n].s_s[number_s_s].comps =
	  (struct s_s_comp *) PHRQ_realloc (s_s_assemblage[n].s_s[number_s_s].
					    comps,
					    (size_t) (count_comps) *
					    sizeof (struct s_s_comp));
	if (s_s_assemblage[n].s_s[number_s_s].comps == NULL)
	  malloc_error ();
      }
      /*
       *   Read phase name of component
       */
      ptr = next_char;
      copy_token (token, &ptr, &l);
      s_s_assemblage[n].s_s[number_s_s].comps[0].name = string_hsave (token);
      /*
       *   Read moles of component
       */
      if ((j = copy_token (token, &ptr, &l)) == EMPTY)
      {
	s_s_assemblage[n].s_s[number_s_s].comps[0].moles = NAN;
      }
      else
      {
	j = sscanf (token, SCANFORMAT, &dummy);
	s_s_assemblage[n].s_s[number_s_s].comps[0].moles = (LDBLE) dummy;
	if (j != 1)
	{
	  error_msg ("Expected moles of solid solution.", CONTINUE);
	  error_msg (line_save, CONTINUE);
	  input_error++;
	}
      }
      break;
    case 17:			/* comp2 */
      if (count_comps < 2)
      {
	s_s_assemblage[n].s_s[number_s_s].count_comps = 2;
	count_comps = 2;
	s_s_assemblage[n].s_s[number_s_s].comps =
	  (struct s_s_comp *) PHRQ_realloc (s_s_assemblage[n].s_s[number_s_s].
					    comps,
					    (size_t) (count_comps) *
					    sizeof (struct s_s_comp));
	if (s_s_assemblage[n].s_s[number_s_s].comps == NULL)
	  malloc_error ();
      }
      /*
       *   Read phase name of component
       */
      ptr = next_char;
      copy_token (token, &ptr, &l);
      s_s_assemblage[n].s_s[number_s_s].comps[1].name = string_hsave (token);
      /*
       *   Read moles of component
       */
      if ((j = copy_token (token, &ptr, &l)) == EMPTY)
      {
	s_s_assemblage[n].s_s[number_s_s].comps[1].moles = NAN;
      }
      else
      {
	j = sscanf (token, SCANFORMAT, &dummy);
	s_s_assemblage[n].s_s[number_s_s].comps[1].moles = (LDBLE) dummy;
	if (j != 1)
	{
	  error_msg ("Expected moles of solid solution.", CONTINUE);
	  error_msg (line_save, CONTINUE);
	  input_error++;
	}
      }
      break;
/* 
 * New solid solution
 */
    case OPTION_DEFAULT:
      number_s_s = count_s_s++;
      /*
       *   Make space, set default
       */

      /* realloc space for one s_s */
      s_s_assemblage[n].s_s =
	(struct s_s *) PHRQ_realloc (s_s_assemblage[n].s_s,
				     (size_t) (count_s_s +
					       1) * sizeof (struct s_s));
      if (s_s_assemblage[n].s_s == NULL)
	malloc_error ();

      /* malloc space for one component */
      s_s_assemblage[n].s_s[number_s_s].comps =
	(struct s_s_comp *) PHRQ_malloc ((size_t) sizeof (struct s_s_comp));
      if (s_s_assemblage[n].s_s[number_s_s].comps == NULL)
	malloc_error ();
      count_comps = 0;
      s_s_assemblage[n].s_s[number_s_s].a0 = 0.0;
      s_s_assemblage[n].s_s[number_s_s].a1 = 0.0;
      s_s_assemblage[n].s_s[number_s_s].count_comps = 0;
      s_s_assemblage[n].s_s[number_s_s].input_case = 0;
      s_s_assemblage[n].s_s[number_s_s].miscibility = FALSE;
      s_s_assemblage[n].s_s[number_s_s].p[0] = 0.0;
      s_s_assemblage[n].s_s[number_s_s].p[1] = 0.0;
      s_s_assemblage[n].s_s[number_s_s].tk = 298.15;
      s_s_assemblage[n].s_s[number_s_s].comps->name = NULL;
      s_s_assemblage[n].s_s[number_s_s].comps->phase = NULL;
      /*
       *   Read solid solution name
       */
      ptr = line;
      copy_token (token, &ptr, &l);
      s_s_assemblage[n].s_s[number_s_s].name = string_hsave (token);
      s_s_assemblage[n].s_s[number_s_s].total_moles = NAN;
      break;
    }
    if (return_value == EOF || return_value == KEYWORD)
      break;
  }
  for (i = 0; i < s_s_assemblage[n].count_s_s; i++)
  {
    if (s_s_assemblage[n].s_s[i].p[0] != 0.0 ||
	s_s_assemblage[n].s_s[i].p[1] != 0.0)
    {
      if (s_s_assemblage[n].s_s[number_s_s].count_comps != 2)
      {
	sprintf (error_string,
		 "Solid solution, %s, is nonideal. Must define exactly two components (-comp1 and -comp2).",
		 s_s_assemblage[n].s_s[number_s_s].name);
	error_msg (error_string, CONTINUE);
	input_error++;
      }
    }
  }
/*
 *   Sort components by name (lowercase)
 */

  s_s_assemblage[n].count_s_s = count_s_s;
  qsort (s_s_assemblage[n].s_s,
	 (size_t) count_s_s, (size_t) sizeof (struct s_s), s_s_compare);
  return (return_value);
}

/* ---------------------------------------------------------------------- */
int
read_llnl_aqueous_model_parameters (void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads aqueous model parameters
 *
 *      Arguments:
 *         none
 *
 *      Returns:
 *         KEYWORD if keyword encountered, input_error may be incremented if
 *                    a keyword is encountered in an unexpected position
 *         EOF     if eof encountered while reading mass balance concentrations
 *         ERROR   if error occurred reading data
 *
 */
  extern int check_line_return;	/* input.c */
  int i, count_alloc;
  char token[MAX_LENGTH];


  int return_value, opt;
  char *next_char;
  const char *opt_list[] = {
    "temperatures",		/* 0 */
    "temperature",		/* 1 */
    "temp",			/* 2 */
    "adh",			/* 3 */
    "debye_huckel_a",		/* 4 */
    "dh_a",			/* 5 */
    "bdh",			/* 6 */
    "debye_huckel_b",		/* 7 */
    "dh_b",			/* 8 */
    "bdot",			/* 9 */
    "b_dot",			/* 10 */
    "c_co2",			/* 11 */
    "co2_coefs"			/* 12 */
  };
  int count_opt_list = 13;
/*
 *   Initialize
 */
/*
 *   Read aqueous model parameters
 */
  return_value = UNKNOWN;
  opt = get_option (opt_list, count_opt_list, &next_char);
  for (;;)
  {
    next_char = line;
    if (opt >= 0)
    {
      copy_token (token, &next_char, &i);
    }
    switch (opt)
    {
    case OPTION_EOF:		/* end of file */
      return_value = EOF;
      break;
    case OPTION_KEYWORD:	/* keyword */
      return_value = KEYWORD;
      break;
    case OPTION_DEFAULT:
    case OPTION_ERROR:
      input_error++;
      error_msg ("Unknown input in LLNL_AQUEOUS_MODEL_PARAMETERS keyword.",
		 CONTINUE);
      error_msg (line_save, CONTINUE);
      break;

/* 
 * New component
 */
    case 0:			/* temperatures */
    case 1:			/* temperature */
    case 2:			/* temp */
      count_alloc = 1;
      llnl_count_temp = 0;
      i =
	read_lines_doubles (next_char, &(llnl_temp), &(llnl_count_temp),
			    &(count_alloc), opt_list, count_opt_list, &opt);
      /*
         ptr = next_char;
         llnl_temp = read_list_doubles(&ptr, &count);
         llnl_count_temp = count;
       */
      break;
    case 3:			/* adh */
    case 4:			/* debye_huckel_a */
    case 5:			/* dh_a */
      count_alloc = 1;
      llnl_count_adh = 0;
      i =
	read_lines_doubles (next_char, &(llnl_adh), &(llnl_count_adh),
			    &(count_alloc), opt_list, count_opt_list, &opt);
      /*
         ptr = next_char;
         llnl_adh = read_list_doubles(&ptr, &count);
         llnl_count_adh = count;
       */
      break;
    case 6:			/* bdh */
    case 7:			/* debye_huckel_b */
    case 8:			/* dh_b */
      count_alloc = 1;
      llnl_count_bdh = 0;
      i =
	read_lines_doubles (next_char, &(llnl_bdh), &(llnl_count_bdh),
			    &(count_alloc), opt_list, count_opt_list, &opt);
      /*
         ptr = next_char;
         llnl_bdh = read_list_doubles(&ptr, &count);
         llnl_count_bdh = count;
       */
      break;
    case 9:			/* bdot */
    case 10:			/* b_dot */
      count_alloc = 1;
      llnl_count_bdot = 0;
      i =
	read_lines_doubles (next_char, &(llnl_bdot), &(llnl_count_bdot),
			    &(count_alloc), opt_list, count_opt_list, &opt);
      /*
         ptr = next_char;
         llnl_bdot = read_list_doubles(&ptr, &count);
         llnl_count_bdot = count;
       */
      break;
    case 11:			/* c_co2 */
    case 12:			/* co2_coefs */
      count_alloc = 1;
      llnl_count_co2_coefs = 0;
      i =
	read_lines_doubles (next_char, &(llnl_co2_coefs),
			    &(llnl_count_co2_coefs), &(count_alloc), opt_list,
			    count_opt_list, &opt);
      /*
         ptr = next_char;
         llnl_co2_coefs = read_list_doubles(&ptr, &count);
         llnl_count_co2_coefs = count;
       */
      break;
    }
    return_value = check_line_return;
    if (return_value == EOF || return_value == KEYWORD)
      break;
  }
  /* check consistency */
  if ((llnl_count_temp <= 0) ||
      (llnl_count_temp != llnl_count_adh) ||
      (llnl_count_temp != llnl_count_bdh) ||
      (llnl_count_temp != llnl_count_bdot))
  {
    error_msg
      ("Must define equal number (>0) of temperatures, dh_a, dh_b, and bdot parameters\nin LLNL_AQUEOUS_MODEL",
       CONTINUE);
    input_error++;
  }
  if (llnl_count_co2_coefs != 5)
  {
    error_msg
      ("Must define 5 CO2 activity coefficient parameters in LLNL_AQUEOUS_MODEL",
       CONTINUE);
    input_error++;
  }
  for (i = 1; i < llnl_count_temp; i++)
  {
    if (llnl_temp[i - 1] > llnl_temp[i])
    {
      error_msg
	("Temperatures must be in ascending order in LLNL_AQUEOUS_MODEL",
	 CONTINUE);
      input_error++;
    }
  }

  return (return_value);
}

/* ---------------------------------------------------------------------- */
int
read_lines_doubles (char *next_char, LDBLE ** d, int *count_d,
		    int *count_alloc, const char **opt_list,
		    int count_opt_list, int *opt)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads LDBLEs on line starting at next_char
 *      and on succeeding lines. Appends to d.
 *      Stops at KEYWORD, OPTION, and EOF
 *
 *      Input Arguments:
 *         next_char    points to line to read from
 *         d            points to array of LDBLEs, must be malloced
 *         count_d      number of elements in array
 *         count_alloc  number of elements malloced
 *
 *      Output Arguments:
 *         d            points to array of LDBLEs, may have been
 *                          realloced
 *         count_d      updated number of elements in array
 *         count_alloc  updated of elements malloced
 *
 *      Returns:
 *         KEYWORD
 *         OPTION
 *         EOF
 *         ERROR if any errors reading LDBLEs
 */

  if (read_line_doubles (next_char, d, count_d, count_alloc) == ERROR)
  {
    return (ERROR);
  }
  for (;;)
  {
    *opt = get_option (opt_list, count_opt_list, &next_char);
    if (*opt == OPTION_KEYWORD || *opt == OPTION_EOF || *opt == OPTION_ERROR)
    {
      break;
    }
    else if (*opt >= 0)
    {
      break;
    }
    next_char = line;
    if (read_line_doubles (next_char, d, count_d, count_alloc) == ERROR)
    {
      return (ERROR);
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
read_line_doubles (char *next_char, LDBLE ** d, int *count_d,
		   int *count_alloc)
/* ---------------------------------------------------------------------- */
{
  int i, j, l, n;
  LDBLE value;
  char token[MAX_LENGTH];

  for (;;)
  {
    j = copy_token (token, &next_char, &l);
    if (j == EMPTY)
    {
      break;
    }
    if (j != DIGIT)
    {
      return (ERROR);
    }
    if (replace ("*", " ", token) == TRUE)
    {
      if (sscanf (token, "%d" SCANFORMAT, &n, &value) != 2)
      {
	return (ERROR);
      }
    }
    else
    {
      sscanf (token, SCANFORMAT, &value);
      n = 1;
    }
    for (;;)
    {
      if ((*count_d) + n > (*count_alloc))
      {
	*count_alloc *= 2;
	*d =
	  (LDBLE *) PHRQ_realloc (*d,
				  (size_t) (*count_alloc) * sizeof (LDBLE));
	if (*d == NULL)
	  malloc_error ();
      }
      else
      {
	break;
      }
    }
    for (i = 0; i < n; i++)
    {
      (*d)[(*count_d) + i] = value;
    }
    *count_d += n;
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
next_keyword_or_option (const char **opt_list, int count_opt_list)
/* ---------------------------------------------------------------------- */
{
/*
 *   Reads to next keyword or option or eof
 *
 *   Returns:
 *       KEYWORD
 *       OPTION
 *       EOF
 */
  int opt;
  char *next_char;

  for (;;)
  {
    opt = get_option (opt_list, count_opt_list, &next_char);
    if (opt == OPTION_EOF)
    {				/* end of file */
      break;
    }
    else if (opt == OPTION_KEYWORD)
    {				/* keyword */
      break;
    }
    else if (opt >= 0 && opt < count_opt_list)
    {
      break;
    }
    else
    {
      error_msg ("Expected a keyword or option.", CONTINUE);
      error_msg (line_save, CONTINUE);
      input_error++;
    }
  }
  return (opt);
}

/* ---------------------------------------------------------------------- */
int
read_named_logk (void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads K's that can be used to calculate K's for species
 *
 *      Arguments:
 *         none
 *
 *      Returns:
 *         KEYWORD if keyword encountered, input_error may be incremented if
 *                    a keyword is encountered in an unexpected position
 *         EOF     if eof encountered while reading mass balance concentrations
 *         ERROR   if error occurred reading data
 *
 */

  int j, l;
  int i, empty;
  struct logk *logk_ptr;
  char token[MAX_LENGTH];

  int return_value, opt, opt_save;
  char *next_char;
  const char *opt_list[] = {
    "log_k",			/* 0 */
    "logk",			/* 1 */
    "delta_h",			/* 2 */
    "deltah",			/* 3 */
    "analytical_expression",	/* 4 */
    "a_e",			/* 5 */
    "ae",			/* 6 */
    "ln_alpha1000",		/* 7 */
    "add_logk",			/* 8 */
    "add_log_k"			/* 9 */
  };
  int count_opt_list = 10;
  logk_ptr = NULL;
/*
 *   Read name followed by options
 */
  opt_save = OPTION_DEFAULT;
  return_value = UNKNOWN;
  for (;;)
  {
    opt = get_option (opt_list, count_opt_list, &next_char);
    if (opt == OPTION_DEFAULT)
    {
      opt = opt_save;
    }
    switch (opt)
    {
    case OPTION_EOF:		/* end of file */
      return_value = EOF;
      break;
    case OPTION_KEYWORD:	/* keyword */
      return_value = KEYWORD;
      break;
    case OPTION_ERROR:
      input_error++;
      error_msg ("Unknown input in SPECIES keyword.", CONTINUE);
      error_msg (line_save, CONTINUE);
      break;
    case 0:			/* log_k */
    case 1:			/* logk */
      if (logk_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      read_log_k_only (next_char, &logk_ptr->log_k[0]);
      logk_copy2orig (logk_ptr);
      opt_save = OPTION_DEFAULT;
      break;
    case 2:			/* delta_h */
    case 3:			/* deltah */
      if (logk_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      read_delta_h_only (next_char, &logk_ptr->log_k[1],
			 &logk_ptr->original_units);
      logk_copy2orig (logk_ptr);
      opt_save = OPTION_DEFAULT;
      break;
    case 4:			/* analytical_expression */
    case 5:			/* a_e */
    case 6:			/* ae */
      if (logk_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      read_analytical_expression_only (next_char, &(logk_ptr->log_k[2]));
      logk_copy2orig (logk_ptr);
      opt_save = OPTION_DEFAULT;
      break;
    case 7:			/* ln_alpha1000 */
      if (logk_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      empty = TRUE;
      for (i = 2; i < 7; i++)
      {
	if (logk_ptr->log_k[i] != 0.0)
	{
	  empty = FALSE;
	  logk_ptr->log_k[i] = 0.0;
	}
      }
      if (empty == FALSE)
      {
	sprintf (error_string,
		 "Analytical expression previously defined for %s in NAMED_EXPRESSIONS\nAnalytical expression will be overwritten.",
		 logk_ptr->name);
	warning_msg (error_string);
      }
      read_analytical_expression_only (next_char, &(logk_ptr->log_k[2]));
      for (i = 2; i < 7; i++)
      {
	logk_ptr->log_k[i] /= 1000. * LOG_10;
      }
      logk_copy2orig (logk_ptr);
      opt_save = OPTION_DEFAULT;
      break;
    case 8:			/* add_logk */
    case 9:			/* add_log_k */
      if (logk_ptr == NULL)
      {
	sprintf (error_string, "No reaction defined before option, %s.",
		 opt_list[opt]);
	error_msg (error_string, CONTINUE);
	input_error++;
	break;
      }
      if (logk_ptr->count_add_logk == 0)
      {
	logk_ptr->add_logk =
	  (struct name_coef *) PHRQ_malloc (sizeof (struct name_coef));
	if (logk_ptr->add_logk == NULL)
	  malloc_error ();
      }
      else
      {
	logk_ptr->add_logk =
	  (struct name_coef *) PHRQ_realloc (logk_ptr->add_logk,
					     (size_t) ((logk_ptr->
							count_add_logk +
							1) *
						       sizeof (struct
							       name_coef)));
	if (logk_ptr->add_logk == NULL)
	  malloc_error ();
      }
      /* read name */
      if (copy_token (token, &next_char, &i) == EMPTY)
      {
	input_error++;
	sprintf (error_string, "Expected the name of a NAMED_EXPRESSION.");
	error_msg (error_string, CONTINUE);
	break;
      }
      logk_ptr->add_logk[logk_ptr->count_add_logk].name =
	string_hsave (token);
      /* read coef */
      i =
	sscanf (next_char, SCANFORMAT,
		&logk_ptr->add_logk[logk_ptr->count_add_logk].coef);
      if (i <= 0)
      {
	logk_ptr->add_logk[logk_ptr->count_add_logk].coef = 1;
      }
      logk_ptr->count_add_logk++;
      opt_save = OPTION_DEFAULT;
      break;
    case OPTION_DEFAULT:
/*
 *   Get space for logk information
 */
      logk_ptr = NULL;
      j = copy_token (token, &next_char, &l);

      logk_ptr = logk_store (token, TRUE);
/*
 *   Get pointer to each species in the reaction, store new species if necessary
 */
      opt_save = OPTION_DEFAULT;
      break;
    }
    if (return_value == EOF || return_value == KEYWORD)
      break;
  }
  return (return_value);
}

/* ---------------------------------------------------------------------- */
int
read_copy (void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Reads solution, 
 *         equilibrium_phases,
 *         exchange,
 *         surface,
 *         solid_solution,
 *         gas_phase,
 *         kinetics,
 *         mix, 
 *         reaction, 
 *         reaction_temperature
 *
 */
  int i, l, n, n_user, n_user_start, n_user_end, return_value;
  char *ptr;
  char token[MAX_LENGTH], token1[MAX_LENGTH];;
/*
 *   Read "copy"
 */
  ptr = line;
  copy_token (token, &ptr, &l);
/*
 *   Read keyword
 */
  copy_token (token, &ptr, &l);
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
    input_error++;
    error_msg
      ("Expecting keyword solution, mix, kinetics, reaction, reaction_temperature, equilibrium_phases, exchange, surface, gas_phase, or solid_solutions.",
       CONTINUE);
    error_msg (line_save, CONTINUE);
    check_line ("End of use", FALSE, TRUE, TRUE, TRUE);
    /* empty, eof, keyword, print */
    return (ERROR);
  }
/*
 *   Read source index
 */
  strcpy (token1, token);
  i = copy_token (token, &ptr, &l);
  if (i == DIGIT)
  {
    sscanf (token, "%d", &n_user);
    if (n_user < 0)
    {
      error_msg ("Source index number must be a positive integer.", CONTINUE);
      error_msg (line_save, CONTINUE);
      input_error++;
      return (ERROR);
    }
    if (strstr (token, "-") != NULL)
    {
      error_msg ("COPY does not accept a range of numbers for source index",
		 CONTINUE);
      error_msg (line_save, CONTINUE);
      input_error++;
      return (ERROR);
    }
  }
  else
  {
    error_msg ("Source index number must be a positive integer.", CONTINUE);
    error_msg (line_save, CONTINUE);
    input_error++;
    return (ERROR);
  }
/*
 *   Read target index or range of indices
 */
  i = copy_token (token, &ptr, &l);
  if (i == DIGIT)
  {
    replace ("-", " ", token);
    n = sscanf (token, "%d%d", &n_user_start, &n_user_end);
    if (n == 1)
    {
      n_user_end = n_user_start;
    }
    if (n_user_start < 0)
    {
      error_msg ("Target index number must be a positive integer.", CONTINUE);
      error_msg (line_save, CONTINUE);
      input_error++;
      return (ERROR);
    }
  }
  else
  {
    error_msg ("Target index number must be a positive integer.", CONTINUE);
    error_msg (line_save, CONTINUE);
    input_error++;
    return (ERROR);
  }

  switch (next_keyword)
  {
  case 4:			/* Solution */
    copier_add (&copy_solution, n_user, n_user_start, n_user_end);
    break;
  case 6:			/* Pure phases */
  case 26:
  case 27:
  case 28:
  case 29:
    copier_add (&copy_pp_assemblage, n_user, n_user_start, n_user_end);
    break;
  case 7:			/* Reaction */
    copier_add (&copy_irrev, n_user, n_user_start, n_user_end);
    break;
  case 8:			/* Mix */
    copier_add (&copy_mix, n_user, n_user_start, n_user_end);
    break;
  case 13:			/* Ex */
    copier_add (&copy_exchange, n_user, n_user_start, n_user_end);
    break;
  case 16:			/* Surface */
    copier_add (&copy_surface, n_user, n_user_start, n_user_end);
    break;
  case 17:			/* Temperature */
    copier_add (&copy_temperature, n_user, n_user_start, n_user_end);
    break;
  case 19:			/* Gas */
    copier_add (&copy_gas_phase, n_user, n_user_start, n_user_end);
    break;
  case 33:			/* Kinetics */
    copier_add (&copy_kinetics, n_user, n_user_start, n_user_end);
    break;
  case 40:			/* solid_solutions */
  case 41:			/* solid_solution */
    copier_add (&copy_s_s_assemblage, n_user, n_user_start, n_user_end);
    break;
  }
  return_value = check_line ("End of COPY", FALSE, TRUE, TRUE, TRUE);
  /* empty, eof, keyword, print */
  return (return_value);
}
