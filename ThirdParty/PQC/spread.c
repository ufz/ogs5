#define EXTERNAL extern
#include "global.h"
#include "phqalloc.h"
#include "output.h"
#include "phrqproto.h"

static char const svnid[] =
  "$Id: spread.c 4 2009-04-21 17:29:29Z delucia $";

#define STATIC static

#define STRING 11
#define NUMBER 12
#define MIXED 13
#define EOL 14

#define OPTION_EOF -1
#define OPTION_KEYWORD -2
#define OPTION_ERROR -3
#define OPTION_DEFAULT -4
#define OPT_1 -5

static int copy_token_tab (char *token_ptr, char **ptr, int *length);
static int get_option_string (const char **opt_list, int count_opt_list,
			      char **next_char);
static int spread_row_free (struct spread_row *spread_row_ptr);
static int spread_row_to_solution (struct spread_row *heading,
				   struct spread_row *units,
				   struct spread_row *data,
				   struct defaults defaults);
static struct spread_row *string_to_spread_row (char *string);
#ifdef PHREEQCI_GUI
static void add_row (struct spread_row *spread_row_ptr);
static void copy_defaults (struct defaults *dest_ptr,
			   struct defaults *src_ptr);
void free_spread (void);
static struct spread_row *copy_row (struct spread_row *spread_row_ptr);
#endif

/* ---------------------------------------------------------------------- */
int
read_solution_spread (void)
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
  struct spread_row *heading, *row_ptr, *units;
  int i, j, l, j1, num, count;
  int strings, numbers;
  int spread_lines;
  char token[MAX_LENGTH], token1[MAX_LENGTH];
  char *ptr;
  struct defaults defaults = {
    25,
    1,
    "mmol/kgw",
    "pe",
    7,
    4,
    1,
    1,
    iso_defaults,
  };
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
    "water",			/* 10 */
    "isotope_uncertainty",	/* 11 */
    "uncertainty",		/* 12 */
    "uncertainties"		/* 13 */
  };
  int count_opt_list = 14;
  if (svnid == NULL)
    fprintf (stderr, " ");

  heading = NULL;
  units = NULL;
  defaults.count_iso = count_iso_defaults;
  defaults.iso =
    (struct iso *) PHRQ_malloc ((size_t) defaults.count_iso *
				sizeof (struct iso));
  if (defaults.iso == NULL)
    malloc_error ();
  memcpy (defaults.iso, iso_defaults,
	  (size_t) defaults.count_iso * sizeof (struct iso));
  return_value = UNKNOWN;
  spread_lines = 0;
/*
 *   Loop on solutions
 */
  for (;;)
  {
    opt = get_option (opt_list, count_opt_list, &next_char);
    if (spread_lines == 0 && opt != OPTION_DEFAULT)
    {
      row_ptr = string_to_spread_row (line);
      count = 0;
      ptr = line;
      numbers = 0;
      strings = 0;
      while (((j = copy_token (token, &ptr, &l)) != EMPTY))
      {
	count++;
	if (j == UPPER || j == LOWER)
	  strings++;
	if (j == DIGIT)
	  numbers++;
      }
#ifdef SKIP
      for (i = 0; i < row_ptr->count; i++)
      {
	if (row_ptr->type_vector[i] == STRING)
	{
	  strings++;
	}
	else if (row_ptr->type_vector[i] == NUMBER)
	{
	  numbers++;
	}
      }
#endif
      /*
       * Is 2nd token all number
       */
      ptr = line;
      copy_token (token, &ptr, &l);
      j = copy_token (token, &ptr, &l);
      num = FALSE;
      if (j == DIGIT)
      {
	strtod (token, &ptr);
	j1 = copy_token (token1, &ptr, &l);
	if (j1 != EMPTY)
	{
	  num = FALSE;
	}
	else
	{
	  num = TRUE;
	}
      }

      /*
       *   Starts with hyphen
       */
      ptr = line;
      copy_token (token, &ptr, &l);
      if (token[0] == '-')
      {
	opt = opt;
      }
      else
      {
	switch (opt)
	{
	case 0:		/* temp */
	case 1:		/* temperature */
	case 2:		/* dens */
	case 3:		/* density */
	case 10:		/* water */
	  if (count == 2 && num == TRUE)
	  {
	    opt = opt;
	  }
	  else
	  {
	    opt = OPTION_DEFAULT;
	  }
	  break;
	case 6:		/* ph */
	case 7:		/* pe */
	  if ((count == 2 || count == 3 || count == 4) && num == TRUE)
	  {
	    opt = opt;
	  }
	  else
	  {
	    opt = OPTION_DEFAULT;
	  }
	  break;
	case 5:		/* redox */
	case 4:		/* units */
	case 8:		/* unit */
	  if (count == 2)
	  {
	    opt = opt;
	  }
	  else
	  {
	    opt = OPTION_DEFAULT;
	  }
	  break;
	case 9:		/* isotope */
	  if (row_ptr->count > 4)
	  {
	    opt = OPTION_DEFAULT;
	  }
	  else
	  {
	    opt = opt;
	  }
	  break;
	case 11:		/* isotope_uncertainty */
	case 12:		/* uncertainty */
	case 13:		/* uncertainties */
	  if (row_ptr->count > 3)
	  {
	    opt = OPTION_DEFAULT;
	  }
	  else
	  {
	    opt = opt;
	  }
	  break;
	}
      }
      spread_row_free (row_ptr);
    }
    if (opt == OPTION_DEFAULT)
    {
      if (spread_lines == 0)
      {
	opt = 100;
      }
      spread_lines++;
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
    case OPTION_DEFAULT:	/* solution definition */
      row_ptr = string_to_spread_row (line);
      if (spread_lines == 2)
      {
	numbers = 0;
	strings = 0;
	for (i = 0; i < heading->count; i++)
	{
	  if (row_ptr->type_vector[i] == STRING)
	  {
	    strings++;
	  }
	  else if (row_ptr->type_vector[i] == NUMBER)
	  {
	    numbers++;
	  }
#ifdef SKIP
	  if (row_ptr->type_vector[i] == STRING &&
	      (strcmp_nocase (heading->char_vector[i], "units") != 0) &&
	      (strcmp_nocase (heading->char_vector[i], "unit") != 0) &&
	      (strcmp_nocase (heading->char_vector[i], "description") != 0) &&
	      (strcmp_nocase (heading->char_vector[i], "desc") != 0) &&
	      (strcmp_nocase (heading->char_vector[i], "descriptor") != 0) &&
	      (strcmp_nocase (heading->char_vector[i], "redox") != 0))
	  {
	    break;
	  }
#endif
	}
#ifdef SKIP
	if (i < heading->count)
	{
	  units = row_ptr;
	  break;
	}
#endif
	if (numbers == 0)
	{
	  units = row_ptr;
	  break;
	}
      }
      spread_row_to_solution (heading, units, row_ptr, defaults);
#ifdef PHREEQCI_GUI
      add_row (row_ptr);
#endif
      spread_row_free (row_ptr);
      break;
    case 0:			/* temperature */
    case 1:
      sscanf (next_char, SCANFORMAT, &(defaults.temp));
      break;
    case 2:			/* density */
    case 3:
      sscanf (next_char, SCANFORMAT, &(defaults.density));
      break;
    case 4:			/* units */
    case 8:			/* unit */
      if (copy_token (token, &next_char, &l) == EMPTY)
	break;
      if (check_units (token, FALSE, FALSE, NULL, TRUE) == OK)
      {
	defaults.units = string_hsave (token);
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
	defaults.redox = string_hsave (token);
      }
      else
      {
	input_error++;
      }
      break;
    case 6:			/* ph */
      copy_token (token, &next_char, &l);
      sscanf (token, SCANFORMAT, &(defaults.ph));
      if (copy_token (token, &next_char, &l) != EMPTY)
      {
	warning_msg
	  ("Not possible to use phase name or saturation index in definition of default pH in SOLUTION_SPREAD.");
      }
      break;
    case 7:			/* pe */
      copy_token (token, &next_char, &l);
      sscanf (token, SCANFORMAT, &(defaults.pe));
      if (copy_token (token, &next_char, &l) != EMPTY)
      {
	warning_msg
	  ("Not possible to use phase name or saturation index in definition of default pe in SOLUTION_SPREAD.");
      }
      break;
    case 11:			/* isotope_uncertainty */
    case 12:			/* uncertainty */
    case 13:			/* uncertainties */
      if (copy_token (token, &next_char, &l) != DIGIT)
      {
	input_error++;
	sprintf (error_string, "Expected isotope name to"
		 " begin with an isotopic number.");
	error_msg (error_string, CONTINUE);
	continue;
      }
      for (i = 0; i < defaults.count_iso; i++)
      {
	if (strcmp (token, defaults.iso[i].name) == 0)
	{
	  break;
	}
      }
      if (i == defaults.count_iso)
      {
	defaults.iso =
	  (struct iso *) PHRQ_realloc (defaults.iso,
				       (size_t) (i +
						 1) * sizeof (struct iso));
	if (defaults.iso == NULL)
	  malloc_error ();
	defaults.iso[i].name = string_duplicate (token);
	defaults.iso[i].value = NAN;
	defaults.iso[i].uncertainty = NAN;
	defaults.count_iso++;
      }

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
	else
	{
	  sscanf (token, SCANFORMAT, &(defaults.iso[i].uncertainty));
	}
      }
      else
      {
	defaults.iso[i].uncertainty = NAN;
      }
      break;
    case 10:			/* water */
      j = copy_token (token, &next_char, &l);
      if (j != DIGIT)
      {
	input_error++;
	sprintf (error_string,
		 "Expected numeric value for mass of water in solution.");
	error_msg (error_string, CONTINUE);
      }
      else
      {
	sscanf (token, SCANFORMAT, &(defaults.water));
      }
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
      for (i = 0; i < defaults.count_iso; i++)
      {
	if (strcmp (token, defaults.iso[i].name) == 0)
	{
	  break;
	}
      }
      if (i == defaults.count_iso)
      {
	defaults.iso =
	  (struct iso *) PHRQ_realloc (defaults.iso,
				       (size_t) (i +
						 1) * sizeof (struct iso));
	if (defaults.iso == NULL)
	  malloc_error ();
	defaults.iso[i].name = string_duplicate (token);
	defaults.iso[i].value = NAN;
	defaults.iso[i].uncertainty = NAN;
	defaults.count_iso++;
      }
      /* read and store isotope ratio */
      if (copy_token (token, &next_char, &l) != DIGIT)
      {
	input_error++;
	sprintf (error_string,
		 "Expected numeric value for default isotope ratio.");
	error_msg (error_string, CONTINUE);
	break;
      }
      sscanf (token, SCANFORMAT, &(defaults.iso[i].value));
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
	else
	{
	  sscanf (token, SCANFORMAT, &(defaults.iso[i].uncertainty));
	}
      }
      break;
    case 100:			/* read headings */
      heading = string_to_spread_row (line);
      for (i = 0; i < heading->count; i++)
      {
	while (replace (" ", "", heading->char_vector[i]) == TRUE);
	while (replace (",", "_", heading->char_vector[i]) == TRUE);
      }

      break;
    }
    if (return_value == EOF || return_value == KEYWORD)
      break;
  }
#ifdef PHREEQCI_GUI
  if (heading)
    g_spread_sheet.heading = copy_row (heading);
  if (units)
    g_spread_sheet.units = copy_row (units);
  copy_defaults (&g_spread_sheet.defaults, &defaults);
#endif
  spread_row_free (heading);
  spread_row_free (units);
  /* free non-default iso names */
  for (i = count_iso_defaults; i < defaults.count_iso; i++)
  {
    defaults.iso[i].name = (char *) free_check_null (defaults.iso[i].name);
  }
  defaults.iso = (struct iso *) free_check_null (defaults.iso);
  return (return_value);
}

/* ---------------------------------------------------------------------- */
int
spread_row_to_solution (struct spread_row *heading, struct spread_row *units,
			struct spread_row *data, struct defaults defaults)
/* ---------------------------------------------------------------------- */
{
  int i, j, n, l, next_keyword_save;
  int n_user, n_user_end;
  int default_pe, alk;
  int count_isotopes;
  int max_mass_balance, count_mass_balance;
  char *ptr, *ptr1;
  char *description;
  char token[MAX_LENGTH], token1[MAX_LENGTH];
  char string[2 * MAX_LENGTH];
  LDBLE dummy;

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
    "water",			/* 10 */
    "description",		/* 11 */
    "desc",			/* 12 */
    "descriptor"		/* 13 */
  };
  int count_opt_list = 14;

/*
 *      look for solution number
 */
  n_user = -1;
  n_user_end = -1;
  description = string_duplicate ("");
  for (i = 0; i < heading->count; i++)
  {
    if (strcmp_nocase (heading->char_vector[i], "number") == 0)
    {
      break;
    }
  }
  if (i == heading->count || data->type_vector[i] == EMPTY
      || data->count <= i)
  {
    n_user = -1;
#ifdef SKIP
    for (i = 0; i < count_solution; i++)
    {
      if (n_user <= solution[i]->n_user)
      {
	n_user = solution[i]->n_user + 1;
      }
    }
#endif
  }
  else if (data->type_vector[i] == STRING)
  {
    input_error++;
    sprintf (error_string,
	     "Expected solution number or number range in 'number' column, found:  %s.",
	     data->char_vector[i]);
    error_msg (error_string, CONTINUE);
  }
  else
  {
    strcpy (string, "solution_s ");
    strcat (string, data->char_vector[i]);
    ptr = string;
    description = (char *) free_check_null (description);
    next_keyword_save = next_keyword;
    next_keyword = 42;
    read_number_description (ptr, &n_user, &n_user_end, &description);
    next_keyword = next_keyword_save;
  }
/*
 *   set up solution
 */

  if (n_user >= 0 && solution_bsearch (n_user, &n, FALSE) != NULL)
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
  solution[n]->tc = defaults.temp;
  solution[n]->ph = defaults.ph;
  solution[n]->density = defaults.density;
  solution[n]->solution_pe = defaults.pe;
  solution[n]->mass_water = defaults.water;
  solution[n]->ah2o = 1.0;
  solution[n]->mu = 1e-7;
  solution[n]->cb = 0.0;
  default_pe = 0;
  solution[n]->units = defaults.units;
  solution[n]->totals[0].description = NULL;
  count_mass_balance = 0;
  count_isotopes = 0;
  default_pe = pe_data_store (&(solution[n]->pe), defaults.redox);
/*
 *   Read concentration data
 */
  return_value = UNKNOWN;
  for (i = 0; i < heading->count; i++)
  {
    if (strcmp_nocase (heading->char_vector[i], "number") == 0)
      continue;
    if (strcmp_nocase (heading->char_vector[i], "uncertainty") == 0)
      continue;
    if (strcmp_nocase (heading->char_vector[i], "uncertainties") == 0)
      continue;
    if (strcmp_nocase (heading->char_vector[i], "isotope_uncertainty") == 0)
      continue;
    /*
     *  Copy in element name
     */
    if (heading->type_vector[i] == EMPTY)
      continue;
    strcpy (string, heading->char_vector[i]);
    strcat (string, " ");
    /*
     *  Copy in concentration data
     */
    if (i >= data->count || data->type_vector[i] == EMPTY)
      continue;
    strcat (string, data->char_vector[i]);
    strcat (string, " ");
    /*
     *  Copy in concentration data
     */
    if (units != NULL && i < units->count && units->type_vector[i] != EMPTY)
    {
      strcat (string, units->char_vector[i]);
    }
/*
 *   Parse string just like read_solution input 
 */
    next_char = string;
    opt = get_option_string (opt_list, count_opt_list, &next_char);
    if (opt == OPTION_DEFAULT && heading->type_vector[i] == NUMBER)
    {
      opt = 9;
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
      sscanf (next_char, SCANFORMAT, &(solution[n]->tc));
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
      next_char = string;
      if (read_conc (n, count_mass_balance, next_char) == ERROR)
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
      next_char = string;
      if (read_conc (n, count_mass_balance, next_char) == ERROR)
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
      next_char = string;
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

      /* read and store isotope ratio uncertainty */
      /* first choice is next column */
      if ((i + 1) < heading->count &&
	  (strcmp_nocase (heading->char_vector[i + 1], "uncertainty") == 0 ||
	   strcmp_nocase (heading->char_vector[i + 1],
			  "isotope_uncertainty") == 0
	   || strcmp_nocase (heading->char_vector[i + 1],
			     "uncertainties") == 0) && (i + 1) < data->count
	  && data->type_vector[i + 1] == NUMBER)
      {
	solution[n]->isotopes[count_isotopes].ratio_uncertainty =
	  data->d_vector[i + 1];
      }
      else
      {
	next_char = string;
	copy_token (token, &next_char, &l);
	for (j = 0; j < defaults.count_iso; j++)
	{
	  if (strcmp (token, defaults.iso[j].name) == 0)
	  {
	    solution[n]->isotopes[count_isotopes].ratio_uncertainty =
	      defaults.iso[j].uncertainty;
	    break;
	  }
	}
	if (j == defaults.count_iso)
	{
	  solution[n]->isotopes[count_isotopes].ratio_uncertainty = NAN;
	}
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
    case 11:			/* description */
    case 12:			/* desc */
    case 13:			/* descriptor */
      solution[n]->description =
	(char *) free_check_null (solution[n]->description);
      solution[n]->description = string_duplicate (next_char);
      break;
    case OPTION_DEFAULT:
/*
 *   Read concentration
 */
      next_char = string;
      if (copy_token (token, &next_char, &l) == LOWER)
	continue;
      next_char = string;
      if (read_conc (n, count_mass_balance, next_char) == ERROR)
      {
#ifdef SKIP
	input_error++;
	break;
#endif
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
/*
 *   Sort totals by description
 */
  qsort (solution[n]->totals,
	 (size_t) count_mass_balance,
	 (size_t) sizeof (struct conc), conc_compare);
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
  if (count_isotopes > 0)
  {
    qsort (solution[n]->isotopes,
	   (size_t) count_isotopes,
	   (size_t) sizeof (struct isotope), isotope_compare);
  }
  else
  {
    solution[n]->isotopes =
      (struct isotope *) free_check_null (solution[n]->isotopes);
  }
  return (return_value);
}

/* ---------------------------------------------------------------------- */
struct spread_row *
string_to_spread_row (char *string)
/* ---------------------------------------------------------------------- */
{
  int j, l;
  static int length = 10; /* possible memory error if length of line is smaller than previous line */
  char *token;
  char *ptr;
  struct spread_row *spread_row_ptr;
/*
 *   Allocate space
 */
  token =
    (char *) PHRQ_malloc (strlen(line)+1);
  if (token == NULL)
    malloc_error ();
  spread_row_ptr =
    (struct spread_row *) PHRQ_malloc ((size_t) sizeof (struct spread_row));
  if (spread_row_ptr == NULL)
    malloc_error ();
  spread_row_ptr->char_vector =
    (char **) PHRQ_malloc ((size_t) length * sizeof (char *));
  if (spread_row_ptr->char_vector == NULL)
    malloc_error ();
  spread_row_ptr->d_vector =
    (LDBLE *) PHRQ_malloc ((size_t) length * sizeof (LDBLE));
  if (spread_row_ptr->d_vector == NULL)
    malloc_error ();
  spread_row_ptr->type_vector =
    (int *) PHRQ_malloc ((size_t) length * sizeof (int));
  if (spread_row_ptr->type_vector == NULL)
    malloc_error ();
  spread_row_ptr->count = 0;
  spread_row_ptr->empty = 0;
  spread_row_ptr->string = 0;
  spread_row_ptr->number = 0;
  ptr = string;
/*
 *   Split by tabs, reallocate space
 */
  for (;;)
  {
    if (spread_row_ptr->count + 1 > length)
    {
      length *= 2;

      spread_row_ptr->char_vector =
	(char **) PHRQ_realloc (spread_row_ptr->char_vector,
				(size_t) length * sizeof (char *));
      if (spread_row_ptr->char_vector == NULL)
	malloc_error ();

      spread_row_ptr->d_vector =
	(LDBLE *) PHRQ_realloc (spread_row_ptr->d_vector,
				(size_t) length * sizeof (LDBLE));
      if (spread_row_ptr->d_vector == NULL)
	malloc_error ();

      spread_row_ptr->type_vector =
	(int *) PHRQ_realloc (spread_row_ptr->type_vector,
			      (size_t) length * sizeof (int));
      if (spread_row_ptr->type_vector == NULL)
	malloc_error ();
    }
    j = copy_token_tab (token, &ptr, &l);
    if (j == EOL)
      break;
    spread_row_ptr->char_vector[spread_row_ptr->count] =
      string_duplicate (token);
    spread_row_ptr->d_vector[spread_row_ptr->count] = NAN;
    if (j == EMPTY || l == 0)
    {
      spread_row_ptr->empty++;
      spread_row_ptr->type_vector[spread_row_ptr->count] = EMPTY;
    }
    else if (j == UPPER || j == LOWER)
    {
      spread_row_ptr->string++;
      spread_row_ptr->type_vector[spread_row_ptr->count] = STRING;
    }
    else if (j == DIGIT)
    {
      spread_row_ptr->number++;
      spread_row_ptr->d_vector[spread_row_ptr->count] = strtod (token, NULL);
      spread_row_ptr->type_vector[spread_row_ptr->count] = NUMBER;
    }
    spread_row_ptr->count++;
  }
/*
 *   Clean up and return
 */
  if (spread_row_ptr->count == 0)
  {
    spread_row_ptr->char_vector =
      (char **) free_check_null (spread_row_ptr->char_vector);
    spread_row_ptr->d_vector =
      (LDBLE *) free_check_null (spread_row_ptr->d_vector);
    spread_row_ptr->type_vector =
      (int *) free_check_null (spread_row_ptr->type_vector);
  }
  else
  {
/*  Do not realloc to smaller size, memory error */
/*
    spread_row_ptr->char_vector =
      (char **) PHRQ_realloc (spread_row_ptr->char_vector,
			      (size_t) spread_row_ptr->count *
			      sizeof (char *));
    if (spread_row_ptr->char_vector == NULL)
      malloc_error ();
    spread_row_ptr->d_vector =
      (LDBLE *) PHRQ_realloc (spread_row_ptr->d_vector,
			      (size_t) spread_row_ptr->count *
			      sizeof (LDBLE));
    if (spread_row_ptr->d_vector == NULL)
      malloc_error ();
    spread_row_ptr->type_vector =
      (int *) PHRQ_realloc (spread_row_ptr->type_vector,
			    (size_t) spread_row_ptr->count * sizeof (int));
    if (spread_row_ptr->type_vector == NULL)
      malloc_error ();
*/
  }
  token = (char *) free_check_null (token);
  return (spread_row_ptr);
}

/* ---------------------------------------------------------------------- */
int
spread_row_free (struct spread_row *spread_row_ptr)
/* ---------------------------------------------------------------------- */
{
  int i;

  if (spread_row_ptr == NULL)
    return (OK);
  for (i = 0; i < spread_row_ptr->count; i++)
  {
    spread_row_ptr->char_vector[i] =
      (char *) free_check_null (spread_row_ptr->char_vector[i]);
  }

  spread_row_ptr->char_vector =
    (char **) free_check_null (spread_row_ptr->char_vector);
  spread_row_ptr->d_vector =
    (double *) free_check_null (spread_row_ptr->d_vector);
  spread_row_ptr->type_vector =
    (int *) free_check_null (spread_row_ptr->type_vector);
  spread_row_ptr = (struct spread_row *) free_check_null (spread_row_ptr);
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
copy_token_tab (char *token_ptr, char **ptr, int *length)
/* ---------------------------------------------------------------------- */
{
/*
 *   Copies from **ptr to *token_ptr until first tab is encountered.
 *
 *   Arguments:
 *      *token_ptr  output, place to store token
 *
 *     **ptr        input, character string to read token from
 *                  output, next position after token
 *
 *       length     output, length of token
 *
 *   Returns:
 *      UPPER,
 *      LOWER,
 *      DIGIT,
 *      EMPTY,
 *      EOL,
 *      UNKNOWN.
 */
  int i, j, return_value;
  char c;
/*
 *   Strip leading spaces
 */
  while ((c = **ptr) == ' ')
    (*ptr)++;
/*
 *   Check what we have
 */
  if (isupper ((int) c))
  {
    return_value = UPPER;
  }
  else if (islower ((int) c))
  {
    return_value = LOWER;
  }
  else if (isdigit ((int) c) || c == '.' || c == '-')
  {
    return_value = DIGIT;
  }
  else if (c == '\0')
  {
    return_value = EOL;
    return (return_value);
  }
  else if (c == '\t')
  {
    return_value = EMPTY;
  }
  else
  {
    return_value = UNKNOWN;
  }
/*
 *   Begin copying to token
 */
  i = 0;
#ifdef SKIP
  while ((c = **ptr) != '\t' && c != '\0')
  {
    token_ptr[i] = c;
    (*ptr)++;
    i++;
  }
#endif
  for (;;)
  {
    c = **ptr;
    if (c == '\t')
    {
      (*ptr)++;
      break;
    }
    else if (c == '\0')
    {
      break;
    }
    else
    {
      token_ptr[i] = c;
      (*ptr)++;
      i++;
    }
  }
  token_ptr[i] = '\0';
  *length = i;
/*
 *   Strip trailing spaces
 */
  for (j = i - 1; j >= 0; j--)
  {
    if (j != ' ')
      break;
  }
  if (j != i - 1)
  {
    token_ptr[j + 1] = '\0';
    *length = j + 1;
  }
  return (return_value);
}

/* ---------------------------------------------------------------------- */
static int
get_option_string (const char **opt_list, int count_opt_list,
		   char **next_char)
/* ---------------------------------------------------------------------- */
{
/*
 *   Read a line and check for options
 */
  int j;
  int opt_l, opt;
  char *opt_ptr;
  char option[MAX_LENGTH];

  opt_ptr = *next_char;
  if (opt_ptr[0] == '-')
  {
    opt_ptr++;
    copy_token (option, &opt_ptr, &opt_l);
    if (find_option (&(option[1]), &opt, opt_list, count_opt_list, FALSE) ==
	OK)
    {
      j = opt;
      *next_char = opt_ptr;
    }
    else
    {
      error_msg ("Unknown option.", CONTINUE);
      error_msg (*next_char, CONTINUE);
      input_error++;
      j = OPTION_ERROR;
    }
  }
  else
  {
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
    }
  }
  return (j);
}

#ifdef PHREEQCI_GUI
/* ---------------------------------------------------------------------- */
void
free_spread (void)
/* ---------------------------------------------------------------------- */
{
  int i;
  spread_row_free (g_spread_sheet.heading);
  spread_row_free (g_spread_sheet.units);
  for (i = 0; i < g_spread_sheet.count_rows; i++)
  {
    spread_row_free (g_spread_sheet.rows[i]);
  }
  g_spread_sheet.rows = free_check_null (g_spread_sheet.rows);

  for (i = 0; i < g_spread_sheet.defaults.count_iso; i++)
  {
    g_spread_sheet.defaults.iso[i].name =
      free_check_null (g_spread_sheet.defaults.iso[i].name);
  }
  g_spread_sheet.defaults.iso = free_check_null (g_spread_sheet.defaults.iso);

  g_spread_sheet.defaults.redox =
    free_check_null (g_spread_sheet.defaults.redox);
  g_spread_sheet.defaults.units =
    free_check_null (g_spread_sheet.defaults.units);
}

/* ---------------------------------------------------------------------- */
void
add_row (struct spread_row *spread_row_ptr)
/* ---------------------------------------------------------------------- */
{
  g_spread_sheet.rows =
    (struct spread_row **) PHRQ_realloc (g_spread_sheet.rows,
					 sizeof (struct spread_row *) *
					 (g_spread_sheet.count_rows + 1));
  if (g_spread_sheet.rows == NULL)
  {
    malloc_error ();
  }
  else
  {
    g_spread_sheet.rows[g_spread_sheet.count_rows++] =
      copy_row (spread_row_ptr);
  }
}

/* ---------------------------------------------------------------------- */
struct spread_row *
copy_row (struct spread_row *spread_row_ptr)
/* ---------------------------------------------------------------------- */
{
  int i;
  struct spread_row *new_spread_row_ptr;
/*
 *   Allocate space
 */
  new_spread_row_ptr =
    (struct spread_row *) PHRQ_malloc ((size_t) sizeof (struct spread_row));
  if (new_spread_row_ptr == NULL)
    malloc_error ();
  new_spread_row_ptr->char_vector =
    (char **) PHRQ_malloc ((size_t) spread_row_ptr->count * sizeof (char *));
  if (new_spread_row_ptr->char_vector == NULL)
    malloc_error ();
  new_spread_row_ptr->d_vector =
    (LDBLE *) PHRQ_malloc ((size_t) spread_row_ptr->count * sizeof (LDBLE));
  if (new_spread_row_ptr->d_vector == NULL)
    malloc_error ();
  new_spread_row_ptr->type_vector =
    (int *) PHRQ_malloc ((size_t) spread_row_ptr->count * sizeof (int));
  if (new_spread_row_ptr->type_vector == NULL)
    malloc_error ();

  for (i = 0; i < spread_row_ptr->count; i++)
  {
    new_spread_row_ptr->char_vector[i] =
      string_duplicate (spread_row_ptr->char_vector[i]);
    new_spread_row_ptr->d_vector[i] = spread_row_ptr->d_vector[i];
    new_spread_row_ptr->type_vector[i] = spread_row_ptr->type_vector[i];
  }
  new_spread_row_ptr->count = spread_row_ptr->count;
  new_spread_row_ptr->empty = spread_row_ptr->empty;
  new_spread_row_ptr->number = spread_row_ptr->number;
  new_spread_row_ptr->string = spread_row_ptr->string;

  return new_spread_row_ptr;
}

/* ---------------------------------------------------------------------- */
void
copy_defaults (struct defaults *dest_ptr, struct defaults *src_ptr)
/* ---------------------------------------------------------------------- */
{
  int i;
  dest_ptr->count_iso = src_ptr->count_iso;
  dest_ptr->density = src_ptr->density;
  dest_ptr->iso =
    (struct iso *) PHRQ_malloc (sizeof (struct iso) * src_ptr->count_iso);
  if (dest_ptr->iso == NULL)
  {
    malloc_error ();
  }
  else
  {
    for (i = 0; i < src_ptr->count_iso; i++)
    {
      dest_ptr->iso[i] = src_ptr->iso[i];
      dest_ptr->iso[i].name = string_duplicate (src_ptr->iso[i].name);
    }
  }

  dest_ptr->pe = src_ptr->pe;
  dest_ptr->ph = src_ptr->ph;
  dest_ptr->redox = string_duplicate (src_ptr->redox);
  dest_ptr->temp = src_ptr->temp;
  dest_ptr->units = string_duplicate (src_ptr->units);
  dest_ptr->water = src_ptr->water;
}

#endif
