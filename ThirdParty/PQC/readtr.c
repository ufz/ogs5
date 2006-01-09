#define EXTERNAL extern
#include "global.h"
#include "phqalloc.h"
#include "output.h"
#include "phrqproto.h"

static char const svnid[] =
  "$Id: readtr.c 4 2009-04-21 17:29:29Z delucia $";

static int read_line_LDBLEs (char *next_char, LDBLE ** d, int *count_d,
			     int *count_alloc);

#define OPTION_EOF -1
#define OPTION_KEYWORD -2
#define OPTION_ERROR -3
#define OPTION_DEFAULT -4
#define OPTION_DEFAULT2 -5

/* ---------------------------------------------------------------------- */
int
read_transport (void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads advection and column information
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
  int i, j, l;
  int old_cells, max;
  int count_length, count_disp, count_punch, count_print;
  int count_length_alloc, count_disp_alloc;
  char token[MAX_LENGTH];
  char *description;
  int n_user, n_user_end;
  LDBLE *length, *disp;
  int *punch_temp, *print_temp;
  int return_value, opt, opt_save;
  char *next_char, *next_char_save;
  char file_name[MAX_LENGTH];

  const char *opt_list[] = {
    "cells",			/* 0 */
    "shifts",			/* 1 */
    "print",			/* 2 */
    "selected_output",		/* 3 */
    "bcond",			/* 4 */
    "timest",			/* 5 */
    "diffc",			/* 6 */
    "tempr",			/* 7 */
    "length",			/* 8 */
    "disp",			/* 9 */
    "punch",			/* 10 */
    "stagnant",			/* 11 */
    "bc",			/* 12 */
    "boundary_conditions",	/* 13 */
    "time_step",		/* 14 */
    "temp_retardation_factor",	/* 15 */
    "diffusion_coefficient",	/* 16 */
    "dispersivity",		/* 17 */
    "direction",		/* 18 */
    "temperature_retardation_factor",	/* 19 */
    "print_cells",		/* 20 */
    "selected_cells",		/* 21 */
    "flow_direction",		/* 22 */
    "flow",			/* 23 */
    "lengths",			/* 24 */
    "dispersivities",		/* 25 */
    "dump",			/* 26 */
    "output",			/* 27 */
    "output_frequency",		/* 28 */
    "selected_output_frequency",	/* 29 */
    "punch_cells",		/* 30 */
    "dump_frequency",		/* 31 */
    "dump_restart",		/* 32 */
    "punch_frequency",		/* 33 */
    "print_frequency",		/* 34 */
    "correct_disp",		/* 35 */
    "initial_time",		/* 36 */
    "warning",			/* 37 */
    "warnings",			/* 38 */
    "thermal_diffusion",	/* 39 */
    "multi_d"			/* 40 */
  };
  int count_opt_list = 41;
  if (svnid == NULL)
    fprintf (stderr, " ");

  strcpy (file_name, "phreeqc.dmp");
/*
 *   Initialize
 */
  simul_tr++;
  if (simul_tr == 1)
  {
    correct_disp = FALSE;
    old_cells = 0;
  }
  else
    old_cells = count_cells;
  count_length = count_disp = count_punch = count_print = 0;

  length = (LDBLE *) PHRQ_malloc (sizeof (LDBLE));
  if (length == NULL)
    malloc_error ();

  disp = (LDBLE *) PHRQ_malloc (sizeof (LDBLE));
  if (disp == NULL)
    malloc_error ();

  punch_temp = (int *) PHRQ_malloc (sizeof (int));
  if (punch_temp == NULL)
    malloc_error ();

  print_temp = (int *) PHRQ_malloc (sizeof (int));
  if (print_temp == NULL)
    malloc_error ();

  count_length_alloc = count_disp_alloc = 1;
  transport_start = 1;
/*
 *   Read transport number (not currently used)
 */
  ptr = line;
  read_number_description (ptr, &n_user, &n_user_end, &description);
  description = (char *) free_check_null (description);
/*
 *   Set use data to last read
 */
  use.trans_in = TRUE;
/*
 *   Read lines
 */
  opt_save = OPTION_DEFAULT;
  return_value = UNKNOWN;
  for (;;)
  {
    opt = get_option (opt_list, count_opt_list, &next_char);
    if (opt == OPTION_DEFAULT)
      opt = opt_save;
    switch (opt)
    {
    case OPTION_EOF:		/* end of file */
      return_value = EOF;
      break;
    case OPTION_KEYWORD:	/* keyword */
      return_value = KEYWORD;
      break;
    case OPTION_ERROR:
    case OPTION_DEFAULT:
      input_error++;
      error_msg ("Unknown input in TRANSPORT keyword.", CONTINUE);
      error_msg (line_save, CONTINUE);
      break;
    case 0:			/* cells */
      sscanf (next_char, "%d", &count_cells);
      opt_save = OPTION_DEFAULT;
      break;
    case 1:			/* shifts */
      if (copy_token (token, &next_char, &l) == DIGIT)
	sscanf (token, "%d", &count_shifts);
      else
      {
	warning_msg ("Expected the number of shifts. One shift is assumed.");
	count_shifts = 1;
      }
      j = copy_token (token, &next_char, &l);
      if (j != EMPTY)
      {
	if (j == DIGIT)
	  sscanf (token, "%d", &ishift);
	else
	{
	  input_error++;
	  error_msg
	    ("Expected shift direction, -1, 0, 1. Use -direction instead.",
	     CONTINUE);
	  ishift = 1;
	}
      }
      opt_save = OPTION_DEFAULT;
      break;
    case 2:			/* print */
    case 20:			/* print_cells */
      print_temp =
	read_list_ints_range (&next_char, &count_print, TRUE, print_temp);
#ifdef SKIP
      while (copy_token (token, &next_char, &l) == DIGIT)
      {
	sscanf (token, "%d", &l);
	print_temp =
	  PHRQ_realloc (print_temp,
			(size_t) (count_print + 1) * sizeof (LDBLE));
	if (print_temp == NULL)
	  malloc_error ();
	print_temp[count_print] = l;
	count_print++;
      }
#endif
      opt_save = 2;
      break;
    case 3:			/* selected_output */
    case 29:			/* selected_output_frequency */
    case 33:			/* punch_frequency */
      sscanf (next_char, "%d", &punch_modulus);
      opt_save = OPTION_DEFAULT;
      if (punch_modulus <= 0) 
      {
	sprintf (error_string,
		 "Punch frequency must be greater than 0. Frequency set to 1000.");
	warning_msg (error_string);
	punch_modulus = 1000;
      }

      break;
    case 4:			/* bcond */
    case 12:			/* bc */
    case 13:			/* boundary_conditions */
      /* first cell boundary condition */
      i = copy_token (token, &next_char, &l);
      str_tolower (token);
      if (i == DIGIT)
      {
	sscanf (token, "%d", &bcon_first);
	if (bcon_first < 1 || bcon_first > 3)
	{
	  input_error++;
	  error_msg
	    ("Expected boundary condition to be 'constant' (1), 'closed' (2) , or 'flux' (3).",
	     CONTINUE);
	}
      }
      else if (i == EMPTY)
	bcon_first = 3;
      else if (strstr (token, "co") == token)
	bcon_first = 1;
      else if (strstr (token, "cl") == token)
	bcon_first = 2;
      else if (strstr (token, "f") == token)
	bcon_first = 3;
      else
      {
	input_error++;
	error_msg
	  ("Expected boundary condition to be 'constant', 'closed', or 'flux'.",
	   CONTINUE);
      }

      /* last cell boundary condition */
      i = copy_token (token, &next_char, &l);
      str_tolower (token);
      if (i == DIGIT)
      {
	sscanf (token, "%d", &bcon_last);
	if (bcon_last < 1 || bcon_last > 3)
	{
	  input_error++;
	  error_msg
	    ("Expected boundary condition to be 'constant' (1), 'closed' (2) , or 'flux' (3).",
	     CONTINUE);
	}
      }
      else if (i == EMPTY)
	bcon_last = 3;
      else if (strstr (token, "co") == token)
	bcon_last = 1;
      else if (strstr (token, "cl") == token)
	bcon_last = 2;
      else if (strstr (token, "f") == token)
	bcon_last = 3;
      else
      {
	input_error++;
	error_msg
	  ("Expected boundary condition to be 'constant', 'closed', or 'flux'.",
	   CONTINUE);
      }
      opt_save = OPTION_DEFAULT;
      break;
    case 5:			/* timest */
    case 14:			/* time_step */
      if (copy_token (token, &next_char, &l) == DIGIT)
	sscanf (token, SCANFORMAT, &timest);
      if (copy_token (token, &next_char, &l) == DIGIT)
	sscanf (token, SCANFORMAT, &mcd_substeps);
      if (mcd_substeps < 1)
      {
	mcd_substeps = 1.0;
	warning_msg ("Substep factor in MCD must be >= 1.0\n"
		     "mcd_substeps = 1.0 assumed.");
      }
      opt_save = OPTION_DEFAULT;
      break;
    case 6:			/* diffc */
    case 16:			/* diffusion_coefficient */
      sscanf (next_char, SCANFORMAT, &diffc);
      opt_save = OPTION_DEFAULT;
      break;
    case 7:			/* tempr */
    case 15:			/* temp_retardation_factor */
    case 19:			/* temperature_retardation_factor */
    case 39:			/* thermal_diffusion */
      if (copy_token (token, &next_char, &l) == DIGIT)
	sscanf (token, SCANFORMAT, &tempr);
      if (tempr < 1)
      {
	tempr = 1;
	warning_msg ("Temperature retardation factor < 1 is not possible.\n"
		     "Temperature retardation factor = 1 assumed.");
      }
      j = copy_token (token, &next_char, &l);
      if (j == DIGIT)
	sscanf (token, SCANFORMAT, &heat_diffc);
      opt_save = OPTION_DEFAULT;
      break;
    case 8:			/* length */
    case 24:			/* lengths */
      if (read_line_LDBLEs
	  (next_char, &length, &count_length, &count_length_alloc) == ERROR)
      {
	input_error++;
	error_msg ("Reading lengths in TRANSPORT keyword.\n", CONTINUE);
      }
      opt_save = 8;
      break;
    case 9:			/* disp */
    case 17:			/* dispersivity */
    case 25:			/* dispersivities */
      if (read_line_LDBLEs (next_char, &disp, &count_disp, &count_disp_alloc)
	  == ERROR)
      {
	input_error++;
	error_msg ("Reading dispersivities in TRANSPORT keyword.\n",
		   CONTINUE);
      }
      opt_save = 9;
      break;
    case 10:			/* punch */
    case 21:			/* selected_cells */
    case 30:			/* punch_cells */
      punch_temp =
	read_list_ints_range (&next_char, &count_punch, TRUE, punch_temp);
#ifdef SKIP
      while (copy_token (token, &next_char, &l) == DIGIT)
      {
	sscanf (token, "%d", &l);
	punch_temp =
	  PHRQ_realloc (punch_temp,
			(size_t) (count_punch + 1) * sizeof (LDBLE));
	if (punch_temp == NULL)
	  malloc_error ();
	punch_temp[count_punch] = l;
	count_punch++;
      }
#endif
      opt_save = 10;
      break;
    case 11:			/* stagnant */
      if (copy_token (token, &next_char, &l) != EMPTY)
      {
	/* exchange factor */
	if (sscanf (token, "%d", &(stag_data->count_stag)) != 1)
	{
	  input_error++;
	  sprintf (error_string, "Expecting number of stagnant layers.");
	  error_msg (error_string, CONTINUE);
	  break;
	}

	/* exchange factor */
	j = copy_token (token, &next_char, &l);
	if (j != EMPTY)
	{
	  if (sscanf (token, SCANFORMAT, &(stag_data->exch_f)) != 1)
	  {
	    input_error++;
	    sprintf (error_string,
		     "Expecting exchange factor for stagnant layers.");
	    error_msg (error_string, CONTINUE);
	    break;
	  }
	  copy_token (token, &next_char, &l);
	  if (sscanf (token, SCANFORMAT, &(stag_data->th_m)) != 1)
	  {
	    input_error++;
	    sprintf (error_string, "Expecting porosity in the mobile zone.");
	    error_msg (error_string, CONTINUE);
	    break;
	  }
	  copy_token (token, &next_char, &l);
	  if (sscanf (token, SCANFORMAT, &(stag_data->th_im)) != 1)
	  {
	    input_error++;
	    sprintf (error_string,
		     "Expecting porosity in the immobile zone.");
	    error_msg (error_string, CONTINUE);
	    break;
	  }
	}
      }
      opt_save = OPTION_DEFAULT;
      break;
    case 18:			/* direction */
    case 22:			/* flow_direction */
    case 23:			/* flow */
      copy_token (token, &next_char, &l);
      str_tolower (token);
      if (strstr (token, "f") == token)
	ishift = 1;
      else if (strstr (token, "b") == token)
	ishift = -1;
      else if (strstr (token, "d") == token)
	ishift = 0;
      else if (strstr (token, "n") == token)
	ishift = 0;
      else
      {
	input_error++;
	error_msg
	  ("Expected flow direction to be 'forward', 'back', or 'no_flow'.",
	   CONTINUE);
      }
      opt_save = OPTION_DEFAULT;
      break;
    case 26:			/* dump */
      dump_in = TRUE;
      next_char_save = next_char;
      if (copy_token (file_name, &next_char, &l) == EMPTY)
	strcpy (file_name, "phreeqc.dmp");
      else
      {
	string_trim (next_char_save);
	strcpy (file_name, next_char_save);
      }
#ifdef SKIP
      /* Can not define dump_modulus an transport_start here */
      if (copy_token (token, &next_char, &l) == DIGIT)
      {
	sscanf (token, "%d", &dump_modulus);
      }
      if (copy_token (token, &next_char, &l) == DIGIT)
      {
	sscanf (token, "%d", &transport_start);
      }
#endif
      opt_save = OPTION_DEFAULT;
      break;
    case 27:			/* output */
    case 28:			/* output_frequency */
    case 34:			/* print_frequency */
      sscanf (next_char, "%d", &print_modulus);
      opt_save = OPTION_DEFAULT;
      if (print_modulus <= 0) 
      {
	sprintf (error_string,
		 "Print frequency must be greater than 0. Frequency set to 1000.");
	warning_msg (error_string);
	print_modulus = 1000;
      }
      break;
    case 31:			/* dump_frequency */
      dump_in = TRUE;
      if (copy_token (token, &next_char, &l) == DIGIT)
	sscanf (token, "%d", &dump_modulus);
      else
      {
	warning_msg ("Expected integer value for dump_frequency.");
	dump_modulus = 0;
      }
      opt_save = OPTION_DEFAULT;
      break;
    case 32:			/* dump_restart */
      dump_in = TRUE;
      if (copy_token (token, &next_char, &l) == DIGIT)
	sscanf (token, "%d", &transport_start);
      else
      {
	warning_msg
	  ("Expected shift number to start calculations, 1 will be used.");
	transport_start = 1;
      }
      opt_save = OPTION_DEFAULT;
      break;
    case 35:			/* correct_dispersion */
      correct_disp = get_true_false (next_char, TRUE);
      opt_save = OPTION_DEFAULT;
      break;
    case 36:			/* initial_time */
      sscanf (next_char, SCANFORMAT, &initial_total_time);
      opt_save = OPTION_DEFAULT;
      break;
    case 37:			/* warning */
    case 38:			/* warnings */
      transport_warnings = get_true_false (next_char, TRUE);
      break;
    case 40:			/* multicomponent diffusion */
      copy_token (token, &next_char, &l);
      str_tolower (token);
      if (strstr (token, "f") == token)
	multi_Dflag = 0;
      else if (strstr (token, "t") == token)
	multi_Dflag = 1;
      else
      {
	input_error++;
	error_msg
	  ("Expected multicomponent diffusion flag: 'true' or 'false'.",
	   CONTINUE);
      }
      default_Dw = 1e-9;
      multi_Dpor = 0.3;
      multi_Dpor_lim = 0.0;
      multi_Dn = 1.0;
      if (copy_token (token, &next_char, &l) == EMPTY)
	break;
      else
      {
	/* default species diffusion coeff */
	if (sscanf (token, SCANFORMAT, &default_Dw) != 1)
	{
	  input_error++;
	  sprintf (error_string,
		   "Expected default species diffusion coefficient in water at 25oC, m2/s.");
	  error_msg (error_string, CONTINUE);
	  break;
	}
      }
      if (copy_token (token, &next_char, &l) == EMPTY)
	break;
      else
      {
	/* porosity */
	if (sscanf (token, SCANFORMAT, &multi_Dpor) != 1)
	{
	  input_error++;
	  sprintf (error_string,
		   "Expected porosity to calculate diffusion coefficient.");
	  error_msg (error_string, CONTINUE);
	  break;
	}
      }
      if (copy_token (token, &next_char, &l) == EMPTY)
	break;
      else
      {
	/* porosity */
	if (sscanf (token, SCANFORMAT, &multi_Dpor_lim) != 1)
	{
	  input_error++;
	  sprintf (error_string,
		   "Expected porosity limit for diffusive transport.");
	  error_msg (error_string, CONTINUE);
	  break;
	}
      }
      if (copy_token (token, &next_char, &l) == EMPTY)
	break;
      else
      {
	if (sscanf (token, SCANFORMAT, &multi_Dn) != 1)
	{
	  input_error++;
	  sprintf (error_string,
		   "Expected exponent for porosity reduction of diffusion coefficient (Dp = Dw * (por)^n).");
	  error_msg (error_string, CONTINUE);
	  break;
	}
      }
      opt_save = OPTION_DEFAULT;
      break;
    }
    if (return_value == EOF || return_value == KEYWORD)
      break;
  }
/*
 *   Determine number of cells
 */
  max = count_cells;
  if (count_length > max)
    max = count_length;
  if (count_disp > max)
    max = count_disp;
  if (max > count_cells)
  {
    if (max == count_length)
    {
      sprintf (token,
	       "Number of cells is increased to number of 'lengths' %d.",
	       count_length);
      warning_msg (token);
    }
    else
    {
      sprintf (token,
	       "Number of cells is increased to number of dispersivities %d.",
	       count_disp);
      warning_msg (token);
    }
  }
/*
 *   Allocate space for cell_data
 */
  cell_data =
    (struct cell_data *) PHRQ_realloc (cell_data,
				       (size_t) (max *
						 (1 + stag_data->count_stag) +
						 1) *
				       sizeof (struct cell_data));
  if (cell_data == NULL)
    malloc_error ();
/*
 *   Fill in data for lengths
 */
  if (count_length == 0)
  {
    if (old_cells < max)
    {
      sprintf (error_string,
	       "No cell-lengths were read; length = 1 m assumed.");
      warning_msg (error_string);
      for (i = 0; i < max; i++)
	cell_data[i].length = 1.0;
    }
  }
  else
  {
    for (i = 0; i < count_length; i++)
    {
      cell_data[i].length = length[i];
    }
    if (max > count_length)
    {
      sprintf (error_string,
	       "Cell-lengths were read for %d cells. Last value is used till cell %d.",
	       count_length, max);
      warning_msg (error_string);
      for (i = count_length - 1; i < max; i++)
	cell_data[i].length = length[count_length - 1];
    }
  }
  cell_data[0].mid_cell_x = cell_data[0].length / 2;
  for (i = 1; i < max; i++)
  {
    cell_data[i].mid_cell_x = cell_data[i - 1].mid_cell_x +
      (cell_data[i - 1].length + cell_data[i].length) / 2;
  }
  cell_data[max].mid_cell_x =
    cell_data[max - 1].mid_cell_x + cell_data[max - 1].length;
/*
 *   Fill in data for dispersivities
 */
  if (count_disp == 0)
  {
    if (old_cells < max)
    {
      sprintf (error_string,
	       "No dispersivities were read; disp = 0 assumed.");
      warning_msg (error_string);
      for (i = 0; i < max; i++)
	cell_data[i].disp = 0.0;
    }
  }
  else
  {
    for (i = 0; i < count_disp; i++)
      cell_data[i].disp = disp[i];
    if (max > count_disp)
    {
      sprintf (error_string,
	       "Dispersivities were read for %d cells. Last value is used till cell %d.",
	       count_disp, max);
      warning_msg (error_string);
      for (i = count_disp - 1; i < max; i++)
	cell_data[i].disp = disp[count_disp - 1];
    }
  }
  count_cells = max;
/*
 *  Account for stagnant cells
 */
  if (stag_data->count_stag > 0)
  {
    max = count_cells * (1 + stag_data->count_stag) + 1;
    for (i = 0; i < count_cells; i++)
    {
      for (l = 1; l <= stag_data->count_stag; l++)
	cell_data[i + 1 + l * count_cells].mid_cell_x =
	  cell_data[i].mid_cell_x;
    }
  }
/*
 *   Fill in data for punch
 */
  if (count_punch != 0)
  {
    for (i = 0; i < max; i++)
      cell_data[i].punch = FALSE;
    for (i = 0; i < count_punch; i++)
    {
      if (punch_temp[i] > max || punch_temp[i] < 1)
      {
	sprintf (error_string,
		 "Cell number for punch is out of range, %d. Request ignored.",
		 punch_temp[i]);
	warning_msg (error_string);
      }
      else
	cell_data[punch_temp[i] - 1].punch = TRUE;
    }
  }
  else if (simul_tr == 1)
    for (i = 0; i < max; i++)
      cell_data[i].punch = TRUE;
/*
 *   Fill in data for print
 */
  if (count_print != 0)
  {
    for (i = 0; i < max; i++)
      cell_data[i].print = FALSE;
    for (i = 0; i < count_print; i++)
    {
      if (print_temp[i] > max || print_temp[i] < 1)
      {
	sprintf (error_string,
		 "Cell number for print is out of range, %d. Request ignored.",
		 print_temp[i]);
	warning_msg (error_string);
      }
      else
	cell_data[print_temp[i] - 1].print = TRUE;
    }
  }
  else if (simul_tr == 1)
    for (i = 0; i < max; i++)
      cell_data[i].print = TRUE;
/*
 *   Fill in porosities
 */
  for (i = 0; i < max; i++)
    cell_data[i].por = multi_Dpor;
/*
 *   Calculate dump_modulus
 */
  if (dump_in == TRUE)
  {
    if (dump_modulus == 0)
    {
      warning_msg
	("Expected dump_modulus. Value of 'shifts/2' will be used.");
      dump_modulus = count_shifts / 2;
      if (dump_modulus == 0)
	dump_modulus = 1;
    }
    if (transport_start > count_shifts)
    {
      input_error++;
      sprintf (error_string,
	       "Starting shift for transport, %d, is greater than number of shifts, %d.",
	       transport_start, count_shifts);
      error_msg (error_string, CONTINUE);
    }
  }
/*
 *  Check boundary conditions
 */
  if ((ishift != 0) && ((bcon_first == 2) || (bcon_last == 2)))
  {
    warning_msg
      ("Boundary condition = 'closed' not possible with advective transport.\n\t Boundary condition = 'flux' assumed.");
    if (bcon_first == 2)
      bcon_first = 3;
    if (bcon_last == 2)
      bcon_last = 3;
  }
/*
 *  Warn for multicomponent diffusion (the effect is small)...
 */
#ifdef SKIP
  if (multi_Dflag == TRUE && count_disp > 0)
  {
    sprintf (error_string,
	     "A model with multicomponent diffusion and dispersivity > 0\n\t may show ill-defined, too high dispersion...");
    warning_msg (error_string);
  }
#endif
/*
 *  Retain data from previous run
 */
  if (simul_tr > 1)
  {
    if ((count_length == 0) && (count_disp == 0))
      dup_print ("Column data retained from former run", TRUE);
  }
/*
 *  Check heat_diffc
 */
  if (heat_diffc < 0)
    heat_diffc = diffc;
  else if (stag_data->count_stag == 1)
  {
    if (stag_data->exch_f > 0)
    {
      if (diffc <= 0 && heat_diffc > 0)
      {
	input_error++;
	sprintf (error_string,
		 "Must enter diffusion coefficient (-diffc) when modeling thermal diffusion.");
	error_msg (error_string, CONTINUE);
      }
      else if (heat_diffc > diffc)
      {
	sprintf (error_string,
		 "Thermal diffusion is calculated assuming exchange factor was for\n\t effective (non-thermal) diffusion coefficient = %e.",
		 (double) diffc);
	warning_msg (error_string);
      }
    }
    else
    {
      if (heat_diffc > diffc)
      {
	input_error++;
	sprintf (error_string,
		 "Must enter value for mobile/stagnant exchange factor when modeling thermal diffusion.");
	error_msg (error_string, CONTINUE);
      }
    }
  }
  else if (stag_data->count_stag > 1 && heat_diffc > diffc)
  {
    input_error++;
    sprintf (error_string,
	     "Only one stagant layer permitted (-stag) when modeling thermal diffusion.");
    error_msg (error_string, CONTINUE);
  }
/*
 *   free storage for length, disp, punch
 */
  length = (LDBLE *) free_check_null (length);
  disp = (LDBLE *) free_check_null (disp);
  punch_temp = (int *) free_check_null (punch_temp);
  print_temp = (int *) free_check_null (print_temp);

  if (dump_in == TRUE)
  {
    if (output_open (OUTPUT_DUMP, file_name) != OK)
    {
      sprintf (error_string, "Can't open file, %s.", file_name);
      error_msg (error_string, CONTINUE);
      input_error++;
    }
  }
  return (return_value);
}

/* ---------------------------------------------------------------------- */
int
read_line_LDBLEs (char *next_char, LDBLE ** d, int *count_d, int *count_alloc)
/* ---------------------------------------------------------------------- */
{
  int i, j, l, n;
  LDBLE value;
  char token[MAX_LENGTH];

  for (;;)
  {
    j = copy_token (token, &next_char, &l);
    if (j == EMPTY)
      break;
    if (j != DIGIT)
      return (ERROR);
    if (replace ("*", " ", token) == TRUE)
    {
      if (sscanf (token, "%d" SCANFORMAT, &n, &value) != 2)
	return (ERROR);
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
	break;
    }
    for (i = 0; i < n; i++)
      (*d)[(*count_d) + i] = value;
    *count_d += n;
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
dump (void)
/* ---------------------------------------------------------------------- */
{
/*
 * dumps solution compositions to file
 */

  int i, j, k, l;

  if (dump_in == FALSE || pr.dump == FALSE)
    return (OK);
  output_rewind (OUTPUT_DUMP);
  output_msg (OUTPUT_DUMP,
	      "# Dumpfile\n# Transport simulation %d.  Shift %d.\n#\n",
	      simul_tr, transport_step);

  for (k = 0; k <= 1 + (1 + stag_data->count_stag) * count_cells; k++)
  {
    dump_reaction (k);
    dump_kinetics (k);
    output_msg (OUTPUT_DUMP, "END\n");
    dump_solution (k);
    dump_pp_assemblage (k);
    dump_exchange (k);
    dump_surface (k);
    dump_gas_phase (k);
    dump_s_s_assemblage (k);
    output_msg (OUTPUT_DUMP, "END\n");
  }

  output_msg (OUTPUT_DUMP, "KNOBS\n");
  output_msg (OUTPUT_DUMP, "\t-iter%15d\n", itmax);
  output_msg (OUTPUT_DUMP, "\t-tol %15.3e\n", (double) ineq_tol);
  output_msg (OUTPUT_DUMP, "\t-step%15.3e\n", (double) step_size);
  output_msg (OUTPUT_DUMP, "\t-pe_s%15.3e\n", (double) pe_step_size);
  output_msg (OUTPUT_DUMP, "\t-diag      ");
  if (diagonal_scale == TRUE)
    output_msg (OUTPUT_DUMP, "true\n");
  else
    output_msg (OUTPUT_DUMP, "false\n");

  output_msg (OUTPUT_DUMP, "SELECTED_OUTPUT\n");
  output_msg (OUTPUT_DUMP, "\t-file  %-15s\n", "sel_o$$$.prn");
  if (punch.count_totals != 0)
  {
    output_msg (OUTPUT_DUMP, "\t-tot ");
    for (i = 0; i < punch.count_totals; i++)
      output_msg (OUTPUT_DUMP, "  %s", punch.totals[i].name);
    output_msg (OUTPUT_DUMP, "\n");
  }
  if (punch.count_molalities != 0)
  {
    output_msg (OUTPUT_DUMP, "\t-mol ");
    for (i = 0; i < punch.count_molalities; i++)
      output_msg (OUTPUT_DUMP, "  %s", punch.molalities[i].name);
    output_msg (OUTPUT_DUMP, "\n");
  }
  if (punch.count_activities != 0)
  {
    output_msg (OUTPUT_DUMP, "\t-act ");
    for (i = 0; i < punch.count_activities; i++)
      output_msg (OUTPUT_DUMP, "  %s", punch.activities[i].name);
    output_msg (OUTPUT_DUMP, "\n");
  }
  if (punch.count_pure_phases != 0)
  {
    output_msg (OUTPUT_DUMP, "\t-equ ");
    for (i = 0; i < punch.count_pure_phases; i++)
      output_msg (OUTPUT_DUMP, "  %s", punch.pure_phases[i].name);
    output_msg (OUTPUT_DUMP, "\n");
  }
  if (punch.count_si != 0)
  {
    output_msg (OUTPUT_DUMP, "\t-si ");
    for (i = 0; i < punch.count_si; i++)
      output_msg (OUTPUT_DUMP, "  %s", punch.si[i].name);
    output_msg (OUTPUT_DUMP, "\n");
  }
  if (punch.count_gases != 0)
  {
    output_msg (OUTPUT_DUMP, "\t-gas ");
    for (i = 0; i < punch.count_gases; i++)
      output_msg (OUTPUT_DUMP, "  %s", punch.gases[i].name);
    output_msg (OUTPUT_DUMP, "\n");
  }
  if (punch.count_s_s != 0)
  {
    output_msg (OUTPUT_DUMP, "\t-solid_solutions ");
    for (i = 0; i < punch.count_s_s; i++)
      output_msg (OUTPUT_DUMP, "  %s", punch.s_s[i].name);
    output_msg (OUTPUT_DUMP, "\n");
  }
  if (punch.count_kinetics != 0)
  {
    output_msg (OUTPUT_DUMP, "\t-kin ");
    for (i = 0; i < punch.count_kinetics; i++)
      output_msg (OUTPUT_DUMP, "  %s", punch.kinetics[i].name);
    output_msg (OUTPUT_DUMP, "\n");
  }
  output_msg (OUTPUT_DUMP, "TRANSPORT\n");
  output_msg (OUTPUT_DUMP, "\t-cells %6d\n", count_cells);
  output_msg (OUTPUT_DUMP, "\t-shifts%6d%6d\n", count_shifts, ishift);
  output_msg (OUTPUT_DUMP, "\t-output_frequency %6d\n", print_modulus);
  output_msg (OUTPUT_DUMP, "\t-selected_output_frequency %6d\n",
	      punch_modulus);
  output_msg (OUTPUT_DUMP, "\t-bcon  %6d%6d\n", bcon_first, bcon_last);
  output_msg (OUTPUT_DUMP, "\t-timest %13.5e\n", (double) timest);
  output_msg (OUTPUT_DUMP, "\t-diffc  %13.5e\n", (double) diffc);
  output_msg (OUTPUT_DUMP, "\t-tempr  %13.5e\n", (double) tempr);
  if (correct_disp == TRUE)
    output_msg (OUTPUT_DUMP, "\t-correct_disp %s\n", "True");
  else
    output_msg (OUTPUT_DUMP, "\t-correct_disp %s\n", "False");
  output_msg (OUTPUT_DUMP, "\t-length\n");
  for (i = 0; i < count_cells; i++)
  {
    output_msg (OUTPUT_DUMP, "%12.3e", (double) cell_data[i].length);
    if (i > 0 && (i % 8) == 0)
      output_msg (OUTPUT_DUMP, "\n");
  }
  output_msg (OUTPUT_DUMP, "\n");
  output_msg (OUTPUT_DUMP, "\t-disp\n");
  for (i = 0; i < count_cells; i++)
  {
    output_msg (OUTPUT_DUMP, "%12.3e", (double) cell_data[i].disp);
    if (i > 0 && (i % 8) == 0)
      output_msg (OUTPUT_DUMP, "\n");
  }
  output_msg (OUTPUT_DUMP, "\n");
  output_msg (OUTPUT_DUMP, "\t-punch_cells");
  if (stag_data->count_stag > 0)
    j = 1 + (1 + stag_data->count_stag) * count_cells;
  else
    j = count_cells;
  l = 0;
  for (i = 0; i < j; i++)
  {
    if (cell_data[i].punch != TRUE)
      continue;
    output_msg (OUTPUT_DUMP, "  %d", i + 1);
    l++;
    if ((l % 20) == 0)
      output_msg (OUTPUT_DUMP, "\n");
  }
  output_msg (OUTPUT_DUMP, "\n");
  output_msg (OUTPUT_DUMP, "\t-print_cells");
  if (stag_data->count_stag > 0)
    j = 1 + (1 + stag_data->count_stag) * count_cells;
  else
    j = count_cells;
  l = 0;
  for (i = 0; i < j; i++)
  {
    if (cell_data[i].print != TRUE)
      continue;
    output_msg (OUTPUT_DUMP, "  %d", i + 1);
    l++;
    if ((l % 20) == 0)
      output_msg (OUTPUT_DUMP, "\n");
  }
  output_msg (OUTPUT_DUMP, "\n");
  output_msg (OUTPUT_DUMP, "\t-dump            $$$.dmp\n");
  output_msg (OUTPUT_DUMP, "\t-dump_frequency  %d\n", dump_modulus);
  output_msg (OUTPUT_DUMP, "\t-dump_restart    %d\n", transport_step + 1);

  output_msg (OUTPUT_DUMP, "END\n");
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
dump_exchange (int k)
/* ---------------------------------------------------------------------- */
{
/*
 *   Print moles of each exchange species
 */
  int i, j, n;
  struct exchange *exchange_ptr;

  exchange_ptr = exchange_bsearch (k, &n);
  if (exchange_ptr == NULL)
    return (OK);

  output_msg (OUTPUT_DUMP, "EXCHANGE %d\n", k);
#ifdef SKIP
  for (i = 0; i < exchange[n].count_comps; i++)
  {
    output_msg (OUTPUT_DUMP, "\t%-15s", exchange[n].comps[i].formula);
    if (exchange[n].comps[i].phase_name != NULL)
    {
      output_msg (OUTPUT_DUMP, "\t%-15s%15.6e\n",
		  exchange[n].comps[i].phase_name,
		  exchange[n].comps[i].phase_proportion);
    }
    else
    {
      for (j = 0; exchange[n].comps[i].totals[j].elt != NULL; j++)
      {
	if (strcmp
	    (exchange_ptr->comps[i].totals[j].elt->name,
	     exchange[n].comps[i].formula) == NULL)
	{
	  output_msg (OUTPUT_DUMP, "%15.6e\n",
		      exchange[n].comps[i].totals[j].coef);
	}
      }
    }
  }
#endif
  for (i = 0; i < exchange[n].count_comps; i++)
  {
    if (exchange[n].comps[i].phase_name != NULL)
    {
      output_msg (OUTPUT_DUMP, "\t%-15s", exchange[n].comps[i].formula);
      output_msg (OUTPUT_DUMP, "\t%-15s%15.6e\n",
		  exchange[n].comps[i].phase_name,
		  (double) exchange[n].comps[i].phase_proportion);
    }
    else
    {
      for (j = 0; exchange[n].comps[i].totals[j].elt != NULL; j++)
      {
/*				if (exchange_ptr->comps[i].totals[j].elt->master->type != EX) continue; */
	output_msg (OUTPUT_DUMP, "\t%-15s",
		    exchange[n].comps[i].totals[j].elt->name);
	output_msg (OUTPUT_DUMP, "%15.6e",
		    (double) exchange[n].comps[i].totals[j].coef);

	if (exchange[n].comps[i].rate_name != NULL)
	{
	  output_msg (OUTPUT_DUMP, "\t%-15s%15.6e",
		      exchange[n].comps[i].rate_name,
		      (double) exchange[n].comps[i].phase_proportion);
	}
#ifdef SKIP
	for (l = 0; l < count_kin_exch; l++)
	{
	  if (strcmp_nocase
	      (kin_exch[l].exch_name,
	       exchange[n].comps[i].totals[j].elt->name) != 0)
	    continue;
	  output_msg (OUTPUT_DUMP, "\t%-15s%15.6e", kin_exch[l].phase_name,
		      kin_exch[l].phase_proportion);
	}
#endif
	output_msg (OUTPUT_DUMP, "\n");
      }
    }
  }
  output_msg (OUTPUT_DUMP, "\t-equil  %d\n", k);
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
dump_pp_assemblage (int k)
/* ---------------------------------------------------------------------- */
{
/*
 *   Prints saturation indices and masses of pure_phases in pp_assemblage
 */
  int j, n;

  if (pp_assemblage_search (k, &n) == NULL)
    return (OK);

  output_msg (OUTPUT_DUMP, "EQUILIBRIUM_PHASES %d\n", k);
  for (j = 0; j < pp_assemblage[n].count_comps; j++)
  {
    output_msg (OUTPUT_DUMP, "\t%-15s%15.6e",
		pp_assemblage[n].pure_phases[j].name,
		(double) pp_assemblage[n].pure_phases[j].si);

    if (pp_assemblage[n].pure_phases[j].add_formula == NULL)
    {
      output_msg (OUTPUT_DUMP, "%15.6e",
		  (double) pp_assemblage[n].pure_phases[j].moles);
    }
    else
    {
      output_msg (OUTPUT_DUMP, "\t%-15s%15.6e",
		  pp_assemblage[n].pure_phases[j].add_formula,
		  (double) pp_assemblage[n].pure_phases[j].moles);
    }
    output_msg (OUTPUT_DUMP, "\n");
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
dump_reaction (int k)
/* ---------------------------------------------------------------------- */
{
/*
 *      Dumps reaction data
 */
  int i;
  struct irrev *irrev_ptr;

  if ((irrev_ptr = irrev_search (k, &i)) == NULL)
    return (OK);

  output_msg (OUTPUT_DUMP, "REACTION %d\n", k);
  for (i = 0; i < irrev_ptr->count_list; i++)
  {
    output_msg (OUTPUT_DUMP, "\t%-15s%15.6e\n",
		irrev_ptr->list[i].name, (double) irrev_ptr->list[i].coef);
  }
  if (irrev_ptr->count_steps >= 1)
  {
    for (i = 0; i < irrev_ptr->count_steps; i++)
      output_msg (OUTPUT_DUMP, "%15.6e", (double) irrev_ptr->steps[i]);
  }
  else
  {
    output_msg (OUTPUT_DUMP, "\t%15.6e in %d step",
		(double) irrev_ptr->steps[0], -irrev_ptr->count_steps);
  }
  output_msg (OUTPUT_DUMP, "\n");
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
dump_surface (int k)
/* ---------------------------------------------------------------------- */
{
/*
 *   Print moles of each surface master species
 */
  int i, j, l, n;
  struct surface *surface_ptr;

  if ((surface_ptr = surface_bsearch (k, &n)) == NULL)
    return (OK);

  output_msg (OUTPUT_DUMP, "SURFACE %d\n", k);
  for (i = 0; i < surface[n].count_comps; i++)
  {
    for (j = 0; surface[n].comps[i].totals[j].elt != NULL; j++)
    {
      if (surface_ptr->comps[i].totals[j].elt->master->type != SURF)
	continue;
      output_msg (OUTPUT_DUMP, "\t%-13s",
		  surface[n].comps[i].totals[j].elt->name);

      if (surface[n].comps[i].phase_name != NULL)
      {
	output_msg (OUTPUT_DUMP, "\t%-13s  equi  %13.5e",
		    surface[n].comps[i].phase_name,
		    (double) surface[n].comps[i].phase_proportion);
      }
      else if (surface[n].comps[i].rate_name != NULL)
      {
	output_msg (OUTPUT_DUMP, "\t%-13s  kine  %13.5e",
		    surface[n].comps[i].rate_name,
		    (double) surface[n].comps[i].phase_proportion);
      }
      else
      {
	output_msg (OUTPUT_DUMP, "\t%13.5e",
		    (double) surface[n].comps[i].totals[j].coef);
#ifdef SKIP
	for (l = 0; l < count_kin_surf; l++)
	{
	  if (strcmp_nocase
	      (kin_surf[l].surf_name,
	       surface[n].comps[i].totals[j].elt->name) != 0)
	    continue;
	  output_msg (OUTPUT_DUMP, "\t%-13s %13.5e", kin_surf[l].phase_name,
		      kin_surf[l].phase_proportion);
	}
#endif
      }
      /*if (surface[n].edl == TRUE) { */
      if (surface[n].type == DDL)
      {
	l = surface[n].comps[i].charge;
	output_msg (OUTPUT_DUMP, " %13.5e",
		    (double) surface[n].charge[l].specific_area);
	if (surface[n].comps[i].phase_name == NULL
	    && surface[n].comps[i].rate_name == NULL)
	{
	  output_msg (OUTPUT_DUMP, " %13.5e\n",
		      (double) surface[n].charge[l].grams);
	}
	else
	  output_msg (OUTPUT_DUMP, "\n");
      }
      else
	output_msg (OUTPUT_DUMP, "\n");
    }
  }
  if (surface[n].type == NO_EDL)
  {
    output_msg (OUTPUT_DUMP, "\t-no_edl\n");
  }
  else if (surface[n].type == DDL)
  {
  }
  else if (surface[n].type == CD_MUSIC)
  {
    output_msg (OUTPUT_DUMP, "\t-cd_music\n");
  }
  if (surface[n].dl_type == DONNAN_DL)
  {
    output_msg (OUTPUT_DUMP, "\t-donnan\t%13.5e\n",
		(double) surface[n].thickness);
    /*} else if (surface[n].edl == FALSE) { */
  }
  else if (surface[n].dl_type == BORKOVEK_DL)
  {
    output_msg (OUTPUT_DUMP, "\t-diffuse_layer\t%13.5e\n",
		(double) surface[n].thickness);
  }


  output_msg (OUTPUT_DUMP, "\t-equil  %d\n", k);
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
dump_gas_phase (int k)
/* ---------------------------------------------------------------------- */
{
/*
 *      Dumps gas phase data
 */
  int i, n;
  LDBLE lp;
  struct rxn_token *rxn_ptr;
  struct gas_phase *gas_phase_ptr;

  if ((gas_phase_ptr = gas_phase_bsearch (k, &n)) == NULL)
    return (OK);

  output_msg (OUTPUT_DUMP, "GAS_PHASE  %d\n", k);
  output_msg (OUTPUT_DUMP, "\t-pressure%15.6e\n",
	      (double) gas_phase[n].total_p);
  output_msg (OUTPUT_DUMP, "\t-volume%15.6e\n", (double) gas_phase[n].volume);
  output_msg (OUTPUT_DUMP, "\t-temperature%15.6e\n",
	      (double) gas_phase[n].temperature);
  for (i = 0; i < gas_phase[n].count_comps; i++)
  {
/*
 *   Calculate partial pressure
 */
    lp = -gas_phase_ptr->comps[i].phase->lk;
    if (use.gas_phase_ptr->comps[i].phase->rxn_x != NULL)
    {
      for (rxn_ptr = use.gas_phase_ptr->comps[i].phase->rxn_x->token + 1;
	   rxn_ptr->s != NULL; rxn_ptr++)
      {
	lp += rxn_ptr->s->la * rxn_ptr->coef;
      }
    }
    else
    {
      lp = -99.9;
    }
    output_msg (OUTPUT_DUMP, "\t%-15s%15.6e\n",
		gas_phase[n].comps[i].name, (double) lp);
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
dump_s_s_assemblage (int k)
/* ---------------------------------------------------------------------- */
{
/*
 *      Dumps solid solution data
 */
  int i, j, n;
  struct s_s_assemblage *s_s_assemblage_ptr;

  s_s_assemblage_ptr = s_s_assemblage_bsearch (k, &n);
  if (s_s_assemblage_ptr == NULL)
    return (OK);

  output_msg (OUTPUT_DUMP, "SOLID_SOLUTIONS  %d\n", k);
  for (i = 0; i < s_s_assemblage[n].count_s_s; i++)
  {
    output_msg (OUTPUT_DUMP, "\t%-15s\n", s_s_assemblage[n].s_s[i].name);
    for (j = 0; j < s_s_assemblage[n].s_s[i].count_comps; j++)
    {
      output_msg (OUTPUT_DUMP, "\t\t-comp\t%-15s%15.6e\n",
		  s_s_assemblage[n].s_s[i].comps[j].name,
		  (double) s_s_assemblage[n].s_s[i].comps[j].moles);
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
dump_kinetics (int k)
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
  output_msg (OUTPUT_DUMP, "KINETICS  %d\n", k);
  for (i = 0; i < kinetics[n].count_comps; i++)
  {
    output_msg (OUTPUT_DUMP, "%-15s\n", kinetics_ptr->comps[i].rate_name);

    kinetics_comp_ptr = &kinetics_ptr->comps[i];
    output_msg (OUTPUT_DUMP, "\t-formula ");
    for (j = 0; j < kinetics_comp_ptr->count_list; j++)
    {
      output_msg (OUTPUT_DUMP, "   %s  %12.3e",
		  kinetics_comp_ptr->list[j].name,
		  (double) kinetics_comp_ptr->list[j].coef);
    }
    output_msg (OUTPUT_DUMP, "\n");

    output_msg (OUTPUT_DUMP, "\t-tol %15.2e\n",
		(double) kinetics_ptr->comps[i].tol);
    output_msg (OUTPUT_DUMP, "\t-m0  %15.6e\n",
		(double) kinetics_ptr->comps[i].m0);
    output_msg (OUTPUT_DUMP, "\t-m   %15.6e\n",
		(double) kinetics_ptr->comps[i].m);

    if (kinetics_comp_ptr->count_d_params != 0)
    {
      output_msg (OUTPUT_DUMP, "\t-parm");
      for (j = 0; j < kinetics_comp_ptr->count_d_params; j++)
      {
	output_msg (OUTPUT_DUMP, "%15.6e",
		    (double) kinetics_comp_ptr->d_params[j]);
      }
      output_msg (OUTPUT_DUMP, "\n");
    }

/* not dumped:		kinetics_comp_ptr->count_c_params = 0 */
  }
  output_msg (OUTPUT_DUMP, "\t-step_divide   %15.6e\n",
	      (double) kinetics[n].step_divide);
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
dump_solution (int k)
/* ---------------------------------------------------------------------- */
{
/*
 *      Dumps solution data
 */
  int i, j;
  struct solution *solution_ptr;

  solution_ptr = solution_bsearch (k, &j, TRUE);
  if (solution_ptr != NULL)
  {
    output_msg (OUTPUT_DUMP, "SOLUTION %d\n", k);
    output_msg (OUTPUT_DUMP, "\tunits  mol/kgw\n");
    output_msg (OUTPUT_DUMP, "\ttemp  %5.1f\n", (double) solution_ptr->tc);
    output_msg (OUTPUT_DUMP, "\tpH    %9.5f\n", (double) solution_ptr->ph);
    output_msg (OUTPUT_DUMP, "\tpe    %9.5f\n",
		(double) solution_ptr->solution_pe);
    output_msg (OUTPUT_DUMP, "\twater %9.5f\n",
		(double) solution_ptr->mass_water);
    for (i = 0; solution_ptr->totals[i].description != NULL; i++)
    {
      /*
       *                              write out totals...
       */
      output_msg (OUTPUT_DUMP, "\t%-6s", solution_ptr->totals[i].description);
      output_msg (OUTPUT_DUMP, " %13.5e\n",
		  (double) (solution_ptr->totals[i].moles /
			    solution_ptr->mass_water));
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
dump_mix (int k)
/* ---------------------------------------------------------------------- */
{
/*
 *      Dumps solution data
 */
  int j;
  struct mix *mix_ptr;

  mix_ptr = mix_bsearch (k, &j);
  if (mix_ptr == NULL)
    return (ERROR);
  for (j = 0; j < mix_ptr->count_comps; j++)
  {
    dump_solution (mix_ptr->comps[j].n_solution);
  }
  output_msg (OUTPUT_DUMP, "MIX %d\n", k);
  for (j = 0; j < mix_ptr->count_comps; j++)
  {
    output_msg (OUTPUT_DUMP, "\t5%d\t%13.5e\n", mix_ptr->comps[j].n_solution,
		(double) mix_ptr->comps[j].fraction);
  }
  return (OK);
}
