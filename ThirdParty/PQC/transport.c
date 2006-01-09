#define  EXTERNAL extern
#include "global.h"
#include "phqalloc.h"
#include "output.h"
#include "phrqproto.h"
static char const svnid[] =
  "$Id: transport.c 4 2009-04-21 17:29:29Z delucia $";

struct spec
{
  char *name;			/* name of species */
  LDBLE a;			/* activity */
  LDBLE lm;			/* log(concentration) */
  LDBLE lg;			/* log(gamma) */
  LDBLE c;			/* concentration */
  LDBLE z;			/* charge number */
  LDBLE Dp;			/* pore water diffusion coefficient, m2/s */
};
struct sol_D
{
  int count_spec;		/* number of species */
  struct spec *spec;
} *sol_D;

struct J_ij
{
  char *name;
  LDBLE tot;
} *J_ij;

static int multi_D (LDBLE stagkin_time);
static int find_J (int cell_no);
static int fill_spec (int cell_no);
static int sort_species_name (const void *ptr1, const void *ptr2);
static int multi_Dstag (int mobile_cell);
static int find_Jstag (int icell, int jcell, LDBLE mixf);

LDBLE diffc_max, J_ij_sum;
int J_ij_count_spec;
int transp_surf;
static int disp_surf (LDBLE stagkin_time);
static int diff_stag_surf (int mobile_cell);
static int check_surfaces (struct surface *surface_ptr1,
			   struct surface *surface_ptr2);
int sum_surface_comp (struct surface *source1, LDBLE f1,
		      struct surface *source2, int k, LDBLE f2,
		      struct surface *target, LDBLE new_Dw);
static int mobile_surface_copy (struct surface *surface_old_ptr,
				struct surface *surf_ptr1, int n_user_new,
				int move_old);

static int init_mix (void);
static int init_heat_mix (int nmix);
static int heat_mix (int heat_nmix);
static int mix_stag (int i, LDBLE stagkin_time, int punch,
		     LDBLE step_fraction_kin);


LDBLE *heat_mix_array;
LDBLE *temp1, *temp2;
int nmix, heat_nmix;
LDBLE heat_mix_f_imm, heat_mix_f_m;

/* ---------------------------------------------------------------------- */
int
transport (void)
/* ---------------------------------------------------------------------- */
{
  int i, j, k, n;
  int j_imm, n_m, n_imm;
  LDBLE b, f, mix_f_m, mix_f_imm;
  LDBLE water_m, water_imm;
  int first_c, last_c, b_c;
  int max_iter;
  char token[MAX_LENGTH];
  LDBLE kin_time, stagkin_time, kin_time_save;
  struct mix *mix_ptr;
  struct surface *surf_ptr, *surf_ptr1;
  int punch_boolean;
  LDBLE step_fraction;
  LDBLE cb_tol;
  if (svnid == NULL)
    fprintf (stderr, " ");

  state = TRANSPORT;
  diffc_max = 0.0;
  cb_tol = 1e-9;
  transp_surf = 0;

/*        mass_water_switch = TRUE; */
/*
 *   Check existence of solutions
 */
  j = -1;
  /* check column solutions */
  for (i = 1; i <= count_cells; i++)
  {
    use.solution_ptr = solution_bsearch (i, &n, TRUE);
    if (use.solution_ptr == NULL)
    {
      input_error++;
      sprintf (error_string,
	       "Solution %d is needed for transport, but is not defined.", i);
      error_msg (error_string, CONTINUE);
    }
    else
      cell_data[i - 1].temp = use.solution_ptr->tc;
  }

  if (multi_Dflag == TRUE)
  {
    sol_D =
      (struct sol_D *)
      PHRQ_malloc ((size_t)
		   (count_cells + 2 +
		    stag_data->count_stag * count_cells) *
		   sizeof (struct sol_D));
    if (sol_D == NULL)
      malloc_error ();
    for (i = 0; i < count_cells + 2 + stag_data->count_stag * count_cells;
	 i++)
    {
      sol_D[i].count_spec = 0;
      sol_D[i].spec = NULL;
    }
  }
  /* check solution 0 */
  if (solution_bsearch (0, &n, FALSE) == NULL)
  {
    if (ishift == 1)
    {
      input_error++;
      sprintf (error_string,
	       "Solution 0 is needed for transport, but is not defined.");
      error_msg (error_string, CONTINUE);
    }
    else
      solution_duplicate (1, 0);
  }

  /* check solution count_cells */
  if (solution_bsearch (count_cells + 1, &n, FALSE) == NULL)
  {
    if (ishift == -1)
    {
      input_error++;
      sprintf (error_string,
	       "Solution %d is needed for transport, but is not defined.",
	       count_cells + 1);
      error_msg (error_string, CONTINUE);
    }
    else
      solution_duplicate (count_cells, count_cells + 1);
  }
/*
 *   Initialize temperature in stagnant cells ...
 */
  for (n = 1; n <= stag_data->count_stag; n++)
  {
    for (i = 1; i <= count_cells; i++)
    {
      k = i + 1 + n * count_cells;
      use.solution_ptr = solution_bsearch (k, &use.n_solution, FALSE);
      if (use.solution_ptr != NULL)
	cell_data[k - 1].temp = use.solution_ptr->tc;
    }
  }
/*
 * First equilibrate solutions
 */
  dup_print ("Equilibrating initial solutions", TRUE);
  transport_step = 0;
  for (i = 0; i <= count_cells + 1; i++)
  {
    set_initial_moles (i);
    cell_no = i;
    set_and_run_wrapper (i, NOMIX, FALSE, i, 0.0);
    if (use.surface_ptr != NULL && use.surface_ptr->transport == TRUE)
      transp_surf = TRUE;
    if (multi_Dflag == TRUE)
    {
#ifdef SKIP
      /* check for charge balance of solutions... */
      if (fabs (cb_x) > (cb_tol * mu_x))
      {
	/* make it an error... */
	input_error++;
	sprintf (error_string,
		 "Solution %d must be charge-balanced for multicomponent diffusion.",
		 i);
	error_msg (error_string, CONTINUE);
	/* or just a warning ... */
	sprintf (token,
		 "Solution %d has %.2e charge imbalance in multicomponent diffusion.",
		 i, cb_x);
	warning_msg (token);
      }
#endif
      fill_spec (cell_no);
    }
    if (cell_no > 0 && cell_no <= count_cells) 
    {
      if ((cell_data[i - 1].punch == TRUE) && (cell_no != count_cells + 1))
	punch_all ();
      if ((cell_data[i - 1].print == TRUE) && (cell_no != count_cells + 1))
	print_all ();
    }
/*    if (i > 0 && i <= count_cells)*/
      saver ();
  }
/*
 * Also stagnant cells
 */
  for (n = 1; n <= stag_data->count_stag; n++)
  {
    for (i = 1; i <= count_cells; i++)
    {
      k = i + 1 + n * count_cells;
      cell_no = k;
      if (solution_bsearch (k, &use.n_solution, FALSE) != 0)
      {
	set_initial_moles (k);
	set_and_run_wrapper (k, NOMIX, FALSE, k, 0.0);
	if (multi_Dflag == TRUE)
	{
	  fill_spec (cell_no);
	}
	if ((cell_data[k - 1].punch == TRUE))
	  punch_all ();
	if ((cell_data[k - 1].print == TRUE)
	    && (transport_step % print_modulus == 0))
	  print_all ();
	saver ();
      }
    }
  }
/*
 *  Initialize mixing factors, define kinetics times
 *  for multicomponent diffusion, limit mixing by diffc_max (usually from H+)
 */
  if (multi_Dflag == TRUE)
    diffc = diffc_max;
  if ((stag_data->exch_f > 0) && (stag_data->count_stag == 1))
  {
    /* multi_D calc's are OK if all cells have the same amount of water */
/*              if (multi_Dflag == TRUE) {
                        sprintf(token, "multi_D calc's and stagnant: define MIXing factors explicitly, or \n\t give all cells the same amount of water.");
                        warning_msg(token);
                }
 */ mix_ptr = &mix[0];
    for (i = 0; i < count_mix; i++)
      mix_free (mix_ptr++);
    count_mix = 2 * count_cells;
/*
 * stagnant mix factors go in mix[0 .. count_cells]
 */
    mix =
      (struct mix *) PHRQ_realloc (mix,
				   (size_t) count_mix * sizeof (struct mix));
    if (mix == NULL)
      malloc_error ();
    memset (mix, 0, sizeof (struct mix) * count_mix);
  }
/*
 * mix[] is extended in init_mix(), to accommodate column mix factors
 */
  nmix = init_mix ();
  heat_nmix = init_heat_mix (nmix);
  if (nmix < 2)
    stagkin_time = timest;
  else
    stagkin_time = timest / nmix;
  if (ishift != 0)
    kin_time = timest / (1 + nmix);
  else
    kin_time = stagkin_time;
  kin_time_save = kin_time;

/* Reaction defined for a shift... */
  step_fraction = 1.0 / (1.0 + nmix);
/*
 *   Set boundary conditions, transport direction
 */
  last_model.force_prep = TRUE;
  if ((ishift == 0) || (bcon_first == 1) || (bcon_last == 1))
    b_c = 1;
  else
    b_c = 0;
  if (ishift >= 0)
  {
    last_c = count_cells;
    first_c = 1;
  }
  else
  {
    last_c = 1;
    first_c = count_cells;
  }
/*
 * Define stagnant/mobile mix structure, if not read explicitly.
 *
 * With count_stag = 1, mix factors are calculated from exchange factor à
 * (= exch_f), mobile é_m (= th_m) and immobile é_im (= th_im) porosity.
 * These variables are read under keyword TRANSPORT, after stagnant, in
 * structure stag_data.
 * MIX 'cell_no' in input file can be an alternative for the calculation here.
 */

  if ((stag_data->exch_f > 0) && (stag_data->count_stag == 1))
  {

    b = stag_data->th_m / (stag_data->th_m + stag_data->th_im);
    f = exp (-stag_data->exch_f * stagkin_time / (b * stag_data->th_im));
    mix_f_imm = b - b * f;
    mix_f_m = mix_f_imm * stag_data->th_im / stag_data->th_m;

    n = 0;
    for (i = 1; i <= count_cells; i++)
    {
      j = i;
      j_imm = j + (1 + count_cells);
      if (solution_bsearch (j, &n_m, TRUE) == NULL)
	error_msg ("Could not find mobile cell solution in TRANSPORT.", STOP);
      if (solution_bsearch (j_imm, &n_imm, TRUE) == NULL)
	error_msg ("Could not find immobile cell solution in TRANSPORT.",
		   STOP);
      water_m = solution[n_m]->mass_water;
      water_imm = solution[n_imm]->mass_water;
/*
 * Define C_m = (1 - mix_f_m) * C_m0  +  mix_f_m) * C_im0
 */
      mix[n].comps =
	(struct mix_comp *) PHRQ_malloc ((size_t) 2 *
					 sizeof (struct mix_comp));
      if (mix[n].comps == NULL)
	malloc_error ();
      mix[n].count_comps = 2;
      mix[n].description = string_duplicate (" ");
      mix[n].n_user = j;
      mix[n].n_user_end = j;
      mix[n].comps[0].n_solution = j;
      mix[n].comps[0].fraction = 1 - mix_f_m;
      mix[n].comps[1].n_solution = j_imm;
      mix[n].comps[1].fraction = mix_f_m * water_m / water_imm;
      n++;
/*
 * Define C_im = mix_f_imm * C_m0  +  (1 - mix_f_imm) * C_im0,  or...
 */
      mix[n].comps =
	(struct mix_comp *) PHRQ_malloc ((size_t) 2 *
					 sizeof (struct mix_comp));
      if (mix[n].comps == NULL)
	malloc_error ();
      mix[n].count_comps = 2;
      mix[n].description = string_duplicate (" ");
      mix[n].n_user = j_imm;
      mix[n].n_user_end = j_imm;
      mix[n].comps[0].n_solution = j_imm;
      mix[n].comps[0].fraction = 1 - mix_f_imm;
      mix[n].comps[1].n_solution = j;
      mix[n].comps[1].fraction = mix_f_imm * water_imm / water_m;
      n++;
    }

    if (heat_nmix > 0)
    {
/*
 * Assumption: D_e used for calculating exch_f in input file equals diffc
 */
      f = stag_data->exch_f * (heat_diffc - diffc) / diffc / tempr;
      f = exp (-f * stagkin_time / (b * stag_data->th_im));
      heat_mix_f_imm = b - b * f;
      heat_mix_f_m = heat_mix_f_imm * stag_data->th_im / stag_data->th_m;
    }
  }
/*
 *   Stop if error
 */
  if (input_error > 0)
  {
    error_msg ("Program terminating due to input errors.", STOP);
  }
/*
 * Now transport
 */
  max_iter = 0;
  for (transport_step = transport_start; transport_step <= count_shifts;
       transport_step++)
  {
    /*
     *  Set initial moles of phases
     */
    for (i = 1; i <= count_cells; i++)
      set_initial_moles (i);
    /*
     * Also stagnant cells
     */
    for (n = 1; n <= stag_data->count_stag; n++)
    {
      for (i = 1; i <= count_cells; i++)
      {
	k = i + 1 + n * count_cells;
	cell_no = k;
	if (solution_bsearch (k, &use.n_solution, FALSE) != 0)
	{
	  set_initial_moles (k);
	}
      }
    }
/*
 * Start diffusing if boundary cond = 1, (fixed c, or closed)
 */
    if (b_c == 1)
    {

      /* For half of mixing steps */
      for (j = 1; j <= floor ((double) nmix / 2); j++)
      {
	rate_sim_time_start =
	  (transport_step - 1) * timest + (j - 1) * kin_time;
	rate_sim_time = rate_sim_time_start + kin_time;

	if (multi_Dflag)
	  sprintf (token,
		   "Transport step %3d. Multicomponent diffusion run %3d.",
		   transport_step, j);
	else
	  sprintf (token, "Transport step %3d. Mixrun %3d.", transport_step,
		   j);
	dup_print (token, FALSE);

	if (heat_nmix > 0)
	{
	  heat_mix (heat_nmix);
	  /* equilibrate again ... */
	  for (i = 1; i <= count_cells; i++)
	  {
	    cell_no = i;
	    set_and_run_wrapper (i, NOMIX, FALSE, i, 0.0);
	    if (multi_Dflag)
	      fill_spec (i);
	    saver ();
	  }
	}
	/* Go through cells */
	if (transp_surf)
	{
	  if (disp_surf (stagkin_time) == ERROR)
	    error_msg ("Error in surface transport, stopping.", STOP);
	}
	if (multi_Dflag)
	  multi_D (stagkin_time);

	for (i = 1; i <= count_cells; i++)
	{

	  if (iterations > max_iter)
	    max_iter = iterations;
	  if (multi_Dflag)
	    sprintf (token,
		     "Transport step %3d. MCDrun %3d. Cell %3d. (Max. iter %3d)",
		     transport_step, j, i, max_iter);
	  else
	    sprintf (token,
		     "Transport step %3d. Mixrun %3d. Cell %3d. (Max. iter %3d)",
		     transport_step, j, i, max_iter);
	  status (0, token);

	  cell_no = i;
	  run_reactions (i, kin_time, DISP, step_fraction);
	  if (multi_Dflag)
	    fill_spec (i);

	  /* punch and output file */
	  if ((ishift == 0) && (j == nmix) && ((stag_data->count_stag == 0) ||
					       solution_bsearch (i + 1 +
								 count_cells,
								 &use.
								 n_solution,
								 FALSE) == 0))
	  {
	    if ((cell_data[i - 1].punch == TRUE)
		&& (transport_step % punch_modulus == 0))
	      punch_all ();
	    if ((cell_data[i - 1].print == TRUE)
		&& (transport_step % print_modulus == 0))
	      print_all ();
	  }
	  if (i > 1)
	    solution_duplicate (-2, i - 1);
	  saver ();

	  /* maybe sorb a surface component... */
	  if ((ishift == 0) && (j == nmix) && ((stag_data->count_stag == 0) ||
					       solution_bsearch (i + 1 +
								 count_cells,
								 &use.
								 n_solution,
								 FALSE) == 0))
	  {
	    if (change_surf_count > 0)
	    {
	      for (k = 0; k < change_surf_count; k++)
	      {
		if (change_surf[k].cell_no != i)
		  break;
		reformat_surf (change_surf[k].comp_name,
			       change_surf[k].fraction,
			       change_surf[k].new_comp_name,
			       change_surf[k].new_Dw, change_surf[k].cell_no);
		change_surf[k].cell_no = -99;
	      }
	      change_surf_count = 0;
	    }
	  }
	}
	solution_duplicate (-2, count_cells);

	/* Stagnant zone mixing after completion of each
	   diffusive/dispersive step ...  */
	rate_sim_time_start =
	  (transport_step - 1) * timest + (j - 1) * stagkin_time;
	rate_sim_time = rate_sim_time_start + stagkin_time;

	if (stag_data->count_stag > 0)
	{
	  if ((ishift == 0) && (j == nmix))
	    punch_boolean = TRUE;
	  else
	    punch_boolean = FALSE;
	  for (i = 1; i <= count_cells; i++)
	    mix_stag (i, stagkin_time, punch_boolean, step_fraction);
	}
      }
    }
/*
 * Advective transport
 */
    if (ishift != 0)
    {
      sprintf (token, "Transport step %3d.", transport_step);
      dup_print (token, FALSE);
      if (b_c == 1)
	rate_sim_time_start =
	  (transport_step - 1) * timest + (j - 1) * kin_time;
      else
	rate_sim_time_start = (transport_step - 1) * timest;
      rate_sim_time = rate_sim_time_start + kin_time;

      /* halftime kinetics for resident water in first cell ... */
      if (kinetics_bsearch (first_c, &i) != NULL && count_cells > 1)
      {
	cell_no = first_c;
	kin_time = kin_time_save / 2;
	run_reactions (first_c, kin_time, NOMIX, 0.0);
	saver ();
	kin_time = kin_time_save;
      }

      /* for each cell in column */
      for (i = last_c; i != (first_c - ishift); i -= ishift)
	solution_duplicate (i - ishift, i);

/* if boundary_solutions must be flushed by the flux from the column...
      if (ishift == 1 && bcon_last == 3)
	solution_duplicate (last_c, last_c + 1);
      else if (ishift == -1 && bcon_first == 3)
	solution_duplicate (last_c, last_c - 1);
 */
      if (transp_surf)
      {
	for (i = last_c + ishift; i != (first_c - ishift); i -= ishift)
	{
	  if ((ishift == 1 && i == last_c + 1 && bcon_last != 3) ||
	      (ishift == -1 && i == last_c - 1 && bcon_first != 3))
	    continue;
	  if ((surf_ptr =
	       surface_bsearch (i - ishift, &use.n_surface)) == NULL)
	  {
	    if ((surface_bsearch (i, &use.n_surface) != NULL) &&
		((i == 0 && bcon_first == 3)
		 || (i == count_cells + 1 && bcon_last == 3)))
	      surface_delete (i);
	    continue;
	  }
	  if (surf_ptr->transport)
	  {
	    if ((surf_ptr1 = surface_bsearch (i, &use.n_surface)) == NULL)
	    {
	      n = count_surface++;
	      space ((void **) ((void *) &surface), count_surface,
		     &max_surface, sizeof (struct surface));
	      surf_ptr1 = &surface[n];
	      surf_ptr1->count_comps = 0;
	    }
	    if (i == first_c)
	      mobile_surface_copy (surf_ptr, surf_ptr1, i, FALSE);
	    else
	      mobile_surface_copy (surf_ptr, surf_ptr1, i, TRUE);
	    if (surface[n].n_user < surface[n - 1].n_user)
	      surface_sort ();
	  }
	}
      }

/*
 * thermal diffusion when nmix = 0...
 */
      if ((nmix == 0) && (heat_nmix > 0))
      {
	heat_mix (heat_nmix);
	/* equilibrate again ... */
	for (i = 1; i <= count_cells; i++)
	{
	  cell_no = i;
	  set_and_run_wrapper (i, NOMIX, FALSE, i, 0.0);
	  if (multi_Dflag)
	    fill_spec (i);
	  saver ();
	}
      }

      for (i = 1; i <= count_cells; i++)
      {
	if (i == first_c && count_cells > 1)
	  kin_time /= 2;
	if (multi_Dflag)
	  sprintf (token,
		   "Transport step %3d. MCDrun %3d. Cell %3d. (Max. iter %3d)",
		   transport_step, 0, i, max_iter);
	else
	  sprintf (token,
		   "Transport step %3d. Mixrun %3d. Cell %3d. (Max. iter %3d)",
		   transport_step, 0, i, max_iter);
	status (0, token);
	cell_no = i;
	run_reactions (i, kin_time, NOMIX, step_fraction);
	if (multi_Dflag == TRUE)
	  fill_spec (i);
	if (iterations > max_iter)
	  max_iter = iterations;

	if ((nmix == 0) && ((stag_data->count_stag == 0) ||
			    (solution_bsearch
			     (i + 1 + count_cells, &use.n_solution,
			      FALSE) == 0)))
	{
	  if ((cell_data[i - 1].punch == TRUE)
	      && (transport_step % punch_modulus == 0))
	    punch_all ();
	  if ((cell_data[i - 1].print == TRUE)
	      && (transport_step % print_modulus == 0))
	    print_all ();
	}
	if (i == first_c && count_cells > 1)
	  kin_time = kin_time_save;
	saver ();

	/* maybe sorb a surface component... */
	if ((nmix == 0) && ((stag_data->count_stag == 0) ||
			    (solution_bsearch
			     (i + 1 + count_cells, &use.n_solution,
			      FALSE) == 0)))
	{
	  if (change_surf_count > 0)
	  {
	    for (k = 0; k < change_surf_count; k++)
	    {
	      if (change_surf[k].cell_no != i)
		break;
	      reformat_surf (change_surf[k].comp_name,
			     change_surf[k].fraction,
			     change_surf[k].new_comp_name,
			     change_surf[k].new_Dw, change_surf[k].cell_no);
	      change_surf[k].cell_no = -99;
	    }
	    change_surf_count = 0;
	  }
	}

	/* If nmix is zero, stagnant zone mixing after
	   advective step ... */
	if ((nmix == 0) && (stag_data->count_stag > 0))
	{
	  mix_stag (i, stagkin_time, TRUE, step_fraction);
	}
      }
    }
/*
 * Further dispersive and diffusive transport
 */
    if (b_c != 1)
      j = 1;
    for (j = j; j <= nmix; j++)
    {
      if (multi_Dflag)
	sprintf (token,
		 "Transport step %3d. Multicomponent diffusion run %3d.",
		 transport_step, j);
      else
	sprintf (token, "Transport step %3d. Mixrun %3d.", transport_step, j);
      dup_print (token, FALSE);
      rate_sim_time_start =
	(transport_step - 1) * timest + (j - 1) * kin_time;
      if (ishift != 0)
	rate_sim_time_start += kin_time;
      rate_sim_time = rate_sim_time_start + kin_time;

      if (heat_nmix > 0)
      {
	heat_mix (heat_nmix);
	/* equilibrate again ... */
	for (i = 1; i <= count_cells; i++)
	{
	  cell_no = i;
	  set_and_run_wrapper (i, NOMIX, FALSE, i, 0.0);
	  if (multi_Dflag)
	    fill_spec (i);
	  saver ();
	}
      }
      if (transp_surf)
      {
	if (disp_surf (stagkin_time) == ERROR)
	  error_msg ("Error in surface transport, stopping.", STOP);
      }

      if (multi_Dflag == TRUE)
	multi_D (stagkin_time);

      /* for each cell in column */
      for (i = 1; i <= count_cells; i++)
      {
	if (iterations > max_iter)
	  max_iter = iterations;
	if (multi_Dflag)
	  sprintf (token,
		   "Transport step %3d. MCDrun %3d. Cell %3d. (Max. iter %3d)",
		   transport_step, j, i, max_iter);
	else
	  sprintf (token,
		   "Transport step %3d. Mixrun %3d. Cell %3d. (Max. iter %3d)",
		   transport_step, j, i, max_iter);
	status (0, token);
	cell_no = i;
	run_reactions (i, kin_time, DISP, step_fraction);
	if (multi_Dflag == TRUE)
	  fill_spec (i);

	if ((j == nmix) && ((stag_data->count_stag == 0) ||
			    (solution_bsearch
			     (i + 1 + count_cells, &use.n_solution,
			      FALSE) == 0)))
	{
	  if ((cell_data[i - 1].punch == TRUE)
	      && (transport_step % punch_modulus == 0))
	    punch_all ();
	  if ((cell_data[i - 1].print == TRUE)
	      && (transport_step % print_modulus == 0))
	    print_all ();
	}
	if (i > 1)
	  solution_duplicate (-2, i - 1);
	saver ();

	/* maybe sorb a surface component... */
	if ((j == nmix) && ((stag_data->count_stag == 0) ||
			    (solution_bsearch
			     (i + 1 + count_cells, &use.n_solution,
			      FALSE) == 0)))
	{
	  if (change_surf_count > 0)
	  {
	    for (k = 0; k < change_surf_count; k++)
	    {
	      if (change_surf[k].cell_no != i)
		break;
	      reformat_surf (change_surf[k].comp_name,
			     change_surf[k].fraction,
			     change_surf[k].new_comp_name,
			     change_surf[k].new_Dw, change_surf[k].cell_no);
	      change_surf[k].cell_no = -99;
	    }
	    change_surf_count = 0;
	  }
	}
      }
      solution_duplicate (-2, count_cells);

      /* Stagnant zone mixing after completion of each
         diffusive/dispersive step ... */
      rate_sim_time_start =
	(transport_step - 1) * timest + (j - 1) * stagkin_time;
      rate_sim_time = rate_sim_time_start + stagkin_time;

      if (stag_data->count_stag > 0)
      {
	if (j == nmix)
	  punch_boolean = TRUE;
	else
	  punch_boolean = FALSE;
	for (i = 1; i <= count_cells; i++)
	  mix_stag (i, stagkin_time, punch_boolean, step_fraction);
      }
    }
    if (dump_modulus != 0 && (transport_step % dump_modulus) == 0)
      dump ();
  }
#ifdef DOS
  output_msg (OUTPUT_SCREEN, "\n");
#else
  output_msg (OUTPUT_SCREEN, "%s%-80s", "\n", " ");
#endif
  /* free_model_allocs(); */
/*
 * free mix structures
 */
  if ((stag_data->exch_f > 0) && (stag_data->count_stag == 1))
  {
    mix_ptr = &mix[0];
    for (i = 0; i < count_mix; i++)
      mix_free (mix_ptr++);
    count_mix = 0;
  }
  else
  {
    if (nmix > 0)
    {
      mix_ptr = &mix[count_mix - count_cells];
      for (i = count_mix - count_cells; i < count_mix; i++)
	mix_free (mix_ptr++);

      count_mix -= count_cells;
      mix =
	(struct mix *) PHRQ_realloc (mix,
				     (size_t) (count_mix +
					       1) * sizeof (struct mix));
      if (mix == NULL)
	malloc_error ();
    }
  }
  if (heat_nmix > 0)
  {
    heat_mix_array = (double *) free_check_null (heat_mix_array);
    temp1 = (double *) free_check_null (temp1);
    temp2 = (double *) free_check_null (temp2);
  }
  if (multi_Dflag == TRUE)
  {
    for (i = 0; i < count_cells + 2 + stag_data->count_stag * count_cells;
	 i++)
      sol_D[i].spec = (struct spec *) free_check_null (sol_D[i].spec);
    sol_D = (struct sol_D *) free_check_null (sol_D);
  }

  initial_total_time += rate_sim_time;
  rate_sim_time = 0;
  mass_water_switch = FALSE;
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
init_mix (void)
/* ---------------------------------------------------------------------- */
{
  LDBLE dav, lav, mixf, maxmix, corr_disp, diffc_here, mD;
  int i, n, nmix, count_comps, max_mix;
  LDBLE *m;

  m = (LDBLE *) PHRQ_malloc ((count_cells + 1) * sizeof (LDBLE));
  if (m == NULL)
    malloc_error ();
  if (multi_Dflag == TRUE)
    diffc_here = 0.0;
  else
    diffc_here = diffc;
/*
 * Define mixing factors among inner cells
 */
  corr_disp = 1.;
  if (correct_disp == TRUE && ishift != 0)
  {
    if (bcon_first == 3)
      corr_disp += 1. / count_cells;
    if (bcon_last == 3)
      corr_disp += 1. / count_cells;
  }
  maxmix = 0.0;
  for (i = 1; i < count_cells; i++)
  {
    lav = (cell_data[i - 1].length + cell_data[i].length) / 2;
    if (ishift != 0)
      dav = (cell_data[i - 1].disp + cell_data[i].disp) / 2;
    else
      dav = 0;

    mixf =
      (diffc_here * timest / lav + dav) * corr_disp / cell_data[i].length;
    if (mixf > maxmix)
      maxmix = mixf;
    m[i] = mixf;		/* m[i] has mixf with lower cell */
    if (multi_Dflag == TRUE)
    {
      mD = diffc_max * timest / (lav * lav);
      if (mD > maxmix)
	maxmix = mD;
    }
  }
/*
 * Also for boundary cells
 */
  if (bcon_first == 1)
  {
    lav = cell_data[0].length;
    if (ishift != 0)
      dav = cell_data[0].disp;
    else
      dav = 0;

    mixf = (diffc_here * timest / lav + dav) / lav;
    if (mixf > maxmix)
      maxmix = mixf;
    m[0] = 2 * mixf;
    if (multi_Dflag == TRUE)
    {
      mD = diffc_max * timest / (lav * lav);
      if (mD > maxmix)
	maxmix = mD;
    }
  }
  else
    m[0] = 0;

  if (bcon_last == 1)
  {
    lav = cell_data[count_cells - 1].length;
    if (ishift != 0)
      dav = cell_data[count_cells - 1].disp;
    else
      dav = 0;

    mixf = (diffc_here * timest / lav + dav) / lav;
    if (mixf > maxmix)
      maxmix = mixf;
    m[count_cells] = 2 * mixf;
    if (multi_Dflag == TRUE)
    {
      mD = diffc_max * timest / (lav * lav);
      if (mD > maxmix)
	maxmix = mD;
    }
  }
  else
    m[count_cells] = 0;

/*
 * Find number of mixes
 */
  if (maxmix == 0)
  {
    nmix = 0;
    if (multi_Dflag == TRUE && mcd_substeps > 1 && stag_data->count_stag > 0)
      nmix = (int) ceil (mcd_substeps);
  }
  else
  {
    if ((bcon_first == 1) || (bcon_last == 1))
      nmix = 1 + (int) floor (4.5 * maxmix);
    else
      nmix = 1 + (int) floor (3.0 * maxmix);

    if ((ishift != 0) && ((bcon_first == 1) || (bcon_last == 1)))
    {
      if (nmix < 2)
	nmix = 2;
    }
    if (multi_Dflag == TRUE && mcd_substeps > 1)
      nmix = (int) ceil (nmix * mcd_substeps);

    for (i = 0; i <= count_cells; i++)
      m[i] /= nmix;
  }
  /*
   * Fill mix structure
   */
  if (nmix != 0)
  {
    mix =
      (struct mix *) PHRQ_realloc (mix,
				   (size_t) (count_mix +
					     count_cells) *
				   sizeof (struct mix));
    if (mix == NULL)
      malloc_error ();
    count_mix += count_cells;
    for (n = count_mix - count_cells; n < count_mix; n++)
    {
      mix[n].description = NULL;
      mix[n].count_comps = 3;
      mix[n].comps =
	(struct mix_comp *) PHRQ_malloc ((size_t) 3 *
					 sizeof (struct mix_comp));
      if (mix[n].comps == NULL)
	malloc_error ();
    }

    n = count_mix - count_cells;
/*
 * max_mix brings n_user outside range of active cells
 * mix[n].n_user = mix[n].n_user_end = -999 has same effect
 * but max_mix keeps mix in sort order in case mix_bsearch
 * is used
 */
    if (n - 1 <= 0)
      max_mix = 1;
    else
      max_mix = mix[n - 1].n_user + 1;

    if (max_mix < count_cells * (stag_data->count_stag + 1) + 1)
      max_mix = count_cells * (stag_data->count_stag + 1) + 1;

    for (i = 1; i <= count_cells; i++)
    {
      dav = 0;
      count_comps = 0;
      mix[n].description = (char *) free_check_null (mix[n].description);
      mix[n].description = string_duplicate (" ");
/*
 * again, max_mix brings n_user outside range of active cells, etc...
 */
      mix[n].n_user = max_mix + i;
      mix[n].n_user_end = max_mix + i;

      mix[n].comps[count_comps].n_solution = i - 1;
      mix[n].comps[count_comps].fraction = m[i - 1];
      dav += m[i - 1];
      count_comps++;
      mix[n].comps[count_comps].n_solution = i + 1;
      mix[n].comps[count_comps].fraction = m[i];
      dav += m[i];
      count_comps++;
      mix[n].comps[count_comps].n_solution = i;
      mix[n].comps[count_comps].fraction = 1.0 - dav;

      n++;
    }
  }
  m = (LDBLE *) free_check_null (m);
  return (nmix);
}

/* ---------------------------------------------------------------------- */
int
mix_stag (int i, LDBLE kin_time, int punch, LDBLE step_fraction)
/* ---------------------------------------------------------------------- */
{
  int j, n, k, l;
  LDBLE t_imm;
  struct solution *ptr_imm, *ptr_m;
/*
 * Kinetics in transport cell is done while transporting
 */
  for (n = 1; n <= stag_data->count_stag; n++)
  {
    k = i + 1 + n * count_cells;
    if ((ptr_imm = solution_bsearch (k, &use.n_solution, FALSE)) != NULL)
    {
      if (n == 1)
      {
	if (heat_nmix > 0)
	{
	  ptr_m = solution_bsearch (i, &use.n_solution, FALSE);
	  t_imm =
	    heat_mix_f_imm * ptr_m->tc + (1 - heat_mix_f_imm) * ptr_imm->tc;
	  ptr_m->tc =
	    heat_mix_f_m * ptr_imm->tc + (1 - heat_mix_f_m) * ptr_m->tc;
	  cell_data[i - 1].temp = ptr_m->tc;
	  cell_data[k - 1].temp = ptr_imm->tc = t_imm;
	  /* equilibrate again ... */
	  cell_no = i;
	  set_and_run_wrapper (i, NOMIX, FALSE, i, 0.0);
	  if (multi_Dflag == TRUE)
	    fill_spec (cell_no);
	  saver ();
	  cell_no = k;
	  set_and_run_wrapper (k, NOMIX, FALSE, k, 0.0);
	  if (multi_Dflag == TRUE)
	    fill_spec (cell_no);
	  saver ();
	}
/*
 * Mobile cell, kinetics already done ...
 */
	cell_no = i;
	if (transp_surf)
	{
	  if (diff_stag_surf (i) == ERROR)
	    error_msg ("Error in surface transport, stopping.", STOP);
	}
	if (multi_Dflag == TRUE)
	  multi_Dstag (i);
	set_and_run_wrapper (i, STAG, FALSE, -2, 0.0);
	if (multi_Dflag == TRUE)
	  fill_spec (cell_no);
	if ((use.kinetics_ptr = kinetics_bsearch (i, &l)) != NULL)
	{
	  use.n_kinetics_user = i;
	  use.kinetics_in = TRUE;
	}

	if (punch && (cell_data[i - 1].print == TRUE) &&
	    (transport_step % print_modulus == 0))
	  print_all ();
	if (punch && (cell_data[i - 1].punch == TRUE) &&
	    (transport_step % punch_modulus == 0))
	  punch_all ();
	saver ();

	/* maybe sorb a surface component... */
	if (punch && change_surf_count)
	{
	  for (j = 0; j < change_surf_count; j++)
	  {
	    if (change_surf[j].cell_no != i)
	      break;
	    reformat_surf (change_surf[j].comp_name, change_surf[j].fraction,
			   change_surf[j].new_comp_name,
			   change_surf[j].new_Dw, change_surf[j].cell_no);
	    change_surf[j].cell_no = -99;
	  }
	  change_surf_count = 0;
	}
      }

      cell_no = k;
      run_reactions (k, kin_time, STAG, step_fraction);
      if (multi_Dflag == TRUE)
	fill_spec (cell_no);

      if ((cell_data[k - 1].print == TRUE) && (punch == TRUE) &&
	  (transport_step % print_modulus == 0))
	print_all ();
      if ((cell_data[k - 1].punch == TRUE) && (punch == TRUE) &&
	  (transport_step % punch_modulus == 0))
	punch_all ();
      saver ();

      /* maybe sorb a surface component... */
      if (punch && change_surf_count)
      {
	for (j = 0; j < change_surf_count; j++)
	{
	  if (change_surf[j].cell_no != k)
	    break;
	  reformat_surf (change_surf[j].comp_name, change_surf[j].fraction,
			 change_surf[j].new_comp_name, change_surf[j].new_Dw,
			 change_surf[j].cell_no);
	  change_surf[j].cell_no = -99;
	}
	change_surf_count = 0;
      }
    }
  }
  for (n = 1; n <= stag_data->count_stag; n++)
  {
    k = i + 1 + n * count_cells;
    if (solution_bsearch (k, &use.n_solution, FALSE) != 0)
    {
      solution_duplicate (-2 - k, k);
      if (n == 1)
	solution_duplicate (-2, i);
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
init_heat_mix (int nmix)
/* ---------------------------------------------------------------------- */
{
  LDBLE lav, mixf, maxmix, corr_disp;
  int i, j, k, n;
  int heat_nmix;
  LDBLE t0;
/*
 * Check for need to model thermal diffusion...
 */
  if (heat_diffc <= diffc)
    return (0);
  if (count_cells < 2)
    return (0);

  heat_nmix = 0;
  t0 = solution_bsearch (0, &n, FALSE)->tc;
  for (i = 0; i < count_cells; i++)
  {
    if (fabs (cell_data[i].temp - t0) > 1.0)
    {
      heat_nmix = 1;
      break;
    }
  }
  if (heat_nmix == 0)
  {
    if (fabs (solution_bsearch (count_cells + 1, &n, FALSE)->tc - t0) > 1.0)
      heat_nmix = 1;
    for (n = 1; n <= stag_data->count_stag; n++)
    {
      for (i = 1; i < count_cells; i++)
      {
	k = i + 1 + n * count_cells;
	if (solution_bsearch (k, &j, FALSE) != 0)
	{
	  if (fabs (cell_data[k - 1].temp - t0) > 1.0)
	  {
	    heat_nmix = 1;
	    break;
	  }
	}
      }
    }
  }
  if (heat_nmix == 0)
    return (0);
/*
 * Initialize arrays...
 */
  heat_mix_array = (LDBLE *) PHRQ_malloc ((count_cells + 2) * sizeof (LDBLE));
  if (heat_mix_array == NULL)
    malloc_error ();

  temp1 = (LDBLE *) PHRQ_malloc ((count_cells + 2) * sizeof (LDBLE));
  if (temp1 == NULL)
    malloc_error ();

  temp2 = (LDBLE *) PHRQ_malloc ((count_cells + 2) * sizeof (LDBLE));
  if (temp2 == NULL)
    malloc_error ();
/*
 * Define mixing factors among inner cells...
 */
  corr_disp = 1.;
  if (correct_disp == TRUE && ishift != 0)
  {
    if (bcon_first == 3)
      corr_disp += 1. / count_cells;
    if (bcon_last == 3)
      corr_disp += 1. / count_cells;
  }
  if (nmix > 0)
    corr_disp /= nmix;
  maxmix = 0.0;
  for (i = 1; i < count_cells; i++)
  {
    lav = (cell_data[i - 1].length + cell_data[i].length) / 2;
    mixf = (heat_diffc - diffc) * timest * corr_disp / tempr / (lav * lav);
    if (mixf > maxmix)
      maxmix = mixf;
    heat_mix_array[i + 1] = mixf;	/* m[i] has mixf with lower cell */
  }
/*
 * Also for boundary cells
 */
  if (bcon_first == 1)
  {
    lav = cell_data[0].length;
    mixf = (heat_diffc - diffc) * timest * corr_disp / tempr / (lav * lav);
    if (2 * mixf > maxmix)
      maxmix = 2 * mixf;
    heat_mix_array[1] = 2 * mixf;
  }
  else
    heat_mix_array[1] = 0;

  if (bcon_last == 1)
  {
    lav = cell_data[count_cells - 1].length;
    mixf = (heat_diffc - diffc) * timest * corr_disp / tempr / (lav * lav);
    if (2 * mixf > maxmix)
      maxmix = 2 * mixf;
    heat_mix_array[count_cells + 1] = 2 * mixf;
  }
  else
    heat_mix_array[count_cells + 1] = 0;
/*
 * Find number of mixes
 */
  if (maxmix == 0)
    heat_nmix = 0;
  else
  {
    heat_nmix = 1 + (int) floor (3.0 * maxmix);
    for (i = 1; i <= count_cells + 1; i++)
      heat_mix_array[i] /= heat_nmix;
  }

  return (heat_nmix);
}

/* ---------------------------------------------------------------------- */
int
heat_mix (int heat_nmix)
/* ---------------------------------------------------------------------- */
{
  int i, j;

  for (i = 1; i <= count_cells; i++)
    temp1[i] = solution_bsearch (i, &j, FALSE)->tc;
  temp1[0] = solution_bsearch (0, &j, FALSE)->tc;
  temp1[count_cells + 1] =
    solution_bsearch ((count_cells + 1), &j, FALSE)->tc;

  for (i = 1; i <= heat_nmix; i++)
  {
    for (j = 1; j <= count_cells; j++)
      temp2[j] =
	heat_mix_array[j] * temp1[j - 1] + heat_mix_array[j + 1] * temp1[j +
									 1] +
	(1 - heat_mix_array[j] - heat_mix_array[j + 1]) * temp1[j];
    for (j = 1; j <= count_cells; j++)
      temp1[j] = temp2[j];
  }

  for (i = 1; i <= count_cells; i++)
  {
    cell_data[i - 1].temp = temp1[i];
    solution_bsearch (i, &j, FALSE)->tc = temp1[i];
  }

  return (OK);
}

/* ---------------------------------------------------------------------- */
int
set_initial_moles (int i)
/* ---------------------------------------------------------------------- */
{
  struct pp_assemblage *pp_assemblage_ptr;
  struct gas_phase *gas_phase_ptr;
  struct kinetics *kinetics_ptr;
  struct s_s_assemblage *s_s_assemblage_ptr;
  int j, k, n;
  /*
   *   Pure phase assemblage
   */
  pp_assemblage_ptr = pp_assemblage_bsearch (i, &n);
  if (pp_assemblage_ptr != NULL)
  {
    for (j = 0; j < pp_assemblage_ptr->count_comps; j++)
    {
      pp_assemblage_ptr->pure_phases[j].initial_moles =
	pp_assemblage_ptr->pure_phases[j].moles;
      if (pp_assemblage_ptr->pure_phases[j].initial_moles < 0)
	pp_assemblage_ptr->pure_phases[j].initial_moles = 0;
    }
  }
  /*
   *   Gas phase
   */
  gas_phase_ptr = gas_phase_bsearch (i, &n);
  if (gas_phase_ptr != NULL)
  {
    for (j = 0; j < gas_phase_ptr->count_comps; j++)
      gas_phase_ptr->comps[j].initial_moles = gas_phase_ptr->comps[j].moles;
  }
  /*
   *   Kinetics
   */
  kinetics_ptr = kinetics_bsearch (i, &n);
  if (kinetics_ptr != NULL)
  {
    for (j = 0; j < kinetics_ptr->count_comps; j++)
      kinetics_ptr->comps[j].initial_moles = kinetics_ptr->comps[j].m;
  }
  /*
   *   Solid solutions
   */
  s_s_assemblage_ptr = s_s_assemblage_bsearch (i, &n);
  if (s_s_assemblage_ptr != NULL)
  {
    for (k = 0; k < s_s_assemblage_ptr->count_s_s; k++)
    {
      for (j = 0; j < s_s_assemblage_ptr->s_s[k].count_comps; j++)
	s_s_assemblage_ptr->s_s[k].comps[j].init_moles =
	  s_s_assemblage_ptr->s_s[k].comps[j].moles;
    }
  }
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
multi_D (LDBLE DDt)
/* ---------------------------------------------------------------------- */
{
/* basic scheme:
 * 1. determine mole transfer (mol/s) of solute species with concentrations above
      the limit set in routine fill_spec.
 * 2. sum up as mole transfer of master_species
 * 3. add moles of master_species to solutions for mixing timestep DDt
 */
  int i, j, k, l, length, length2;
  int first_c, last_c;
  char *ptr;
  char token[MAX_LENGTH];
  struct m_s
  {
    char *name;
    LDBLE tot;
  } *m_s;
  int count_m_s;
  LDBLE tot_h, tot_o, temp;

  m_s =
    (struct m_s *) PHRQ_malloc ((size_t) count_elements *
				sizeof (struct m_s));
  if (m_s == NULL)
    malloc_error ();

  if (bcon_first == 1)
    first_c = 0;
  else
    first_c = 1;
  if (bcon_last == 1)
    last_c = count_cells;
  else
    last_c = count_cells - 1;
  for (i = first_c; i <= last_c; i++)
  {
/*
 * 1. obtain J_ij in mol/s...
 */
    find_J (i);
/*
 * 2. sum up the primary or secondary master_species from all the solute species
 *      H and O go in total_h and total_o
 */
    tot_h = tot_o = 0.0;
    count_m_s = 0;
    for (j = 0; j < J_ij_count_spec; j++)
    {
      ptr = J_ij[j].name;
      count_elts = 0;
      get_elts_in_species (&ptr, 1);
      for (k = 0; k < count_elts; k++)
      {
	if (strcmp (elt_list[k].elt->name, "H") == 0)
	  tot_h += elt_list[k].coef * J_ij[j].tot;
	else if (strcmp (elt_list[k].elt->name, "O") == 0)
	  tot_o += elt_list[k].coef * J_ij[j].tot;
	else
	{
	  for (l = 0; l < count_m_s; l++)
	  {
	    length = (int) strlen (elt_list[k].elt->name);
	    if (strncmp (m_s[l].name, elt_list[k].elt->name, length) == 0)
	    {
	      m_s[l].tot += elt_list[k].coef * J_ij[j].tot;
	      break;
	    }
	  }
	  if (l == count_m_s)
	  {
	    m_s[l].name = string_hsave (elt_list[k].elt->name);
	    m_s[l].tot = elt_list[k].coef * J_ij[j].tot;
	    count_m_s++;
	  }
	}
      }
    }
/*
 * multiply with timestep...
*/
    for (l = 0; l < count_m_s; l++)
      m_s[l].tot *= DDt;
    tot_h *= DDt;
    tot_o *= DDt;
    J_ij_sum *= DDt;
/*
 * 3. find the solutions and add or subtract mol/L...
 */
    if (i > 0)
    {
      use.solution_ptr = solution_bsearch (i, &use.n_solution, FALSE);
      use.solution_ptr->total_h -= tot_h;
      use.solution_ptr->total_o -= tot_o;
      use.solution_ptr->cb -= J_ij_sum;
      for (l = 0; l < count_m_s; l++)
      {
	temp = 0.0;
	length = (int) strlen (m_s[l].name);
	for (j = 0; use.solution_ptr->totals[j].description != NULL; j++)
	{
          length2 = (int) (size_t) strcspn(use.solution_ptr->totals[j].description, "(");
          if (strncmp
	      (m_s[l].name, use.solution_ptr->totals[j].description,
	       length) == 0 && length == length2)
	  {
	    if (use.solution_ptr->totals[j].moles < m_s[l].tot)
	    {
	      temp = use.solution_ptr->totals[j].moles;
	      use.solution_ptr->totals[j].moles = 0;
	      /* see if other redox states have more moles... */
	      for (k = 1; use.solution_ptr->totals[j + k].description != NULL;
		   k++)
	      {
                length2 = (int) (size_t) strcspn(use.solution_ptr->totals[j + k].description, "(");
		if (strncmp
		    (m_s[l].name, use.solution_ptr->totals[j + k].description,
		     length) == 0 && length == length2)
		{
		  temp += use.solution_ptr->totals[j + k].moles;
		  if (temp < m_s[l].tot)
		  {
		    use.solution_ptr->totals[j + k].moles = 0;
		  }
		  else
		  {
		    use.solution_ptr->totals[j + k].moles = temp - m_s[l].tot;
		    temp = 0.0;
		    break;
		  }
		}
	      }
	      if (temp != 0.0 && m_s[l].tot - temp > 1e-12)
	      {
		sprintf (token,
			 "Negative concentration in MCD: added %.1e moles %s in cell %d.",
			 m_s[l].tot - temp, m_s[l].name, i);
		warning_msg (token);
	      }
	    }
	    else
	      use.solution_ptr->totals[j].moles -= m_s[l].tot;
	    break;
	  }
	}
	if (use.solution_ptr->totals[j].description == NULL)
	{
	  use.solution_ptr->totals =
	    (struct conc *) PHRQ_realloc (use.solution_ptr->totals,
					  (size_t) (j +
						    2) *
					  sizeof (struct conc));
	  use.solution_ptr->totals[j].description =
	    string_hsave (m_s[l].name);
	  use.solution_ptr->totals[j].moles = -m_s[l].tot;
	  if (use.solution_ptr->totals[j].moles < 0)
	  {
	    if (use.solution_ptr->totals[j].moles < -1e-12)
	    {
	      sprintf (token,
		     "Negative concentration in MCD: added %.2e moles %s in cell %d.",
		     -use.solution_ptr->totals[j].moles, m_s[l].name, i);
	      warning_msg (token);
	    }
	    use.solution_ptr->totals[j].moles = 0;
	  }
	  use.solution_ptr->totals[j + 1].description = NULL;
	}
      }
    }
    if (i < count_cells)
    {
      use.solution_ptr = solution_bsearch (i + 1, &use.n_solution, FALSE);
      use.solution_ptr->total_h += tot_h;
      use.solution_ptr->total_o += tot_o;
      use.solution_ptr->cb += J_ij_sum;
      for (l = 0; l < count_m_s; l++)
      {
	temp = 0.0;
	length = (int) strlen (m_s[l].name);
	for (j = 0; use.solution_ptr->totals[j].description != NULL; j++)
	{
          length2 = (int) (size_t) strcspn(use.solution_ptr->totals[j].description, "(");
	  if (strncmp
	      (m_s[l].name, use.solution_ptr->totals[j].description,
	       length) == 0 && length == length2)
	  {
	    if (use.solution_ptr->totals[j].moles < -m_s[l].tot)
	    {
	      temp = use.solution_ptr->totals[j].moles;
	      use.solution_ptr->totals[j].moles = 0;
	      /* see if other redox states have more moles... */
	      for (k = 1; use.solution_ptr->totals[j + k].description != NULL;
		   k++)
	      {
                length2 = (int) (size_t) strcspn(use.solution_ptr->totals[j + k].description, "(");
		if (strncmp
		    (m_s[l].name, use.solution_ptr->totals[j + k].description,
		     length) == 0 && length == length2)
		{
		  temp += use.solution_ptr->totals[j + k].moles;
		  if (temp < -m_s[l].tot)
		  {
		    use.solution_ptr->totals[j + k].moles = 0;
		  }
		  else
		  {
		    use.solution_ptr->totals[j + k].moles = temp + m_s[l].tot;
		    temp = 0.0;
		    break;
		  }
		}
	      }
	      if (temp != 0.0 && -m_s[l].tot - temp > 1e-12)
	      {
		sprintf (token,
			 "Negative concentration in MCD: added %.3e moles %s in cell %d",
			 -m_s[l].tot - temp, m_s[l].name, i + 1);
		warning_msg (token);
	      }
	    }
	    else
	      use.solution_ptr->totals[j].moles += m_s[l].tot;
	    break;
	  }
	}
	if (use.solution_ptr->totals[j].description == NULL)
	{
	  use.solution_ptr->totals =
	    (struct conc *) PHRQ_realloc (use.solution_ptr->totals,
					  (size_t) (j +
						    2) *
					  sizeof (struct conc));
	  use.solution_ptr->totals[j].description =
	    string_hsave (m_s[l].name);
	  use.solution_ptr->totals[j].moles = m_s[l].tot;
	  if (use.solution_ptr->totals[j].moles < 0)
	  {
	    if (use.solution_ptr->totals[j].moles < -1e-12)
            {
	      sprintf (token,
		     "Negative concentration in MCD: added %.4e moles %s in cell %d.",
		     -use.solution_ptr->totals[j].moles, m_s[l].name, i + 1);
	      warning_msg (token);
            }
	    use.solution_ptr->totals[j].moles = 0;
	  }
	  use.solution_ptr->totals[j + 1].description = NULL;
	}
      }
    }
  }
  m_s = (struct m_s *) free_check_null (m_s);
  J_ij = (struct J_ij *) free_check_null (J_ij);
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
find_J (int cell_no)
/* ---------------------------------------------------------------------- */
{
/* mole transfer (mol/s) of the individual master_species:
 * Eqn 1:
 * J_ij = A_ij * (-D_i*grad(a) + D_i*z_i*c_i * SUM(D_i*z_i*grad(a)) / SUM(D_i*(z_i)^2*c_i))
 */
  int i, i_max, j, j_max, k, l, n;
  LDBLE lav, ddlm, aq1, aq2;
  LDBLE *grad, *D, *z, *Dz, *Dzc, *Dzc_dl, Dz2c, Dz2c_dl, c, A_ij;
  int *o_c, only_counter;
  struct surface *s_ptr1, *s_ptr2;
  struct surface_charge *s_charge_ptr, *s_charge_ptr1, *s_charge_ptr2;
  LDBLE dl_s, dl_aq1, dl_aq2, visc1, visc2, c_dl, temp;
  char token[MAX_LENGTH], token1[MAX_LENGTH];

  if (cell_no == 0)
  {
    lav = cell_data[0].length / 2;
    aq1 = solution_bsearch (0, &n, TRUE)->mass_water;
    aq2 = solution_bsearch (1, &n, TRUE)->mass_water;
    A_ij = aq2 / cell_data[0].length;
  }
  else if (cell_no == count_cells)
  {
    lav = cell_data[count_cells - 1].length / 2;
    aq1 = solution_bsearch (count_cells, &n, TRUE)->mass_water;
    aq2 = solution_bsearch (count_cells + 1, &n, TRUE)->mass_water;
    A_ij = aq1 / cell_data[count_cells - 1].length;
  }
  else
  {
    lav = (cell_data[cell_no - 1].length + cell_data[cell_no].length) / 2;
    aq1 = solution_bsearch (cell_no, &n, TRUE)->mass_water;
    aq2 = solution_bsearch (cell_no + 1, &n, TRUE)->mass_water;
    A_ij = aq1 / (cell_data[cell_no - 1].length * cell_data[cell_no - 1].por);
    A_ij += aq2 / (cell_data[cell_no].length * cell_data[cell_no].por);
    A_ij /= 2;
    A_ij *=
      (cell_data[cell_no - 1].por <
       cell_data[cell_no].por ? cell_data[cell_no -
					  1].por : cell_data[cell_no].por);
  }
/*
 * check if DL calculations must be made...
 */
  s_charge_ptr1 = s_charge_ptr2 = NULL;
  dl_s = dl_aq1 = dl_aq2 = temp = 0.0;
  visc1 = visc2 = 1.0;
  only_counter = FALSE;

  s_ptr1 = surface_bsearch (cell_no, &i);
  if (s_ptr1 != NULL)
  {
    if (s_ptr1->dl_type != NO_DL)
    {
      if (s_ptr1->only_counter_ions)
	only_counter = TRUE;
      /* find the one (and only one...) immobile surface comp with DL... */
      for (i = 0; i < s_ptr1->count_comps; i++)
      {
	if (s_ptr1->comps[i].Dw == 0)
	{
	  s_charge_ptr1 = &s_ptr1->charge[s_ptr1->comps[i].charge];
	  dl_aq1 = s_charge_ptr1->mass_water;
	  visc1 = s_ptr1->DDL_viscosity;
	  /* check for more comps with Dw = 0 */
	  for (j = i + 1; j < s_ptr1->count_comps; j++)
	  {
	    if (s_ptr1->comps[j].Dw == 0
		&& s_ptr1->comps[j].charge != s_ptr1->comps[i].charge)
	    {
	      k = (int) strcspn (s_ptr1->comps[i].formula, "_");
	      strncpy (token1, s_ptr1->comps[i].formula, k);
	      token1[k] = '\0';
	      sprintf (token,
		       "MCD found more than 1 fixed surface with a DDL, used %s.",
		       token1);
	      warning_msg (token);
	      break;
	    }
	  }
	  break;
	}
      }
    }
  }
  s_ptr2 = surface_bsearch (cell_no + 1, &i);
  if (s_ptr2 != NULL)
  {
    if (s_ptr2->dl_type != NO_DL)
    {
      if (s_ptr2->only_counter_ions)
	only_counter = TRUE;
      for (i = 0; i < s_ptr2->count_comps; i++)
      {
	if (s_ptr2->comps[i].Dw == 0)
	{
	  s_charge_ptr2 = &s_ptr2->charge[s_ptr2->comps[i].charge];
	  dl_aq2 = s_charge_ptr2->mass_water;
	  visc2 = s_ptr2->DDL_viscosity;
	  /* check for more comps with Dw = 0 */
	  for (j = i + 1; j < s_ptr2->count_comps; j++)
	  {
	    if (s_ptr2->comps[j].Dw == 0
		&& s_ptr2->comps[j].charge != s_ptr2->comps[i].charge)
	    {
	      k = (int) strcspn (s_ptr2->comps[i].formula, "_");
	      strncpy (token1, s_ptr2->comps[i].formula, k);
	      token1[k] = '\0';
	      sprintf (token,
		       "MCD found more than 1 fixed surface with a DDL, used %s.",
		       token1);
	      warning_msg (token);
	      break;
	    }
	  }
	  break;
	}
      }
    }
  }
  if (cell_no == 0)
    visc1 = visc2;
  else if (cell_no == count_cells)
    visc2 = visc1;

  /* in each cell: DL surface = mass_water_DL / (cell_length * tortuosity)
     free pore surface = mass_water_free / (cell_length * tortuosity)
     determine DL surface as a fraction of the total pore surface... */
  if (dl_aq1 > 0)
    dl_s = dl_aq1 / (dl_aq1 + aq1);
  if (dl_aq2 > 0)
    temp = dl_aq2 / (dl_aq2 + aq2);
  if (dl_aq1 > 0 && dl_aq2 > 0)
    /* average the 2... */
    dl_s = 0.5 * (dl_s + temp);
  else if (dl_aq2 > 0)
    /* there is one DL surface... */
    dl_s = temp;
  /* total pore surface becomes... */
  if (dl_s > 0)
  {
    if (cell_no == 0)
    {
      if (dl_aq1 > 0 && dl_aq2 > 0)
	temp = 0.5 * (dl_aq1 + dl_aq2);
      else if (dl_aq2 > 0)
	temp = dl_aq2;
      else
	temp = dl_aq1;
      A_ij += temp / cell_data[0].length;
    }
    else if (cell_no == count_cells)
    {
      if (dl_aq1 > 0 && dl_aq2 > 0)
	temp = 0.5 * (dl_aq1 + dl_aq2);
      else if (dl_aq2 > 0)
	temp = dl_aq2;
      else
	temp = dl_aq1;
      A_ij += temp / cell_data[count_cells - 1].length;
    }
    else
    {
      temp =
	dl_aq1 / (cell_data[cell_no - 1].length * cell_data[cell_no - 1].por);
      temp += dl_aq2 / (cell_data[cell_no].length * cell_data[cell_no].por);
      temp *=
	0.5 * (cell_data[cell_no - 1].por <
	       cell_data[cell_no].por ? cell_data[cell_no -
						  1].por : cell_data[cell_no].
	       por);
      A_ij += temp;
    }
  }

/*
 * malloc sufficient space...
 */
  k = sol_D[cell_no].count_spec + sol_D[cell_no + 1].count_spec;

  J_ij = (struct J_ij *) free_check_null (J_ij);
  J_ij = (struct J_ij *) PHRQ_malloc ((size_t) k * sizeof (struct J_ij));
  if (J_ij == NULL)
    malloc_error ();

  grad = (LDBLE *) PHRQ_malloc ((size_t) k * sizeof (LDBLE));
  if (grad == NULL)
    malloc_error ();

  D = (LDBLE *) PHRQ_malloc ((size_t) k * sizeof (LDBLE));
  if (D == NULL)
    malloc_error ();

  z = (LDBLE *) PHRQ_malloc ((size_t) k * sizeof (LDBLE));
  if (z == NULL)
    malloc_error ();

  Dz = (LDBLE *) PHRQ_malloc ((size_t) k * sizeof (LDBLE));
  if (Dz == NULL)
    malloc_error ();

  Dzc = (LDBLE *) PHRQ_malloc ((size_t) k * sizeof (LDBLE));
  if (Dzc == NULL)
    malloc_error ();

  Dzc_dl = (LDBLE *) PHRQ_malloc ((size_t) k * sizeof (LDBLE));
  if (Dzc_dl == NULL)
    malloc_error ();

  o_c = (int *) PHRQ_malloc ((size_t) k * sizeof (int));
  if (o_c == NULL)
    malloc_error ();

  for (i = 0; i < k; i++)
  {
    J_ij[i].name = NULL;
    J_ij[i].tot = 0.0;
    grad[i] = 0.0;
    D[i] = 0.0;
    z[i] = 0.0;
    Dz[i] = 0.0;
    Dzc[i] = 0.0;
    Dzc_dl[i] = 0.0;
    o_c[i] = 0;
  }

  Dz2c = Dz2c_dl = 0.0;

/*
 * coefficients in Eqn (1)...
 */
  i = j = k = 0;
  i_max = sol_D[cell_no].count_spec;
  j_max = sol_D[cell_no + 1].count_spec;

  while (i < i_max || j < j_max)
  {
    if (j == j_max
	|| (i < i_max
	    && strcmp (sol_D[cell_no].spec[i].name,
		       sol_D[cell_no + 1].spec[j].name) < 0))
    {
      /* species 'name' is only in cell_no */
      J_ij[k].name = string_hsave (sol_D[cell_no].spec[i].name);
      D[k] = sol_D[cell_no].spec[i].Dp;
      z[k] = sol_D[cell_no].spec[i].z;
      Dz[k] = D[k] * z[k];
      Dzc[k] = Dz[k] * sol_D[cell_no].spec[i].c / 2;
      if (dl_s > 0)
      {
	o_c[k] = 1;
	s_charge_ptr = (dl_aq1 > 0) ? s_charge_ptr1 : s_charge_ptr2;
	for (l = 0; l < s_charge_ptr->count_g; l++)
	{
	  if (equal (s_charge_ptr->g[l].charge, z[k], G_TOL) == TRUE)
	  {
	    if (only_counter &&
		((s_charge_ptr->la_psi < 0 && z[k] < 0)
		 || (s_charge_ptr->la_psi > 0 && z[k] > 0)))
	    {
	      o_c[k] = 0;
	      Dzc_dl[k] = 0;
	    }
	    else
	    {
	      if (dl_aq1 > 0)
		Dzc_dl[k] =
		  Dz[k] * sol_D[cell_no].spec[i].c / 2 * (1 +
							  s_charge_ptr->g[l].
							  g * aq1 / dl_aq1);
	      else
		Dzc_dl[k] = Dz[k] * sol_D[cell_no].spec[i].c / 2;
	    }
	    break;
	  }
	}
	Dz2c_dl += Dzc_dl[k] * z[k];
      }
      Dz2c += Dzc[k] * z[k];
      grad[k] = -sol_D[cell_no].spec[i].a / lav;
      if (i < i_max)
	i++;
    }
    else if (i == i_max
	     ||
	     (strcmp
	      (sol_D[cell_no].spec[i].name,
	       sol_D[cell_no + 1].spec[j].name) > 0 && j < j_max))
    {
      /* species 'name' is only in (cell_no + 1) */
      J_ij[k].name = string_hsave (sol_D[cell_no + 1].spec[j].name);
      D[k] = sol_D[cell_no + 1].spec[j].Dp;
      z[k] = sol_D[cell_no + 1].spec[j].z;
      Dz[k] = D[k] * z[k];
      Dzc[k] = Dz[k] * sol_D[cell_no + 1].spec[j].c / 2;
      if (dl_s > 0)
      {
	o_c[k] = 1;
	s_charge_ptr = (dl_aq2 > 0) ? s_charge_ptr2 : s_charge_ptr1;
	for (l = 0; l < s_charge_ptr->count_g; l++)
	{
	  if (equal (s_charge_ptr->g[l].charge, z[k], G_TOL) == TRUE)
	  {
	    if (only_counter &&
		((s_charge_ptr->la_psi < 0 && z[k] < 0)
		 || (s_charge_ptr->la_psi > 0 && z[k] > 0)))
	    {
	      o_c[k] = 0;
	      Dzc_dl[k] = 0;
	    }
	    else
	    {
	      if (dl_aq2 > 0)
		Dzc_dl[k] =
		  Dz[k] * sol_D[cell_no + 1].spec[j].c / 2 * (1 +
							      s_charge_ptr->
							      g[l].g * aq2 /
							      dl_aq2);
	      else
		Dzc_dl[k] = Dz[k] * sol_D[cell_no + 1].spec[j].c / 2;
	    }
	    break;
	  }
	}
	Dz2c_dl += Dzc_dl[k] * z[k];
      }
      Dz2c += Dzc[k] * z[k];
      grad[k] = sol_D[cell_no + 1].spec[j].a / lav;
      if (j < j_max)
	j++;
    }
    else
      if (strcmp
	  (sol_D[cell_no].spec[i].name, sol_D[cell_no + 1].spec[j].name) == 0)
    {
      /* species 'name' is in both cells */
      J_ij[k].name = string_hsave (sol_D[cell_no].spec[i].name);
      if (sol_D[cell_no].spec[i].Dp == 0
	  || sol_D[cell_no + 1].spec[j].Dp == 0)
	D[k] = 0.0;
      else
	D[k] =
	  (sol_D[cell_no].spec[i].Dp + sol_D[cell_no + 1].spec[j].Dp) / 2;
      z[k] = sol_D[cell_no].spec[i].z;
      Dz[k] = D[k] * z[k];
      Dzc[k] =
	Dz[k] * (sol_D[cell_no].spec[i].c + sol_D[cell_no + 1].spec[j].c) / 2;
/*                      Dzc[k] = Dz[k] * (sol_D[cell_no].spec[i].c > sol_D[cell_no + 1].spec[j].c ? sol_D[cell_no].spec[i].c : sol_D[cell_no + 1].spec[j].c);
 */
      if (dl_s > 0)
      {
	c_dl = 0.0;
	o_c[k] = 1;
	if (dl_aq1 > 0)
	{
	  for (l = 0; l < s_charge_ptr1->count_g; l++)
	  {
	    if (equal (s_charge_ptr1->g[l].charge, z[k], G_TOL) == TRUE)
	    {
	      if (only_counter &&
		  ((s_charge_ptr1->la_psi < 0 && z[k] < 0)
		   || (s_charge_ptr1->la_psi > 0 && z[k] > 0)))
	      {
		o_c[k] = 0;
		c_dl = 0;
	      }
	      else
		c_dl =
		  sol_D[cell_no].spec[i].c / 2 * (1 +
						  s_charge_ptr1->g[l].g *
						  aq1 / dl_aq1);
	      break;
	    }
	  }
	}
	else
	  c_dl = sol_D[cell_no].spec[i].c / 2;

	if (dl_aq2 > 0)
	{
	  for (l = 0; l < s_charge_ptr2->count_g; l++)
	  {
	    if (equal (s_charge_ptr2->g[l].charge, z[k], G_TOL) == TRUE)
	    {
	      if (only_counter &&
		  ((s_charge_ptr2->la_psi < 0 && z[k] < 0)
		   || (s_charge_ptr2->la_psi > 0 && z[k] > 0)))
	      {
		o_c[k] = 0;
		c_dl = 0;
	      }
	      if (o_c[k] == 1)
		c_dl +=
		  sol_D[cell_no + 1].spec[j].c / 2 * (1 +
						      s_charge_ptr2->g[l].g *
						      aq2 / dl_aq2);
	      break;
	    }
	  }
	}
	else if (o_c[k] == 1)
	  c_dl += sol_D[cell_no + 1].spec[j].c / 2;

	Dzc_dl[k] = Dz[k] * c_dl;
	Dz2c_dl += Dzc_dl[k] * z[k];
      }
      Dz2c += Dzc[k] * z[k];
      grad[k] =
	(sol_D[cell_no + 1].spec[j].c - sol_D[cell_no].spec[i].c) / lav;
      ddlm = sol_D[cell_no + 1].spec[j].lm - sol_D[cell_no].spec[i].lm;
      if (fabs (ddlm) > 1e-10)
	grad[k] *=
	  (1 +
	   (sol_D[cell_no + 1].spec[j].lg -
	    sol_D[cell_no].spec[i].lg) / ddlm);
      if (i < i_max)
	i++;
      if (j < j_max)
	j++;
    }
    k++;
  }
/*
 * fill in J_ij...
 */
  if (Dz2c == 0)
    k = 0;
  J_ij_count_spec = k;
  J_ij_sum = 0;
  c = c_dl = 0.0;
  for (j = 0; j < k; j++)
  {
    c += Dz[j] * grad[j];
    c_dl += o_c[j] * Dz[j] * grad[j];
  }
  for (i = 0; i < k; i++)
  {
    J_ij[i].tot = -D[i] * grad[i] + c * Dzc[i] / Dz2c;
    J_ij[i].tot *= (1 - dl_s);
    if (Dz2c_dl > 0)
    {
      temp =
	(-D[i] * grad[i] + c_dl * Dzc_dl[i] / Dz2c_dl) * (2 / (visc1 + visc2));
      if ((J_ij[i].tot <= 0 && temp <= 0) || (J_ij[i].tot > 0 && temp > 0))
      {
	J_ij[i].tot += o_c[i] * temp * dl_s;
      }
    }
    J_ij[i].tot *= A_ij;
    J_ij_sum += z[i] * J_ij[i].tot;
  }
/* appt 3 July 07, improved convergence without transporting charge imbalance */
  J_ij_sum = 0;
  D = (LDBLE *) free_check_null (D);
  z = (LDBLE *) free_check_null (z);
  Dz = (LDBLE *) free_check_null (Dz);
  Dzc = (LDBLE *) free_check_null (Dzc);
  Dzc_dl = (LDBLE *) free_check_null (Dzc_dl);
  grad = (LDBLE *) free_check_null (grad);
  o_c = (int *) free_check_null (o_c);

  return (OK);
}

/* ---------------------------------------------------------------------- */
int
fill_spec (int cell_no)
/* ---------------------------------------------------------------------- */
{
/* copy species activities into sol_D.spec... */

  int i, count_spec;
  LDBLE lm;
  LDBLE por, por_factor, viscos;

  sol_D[cell_no].spec = (struct spec *) free_check_null (sol_D[cell_no].spec);
  sol_D[cell_no].spec =
    (struct spec *) PHRQ_malloc ((size_t) count_species_list *
				 sizeof (struct spec));
  if (sol_D[cell_no].spec == NULL)
    malloc_error ();

  if (cell_no == 0)
    por = cell_data[0].por;
  else if (cell_no == count_cells + 1)
    por = cell_data[count_cells - 1].por;
  else
    por = cell_data[cell_no - 1].por;

  if (por < multi_Dpor_lim)
    por = por_factor = 0.0;
  else
    por_factor = pow (por, multi_Dn);
/*
 * correct diffusion coefficient for temperature and viscosity, D_T = D_298 * Tk * viscos_298 / (298 * viscos)
 */
  viscos =
    pow (10.,
	 -(1.37023 * (tc_x - 20) +
	   0.000836 * (tc_x - 20) * (tc_x - 20)) / (109 + tc_x));
/*
 * put temperature factor in por_factor which corrects for porous medium...
 */
  por_factor *= tk_x * 0.88862 / (298.15 * viscos);

  count_spec = 0;
/*
 * sort species by name...
 */
  if (count_species_list > 0)
    qsort (&species_list[0], (size_t) count_species_list,
	   (size_t) sizeof (struct species_list), sort_species_name);

  for (i = 0; i < count_species_list; i++)
  {
    if (species_list[i].s->type == EX)
      continue;
    if (species_list[i].s->type == SURF)
      continue;
/*
 *   copy species data
 */
    if (i > 0
	&& strcmp (species_list[i].s->name, species_list[i - 1].s->name) == 0)
      continue;
    if (species_list[i].s == s_h2o)
      continue;
/* appt                 lm = s_h2o->la - species_list[i].s->lg;
                else
 */ lm = species_list[i].s->lm;

    if (lm > -20)
    {
      sol_D[cell_no].spec[count_spec].name =
	string_hsave (species_list[i].s->name);
      sol_D[cell_no].spec[count_spec].c =
	species_list[i].s->moles / mass_water_aq_x;
      sol_D[cell_no].spec[count_spec].a = under (lm + species_list[i].s->lg);
      sol_D[cell_no].spec[count_spec].lm = lm;
      sol_D[cell_no].spec[count_spec].lg = species_list[i].s->lg;
      sol_D[cell_no].spec[count_spec].z = species_list[i].s->z;
      if (species_list[i].s->dw == 0)
	sol_D[cell_no].spec[count_spec].Dp = default_Dw * por_factor;
      else
	sol_D[cell_no].spec[count_spec].Dp =
	  species_list[i].s->dw * por_factor;
      if (sol_D[cell_no].spec[count_spec].Dp > diffc_max)
	diffc_max = sol_D[cell_no].spec[count_spec].Dp;
      count_spec++;
    }
  }
  sol_D[cell_no].spec =
    (struct spec *) PHRQ_realloc (sol_D[cell_no].spec,
				  (size_t) count_spec * sizeof (struct spec));
  if (sol_D[cell_no].spec == NULL)
    malloc_error ();

  sol_D[cell_no].count_spec = count_spec;

  return (OK);
}

/* ---------------------------------------------------------------------- */
int
sort_species_name (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
  const struct species_list *nptr1, *nptr2;

  nptr1 = (const struct species_list *) ptr1;
  nptr2 = (const struct species_list *) ptr2;

  return (strcmp (nptr1->s->name, nptr2->s->name));
}

/* ---------------------------------------------------------------------- */
int
multi_Dstag (int mobile_cell)
/* ---------------------------------------------------------------------- */
{
/* basic scheme follows multi_D, but uses the mix factors predefined with MIX
 * 1. determine mole transfer (mol/s) of solute species for the interface between 2 cells.
 * 2. sum up as mole transfer of master_species
 * 3. add moles of master_species to the 2 cells
 *      NOTE. Define the water content of stagnant cells relative to the
 *      mobile cell (with, for example, 1 kg water)
 *      Define properties of each interface only 1 time with MIX.
 */
  int icell, jcell, i, j, k, l, n, length, length2;
  char *ptr;
  char token[MAX_LENGTH];
  struct m_s
  {
    char *name;
    LDBLE tot;
  } *m_s;
  int count_m_s;
  LDBLE tot_h, tot_o, mixf, temp;

  m_s =
    (struct m_s *) PHRQ_malloc ((size_t) count_elements *
				sizeof (struct m_s));
  if (m_s == NULL)
    malloc_error ();

  for (n = 0; n < stag_data->count_stag; n++)
  {
    icell = mobile_cell + 1 + n * count_cells;
    if (n == 0)
      icell -= 1;
/*
 *      find the mix ptr for icell and go along the cells that mix with it
 */
    use.mix_ptr = mix_search (icell, &use.n_mix, FALSE);
    if (use.mix_ptr == NULL)
      continue;
    for (i = 0; i < use.mix_ptr->count_comps; i++)
    {
      if ((jcell = use.mix_ptr->comps[i].n_solution) <= icell)
	continue;
      mixf = use.mix_ptr->comps[i].fraction;
      if (mcd_substeps > 1)
	mixf /= nmix;
/*
 * 1. obtain J_ij...
 */
      find_Jstag (icell, jcell, mixf);
/*
 * 2. sum up the primary or secondary master_species from all the solute species
 *      H and O go in total_h and total_o
 */
      tot_h = tot_o = 0.0;
      count_m_s = 0;
      for (j = 0; j < J_ij_count_spec; j++)
      {
	ptr = J_ij[j].name;
	count_elts = 0;
	get_elts_in_species (&ptr, 1);
	for (k = 0; k < count_elts; k++)
	{
	  if (strcmp (elt_list[k].elt->name, "H") == 0)
	    tot_h += elt_list[k].coef * J_ij[j].tot;
	  else if (strcmp (elt_list[k].elt->name, "O") == 0)
	    tot_o += elt_list[k].coef * J_ij[j].tot;
	  else
	  {
	    for (l = 0; l < count_m_s; l++)
	    {
	      length = (int) strlen (elt_list[k].elt->name);
	      if (strncmp (m_s[l].name, elt_list[k].elt->name, length) == 0)
	      {
		m_s[l].tot += elt_list[k].coef * J_ij[j].tot;
		break;
	      }
	    }
	    if (l == count_m_s)
	    {
	      m_s[l].name = string_hsave (elt_list[k].elt->name);
	      m_s[l].tot = elt_list[k].coef * J_ij[j].tot;
	      count_m_s++;
	    }
	  }
	}
      }
/*
 * timestep is in mixf.
 * NOTE. The timestep calculated in init_mix for MCD (by PHREEQC) must be equal or smaller than
 *        the timestep taken (by the user) for calculating mixf in MIX.
 *        Make this timestep small enough, consider the largest Dw in phreeqd.dat (usually H+).
 *        Dw used for calculating mixf must be given as default_Dw in the input file.
 */
/*
 * 3. find the solutions, add or subtract the moles...
 */
      use.solution_ptr = solution_bsearch (icell, &use.n_solution, FALSE);
      use.solution_ptr->total_h -= tot_h;
      use.solution_ptr->total_o -= tot_o;
      use.solution_ptr->cb -= J_ij_sum;
      for (l = 0; l < count_m_s; l++)
      {
	temp = 0.0;
	length = (int) strlen (m_s[l].name);
	for (j = 0; use.solution_ptr->totals[j].description != NULL; j++)
	{
          length2 = (int) (size_t) strcspn(use.solution_ptr->totals[j].description, "(");
          if (strncmp
	      (m_s[l].name, use.solution_ptr->totals[j].description,
	       length) == 0 && length == length2)
	  {
	    if (use.solution_ptr->totals[j].moles < m_s[l].tot)
	    {
	      temp = use.solution_ptr->totals[j].moles;
	      use.solution_ptr->totals[j].moles = 0;
	      /* see if other redox states have more moles... */
	      for (k = 1; use.solution_ptr->totals[j + k].description != NULL;
		   k++)
	      {
                length2 = (int) (size_t) strcspn(use.solution_ptr->totals[j + k].description, "(");
		if (strncmp
		    (m_s[l].name, use.solution_ptr->totals[j + k].description,
		     length) == 0 && length == length2)
		{
		  temp += use.solution_ptr->totals[j + k].moles;
		  if (temp < m_s[l].tot)
		  {
		    use.solution_ptr->totals[j + k].moles = 0;
		  }
		  else
		  {
		    use.solution_ptr->totals[j + k].moles = temp - m_s[l].tot;
		    temp = 0.0;
		    break;
		  }
		}
	      }
	      if (temp != 0.0 && m_s[l].tot - temp > 1e-12)
	      {
		sprintf (token,
			 "Negative concentration in MCD: added %.1e moles %s in cell %d.",
			 m_s[l].tot - temp, m_s[l].name, icell);
		warning_msg (token);
	      }
	    }
	    else
	      use.solution_ptr->totals[j].moles -= m_s[l].tot;
	    break;
	  }
	}
	if (use.solution_ptr->totals[j].description == NULL)
	{
	  use.solution_ptr->totals =
	    (struct conc *) PHRQ_realloc (use.solution_ptr->totals,
					  (size_t) (j +
						    2) *
					  sizeof (struct conc));
	  use.solution_ptr->totals[j].description =
	    string_hsave (m_s[l].name);
	  use.solution_ptr->totals[j].moles = -m_s[l].tot;
	  if (use.solution_ptr->totals[j].moles < 0)
	  {
	    if (use.solution_ptr->totals[j].moles < -1e-12)
	    {
	      sprintf (token,
		     "Negative concentration in MCD: added %.2e moles %s in cell %d",
		     -use.solution_ptr->totals[j].moles, m_s[l].name, icell);
	      warning_msg (token);
	    }
	    use.solution_ptr->totals[j].moles = 0;
	  }
	  use.solution_ptr->totals[j + 1].description = NULL;
	}
      }
      use.solution_ptr = solution_bsearch (jcell, &use.n_solution, FALSE);
      use.solution_ptr->total_h += tot_h;
      use.solution_ptr->total_o += tot_o;
      use.solution_ptr->cb += J_ij_sum;
      for (l = 0; l < count_m_s; l++)
      {
	temp = 0.0;
	length = (int) strlen (m_s[l].name);
	for (j = 0; use.solution_ptr->totals[j].description != NULL; j++)
	{
          length2 = (int) (size_t) strcspn(use.solution_ptr->totals[j].description, "(");
	  if (strncmp
	      (m_s[l].name, use.solution_ptr->totals[j].description,
	       length) == 0 && length == length2)
	  {
	    if (use.solution_ptr->totals[j].moles < -m_s[l].tot)
	    {
	      temp = use.solution_ptr->totals[j].moles;
	      use.solution_ptr->totals[j].moles = 0;
	      /* see if other redox states have more moles... */
	      for (k = 1; use.solution_ptr->totals[j + k].description != NULL;
		   k++)
	      {
                length2 = (int) (size_t) strcspn(use.solution_ptr->totals[j + k].description, "(");
		if (strncmp
		    (m_s[l].name, use.solution_ptr->totals[j + k].description,
		     length) == 0 && length == length2)
		{
		  temp += use.solution_ptr->totals[j + k].moles;
		  if (temp < -m_s[l].tot)
		  {
		    use.solution_ptr->totals[j + k].moles = 0;
		  }
		  else
		  {
		    use.solution_ptr->totals[j + k].moles = temp + m_s[l].tot;
		    temp = 0.0;
		    break;
		  }
		}
	      }
	      if (temp != 0.0 && -m_s[l].tot - temp > 1e-12)
	      {
		sprintf (token,
			 "Negative concentration in MCD: added %.3e moles %s in cell %d",
			 -m_s[l].tot - temp, m_s[l].name, jcell);
		warning_msg (token);
	      }
	    }
	    else
	      use.solution_ptr->totals[j].moles += m_s[l].tot;
	    break;
	  }
	}
	if (use.solution_ptr->totals[j].description == NULL)
	{
	  use.solution_ptr->totals =
	    (struct conc *) PHRQ_realloc (use.solution_ptr->totals,
					  (size_t) (j +
						    2) *
					  sizeof (struct conc));
	  use.solution_ptr->totals[j].description =
	    string_hsave (m_s[l].name);
	  use.solution_ptr->totals[j].moles = m_s[l].tot;
	  if (use.solution_ptr->totals[j].moles < 0)
	  {
	    if (use.solution_ptr->totals[j].moles < -1e-12)
	    {
	      sprintf (token,
		"Negative concentration in MCD: added %.4e moles %s in cell %d",
		-use.solution_ptr->totals[j].moles, m_s[l].name, jcell);
	      warning_msg (token);
	    }
	    use.solution_ptr->totals[j].moles = 0;
	  }
	  use.solution_ptr->totals[j + 1].description = NULL;
	}
      }
    }
  }
  m_s = (struct m_s *) free_check_null (m_s);
  J_ij = (struct J_ij *) free_check_null (J_ij);
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
find_Jstag (int icell, int jcell, LDBLE mixf)
/* ---------------------------------------------------------------------- */
{
/* mole transfer (mol/s) of the individual master_species:
 * Eqn 1:
 * J_ij = mixf_ij * (-D_i*grad(a) + D_i*z_i*c_i * SUM(D_i*z_i*grad(a)) / SUM(D_i*(z_i)^2*c_i))
 *              mixf_ij = mixf / (Dw * init_pf) * new_por / init_por
 *              mixf is defined in MIX; Dw is default multicomponent diffusion coefficient;
 *              init_pf equals multi_Dpor^multi_Dn;
 */
  int i, i_max, j, j_max, k, l;
  LDBLE ddlm, aq1, aq2;
  LDBLE *grad, *D, *z, *Dz, *Dzc, *Dzc_dl, Dz2c, Dz2c_dl, c;
  int *o_c, only_counter;
  struct surface *s_ptr1, *s_ptr2;
  struct surface_charge *s_charge_ptr, *s_charge_ptr1, *s_charge_ptr2;
  LDBLE dl_s, dl_aq1, dl_aq2, c_dl, visc1, visc2, temp;
  char token[MAX_LENGTH], token1[MAX_LENGTH];

  if (cell_data[icell - 1].por < multi_Dpor_lim
      || cell_data[jcell - 1].por < multi_Dpor_lim)
  {
    J_ij_count_spec = 0;
    J_ij_sum = 0;
    return (OK);
  }

  mixf *=
    (cell_data[icell - 1].por <=
     cell_data[jcell - 1].por ? cell_data[icell - 1].por : cell_data[jcell -
								     1].por);
  mixf /= (default_Dw * pow (multi_Dpor, multi_Dn) * multi_Dpor);
  aq1 = solution_bsearch (icell, &i, TRUE)->mass_water;
  aq2 = solution_bsearch (jcell, &i, TRUE)->mass_water;

/*
 * check if DL calculations must be made...
 */
  s_charge_ptr1 = s_charge_ptr2 = NULL;
  s_ptr1 = s_ptr2 = NULL;
  dl_s = dl_aq1 = dl_aq2 = temp = 0.0;
  visc1 = visc2 = 1.0;
  only_counter = FALSE;

  s_ptr1 = surface_bsearch (icell, &i);
  if (s_ptr1 != NULL)
  {
    if (s_ptr1->dl_type != NO_DL)
    {
      if (s_ptr1->only_counter_ions)
	only_counter = TRUE;
      /* find the one (and only one...) immobile surface comp with DL... */
      for (i = 0; i < s_ptr1->count_comps; i++)
      {
	if (s_ptr1->comps[i].Dw == 0)
	{
	  s_charge_ptr1 = &s_ptr1->charge[s_ptr1->comps[i].charge];
	  dl_aq1 = s_charge_ptr1->mass_water;
	  visc1 = s_ptr1->DDL_viscosity;
	  /* check for more comps with Dw = 0 */
	  for (j = i + 1; j < s_ptr1->count_comps; j++)
	  {
	    if (s_ptr1->comps[j].Dw == 0
		&& s_ptr1->comps[j].charge != s_ptr1->comps[i].charge)
	    {
	      k = (int) strcspn (s_ptr1->comps[i].formula, "_");
	      strncpy (token1, s_ptr1->comps[i].formula, k);
	      token1[k] = '\0';
	      sprintf (token,
		       "MCD found more than 1 fixed surface with a DDL, used %s.",
		       token1);
	      warning_msg (token);
	      break;
	    }
	  }
	  break;
	}
      }
    }
  }
  s_ptr2 = surface_bsearch (jcell, &i);
  if (s_ptr2 != NULL)
  {
    if (s_ptr2->dl_type != NO_DL)
    {
      if (s_ptr2->only_counter_ions)
	only_counter = TRUE;
      for (i = 0; i < s_ptr2->count_comps; i++)
      {
	if (s_ptr2->comps[i].Dw == 0)
	{
	  s_charge_ptr2 = &s_ptr2->charge[s_ptr2->comps[i].charge];
	  dl_aq2 = s_charge_ptr2->mass_water;
	  visc2 = s_ptr2->DDL_viscosity;
	  /* check for more comps with Dw = 0 */
	  for (j = i + 1; j < s_ptr2->count_comps; j++)
	  {
	    if (s_ptr2->comps[j].Dw == 0
		&& s_ptr2->comps[j].charge != s_ptr2->comps[i].charge)
	    {
	      k = (int) strcspn (s_ptr2->comps[i].formula, "_");
	      strncpy (token1, s_ptr2->comps[i].formula, k);
	      token1[k] = '\0';
	      sprintf (token,
		       "MCD found more than 1 fixed surface with a DDL, used %s.",
		       token1);
	      warning_msg (token);
	      break;
	    }
	  }
	  break;
	}
      }
    }
  }
  /* in each cell: DL surface = mass_water_DL / (cell_length * tortuosity)
     free pore surface = mass_water_free / (cell_length * tortuosity)
     determine DL surface as a fraction of the total pore surface... */
  if (dl_aq1 > 0)
    dl_s = dl_aq1 / (dl_aq1 + aq1);
  if (dl_aq2 > 0)
    temp = dl_aq2 / (dl_aq2 + aq2);
  if (dl_aq1 > 0 && dl_aq2 > 0)
    /* average the 2... */
    dl_s = (dl_s + temp) / 2;
  else if (dl_aq2 > 0)
    /* there is one DL surface... */
    dl_s = temp;

/*
 * malloc sufficient space...
 */
  k = sol_D[icell].count_spec + sol_D[jcell].count_spec;

  J_ij = (struct J_ij *) free_check_null (J_ij);
  J_ij = (struct J_ij *) PHRQ_malloc ((size_t) k * sizeof (struct J_ij));
  if (J_ij == NULL)
    malloc_error ();

  grad = (LDBLE *) PHRQ_malloc ((size_t) k * sizeof (LDBLE));
  if (grad == NULL)
    malloc_error ();

  D = (LDBLE *) PHRQ_malloc ((size_t) k * sizeof (LDBLE));
  if (D == NULL)
    malloc_error ();

  z = (LDBLE *) PHRQ_malloc ((size_t) k * sizeof (LDBLE));
  if (z == NULL)
    malloc_error ();

  Dz = (LDBLE *) PHRQ_malloc ((size_t) k * sizeof (LDBLE));
  if (Dz == NULL)
    malloc_error ();

  Dzc = (LDBLE *) PHRQ_malloc ((size_t) k * sizeof (LDBLE));
  if (Dzc == NULL)
    malloc_error ();

  Dzc_dl = (LDBLE *) PHRQ_malloc ((size_t) k * sizeof (LDBLE));
  if (Dzc_dl == NULL)
    malloc_error ();

  o_c = (int *) PHRQ_malloc ((size_t) k * sizeof (int));
  if (o_c == NULL)
    malloc_error ();
  Dz2c = Dz2c_dl = 0.0;

/*
 * coefficients in Eqn (1)...
 */
  i = j = k = 0;
  i_max = sol_D[icell].count_spec;
  j_max = sol_D[jcell].count_spec;

  while (i < i_max || j < j_max)
  {
    if (j == j_max
	|| (strcmp (sol_D[icell].spec[i].name, sol_D[jcell].spec[j].name) < 0
	    && i < i_max))
    {
      /* species 'name' is only in icell */
      J_ij[k].name = string_hsave (sol_D[icell].spec[i].name);
      D[k] = sol_D[icell].spec[i].Dp;
      z[k] = sol_D[icell].spec[i].z;
      Dz[k] = D[k] * z[k];
      Dzc[k] = Dz[k] * sol_D[icell].spec[i].c / 2;
      if (dl_s > 0)
      {
	o_c[k] = 1;
	s_charge_ptr = (dl_aq1 > 0) ? s_charge_ptr1 : s_charge_ptr2;
	for (l = 0; l < s_charge_ptr->count_g; l++)
	{
	  if (equal (s_charge_ptr->g[l].charge, z[k], G_TOL) == TRUE)
	  {
	    if (only_counter &&
		((s_charge_ptr->la_psi < 0 && z[k] < 0)
		 || (s_charge_ptr->la_psi > 0 && z[k] > 0)))
	    {
	      o_c[k] = 0;
	      Dzc_dl[k] = 0;
	    }
	    else
	    {
	      if (dl_aq1 > 0)
		Dzc_dl[k] =
		  Dz[k] * sol_D[icell].spec[i].c / 2 * (1 +
							  s_charge_ptr->g[l].
							  g * aq1 / dl_aq1);
	      else
		Dzc_dl[k] = Dz[k] * sol_D[icell].spec[i].c / 2;
	    }
	    break;
	  }
	}
	Dz2c_dl += Dzc_dl[k] * z[k];
      }
      Dz2c += Dzc[k] * z[k];
      grad[k] = -sol_D[icell].spec[i].a;
      if (i < i_max)
	i++;
    }
    else if (i == i_max
	     || (strcmp (sol_D[icell].spec[i].name, sol_D[jcell].spec[j].name)
		 > 0 && j < j_max))
    {
      /* species 'name' is only in jcell */
      J_ij[k].name = string_hsave (sol_D[jcell].spec[j].name);
      D[k] = sol_D[jcell].spec[j].Dp;
      z[k] = sol_D[jcell].spec[j].z;
      Dz[k] = D[k] * z[k];
      Dzc[k] = Dz[k] * sol_D[jcell].spec[j].c / 2;
      if (dl_s > 0)
      {
	o_c[k] = 1;
	s_charge_ptr = (dl_aq2 > 0) ? s_charge_ptr2 : s_charge_ptr1;
	for (l = 0; l < s_charge_ptr->count_g; l++)
	{
	  if (equal (s_charge_ptr->g[l].charge, z[k], G_TOL) == TRUE)
	  {
	    if (only_counter &&
		((s_charge_ptr->la_psi < 0 && z[k] < 0)
		 || (s_charge_ptr->la_psi > 0 && z[k] > 0)))
	    {
	      o_c[k] = 0;
	      Dzc_dl[k] = 0;
	    }
	    else
	    {
	      if (dl_aq2 > 0)
		Dzc_dl[k] =
		  Dz[k] * sol_D[jcell].spec[j].c / 2 * (1 +
							      s_charge_ptr->
							      g[l].g * aq2 /
							      dl_aq2);
	      else
		Dzc_dl[k] = Dz[k] * sol_D[jcell].spec[j].c / 2;
	    }
	    break;
	  }
	}
	Dz2c_dl += Dzc_dl[k] * z[k];
      }
      Dz2c += Dzc[k] * z[k];
      grad[k] = sol_D[jcell].spec[j].a;
      if (j < j_max)
	j++;
    }
    else if (strcmp (sol_D[icell].spec[i].name, sol_D[jcell].spec[j].name) ==
	     0)
    {
      /* species 'name' is in both cells */
      J_ij[k].name = string_hsave (sol_D[icell].spec[i].name);
      if (sol_D[icell].spec[i].Dp == 0 || sol_D[jcell].spec[j].Dp == 0)
	D[k] = 0.0;
      else
	D[k] = (sol_D[icell].spec[i].Dp + sol_D[jcell].spec[j].Dp) / 2;
      z[k] = sol_D[icell].spec[i].z;
      Dz[k] = D[k] * z[k];
      Dzc[k] = Dz[k] * (sol_D[icell].spec[i].c + sol_D[jcell].spec[j].c) / 2;
/*                      Dzc[k] = Dz[k] * (sol_D[icell].spec[i].c > sol_D[jcell].spec[j].c ? sol_D[icell].spec[i].c : sol_D[jcell].spec[j].c);
 */
      if (dl_s > 0)
      {
	c_dl = 0.0;
	o_c[k] = 1;
	if (dl_aq1 > 0)
	{
	  for (l = 0; l < s_charge_ptr1->count_g; l++)
	  {
	    if (equal (s_charge_ptr1->g[l].charge, z[k], G_TOL) == TRUE)
	    {
	      if (only_counter &&
		  ((s_charge_ptr1->la_psi < 0 && z[k] < 0)
		   || (s_charge_ptr1->la_psi > 0 && z[k] > 0)))
	      {
		o_c[k] = 0;
		c_dl = 0;
	      }
	      else
		c_dl =
		  sol_D[icell].spec[i].c / 2 * (1 +
						s_charge_ptr1->g[l].g * aq1 /
						dl_aq1);
	      break;
	    }
	  }
	}
	else
	  c_dl = sol_D[icell].spec[i].c / 2;

	if (dl_aq2 > 0)
	{
	  for (l = 0; l < s_charge_ptr2->count_g; l++)
	  {
	    if (equal (s_charge_ptr2->g[l].charge, z[k], G_TOL) == TRUE)
	    {
	      if (only_counter &&
		  ((s_charge_ptr2->la_psi < 0 && z[k] < 0)
		   || (s_charge_ptr2->la_psi > 0 && z[k] > 0)))
	      {
		o_c[k] = 0;
		c_dl = 0;
	      }
	      if (o_c[k] == 1)
		c_dl +=
		  sol_D[jcell].spec[j].c / 2 * (1 +
						s_charge_ptr2->g[l].g * aq2 /
						dl_aq2);
	      break;
	    }
	  }
	}
	else if (o_c[k] == 1)
	  c_dl += sol_D[jcell].spec[j].c / 2;

	Dzc_dl[k] = Dz[k] * c_dl;
	Dz2c_dl += Dzc_dl[k] * z[k];
      }
      Dz2c += Dzc[k] * z[k];
      grad[k] = (sol_D[jcell].spec[j].c - sol_D[icell].spec[i].c);
      ddlm = sol_D[jcell].spec[j].lm - sol_D[icell].spec[i].lm;
      if (fabs (ddlm) > 1e-10)
	grad[k] *=
	  (1 + (sol_D[jcell].spec[j].lg - sol_D[icell].spec[i].lg) / ddlm);
      if (i < i_max)
	i++;
      if (j < j_max)
	j++;
    }
    k++;
  }
/*
 * fill in J_ij...
 */
  if (Dz2c == 0)
    k = 0;
  J_ij_count_spec = k;
  J_ij_sum = 0;
  c = c_dl = 0.0;
  for (j = 0; j < k; j++)
  {
    c += Dz[j] * grad[j];
    c_dl += o_c[j] * Dz[j] * grad[j];
  }
  for (i = 0; i < k; i++)
  {
    J_ij[i].tot = -D[i] * grad[i] + c * Dzc[i] / Dz2c;
    J_ij[i].tot *= (1 - dl_s);
    if (Dz2c_dl > 0)
    {
      temp =
	(-D[i] * grad[i] + c_dl * Dzc_dl[i] / Dz2c_dl) * (2 / (visc1 + visc2));
      if ((J_ij[i].tot <= 0 && temp <= 0) || (J_ij[i].tot > 0 && temp > 0))
      {
	J_ij[i].tot += o_c[i] * temp * dl_s;
      }
    }
    J_ij[i].tot *= mixf;
    J_ij_sum += z[i] * J_ij[i].tot;
  }

/* appt 3 July 07, improved convergence without transporting charge imbalance */
  J_ij_sum = 0;
  D = (LDBLE *) free_check_null (D);
  z = (LDBLE *) free_check_null (z);
  Dz = (LDBLE *) free_check_null (Dz);
  Dzc = (LDBLE *) free_check_null (Dzc);
  Dzc_dl = (LDBLE *) free_check_null (Dzc_dl);
  grad = (LDBLE *) free_check_null (grad);
  o_c = (int *) free_check_null (o_c);

  return (OK);
}

/* ---------------------------------------------------------------------- */
int
disp_surf (LDBLE DDt)
/* ---------------------------------------------------------------------- */
/*
 *  Disperse/diffuse surfaces:
 *       surface[n1] = SS(mixf * surface[n2]) + (1 - SS(mixf) * surface[i1])
 *  where SS is Sum of, f2 is a mixing factor.
 *  Implementation:
 *  Mobile surface comps and charges are mixed face by face and 1 by 1 in surface[n1]:
  Step (from cell 1 to count_cells + 1):
 *  1. Surface n1 is made a copy of cell[i1]'s surface, if it exists, or else
 *       b. a copy of the first encountered mobile surface[n2] from cell i2 = i1 - 1.
 *  2  a. Column surfaces are mixed by dispersion/diffusion
 *       b. Column surfaces are mixed by MCD (if multi_Dflag is true)
 *  3. Surfaces n1 and n2 are stored in a temporary surface, with n_user = i1 or i2.
 *  4. When done for all cells, new surfaces are copied into the cells.
 *       NOTE... The surfaces' diffuse_layer, edl and phases/kinetics relations must be identical,
 *                       but only mobile surface_comp's (Dw > 0) and their surface_charge's are transported.
 */
{
  int i, i1, i2, j, k, k1, n, n1, n2, temp_count_surface;
  int charge_done, surf1, surf2;
  LDBLE f1, f2, mixf, mixf_store, mcd_mixf;
  LDBLE lav, A_ij, por, Dp1, Dp2;
  struct mix *mix_ptr;
  struct surface *surface_ptr1, *surface_ptr2;
  LDBLE viscos_f;
/*
 * temperature and viscosity correction for MCD coefficient, D_T = D_298 * Tk * viscos_298 / (298 * viscos)
 */
  viscos_f =
    pow (10,
	 -(1.37023 * (tc_x - 20) +
	   0.000836 * (tc_x - 20) * (tc_x - 20)) / (109 + tc_x));
  viscos_f = tk_x * 0.88862 / (298.15 * viscos_f);

  n1 = count_surface;
  n2 = n1 + 1;
  temp_count_surface = count_surface + 2;
  space ((void **) ((void *) &surface), temp_count_surface, &max_surface,
	 sizeof (struct surface));
  surface[n1].count_comps = surface[n1 + 1].count_comps = 0;

  for (i1 = 1; i1 <= count_cells + 1; i1++)
  {

    if (i1 <= count_cells && cell_data[i1 - 1].por < multi_Dpor_lim)
      continue;

    if (i1 == 1 && bcon_first != 1)
      continue;
    if (i1 == count_cells + 1 && bcon_last != 1)
      continue;

    i2 = i1 - 1;
    if (i2 > 0 && cell_data[i2 - 1].por < multi_Dpor_lim)
      continue;
/*
 * step 1. define surface n1 from cell i1, if it exists...
 */
    surface[n1].n_user = surface[n1 + 1].n_user = -99;
    if ((surface_ptr1 = surface_bsearch (i1, &n)) != NULL)
      surface_copy (surface_ptr1, &surface[n1], i1);

    surface_ptr2 = surface_bsearch (i2, &n);

    surf1 = surf2 = 0;
    if (surface_ptr2 != NULL)
    {
      if (surface_ptr2->transport == TRUE)
	surf2 = 1;
    }
    if (surface_ptr1 != NULL)
    {
      if (surface_ptr1->transport == TRUE)
	surf1 = 1;
    }

    mixf = mixf_store = 0;
    if (surf1 || surf2)
    {
/*
 * step 2a. Dispersive mixing of mobile surfaces from column cells i2 < i1...
 */
      if (i1 <= count_cells)
	mix_ptr = &mix[count_mix - count_cells + i1 - 1];
      else
	mix_ptr = &mix[count_mix - 1];

      for (j = 0; j < mix_ptr->count_comps; j++)
      {
	if (i1 <= count_cells)
	{
	  if (mix_ptr->comps[j].n_solution == i2)
	  {
	    mixf = mixf_store = mix_ptr->comps[j].fraction;
	    break;
	  }
	}
	else if (mix_ptr->comps[j].n_solution == i1)
	{
	  mixf = mixf_store = mix_ptr->comps[j].fraction;
	  break;
	}
      }

      /* step 1b. If cell i1 has no surface, get the mobile comps from surface i2... */
      if (surface[n1].n_user == -99)
      {
	mobile_surface_copy (surface_ptr2, &surface[n1], i1, FALSE);
	/* limit charges to 1... */
	if (surface[n1].count_charge > 1)
	{
	  if (sum_surface_comp
	      (&surface[n1], 0, &surface[n1], 0, 1, &surface[n1],
	       surface[n1].comps[0].Dw) == ERROR)
	    return (ERROR);
	}
	f1 = 0;
      }
      else
	f1 = 1;
      /* find the (possibly modified) surface in the previous cell i2... */
      f2 = 1;
      if (i2 > 0 || bcon_first == 1)
      {
	for (k = count_surface + 2; k < temp_count_surface; k++)
	{
	  if (surface[k].n_user == i2)
	  {
	    n2 = k;
	    break;
	  }
	}
	/* if not found... */
	if (k == temp_count_surface)
	{
	  n2 = n1 + 1;
	  /* copy it from surface_ptr2... */
	  if (surface_ptr2 != NULL)
	    surface_copy (surface_ptr2, &surface[n2], i2);
	  else
	  {
	    /* or make it a mobile copy of the surface in cell i1... */
	    mobile_surface_copy (surface_ptr1, &surface[n2], i2, FALSE);
	    /* limit charges to 1... */
	    if (surface[n2].count_charge > 1)
	    {
	      if (sum_surface_comp
		  (&surface[n2], 0, &surface[n2], 0, 1, &surface[n2],
		   surface[n2].comps[0].Dw) == ERROR)
		return (ERROR);
	    }
	    f2 = 0;
	  }
	}
      }

/*
 * For MCD, step 2b. Find physical dimensions and porosity of the cells...
 */
      if (i1 == 1)
      {
	por = cell_data[0].por;
	lav = cell_data[0].length / 2;
	A_ij =
	  solution_bsearch (1, &n, TRUE)->mass_water / cell_data[0].length;
      }
      else if (i1 == count_cells + 1)
      {
	por = cell_data[count_cells - 1].por;
	lav = cell_data[count_cells - 1].length / 2;
	A_ij =
	  solution_bsearch (count_cells, &n,
			    TRUE)->mass_water / cell_data[count_cells -
							  1].length;
      }
      else
      {
	por = cell_data[i2 - 1].por;
	lav = (cell_data[i1 - 1].length + cell_data[i2 - 1].length) / 2;
	A_ij =
	  solution_bsearch (i1, &n,
			    TRUE)->mass_water / (cell_data[i1 -
							   1].length *
						 cell_data[i1 - 1].por);
	A_ij +=
	  solution_bsearch (i2, &n,
			    TRUE)->mass_water / (cell_data[i2 -
							   1].length *
						 cell_data[i2 - 1].por);
	A_ij /= 2;
	A_ij *=
	  (cell_data[i1 - 1].por <
	   cell_data[i2 - 1].por ? cell_data[i1 - 1].por : cell_data[i2 -
								     1].por);
      }

      /* mix in comps with the same charge structure... */
      if (surf2)
      {
	for (k = 0; k < surface_ptr2->count_comps; k++)
	{
	  if (surface_ptr2->comps[k].Dw == 0)
	    continue;
	  charge_done = FALSE;
	  for (k1 = 0; k1 < k; k1++)
	  {
	    if (surface_ptr2->comps[k1].charge ==
		surface_ptr2->comps[k].charge)
	    {
	      charge_done = TRUE;
	      break;
	    }
	  }
	  if (charge_done)
	    continue;

	  /* Define the MCD diffusion coefficient... */
	  mcd_mixf = 0;
	  if (multi_Dflag)
	  {
	    Dp2 = surface_ptr2->comps[k].Dw * pow (por, multi_Dn);
	    Dp1 = 0;
	    if (surface_ptr1 != NULL && surface_ptr1->transport)
	    {
	      for (k1 = 0; k1 < surface_ptr1->count_comps; k1++)
	      {
		if (elt_list_compare
		    (surface_ptr1->comps[k1].totals,
		     surface_ptr2->comps[k].totals) != 0)
		  continue;
		Dp1 =
		  surface_ptr1->comps[k].Dw * pow (cell_data[i1 - 1].por,
						   multi_Dn);
		break;
	      }
	    }
	    if (Dp1 > 0)
	      Dp2 = (Dp2 + Dp1) / 2;

	    /* and the mixing factor... */
	    mcd_mixf = Dp2 * viscos_f * A_ij * DDt / lav;
	  }
	  mixf = mixf_store + mcd_mixf;

	  if (mixf < 1e-8)
	    mixf = 0;
	  if (mixf > 0.99999999)
	    mixf = 0.99999999;
	  if (i1 <= count_cells)
	  {
	    if (sum_surface_comp
		(&surface[n1], f1, surface_ptr2, k, mixf, &surface[n1],
		 surface_ptr2->comps[k].Dw) == ERROR)
	      return (ERROR);
	    f1 = 1;
	  }
	  if (i2 > 0)
	  {
	    if (sum_surface_comp
		(&surface[n2], f2, surface_ptr2, k, -mixf, &surface[n2],
		 surface_ptr2->comps[k].Dw) == ERROR)
	      return (ERROR);
	    f2 = 1;
	  }
	  surface[n1].n_user = i1;
	}
      }
      if (surf1)
      {
	for (k = 0; k < surface_ptr1->count_comps; k++)
	{
	  if (surface_ptr1->comps[k].Dw == 0)
	    continue;
	  charge_done = FALSE;
	  for (k1 = 0; k1 < k; k1++)
	  {
	    if (surface_ptr1->comps[k1].charge ==
		surface_ptr1->comps[k].charge)
	    {
	      charge_done = TRUE;
	      break;
	    }
	  }
	  if (charge_done)
	    continue;

	  /* Define the MCD diffusion coefficient... */
	  mcd_mixf = 0;
	  if (multi_Dflag)
	  {
	    if (i1 <= count_cells)
	      Dp1 =
		surface_ptr1->comps[k].Dw * pow (cell_data[i1 - 1].por,
						 multi_Dn);
	    else
	      Dp1 = surface_ptr1->comps[k].Dw * pow (por, multi_Dn);
	    Dp2 = 0;
	    if (surface_ptr2 != NULL && surface_ptr2->transport)
	    {
	      for (k1 = 0; k1 < surface_ptr2->count_comps; k1++)
	      {
		if (elt_list_compare
		    (surface_ptr2->comps[k1].totals,
		     surface_ptr1->comps[k].totals) != 0)
		  continue;
		Dp2 = surface_ptr2->comps[k].Dw * pow (por, multi_Dn);
		break;
	      }
	    }
	    if (Dp2 > 0)
	      Dp1 = (Dp1 + Dp2) / 2;

	    /* and the mixing factor... */
	    mcd_mixf = Dp1 * viscos_f * A_ij * DDt / lav;
	  }
	  mixf = mixf_store + mcd_mixf;

	  if (mixf < 1e-8)
	    mixf = 0;
	  if (mixf > 0.99999999)
	    mixf = 0.99999999;
	  if (i2 > 0)
	  {
	    if (sum_surface_comp
		(&surface[n2], f2, surface_ptr1, k, mixf, &surface[n2],
		 surface_ptr1->comps[k].Dw) == ERROR)
	      return (ERROR);
	    f2 = 1;
	  }
	  if (i1 <= count_cells)
	  {
	    if (sum_surface_comp
		(&surface[n1], f1, surface_ptr1, k, -mixf, &surface[n1],
		 surface_ptr1->comps[k].Dw) == ERROR)
	      return (ERROR);
	    f1 = 1;
	  }
	  surface[n2].n_user = i2;
	}
      }
    }

/*
 *  Step 3. copy surface[n1] and [n2] in a new temporary surface...
 */
    if (surface[n1].n_user == -99)
      continue;

    n = temp_count_surface++;

    space ((void **) ((void *) &surface), temp_count_surface, &max_surface,
	   sizeof (struct surface));
    surface_copy (&surface[n1], &surface[n], i1);

    surface[n1].n_user = -99;
    surface[n].n_user = surface[n].n_user_end = i1;
    surface_free (&surface[n1]);
    surface[n1].count_comps = 0;
    if (n2 == n1 + 1 && surface[n2].n_user == i2)
    {
      n = temp_count_surface++;

      space ((void **) ((void *) &surface), temp_count_surface, &max_surface,
	     sizeof (struct surface));
      surface_copy (&surface[n2], &surface[n], i2);

      surface[n2].n_user = -99;
      surface[n].n_user = surface[n].n_user_end = i2;
      surface_free (&surface[n1 + 1]);
      surface[n1 + 1].count_comps = 0;
    }
  }
/*
 * Step 4. Dispersion/diffusion is done. New surfaces can be copied in the cell's surface...
 */
  n2 = 0;

  for (j = count_surface + 2; j < temp_count_surface; j++)
  {
    i = surface[j].n_user;
    if ((i == 0 && bcon_first == 1)
	|| (i == count_cells + 1 && bcon_last == 1))
    {
      surface_free (&surface[j]);
      surface[j].count_comps = 0;
      continue;
    }
    if (i >= 0 && i <= 1 + count_cells * (1 + stag_data->count_stag))
    {
      surface_ptr1 = surface_bsearch (i, &n);
      if (surface_ptr1 != NULL)
      {
	surface_free (surface_ptr1);
	surface_copy (&surface[j], surface_ptr1, i);
      }
      else
      {
	surface_copy (&surface[j], &surface[count_surface + n2], i);
	n2++;
      }
    }
    surface_free (&surface[j]);
    surface[j].count_comps = 0;
  }
  count_surface += n2;
  if (n2 > 0)
    qsort (surface, (size_t) count_surface, (size_t) sizeof (struct surface),
	   surface_compare);

  return (OK);
}

/* ---------------------------------------------------------------------- */
int
sum_surface_comp (struct surface *source1, LDBLE f1, struct surface *source2,
		  int k, LDBLE f2, struct surface *target, LDBLE new_Dw)
/* ---------------------------------------------------------------------- */
{
/*
 *   Takes fraction f1 of the 1st surface, adds fraction f2 of the 2nd surface's comps[k] and its charge.
 *   The result goes in target
 */
  int i, i1, i2;
  int new_n_user, found, charge_in;
  struct surface temp_surface, *surface_ptr;

  struct surface *surface_ptr1, *surface_ptr2;
  char token[MAX_LENGTH];
  int count_comps, count_comps1;
  int count_charge, count_charge1, charge1, charge2;

/*
 *   Find surfaces
 */
  surface_ptr1 = source1;
  if (surface_ptr1 == NULL)
  {
    sprintf (error_string, "Null pointer for surface 1 in sum_surface.");
    error_msg (error_string, CONTINUE);
    input_error++;
    return (ERROR);
  }
  surface_ptr2 = source2;

  if (check_surfaces (source1, source2) == ERROR)
    return (ERROR);
/*
 *   Store data for structure surface
 */
  new_n_user = surface_ptr1->n_user;
  surface_copy (surface_ptr1, &temp_surface, new_n_user);
  sprintf (token, "Copy");
  temp_surface.description =
    (char *) free_check_null (temp_surface.description);
  temp_surface.description = string_duplicate (token);
  temp_surface.solution_equilibria = FALSE;
  temp_surface.n_solution = -99;
/*
 *   Multiply component compositions by f1
 */
  for (i = 0; i < surface_ptr1->count_comps; i++)
  {
    temp_surface.comps[i].moles *= f1;
    temp_surface.comps[i].cb *= f1;
    count_elts = 0;
    add_elt_list (surface_ptr1->comps[i].totals, f1);
    temp_surface.comps[i].totals =
      (struct elt_list *) free_check_null (temp_surface.comps[i].totals);
    temp_surface.comps[i].totals = elt_list_save ();
  }
  if (temp_surface.type != NO_EDL)
  {
    for (i = 0; i < surface_ptr1->count_charge; i++)
    {
      temp_surface.charge[i].grams *= f1;
      temp_surface.charge[i].charge_balance *= f1;
      temp_surface.charge[i].mass_water *= f1;
      temp_surface.charge[i].g = NULL;
      temp_surface.charge[i].count_g = 0;
      count_elts = 0;
      if (surface_ptr1->charge[i].diffuse_layer_totals != NULL)
      {
	add_elt_list (surface_ptr1->charge[i].diffuse_layer_totals, f1);
	temp_surface.charge[i].diffuse_layer_totals =
	  (struct elt_list *) free_check_null (temp_surface.charge[i].
					       diffuse_layer_totals);
	temp_surface.charge[i].g =
	  (struct surface_diff_layer *) free_check_null (temp_surface.
							 charge[i].g);
	temp_surface.charge[i].diffuse_layer_totals = elt_list_save ();
      }
      else
      {
	temp_surface.charge[i].diffuse_layer_totals = NULL;
      }
    }
  }
/*
 *   Add in surface_ptr2
 */

  count_comps = count_comps1 = surface_ptr1->count_comps;
  count_charge = count_charge1 = surface_ptr1->count_charge;

  charge2 = surface_ptr2->comps[k].charge;
  charge_in = FALSE;
  for (i2 = 0; i2 < surface_ptr2->count_comps; i2++)
  {

    if (surface_ptr2->comps[i2].charge != charge2)
      continue;
    /*
     *  handle comps...
     */
    found = FALSE;
    for (i1 = 0; i1 < count_comps1; i1++)
    {

      /* see if first surface contains the 2nd surface's component */
      if (surface_ptr2->comps[i2].formula == surface_ptr1->comps[i1].formula)
      {
	found = TRUE;
	if ((surface_ptr1->comps[i1].phase_name != NULL
	     && surface_ptr2->comps[i2].phase_name == NULL)
	    || (surface_ptr1->comps[i1].phase_name == NULL
		&& surface_ptr2->comps[i2].phase_name != NULL))
	{
	  sprintf (error_string,
		   "Surfaces differ in use of related phases (sites proportional to moles of an equilibrium phase). Can not mix.");
	  error_msg (error_string, CONTINUE);
	  input_error++;
	  return (ERROR);
	}
	else if (surface_ptr1->comps[i1].phase_name != NULL
		 && surface_ptr2->comps[i2].phase_name != NULL
		 && strcmp_nocase (surface_ptr1->comps[i1].phase_name,
				   surface_ptr2->comps[i2].phase_name) != 0)
	{
	  sprintf (error_string,
		   "Surfaces differ in use of related phases (sites proportional to moles of an equilibrium phase). Can not mix.");
	  error_msg (error_string, CONTINUE);
	  input_error++;
	  return (ERROR);
	}
	if ((surface_ptr1->comps[i1].rate_name != NULL
	     && surface_ptr2->comps[i2].rate_name == NULL)
	    || (surface_ptr1->comps[i1].rate_name == NULL
		&& surface_ptr2->comps[i2].rate_name != NULL))
	{
	  sprintf (error_string,
		   "Surfaces differ in use of related rate (sites proportional to moles of a kinetic reactant). Can not mix.");
	  error_msg (error_string, CONTINUE);
	  input_error++;
	  return (ERROR);
	}
	else if (surface_ptr1->comps[i1].rate_name != NULL
		 && surface_ptr2->comps[i2].rate_name != NULL
		 && strcmp_nocase (surface_ptr1->comps[i1].rate_name,
				   surface_ptr2->comps[i2].rate_name) != 0)
	{
	  sprintf (error_string,
		   "Surfaces differ in use of related rates (sites proportional to moles of a kinetic reactant). Can not mix.");
	  error_msg (error_string, CONTINUE);
	  input_error++;
	  return (ERROR);
	}
	temp_surface.comps[i1].moles += surface_ptr2->comps[i2].moles * f2;
	count_elts = 0;
	add_elt_list (temp_surface.comps[i1].totals, 1.0);
	add_elt_list (surface_ptr2->comps[i2].totals, f2);
	if (count_elts > 0)
	{
	  qsort (elt_list, (size_t) count_elts,
		 (size_t) sizeof (struct elt_list), elt_list_compare);
	  elt_list_combine ();
	}
	temp_surface.comps[i1].totals =
	  (struct elt_list *) free_check_null (temp_surface.comps[i1].totals);
	temp_surface.comps[i1].totals = elt_list_save ();
	temp_surface.comps[i1].formula_totals =
	  (struct elt_list *) free_check_null (temp_surface.comps[i1].
					       formula_totals);
	temp_surface.comps[i1].formula_totals = elt_list_save ();
	temp_surface.comps[i1].Dw = new_Dw;
	charge1 = temp_surface.comps[i1].charge;
	break;
      }
    }
    /* sum charge structures if not yet done... */
    if (found && temp_surface.type != NO_EDL && !charge_in)
    {
      temp_surface.charge[charge1].grams +=
	surface_ptr2->charge[charge2].grams * f2;
      temp_surface.charge[charge1].charge_balance +=
	surface_ptr2->charge[charge2].charge_balance * f2;
      temp_surface.charge[charge1].mass_water +=
	surface_ptr2->charge[charge2].mass_water * f2;
      count_elts = 0;
      add_elt_list (temp_surface.charge[charge1].diffuse_layer_totals, 1.0);
      add_elt_list (surface_ptr2->charge[charge2].diffuse_layer_totals, f2);
      if (count_elts > 0)
      {
	qsort (elt_list, (size_t) count_elts,
	       (size_t) sizeof (struct elt_list), elt_list_compare);
	elt_list_combine ();
      }
      temp_surface.charge[charge1].diffuse_layer_totals =
	(struct elt_list *) free_check_null (temp_surface.charge[charge1].
					     diffuse_layer_totals);
      temp_surface.charge[charge1].diffuse_layer_totals = elt_list_save ();
      charge_in = TRUE;
    }

    if (found == FALSE)
    {

      /* 2nd surface's component not in 1st surface */
      temp_surface.comps =
	(struct surface_comp *) PHRQ_realloc (temp_surface.comps,
					      (size_t) (count_comps +
							1) *
					      sizeof (struct surface_comp));
      if (temp_surface.comps == NULL)
	malloc_error ();
      memcpy (&temp_surface.comps[count_comps], &surface_ptr2->comps[i2],
	      sizeof (struct surface_comp));
      temp_surface.comps[count_comps].moles *= f2;
      temp_surface.comps[count_comps].cb *= f2;
      count_elts = 0;
      add_elt_list (surface_ptr2->comps[i2].totals, f2);
      temp_surface.comps[count_comps].totals = elt_list_save ();
      temp_surface.comps[count_comps].formula_totals = elt_list_save ();
      temp_surface.comps[i1].Dw = new_Dw;
      count_comps++;

      /* sum charge structures if not yet done... */
      if (temp_surface.type != NO_EDL && !charge_in)
      {

	temp_surface.charge =
	  (struct surface_charge *) PHRQ_realloc (temp_surface.charge,
						  (size_t) (count_charge +
							    1) *
						  sizeof (struct
							  surface_charge));
	if (temp_surface.charge == NULL)
	  malloc_error ();
	memcpy (&temp_surface.charge[count_charge],
		&surface_ptr2->charge[charge2],
		sizeof (struct surface_charge));
	temp_surface.charge[count_charge].grams *= f2;
	temp_surface.charge[count_charge].charge_balance *= f2;
	temp_surface.charge[count_charge].mass_water *= f2;
	temp_surface.charge[count_charge].g = NULL;
	temp_surface.charge[count_charge].count_g = 0;
	count_elts = 0;
	add_elt_list (surface_ptr2->charge[surface_ptr2->comps[k].charge].
		      diffuse_layer_totals, f2);
	temp_surface.charge[count_charge].diffuse_layer_totals =
	  elt_list_save ();
	count_charge++;
	charge_in = TRUE;
      }
      temp_surface.comps[count_comps - 1].charge = count_charge - 1;
    }
  }

  temp_surface.count_comps = count_comps;
  temp_surface.count_charge = count_charge;
  temp_surface.transport = FALSE;
  for (i = 0; i < count_comps; i++)
  {
    if (temp_surface.comps[i].Dw > 0)
    {
      temp_surface.transport = TRUE;
      break;
    }
  }

/*
 *   Finish up
 */
  surface_ptr = target;
  if (surface_ptr == NULL)
  {
    sprintf (error_string, "Target surface is NULL in sum_surface_comp");
    error_msg (error_string, CONTINUE);
    input_error++;
    return (ERROR);
  }
  surface_free (surface_ptr);
  surface_copy (&temp_surface, surface_ptr, new_n_user);
  surface_free (&temp_surface);
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
check_surfaces (struct surface *surface_ptr1, struct surface *surface_ptr2)
/* ---------------------------------------------------------------------- */
{
/*  checks if surfaces can be mixed...
 */
  int n_user1, n_user2, return_code;

  return_code = OK;
  n_user1 = surface_ptr1->n_user;
  n_user2 = surface_ptr2->n_user;

  /*if (surface_ptr1->diffuse_layer != surface_ptr2->diffuse_layer) { */
  if (surface_ptr1->dl_type != surface_ptr2->dl_type)
  {
    sprintf (error_string,
	     "Surfaces %d and %d differ in definition of diffuse layer. Can not mix.",
	     n_user1, n_user2);
    error_msg (error_string, CONTINUE);
    return_code = ERROR;
    input_error++;
  }
  /*if (surface_ptr1->edl != surface_ptr2->edl) { */
  if (surface_ptr1->type != surface_ptr2->type)
  {
    sprintf (error_string,
	     "Surfaces %d and %d differ in use of electrical double layer. Can not mix.",
	     n_user1, n_user2);
    error_msg (error_string, CONTINUE);
    return_code = ERROR;
    input_error++;
  }
  if (surface_ptr1->only_counter_ions != surface_ptr2->only_counter_ions)
  {
    sprintf (error_string,
	     "Surfaces %d and %d differ in use of only counter ions in the diffuse layer. Can not mix.",
	     n_user1, n_user2);
    error_msg (error_string, CONTINUE);
    return_code = ERROR;
    input_error++;
  }
  if (surface_ptr1->related_phases != surface_ptr2->related_phases)
  {
    sprintf (error_string,
	     "Surfaces %d and %d differ in use of related phases (sites proportional to moles of an equilibrium phase). Can not mix.",
	     n_user1, n_user2);
    error_msg (error_string, CONTINUE);
    return_code = ERROR;
    input_error++;
  }
  if (surface_ptr1->related_rate != surface_ptr2->related_rate)
  {
    sprintf (error_string,
	     "Surfaces %d and %d differ in use of related rate (sites proportional to moles of a kinetic reactant). Can not mix.",
	     n_user1, n_user2);
    error_msg (error_string, CONTINUE);
    return_code = ERROR;
    input_error++;
  }

  return (return_code);
}

/* ---------------------------------------------------------------------- */
int
mobile_surface_copy (struct surface *surface_old_ptr,
		     struct surface *surf_ptr1, int n_user_new, int move_old)
/* ---------------------------------------------------------------------- */
{
/*
 *   Copies mobile comps from surface_old_ptr to surf_ptr1,
 *   comps and charges with Dw > 0 are moved if move_old == TRUE, else copied.
 *   NOTE... when all comps are moved, the old surface is deleted and surfaces are sorted again,
 *                 which will modify pointers and surface numbers.
 *   User number of new surface structure is n_user_new, structure is freed when n_user_new is already defined
 */
  int count_comps, count_charge, new_charge, n_user_old, charge_done;
  int i, i1, i2, j, k, k1, n;
  char token[MAX_LENGTH];
  struct surface temp_surface, *surf_ptr;
/*
 *   Store moving surface's properties in temp_surface
 */
  temp_surface.n_user = temp_surface.n_user_end = n_user_new;
  temp_surface.new_def = surface_old_ptr->new_def;
  sprintf (token, "Surface defined in simulation %d.", simulation);
  temp_surface.description = string_duplicate (token);
  temp_surface.dl_type = surface_old_ptr->dl_type;
  temp_surface.type = surface_old_ptr->type;
  temp_surface.only_counter_ions = surface_old_ptr->only_counter_ions;
  temp_surface.thickness = surface_old_ptr->thickness;
  temp_surface.debye_lengths = surface_old_ptr->debye_lengths;
  temp_surface.DDL_viscosity = surface_old_ptr->DDL_viscosity;
  temp_surface.DDL_limit = surface_old_ptr->DDL_limit;
  temp_surface.solution_equilibria = FALSE;
  temp_surface.n_solution = -99;
  temp_surface.related_phases = surface_old_ptr->related_phases;
  temp_surface.related_rate = surface_old_ptr->related_rate;
  temp_surface.transport = TRUE;

  temp_surface.comps =
    (struct surface_comp *) PHRQ_malloc (sizeof (struct surface_comp));
  if (temp_surface.comps == NULL)
    malloc_error ();
  temp_surface.charge =
    (struct surface_charge *) PHRQ_malloc (sizeof (struct surface_charge));
  if (temp_surface.charge == NULL)
    malloc_error ();

  count_comps = surface_old_ptr->count_comps;
  count_charge = surface_old_ptr->count_charge;
  i1 = i2 = 0;
  /* see if comps must be moved, Dw > 0 */
  for (i = 0; i < count_comps; i++)
  {
    if (surface_old_ptr->comps[i].Dw > 0)
    {
      i1++;
      if (i1 > 1)
      {
	temp_surface.comps =
	  (struct surface_comp *) PHRQ_realloc (temp_surface.comps,
						(size_t) i1 *
						sizeof (struct surface_comp));
	if (temp_surface.comps == NULL)
	  malloc_error ();
      }
      memcpy (&temp_surface.comps[i1 - 1], &surface_old_ptr->comps[i],
	      sizeof (struct surface_comp));
      temp_surface.comps[i1 - 1].totals =
	elt_list_dup (surface_old_ptr->comps[i].totals);
      temp_surface.comps[i1 - 1].formula_totals =
	elt_list_dup (surface_old_ptr->comps[i].formula_totals);

      /* charge structure */
      new_charge = TRUE;
      for (j = 0; j < i; j++)
      {
	if (surface_old_ptr->comps[j].charge ==
	    surface_old_ptr->comps[i].charge)
	{
	  temp_surface.comps[i1 - 1].charge = i2 - 1;
	  new_charge = FALSE;
	  break;
	}
      }
      if (new_charge)
      {
	i2++;
	if (i2 > 1)
	{
	  temp_surface.charge =
	    (struct surface_charge *) PHRQ_realloc (temp_surface.charge,
						    (size_t) i2 *
						    sizeof (struct
							    surface_charge));
	  if (temp_surface.charge == NULL)
	    malloc_error ();
	}
	memcpy (&temp_surface.charge[i2 - 1],
		&surface_old_ptr->charge[surface_old_ptr->comps[i].charge],
		(size_t) sizeof (struct surface_charge));
	temp_surface.charge[i2 - 1].diffuse_layer_totals =
	  elt_list_dup (surface_old_ptr->
			charge[surface_old_ptr->comps[i].charge].
			diffuse_layer_totals);
	temp_surface.charge[i2 - 1].g = NULL;

	temp_surface.comps[i1 - 1].charge = i2 - 1;
      }
    }
  }
  temp_surface.count_comps = i1;
  temp_surface.count_charge = i2;
  if (i1 > 0)
  {
    /* OK, store moved parts from old surface, but first:
       get immobile surface comps from new surface... */
    if ((surf_ptr = surface_bsearch (n_user_new, &n)) != NULL)
    {
      for (k = 0; k < surf_ptr->count_comps; k++)
      {
	if (surf_ptr->comps[k].Dw > 0)
	  continue;
	charge_done = FALSE;
	for (k1 = 0; k1 < k; k1++)
	{
	  if (surf_ptr->comps[k1].charge == surf_ptr->comps[k].charge)
	  {
	    charge_done = TRUE;
	    break;
	  }
	}
	if (charge_done)
	  continue;

	if (sum_surface_comp
	    (&temp_surface, 1, surf_ptr, k, 1, &temp_surface, 0) == ERROR)
	  return (ERROR);
      }
    }
    if (surf_ptr1->count_comps > 0)
      surface_free (surf_ptr1);
    surface_copy (&temp_surface, surf_ptr1, n_user_new);
  }

  /* delete moved parts from old surface */
  if (move_old && i1 > 0)
  {
    n_user_old = surface_old_ptr->n_user;
    if (i1 != count_comps)
    {
      /* redefine old surface with only unmovable comps (Dw = 0) */
      /* other temp_surface members were set above */
      temp_surface.n_user = temp_surface.n_user_end = n_user_old;
      temp_surface.transport = FALSE;

      for (i = 0; i < temp_surface.count_comps; i++)
      {
	temp_surface.comps[i].totals =
	  (struct elt_list *) free_check_null (temp_surface.comps[i].totals);
	temp_surface.comps[i].formula_totals =
	  (struct elt_list *) free_check_null (temp_surface.comps[i].
					       formula_totals);
      }
      for (i = 0; i < temp_surface.count_charge; i++)
      {
	temp_surface.charge[i].diffuse_layer_totals =
	  (struct elt_list *) free_check_null (temp_surface.charge[i].
					       diffuse_layer_totals);
	temp_surface.charge[i].g =
	  (struct surface_diff_layer *) free_check_null (temp_surface.
							 charge[i].g);
      }
      temp_surface.comps =
	(struct surface_comp *) PHRQ_realloc (temp_surface.comps,
					      sizeof (struct surface_comp));
      if (temp_surface.comps == NULL)
	malloc_error ();
      temp_surface.charge =
	(struct surface_charge *) PHRQ_realloc (temp_surface.charge,
						sizeof (struct
							surface_charge));
      if (temp_surface.charge == NULL)
	malloc_error ();

      i1 = i2 = 0;
      for (i = 0; i < count_comps; i++)
      {
	if (surface_old_ptr->comps[i].Dw == 0)
	{
	  i1++;
	  if (i1 > 1)
	  {
	    temp_surface.comps =
	      (struct surface_comp *) PHRQ_realloc (temp_surface.comps,
						    (size_t) i1 *
						    sizeof (struct
							    surface_comp));
	    if (temp_surface.comps == NULL)
	      malloc_error ();
	  }
	  memcpy (&temp_surface.comps[i1 - 1], &surface_old_ptr->comps[i],
		  (size_t) sizeof (struct surface_comp));
	  temp_surface.comps[i1 - 1].totals =
	    elt_list_dup (surface_old_ptr->comps[i].totals);
	  temp_surface.comps[i1 - 1].formula_totals =
	    elt_list_dup (surface_old_ptr->comps[i].formula_totals);

	  /* charge structure */
	  new_charge = TRUE;
	  for (j = 0; j < i; j++)
	  {
	    if (surface_old_ptr->comps[j].charge ==
		surface_old_ptr->comps[i].charge)
	    {
	      temp_surface.comps[i1 - 1].charge = i2 - 1;
	      new_charge = FALSE;
	      break;
	    }
	  }
	  if (new_charge)
	  {
	    i2++;
	    if (i2 > 1)
	    {
	      temp_surface.charge =
		(struct surface_charge *) PHRQ_realloc (temp_surface.charge,
							(size_t) i2 *
							sizeof (struct
								surface_charge));
	      if (temp_surface.charge == NULL)
		malloc_error ();
	    }
	    memcpy (&temp_surface.charge[i2 - 1],
		    &surface_old_ptr->charge[surface_old_ptr->comps[i].
					     charge],
		    (size_t) sizeof (struct surface_charge));
	    temp_surface.charge[i2 - 1].diffuse_layer_totals =
	      elt_list_dup (surface_old_ptr->
			    charge[surface_old_ptr->comps[i].charge].
			    diffuse_layer_totals);
	    temp_surface.charge[i2 - 1].g = NULL;

	    temp_surface.comps[i1 - 1].charge = i2 - 1;
	  }
	}
      }
      temp_surface.count_comps = i1;
      temp_surface.count_charge = i2;
      if (i1 > 0)
      {
	temp_surface.related_phases = temp_surface.related_rate = FALSE;
	for (i = 0; i < i1; i++)
	{
	  if (surface_old_ptr->comps[i].phase_name != NULL)
	    temp_surface.related_phases = TRUE;
	  if (surface_old_ptr->comps[i].rate_name != NULL)
	    temp_surface.related_rate = TRUE;
	}
      }
      surface_free (surface_old_ptr);
      surface_copy (&temp_surface, surface_old_ptr, n_user_old);
    }
    else
    {
      surface_sort ();
      surface_delete (n_user_old);
    }
  }

  surface_free (&temp_surface);
  return (OK);
}

/* ---------------------------------------------------------------------- */
int
diff_stag_surf (int mobile_cell)
/* ---------------------------------------------------------------------- */
/*
 *  Diffuse stagnant and mobile surfaces, following the steps of disp_surf.
 *  First the mobile/stagnant surfaces are mixed, then the stagnant surfaces
 *  when not already done.
 *  If mixing factors among the cells are defined expicitly, it is assumed that
 *  mixing with a lower numbered cell was done when that cell was processed:
 *  for any cell in MCD, need only include the mixing factors for higher numbered cells.
 */
{
  int i, i1, i2, j, k, k1, n, n1, n2, ns, temp_count_surface;
  int charge_done, surf1, surf2;
  LDBLE f1, f2, mixf, mixf_store;
  LDBLE Dp1, Dp2;
  struct mix *mix_ptr;
  struct surface *surface_ptr1, *surface_ptr2;
  LDBLE viscos_f;
/*
 * temperature and viscosity correction for MCD coefficient, D_T = D_298 * Tk * viscos_298 / (298 * viscos)
 */
  viscos_f =
    pow (10,
	 -(1.37023 * (tc_x - 20) +
	   0.000836 * (tc_x - 20) * (tc_x - 20)) / (109 + tc_x));
  viscos_f = tk_x * 0.88862 / (298.15 * viscos_f);

  n1 = count_surface;
  n2 = n1 + 1;
  temp_count_surface = count_surface + 2;
  space ((void **) ((void *) &surface), temp_count_surface, &max_surface,
	 sizeof (struct surface));
  surface[n1].count_comps = surface[n1 + 1].count_comps = 0;

  for (ns = 0; ns < stag_data->count_stag; ns++)
  {

    i1 = mobile_cell + 1 + ns * count_cells;
    if (ns == 0)
      i1--;

    if (cell_data[i1 - 1].por < multi_Dpor_lim)
      continue;
    surface[n1].n_user = surface[n1 + 1].n_user = -99;
    surface_ptr1 = surface_bsearch (i1, &n);

/*
 * step 2a. mix surfaces...
 */
    mix_ptr = mix_search (i1, &n, FALSE);
    if (mix_ptr == NULL)
      continue;

    for (j = 0; j < mix_ptr->count_comps; j++)
    {

      if ((i2 = mix_ptr->comps[j].n_solution) <= i1)
	continue;
      if (cell_data[i2 - 1].por < multi_Dpor_lim)
	continue;
      surface_ptr2 = surface_bsearch (i2, &n);

      surf1 = surf2 = 0;
      if (surface_ptr2 != NULL)
      {
	if (surface_ptr2->transport == TRUE)
	  surf2 = 1;
      }
      if (surface_ptr1 != NULL)
      {
	if (surface_ptr1->transport == TRUE)
	  surf1 = 1;
      }
      if (!surf1 && !surf2)
	continue;

      mixf = mixf_store = mix_ptr->comps[j].fraction;

      /* find the (possibly modified) surface in cell i1... */
      f1 = 1;
      for (k = count_surface + 2; k < temp_count_surface; k++)
      {
	if (surface[k].n_user == i1)
	{
	  n1 = k;
	  break;
	}
      }
      /* if not found... */
      if (k == temp_count_surface)
      {
	n1 = count_surface;
	/* copy it from surface_ptr1... */
	if (surface_ptr1 != NULL)
	  surface_copy (surface_ptr1, &surface[n1], i1);
	else
	{
	  /* or make it a mobile copy of the surface in cell i2... */
	  mobile_surface_copy (surface_ptr2, &surface[n1], i1, FALSE);
	  /* limit comps to 1... */
	  if (surface[n1].count_comps > 1)
	  {
	    if (sum_surface_comp
		(&surface[n1], 0, &surface[n1], 0, 1, &surface[n1],
		 surface[n1].comps[0].Dw) == ERROR)
	      return (ERROR);
	  }
	  f1 = 0;
	}
      }
      /* find the (possibly modified) surface in cell i2... */
      f2 = 1;
      for (k = count_surface + 2; k < temp_count_surface; k++)
      {
	if (surface[k].n_user == i2)
	{
	  n2 = k;
	  break;
	}
      }
      /* if not found... */
      if (k == temp_count_surface)
      {
	n2 = count_surface + 1;
	/* copy it from surface_ptr2... */
	if (surface_ptr2 != NULL)
	  surface_copy (surface_ptr2, &surface[n2], i2);
	else
	{
	  /* or make it a mobile copy of the surface in cell i1... */
	  mobile_surface_copy (surface_ptr1, &surface[n2], i2, FALSE);
	  /* limit comps to 1... */
	  if (surface[n2].count_comps > 1)
	  {
	    if (sum_surface_comp
		(&surface[n2], 0, &surface[n2], 0, 1, &surface[n2],
		 surface[n2].comps[0].Dw) == ERROR)
	      return (ERROR);
	  }
	  f2 = 0;
	}
      }

      /* For MCD, step 2b. Adapt mixf to default values... */
      if (multi_Dflag)
      {
	mixf_store *=
	  (cell_data[i1 - 1].por <=
	   cell_data[i2 - 1].por ? cell_data[i1 - 1].por : cell_data[i2 -
								     1].por);
	mixf_store /= (default_Dw * pow (multi_Dpor, multi_Dn) * multi_Dpor);
      }

      /* mix in comps with the same charge structure... */
      if (surf2)
      {
	for (k = 0; k < surface_ptr2->count_comps; k++)
	{
	  if (surface_ptr2->comps[k].Dw == 0)
	    continue;
	  charge_done = FALSE;
	  for (k1 = 0; k1 < k; k1++)
	  {
	    if (surface_ptr2->comps[k1].charge ==
		surface_ptr2->comps[k].charge)
	    {
	      charge_done = TRUE;
	      break;
	    }
	  }
	  if (charge_done)
	    continue;

	  /* find diffusion coefficients of surfaces... */
	  if (multi_Dflag)
	  {
	    Dp2 =
	      surface_ptr2->comps[k].Dw * pow (cell_data[i2 - 1].por,
					       multi_Dn) * viscos_f;
	    Dp1 = 0;
	    if (surf1)
	    {
	      for (k1 = 0; k1 < surface_ptr1->count_comps; k1++)
	      {
		if (elt_list_compare
		    (surface_ptr1->comps[k1].totals,
		     surface_ptr2->comps[k].totals) != 0)
		  continue;
		Dp1 =
		  surface_ptr1->comps[k1].Dw * pow (cell_data[i1 - 1].por,
						    multi_Dn) * viscos_f;
		break;
	      }
	    }
	    if (Dp1 > 0)
	      Dp2 = (Dp2 + Dp1) / 2;

	    /* and adapt the mixing factor... */
	    mixf = mixf_store * Dp2;
	    mixf /= solution_bsearch (i2, &n, FALSE)->mass_water;
	  }

	  if (mixf < 1e-8)
	    mixf = 0;
	  if (mixf > 0.99999999)
	    mixf = 0.99999999;
	  if (sum_surface_comp
	      (&surface[n1], f1, surface_ptr2, k, mixf, &surface[n1],
	       surface_ptr2->comps[k].Dw) == ERROR)
	    return (ERROR);
	  f1 = 1;

	  if (sum_surface_comp
	      (&surface[n2], f2, surface_ptr2, k, -mixf, &surface[n2],
	       surface_ptr2->comps[k].Dw) == ERROR)
	    return (ERROR);
	}
      }

      if (surf1)
      {
	for (k = 0; k < surface_ptr1->count_comps; k++)
	{
	  if (surface_ptr1->comps[k].Dw == 0)
	    continue;
	  charge_done = FALSE;
	  for (k1 = 0; k1 < k; k1++)
	  {
	    if (surface_ptr1->comps[k1].charge ==
		surface_ptr1->comps[k].charge)
	    {
	      charge_done = TRUE;
	      break;
	    }
	  }
	  if (charge_done)
	    continue;

	  /* find diffusion coefficients of surfaces... */
	  if (multi_Dflag)
	  {
	    Dp1 =
	      surface_ptr1->comps[k].Dw * pow (cell_data[i1 - 1].por,
					       multi_Dn) * viscos_f;

	    Dp2 = 0;
	    if (surf2)
	    {
	      for (k1 = 0; k1 < surface_ptr2->count_comps; k1++)
	      {
		if (elt_list_compare
		    (surface_ptr2->comps[k1].totals,
		     surface_ptr1->comps[k].totals) != 0)
		  continue;
		Dp2 =
		  surface_ptr2->comps[k1].Dw * pow (cell_data[i2 - 1].por,
						    multi_Dn) * viscos_f;
		break;
	      }
	    }
	    if (Dp2 > 0)
	      Dp1 = (Dp1 + Dp2) / 2;

	    /* and adapt the mixing factor... */
	    mixf = mixf_store * Dp1;
	    mixf /= solution_bsearch (i1, &n, FALSE)->mass_water;
	  }

	  if (mixf < 1e-8)
	    mixf = 0;
	  if (mixf > 0.99999999)
	    mixf = 0.99999999;
	  if (sum_surface_comp
	      (&surface[n2], f2, surface_ptr1, k, mixf, &surface[n2],
	       surface_ptr1->comps[k].Dw) == ERROR)
	    return (ERROR);
	  f2 = 1;

	  if (sum_surface_comp
	      (&surface[n1], f1, surface_ptr1, k, -mixf, &surface[n1],
	       surface_ptr1->comps[k].Dw) == ERROR)
	    return (ERROR);
	}
      }

/*
 *  Step 3. copy surface[n1] and [n2] in a new temporary surface...
 */
      if (surface[n1].n_user == -99)
	continue;

      if (n1 == count_surface)
      {
	n = temp_count_surface++;

	space ((void **) ((void *) &surface), temp_count_surface,
	       &max_surface, sizeof (struct surface));
	surface_copy (&surface[n1], &surface[n], i1);

	surface_free (&surface[n1]);
	surface[n1].count_comps = 0;
      }
      else
	n1 = count_surface;

      if (n2 == count_surface + 1)
      {
	n = temp_count_surface++;

	space ((void **) ((void *) &surface), temp_count_surface,
	       &max_surface, sizeof (struct surface));
	surface_copy (&surface[n2], &surface[n], i2);

	surface_free (&surface[n2]);
	surface[n2].count_comps = 0;
      }
      else
	n2 = count_surface + 1;
    }
  }
/*
 * Step 4. Diffusion is done. New surfaces can be copied in the cells...
 */

  n2 = 0;

  for (j = count_surface + 2; j < temp_count_surface; j++)
  {
    i = surface[j].n_user;
    if ((i == 0 && bcon_first == 1)
	|| (i == count_cells + 1 && bcon_last == 1))
    {
      surface_free (&surface[j]);
      surface[j].count_comps = 0;
      continue;
    }
    if (i >= 0 && i <= 1 + count_cells * (1 + stag_data->count_stag))
    {
      surface_ptr1 = surface_bsearch (i, &n);
      if (surface_ptr1 != NULL)
      {
	surface_free (surface_ptr1);
	surface_copy (&surface[j], surface_ptr1, i);
      }
      else
      {
	surface_copy (&surface[j], &surface[count_surface + n2], i);
	n2++;
      }
    }
    surface_free (&surface[j]);
    surface[j].count_comps = 0;
  }
  count_surface += n2;
  if (n2 > 0)
    qsort (surface, (size_t) count_surface, (size_t) sizeof (struct surface),
	   surface_compare);

  return (OK);
}

/* ---------------------------------------------------------------------- */
int
reformat_surf (char *comp_name, LDBLE fraction, char *new_comp_name,
	       LDBLE new_Dw, int cell)
/* ---------------------------------------------------------------------- */
{
  int i, i1, j, k, k1, n, length, length1, cn, cn1, charge_in;
  char string[MAX_LENGTH];
  struct surface temp_surface, *surf_ptr;

  if ((surf_ptr = surface_bsearch (cell, &n)) == ERROR)
    return (OK);

  if (fraction > 0.99999999)
    fraction = 0.99999999;
  length = (int) strlen (comp_name);
  length1 = (int) strlen (new_comp_name);
  for (k = 0; k < surf_ptr->count_comps; k++)
  {
    if (strncmp (comp_name, surf_ptr->comps[k].formula, length) == 0)
      break;
  }
  if (k == surf_ptr->count_comps)
    return (OK);

  surface_copy (surf_ptr, &temp_surface, cell);

  for (k1 = 0; k1 < temp_surface.count_comps; k1++)
  {
    if (strncmp (new_comp_name, temp_surface.comps[k1].formula, length1) == 0)
      break;
  }

  charge_in = FALSE;
  if (k1 == temp_surface.count_comps)
  {
    /* new_comp_name absent */
    /* rename the comps... */
    for (i = 0; i < temp_surface.count_comps; i++)
    {
      if (temp_surface.comps[i].charge != temp_surface.comps[k].charge)
	continue;

      strcpy (string, temp_surface.comps[i].formula);
      replace (comp_name, new_comp_name, string);
      temp_surface.comps[i].formula = string_hsave (string);

      temp_surface.comps[i].moles *= fraction;

      count_elts = 0;
      add_elt_list (temp_surface.comps[i].totals, fraction);
      temp_surface.comps[i].totals =
	(struct elt_list *) free_check_null (temp_surface.comps[i].totals);
      temp_surface.comps[i].totals = elt_list_save ();
/* appt, is this correct ? */
      temp_surface.comps[i].formula_totals =
	(struct elt_list *) free_check_null (temp_surface.comps[i].
					     formula_totals);
      temp_surface.comps[i].formula_totals = elt_list_save ();

      /* rename surface comp in element list */
      for (j = 0; temp_surface.comps[i].totals[j].elt != NULL; j++)
      {
	if (strncmp
	    (comp_name, temp_surface.comps[i].totals[j].elt->name,
	     length) == 0)
	{
	  strcpy (string, temp_surface.comps[i].totals[j].elt->name);
	  replace (comp_name, new_comp_name, string);
	  temp_surface.comps[i].totals[j].elt = element_store (string);
	}
      }
      for (j = 0; temp_surface.comps[i].formula_totals[j].elt != NULL; j++)
      {
	if (strncmp
	    (comp_name, temp_surface.comps[i].formula_totals[j].elt->name,
	     length) == 0)
	{
	  strcpy (string, temp_surface.comps[i].formula_totals[j].elt->name);
	  replace (comp_name, new_comp_name, string);
	  temp_surface.comps[i].formula_totals[j].elt =
	    element_store (string);
	}
      }
      /* add charge */
      if (!charge_in)
      {
	cn = temp_surface.comps[i].charge;

	temp_surface.charge[cn].grams *= fraction;
	temp_surface.charge[cn].charge_balance *= fraction;
	temp_surface.charge[cn].mass_water *= fraction;

	count_elts = 0;
	add_elt_list (temp_surface.charge[cn].diffuse_layer_totals, fraction);
	temp_surface.charge[cn].diffuse_layer_totals =
	  (struct elt_list *) free_check_null (temp_surface.charge[cn].
					       diffuse_layer_totals);
	temp_surface.charge[cn].diffuse_layer_totals = elt_list_save ();

	strcpy (string, temp_surface.charge[cn].name);
	replace (comp_name, new_comp_name, string);
	temp_surface.charge[cn].name = string_hsave (string);
	/*  renaming charge psi_master seems unnecessary... */
#ifdef SKIP
	strcpy (string, temp_surface.charge[cn].psi_master[0].elt->name);
	replace (comp_name, new_comp_name, string);
	temp_surface.charge[cn].psi_master[0].elt->name =
	  string_hsave (string);
#endif

	charge_in = TRUE;
      }
      temp_surface.comps[i].Dw = new_Dw;
    }

    /* add remainder of old comp and charge */
    if (fraction < 1)
      sum_surface_comp (&temp_surface, 1, surf_ptr, k, 1 - fraction,
			&temp_surface, surf_ptr->comps[k].Dw);
  }
  else
  {
    /* new_comp_name is already in surface */
    cn = surf_ptr->comps[k].charge;
    for (i = 0; i < surf_ptr->count_comps; i++)
    {

      if (surf_ptr->comps[i].charge != cn)
	continue;

      strcpy (string, surf_ptr->comps[i].formula);
      replace (comp_name, new_comp_name, string);

      /* add comps to new_comps */
      for (i1 = 0; i1 < temp_surface.count_comps; i1++)
      {

	if (strcmp (string, temp_surface.comps[i1].formula) != 0)
	  continue;

	temp_surface.comps[i1].moles += surf_ptr->comps[i].moles * fraction;
	temp_surface.comps[i].moles -= surf_ptr->comps[i].moles * fraction;

	/* rename surface comp in element list */
	for (j = 0; temp_surface.comps[i].totals[j].elt != NULL; j++)
	{
	  if (strncmp
	      (comp_name, temp_surface.comps[i].totals[j].elt->name,
	       length) == 0)
	  {
	    strcpy (string, temp_surface.comps[i].totals[j].elt->name);
	    replace (comp_name, new_comp_name, string);
	    temp_surface.comps[i].totals[j].elt = element_store (string);
	  }
	}
	for (j = 0; temp_surface.comps[i].formula_totals[j].elt != NULL; j++)
	{
	  if (strncmp
	      (comp_name, temp_surface.comps[i].formula_totals[j].elt->name,
	       length) == 0)
	  {
	    strcpy (string,
		    temp_surface.comps[i].formula_totals[j].elt->name);
	    replace (comp_name, new_comp_name, string);
	    temp_surface.comps[i].formula_totals[j].elt =
	      element_store (string);
	  }
	}
	count_elts = 0;
	add_elt_list (temp_surface.comps[i1].totals, 1.0);
	add_elt_list (temp_surface.comps[i].totals, fraction);
	if (count_elts > 0)
	{
	  qsort (elt_list, (size_t) count_elts,
		 (size_t) sizeof (struct elt_list), elt_list_compare);
	  elt_list_combine ();
	}
	temp_surface.comps[i1].totals =
	  (struct elt_list *) free_check_null (temp_surface.comps[i1].totals);
	temp_surface.comps[i1].totals = elt_list_save ();
	temp_surface.comps[i1].formula_totals =
	  (struct elt_list *) free_check_null (temp_surface.comps[i1].
					       formula_totals);
	temp_surface.comps[i1].formula_totals = elt_list_save ();
	temp_surface.comps[i1].Dw = new_Dw;

	count_elts = 0;
	add_elt_list (surf_ptr->comps[i].totals, 1.0 - fraction);
	temp_surface.comps[i].totals =
	  (struct elt_list *) free_check_null (temp_surface.comps[i].totals);
	temp_surface.comps[i].totals = elt_list_save ();
	temp_surface.comps[i].formula_totals =
	  (struct elt_list *) free_check_null (temp_surface.comps[i].
					       formula_totals);
	temp_surface.comps[i].formula_totals = elt_list_save ();

	/* add charge */
	if (!charge_in)
	{
	  cn1 = temp_surface.comps[i1].charge;
	  temp_surface.charge[cn1].grams +=
	    surf_ptr->charge[cn].grams * fraction;
	  temp_surface.charge[cn1].charge_balance +=
	    surf_ptr->charge[cn].charge_balance * fraction;
	  temp_surface.charge[cn1].mass_water +=
	    surf_ptr->charge[cn].mass_water * fraction;

	  temp_surface.charge[cn].grams -=
	    surf_ptr->charge[cn].grams * fraction;
	  temp_surface.charge[cn].charge_balance -=
	    surf_ptr->charge[cn].charge_balance * fraction;
	  temp_surface.charge[cn].mass_water -=
	    surf_ptr->charge[cn].mass_water * fraction;

	  count_elts = 0;
	  add_elt_list (temp_surface.charge[cn1].diffuse_layer_totals, 1.0);
	  add_elt_list (surf_ptr->charge[cn].diffuse_layer_totals, fraction);
	  if (count_elts > 0)
	  {
	    qsort (elt_list, (size_t) count_elts,
		   (size_t) sizeof (struct elt_list), elt_list_compare);
	    elt_list_combine ();
	  }
	  temp_surface.charge[cn1].diffuse_layer_totals =
	    (struct elt_list *) free_check_null (temp_surface.charge[cn1].
						 diffuse_layer_totals);
	  temp_surface.charge[cn1].diffuse_layer_totals = elt_list_save ();

	  count_elts = 0;
	  add_elt_list (surf_ptr->charge[cn].diffuse_layer_totals,
			1.0 - fraction);
	  temp_surface.charge[cn].diffuse_layer_totals =
	    (struct elt_list *) free_check_null (temp_surface.charge[cn].
						 diffuse_layer_totals);
	  temp_surface.charge[cn].diffuse_layer_totals = elt_list_save ();

	  charge_in = TRUE;
	}
      }
    }
  }

  temp_surface.transport = FALSE;
  for (i = 0; i < temp_surface.count_comps; i++)
  {
    if (temp_surface.comps[i].Dw > 0)
      temp_surface.transport = TRUE;
  }
  surface_free (surf_ptr);
  surface_copy (&temp_surface, surf_ptr, cell);
  surface_free (&temp_surface);

  return (OK);
}
