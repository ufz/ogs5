#ifdef PHREEQC_IDENT
static char const svnidkinetics[] =
  "$Id: kinetics.h 4 2009-04-21 17:29:29Z delucia $";
#endif
KINETICS_EXTERNAL void *cvode_kinetics_ptr;
KINETICS_EXTERNAL int cvode_test;
KINETICS_EXTERNAL int cvode_error;
KINETICS_EXTERNAL int cvode_n_user;
KINETICS_EXTERNAL int cvode_n_reactions;
KINETICS_EXTERNAL realtype cvode_step_fraction;
KINETICS_EXTERNAL realtype cvode_rate_sim_time;
KINETICS_EXTERNAL realtype cvode_rate_sim_time_start;
KINETICS_EXTERNAL realtype cvode_last_good_time;
KINETICS_EXTERNAL realtype cvode_prev_good_time;
KINETICS_EXTERNAL N_Vector cvode_last_good_y;
KINETICS_EXTERNAL N_Vector cvode_prev_good_y;
KINETICS_EXTERNAL M_Env kinetics_machEnv;
KINETICS_EXTERNAL N_Vector kinetics_y, kinetics_abstol;
KINETICS_EXTERNAL void *kinetics_cvode_mem;
KINETICS_EXTERNAL struct pp_assemblage *cvode_pp_assemblage_save;
KINETICS_EXTERNAL struct s_s_assemblage *cvode_s_s_assemblage_save;
