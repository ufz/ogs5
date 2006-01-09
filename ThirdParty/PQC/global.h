
#ifdef PHREEQC_IDENT
static char const svnidglobal[] =
  "$Id: global.h 4 2009-04-21 17:29:29Z delucia $";
#endif
#ifndef _INC_GLOBAL_H
#define _INC_GLOBAL_H

#define NO_DOS
/* #define PHREEQ98 */ /* PHREEQ98: code for graphical user interface */
#ifdef PHREEQ98
#define isnan _isnan
#endif
/*
 * uncomment following line, to use default DOS file name for
 * output file
 */
/*#define DOS*/
/*
 *   BUG FIX FOR DGs
 */
#ifndef __OPEN_NAMESPACE__
#define __OPEN_NAMESPACE__
#endif
/* ----------------------------------------------------------------------
 *   INCLUDE FILES
 * ---------------------------------------------------------------------- */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <errno.h>
#include <float.h>
#include <assert.h>
#include <setjmp.h>
#include "phrqtype.h"

/* must be defined here and in cl.c */
/* #include <nan.h> */
#ifndef NAN
#   define NAN -99999999
#endif
#define MISSING -9999.999
/* search.h -- declarations for POSIX/SVID-compatible search functions */

/* HSEARCH(3C) */
typedef struct entry
{
  char *key;
  void *data;
} ENTRY;
typedef enum
{ FIND, ENTER } ACTION;

/* TSEARCH(3C) */
typedef enum
{ preorder, postorder, endorder, leaf } VISIT;

/* ----------------------------------------------------------------------
 *   DEFINITIONS
 * ---------------------------------------------------------------------- */
#define F_C_MOL 96493.5		/* C/mol or joule/volt-eq */
#define F_KJ_V_EQ  96.4935	/* kJ/volt-eq */
#define F_KCAL_V_EQ 23.0623	/* kcal/volt-eq */
#define R_LITER_ATM 0.0820597	/* L-atm/deg-mol */
#define R_KCAL_DEG_MOL 0.00198726	/* kcal/deg-mol */
#define R_KJ_DEG_MOL 0.00831470	/* kJ/deg-mol */
#define EPSILON 78.5		/* dialectric constant, dimensionless */
#define EPSILON_ZERO 8.854e-12	/* permittivity of free space, C/V-m = C**2/m-J */
#define JOULES_PER_CALORIE 4.1840
#define AVOGADRO 6.02252e23	/* atoms / mole */
typedef enum
{ kcal, cal, kjoules, joules } DELTA_H_UNIT;

#define TRUE 1
#define FALSE 0
#define OK 1
#define ERROR 0
#define STOP 1
#define CONTINUE 0

#define DISP 2
#define STAG 3
#define NOMIX 4

#define CONVERGED 2
#define MASS_BALANCE 3
/*
  #define OSCILLATE 4
  #define H2O_LIMITS 5
*/
#define REWRITE 2
#define INIT -1

/* check_line values, plus EMPTY, EOF, OK */
#define KEYWORD 3

/* copy_token values */
#define EMPTY 2
#define UPPER 4
#define LOWER 5
#define DIGIT 6
#define UNKNOWN 7
#define OPTION 8

/* species types */
#define AQ 0
#define HPLUS 1
#define H2O 2
#define EMINUS 3
#define SOLID 4
#define EX 5
#define SURF 6
#define SURF_PSI 7
#define SURF_PSI1 8
#define SURF_PSI2 9

/* unknown types */
#define MB 10
#define ALK 11
#define CB 12
#define SOLUTION_PHASE_BOUNDARY 13
#define MU 14
#define AH2O 15
#define MH 16
#define MH2O 17
#define PP 18
#define EXCH 19
#define SURFACE 20
#define SURFACE_CB 21
#define SURFACE_CB1 22
#define SURFACE_CB2 23
#define GAS_MOLES 24
#define S_S_MOLES 25
#define PITZER_GAMMA 26
/* state */
#define INITIALIZE               0
#define INITIAL_SOLUTION   1
#define INITIAL_EXCHANGE   2
#define INITIAL_SURFACE 3
#define INITIAL_GAS_PHASE  4
#define REACTION                   5
#define INVERSE                 6
#define ADVECTION                 7
#define TRANSPORT                 8
#define PHAST                     9

/* constaints in mass balance */
#define EITHER 0
#define DISSOLVE 1
#define PRECIPITATE -1

/* gas phase type */
#define PRESSURE 1
#define VOLUME 2

#define MAX_PP_ASSEMBLAGE 10	/* default estimate of the number of phase assemblages */
#define MAX_ADD_EQUATIONS 20	/* maximum number of equations added together to reduce eqn to
				   master species */
#define MAX_ELEMENTS 50		/* default estimate of the number of elements */
#define MAX_LENGTH 256		/* maximum number of characters component name */
#define MAX_LINE 400		/* estimate of maximum line length */
#define MAX_LM 3.0		/* maximum log molality allowed in intermediate iterations */
#define MIN_LM -30.0		/* minimum log molality allowed before molality set to zero */
#define MAX_MASS_BALANCE 10	/* initial guess of number mass balance equations for a solution */
#define MAX_MASTER 50		/* default estimate of the number of master species */
#define MAX_ELTS 15		/* default estimate for maximum number of times elements occur in
				   an equation */
#define MAX_PHASES 500		/* initial guess of number of phases defined */
#define MAX_SOLUTION 10		/* The maximum number of solutions allowed */
#define MAX_S 500		/* default estimate for maximum number of species in aqueous model */
#define MAX_STRINGS 3000
#define MAX_SUM_JACOB0 50	/* list used to calculate jacobian */
#define MAX_SUM_JACOB1 500	/* list used to calculate jacobian */
#define MAX_SUM_JACOB2 500	/* list used to calculate jacobian */
#define MAX_SUM_MB 500		/* list used to calculate mass balance sums */
#define MAX_TRXN 16		/* default estimate for maximum number of components in an eqn */
#define MAX_UNKNOWNS 15		/* default estimate for maximum number of unknowns in model */
#define TOL 1e-9		/* tolerance for comparisons of double numbers */
#define LOG_ZERO_MOLALITY -30	/* molalities <= LOG_ZERO_MOLALITY are considered equal to zero */
#define MIN_TOTAL 1e-25
#define MIN_TOTAL_SS MIN_TOTAL
#define MIN_RELATED_SURFACE MIN_TOTAL*100
#define MIN_RELATED_LOG_ACTIVITY -30
/* ----------------------------------------------------------------------
 *   STRUCTURES
 * ---------------------------------------------------------------------- */
enum SURFACE_TYPE
{ UNKNOWN_DL, NO_EDL, DDL, CD_MUSIC };
enum DIFFUSE_LAYER_TYPE
{ NO_DL, BORKOVEK_DL, DONNAN_DL };
enum SITES_UNITS
{ SITES_ABSOLUTE, SITES_DENSITY };
struct model
{
  int force_prep;
  LDBLE temperature;
  int count_exchange;
  struct master **exchange;

  int count_kinetics;
  struct kinetics *kinetics;

  int count_gas_phase;
  struct phase **gas_phase;

  int count_s_s_assemblage;
  char **s_s_assemblage;

  int count_pp_assemblage;
  struct phase **pp_assemblage;
  char **add_formula;
  LDBLE *si;

  /*int diffuse_layer; */
  /*int edl; */
  enum DIFFUSE_LAYER_TYPE dl_type;
  enum SURFACE_TYPE surface_type;
  int only_counter_ions;
  /*int donnan; */
  LDBLE thickness;
  int count_surface_comp;
  char **surface_comp;
  int count_surface_charge;
  char **surface_charge;
};
EXTERNAL struct model last_model;
EXTERNAL int same_model;
EXTERNAL int same_temperature;

struct name_master
{
  char *name;
  struct master *master;
};
struct name_species
{
  char *name;
  struct species *s;
};
struct name_phase
{
  char *name;
  struct phase *phase;
};
struct punch
{
  int in;
  int new_def;
  struct name_master *totals;
  int count_totals;
  struct name_species *molalities;
  int count_molalities;
  struct name_species *activities;
  int count_activities;
  struct name_phase *pure_phases;
  int count_pure_phases;
  struct name_phase *si;
  int count_si;
  struct name_phase *gases;
  int count_gases;
  struct name_phase *s_s;
  int count_s_s;
  struct name_phase *kinetics;
  int count_kinetics;
  struct name_master *isotopes;
  int count_isotopes;
  struct name_master *calculate_values;
  int count_calculate_values;
  int inverse;
  int sim;
  int state;
  int soln;
  int dist;
  int time;
  int step;
  int rxn;
  int temp;
  int ph;
  int pe;
  int alk;
  int mu;
  int water;
  int high_precision;
  int user_punch;
  int charge_balance;
  int percent_error;
};
EXTERNAL struct punch punch;
/* ----------------------------------------------------------------------
 *   Temperatures
 * ---------------------------------------------------------------------- */
struct temperature
{
  int n_user;
  int n_user_end;
  char *description;
  LDBLE *t;
  int count_t;
};
EXTERNAL struct temperature *temperature;
EXTERNAL int count_temperature;
/* ----------------------------------------------------------------------
 *   Surface
 * --------------------------------------------------------------------- */
struct surface
{
  int n_user;
  int n_user_end;
  int new_def;
  /*int diffuse_layer; */
  /*int edl; */
  int only_counter_ions;
  /*int donnan; */
  enum DIFFUSE_LAYER_TYPE dl_type;
  enum SURFACE_TYPE type;
  enum SITES_UNITS sites_units;
  LDBLE thickness;
  LDBLE debye_lengths;
  LDBLE DDL_viscosity;		/* viscosity relative to pure water */
  LDBLE DDL_limit;		/* limits DDL water to this fraction of bulk water */
  char *description;
  int solution_equilibria;
  int n_solution;
  int count_comps;
  struct surface_comp *comps;
  int count_charge;
  struct surface_charge *charge;
  int related_phases;
  int related_rate;
  int transport;		/* transports comp's and charges if true */
};
struct surface_comp
{
  char *formula;
  struct elt_list *formula_totals;
  LDBLE formula_z;
  LDBLE moles;
  struct master *master;
  struct elt_list *totals;
  LDBLE la;
  int charge;
  LDBLE cb;
  char *phase_name;
  LDBLE phase_proportion;
  char *rate_name;
  LDBLE Dw;			/* diffusion coefficient in water, used in MCD. No transport if 0 */
};
struct surface_charge
{
  char *name;
  LDBLE specific_area;
  LDBLE grams;
  LDBLE charge_balance;
  LDBLE mass_water;
  struct elt_list *diffuse_layer_totals;
  int count_g;
  struct surface_diff_layer *g;	/* stores g and dg/dXd for each ionic charge */
  LDBLE la_psi, la_psi1, la_psi2;
  LDBLE psi, psi1, psi2;
  double capacitance[2];
  double sigma0, sigma1, sigma2, sigmaddl;
};
struct surface_diff_layer
{
  LDBLE charge;
  LDBLE g;
  LDBLE dg;
  LDBLE psi_to_z;
};
EXTERNAL int g_iterations;
EXTERNAL LDBLE G_TOL;
EXTERNAL struct surface *surface;
EXTERNAL struct surface *dbg_surface;
EXTERNAL int count_surface;
EXTERNAL int max_surface;
EXTERNAL struct Charge_Group
{
  LDBLE z;
  LDBLE eq;
} *charge_group;
EXTERNAL int change_surf_count;
EXTERNAL struct Change_Surf
{
  char *comp_name;
  LDBLE fraction;
  char *new_comp_name;
  LDBLE new_Dw;
  int cell_no;
  int next;
} *change_surf;
/* ----------------------------------------------------------------------
 *   Exchange
 * ---------------------------------------------------------------------- */
struct exchange
{
  int n_user;
  int n_user_end;
  int new_def;
  char *description;
  int solution_equilibria;
  int n_solution;
  int count_comps;
  struct exch_comp *comps;
  int related_phases;
  int related_rate;
  int pitzer_exchange_gammas;
};
struct exch_comp
{
  char *formula;
  LDBLE formula_z;
  struct elt_list *formula_totals;
  LDBLE moles;
  struct master *master;
  struct elt_list *totals;
  LDBLE la;
  LDBLE charge_balance;
  char *phase_name;
  LDBLE phase_proportion;
  char *rate_name;
};
EXTERNAL struct exchange *exchange;
EXTERNAL struct exchange *dbg_exchange;
EXTERNAL int count_exchange;
EXTERNAL int max_exchange;
/* ----------------------------------------------------------------------
 *   Kinetics
 * ---------------------------------------------------------------------- */
struct kinetics
{
  int n_user;
  int n_user_end;
  char *description;
  int count_comps;
  struct kinetics_comp *comps;
  int count_steps;
  LDBLE *steps;
  LDBLE step_divide;
  /*char *units; */
  struct elt_list *totals;
  int rk;
  int bad_step_max;
  int use_cvode;
  int cvode_order;
  int cvode_steps;
};
struct kinetics_comp
{
  char *rate_name;
#ifdef SKIP
  char *formula;
#endif
  struct name_coef *list;
  int count_list;
  /*        struct phase *phase; */
  LDBLE tol;
  LDBLE m;
  LDBLE initial_moles;
  LDBLE m0;
  LDBLE moles;
  int count_c_params;
  char **c_params;
  int count_d_params;
  LDBLE *d_params;
};
EXTERNAL struct kinetics *kinetics;
EXTERNAL struct kinetics *dbg_kinetics;
EXTERNAL int count_kinetics;
EXTERNAL int max_kinetics;

struct save_values
{
  LDBLE value;
  int count_subscripts;
  int *subscripts;
};
EXTERNAL int count_save_values;
EXTERNAL struct save_values *save_values;

#ifdef SKIP
struct kin_exch
{
  char *exch_name;
  char *phase_name;
  LDBLE phase_proportion;
};
EXTERNAL struct kin_exch *kin_exch;
EXTERNAL int count_kin_exch;
struct kin_surf
{
  char *surf_name;
  char *phase_name;
  LDBLE phase_proportion;
};
EXTERNAL struct kin_surf *kin_surf;
EXTERNAL int count_kin_surf;
#endif
/*----------------------------------------------------------------------
 *   Save
 *---------------------------------------------------------------------- */
struct save
{
  int solution;
  int n_solution_user;
  int n_solution_user_end;
  int mix;
  int n_mix_user;
  int n_mix_user_end;
  int irrev;
  int n_irrev_user;
  int n_irrev_user_end;
  int pp_assemblage;
  int n_pp_assemblage_user;
  int n_pp_assemblage_user_end;
  int exchange;
  int n_exchange_user;
  int n_exchange_user_end;
  int kinetics;
  int n_kinetics_user;
  int n_kinetics_user_end;
  int surface;
  int n_surface_user;
  int n_surface_user_end;
  int gas_phase;
  int n_gas_phase_user;
  int n_gas_phase_user_end;
  int s_s_assemblage;
  int n_s_s_assemblage_user;
  int n_s_s_assemblage_user_end;
};
EXTERNAL struct save save;
/*----------------------------------------------------------------------
 *   Use
 *---------------------------------------------------------------------- */
struct Use
{
  int solution_in;
  int n_solution_user;
  int n_solution;
  struct solution *solution_ptr;

  int pp_assemblage_in;
  int n_pp_assemblage_user;
  int n_pp_assemblage;
  struct pp_assemblage *pp_assemblage_ptr;

  int mix_in;
  int n_mix_user;
  int n_mix;
  struct mix *mix_ptr;
  int n_mix_user_orig;

  int irrev_in;
  int n_irrev_user;
  int n_irrev;
  struct irrev *irrev_ptr;

  int exchange_in;
  int n_exchange_user;
  int n_exchange;
  struct exchange *exchange_ptr;

  int kinetics_in;
  int n_kinetics_user;
  int n_kinetics;
  struct kinetics *kinetics_ptr;

  int surface_in;
  int n_surface_user;
  int n_surface;
  struct surface *surface_ptr;

  int temperature_in;
  int n_temperature_user;
  int n_temperature;
  struct temperature *temperature_ptr;

  int inverse_in;
  int n_inverse_user;
  int n_inverse;
  struct inverse *inverse_ptr;

  int gas_phase_in;
  int n_gas_phase_user;
  int n_gas_phase;
  struct gas_phase *gas_phase_ptr;

  int s_s_assemblage_in;
  int n_s_s_assemblage_user;
  int n_s_s_assemblage;
  struct s_s_assemblage *s_s_assemblage_ptr;

  int trans_in;
  int advect_in;
};
EXTERNAL struct Use use;
EXTERNAL struct Use *dbg_use;
/*----------------------------------------------------------------------
 *   Copy
 *---------------------------------------------------------------------- */
struct copier
{
  int count;
  int max;
  int *n_user;
  int *start;
  int *end;
};
EXTERNAL struct copier copy_solution;
EXTERNAL struct copier copy_pp_assemblage;
EXTERNAL struct copier copy_exchange;
EXTERNAL struct copier copy_surface;
EXTERNAL struct copier copy_s_s_assemblage;
EXTERNAL struct copier copy_gas_phase;
EXTERNAL struct copier copy_kinetics;
EXTERNAL struct copier copy_mix;
EXTERNAL struct copier copy_irrev;
EXTERNAL struct copier copy_temperature;


/*----------------------------------------------------------------------
 *   Inverse
 *---------------------------------------------------------------------- */
struct inverse
{
  int n_user;
  char *description;
  int new_def;
  int minimal;
  int range;
  int mp;
  LDBLE mp_censor;
  LDBLE range_max;
  LDBLE tolerance;
  LDBLE mp_tolerance;
  int count_uncertainties;
  LDBLE *uncertainties;
  int count_ph_uncertainties;
  LDBLE *ph_uncertainties;
#ifdef SKIP
  LDBLE *alk_uncertainties;
#endif
  LDBLE water_uncertainty;
  int mineral_water;
  int carbon;
  LDBLE *dalk_dph;
  LDBLE *dalk_dc;
  int count_solns;
  int *solns;
  int count_force_solns;
  int *force_solns;
  int count_elts;
  struct inv_elts *elts;
  int count_phases;
  struct inv_phases *phases;
  int count_master_list;
  struct master **master_list;
  int count_redox_rxns;
  int count_isotopes;
  struct inv_isotope *isotopes;
  int count_i_u;
  struct inv_isotope *i_u;
  int count_isotope_unknowns;
  struct isotope *isotope_unknowns;
  char *netpath;
  char *pat;
};
struct inv_elts
{
  char *name;
  struct master *master;
  int row;
  int count_uncertainties;
  LDBLE *uncertainties;
};
struct inv_isotope
{
  char *isotope_name;
  LDBLE isotope_number;
  char *elt_name;
  int count_uncertainties;
  LDBLE *uncertainties;
};
struct inv_phases
{
  char *name;
  struct phase *phase;
  int column;
  int constraint;
  int force;
  int count_isotopes;
  struct isotope *isotopes;
};
EXTERNAL struct inverse *inverse;
EXTERNAL int count_inverse;

/*----------------------------------------------------------------------
 *   Mix
 *---------------------------------------------------------------------- */
struct mix
{
  int n_user;
  int n_user_end;
  char *description;
  int count_comps;
  struct mix_comp *comps;
};
struct mix_comp
{
  int n_solution;
  LDBLE fraction;
};
EXTERNAL struct mix *mix;
EXTERNAL struct mix *dbg_mix;
EXTERNAL int count_mix;
/*----------------------------------------------------------------------
 *   Irreversible reaction
 *---------------------------------------------------------------------- */
struct irrev
{
  int n_user;
  int n_user_end;
  char *description;
  struct name_coef *list;
  struct elt_list *elts;
  LDBLE *steps;
  char *units;
  int count_steps;
  int count_list;
};
struct name_coef
{
  char *name;
  LDBLE coef;
};
EXTERNAL struct irrev *irrev;
EXTERNAL struct irrev *dbg_irrev;
EXTERNAL int count_irrev;
/*----------------------------------------------------------------------
 *   Gas phase
 *---------------------------------------------------------------------- */
struct gas_phase
{
  int n_user;
  int n_user_end;
  char *description;
  int new_def;
  int solution_equilibria;
  int n_solution;
  int type;
  LDBLE total_p;
  LDBLE total_moles;
  LDBLE volume;
  LDBLE temperature;
  int count_comps;
  struct gas_comp *comps;
};
struct gas_comp
{
  struct phase *phase;
  char *name;
  LDBLE p_read;
  LDBLE moles;
  LDBLE initial_moles;
};
EXTERNAL int count_gas_phase;
EXTERNAL int max_gas_phase;
EXTERNAL struct gas_phase *gas_phase;
/*----------------------------------------------------------------------
 *   Solid solution
 *---------------------------------------------------------------------- */
struct s_s_assemblage
{
  int n_user;
  int n_user_end;
  char *description;
  int new_def;
  /*        int type; */
  /*        int solution_equilibria; */
  /*        int n_solution; */
  int count_s_s;
  struct s_s *s_s;
};
struct s_s
{
  char *name;
  struct s_s_comp *comps;
  int count_comps;
  LDBLE total_moles;
  LDBLE dn;
  LDBLE a0, a1;
  LDBLE ag0, ag1;
  int s_s_in;
  int miscibility;
  int spinodal;
  LDBLE tk, xb1, xb2;
  int input_case;
  LDBLE p[4];
};
struct s_s_comp
{
  char *name;
  struct phase *phase;
  LDBLE initial_moles;
  LDBLE moles;
  LDBLE init_moles;
  LDBLE delta;
  LDBLE fraction_x;
  LDBLE log10_lambda;
  LDBLE log10_fraction_x;
  LDBLE dn, dnc, dnb;
};
EXTERNAL int count_s_s_assemblage;
EXTERNAL int max_s_s_assemblage;
EXTERNAL struct s_s_assemblage *s_s_assemblage;
/*----------------------------------------------------------------------
 *   Pure-phase assemblage
 *---------------------------------------------------------------------- */
struct pp_assemblage
{
  int n_user;
  int n_user_end;
  char *description;
  int new_def;
  struct elt_list *next_elt;
  int count_comps;
  struct pure_phase *pure_phases;
};
struct pure_phase
{
  struct phase *phase;
  char *name;
  char *add_formula;
  LDBLE si;
  LDBLE moles;
  LDBLE delta;
  LDBLE initial_moles;
  int force_equality;
  int dissolve_only;
};
EXTERNAL int count_pp_assemblage;
EXTERNAL int max_pp_assemblage;
EXTERNAL struct pp_assemblage *pp_assemblage;
EXTERNAL struct pp_assemblage *dbg_pp_assemblage;
/*----------------------------------------------------------------------
 *   Species_list
 *---------------------------------------------------------------------- */
struct species_list
{
  struct species *master_s;
  struct species *s;
  LDBLE coef;
};
EXTERNAL int count_species_list;
EXTERNAL int max_species_list;
EXTERNAL struct species_list *species_list;
/*----------------------------------------------------------------------
 *   Jacobian and Mass balance lists
 *---------------------------------------------------------------------- */
struct list0
{
  LDBLE *target;
  LDBLE coef;
};
EXTERNAL int count_sum_jacob0;	/* number of elements in sum_jacob0 */
EXTERNAL int max_sum_jacob0;	/* calculated maximum number of elements in sum_jacob0 */
EXTERNAL struct list0 *sum_jacob0;	/* array of pointers to targets and coefficients for array */

struct list1
{
  LDBLE *source;
  LDBLE *target;
};
EXTERNAL int count_sum_mb1;	/* number of elements in sum_mb1 */
EXTERNAL int max_sum_mb1;	/* calculated maximum number of elements in sum_mb1 */
EXTERNAL struct list1 *sum_mb1;	/* array of pointers to sources and targets for mass
				   balance summations with coef = 1.0 */
EXTERNAL int count_sum_jacob1;	/* number of elements in sum_jacob1 */
EXTERNAL int max_sum_jacob1;	/* calculated maximum number of elements in sum_jacob1 */
EXTERNAL struct list1 *sum_jacob1;	/* array of pointers to sources and targets for array
					   equations with coef = 1.0 */
struct list2
{
  LDBLE *source;
  LDBLE *target;
  LDBLE coef;
};
EXTERNAL int count_sum_mb2;	/* number of elements in sum_mb2 */
EXTERNAL int max_sum_mb2;	/* calculated maximum number of elements in sum_mb2 */
EXTERNAL struct list2 *sum_mb2;	/* array of coefficients and pointers to sources and
				   targets for mass balance summations with coef != 1.0 */
EXTERNAL int count_sum_jacob2;	/* number of elements in sum_jacob2 */
EXTERNAL int max_sum_jacob2;	/* calculated maximum number of elements in sum_jacob2 */
EXTERNAL struct list2 *sum_jacob2;	/* array of coefficients and pointers to sources and
					   targets, coef != 1.0 */
EXTERNAL int count_sum_delta;	/* number of elements in sum_delta */
EXTERNAL int max_sum_delta;	/* calculated maximum number of elements in sum_delta */
EXTERNAL struct list2 *sum_delta;	/* array of pointers to sources, targets and coefficients for
					   summing deltas for mass balance equations */
/*----------------------------------------------------------------------
 *   Solution
 *---------------------------------------------------------------------- */
struct solution
{
  int new_def;
  int n_user;
  int n_user_end;
  char *description;
  LDBLE tc;
  LDBLE ph;
  LDBLE solution_pe;
  LDBLE mu;
  LDBLE ah2o;
  LDBLE density;
  LDBLE total_h;
  LDBLE total_o;
  LDBLE cb;
  LDBLE mass_water;
  LDBLE total_alkalinity;
  char *units;
  struct pe_data *pe;
  int default_pe;
  struct conc *totals;
  struct master_activity *master_activity;
  int count_master_activity;
  int count_isotopes;
  struct isotope *isotopes;
  struct master_activity *species_gamma;
  int count_species_gamma;
};
struct master_activity
{
  char *description;
  LDBLE la;
};
struct conc
{
  char *description;
  /*int skip; */
  LDBLE moles;
  LDBLE input_conc;
  char *units;
  char *equation_name;
  struct phase *phase;
  LDBLE phase_si;
  int n_pe;
  char *as;
  LDBLE gfw;
};
struct pe_data
{
  char *name;
  struct reaction *rxn;
};
struct isotope
{
  LDBLE isotope_number;
  char *elt_name;
  char *isotope_name;
  LDBLE total;
  LDBLE ratio;
  LDBLE ratio_uncertainty;
  LDBLE x_ratio_uncertainty;
  struct master *master;
  struct master *primary;
  LDBLE coef;			/* coefficient of element in phase */
};
EXTERNAL struct solution **solution;
EXTERNAL struct solution **dbg_solution;
EXTERNAL int count_solution;
EXTERNAL int max_solution;
struct iso
{
  char *name;
  LDBLE value;
  LDBLE uncertainty;
};
#ifdef MAINSUBS
struct iso iso_defaults[] = {
  {"13C", -10, 1},
  {"13C(4)", -10, 1},
  {"13C(-4)", -50, 5},
  {"34S", 10, 1},
  {"34S(6)", 10, 1},
  {"34S(-2)", -30, 5},
  {"2H", -28, 1},
  {"18O", -5, .1},
  {"87Sr", .71, .01},
  {"11B", 20, 5}
};
int count_iso_defaults = (sizeof (iso_defaults) / sizeof (struct iso));
#else
extern struct iso iso_defaults[];
extern int count_iso_defaults;
#endif
/*----------------------------------------------------------------------
 *   Global solution
 *---------------------------------------------------------------------- */
EXTERNAL char *title_x;
EXTERNAL int new_x;
EXTERNAL char *description_x;
EXTERNAL LDBLE tc_x;
EXTERNAL LDBLE tk_x;
EXTERNAL LDBLE ph_x;
EXTERNAL LDBLE solution_pe_x;
EXTERNAL LDBLE mu_x;
EXTERNAL LDBLE ah2o_x;
EXTERNAL LDBLE density_x;
EXTERNAL LDBLE total_h_x;
EXTERNAL LDBLE total_o_x;
EXTERNAL LDBLE cb_x;
EXTERNAL LDBLE total_ions_x;
EXTERNAL LDBLE mass_water_aq_x;
EXTERNAL LDBLE mass_water_surfaces_x;
EXTERNAL LDBLE mass_water_bulk_x;
EXTERNAL char *units_x;
EXTERNAL struct pe_data *pe_x;
EXTERNAL int count_isotopes_x;
EXTERNAL struct isotope *isotopes_x;
EXTERNAL int default_pe_x;
/*EXTERNAL int diffuse_layer_x;*/
EXTERNAL enum DIFFUSE_LAYER_TYPE dl_type_x;
EXTERNAL LDBLE total_carbon;
EXTERNAL LDBLE total_co2;
EXTERNAL LDBLE total_alkalinity;
EXTERNAL LDBLE gfw_water;
EXTERNAL LDBLE step_x;
EXTERNAL LDBLE kin_time_x;
/*----------------------------------------------------------------------
 *   Transport data
 *---------------------------------------------------------------------- */
EXTERNAL int count_cells;
EXTERNAL int count_shifts;
EXTERNAL int ishift;
EXTERNAL int bcon_first;
EXTERNAL int bcon_last;
EXTERNAL int correct_disp;
EXTERNAL LDBLE tempr;
EXTERNAL LDBLE timest;
EXTERNAL int simul_tr;
EXTERNAL LDBLE diffc;
EXTERNAL LDBLE heat_diffc;
EXTERNAL int cell;
EXTERNAL LDBLE mcd_substeps;
EXTERNAL struct stag_data
{
  int count_stag;
  LDBLE exch_f;
  LDBLE th_m;
  LDBLE th_im;
} *stag_data;
EXTERNAL int print_modulus;
EXTERNAL int punch_modulus;
EXTERNAL int dump_in;
EXTERNAL int dump_modulus;
EXTERNAL int transport_warnings;
EXTERNAL struct cell_data
{
  LDBLE length;
  LDBLE mid_cell_x;
  LDBLE disp;
  LDBLE temp;
  LDBLE por;			/* cell porosities */
  int punch;
  int print;
} *cell_data;
EXTERNAL int multi_Dflag;	/* signals calc'n of multicomponent diffusion */
EXTERNAL LDBLE default_Dw;	/* default species diffusion coefficient in water at 25oC, m2/s */
EXTERNAL LDBLE multi_Dpor;	/* uniform porosity of solid medium */
EXTERNAL LDBLE multi_Dpor_lim;	/* limiting porosity of solid medium where transport stops */
EXTERNAL LDBLE multi_Dn;	/* exponent to calculate pore water diffusion coefficient,
				   Dp = Dw * (multi_Dpor - multi_Dpor_lim)^multi_Dn */
EXTERNAL int cell_no;
/*----------------------------------------------------------------------
 *   Advection data
 *---------------------------------------------------------------------- */
EXTERNAL int count_ad_cells;
EXTERNAL int count_ad_shifts;
EXTERNAL int print_ad_modulus;
EXTERNAL int punch_ad_modulus;
EXTERNAL int *advection_punch, *advection_print;
EXTERNAL LDBLE advection_kin_time;
EXTERNAL LDBLE advection_kin_time_defined;
EXTERNAL int advection_warnings;

/*----------------------------------------------------------------------
 *   Keywords
 *---------------------------------------------------------------------- */
struct key
{
  char *name;
  int keycount;
};
#ifdef MAINSUBS
/* list of valid keywords */
struct key keyword[] = {
  {"eof", 0},
  {"end", 0},
  {"solution_species", 0},
  {"solution_master_species", 0},
  {"solution", 0},
  {"phases", 0},
  {"pure_phases", 0},
  {"reaction", 0},
  {"mix", 0},
  {"use", 0},
  {"save", 0},
  {"exchange_species", 0},
  {"exchange_master_species", 0},
  {"exchange", 0},
  {"surface_species", 0},
  {"surface_master_species", 0},
  {"surface", 0},
  {"reaction_temperature", 0},
  {"inverse_modeling", 0},
  {"gas_phase", 0},
  {"transport", 0},
  {"debug", 0},
  {"selected_output", 0},
  {"select_output", 0},
  {"knobs", 0},
  {"print", 0},
  {"equilibrium_phases", 0},
  {"equilibria", 0},
  {"equilibrium", 0},
  {"pure", 0},
  {"title", 0},
  {"comment", 0},
  {"advection", 0},
  {"kinetics", 0},
  {"incremental_reactions", 0},
  {"incremental", 0},
  {"rates", 0},
  {"solution_s", 0},
  {"user_print", 0},
  {"user_punch", 0},
  {"solid_solutions", 0},
  {"solid_solution", 0},
  {"solution_spread", 0},
  {"spread_solution", 0},
  {"selected_out", 0},
  {"select_out", 0},
  {"user_graph", 0},
  {"llnl_aqueous_model_parameters", 0},
  {"llnl_aqueous_model", 0},
  {"database", 0},
  {"named_analytical_expression", 0},
  {"named_analytical_expressions", 0},
  {"named_expressions", 0},
  {"named_log_k", 0},
  {"isotopes", 0},
  {"calculate_values", 0},
  {"isotope_ratios", 0},
  {"isotope_alphas", 0},
  {"copy", 0},
  {"pitzer", 0}
#ifdef PHREEQC_CPP
  ,
  {"solution_raw", 0},
  {"exchange_raw", 0},
  {"surface_raw", 0},
  {"equilibrium_phases_raw", 0},
  {"kinetics_raw", 0},
  {"solid_solutions_raw", 0},
  {"gas_phase_raw", 0},
  {"reaction_raw", 0},
  {"mix_raw", 0},
  {"reaction_temperature_raw", 0}
#endif /* PHREEQC_CPP */
};
int NKEYS = (sizeof (keyword) / sizeof (struct key));	/* Number of valid keywords */
#else
extern struct key keyword[];
extern int NKEYS;
#endif
EXTERNAL struct key *keyword_hash;
EXTERNAL int new_model, new_exchange, new_pp_assemblage, new_surface,
  new_reaction, new_temperature, new_mix, new_solution, new_gas_phase,
  new_inverse, new_punch, new_s_s_assemblage, new_kinetics, new_copy,
  new_pitzer;
/*----------------------------------------------------------------------
 *   Elements
 *---------------------------------------------------------------------- */
struct element
{
  char *name;			/* element name */
  /*        int in; */
  struct master *master;
  struct master *primary;
  LDBLE gfw;
};
EXTERNAL struct element **elements;
EXTERNAL int count_elements;
EXTERNAL int max_elements;
EXTERNAL struct element *element_h_one;

/*----------------------------------------------------------------------
 *   Element List
 *---------------------------------------------------------------------- */
struct elt_list
{				/* list of name and number of elements in an equation */
  struct element *elt;		/* pointer to element structure */
  LDBLE coef;			/* number of element e's in eqn */
};
EXTERNAL struct elt_list *elt_list;	/* structure array of working space while reading equations
					   names are in "strings", initially in input order */
EXTERNAL int count_elts;	/* number of elements in elt_list = position of next */
EXTERNAL int max_elts;
/*----------------------------------------------------------------------
 *   Reaction
 *---------------------------------------------------------------------- */
struct reaction
{
  LDBLE logk[8];
  struct rxn_token *token;
};
struct rxn_token
{
  struct species *s;
  LDBLE coef;
  char *name;
};
/*----------------------------------------------------------------------
 *   Species
 *---------------------------------------------------------------------- */
struct species
{				/* all data pertinent to an aqueous species */
  char *name;			/* name of species */
  char *mole_balance;		/* formula for mole balance */
  int in;			/* species used in model if TRUE */
  int number;
  struct master *primary;	/* points to master species list, NULL if not primary master */
  struct master *secondary;	/* points to master species list, NULL if not secondary master */
  LDBLE gfw;			/* gram formula wt of species */
  LDBLE z;			/* charge of species */
  LDBLE dw;			/* tracer diffusion coefficient in water at 25oC, m2/s */
  LDBLE equiv;			/* equivalents in exchange species */
  LDBLE alk;			/* alkalinity of species, used for cec in exchange */
  LDBLE carbon;			/* stoichiometric coefficient of carbon in species */
  LDBLE co2;			/* stoichiometric coefficient of C(4) in species */
  LDBLE h;			/* stoichiometric coefficient of H in species */
  LDBLE o;			/* stoichiometric coefficient of O in species */
  LDBLE dha, dhb;		/* WATEQ Debye Huckel a and b-dot */
  LDBLE lk;			/* log10 k at working temperature */
  LDBLE logk[8];		/* log kt0, delh, 6 coefficients analalytical expression */
  DELTA_H_UNIT original_units;	/* enum with original delta H units */
  int count_add_logk;
  struct name_coef *add_logk;
  LDBLE lg;			/* log10 activity coefficient, gamma */
  LDBLE lg_pitzer;		/* log10 activity coefficient, from pitzer calculation */
  LDBLE lm;			/* log10 molality */
  LDBLE la;			/* log10 activity */
  LDBLE dg;			/* gamma term for jacobian */
  LDBLE dg_total_g;
  LDBLE moles;			/* moles in solution; moles/mass_water = molality */
  int type;			/* flag indicating presence in model and types of equations */
  int gflag;			/* flag for preferred activity coef eqn */
  int exch_gflag;		/* flag for preferred activity coef eqn */
  struct elt_list *next_elt;	/* pointer to next element */
  struct elt_list *next_secondary;
  struct elt_list *next_sys_total;
  int check_equation;		/* switch to check equation for charge and element balance */
  struct reaction *rxn;		/* pointer to data base reaction */
  struct reaction *rxn_s;	/* pointer to reaction converted to secondary and primary
				   master species */
  struct reaction *rxn_x;	/* reaction to be used in model */
  LDBLE tot_g_moles;		/* (1 + sum(g)) * moles */
  LDBLE tot_dh2o_moles;		/* sum(moles*g*Ws/Waq) */
  struct species_diff_layer *diff_layer;	/* information related to diffuse layer factors for each
						   surface */
  double cd_music[5];
  double dz[3];
};
struct logk
{				/* Named log K's */
  char *name;			/* name of species */
  LDBLE lk;			/* log10 k at working temperature */
  LDBLE log_k[8];		/* log kt0, delh, 6 coefficients analalytical expression */
  DELTA_H_UNIT original_units;	/* enum with original delta H units */
  int count_add_logk;
  int done;
  struct name_coef *add_logk;
  LDBLE log_k_original[8];	/* log kt0, delh, 5 coefficients analalytical expression */
};
EXTERNAL struct logk **logk;
EXTERNAL int count_logk;
EXTERNAL int max_logk;
struct species_diff_layer
{
  struct surface_charge *charge;
  int count_g;
  LDBLE g_moles;
  LDBLE dg_g_moles;		/* g_moles*dgterm */
  LDBLE dx_moles;
  LDBLE dh2o_moles;		/* moles*g*Ws/Waq */
  LDBLE drelated_moles;		/* for related phase */
};
EXTERNAL char *moles_per_kilogram_string;
EXTERNAL char *pe_string;

EXTERNAL struct species **s;
EXTERNAL int count_s;
EXTERNAL int max_s;

EXTERNAL struct species **s_x;
EXTERNAL int count_s_x;
EXTERNAL int max_s_x;

EXTERNAL struct species *s_h2o;
EXTERNAL struct species *s_hplus;
EXTERNAL struct species *s_h3oplus;
EXTERNAL struct species *s_eminus;
EXTERNAL struct species *s_co3;
EXTERNAL struct species *s_h2;
EXTERNAL struct species *s_o2;
/*----------------------------------------------------------------------
 *   Phases
 *---------------------------------------------------------------------- */
struct phase
{				/* all data pertinent to a pure solid phase */
  char *name;			/* name of species */
  char *formula;		/* chemical formula */
  int in;			/* species used in model if TRUE */
  LDBLE lk;			/* log10 k at working temperature */
  LDBLE logk[8];		/* log kt0, delh, 6 coefficients analalytical expression */
  DELTA_H_UNIT original_units;	/* enum with original delta H units */
  int count_add_logk;
  struct name_coef *add_logk;
  LDBLE moles_x;
  LDBLE p_soln_x;
  LDBLE fraction_x;
  LDBLE log10_lambda, log10_fraction_x;
  LDBLE dn, dnb, dnc;
  LDBLE gn, gntot;
  LDBLE gn_n, gntot_n;

  int type;			/* flag indicating presence in model and types of equations */
  struct elt_list *next_elt;	/* pointer to list of elements in phase */
  struct elt_list *next_sys_total;
  int check_equation;		/* switch to check equation for charge and element balance */
  struct reaction *rxn;		/* pointer to data base reaction */
  struct reaction *rxn_s;	/* pointer to reaction converted to secondary and primary
				   master species */
  struct reaction *rxn_x;	/* reaction to be used in model */
  int in_system;
};
EXTERNAL struct phase **phases;
EXTERNAL int count_phases;
EXTERNAL int max_phases;
/*----------------------------------------------------------------------
 *   Master species
 *---------------------------------------------------------------------- */
struct master
{				/* list of name and number of elements in an equation */
  int in;			/* TRUE if in model, FALSE if out, REWRITE if other mb eq */
  int number;			/* sequence number in list of masters */
  int last_model;		/* saved to determine if model has changed */
  int type;			/* AQ or EX */
  int primary;			/* TRUE if master species is primary */
  LDBLE coef;			/* coefficient of element in master species */
  LDBLE total;			/* total concentration for element or valence state */
  LDBLE isotope_ratio;
  LDBLE isotope_ratio_uncertainty;
  int isotope;
  LDBLE total_primary;
  /*        LDBLE la;  *//* initial guess of master species log activity */
  struct element *elt;		/* element structure */
  LDBLE alk;			/* alkalinity of species */
  LDBLE gfw;			/* default gfw for species */
  char *gfw_formula;		/* formula from which to calcuate gfw */
  struct unknown *unknown;	/* pointer to unknown structure */
  struct species *s;		/* pointer to species structure */
  struct reaction *rxn_primary;	/* reaction writes master species in terms of primary
				   master species */
  struct reaction *rxn_secondary;	/* reaction writes master species in terms of secondary
					   master species */
  struct reaction **pe_rxn;	/* e- written in terms of redox couple (or e-), points
				   to location */
  int minor_isotope;
};
EXTERNAL struct master **master;	/* structure array of master species */
EXTERNAL struct master **dbg_master;
EXTERNAL int count_master;
EXTERNAL int max_master;
/*----------------------------------------------------------------------
 *   Unknowns
 *---------------------------------------------------------------------- */
struct unknown
{
  int type;
  LDBLE moles;
  LDBLE ln_moles;
  LDBLE f;
  LDBLE sum;
  LDBLE delta;
  LDBLE la;
  int number;
  char *description;
  struct master **master;
  struct phase *phase;
  LDBLE si;
  struct gas_phase *gas_phase;
  struct conc *total;
  struct species *s;
  struct exch_comp *exch_comp;
  struct pure_phase *pure_phase;
  struct s_s *s_s;
  struct s_s_comp *s_s_comp;
  int s_s_comp_number;
  int s_s_in;
  struct surface_comp *surface_comp;
  LDBLE related_moles;
  struct unknown *potential_unknown, *potential_unknown1, *potential_unknown2;
  int count_comp_unknowns;
  struct unknown **comp_unknowns;	/* list for CD_MUSIC of comps that contribute to 0 plane mass-balance term */
  struct unknown *phase_unknown;
  struct surface_charge *surface_charge;
  LDBLE mass_water;
  int dissolve_only;
};
EXTERNAL struct unknown **x;
EXTERNAL int count_unknowns;
EXTERNAL int max_unknowns;

EXTERNAL struct unknown *ah2o_unknown;
EXTERNAL struct unknown *alkalinity_unknown;
EXTERNAL struct unknown *carbon_unknown;
EXTERNAL struct unknown *charge_balance_unknown;
EXTERNAL struct unknown *exchange_unknown;
EXTERNAL struct unknown *mass_hydrogen_unknown;
EXTERNAL struct unknown *mass_oxygen_unknown;
EXTERNAL struct unknown *mb_unknown;
EXTERNAL struct unknown *mu_unknown;
EXTERNAL struct unknown *pe_unknown;
EXTERNAL struct unknown *ph_unknown;
EXTERNAL struct unknown *pure_phase_unknown;
EXTERNAL struct unknown *solution_phase_boundary_unknown;
EXTERNAL struct unknown *surface_unknown;
EXTERNAL struct unknown *gas_unknown;
EXTERNAL struct unknown *s_s_unknown;
/*----------------------------------------------------------------------
 *   Reaction work space
 *---------------------------------------------------------------------- */
struct reaction_temp
{
  LDBLE logk[8];
  struct rxn_token_temp *token;
};
struct rxn_token_temp
{				/* data for equations, aq. species or minerals */
  char *name;			/* pointer to a species name (formula) */
  LDBLE z;			/* charge on species */
  struct species *s;
  struct unknown *unknown;
  LDBLE coef;			/* coefficient of species name */
};
EXTERNAL struct reaction_temp trxn;	/* structure array of working space while reading equations
					   species names are in "temp_strings" */
EXTERNAL int count_trxn;	/* number of reactants in trxn = position of next */
EXTERNAL int max_trxn;
struct unknown_list
{
  struct unknown *unknown;
  LDBLE *source;
  LDBLE *gamma_source;
  /*        int row; */
  /*        int col; */
  LDBLE coef;
};
EXTERNAL struct unknown_list *mb_unknowns;
EXTERNAL int count_mb_unknowns;
EXTERNAL int max_mb_unknowns;
/* ----------------------------------------------------------------------
 *   Print
 * ---------------------------------------------------------------------- */
struct prints
{
  int all;
  int initial_solutions;
  int initial_exchangers;
  int reactions;
  int gas_phase;
  int s_s_assemblage;
  int pp_assemblage;
  int surface;
  int exchange;
  int kinetics;
  int totals;
  int eh;
  int species;
  int saturation_indices;
  int irrev;
  int mix;
  int reaction;
  int use;
  int logfile;
  int punch;
  int status;
  int inverse;
  int dump;
  int user_print;
  int headings;
  int user_graph;
  int echo_input;
  int warnings;
  int initial_isotopes;
  int isotope_ratios;
  int isotope_alphas;
  int hdf;
  int alkalinity;
};
EXTERNAL struct prints pr;
EXTERNAL int status_on;
EXTERNAL int count_warnings;

/* ----------------------------------------------------------------------
 *   RATES
 * ---------------------------------------------------------------------- */
struct rate
{
  char *name;
  char *commands;
  int new_def;
  void *linebase;
  void *varbase;
  void *loopbase;
};
EXTERNAL struct rate *rates;
EXTERNAL int count_rates;
EXTERNAL LDBLE rate_m, rate_m0, *rate_p, rate_time, rate_sim_time_start,
  rate_sim_time_end, rate_sim_time, rate_moles, initial_total_time;
EXTERNAL int count_rate_p;
/* ----------------------------------------------------------------------
 *   USER PRINT COMMANDS
 * ---------------------------------------------------------------------- */
EXTERNAL struct rate *user_print;
EXTERNAL struct rate *user_punch;
EXTERNAL char **user_punch_headings;
EXTERNAL int user_punch_count_headings;
#ifdef PHREEQ98
EXTERNAL struct rate *user_graph;
EXTERNAL char **user_graph_headings;
EXTERNAL int user_graph_count_headings;
#endif

/* ----------------------------------------------------------------------
 *   GLOBAL DECLARATIONS
 * ---------------------------------------------------------------------- */
EXTERNAL char error_string[10 * MAX_LENGTH];
EXTERNAL int simulation;
EXTERNAL int state;
EXTERNAL int reaction_step;
EXTERNAL int transport_step;
EXTERNAL int transport_start;
EXTERNAL int advection_step;
EXTERNAL int stop_program;
EXTERNAL int incremental_reactions;

EXTERNAL int count_strings;
EXTERNAL int max_strings;

EXTERNAL LDBLE *array;
EXTERNAL LDBLE *delta;
EXTERNAL LDBLE *residual;

EXTERNAL int input_error;
EXTERNAL int next_keyword;
EXTERNAL int parse_error;
EXTERNAL int paren_count;
EXTERNAL int iterations;
EXTERNAL int gamma_iterations;
EXTERNAL int run_reactions_iterations;

EXTERNAL int max_line;
EXTERNAL char *line;
EXTERNAL char *line_save;

EXTERNAL LDBLE LOG_10;

EXTERNAL int debug_model;
EXTERNAL int debug_prep;
EXTERNAL int debug_set;
EXTERNAL int debug_diffuse_layer;
EXTERNAL int debug_inverse;

EXTERNAL LDBLE inv_tol_default;
EXTERNAL int itmax;
EXTERNAL LDBLE ineq_tol;
EXTERNAL LDBLE convergence_tolerance;
EXTERNAL LDBLE step_size;
EXTERNAL LDBLE pe_step_size;
EXTERNAL LDBLE step_size_now;
EXTERNAL LDBLE pe_step_size_now;
EXTERNAL LDBLE pp_scale;
EXTERNAL LDBLE pp_column_scale;
EXTERNAL int diagonal_scale;	/* 0 not used, 1 used */
EXTERNAL int mass_water_switch;
EXTERNAL int delay_mass_water;
EXTERNAL LDBLE censor;
EXTERNAL int aqueous_only;
EXTERNAL int negative_concentrations;
EXTERNAL int calculating_deriv;
EXTERNAL int numerical_deriv;

EXTERNAL int count_total_steps;
EXTERNAL int phast;
EXTERNAL LDBLE *llnl_temp, *llnl_adh, *llnl_bdh, *llnl_bdot, *llnl_co2_coefs;
EXTERNAL int llnl_count_temp, llnl_count_adh, llnl_count_bdh, llnl_count_bdot,
  llnl_count_co2_coefs;

EXTERNAL char *selected_output_file_name;
EXTERNAL char *dump_file_name;


// MDL: libphree

EXTERNAL int Pfileprint, Pdo_punch; 
EXTERNAL char ** Pbuff; 
EXTERNAL int Pnsim, PCountLine, Pbuff_dim, Pout_dim, Ppunch_dim;
EXTERNAL LDBLE * Ppunch;
EXTERNAL LDBLE * Pout;

struct spread_row
{
  int count;
  int empty, string, number;
  char **char_vector;
  LDBLE *d_vector;
  int *type_vector;
};
struct defaults
{
  LDBLE temp;
  LDBLE density;
  char *units;
  char *redox;
  LDBLE ph;
  LDBLE pe;
  LDBLE water;
  int count_iso;
  struct iso *iso;
};
struct spread_sheet
{
  struct spread_row *heading;
  struct spread_row *units;
  int count_rows;
  struct spread_row **rows;
  struct defaults defaults;
};
#ifdef PHREEQCI_GUI
EXTERNAL struct spread_sheet g_spread_sheet;
#endif

/* ---------------------------------------------------------------------- */
/*
 *   Hash definitions
 */
/*
** Constants
*/

# define SegmentSize                    256
# define SegmentSizeShift          8	/* log2(SegmentSize) */
# define DirectorySize            256
# define DirectorySizeShift      8	/* log2(DirectorySize)  */
# define Prime1                          37
# define Prime2                          1048583
# define DefaultMaxLoadFactor   5


typedef struct Element
{
  /*
   ** The user only sees the first two fields,
   ** as we pretend to pass back only a pointer to ENTRY.
   ** {S}he doesn't know what else is in here.
   */
  char *Key;
  char *Data;
  struct Element *Next;		/* secret from user    */
} Element, *Segment;

typedef struct
{
  short p;			/* Next bucket to be split      */
  short maxp;			/* upper bound on p during expansion */
  long KeyCount;		/* current # keys       */
  short SegmentCount;		/* current # segments   */
  short MinLoadFactor;
  short MaxLoadFactor;
  Segment *Directory[DirectorySize];
} HashTable;

typedef unsigned long Address;

EXTERNAL HashTable *strings_hash_table;
EXTERNAL HashTable *elements_hash_table;
EXTERNAL HashTable *species_hash_table;
EXTERNAL HashTable *phases_hash_table;
EXTERNAL HashTable *keyword_hash_table;
EXTERNAL HashTable *logk_hash_table;
EXTERNAL HashTable *master_isotope_hash_table;

#if defined(PHREEQCI_GUI)
#include "../../phreeqci_gui.h"
#endif /* defined(PHREEQCI_GUI) */

EXTERNAL struct name_coef match_tokens[50];
EXTERNAL int count_match_tokens;
struct master_isotope
{
  char *name;
  struct master *master;
  struct element *elt;
  char *units;
  LDBLE standard;
  LDBLE ratio;
  LDBLE moles;
  int total_is_major;
  int minor_isotope;
};
EXTERNAL int count_master_isotope;
EXTERNAL struct master_isotope **master_isotope;
EXTERNAL int max_master_isotope;
EXTERNAL int initial_solution_isotopes;

#define OPTION_EOF -1
#define OPTION_KEYWORD -2
#define OPTION_ERROR -3
#define OPTION_DEFAULT -4
#define OPT_1 -5

struct calculate_value
{
  char *name;
  LDBLE value;
  char *commands;
  int new_def;
  int calculated;
  void *linebase;
  void *varbase;
  void *loopbase;
};
EXTERNAL int count_calculate_value;
EXTERNAL struct calculate_value **calculate_value;
EXTERNAL int max_calculate_value;
EXTERNAL HashTable *calculate_value_hash_table;

struct isotope_ratio
{
  char *name;
  char *isotope_name;
  LDBLE ratio;
  LDBLE converted_ratio;
};
EXTERNAL int count_isotope_ratio;
EXTERNAL struct isotope_ratio **isotope_ratio;
EXTERNAL int max_isotope_ratio;
EXTERNAL HashTable *isotope_ratio_hash_table;
struct isotope_alpha
{
  char *name;
  char *named_logk;
  LDBLE value;
};
EXTERNAL int count_isotope_alpha;
EXTERNAL struct isotope_alpha **isotope_alpha;
EXTERNAL int max_isotope_alpha;
EXTERNAL HashTable *isotope_alpha_hash_table;

EXTERNAL int phreeqc_mpi_myself;

enum entity_type
{ Solution, Reaction, Exchange, Surface, Gas_phase, Pure_phase, Ss_phase,
    Kinetics, Mix, Temperature, UnKnown };

EXTERNAL int first_read_input;
EXTERNAL char *user_database;
EXTERNAL int pitzer_model, pitzer_pe;
EXTERNAL int full_pitzer, always_full_pitzer, ICON, IC;
EXTERNAL double COSMOT;
EXTERNAL double AW;
EXTERNAL int have_punch_name;

EXTERNAL jmp_buf mark;
EXTERNAL LDBLE *zeros;
EXTERNAL int zeros_max;
#if defined(WIN32)
#include <windows.h>
#endif

#if defined(WIN32_MEMORY_DEBUG)
#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>
#endif

struct system
{
  struct solution *solution;
  struct exchange *exchange;
  struct pp_assemblage *pp_assemblage;
  struct gas_phase *gas_phase;
  struct s_s_assemblage *s_s_assemblage;
  struct kinetics *kinetics;
  struct surface *surface;
};

EXTERNAL double pore_volume;
#endif /* _INC_GLOBAL_H */
