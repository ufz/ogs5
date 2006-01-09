#ifdef PHREEQC_IDENT
static char const svnidphrqproto[] =
  "$Id: phrqproto.h 4 2009-04-21 17:29:29Z delucia $";
#endif
/* advection.c */
int advection (void);
/* basic.c */
LDBLE system_total (char *total_name, LDBLE * count, char ***names,
		    char ***types, LDBLE ** moles);
int basic_main (char *commands);
void cmd_initialize (void);
void cmd_free (void);
int basic_compile (char *commands, void **lnbase, void **vbase,
		   void **lpbase);
int basic_run (char *commands, void *lnbase, void *vbase, void *lpbase);
/* basicsubs.c */
LDBLE activity (char *species_name);
LDBLE calc_SC (void);
LDBLE calc_logk_n (char *name);
LDBLE calc_logk_p (char *name);
LDBLE calc_logk_s (char *name);
LDBLE calc_surface_charge (char *surface_name);
LDBLE diff_layer_total (char *total_name, char *surface_name);
LDBLE equi_phase (char *phase_name);
LDBLE find_gas_comp (char *gas_comp_name);
LDBLE find_misc1 (char *s_s_name);
LDBLE find_misc2 (char *s_s_name);
LDBLE find_s_s_comp (char *s_s_comp_name);
LDBLE get_calculate_value (char *name);
LDBLE kinetics_moles (char *kinetics_name);
LDBLE log_activity (char *species_name);
LDBLE log_molality (char *species_name);
LDBLE molality (char *species_name);
int reformat_surf (char *comp_name, LDBLE fraction, char *new_comp_name,
		   LDBLE new_Dw, int cell_no);
LDBLE saturation_ratio (char *phase_name);
int saturation_index (char *phase_name, LDBLE * iap, LDBLE * si);
LDBLE solution_sum_secondary (char *total_name);
LDBLE sum_match_gases (char *stemplate, char *name);
LDBLE sum_match_species (char *stemplate, char *name);
LDBLE sum_match_s_s (char *stemplate, char *name);
int match_elts_in_species (char *name, char *stemplate);
int extract_bracket (char **string, char *bracket_string);
LDBLE surf_total (char *total_name, char *surface_name);
int system_species_compare (const void *ptr1, const void *ptr2);
LDBLE system_total (char *total_name, LDBLE * count, char ***names,
		    char ***types, LDBLE ** moles);
int system_total_elements (void);
int system_total_si (void);
int system_total_aq (void);
int system_total_ex (void);
int system_total_surf (void);
int system_total_gas (void);
int system_total_s_s (void);
int system_total_elt (char *total_name);
int system_total_elt_secondary (char *total_name);
LDBLE total (char *total_name);

/* cl1.c */
int cl1 (int k, int l, int m, int n,
	 int nklmd, int n2d,
	 LDBLE * q,
	 int *kode, LDBLE toler,
	 int *iter, LDBLE * x, LDBLE * res, LDBLE * error,
	 LDBLE * cu, int *iu, int *s, int check);

/* default.c */
int close_input_files (void);
int close_output_files (void);
int getc_callback (void *cookie);
int process_file_names (int argc, char *argv[], void **db_cookie,
			void **input_cookie, int log);

/* integrate.c */
int calc_all_g (void);
int calc_init_g (void);
int initial_surface_water (void);
int sum_diffuse_layer (struct surface_charge *surface_charge_ptr1);
int calc_all_donnan (void);
int calc_init_donnan (void);

/* inverse.c */
int inverse_models (void);

/* isotopes.c */
int add_isotopes (struct solution *solution_ptr);
int calculate_values (void);
int calculate_isotope_moles (struct element *elt_ptr,
			     struct solution *solution_ptr,
			     LDBLE total_moles);
LDBLE convert_isotope (struct master_isotope *master_isotope_ptr,
		       LDBLE ratio);
int from_pcil (struct master_isotope *master_isotope_ptr);
int from_permil (struct master_isotope *master_isotope_ptr,
		 LDBLE major_total);
int from_pct (struct master_isotope *master_isotope_ptr, LDBLE major_total);
int from_tu (struct master_isotope *master_isotope_ptr);
struct calculate_value *calculate_value_alloc (void);
int calculate_value_free (struct calculate_value *calculate_value_ptr);
struct calculate_value *calculate_value_search (char *name);
struct calculate_value *calculate_value_store (char *name,
					       int replace_if_found);
struct isotope_alpha *isotope_alpha_alloc (void);
struct isotope_alpha *isotope_alpha_search (char *name);
struct isotope_alpha *isotope_alpha_store (char *name, int replace_if_found);
struct isotope_ratio *isotope_ratio_alloc (void);
struct isotope_ratio *isotope_ratio_search (char *name);
struct isotope_ratio *isotope_ratio_store (char *name, int replace_if_found);
struct master_isotope *master_isotope_store (char *name,
					     int replace_if_found);
struct master_isotope *master_isotope_alloc (void);
struct master_isotope *master_isotope_search (char *name);
int print_initial_solution_isotopes (void);
int print_isotope_ratios (void);
int print_isotope_alphas (void);
int punch_isotopes (void);
int punch_calculate_values (void);
int read_calculate_values (void);
int read_isotopes (void);
int read_isotope_ratios (void);
int read_isotope_alphas (void);

/* kinetics.c */

void cvode_init (void);
int run_reactions (int i, LDBLE kin_time, int use_mix, LDBLE step_fraction);
int set_and_run (int i, int use_mix, int use_kinetics, int nsaver,
		 LDBLE step_fraction);
int set_and_run_wrapper (int i, int use_mix, int use_kinetics, int nsaver,
			 LDBLE step_fraction);
int set_advection (int i, int use_mix, int use_kinetics, int nsaver);
int free_cvode (void);

/* main.c */
int main (int argc, char *argv[]);

/* mainsubs.c */
FILE *file_open (char *query, char *default_name, const char *status,
		 int batch);
int copy_entities (void);
void initialize (void);
int initial_exchangers (int print);
int initial_gas_phases (int print);
int initial_solutions (int print);
int step_save_exch (int n_user);
int step_save_surf (int n_user);
int initial_surfaces (int print);
int k_temp (LDBLE tempc);
int reactions (void);
int saver (void);
int xsolution_save (int k_user);
int do_initialize (void);
int do_status (void);


/* model.c */
int check_residuals (void);
int free_model_allocs (void);
int ineq (int kode);
int model (void);
int jacobian_sums (void);
int mb_gases (void);
int mb_s_s (void);
int mb_sums (void);
int molalities (int allow_overflow);
int reset (void);
int residuals (void);
int set (int initial);
int sum_species (void);
int surface_model (void);

/* p2clib.c */

/* parse.c */
int check_eqn (int association);
int get_charge (char *charge, LDBLE * z);
int get_elt (char **t_ptr, char *element, int *i);
int get_elts_in_species (char **t_ptr, LDBLE coef);
int get_num (char **t_ptr, LDBLE * num);
int get_secondary_in_species (char **t_ptr, LDBLE coef);
int parse_eq (char *eqn, struct elt_list **elt_ptr, int association);

/* pitzer.c */
int gammas_pz (void);
int model_pz (void);
int pitzer (void);
int pitzer_clean_up (void);
int pitzer_init (void);
int pitzer_tidy (void);
int read_pitzer (void);
int set_pz (int initial);

/* prep.c */
int check_same_model (void);
int k_temp (LDBLE tc);
LDBLE k_calc (LDBLE * logk, LDBLE tempk);
int prep (void);
int reprep (void);
int rewrite_master_to_secondary (struct master *master_ptr1,
				 struct master *master_ptr2);
int switch_bases (void);
int write_phase_sys_total (int n);

/* print.c */
int fpunchf (const char *name, const char *format, ...);
int fpunchf_user (int user_index, const char *format, ...);
char *sformatf (const char *format, ...);
int array_print (LDBLE * array_l, int row_count, int column_count,
		 int max_column_count);
int print_all (void);
int print_exchange (void);
int print_gas_phase (void);
int print_master_reactions (void);
int print_reaction (struct reaction *rxn_ptr);
int print_species (void);
int print_surface (void);
int print_user_print (void);
int punch_all (void);
// MDL: prototype for punch_user_punch
int punch_user_punch (void);
/* read.c */
int read_input (void);
int read_conc (int n, int count_mass_balance, char *str);
int *read_list_ints_range (char **ptr, int *count_ints, int positive,
			   int *int_list);
int read_log_k_only (char *ptr, LDBLE * log_k);
int read_number_description (char *ptr, int *n_user, int *n_user_end,
			     char **description);
int check_key (char *str);
int check_units (char *tot_units, int alkalinity, int check_compatibility,
		 const char *default_units, int print);
int find_option (char *item, int *n, const char **list, int count_list,
		 int exact);
int get_option (const char **opt_list, int count_opt_list, char **next_char);
int get_true_false (char *string, int default_value);

/* readtr.c */
int read_transport (void);
int dump (void);
int dump_exchange (int k);
int dump_gas_phase (int k);
int dump_kinetics (int k);
int dump_mix (int k);
int dump_pp_assemblage (int k);
int dump_reaction (int k);
int dump_s_s_assemblage (int k);
int dump_solution (int k);
int dump_surface (int k);

/* spread.c */
int read_solution_spread (void);

/* step.c */
int step (LDBLE step_fraction);
int xsolution_zero (void);
int add_exchange (struct exchange *exchange_ptr);
int add_gas_phase (struct gas_phase *gas_phase_ptr);
int add_kinetics (struct kinetics *kinetics_ptr);
int add_mix (struct mix *mix_ptr);
int add_pp_assemblage (struct pp_assemblage *pp_assemblage_ptr);
int add_reaction (struct irrev *irrev_ptr, int step_number,
		  LDBLE step_fraction);
int add_s_s_assemblage (struct s_s_assemblage *s_s_assemblage_ptr);
int add_solution (struct solution *solution_ptr, LDBLE extensive,
		  LDBLE intensive);
int add_surface (struct surface *surface_ptr);
int add_temperature (struct temperature *temperature_ptr, int step_number);

/* structures.c */
int clean_up (void);
int reinitialize (void);
int copier_add (struct copier *copier_ptr, int n_user, int start, int end);
int copier_free (struct copier *copier_ptr);
int copier_init (struct copier *copier_ptr);
int copy_entities (void);
int element_compare (const void *ptr1, const void *ptr2);
struct element *element_store (char *element);
int elt_list_combine (void);
int elt_list_compare (const void *ptr1, const void *ptr2);
struct elt_list *elt_list_dup (struct elt_list *elt_list_ptr_old);
int elt_list_print (struct elt_list *elt_list_ptr);
struct elt_list *elt_list_save (void);

struct exchange *exchange_alloc (void);
struct exchange *exchange_bsearch (int k, int *n);
int exchange_comp_compare (const void *ptr1, const void *ptr2);
int exchange_copy (struct exchange *exchange_old_ptr,
		   struct exchange *exchange_new_ptr, int n_user_new);
int exchange_compare (const void *ptr1, const void *ptr2);
int exchange_copy_to_last (int n, int n_user);
int exchange_delete (int n_user_old);
int exchange_duplicate (int n_user_old, int n_user_new);
int exchange_init (struct exchange *exchange_ptr, int n_user, int n_user_end,
		   char *description);
int exchange_free (struct exchange *exchange_ptr);
int exchange_ptr_to_user (struct exchange *exchange_old_ptr, int n_user_new);
struct exchange *exchange_replicate (struct exchange *exchange_old_ptr,
				     int n_user_new);
struct exchange *exchange_search (int n_user, int *n, int print);
int exchange_sort (void);

int gas_comp_compare (const void *ptr1, const void *ptr2);
struct gas_phase *gas_phase_alloc (void);
struct gas_phase *gas_phase_bsearch (int k, int *n);
int gas_phase_compare (const void *ptr1, const void *ptr2);
int gas_phase_copy (struct gas_phase *gas_phase_old_ptr,
		    struct gas_phase *gas_phase_new_ptr, int n_user_new);
int gas_phase_copy_to_last (int n, int n_user);
int gas_phase_delete (int n_user_old);
int gas_phase_duplicate (int n_user_old, int n_user_new);
int gas_phase_init (struct gas_phase *gas_phase_ptr, int n_user,
		    int n_user_end, char *description);
int gas_phase_free (struct gas_phase *gas_phase_ptr);
int gas_phase_ptr_to_user (struct gas_phase *gas_phase_ptr_old,
			   int n_user_new);
struct gas_phase *gas_phase_replicate (struct gas_phase *gas_phase_old_ptr,
				       int n_user_new);
struct gas_phase *gas_phase_search (int n_user, int *n);
int gas_phase_sort (void);

enum entity_type get_entity_enum (char *name);

struct inverse *inverse_alloc (void);
int inverse_delete (int i);
int inverse_isotope_compare (const void *ptr1, const void *ptr2);
struct inverse *inverse_search (int n_user, int *n);
int inverse_sort (void);

struct irrev *irrev_bsearch (int k, int *n);
int irrev_copy (struct irrev *irrev_old_ptr, struct irrev *irrev_new_ptr,
		int n_user_new);
int irrev_duplicate (int n_user_old, int n_user_new);
int irrev_free (struct irrev *irrev_ptr);
struct irrev *irrev_search (int n_user, int *n);
int irrev_sort (void);

struct kinetics *kinetics_alloc (void);
struct kinetics *kinetics_bsearch (int k, int *n);
int kinetics_delete (int n_user_old);
int kinetics_comp_duplicate (struct kinetics_comp *kinetics_comp_new_ptr,
			     struct kinetics_comp *kinetics_comp_old_ptr);
int kinetics_compare (const void *ptr1, const void *ptr2);
int kinetics_copy (struct kinetics *kinetics_old_ptr,
		   struct kinetics *kinetics_new_ptr, int n_user_new);
int kinetics_copy_to_last (int n, int n_user);
int kinetics_duplicate (int n_user_old, int n_user_new);
int kinetics_init (struct kinetics *kinetics_ptr, int n_user, int n_user_end,
		   char *description);
int kinetics_free (struct kinetics *kinetics_ptr);
int kinetics_ptr_to_user (struct kinetics *kinetics_ptr_old, int n_user_new);
struct kinetics *kinetics_replicate (struct kinetics *kinetics_old_ptr,
				     int n_user_new);
struct kinetics *kinetics_search (int n_user, int *n, int print);
int kinetics_sort (void);

struct logk *logk_alloc (void);
int logk_copy2orig (struct logk *logk_ptr);
struct logk *logk_store (char *name, int replace_if_found);
struct logk *logk_search (char *name);

struct master *master_alloc (void);
int master_compare (const void *ptr1, const void *ptr2);
int master_delete (char *ptr);
struct master *master_bsearch (const char *ptr);
struct master *master_bsearch_primary (char *ptr);
struct master *master_search (char *ptr, int *n);

struct mix *mix_bsearch (int k, int *n);
int mix_copy (struct mix *mix_old_ptr,
	      struct mix *mix_new_ptr, int n_user_new);
int mix_duplicate (int n_user_old, int n_user_new);
int mix_free (struct mix *mix_ptr);
struct mix *mix_search (int n_user, int *n, int print);
int mix_sort (void);

struct pe_data *pe_data_alloc (void);
struct pe_data *pe_data_dup (struct pe_data *pe_ptr_old);
struct pe_data *pe_data_free (struct pe_data *pe_data_ptr);
int pe_data_store (struct pe_data **pe, const char *token);

struct phase *phase_bsearch (const char *ptr, int *j, int print);
int phase_compare (const void *ptr1, const void *ptr2);
int phase_delete (int i);
struct phase *phase_store (char *name);

struct pp_assemblage *pp_assemblage_alloc (void);
struct pp_assemblage *pp_assemblage_bsearch (int k, int *n);
int pp_assemblage_compare (const void *ptr1, const void *ptr2);
int pp_assemblage_copy (struct pp_assemblage *pp_assemblage_old_ptr,
			struct pp_assemblage *pp_assemblage_new_ptr,
			int n_user_new);
int pp_assemblage_copy_to_last (int n, int n_user);
int pp_assemblage_delete (int n_user_old);
int pp_assemblage_duplicate (int n_user_old, int n_user_new);
int pp_assemblage_free (struct pp_assemblage *pp_assemblage_ptr);
int pp_assemblage_init (struct pp_assemblage *pp_assemblage_ptr, int n_user,
			int n_user_end, char *description);
int pp_assemblage_ptr_to_user (struct pp_assemblage *pp_assemblage_ptr_old,
			       int n_user_new);
struct pp_assemblage *pp_assemblage_replicate (struct pp_assemblage
					       *pp_assemblage_old_ptr,
					       int n_user_new);
struct pp_assemblage *pp_assemblage_search (int n_user, int *n);
int pp_assemblage_sort (void);

int pure_phase_compare (const void *ptr1, const void *ptr2);

struct rate *rate_bsearch (char *ptr, int *j);
int rate_free (struct rate *rate_ptr);
struct rate *rate_search (char *name, int *n);
int rate_sort (void);

struct reaction *rxn_alloc (int ntokens);
struct reaction *rxn_dup (struct reaction *rxn_ptr_old);
LDBLE rxn_find_coef (struct reaction *r_ptr, const char *str);
int rxn_free (struct reaction *rxn_ptr);
int rxn_print (struct reaction *rxn_ptr);

int s_compare (const void *ptr1, const void *ptr2);
int s_delete (int i);
struct species *s_search (char *name);
struct species *s_store (char *name, LDBLE z, int replace_if_found);

struct s_s_assemblage *s_s_assemblage_alloc (void);
struct s_s_assemblage *s_s_assemblage_bsearch (int k, int *n);
int s_s_assemblage_compare (const void *ptr1, const void *ptr2);

int s_s_assemblage_copy (struct s_s_assemblage *s_s_assemblage_old_ptr,
			 struct s_s_assemblage *s_s_assemblage_new_ptr,
			 int n_user_new);
int s_s_assemblage_copy_to_last (int n, int n_user);
int s_s_assemblage_duplicate (int n_user_old, int n_user_new);
int s_s_assemblage_delete (int n_user_old);
int s_s_assemblage_free (struct s_s_assemblage *s_s_assemblage_ptr);
int s_s_assemblage_init (struct s_s_assemblage *s_s_assemblage_ptr,
			 int n_user, int n_user_end, char *description);
int s_s_assemblage_ptr_to_user (struct s_s_assemblage *s_s_assemblage_ptr_old,
				int n_user_new);
struct s_s_assemblage *s_s_assemblage_replicate (struct s_s_assemblage
						 *s_s_assemblage_old_ptr,
						 int n_user_new);
struct s_s_assemblage *s_s_assemblage_search (int n_user, int *n);
int s_s_assemblage_sort (void);
int s_s_compare (const void *ptr1, const void *ptr2);

struct save_values *save_values_bsearch (struct save_values *k, int *n);
int save_values_compare (const void *ptr1, const void *ptr2);
int save_values_sort (void);
int save_values_store (struct save_values *s_v);

int conc_compare (const void *ptr1, const void *ptr2);
int conc_init (struct conc *conc_ptr);
int isotope_compare (const void *ptr1, const void *ptr2);
struct solution *solution_alloc (void);
struct solution *solution_bsearch (int k, int *n, int print);
struct solution *solution_copy (struct solution *solution_old_ptr,
				       int n_user_new);
int solution_copy_to_last (int n, int n_user_new);
int solution_duplicate (int n_user_old, int n_user_new);
int solution_delete (int n_user_old);
int solution_free (struct solution *solution_ptr);
int solution_ptr_to_user (struct solution *solution_old_ptr, int n_user_new);
struct solution *solution_replicate (struct solution *solution_old_ptr,
				     int n_user_new);
int solution_sort (void);

int species_list_compare_alk (const void *ptr1, const void *ptr2);
int species_list_compare_master (const void *ptr1, const void *ptr2);
int species_list_sort (void);

struct Change_Surf *change_surf_alloc (int count);
struct surface *surface_alloc (void);
struct surface *surface_bsearch (int k, int *n);
struct master *surface_get_psi_master (char *name, int plane);
int surface_comp_compare (const void *ptr1, const void *ptr2);
int surface_compare (const void *ptr1, const void *ptr2);
int surface_copy (struct surface *surface_old_ptr,
		  struct surface *surface_new_ptr, int n_user_new);
int surface_copy_to_last (int n, int n_user);
int surface_delete (int n_user_old);
int surface_duplicate (int n_user_old, int n_user_new);
int surface_free (struct surface *surface_ptr);
int surface_init (struct surface *surface_ptr, int n_user, int n_user_end,
		  char *description);
int surface_ptr_to_user (struct surface *surface_ptr_old, int n_user_new);
struct surface *surface_replicate (struct surface *surface_old_ptr,
				   int n_user_new);
struct surface *surface_search (int n_user, int *n, int print);
int surface_sort (void);

int system_duplicate (int i, int save_old);

struct temperature *temperature_bsearch (int k, int *n);
int temperature_copy (struct temperature *temperature_old_ptr,
		      struct temperature *temperature_new_ptr,
		      int n_user_new);
int temperature_duplicate (int n_user_old, int n_user_new);
int temperature_free (struct temperature *temperature_ptr);
struct temperature *temperature_search (int n_user, int *n);
int temperature_sort (void);

int trxn_add (struct reaction *r_ptr, LDBLE coef, int combine);
int trxn_add_phase (struct reaction *r_ptr, LDBLE coef, int combine);
int trxn_combine (void);
int trxn_copy (struct reaction *rxn_ptr);
LDBLE trxn_find_coef (const char *str, int start);
int trxn_print (void);
int trxn_reverse_k (void);
int trxn_sort (void);
int trxn_swap (const char *token);

struct unknown *unknown_alloc (void);
int unknown_delete (int i);
int unknown_free (struct unknown *unknown_ptr);

/* tally.c */
void add_all_components_tally (void);
int build_tally_table (void);
int calc_dummy_kinetic_reaction_tally (struct kinetics *kinetics_ptr);
int diff_tally_table (void);
int extend_tally_table (void);
int free_tally_table (void);
int fill_tally_table (int *n_user, int index_conservative, int n_buffer);
int get_tally_table_rows_columns (int *rows, int *columns);
int get_tally_table_column_heading (int column, int *type, char *string);
int get_tally_table_row_heading (int column, char *string);
int store_tally_table (double *array, int row_dim, int col_dim,
		       double fill_factor);
int zero_tally_table (void);

/* tidy.c */
int add_other_logk (LDBLE * source_k, int count_add_logk,
		    struct name_coef *add_logk);
int add_logks (struct logk *logk_ptr, int repeats);
int select_log_k_expression (LDBLE * source_k, LDBLE * target_k);
int s_s_prep (LDBLE t, struct s_s *s_s_ptr, int print);
int tidy_punch (void);
int tidy_model (void);

/* transport.c */
int transport (void);
int set_initial_moles (int i);
int sum_surface_comp (struct surface *source1, LDBLE f1,
		      struct surface *source2, int k, LDBLE f2,
		      struct surface *target, LDBLE new_Dw);
int reformat_surf (char *comp_name, LDBLE fraction, char *new_comp_name,
		   LDBLE new_Dw, int cell);

/* utilities.c */
int add_elt_list (struct elt_list *elt_list_ptr, LDBLE coef);
int backspace_screen (int spaces);
LDBLE calc_alk (struct reaction *rxn_ptr);
int compute_gfw (const char *string, LDBLE * gfw);
int copy_token (char *token_ptr, char **ptr, int *length);
int dup_print (const char *ptr, int emphasis);
int equal (LDBLE a, LDBLE b, LDBLE eps);
void *free_check_null (void *ptr);
void free_hash_strings (HashTable * Table);
int get_token (char **eqnaddr, char *string, LDBLE * z, int *l);
int hcreate_multi (unsigned Count, HashTable ** HashTable_ptr);
void hdestroy_multi (HashTable * HashTable_ptr);
ENTRY *hsearch_multi (HashTable * Table, ENTRY item, ACTION action);
int islegit (const char c);
void malloc_error (void);
int parse_couple (char *token);
int print_centered (const char *string);
int replace (const char *str1, const char *str2, char *str);
void space (void **ptr, int i, int *max, int struct_size);
void squeeze_white (char *s_l);
int status (int count, const char *str);
void str_tolower (char *str);
void str_toupper (char *str);
int strcmp_nocase (const char *str1, const char *str2);
int strcmp_nocase_arg1 (const char *str1, const char *str2);
char *string_duplicate (const char *token);
char *string_hsave (const char *str);
char *string_pad (char *str, int i);
int string_trim (char *str);
int string_trim_right (char *str);
int string_trim_left (char *str);
LDBLE under (LDBLE xval);
int get_elts_in_species (char **t_ptr, LDBLE coef);
void zero_double (LDBLE *target, int n);
