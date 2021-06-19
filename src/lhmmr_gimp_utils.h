#ifndef __SHMMER_GIMP_UTIlS__
#define __SHMMER_GIMP_UTIlS__

#include <vector>
using namespace std;

struct t_annot_region;

enum { VAR_TYPED_TAG, VAR_UNTYPED_TARGET, VAR_UNTYPED_TARGET_IMPUTED};
enum { MATH_MODE_LOG_SPACE, MATH_MODE_LIN_SPACE };

enum { POSTERIOR_MODE_SINGLE_PATH, POSTERIOR_MODE_WEIGHTED_PW_PATH, // ForeBack posterior modes.
		POSTERIOR_MODE_VITERBI_NEIGHBOR_ALLELE_WEIGHTING, POSTERIOR_MODE_VITERBI_PATH_WEIGHT}; // Viterbi posterior heuristic modes.

void run_Imputation_ForeBack_Full_States(char* reference_tag_target_haplocoded_genotypes_fp, char* ref_sample_ids_list_fp,
	char* testing_tag_haplocoded_genotypes_fp, char* testing_sample_ids_list_fp,
	char* known_tag_target_haplocoded_genotypes_fp,
	char* recombination_rate_dir,
	double min_tar2tag_cM,
	char* fore_back_output_fp);

double get_self(double a);

void run_Imputation_ForeBack_State_Reduction_Sliding_Windows_Math_Mode(char* reference_tag_target_haplocoded_genotypes_fp, char* ref_sample_ids_list_fp,
	char* testing_tag_haplocoded_genotypes_fp, char* testing_sample_ids_list_fp,
	char* known_tag_target_haplocoded_genotypes_fp,
	char* recombination_rate_dir,
	double tag_2_tag_distance_cM,
	double N_e,
	char* target_focus_reg_BED_fp,
	bool math_mode,
	double global_scaler_in_log,
	int l_blocks,
	int max_n_tag_vars_per_window,
	int posterior_mode,
	double l_target_2_center_pred_buffer_in_cM,
	char* geno_probs_op_fp);

void run_Imputation_ForeBack_Full_States_Sliding_Windows_Math_Mode(char* reference_tag_target_haplocoded_genotypes_fp, char* ref_sample_ids_list_fp,
	char* testing_tag_haplocoded_genotypes_fp, char* testing_sample_ids_list_fp,
	char* known_tag_target_haplocoded_genotypes_fp,
	char* recombination_rate_dir,
	double tag_2_tag_distance_cM,
	double N_e,
	char* target_focus_reg_BED_fp,
	bool math_mode,
	double global_scaler_in_log,
	char* geno_probs_op_fp);

void run_Imputation_ForeBack_Full_States_Sliding_Windows(char* reference_tag_target_haplocoded_genotypes_fp, char* ref_sample_ids_list_fp,
	char* testing_tag_haplocoded_genotypes_fp, char* testing_sample_ids_list_fp,
	char* known_tag_target_haplocoded_genotypes_fp,
	char* recombination_rate_dir,
	double tag_2_tag_distance_cM,
	double N_e,
	char* target_focus_reg_BED_fp,
	char* geno_probs_op_fp);

void run_Imputation_Viterbi_Full_States_Sliding_Windows(char* reference_tag_target_haplocoded_genotypes_fp, char* ref_sample_ids_list_fp,
	char* testing_tag_haplocoded_genotypes_fp, char* testing_sample_ids_list_fp,
	char* known_tag_target_haplocoded_genotypes_fp,
	char* recombination_rate_dir,
	double tag_2_tag_distance_cM,
	double N_e,
	char* target_focus_reg_BED_fp,
	bool math_mode,
	double global_scaler_in_log,
	char* geno_probs_op_fp);

struct t_var_block
{
	int start_var_i;
	int end_var_i;
	vector<vector<int>*>* haplogroup_haplo_indices;
};

void compare_Viterbi_arrays(char* array1_fp, char* array2_fp);

vector<t_var_block*>* generate_reduced_state_blocks_constant_size_blocks(vector<t_annot_region*>* regs,
	vector<char*>* sample_ids, int selector_info_index, int n_vars_per_block);

vector<t_var_block*>* generate_reduced_state_blocks_constant_size_blocks(char* reference_haplocoded_tag_target_geno_regs_fp,
	char* ref_sample_ids_list_fp, int selector_info_index, int n_vars_per_block);

void run_Imputation_Viterbi_State_Reduction_Sliding_Windows(char* reference_tag_target_haplocoded_genotypes_fp, char* ref_sample_ids_list_fp,
	char* testing_tag_haplocoded_genotypes_fp, char* testing_sample_ids_list_fp,
	char* known_tag_target_haplocoded_genotypes_fp,
	char* recombination_rate_dir,
	double tag_2_tag_distance_cM,
	double N_e,
	char* target_focus_reg_BED_fp,
	bool math_mode,
	double global_scaler_in_log,
	int l_blocks,
	int max_n_tag_vars_per_window,
	int posterior_mode,
	double l_target_2_center_pred_buffer_in_cM,
	char* geno_probs_op_fp);

void run_Imputation_Viterbi_Full_States_Sliding_Windows_State_Pruning(char* reference_tag_target_haplocoded_genotypes_fp, char* ref_sample_ids_list_fp,
	char* testing_tag_haplocoded_genotypes_fp, char* testing_sample_ids_list_fp,
	char* known_tag_target_haplocoded_genotypes_fp,
	char* recombination_rate_dir,
	double tag_2_tag_distance_cM,
	double N_e,
	char* target_focus_reg_BED_fp,
	bool math_mode,
	double global_scaler_in_log,
	double frac_2_max_2_prune,
	char* geno_probs_op_fp);

void get_R2_per_imputed_genotypes(char* imputed_genotypes_fp, char* imputed_sample_ids_list_fp,
	char* known_genotypes_fp, char* known_sample_ids_list_fp);

void get_R2_per_GIMP_4entry_allelic_probs(char* imputed_genotypes_fp, char* imputed_sample_ids_list_fp,
	char* known_genotypes_fp, char* known_sample_ids_list_fp);

void get_n_local_window_tag_variants_per_target(char* testing_tag_haplocoded_genotypes_fp, 
	char* recombination_rate_dir,
	double tag_2_tag_distance_cM,
	char* op_fp);

void get_PR_stats_per_GIMP_4entry_allelic_probs(char* imputed_genotypes_fp, char* imputed_sample_ids_list_fp, 
												char* known_genotypes_fp, char* known_sample_ids_list_fp, 
												bool flag_imputed_probs_are_linear);

void get_PR_stats_per_3entry_genotype_probs(char* imputed_genotypes_fp, char* imputed_sample_ids_list_fp, char* known_genotypes_fp, char* known_sample_ids_list_fp);

#endif // __SHMMER_GIMP_UTIlS__