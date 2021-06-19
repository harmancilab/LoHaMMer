#ifndef __IMPUTATION_UTILS__
#define __IMPUTATION_UTILS__

struct t_IMPUTE2_op_col_info
{
	int name_col_i;
	int posn_col_i;
	int geno_prob_starting_0_based_col_i;
};

enum { MATCH_BY_NAME, MATCH_BY_START_POSN, N_MATCH_BY_IDS };
char** get_match_by_identifiers();

void resample_phased_haplotypes_per_recombination_rates(char* haplocoded_genotype_matrix_fp,
	char* sample_ids_list_fp,
	char* recombination_rate_dir,
	double upsampling_rate,
	double N_e,
	double allele_error_prob,
	char* op_fp);

void resample_phased_haplotypes_per_recombination_rates_multithreaded(char* haplocoded_genotype_matrix_fp,
	char* sample_ids_list_fp,
	char* recombination_rate_dir,
	double upsampling_rate,
	double N_e,
	double allele_error_prob,
	int n_threads,
	char* op_fp);

vector<double*>* get_unique_haplotypes(vector<double*>* lowMAF_containing_haplotypes, vector<int>* n_cnt_per_uniq_haplotypes);

void count_unique_haplotypes(vector<double*>* lowMAF_containing_haplotypes, vector<int>* n_cnt_per_uniq_haplotypes);

vector<t_annot_region*>* load_recombination_rates(char* recombination_rate_fp);

double get_cumulative_recomb_rate_per_variant(t_annot_region* var_reg, vector<t_annot_region*>* recomb_regs);
double get_cumulative_recomb_rate_per_variant_optimized(t_annot_region* var_reg, vector<t_annot_region*>* sorted_recomb_regs);

void get_signal_level_vs_genotype_level_R2(char* imputed_genotype_signal_regs_BED_fp, char* known_genotype_regs_BED_fp, char* sample_ids_list_fp);

void get_unique_haplotype_indices(vector<double*>* haplotype_alleles, vector<vector<int>*>* per_uniq_haplo_indices);

void update_parents_2_child_imputed_genotype_probabilities(char* trio_list_fp, char* genocoded_genotype_matrix_fp, char* sample_ids_list_fp,
															char* imputed_geno_signals_matrix_fp, char* imputed_signal_sample_ids_list_fp);

void assign_haplotype_consistency_scores_per_low_MAF_imputations(char* train_tag_haplocoded_genotype_matrix_fp,
	char* train_target_haplocoded_genotype_matrix_fp,
	char* train_sample_ids_list_fp,
	char* imp_tag_genocoded_genotype_matrix_fp,
	char* imp_target_genocoded_genotype_matrix_fp,
	char* imp_sample_ids_list_fp,
	char* test_tag_genocoded_genotype_matrix_fp,
	char* test_target_genocoded_genotype_matrix_fp,
	char* test_sample_ids_list_fp,
	int n_vic_vars);

void assign_haplotype_consistency_scores_per_low_MAF_imputation_signal(char* train_tag_haplocoded_genotype_matrix_fp,
	char* train_target_haplocoded_genotype_matrix_fp,
	char* train_sample_ids_list_fp,
	char* imp_tag_genocoded_genotype_matrix_fp,
	char* imp_target_genocoded_genotype_matrix_fp,
	char* imp_target_geno_signal_matrix_fp,
	char* imp_sample_ids_list_fp,
	char* test_tag_genocoded_genotype_matrix_fp,
	char* test_target_genocoded_genotype_matrix_fp,
	char* test_sample_ids_list_fp,
	int n_vic_vars);

void assign_posterior_haplotype_consistency_scores_per_imputation_signal(char* train_tag_haplocoded_genotype_matrix_fp,
	char* train_target_haplocoded_genotype_matrix_fp,
	char* train_sample_ids_list_fp,
	char* imp_tag_genocoded_genotype_matrix_fp,
	char* imp_target_geno_signal_matrix_fp,
	char* imp_sample_ids_list_fp,
	char* test_target_genocoded_genotype_matrix_fp,
	char* test_sample_ids_list_fp,
	double exponential_allelic_prior_weight,
	int n_vic_vars,
	char* probs_5col_fp);

void compute_imputation_stats_per_5col_genotype_probs_w_PR_curve_info(char* imputed_chr_id,
	char* imputed_5col_genotype_probs_fp,
	char* imputed_5col_genotype_probs_sample_ids_list_fp,
	char* known_genotype_regs_fp,
	char* known_genotype_regs_sample_ids_list_fp,
	char* match_by_str,
	double min_score_threshold,
	double max_score_threshold,
	double score_delta,
	char* op_stats_fp);

void analyze_low_MAF_haplotype_statistics_per_target_variants(char* train_tag_haplocoded_genotype_matrix_fp,
	char* train_target_haplocoded_genotype_matrix_fp,
	char* train_sample_ids_list_fp,
	int n_vic_vars);

void compute_imputation_stats_per_5col_genotype_probs(char* imputed_chr_id,
	char* imputed_5col_genotype_probs_fp,
	char* imputed_5col_genotype_probs_sample_ids_list_fp,
	char* known_genotype_regs_fp,
	char* known_genotype_regs_sample_ids_list_fp,
	char* match_by_str,
	char* op_stats_fp);

void train_vicinity_based_LMSE_imputation_model_multithreaded(int start_pos, int end_pos,
	int n_threads,
	int n_tags_vars_per_side,

	// The training data.
	char* tag_genotype_matrix_fp,
	char* tag_sample_ids_list_fp,
	char* target_genotype_matrix_fp,
	char* target_sample_ids_list_fp,

	// The testing data.
	char* testing_tag_genotype_matrix_fp,
	char* testing_tag_sample_ids_list_fp,
	char* testing_target_genotype_matrix_fp,
	char* testing_target_sample_ids_list_fp,

	char* op_dir);

void train_vicinity_based_LMSE_imputation_model(int start_pos, int end_pos,
	int n_tags_vars_per_side,

	// The training data.
	char* tag_genotype_matrix_fp,
	char* tag_sample_ids_list_fp,
	char* target_genotype_matrix_fp,
	char* target_sample_ids_list_fp,

	// The testing data.
	char* testing_tag_genotype_matrix_fp,
	char* testing_tag_sample_ids_list_fp,
	char* testing_target_genotype_matrix_fp,
	char* testing_target_sample_ids_list_fp,

	char* op_dir);

void save_train_test_matrices_per_variant(char* roi_target_vars_BED_fp,
	int n_half_vicinity_tags,
	char* training_tag_geno_sig_regs_fp,
	char* training_target_geno_sig_regs_fp,
	char* training_sample_ids_list_fp,
	char* testing_tag_geno_sig_regs_fp,
	char* testing_target_geno_sig_regs_fp,
	char* testing_sample_ids_list_fp,
	char* op_prefix);

void generate_trio_genotyping_imputability_stats(char* genotype_signal_matrix_fp,
												char* sample_ids_fp,
												char* trio_fmc_sample_ids_fp);

void extract_genotype_probability_info_per_GT_entry_id_VCF(char* vcf_fp,							// This is the VCF file from which the genotypes are read.
	char* vcf_sample_ids_list_fp,		// This is the sample ids list file path.
	char* target_geno_sig_regs_fp,
	char* target_sample_ids_list_fp,
	char* GT_entry_id,
	char* op_prefix);

void extract_genotype_probability_info_per_IMPUTE2_probabilities_output(char* IMPUTE2_output_fp,
	t_IMPUTE2_op_col_info* IMPUTE2_op_col_info,
	char* target_geno_sig_regs_fp,
	char* target_sample_ids_list_fp,
	char* chrom,
	char* op_prefix);

void extract_genotyped_values_per_continuous_genotypes_signal_w_preset_thresholds(char* baseline_model_output_file_path,
	char* input_genotype_matrix_sample_ids_list_fp,
	bool is_signal_haplocoded,
	bool save_haplocoded,
	char* output_genotype_matrix_bed_fp);

void extract_genotyped_values_per_continuous_genotypes_signal_w_adaptive_thresholds_per_training_signal(char* testing_geno_signals_file_path,
	char* testing_sample_ids_list_fp,
	char* training_geno_signals_file_path,
	char* training_geno_signals_sample_ids_list_fp,
	char* training_geno_regs_file_path,
	char* training_geno_regs_sample_ids_list_fp,
	bool is_signal_haplocoded,
	bool save_haplocoded,
	char* output_genotype_matrix_bed_fp,
	char* prob_5col_op_fp);

void extract_IMPUTE2_data_input_files_per_haplocoded_genotype_matrix(char* reference_haplotype_genotype_matrix_matbed_fp, char* reference_haplotype_genotype_matrix_sample_ids_list_fp,
	char* input_genotype_matrix_matbed_fp, char* input_genotype_matrix_sample_ids_list_fp,
	char* h_option_fp,
	char* l_option_fp,
	char* g_option_fp,
	char* strand_g_option_fp);

void extract_BEAGLE_data_input_files_per_haplocoded_genotype_matrix(char* reference_haplotype_genotype_matrix_matbed_fp, char* reference_haplotype_genotype_matrix_sample_ids_list_fp,
	char* input_genotype_matrix_matbed_fp, char* input_genotype_matrix_sample_ids_list_fp,
	char* ref_option_fp,
	char* gt_option_fp,
	bool save_phased_gt_option);

void extract_MACH_data_input_files_per_haplocoded_genotype_matrix(char* reference_haplotype_genotype_matrix_matbed_fp, char* reference_haplotype_genotype_matrix_sample_ids_list_fp,
	char* input_genotype_matrix_matbed_fp, char* input_genotype_matrix_sample_ids_list_fp,
	char* d_option_fp,
	char* p_option_fp,
	char* h_option_fp,
	char* s_option_fp);

void compute_imputation_stats_per_IMPUTED_genotypes(char* imputed_genotyped_regs_fp, char* imputed_genotyped_regs_sample_ids_list_fp,
	char* known_genotype_regs_fp, char* known_genotype_regs_sample_ids_list_fp);

void extract_genotyped_values_per_IMPUTE2_probabilities_output(char* IMPUTE2_output_fp, char* input_genotype_matrix_sample_ids_list_fp, char* output_genotype_matrix_bed_fp,
	double min_prob_2_genotype,
	t_IMPUTE2_op_col_info* IMPUTE2_op_col_info,
	char* chrom);

#endif // __IMPUTATION_UTILS__
