#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "lhmmr_variation_tools.h"
#include "lhmmr_genomics_coords.h"
#include "lhmmr_histogram.h"
#include "lhmmr_imputation_utils.h"
#include "lhmmr_ansi_string.h"
#include "lhmmr_annot_region_tools.h"
#include "lhmmr_file_utils.h"
#include "lhmmr_nucleotide.h"
#include "lhmmr_xlog_math.h"
#include "lhmmr_linear_math.h"
#include "lhmmr_gimp_utils.h"

#include <algorithm>
#include <vector>
using namespace std;

bool __DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__ = false;

void generate_reduced_state_blocks_per_max_unique_haplotypes(char* reference_haplocoded_tag_target_geno_regs_fp,
	char* ref_sample_ids_list_fp, int max_unique_haplotypes_per_block)
{
	vector<t_annot_region*>* reference_haplocoded_tag_target_geno_regs = load_variant_signal_regions_wrapper(reference_haplocoded_tag_target_geno_regs_fp,
		ref_sample_ids_list_fp);
	vector<char*>* ref_sample_ids = buffer_file(ref_sample_ids_list_fp);

	fprintf(stderr, "Generating reduced state haplotype blocks for %d variant regions with %d samples.\n",
		reference_haplocoded_tag_target_geno_regs->size(),
		ref_sample_ids->size());
	vector<int*>* block_coordinates = new vector<int*>();

	int min_n_variants_per_block = 0;

	// Setup the tag matrix for the current vicinity.
	int n_ref_haplotypes = 2 * ref_sample_ids->size();

	double** per_haplo_per_var_alleles = new double*[n_ref_haplotypes + 2];
	for (int hap_i = 0; hap_i < n_ref_haplotypes; hap_i++)
	{
		per_haplo_per_var_alleles[hap_i] = new double[2 * reference_haplocoded_tag_target_geno_regs->size() + 5];

		// Reset all entries to -1.
		for (int i_var = 0;
			i_var < reference_haplocoded_tag_target_geno_regs->size();
			i_var++)
		{
			per_haplo_per_var_alleles[hap_i][i_var] = -1;
		} // i_var loop.				
	} // hap_i loop.

	  // Set the haplotypes.
	for (int j_var = 0;
		j_var < reference_haplocoded_tag_target_geno_regs->size();
		j_var++)
	{
		void** cur_var_info = (void**)(reference_haplocoded_tag_target_geno_regs->at(j_var)->data);
		char* cur_train_sample_haplocoded_geno = (char*)(cur_var_info[0]);

		for (int i_s = 0; i_s < ref_sample_ids->size(); i_s++)
		{
			int rel_var_i = j_var;
			per_haplo_per_var_alleles[i_s * 2][rel_var_i] = get_allele_per_haplotype(cur_train_sample_haplocoded_geno[i_s], 0);
			per_haplo_per_var_alleles[i_s * 2 + 1][rel_var_i] = get_allele_per_haplotype(cur_train_sample_haplocoded_geno[i_s], 1);
		} // i_s loop.
	} // j_var loop.

	vector<double*>* all_haplotypes = new vector<double*>();
	for (int i_hap = 0; i_hap < n_ref_haplotypes; i_hap++)
	{
		all_haplotypes->push_back(per_haplo_per_var_alleles[i_hap]);
	} // i_hap loop.	

	fprintf(stderr, "Setup %d haplotypes on %d regions.\n", all_haplotypes->size(), reference_haplocoded_tag_target_geno_regs->size());

	int i_reg = 0;
	while (i_reg < reference_haplocoded_tag_target_geno_regs->size())
	{
		int cur_l = min_n_variants_per_block;

		// Extend until we get a decent number of haplotypes.
		while (1)
		{
			// Setup the current haplotypes.
			for (int i_hap = 0; i_hap < all_haplotypes->size(); i_hap++)
			{
				all_haplotypes->at(i_hap)[i_reg + cur_l] *= -1;
			} // i_hap loop.

			  // Generate the unique haplotype blocks.
			vector<int>* n_cnt_per_uniq_haplotypes = new vector<int>();
			count_unique_haplotypes(all_haplotypes, n_cnt_per_uniq_haplotypes);

			// Reset the haplotypes.
			for (int i_hap = 0; i_hap < all_haplotypes->size(); i_hap++)
			{
				all_haplotypes->at(i_hap)[i_reg + cur_l] *= -1;
			} // i_hap loop.

			if (n_cnt_per_uniq_haplotypes->size() >= max_unique_haplotypes_per_block)
			{
				int* new_block = new int[2];
				new_block[0] = i_reg;
				new_block[1] = i_reg + cur_l - 1;
				block_coordinates->push_back(new_block);
				i_reg = i_reg + cur_l;
				break;
			}

			cur_l++;
		}

	} // i_reg loop.

	//return(block_coordinates);
}

vector<t_var_block*>* generate_reduced_state_blocks_constant_size_blocks(char* reference_haplocoded_tag_target_geno_regs_fp,
	char* ref_sample_ids_list_fp,
	int genotype_sequence_index_in_info,
	int n_vars_per_block)
{
	vector<t_annot_region*>* reference_haplocoded_tag_target_geno_regs = load_variant_signal_regions_wrapper(reference_haplocoded_tag_target_geno_regs_fp,
		ref_sample_ids_list_fp);
	vector<char*>* ref_sample_ids = buffer_file(ref_sample_ids_list_fp);

	vector<t_var_block*>* blocks = generate_reduced_state_blocks_constant_size_blocks(reference_haplocoded_tag_target_geno_regs,
																						ref_sample_ids,
																						genotype_sequence_index_in_info,
																						n_vars_per_block);

	return(blocks);
}

vector<t_var_block*>* generate_reduced_state_blocks_constant_size_blocks(vector<t_annot_region*>* reference_haplocoded_tag_target_geno_regs,
																		vector<char*>* ref_sample_ids, 
																		int genotype_sequence_index_in_info, 
																		int n_vars_per_block)
{
	fprintf(stderr, "Generating reduced state haplotype blocks for %d variant regions with %d samples using blocks of %d variants.\n", 
			reference_haplocoded_tag_target_geno_regs->size(), 			
			ref_sample_ids->size(),
			n_vars_per_block);

	// Setup the tag matrix for the current vicinity.
	int n_ref_haplotypes = 2 * ref_sample_ids->size();

	vector<t_var_block*>* all_haplo_var_blocks = new vector<t_var_block*>();

	double** per_haplo_per_var_alleles = new double*[n_ref_haplotypes + 2];
	for (int hap_i = 0; hap_i < n_ref_haplotypes; hap_i++)
	{
		per_haplo_per_var_alleles[hap_i] = new double[2 * reference_haplocoded_tag_target_geno_regs->size() + 5];

		// Reset all entries to -1.
		for (int i_var = 0;
			i_var < reference_haplocoded_tag_target_geno_regs->size();
			i_var++)
		{
			per_haplo_per_var_alleles[hap_i][i_var] = -1;
		} // i_var loop.				
	} // hap_i loop.

	// Set the haplotypes.
	for (int j_var = 0;
		j_var < reference_haplocoded_tag_target_geno_regs->size();
		j_var++)
	{
		void** cur_var_info = (void**)(reference_haplocoded_tag_target_geno_regs->at(j_var)->data);
		char* cur_train_sample_haplocoded_geno = (char*)(cur_var_info[genotype_sequence_index_in_info]);

		for (int i_s = 0; i_s < ref_sample_ids->size(); i_s++)
		{
			int rel_var_i = j_var;
			per_haplo_per_var_alleles[i_s * 2][rel_var_i] = get_allele_per_haplotype(cur_train_sample_haplocoded_geno[i_s], 0);
			per_haplo_per_var_alleles[i_s * 2 + 1][rel_var_i] = get_allele_per_haplotype(cur_train_sample_haplocoded_geno[i_s], 1);
		} // i_s loop.
	} // j_var loop.

	vector<double*>* all_haplotypes = new vector<double*>();
	for (int i_hap = 0; i_hap < n_ref_haplotypes; i_hap++)
	{
		all_haplotypes->push_back(per_haplo_per_var_alleles[i_hap]);
	} // i_hap loop.	

	fprintf(stderr, "Setup %d haplotypes on %d regions.\n", all_haplotypes->size(), reference_haplocoded_tag_target_geno_regs->size());

	// Copy the current haplotypes.
	vector<double*>* cur_block_haplotypes = new vector<double*>();
	for (int i_hap = 0; i_hap < all_haplotypes->size(); i_hap++)
	{
		double* cur_block = new double[n_vars_per_block+2];
		cur_block_haplotypes->push_back(cur_block);
	} // i_hap loop.

	int i_reg = 0;
	while(i_reg < reference_haplocoded_tag_target_geno_regs->size())
	{
		// Copy the current haplotypes.
		if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
		{
			fprintf(stderr, "Copying the current block haplotypes @ i_reg=%d\n", i_reg);
		}

		int n_var_per_cur_block = 0;
		for (int i_hap = 0; i_hap < all_haplotypes->size(); i_hap++)
		{
			n_var_per_cur_block = 0;
			for (int i_var = i_reg; 
				i_var < MIN(reference_haplocoded_tag_target_geno_regs->size(), i_reg + n_vars_per_block); 
				i_var++)
			{
				cur_block_haplotypes->at(i_hap)[i_var - i_reg] = all_haplotypes->at(i_hap)[i_var];
				n_var_per_cur_block++;
			} // i_var loop.

			cur_block_haplotypes->at(i_hap)[n_var_per_cur_block] = -1;
		} // i_hap loop.

		// Generate the unique haplotype blocks; assign the 
		if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
		{
			fprintf(stderr, "Computing the unique haplotypes @ i_reg=%d\n", i_reg);
		}

		vector<vector<int>*>* per_uniq_haplo_indices = new vector<vector<int>*>();
		get_unique_haplotype_indices(cur_block_haplotypes, per_uniq_haplo_indices);

		t_var_block* cur_var_block = new t_var_block();
		cur_var_block->start_var_i = i_reg + 1; // These are 1-based.
		cur_var_block->end_var_i = i_reg + n_var_per_cur_block; // These are 1-based.
		cur_var_block->haplogroup_haplo_indices = per_uniq_haplo_indices;

		all_haplo_var_blocks->push_back(cur_var_block);

		if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
		{
			fprintf(stderr, "Block %d-%d: %d unique haplotypes.\n", i_reg, i_reg + n_var_per_cur_block - 1, per_uniq_haplo_indices->size());
		}

		// Move to the next block.
		i_reg = i_reg + n_var_per_cur_block;
	} // i_reg loop.

	for (int i_hap = 0; i_hap < all_haplotypes->size(); i_hap++)
	{
		delete[] cur_block_haplotypes->at(i_hap);
	} // i_hap loop.
	delete cur_block_haplotypes;

	for (int hap_i = 0; hap_i < n_ref_haplotypes; hap_i++)
	{
		delete[] per_haplo_per_var_alleles[hap_i];		
	} // hap_i loop.
	delete[] per_haplo_per_var_alleles;

	delete all_haplotypes;

	return(all_haplo_var_blocks);
}

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
	char* geno_probs_op_fp)
{
	// Make sure the target-2-center is not larger than the window length
	// THIS SHOULD BE TURNED BACK ON.
	if ((2 * l_target_2_center_pred_buffer_in_cM) > tag_2_tag_distance_cM)
	{
		//l_target_2_center_pred_buffer_in_cM = tag_2_tag_distance_cM / 2;
	}

	//double N_e = pow(10, 4);
	double allele_err = pow(10, -4);

	double(*XSUM)(double num1, double num2) = NULL;
	double(*XMUL)(double num1, double num2) = NULL;
	double(*XDIV)(double num1, double num2) = NULL;
	double(*XSUB)(double num1, double num2) = NULL;
	bool(*XCOMPARE)(double num1, double num2) = NULL;
	double(*XCONVERT_LIN)(double num1) = NULL;
	double(*XCONVERT_LOG)(double num1) = NULL;
	double(*XCONVERT_2_LOG)(double num1) = NULL;
	double(*XCONVERT_2_LIN)(double num1) = NULL;
	double ZERO = 0.0;

	if (math_mode == MATH_MODE_LOG_SPACE)
	{
		fprintf(stderr, "Math mode log.\n");
		XSUM = xlog_sum;
		XMUL = xlog_mul;
		XDIV = xlog_div;
		XSUB = xlog_sub;
		XCOMPARE = xlog_comp;
		XCONVERT_LIN = xlog;
		XCONVERT_LOG = get_self;
		XCONVERT_2_LOG = get_self;
		XCONVERT_2_LIN = exp;
		ZERO = xlog(0);
	}
	else if (math_mode == MATH_MODE_LIN_SPACE)
	{
		fprintf(stderr, "Math mode linear.\n");
		XSUM = lin_sum;
		XMUL = lin_mul;
		XDIV = lin_div;
		XSUB = lin_sub;
		XCOMPARE = lin_compare;
		XCONVERT_LIN = get_self;
		XCONVERT_LOG = exp;
		XCONVERT_2_LOG = xlog;
		XCONVERT_2_LIN = get_self;
		ZERO = 0.0;
	}

	if (posterior_mode == POSTERIOR_MODE_SINGLE_PATH)
	{
		fprintf(stderr, "Computing single path posterior probabilities.\n");
	}
	else if (posterior_mode == POSTERIOR_MODE_WEIGHTED_PW_PATH)
	{
		fprintf(stderr, "Computing weighted pairwise path posterior probabilities.\n");
	}
	else
	{
		fprintf(stderr, "Could not determine the posterior probability mode.\n");
		exit(0);
	}

	if (max_n_tag_vars_per_window > 0)
	{
		fprintf(stderr, "Using maximum of %d tag variants per window.\n", max_n_tag_vars_per_window);
	}

	double global_lin_scaler = exp(global_scaler_in_log);
	if (math_mode == MATH_MODE_LOG_SPACE)
	{
		global_lin_scaler = xlog(1.0);
	}

	fprintf(stderr, "Running fore-back sliding window with tag2tag distance of %.3f cMs (%.3f target-2-center buffer) with effective pop. size N_e=%d (Math_Mode: %d; Scaler: %.3f)\n", 
			tag_2_tag_distance_cM, l_target_2_center_pred_buffer_in_cM, (int)N_e, math_mode, global_scaler_in_log);

	// This is the tag+target region genotypes of the reference panel.
	vector<t_annot_region*>* reference_haplocoded_tag_target_geno_regs = load_variant_signal_regions_wrapper(reference_tag_target_haplocoded_genotypes_fp, ref_sample_ids_list_fp);
	vector<char*>* ref_sample_ids = buffer_file(ref_sample_ids_list_fp);
	int n_ref_haplotypes = ref_sample_ids->size() * 2;
	fprintf(stderr, "Loaded %d variant regions for %d reference samples, %d haplotypes.\n",
		reference_haplocoded_tag_target_geno_regs->size(),
		ref_sample_ids->size(),
		n_ref_haplotypes);

	// Reset the scores to 0 to indicates all untyped.
	for (int i_reg = 0; i_reg < reference_haplocoded_tag_target_geno_regs->size(); i_reg++)
	{
		reference_haplocoded_tag_target_geno_regs->at(i_reg)->score = VAR_UNTYPED_TARGET;
	} // i_reg loop.

	  // These are the typed genotypes of the testing panel.
	vector<t_annot_region*>* testing_haplocoded_tag_geno_regs = load_variant_signal_regions_wrapper(testing_tag_haplocoded_genotypes_fp, testing_sample_ids_list_fp);
	vector<char*>* testing_sample_ids = buffer_file(testing_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d testing variant regions for %d testing samples.\n", testing_haplocoded_tag_geno_regs->size(), testing_sample_ids->size());

	// Assign the typed variants to the variants of the reference panel.
	fprintf(stderr, "Assigning testing tag regions with the reference variants.\n");
	vector<t_annot_region*>* intersects = intersect_annot_regions(testing_haplocoded_tag_geno_regs, reference_haplocoded_tag_target_geno_regs, true);
	fprintf(stderr, "Processing %d intersections.\n", intersects->size());
	for (int i_int = 0; i_int < intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* testing_reg = int_info->src_reg;
		t_annot_region* ref_reg = int_info->dest_reg;

		// Set this reference variant as "typed". 
		ref_reg->score = VAR_TYPED_TAG;

		void** testing_reg_info = (void**)(testing_reg->data);
		void** ref_reg_info = (void**)(ref_reg->data);
		testing_reg_info[1] = ref_reg_info[0];

		delete int_info;
	} // i_int loop.
	delete_annot_regions(intersects);

	vector<t_annot_region*>* focus_regs = NULL;
	if (check_file(target_focus_reg_BED_fp))
	{
		focus_regs = load_BED(target_focus_reg_BED_fp);
		fprintf(stderr, "Found %d focus regions.\n", focus_regs->size());
	}

	// Assign the haplocoded regions.
	vector<t_annot_region*>* known_haplocoded_tag_target_geno_regs = NULL;
	if (check_file(known_tag_target_haplocoded_genotypes_fp))
	{
		fprintf(stderr, "Loading known genotypes from %s\n", known_tag_target_haplocoded_genotypes_fp);

		known_haplocoded_tag_target_geno_regs = load_variant_signal_regions_wrapper(known_tag_target_haplocoded_genotypes_fp, testing_sample_ids_list_fp);
		fprintf(stderr, "Loaded %d known tag+target variants for %d samples.\n", known_haplocoded_tag_target_geno_regs->size(), testing_sample_ids->size());

		fprintf(stderr, "Assigning the known genotype regions to reference genotypes.\n");
		intersects = intersect_annot_regions(reference_haplocoded_tag_target_geno_regs, known_haplocoded_tag_target_geno_regs, true);
		fprintf(stderr, "Processing %d known genotype reference intersects.\n", intersects->size());
		for (int i_int = 0; i_int < intersects->size(); i_int++)
		{
			t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
			t_annot_region* ref_reg = int_info->src_reg;
			t_annot_region* known_reg = int_info->dest_reg;

			void** ref_reg_info = (void**)(ref_reg->data);
			void** known_reg_info = (void**)(known_reg->data);
			ref_reg_info[1] = known_reg_info[0];

			delete int_info;
		} // i_int loop.

		delete_annot_regions(intersects);
	}

	// Set the imputed genotypes probabilities for untyped variants.
	for (int i_var = 0; i_var < reference_haplocoded_tag_target_geno_regs->size(); i_var++)
	{
		void** cur_ref_var_info = (void**)(reference_haplocoded_tag_target_geno_regs->at(i_var)->data);
		void** new_info = new void*[5];
		new_info[0] = cur_ref_var_info[0];
		new_info[1] = cur_ref_var_info[1];

		double*** imputed_allele_prob_per_sample_per_hap = new double**[testing_sample_ids->size() + 2];
		for (int i_s = 0; i_s < testing_sample_ids->size(); i_s++)
		{
			imputed_allele_prob_per_sample_per_hap[i_s] = new double*[2];

			imputed_allele_prob_per_sample_per_hap[i_s][0] = new double[2];
			imputed_allele_prob_per_sample_per_hap[i_s][0][0] = ZERO;
			imputed_allele_prob_per_sample_per_hap[i_s][0][1] = ZERO;

			imputed_allele_prob_per_sample_per_hap[i_s][1] = new double[2];
			imputed_allele_prob_per_sample_per_hap[i_s][1][0] = ZERO;
			imputed_allele_prob_per_sample_per_hap[i_s][1][1] = ZERO;
		} // i_s loop.
		new_info[2] = imputed_allele_prob_per_sample_per_hap;

		reference_haplocoded_tag_target_geno_regs->at(i_var)->data = new_info;
	} // i_var loop.

	int n_states = n_ref_haplotypes;
	int n_symbols = 2;

	t_restr_annot_region_list* restr_testing_haplocoded_tag_geno_regs = restructure_annot_regions(testing_haplocoded_tag_geno_regs);
	t_restr_annot_region_list* restr_reference_haplocoded_tag_target_geno_regs = restructure_annot_regions(reference_haplocoded_tag_target_geno_regs, restr_testing_haplocoded_tag_geno_regs->chr_ids);

	for (int i_chr = 0; i_chr < restr_testing_haplocoded_tag_geno_regs->chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Running forward-backward on %s\n", restr_testing_haplocoded_tag_geno_regs->chr_ids->at(i_chr));

		char cur_chr_recombination_rate_fp[1000];
		sprintf(cur_chr_recombination_rate_fp, "%s/%s.map", recombination_rate_dir, restr_testing_haplocoded_tag_geno_regs->chr_ids->at(i_chr));
		vector<t_annot_region*>* cur_chrom_recomb_regs = load_recombination_rates(cur_chr_recombination_rate_fp);
		if (cur_chrom_recomb_regs == NULL)
		{
			fprintf(stderr, "Could not load recombination rates from %s\n", cur_chr_recombination_rate_fp);
			exit(0);
		}

		sort_set_sorting_info(cur_chrom_recomb_regs, sort_regions);

		// Add one element to the beginning of the regions for the current chromosome.
		vector<t_annot_region*>* cur_chr_testing_haplocoded_tag_geno_regs = restr_testing_haplocoded_tag_geno_regs->regions_per_chrom[i_chr];
		vector<t_annot_region*>* cur_chrom_ref_tag_target_var_regs = restr_reference_haplocoded_tag_target_geno_regs->regions_per_chrom[i_chr];

		// Assign the recomb rates.
		fprintf(stderr, "Setting recombination rates using %d recombination mapped regions.\n", cur_chrom_recomb_regs->size());
		for (int i_reg = 0; i_reg < cur_chr_testing_haplocoded_tag_geno_regs->size(); i_reg++)
		{
			double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant_optimized(cur_chr_testing_haplocoded_tag_geno_regs->at(i_reg), cur_chrom_recomb_regs);
			cur_chr_testing_haplocoded_tag_geno_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
		} // i_reg loop.

		for (int i_reg = 0; i_reg < cur_chrom_ref_tag_target_var_regs->size(); i_reg++)
		{
			double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant_optimized(cur_chrom_ref_tag_target_var_regs->at(i_reg), cur_chrom_recomb_regs);
			cur_chrom_ref_tag_target_var_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
		} // i_reg loop.


		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Get the Maximum of the closest target-2-window_centers.
		bool GET_MAX_MIN_TARGET_2_CENTER_DISTANCES = false;
		if (GET_MAX_MIN_TARGET_2_CENTER_DISTANCES)
		{
			double* per_target_min_mid_win_dist = new double[cur_chrom_ref_tag_target_var_regs->size() + 2];
			for (int i_win = 0; i_win < cur_chrom_ref_tag_target_var_regs->size(); i_win++)
			{
				per_target_min_mid_win_dist[i_win] = 1000;
			} // i_win loop.

			for (int win_start_var_i = 0; win_start_var_i < cur_chr_testing_haplocoded_tag_geno_regs->size(); win_start_var_i++)
			{
				double start_cM = cur_chr_testing_haplocoded_tag_geno_regs->at(win_start_var_i)->dbl_score;
				int win_end_var_i = win_start_var_i;
				while (win_end_var_i < cur_chr_testing_haplocoded_tag_geno_regs->size() &&
					fabs(start_cM - cur_chr_testing_haplocoded_tag_geno_regs->at(win_end_var_i)->dbl_score) < tag_2_tag_distance_cM)
				{
					win_end_var_i++;
				} // win_end_var_iloop.

				  // Fix end if passed over the end.
				if (win_end_var_i >= cur_chr_testing_haplocoded_tag_geno_regs->size())
				{
					win_end_var_i = cur_chr_testing_haplocoded_tag_geno_regs->size() - 1;
				} // win_end_var_i loop.

				  // Find the closest targets and update their closest window center parameter.
				double win_center_cM = (cur_chr_testing_haplocoded_tag_geno_regs->at(win_start_var_i)->dbl_score +
					cur_chr_testing_haplocoded_tag_geno_regs->at(win_end_var_i)->dbl_score) / 2;

				for (int i_target = 0; i_target < cur_chrom_ref_tag_target_var_regs->size(); i_target++)
				{
					if (fabs(cur_chrom_ref_tag_target_var_regs->at(i_target)->dbl_score - win_center_cM) < per_target_min_mid_win_dist[i_target])
					{
						per_target_min_mid_win_dist[i_target] = fabs(cur_chrom_ref_tag_target_var_regs->at(i_target)->dbl_score - win_center_cM);
					}
				} // i_target loop.
			} // win_start_var_i loop.

			char dist_fp[1000];
			sprintf(dist_fp, "%s_min_center_2_target_dists.txt", restr_testing_haplocoded_tag_geno_regs->chr_ids->at(i_chr));
			FILE* f_target_dists = open_f(dist_fp, "w");
			double maximum_of_minimum_target_2_center_dist = 0;
			for (int i_target = 0; i_target < cur_chrom_ref_tag_target_var_regs->size(); i_target++)
			{
				fprintf(f_target_dists, "%d\t%.3f\n", cur_chrom_ref_tag_target_var_regs->at(i_target)->start, per_target_min_mid_win_dist[i_target]);
				maximum_of_minimum_target_2_center_dist = MAX(maximum_of_minimum_target_2_center_dist, per_target_min_mid_win_dist[i_target]);
			} // i_target loop.
			fclose(f_target_dists);

			fprintf(stderr, "Maximum of the minimum center-2-target over all variants is %.4f (%.4f)\n", maximum_of_minimum_target_2_center_dist,
				l_target_2_center_pred_buffer_in_cM);
		}
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Get the focus regions on this chromosome.
		int focus_start = cur_chrom_ref_tag_target_var_regs->at(0)->start;
		int focus_end = cur_chrom_ref_tag_target_var_regs->back()->start;
		if (focus_regs != NULL)
		{
			vector<t_annot_region*>* cur_chr_focus_regs = get_regions_per_chromosome(focus_regs, restr_testing_haplocoded_tag_geno_regs->chr_ids->at(i_chr));
			sort(cur_chr_focus_regs->begin(), cur_chr_focus_regs->end(), sort_regions);
			focus_start = cur_chr_focus_regs->at(0)->start;
			focus_end = cur_chr_focus_regs->back()->end;
			delete cur_chr_focus_regs;
		}

		fprintf(stderr, "Focusing on variants @ %d-%d\n", focus_start, focus_end);

		// Find all must be set to false.
		vector<t_annot_region*>* intersects = intersect_annot_regions(cur_chrom_ref_tag_target_var_regs, focus_regs, false);
		vector<t_annot_region*>* cur_chrom_focus_ref_target_var_regs = new vector<t_annot_region*>();
		for (int i_int = 0; i_int < intersects->size(); i_int++)
		{
			t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
			t_annot_region* cur_target_reg = int_info->src_reg;
			cur_chrom_focus_ref_target_var_regs->push_back(cur_target_reg);
		} // i_int loop.

		// Sort the targets in the focus region.
		sort(cur_chrom_focus_ref_target_var_regs->begin(), cur_chrom_focus_ref_target_var_regs->end(), sort_regions);

		vector<t_annot_region*>* prev_win_var_regs = NULL;
		for (int test_sample_i = 0; test_sample_i < testing_sample_ids->size(); test_sample_i++)
		{
			prev_win_var_regs = NULL;

			// Reset all the variant types to untyped.
			for (int i_ref_var = 0; i_ref_var < cur_chrom_focus_ref_target_var_regs->size(); i_ref_var++)
			{
				if (cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->score == VAR_UNTYPED_TARGET_IMPUTED)
				{
					cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->score = VAR_UNTYPED_TARGET;
				}
			} // i_ref_var loop.

			// We keep track of the middle window start coordinates.
			for (int win_mid_start_var_i = 0; win_mid_start_var_i < cur_chr_testing_haplocoded_tag_geno_regs->size(); win_mid_start_var_i++)
			{
				// If the mid-window passed the focus region's end completely, we are done.
				if (cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_start_var_i)->start > focus_end)
				{
					break;
				}

				// We are imputing the variants between current and next variants; select window start and end coordinates.
				double mid_start_cM = cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_start_var_i)->dbl_score;
				int win_mid_end_var_i = win_mid_start_var_i;
				while (win_mid_end_var_i < cur_chr_testing_haplocoded_tag_geno_regs->size() &&
					fabs(mid_start_cM - cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_end_var_i)->dbl_score) < l_target_2_center_pred_buffer_in_cM)
				{
					win_mid_end_var_i++;
				} // win_end_var_i loop.

				if (win_mid_end_var_i >= cur_chr_testing_haplocoded_tag_geno_regs->size())
				{
					win_mid_end_var_i = cur_chr_testing_haplocoded_tag_geno_regs->size() - 1;
				}

				// Get the middle cM, this is the middle of both windows.
				double mid_cM = (cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_start_var_i)->dbl_score +
								cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_end_var_i)->dbl_score) / 2;

				double start_cM = MAX(0, mid_cM - (tag_2_tag_distance_cM / 2));
				int win_start_var_i = win_mid_start_var_i;
				while (win_start_var_i < cur_chr_testing_haplocoded_tag_geno_regs->size() &&
						win_start_var_i > 0 &&
						cur_chr_testing_haplocoded_tag_geno_regs->at(win_start_var_i)->dbl_score > start_cM)
				{
					win_start_var_i--;
				} // win_end_var_i loop.

				// GEt the ending coordiantes for the full window.
				double end_cM = mid_cM + (tag_2_tag_distance_cM / 2);
				int win_end_var_i = win_mid_end_var_i;
				while (win_end_var_i < cur_chr_testing_haplocoded_tag_geno_regs->size() &&
					win_end_var_i > 0 &&
					cur_chr_testing_haplocoded_tag_geno_regs->at(win_end_var_i)->dbl_score < end_cM)
				{
					win_end_var_i++;
				} // win_end_var_i loop.

				  // Fix end if passed over the end.
				if (win_end_var_i >= cur_chr_testing_haplocoded_tag_geno_regs->size())
				{
					win_end_var_i = cur_chr_testing_haplocoded_tag_geno_regs->size() - 1;
				}

				fprintf(stderr, "Sample %d: Prep. win. var_i: %d-(%d-%d)-%d ;; posn: %d-(%d-%d)-%d ;; cM: %.3f-(%.3f-%.3f)-%.3f\n",
					test_sample_i,
					win_start_var_i,
					win_mid_start_var_i,
					win_mid_end_var_i,
					win_end_var_i,
					cur_chr_testing_haplocoded_tag_geno_regs->at(win_start_var_i)->start,
					cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_start_var_i)->start,
					cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_end_var_i)->start,
					cur_chr_testing_haplocoded_tag_geno_regs->at(win_end_var_i)->start,
					cur_chr_testing_haplocoded_tag_geno_regs->at(win_start_var_i)->dbl_score,
					cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_start_var_i)->dbl_score,
					cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_end_var_i)->dbl_score,
					cur_chr_testing_haplocoded_tag_geno_regs->at(win_end_var_i)->dbl_score);

				/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				// Count the number of targets in this window.				
				int n_target_vars_in_cur_win = 0;
				for (int i_ref_var = 0; i_ref_var < cur_chrom_focus_ref_target_var_regs->size(); i_ref_var++)
				{
					double target_cM = cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->dbl_score;
					double left_win_cM = cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_start_var_i)->dbl_score;
					double right_win_cM = cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_end_var_i)->dbl_score;

					// Make sure this is an untyped variant and it is covered well by tags.
					if (cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->score == VAR_UNTYPED_TARGET_IMPUTED ||
						cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->score == VAR_TYPED_TAG ||
						target_cM < left_win_cM ||
						target_cM > right_win_cM )
					{
						
					}
					else
					{
						n_target_vars_in_cur_win++;
					}

					if (cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->start > cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_end_var_i)->end)
					{
						break;
					}
				} // i_ref_var loop.

				if (n_target_vars_in_cur_win == 0)
				{
					if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
					{
						fprintf(stderr, "Skipping window %d-%d\n",
							cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_start_var_i)->start,
							cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_end_var_i)->start);
					}

					continue;
				}
				else
				{
					if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
					{
						fprintf(stderr, "%d target variants in window\n",
							n_target_vars_in_cur_win);
					}
				}
				/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

				start_cM = cur_chr_testing_haplocoded_tag_geno_regs->at(win_start_var_i)->dbl_score;
				end_cM = cur_chr_testing_haplocoded_tag_geno_regs->at(win_end_var_i)->dbl_score;
				
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				// Set the current window's tag variant regions.
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				vector<t_annot_region*>* cur_win_var_regs = new vector<t_annot_region*>();
				
				// Subselect windows if requested.
				if (max_n_tag_vars_per_window == 0 ||
					(win_end_var_i - win_start_var_i + 1) < max_n_tag_vars_per_window)
				{
					cur_win_var_regs->insert(cur_win_var_regs->end(),
						cur_chr_testing_haplocoded_tag_geno_regs->begin() + win_start_var_i,
						cur_chr_testing_haplocoded_tag_geno_regs->begin() + win_end_var_i + 1);
				}
				else
				{
					cur_win_var_regs->push_back(cur_chr_testing_haplocoded_tag_geno_regs->at(win_start_var_i));
					cur_win_var_regs->push_back(cur_chr_testing_haplocoded_tag_geno_regs->at(win_end_var_i));

					double n_vars_2_fit = max_n_tag_vars_per_window - 2;
					double n_vars_in_window = win_end_var_i - win_start_var_i + 1 - 2;
					double frac_i_update = n_vars_in_window / n_vars_2_fit;
					double frac_i = win_start_var_i + 1;

					if (frac_i_update < 1)
					{
						fprintf(stderr, "Sanity check failed: frac_i_update < 1: %lf\n", frac_i_update);
						exit(0);
					}

					int last_var_i = -1;
					while(1)
					{
						int var_i = (int)(floor(frac_i));

						if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
						{
							fprintf(stderr, "Selecting var %d (%.4f) in window @ (%d-%d) [%d windows already]\n", var_i, frac_i, win_start_var_i, win_end_var_i,
								cur_win_var_regs->size());
						}

						if (var_i >= win_end_var_i - 1)
						{
							break;
						}

						if (last_var_i == var_i)
						{
							fprintf(stderr, "Sanity check failed while subsampling windows: %d is repeated.\n");
							exit(0);
						}

						cur_win_var_regs->push_back(cur_chr_testing_haplocoded_tag_geno_regs->at(var_i));
						last_var_i = var_i;

						frac_i += frac_i_update;
					} // frac_i loop.

					if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
					{
						fprintf(stderr, "Subsampled %d/%d::%d tag variants in the window.\n",
							win_end_var_i - win_start_var_i + 1,
							cur_win_var_regs->size(), max_n_tag_vars_per_window);
					}

					sort(cur_win_var_regs->begin(), cur_win_var_regs->end(), sort_regions);
				} // var subsampling check.			

				  // Make sure we inserted the correct elements.
				if (!t_string::compare_strings(cur_win_var_regs->at(0)->name, cur_chr_testing_haplocoded_tag_geno_regs->at(win_start_var_i)->name) ||
					!t_string::compare_strings(cur_win_var_regs->back()->name, cur_chr_testing_haplocoded_tag_geno_regs->at(win_end_var_i)->name))
				{
					fprintf(stderr, "Start and end did not copy as expected.\n");
					exit(0);
				}

				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
				{
					fprintf(stderr, "Computing fore-ward scores for sample %d\n", test_sample_i);
				}

				// Generate the blocks and the reduced states: These must be generated before ends are added.
				vector<t_var_block*>* haplo_var_blocks = generate_reduced_state_blocks_constant_size_blocks(cur_win_var_regs,
					ref_sample_ids,
					1,
					l_blocks);

				fprintf(stderr, "Window (%d/%d/%d tags)::%d targets: %s:%d-%d (%d-%d) spanning %.4f cM between %.4f-%.4f cMs [%d haplo-blocks of %d variants]\n",
						cur_win_var_regs->size(), max_n_tag_vars_per_window, win_end_var_i - win_start_var_i + 1,
						n_target_vars_in_cur_win,
						cur_win_var_regs->at(0)->chrom, cur_win_var_regs->at(0)->start, cur_win_var_regs->back()->start,
						win_start_var_i, win_end_var_i,
						fabs(start_cM - end_cM),
						start_cM, end_cM,
						haplo_var_blocks->size(), l_blocks);

				if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
				{
					fprintf(stderr, "%d haplotype blocks:\n", haplo_var_blocks->size());
					for (int i_block = 0; i_block < haplo_var_blocks->size(); i_block++)
					{
						fprintf(stderr, "%d-%d: %d haplotypes\n",
							haplo_var_blocks->at(i_block)->start_var_i, haplo_var_blocks->at(i_block)->end_var_i,
							haplo_var_blocks->at(i_block)->haplogroup_haplo_indices->size());
					} // i_block loop.
				}

				// Set the start and end state regions.
				t_annot_region* start_state_reg = NULL;
				cur_win_var_regs->insert(cur_win_var_regs->begin(), start_state_reg);

				// Set the number of variants for indexing purposes here.
				int n_vars = cur_win_var_regs->size() - 1;

				// Add the end state region.
				t_annot_region* end_state_reg = NULL;
				cur_win_var_regs->push_back(end_state_reg);

				// Allocate and initialize the forward/backward arrays.
				double*** fore_scores_per_hap = new double**[n_vars + 2];
				double*** back_scores_per_hap = new double**[n_vars + 2];
				for (int var_i = 0; var_i <= n_vars + 1; var_i++)
				{
					fore_scores_per_hap[var_i] = new double*[2];
					back_scores_per_hap[var_i] = new double*[2];

					for (int test_hap = 0; test_hap < 2; test_hap++)
					{
						fore_scores_per_hap[var_i][test_hap] = new double[n_states + 2];
						memset(fore_scores_per_hap[var_i][test_hap], 0, sizeof(double) * (n_states + 1));

						back_scores_per_hap[var_i][test_hap] = new double[n_states + 2];
						memset(back_scores_per_hap[var_i][test_hap], 0, sizeof(double) * (n_states + 2));
					}
				} // i loop.

				// Initialize the state probabilities for both haplotypes.
				for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
				{
					fore_scores_per_hap[0][0][state_i] = XCONVERT_LIN((double)1.0 / n_ref_haplotypes);
					fore_scores_per_hap[0][1][state_i] = XCONVERT_LIN((double)1.0 / n_ref_haplotypes);
				} // state_i loop.

				for (int test_hap_i = 0; test_hap_i < 2; test_hap_i++)
				{
					// Process each block sequentially.
					for (int block_i = 0; block_i < haplo_var_blocks->size(); block_i++)
					{
						// Process the current block.
						int cur_block_start_var_i = haplo_var_blocks->at(block_i)->start_var_i;
						int cur_block_end_var_i = haplo_var_blocks->at(block_i)->end_var_i;

						int n_cur_block_haplogrps = haplo_var_blocks->at(block_i)->haplogroup_haplo_indices->size();
						vector<vector<int>*>* per_hplgrp_hap_i = haplo_var_blocks->at(block_i)->haplogroup_haplo_indices;

						if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
						{
							fprintf(stderr, "Processing test_hap_i:%d;; Block:%d-%d with %d haplogroups.\n",
								test_hap_i,
								cur_block_start_var_i, cur_block_end_var_i,
								n_cur_block_haplogrps);
						}

						////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
						// The first haplotype in each haplogroup is the representative. No specific reason.
						vector<int>* per_hplgrp_representative_hap_i = new vector<int>();
						for (int i_grp = 0; i_grp < n_cur_block_haplogrps; i_grp++)
						{
							int cur_grp_representative_hap_i = per_hplgrp_hap_i->at(i_grp)->at(0);

							per_hplgrp_representative_hap_i->push_back(cur_grp_representative_hap_i);

							if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
							{
								fprintf(stderr, "Block %d::Haplogroup %d; representative_hap_i: %d (%.4f)\n",
									block_i,
									i_grp,
									cur_grp_representative_hap_i);
							}
						} // i_grp loop.

						////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
						// 2. Compute the total other transition for all groups in this block.
						for (int var_i = cur_block_start_var_i;
							var_i <= cur_block_end_var_i;
							var_i++)
						{
							if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
							{
								if (var_i % 10 == 0)
								{
									fprintf(stderr, "Forward: sample_i: %d/%d: var_i: %d         \r", test_sample_i, testing_sample_ids->size(), var_i);
								}
							}

							void** cur_tag_var_info = NULL;
							char* test_sample_geno = NULL;
							char* ref_sample_geno = NULL;

							if (var_i <= n_vars)
							{
								void** cur_tag_var_info = (void**)(cur_win_var_regs->at(var_i)->data);
								test_sample_geno = (char*)(cur_tag_var_info[0]);
								ref_sample_geno = (char*)(cur_tag_var_info[1]);
							}

							/////////////////////////////////////////////////////
							// Pre-compute the transition probabilities:
							double prev_var_cM = cur_win_var_regs->at(1)->dbl_score;
							if (var_i > 1)
							{
								prev_var_cM = cur_win_var_regs->at(var_i - 1)->dbl_score;
							}

							double cur_var_cM = cur_win_var_regs->at(n_vars)->dbl_score;
							if (var_i <= n_vars)
							{
								cur_var_cM = cur_win_var_regs->at(var_i)->dbl_score;
							}

							double r_m = fabs(cur_var_cM - prev_var_cM);
							double rho_m = 4 * N_e * r_m;

							double tau_m = 1 - exp(-1 * rho_m / n_ref_haplotypes);

							double other_prob = tau_m / n_ref_haplotypes;
							double self_prob = (1 - tau_m) + (tau_m / n_ref_haplotypes);

							//fprintf(stderr, "Self prob: %.5f, other_prob: %.5f (%.4f, %.4f)\n", self_prob, other_prob, prev_var_cM, cur_var_cM);
							/////////////////////////////////////////////////////
							// Compute the total other-trans probabilities for each haplogroup.
							vector<double>* per_hpl_grp_other_trans_total = new vector<double>();
							for (int i_grp = 0; i_grp < n_cur_block_haplogrps; i_grp++)
							{
								// The transition and emission probabilities are constant for each haplogroup.
								double cur_grp_total_other_trans_prob = ZERO;

								// Transition prob: Other trans.
								double trans_prob = other_prob;

								for (int i_hap_i = 0;
									i_hap_i < per_hplgrp_hap_i->at(i_grp)->size();
									i_hap_i++)
								{
									int cur_grp_state_i = per_hplgrp_hap_i->at(i_grp)->at(i_hap_i);

									cur_grp_total_other_trans_prob = XSUM(cur_grp_total_other_trans_prob,
																		XMUL(XCONVERT_LIN(trans_prob), fore_scores_per_hap[var_i - 1][test_hap_i][cur_grp_state_i]));
								} // i_sample_i loop.

								per_hpl_grp_other_trans_total->push_back(cur_grp_total_other_trans_prob);

								if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
								{
									fprintf(stderr, "Block %d::Haplogroup %d::Total Other Transition: %.4f\n",
										block_i,
										i_grp,
										cur_grp_total_other_trans_prob);
								}
							} // i_grp loop.

							// Loop over all the states.
							for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
							{
								// Recurse over the previous states.
								fore_scores_per_hap[var_i][test_hap_i][state_i] = ZERO;

								//////////////////////////////////////////////////////////////////////////////////////////////////////////
								//////////////////////////////////////////////////////////////////////////////////////////////////////////
								// Do a self transition with a delta probability.
								double self_min_trans_prob = self_prob - other_prob;

								int ref_sample_i = state_i / 2;
								int ref_hap_i = state_i % 2;
								int cur_ref_allele = get_allele_per_haplotype(ref_sample_geno[ref_sample_i], ref_hap_i);
								int cur_test_allele = get_allele_per_haplotype(test_sample_geno[test_sample_i], test_hap_i);

								double emit_prob_per_ref_allele[2];
								emit_prob_per_ref_allele[cur_ref_allele] = 1 - allele_err;
								emit_prob_per_ref_allele[1 - cur_ref_allele] = allele_err;
								double emit_prob = emit_prob_per_ref_allele[cur_test_allele];

								double trans_emit_prob = XMUL(XCONVERT_LIN(self_min_trans_prob), XCONVERT_LIN(emit_prob));
								trans_emit_prob = XMUL(trans_emit_prob, XCONVERT_LIN(global_lin_scaler));

								fore_scores_per_hap[var_i][test_hap_i][state_i] = XSUM(fore_scores_per_hap[var_i][test_hap_i][state_i],
																				XMUL(trans_emit_prob, fore_scores_per_hap[var_i - 1][test_hap_i][state_i]));

								//////////////////////////////////////////////////////////////////////////////////////////////////////////
								//////////////////////////////////////////////////////////////////////////////////////////////////////////
								// Now do a transition to all the top scoring haplotypes.
								for (int hpl_grp_i = 0; hpl_grp_i < n_cur_block_haplogrps; hpl_grp_i++)
								{
									double cur_haplogrp_total_other_trans = per_hpl_grp_other_trans_total->at(hpl_grp_i);

									// Apply the scaler.
									double emit_prob = emit_prob_per_ref_allele[cur_test_allele];
									double trans_emit_prob = XMUL(XCONVERT_LIN(1.0), XCONVERT_LIN(emit_prob));
									trans_emit_prob = XMUL(trans_emit_prob, XCONVERT_LIN(global_lin_scaler));

									fore_scores_per_hap[var_i][test_hap_i][state_i] = XSUM(fore_scores_per_hap[var_i][test_hap_i][state_i], 
																							XMUL(trans_emit_prob, cur_haplogrp_total_other_trans));
								} // hpl_grp_i loop.
							} // state_i loop.

							delete per_hpl_grp_other_trans_total;
						} // var_i loop.

						delete per_hplgrp_representative_hap_i;
					} // block_i loop.
				} // test_hap_i loop.

				// Move the backward scores.
				if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
				{
					fprintf(stderr, "Computing back-ward scores for sample %d\n", test_sample_i);
				}

				for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
				{
					back_scores_per_hap[n_vars + 1][0][state_i] = XCONVERT_LIN((double)1.0 / n_ref_haplotypes);
					back_scores_per_hap[n_vars + 1][1][state_i] = XCONVERT_LIN((double)1.0 / n_ref_haplotypes);
				} // state_i loop.

				for (int test_hap_i = 0; test_hap_i < 2; test_hap_i++)
				{
					// Process each block sequentially.
					for (int block_i = haplo_var_blocks->size()-1; block_i >= 0; block_i--)
					{
						// Process the current block.
						int cur_block_start_var_i = haplo_var_blocks->at(block_i)->start_var_i;
						int cur_block_end_var_i = haplo_var_blocks->at(block_i)->end_var_i;

						int n_cur_block_haplogrps = haplo_var_blocks->at(block_i)->haplogroup_haplo_indices->size();
						vector<vector<int>*>* per_hplgrp_hap_i = haplo_var_blocks->at(block_i)->haplogroup_haplo_indices;

						if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
						{
							fprintf(stderr, "Processing test_hap_i:%d;; Block:%d-%d with %d haplogroups.\n",
								test_hap_i,
								cur_block_start_var_i, cur_block_end_var_i,
								n_cur_block_haplogrps);
						}

						////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
						// The first haplotype in each haplogroup is the representative. No specific reason.
						vector<int>* per_hplgrp_representative_hap_i = new vector<int>();
						for (int i_grp = 0; i_grp < n_cur_block_haplogrps; i_grp++)
						{
							int cur_grp_representative_hap_i = per_hplgrp_hap_i->at(i_grp)->at(0);

							per_hplgrp_representative_hap_i->push_back(cur_grp_representative_hap_i);

							if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
							{
								fprintf(stderr, "Block %d::Haplogroup %d; representative_hap_i: %d (%.4f)\n",
									block_i,
									i_grp,
									cur_grp_representative_hap_i);
							}
						} // i_grp loop.

						////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
						// 2. Compute the total other transition for all groups in this block.
						for (int var_i = cur_block_end_var_i;
							var_i >= cur_block_start_var_i;
							var_i--)
						{
							if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
							{
								if (var_i % 10 == 0)
								{
									fprintf(stderr, "Backward: sample_i: %d/%d: var_i: %d         \r", test_sample_i, testing_sample_ids->size(), var_i);
								}
							}

							// Do not set the genotypes for the after-last index.
							char* test_sample_geno = NULL;
							char* ref_sample_geno = NULL;
							if (var_i > 0)
							{
								void** cur_tag_var_info = (void**)(cur_win_var_regs->at(var_i)->data);
								test_sample_geno = (char*)(cur_tag_var_info[0]);
								ref_sample_geno = (char*)(cur_tag_var_info[1]);
							}

							// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
							// Set the transition probabilities.
							double cur_cM = cur_win_var_regs->at(1)->dbl_score;
							if (var_i > 0)
							{
								cur_cM = cur_win_var_regs->at(var_i)->dbl_score;
							}

							double next_cM = cur_win_var_regs->at(n_vars)->dbl_score;
							if (var_i <= n_vars - 1)
							{
								next_cM = cur_win_var_regs->at(var_i + 1)->dbl_score;
							}

							double r_m = fabs(cur_cM - next_cM);
							double rho_m = 4 * N_e * r_m;

							double tau_m = 1 - exp(-1 * rho_m / n_ref_haplotypes);

							double other_prob = tau_m / n_ref_haplotypes;
							double self_prob = (1 - tau_m) + (tau_m / n_ref_haplotypes);

							// Go over all the next states and compute total transitions.
							double total_next_other_transitions = ZERO;
							for (int next_state_i = 0; next_state_i < n_ref_haplotypes; next_state_i++)
							{
								total_next_other_transitions = XSUM(total_next_other_transitions, XMUL(other_prob, back_scores_per_hap[var_i+1][test_hap_i][next_state_i]));
							} // state_i loop.

							// Go over all the reference haplotypes.
							for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
							{
								back_scores_per_hap[var_i][test_hap_i][state_i] = ZERO;

								////////////////////////////////////////////////////////////////////////////////////////////////////
								// Do the self transition.
								double self_min_trans_prob = self_prob - other_prob;

								int ref_sample_i = state_i / 2;
								int ref_hap_i = state_i % 2;
								int cur_ref_allele = get_allele_per_haplotype(ref_sample_geno[ref_sample_i], ref_hap_i);
								int cur_test_allele = get_allele_per_haplotype(test_sample_geno[test_sample_i], test_hap_i);

								double emit_prob_per_ref_allele[2];
								emit_prob_per_ref_allele[cur_ref_allele] = 1 - allele_err;
								emit_prob_per_ref_allele[1 - cur_ref_allele] = allele_err;

								double emit_prob = emit_prob_per_ref_allele[cur_test_allele];
								double trans_emit_prob = XMUL(XCONVERT_LIN(emit_prob), XCONVERT_LIN(self_min_trans_prob));
								trans_emit_prob = XMUL(trans_emit_prob, XCONVERT_LIN(global_lin_scaler));

								back_scores_per_hap[var_i][test_hap_i][state_i] = XSUM(back_scores_per_hap[var_i][test_hap_i][state_i],
																					XMUL(back_scores_per_hap[var_i + 1][test_hap_i][state_i], trans_emit_prob));

								////////////////////////////////////////////////////////////////////////////////////////////////////
								// Add the aggregate other transitions w/o emits.
								emit_prob = emit_prob_per_ref_allele[cur_test_allele];
								trans_emit_prob = XMUL(XCONVERT_LIN(1.0), XCONVERT_LIN(emit_prob));
								trans_emit_prob = XMUL(trans_emit_prob, XCONVERT_LIN(global_lin_scaler));
								back_scores_per_hap[var_i][test_hap_i][state_i] = XSUM(back_scores_per_hap[var_i][test_hap_i][state_i], 
																						XMUL(trans_emit_prob, total_next_other_transitions));

							} // state_i loop.
						} // var_i loop.

						delete per_hplgrp_representative_hap_i;
					} // block_i loop.
				} // test_hap_i loop.

				// Compute the total log forward and backward probabilities.
				double per_hap_total_log_fore_prob[2];
				per_hap_total_log_fore_prob[0] = XCONVERT_LIN(0.0);
				per_hap_total_log_fore_prob[1] = XCONVERT_LIN(0.0);
				double per_hap_total_log_back_prob[2];
				per_hap_total_log_back_prob[0] = XCONVERT_LIN(0.0);
				per_hap_total_log_back_prob[1] = XCONVERT_LIN(0.0);

				for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
				{
					per_hap_total_log_fore_prob[0] = XSUM(per_hap_total_log_fore_prob[0], fore_scores_per_hap[n_vars][0][state_i]);
					per_hap_total_log_back_prob[0] = XSUM(per_hap_total_log_back_prob[0], back_scores_per_hap[1][0][state_i]);

					per_hap_total_log_fore_prob[1] = XSUM(per_hap_total_log_fore_prob[1], fore_scores_per_hap[n_vars][1][state_i]);
					per_hap_total_log_back_prob[1] = XSUM(per_hap_total_log_back_prob[1], back_scores_per_hap[1][1][state_i]);
				} // state_i loop.

				fprintf(stderr, "Sample %d: Haplotype 0 total probabilities: fore=%.5f ;; back=%.5f\n", test_sample_i, XCONVERT_2_LOG(per_hap_total_log_fore_prob[0]), XCONVERT_2_LOG(per_hap_total_log_back_prob[0]));
				fprintf(stderr, "Sample %d: Haplotype 1 total probabilities: fore=%.5f ;; back=%.5f\n", test_sample_i, XCONVERT_2_LOG(per_hap_total_log_fore_prob[1]), XCONVERT_2_LOG(per_hap_total_log_back_prob[1]));

				// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
				// Start computing the posterior probabilities for the untyped markers.
				// Now go over the untyped genotypes and compute the posterior probability of all the untyped markers.				
				int n_untyped = 0;

				if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
				{
					fprintf(stderr, "Computing the probabilities of the %d untyped markers in %d-%d.\n",
						cur_chrom_focus_ref_target_var_regs->size(),
						n_target_vars_in_cur_win,
						win_mid_start_var_i, win_mid_end_var_i);
				}

				//// +1 is necessary to match the 1-based coordinates.
				//int rel_win_mid_start_i = win_mid_start_var_i - win_start_var_i + 1;
				//int rel_win_mid_end_i = win_mid_end_var_i - win_start_var_i + 1;

				//// Relative coordinates must match exactly to the absolute coordinates.
				//if (cur_win_var_regs->at(rel_win_mid_end_i)->dbl_score != cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_end_var_i)->dbl_score ||
				//	cur_win_var_regs->at(rel_win_mid_start_i)->dbl_score != cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_start_var_i)->dbl_score)
				//{
				//	fprintf(stderr, "Sanuty check failed: Start: %lf, %lf; End: %lf, %lf\n", 
				//			cur_win_var_regs->at(rel_win_mid_start_i)->dbl_score, cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_start_var_i)->dbl_score,
				//			cur_win_var_regs->at(rel_win_mid_end_i)->dbl_score, cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_end_var_i)->dbl_score);
				//	exit(0);
				//}

				if (win_mid_start_var_i > win_end_var_i ||
					win_mid_end_var_i > win_end_var_i)
				{
					fprintf(stderr, "Sanity check failed: %d-(%d-%d)-%d\n", win_start_var_i, win_mid_start_var_i, win_mid_end_var_i, win_end_var_i);
					exit(0);
				}

				// Go over all the reference target variants in the focus region.
				for (int i_ref_var = 0; i_ref_var < cur_chrom_focus_ref_target_var_regs->size(); i_ref_var++)
				{
					double target_cM = cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->dbl_score;
					//double left_mid_cM = cur_win_var_regs->at(rel_win_mid_start_i)->dbl_score;
					//double right_mid_cM = cur_win_var_regs->at(rel_win_mid_end_i)->dbl_score;

					double left_mid_cM = cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_start_var_i)->dbl_score;
					double right_mid_cM = cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_end_var_i)->dbl_score;

					// Make sure this is an untyped variant and it is covered well by tags.
					if (cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->score == VAR_UNTYPED_TARGET_IMPUTED ||
						cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->score == VAR_TYPED_TAG ||
						target_cM < left_mid_cM ||
						target_cM > right_mid_cM)
					{
						if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
						{
							fprintf(stderr, "Skipping %s:%d @ cM: %.4f (%.4f-%.4f) typed: %d (%.4f, %.4f)\n",
								cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->chrom, cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->start,
								target_cM, left_mid_cM, right_mid_cM,
								cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->score,
								fabs(right_mid_cM - target_cM),
								fabs(left_mid_cM - target_cM));
						}
						continue;
					}

					// This is the untyped variant's information.
					void** cur_untyped_var_ref_var_info = (void**)(cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->data);
					char* untyped_var_ref_geno = (char*)(cur_untyped_var_ref_var_info[0]);
					char* known_var_ref_geno = (char*)(cur_untyped_var_ref_var_info[1]);
					double*** imputed_allele_prob_per_sample_per_hap = (double***)(cur_untyped_var_ref_var_info[2]);

					// If there is probability assigned to this variant for this sample, do not process it.
					if (imputed_allele_prob_per_sample_per_hap[test_sample_i][0][0] != ZERO ||
						imputed_allele_prob_per_sample_per_hap[test_sample_i][0][1] != ZERO)
					{
						continue;
					}

					n_untyped++;

					//fprintf(stderr, "@%d. untyped marker         \r", n_untyped);

					// Find the closest tag variants to this target: This is the fore-back coordinates among the middle window coordinates.
					//for (int var_i = rel_win_mid_start_i; var_i < rel_win_mid_end_i; var_i++)
					for (int var_i = 1; var_i < n_vars; var_i++)
					{
						// Make sure the tag variants span the middle region.
						if (cur_win_var_regs->at(var_i+1)->dbl_score < left_mid_cM ||
							cur_win_var_regs->at(var_i)->dbl_score > right_mid_cM)
						{
							continue;
						}

						// Find the untyped target variant in between the current tags.
						if (cur_win_var_regs->at(var_i)->start < cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->start &&
							cur_win_var_regs->at(var_i + 1)->start > cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->start)
						{
							if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
							{
								fprintf(stderr, "Found target %d (%d): %d-%d-%d\n",
									i_ref_var,
									cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->score,
									cur_win_var_regs->at(var_i)->start,
									cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->start,
									cur_win_var_regs->at(var_i + 1)->start);
							}

							void** cur_tag_var_info = (void**)(cur_win_var_regs->at(var_i)->data);
							char* test_sample_geno = (char*)(cur_tag_var_info[0]);
							char* ref_sample_geno = (char*)(cur_tag_var_info[1]);

							// Get the self transition probability.
							double target_cM = cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->dbl_score;
							double cur_cM = cur_win_var_regs->at(var_i)->dbl_score;
							double next_cM = cur_cM;
							if (var_i <= n_vars - 1)
							{
								next_cM = cur_win_var_regs->at(var_i + 1)->dbl_score;
							}
							double r_m = fabs(cur_cM - next_cM);
							double rho_m = 4 * N_e * r_m;

							double tau_m = 1 - exp(-1 * rho_m / n_ref_haplotypes);

							double other_prob = tau_m / n_ref_haplotypes;
							double self_prob = (1 - tau_m) + (tau_m / n_ref_haplotypes);
							//////////////////////////////////////////////////////////////////////////////////////////////////

							for (int hap_i = 0; hap_i < 2; hap_i++)
							{
								double per_allele_probs[2];
								per_allele_probs[0] = ZERO;
								per_allele_probs[1] = ZERO;

								if (posterior_mode == POSTERIOR_MODE_SINGLE_PATH)
								{
									for (int untyped_state_i = 0; untyped_state_i < n_ref_haplotypes; untyped_state_i++)
									{
										// For this untyped state in the refernce, we would like to sum all the paths that pass through here.
										int ref_sample_i = untyped_state_i / 2;
										int ref_hap_i = untyped_state_i % 2;
										int cur_untyped_ref_allele = get_allele_per_haplotype(untyped_var_ref_geno[ref_sample_i], ref_hap_i);

										// Single path: No changes in state.
										double cur_untyped_state_fore_prob = XMUL(self_prob, fore_scores_per_hap[var_i][hap_i][untyped_state_i]);

										double cur_untyped_state_back_prob = back_scores_per_hap[var_i + 1][hap_i][untyped_state_i];

										per_allele_probs[cur_untyped_ref_allele] = XSUM(per_allele_probs[cur_untyped_ref_allele], XMUL(cur_untyped_state_fore_prob, cur_untyped_state_back_prob));
									} // untyped_state_i loop.
								} // single path posteriors.
								else if (posterior_mode == POSTERIOR_MODE_WEIGHTED_PW_PATH)
								{
									for (int untyped_state_i = 0; untyped_state_i < n_ref_haplotypes; untyped_state_i++)
									{
										for (int untyped_state_j = 0; untyped_state_j < n_ref_haplotypes; untyped_state_j++)
										{
											// For this untyped state in the refernce, we would like to sum all the paths that pass through here.
											int ref_sample_i = untyped_state_i / 2;
											int ref_hap_i = untyped_state_i % 2;
											int cur_untyped_ref_allele_i = get_allele_per_haplotype(untyped_var_ref_geno[ref_sample_i], ref_hap_i);

											double allele_i_weight = fabs(next_cM - target_cM) / fabs(next_cM - cur_cM);
											
											int ref_sample_j = untyped_state_j / 2;
											int ref_hap_j = untyped_state_j % 2;
											int cur_untyped_ref_allele_j = get_allele_per_haplotype(untyped_var_ref_geno[ref_sample_j], ref_hap_j);

											double allele_j_weight = 1 - allele_i_weight;

											// Convert.
											allele_i_weight = XCONVERT_LIN(allele_i_weight);
											allele_j_weight = XCONVERT_LIN(allele_j_weight);

											double cur_untyped_state_fore_prob = XMUL(other_prob, fore_scores_per_hap[var_i][hap_i][untyped_state_i]);
											if (untyped_state_i == untyped_state_j)
											{
												cur_untyped_state_fore_prob = XMUL(self_prob, fore_scores_per_hap[var_i][hap_i][untyped_state_i]);
											}

											double cur_untyped_state_back_prob = back_scores_per_hap[var_i + 1][hap_i][untyped_state_j];

											per_allele_probs[cur_untyped_ref_allele_i] = XSUM(per_allele_probs[cur_untyped_ref_allele_i], 
																								XMUL(XMUL(allele_i_weight, cur_untyped_state_fore_prob), cur_untyped_state_back_prob));
											per_allele_probs[cur_untyped_ref_allele_j] = XSUM(per_allele_probs[cur_untyped_ref_allele_j], 
																								XMUL(XMUL(allele_j_weight, cur_untyped_state_fore_prob), cur_untyped_state_back_prob));
										} // untyped_state_j loop/block.
									} // untyped_state_i loop.
								} // weighted pairwise path posterior.

								// We went through all the states, normalize the probabilities.
								per_allele_probs[0] = XDIV(per_allele_probs[0], per_hap_total_log_fore_prob[hap_i]);
								per_allele_probs[1] = XDIV(per_allele_probs[1], per_hap_total_log_fore_prob[hap_i]);

								int known_geno = -1;
								if (known_var_ref_geno)
								{
									known_geno = (int)(known_var_ref_geno[test_sample_i]);
								}

								if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
								{
									fprintf(stderr, "Var %d on sample %d: %s:%d (%s): Hap%d: [0:%.4f, 1:%.4f]; (known_geno=%d)\n",
										i_ref_var,
										test_sample_i,
										cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->chrom,
										cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->start,
										cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->name,
										hap_i,
										XCONVERT_2_LOG(per_allele_probs[0]),
										XCONVERT_2_LOG(per_allele_probs[1]),
										known_geno);
								}

								imputed_allele_prob_per_sample_per_hap[test_sample_i][hap_i][0] = per_allele_probs[0];
								imputed_allele_prob_per_sample_per_hap[test_sample_i][hap_i][1] = per_allele_probs[1];
							} // hap_i loop.

							// Set the inputated variant type.
							if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
							{
								fprintf(stderr, "Setting target ref_var_i=%d (%d @ %.4f)\n", 
										i_ref_var, 
										cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->start,
										cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->dbl_score);
							}

							// Set the inputated variant type.
							cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->score = VAR_UNTYPED_TARGET_IMPUTED;

							break; // Breaks out of position search.
						} // middle check.
					} // i_var loop.
				} // i_ref_var loop.

				// Free haplogroup and block memories.
				for (int i_block = 0; i_block < haplo_var_blocks->size(); i_block++)
				{
					for (int i_grp = 0; i_grp < haplo_var_blocks->at(i_block)->haplogroup_haplo_indices->size(); i_grp++)
					{
						delete haplo_var_blocks->at(i_block)->haplogroup_haplo_indices->at(i_grp);
					}
				} // i_block loop.

				// Allocate and initialize the forward/backward arrays.
				for (int var_i = 0; var_i <= n_vars + 1; var_i++)
				{
					for (int test_hap = 0; test_hap < 2; test_hap++)
					{
						delete[] fore_scores_per_hap[var_i][test_hap];
						delete[] back_scores_per_hap[var_i][test_hap];
					} // test_hap loop.

					delete[] back_scores_per_hap[var_i];
					delete[] fore_scores_per_hap[var_i];
				} // i loop.

				delete[] fore_scores_per_hap;
				delete[] back_scores_per_hap;

				// Copy the current windows.
				if (prev_win_var_regs != NULL)
				{
					delete prev_win_var_regs;
				}

				prev_win_var_regs = cur_win_var_regs;
			} // win_mid_start_var_i loop.
		} // test_sample_i loop.
	} // i_chr loop.

	  // Save the imputed genotype probabilities.
	FILE* f_geno_probs = open_f(geno_probs_op_fp, "w");
	for (int i_ref_var = 0; i_ref_var < reference_haplocoded_tag_target_geno_regs->size(); i_ref_var++)
	{
		// Make sure this is an untyped variant and it is covered well by tags.
		if (reference_haplocoded_tag_target_geno_regs->at(i_ref_var)->score != VAR_UNTYPED_TARGET_IMPUTED)
		{
			continue;
		}

		// This is the untyped variant's information.
		void** cur_untyped_var_ref_var_info = (void**)(reference_haplocoded_tag_target_geno_regs->at(i_ref_var)->data);
		double*** imputed_allele_prob_per_sample_per_hap = (double***)(cur_untyped_var_ref_var_info[2]);

		fprintf(f_geno_probs, "%s\t%d\t%d\t%s", reference_haplocoded_tag_target_geno_regs->at(i_ref_var)->chrom,
			translate_coord(reference_haplocoded_tag_target_geno_regs->at(i_ref_var)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(reference_haplocoded_tag_target_geno_regs->at(i_ref_var)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			reference_haplocoded_tag_target_geno_regs->at(i_ref_var)->name);

		for (int i_s = 0; i_s < testing_sample_ids->size(); i_s++)
		{
			fprintf(f_geno_probs, "\t%.3f;%.3f;%.3f;%.3f",
				XCONVERT_2_LOG(imputed_allele_prob_per_sample_per_hap[i_s][0][0]),
				XCONVERT_2_LOG(imputed_allele_prob_per_sample_per_hap[i_s][0][1]),
				XCONVERT_2_LOG(imputed_allele_prob_per_sample_per_hap[i_s][1][0]),
				XCONVERT_2_LOG(imputed_allele_prob_per_sample_per_hap[i_s][1][1]));
		} // i_s loop.

		fprintf(f_geno_probs, "\n");
	} // i_ref_var loop.
	close_f(f_geno_probs, geno_probs_op_fp);

	fprintf(stderr, "Done!\n");
}

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
	char* geno_probs_op_fp)
{
	fprintf(stderr, "Running Viterbi state reduction using blocks of length %d\n", l_blocks);

	// Make sure the target-2-center is not larger than the window length
	// THIS SHOULD BE TURNED BACK ON.
	if ((2 * l_target_2_center_pred_buffer_in_cM) > tag_2_tag_distance_cM)
	{
		//l_target_2_center_pred_buffer_in_cM = tag_2_tag_distance_cM / 2;
	}

	double allele_err = pow(10, -4);
	double frac_2_max_2_prune = pow(10, 10);

	double(*XSUM)(double num1, double num2) = NULL;
	double(*XMUL)(double num1, double num2) = NULL;
	double(*XDIV)(double num1, double num2) = NULL;
	double(*XSUB)(double num1, double num2) = NULL;
	bool(*XCOMPARE)(double num1, double num2) = NULL;
	double(*XCONVERT_LIN)(double num1) = NULL;
	double(*XCONVERT_LOG)(double num1) = NULL;
	double(*XCONVERT_2_LOG)(double num1) = NULL;
	double(*XCONVERT_2_LIN)(double num1) = NULL;
	double ZERO = 0.0;

	if (math_mode == MATH_MODE_LOG_SPACE)
	{
		fprintf(stderr, "Math mode log.\n");
		XSUM = xlog_sum;
		XMUL = xlog_mul;
		XDIV = xlog_div;
		XSUB = xlog_sub;
		XCOMPARE = xlog_comp;
		XCONVERT_LIN = xlog;
		XCONVERT_LOG = get_self;
		XCONVERT_2_LOG = get_self;
		XCONVERT_2_LIN = exp;
		ZERO = xlog(0);
	}
	else if (math_mode == MATH_MODE_LIN_SPACE)
	{
		fprintf(stderr, "Math mode linear.\n");
		XSUM = lin_sum;
		XMUL = lin_mul;
		XDIV = lin_div;
		XSUB = lin_sub;
		XCOMPARE = lin_compare;
		XCONVERT_LIN = get_self;
		XCONVERT_LOG = exp;
		XCONVERT_2_LOG = xlog;
		XCONVERT_2_LIN = get_self;
		ZERO = 0.0;
	}

	if (max_n_tag_vars_per_window > 0)
	{
		fprintf(stderr, "Using maximum of %d tag variants per window.\n", max_n_tag_vars_per_window);
	}

	double global_lin_scaler = exp(global_scaler_in_log);
	if (math_mode == MATH_MODE_LOG_SPACE)
	{
		global_lin_scaler = xlog(1.0);
	}

	fprintf(stderr, "Running Viterbi-ML sliding window with tag2tag distance of %.3f cMs (center-2-untyped target buffer of %.3f cMs) with effective pop. size N_e=%d (Scaler: %.3f)\n", tag_2_tag_distance_cM, l_target_2_center_pred_buffer_in_cM, (int)N_e, global_scaler_in_log);

	vector<t_annot_region*>* reference_haplocoded_tag_target_geno_regs = load_variant_signal_regions_wrapper(reference_tag_target_haplocoded_genotypes_fp, ref_sample_ids_list_fp);
	vector<char*>* ref_sample_ids = buffer_file(ref_sample_ids_list_fp);
	int n_ref_haplotypes = ref_sample_ids->size() * 2;
	fprintf(stderr, "Loaded %d variant regions for %d reference samples, %d haplotypes.\n",
		reference_haplocoded_tag_target_geno_regs->size(),
		ref_sample_ids->size(),
		n_ref_haplotypes);

	// Reset the scores to 0 to indicates all untyped.
	for (int i_reg = 0; i_reg < reference_haplocoded_tag_target_geno_regs->size(); i_reg++)
	{
		reference_haplocoded_tag_target_geno_regs->at(i_reg)->score = VAR_UNTYPED_TARGET;
	} // i_reg loop.

	vector<t_annot_region*>* testing_haplocoded_tag_geno_regs = load_variant_signal_regions_wrapper(testing_tag_haplocoded_genotypes_fp, testing_sample_ids_list_fp);
	vector<char*>* testing_sample_ids = buffer_file(testing_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d testing variant regions for %d testing samples.\n", testing_haplocoded_tag_geno_regs->size(), testing_sample_ids->size());

	fprintf(stderr, "Assigning testing tag regions with the reference variants.\n");
	vector<t_annot_region*>* intersects = intersect_annot_regions(testing_haplocoded_tag_geno_regs, reference_haplocoded_tag_target_geno_regs, true);
	fprintf(stderr, "Processing %d intersections.\n", intersects->size());
	for (int i_int = 0; i_int < intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* testing_reg = int_info->src_reg;
		t_annot_region* ref_reg = int_info->dest_reg;

		// Set this reference variant as "typed". 
		ref_reg->score = VAR_TYPED_TAG;

		void** testing_reg_info = (void**)(testing_reg->data);
		void** ref_reg_info = (void**)(ref_reg->data);
		testing_reg_info[1] = ref_reg_info[0];

		delete int_info;
	} // i_int loop.
	delete_annot_regions(intersects);

	vector<t_annot_region*>* focus_regs = NULL;
	if (check_file(target_focus_reg_BED_fp))
	{
		focus_regs = load_BED(target_focus_reg_BED_fp);
		fprintf(stderr, "Found %d focus regions.\n", focus_regs->size());
	}

	// Assign the haplocoded regions.
	vector<t_annot_region*>* known_haplocoded_tag_target_geno_regs = NULL;
	if (check_file(known_tag_target_haplocoded_genotypes_fp))
	{
		fprintf(stderr, "Loading known genotypes from %s\n", known_tag_target_haplocoded_genotypes_fp);

		known_haplocoded_tag_target_geno_regs = load_variant_signal_regions_wrapper(known_tag_target_haplocoded_genotypes_fp, testing_sample_ids_list_fp);
		fprintf(stderr, "Loaded %d known tag+target variants for %d samples.\n", known_haplocoded_tag_target_geno_regs->size(), testing_sample_ids->size());

		fprintf(stderr, "Assigning the known genotype regions to reference genotypes.\n");
		intersects = intersect_annot_regions(reference_haplocoded_tag_target_geno_regs, known_haplocoded_tag_target_geno_regs, true);
		fprintf(stderr, "Processing %d known genotype reference intersects.\n", intersects->size());
		for (int i_int = 0; i_int < intersects->size(); i_int++)
		{
			t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
			t_annot_region* ref_reg = int_info->src_reg;
			t_annot_region* known_reg = int_info->dest_reg;

			void** ref_reg_info = (void**)(ref_reg->data);
			void** known_reg_info = (void**)(known_reg->data);
			ref_reg_info[1] = known_reg_info[0];

			delete int_info;
		} // i_int loop.

		delete_annot_regions(intersects);
	} // known target check.

	// Set the imputed genotypes probabilities for untyped variants.
	for (int i_var = 0; i_var < reference_haplocoded_tag_target_geno_regs->size(); i_var++)
	{
		void** cur_ref_var_info = (void**)(reference_haplocoded_tag_target_geno_regs->at(i_var)->data);
		void** new_info = new void*[5];
		new_info[0] = cur_ref_var_info[0];
		new_info[1] = cur_ref_var_info[1];

		double*** imputed_allele_prob_per_sample_per_hap = new double**[testing_sample_ids->size() + 2];
		for (int i_s = 0; i_s < testing_sample_ids->size(); i_s++)
		{
			imputed_allele_prob_per_sample_per_hap[i_s] = new double*[2];

			imputed_allele_prob_per_sample_per_hap[i_s][0] = new double[2];
			imputed_allele_prob_per_sample_per_hap[i_s][0][0] = ZERO;
			imputed_allele_prob_per_sample_per_hap[i_s][0][1] = ZERO;

			imputed_allele_prob_per_sample_per_hap[i_s][1] = new double[2];
			imputed_allele_prob_per_sample_per_hap[i_s][1][0] = ZERO;
			imputed_allele_prob_per_sample_per_hap[i_s][1][1] = ZERO;
		} // i_s loop.
		new_info[2] = imputed_allele_prob_per_sample_per_hap;

		reference_haplocoded_tag_target_geno_regs->at(i_var)->data = new_info;
	} // i_var loop.

	int n_states = n_ref_haplotypes;
	int n_symbols = 2;

	t_restr_annot_region_list* restr_testing_haplocoded_tag_geno_regs = restructure_annot_regions(testing_haplocoded_tag_geno_regs);
	t_restr_annot_region_list* restr_reference_haplocoded_tag_target_geno_regs = restructure_annot_regions(reference_haplocoded_tag_target_geno_regs, restr_testing_haplocoded_tag_geno_regs->chr_ids);

	char arrays_op_fp[] = "viterbi_arrays.txt.gz";
	//FILE* f_viterbi_arrays = open_f(arrays_op_fp, "w");
	FILE* f_viterbi_arrays = NULL;
	for (int i_chr = 0; i_chr < restr_testing_haplocoded_tag_geno_regs->chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Running Viterbi on %s\n", restr_testing_haplocoded_tag_geno_regs->chr_ids->at(i_chr));

		char cur_chr_recombination_rate_fp[1000];
		sprintf(cur_chr_recombination_rate_fp, "%s/%s.map", recombination_rate_dir, restr_testing_haplocoded_tag_geno_regs->chr_ids->at(i_chr));
		vector<t_annot_region*>* cur_chrom_recomb_regs = load_recombination_rates(cur_chr_recombination_rate_fp);
		if (cur_chrom_recomb_regs == NULL)
		{
			fprintf(stderr, "Could not load recombination rates from %s\n", cur_chr_recombination_rate_fp);
			exit(0);
		}

		sort_set_sorting_info(cur_chrom_recomb_regs, sort_regions);

		// Add one element to the beginning of the regions for the current chromosome.
		vector<t_annot_region*>* cur_chr_testing_haplocoded_tag_geno_regs = restr_testing_haplocoded_tag_geno_regs->regions_per_chrom[i_chr];
		vector<t_annot_region*>* cur_chrom_ref_tag_target_var_regs = restr_reference_haplocoded_tag_target_geno_regs->regions_per_chrom[i_chr];

		// Assign the recomb rates.
		fprintf(stderr, "Setting recombination rates using %d recombination mapped regions.\n", cur_chrom_recomb_regs->size());
		for (int i_reg = 0; i_reg < cur_chr_testing_haplocoded_tag_geno_regs->size(); i_reg++)
		{
			double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant_optimized(cur_chr_testing_haplocoded_tag_geno_regs->at(i_reg), cur_chrom_recomb_regs);
			cur_chr_testing_haplocoded_tag_geno_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
		} // i_reg loop.

		for (int i_reg = 0; i_reg < cur_chrom_ref_tag_target_var_regs->size(); i_reg++)
		{
			double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant_optimized(cur_chrom_ref_tag_target_var_regs->at(i_reg), cur_chrom_recomb_regs);
			cur_chrom_ref_tag_target_var_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
		} // i_reg loop.

		  // Get the focus regions on this chromosome.
		int focus_start = cur_chrom_ref_tag_target_var_regs->at(0)->start;
		int focus_end = cur_chrom_ref_tag_target_var_regs->back()->start;
		if (focus_regs != NULL)
		{
			vector<t_annot_region*>* cur_chr_focus_regs = get_regions_per_chromosome(focus_regs, restr_testing_haplocoded_tag_geno_regs->chr_ids->at(i_chr));
			sort(cur_chr_focus_regs->begin(), cur_chr_focus_regs->end(), sort_regions);
			focus_start = cur_chr_focus_regs->at(0)->start;
			focus_end = cur_chr_focus_regs->back()->end;
			delete cur_chr_focus_regs;
		}

		fprintf(stderr, "Focusing on variants @ %d-%d\n", focus_start, focus_end);

		// Find all must be set to false.
		vector<t_annot_region*>* intersects = intersect_annot_regions(cur_chrom_ref_tag_target_var_regs, focus_regs, false);
		vector<t_annot_region*>* cur_chrom_focus_ref_target_var_regs = new vector<t_annot_region*>();
		for (int i_int = 0; i_int < intersects->size(); i_int++)
		{
			t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
			t_annot_region* cur_target_reg = int_info->src_reg;
			cur_chrom_focus_ref_target_var_regs->push_back(cur_target_reg);
		} // i_int loop.

		  // Sort the targets in the focus region.
		sort(cur_chrom_focus_ref_target_var_regs->begin(), cur_chrom_focus_ref_target_var_regs->end(), sort_regions);

		vector<t_annot_region*>* prev_win_var_regs = NULL;
		for (int test_sample_i = 0; test_sample_i < testing_sample_ids->size(); test_sample_i++)
		{
			prev_win_var_regs = NULL;

			// Reset all the variant types to untyped.
			for (int i_ref_var = 0; i_ref_var < cur_chrom_focus_ref_target_var_regs->size(); i_ref_var++)
			{
				if (cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->score == VAR_UNTYPED_TARGET_IMPUTED)
				{
					cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->score = VAR_UNTYPED_TARGET;
				}
			} // i_ref_var loop.

			// We keep track of the middle window start coordinates.
			for (int win_mid_start_var_i = 0; win_mid_start_var_i < cur_chr_testing_haplocoded_tag_geno_regs->size(); win_mid_start_var_i++)
			{
				// If the mid-window passed the focus region's end completely, we are done.
				if (cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_start_var_i)->start > focus_end)
				{
					break;
				}

				// We are imputing the variants between current and next variants; select window start and end coordinates.
				double mid_start_cM = cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_start_var_i)->dbl_score;
				int win_mid_end_var_i = win_mid_start_var_i;
				while (win_mid_end_var_i < cur_chr_testing_haplocoded_tag_geno_regs->size() &&
					fabs(mid_start_cM - cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_end_var_i)->dbl_score) < l_target_2_center_pred_buffer_in_cM)
				{
					win_mid_end_var_i++;
				} // win_end_var_i loop.

				if (win_mid_end_var_i >= cur_chr_testing_haplocoded_tag_geno_regs->size())
				{
					win_mid_end_var_i = cur_chr_testing_haplocoded_tag_geno_regs->size() - 1;
				}

				// Get the middle cM, this is the middle of both windows.
				double mid_cM = (cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_start_var_i)->dbl_score +
					cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_end_var_i)->dbl_score) / 2;

				double start_cM = MAX(0, mid_cM - (tag_2_tag_distance_cM / 2));
				int win_start_var_i = win_mid_start_var_i;
				while (win_start_var_i < cur_chr_testing_haplocoded_tag_geno_regs->size() &&
					win_start_var_i > 0 &&
					cur_chr_testing_haplocoded_tag_geno_regs->at(win_start_var_i)->dbl_score > start_cM)
				{
					win_start_var_i--;
				} // win_end_var_i loop.

				  // GEt the ending coordiantes for the full window.
				double end_cM = mid_cM + (tag_2_tag_distance_cM / 2);
				int win_end_var_i = win_mid_end_var_i;
				while (win_end_var_i < cur_chr_testing_haplocoded_tag_geno_regs->size() &&
					win_end_var_i > 0 &&
					cur_chr_testing_haplocoded_tag_geno_regs->at(win_end_var_i)->dbl_score < end_cM)
				{
					win_end_var_i++;
				} // win_end_var_i loop.

				  // Fix end if passed over the end.
				if (win_end_var_i >= cur_chr_testing_haplocoded_tag_geno_regs->size())
				{
					win_end_var_i = cur_chr_testing_haplocoded_tag_geno_regs->size() - 1;
				}

				fprintf(stderr, "Sample %d: Prep. win. var_i: %d-(%d-%d)-%d ;; posn: %d-(%d-%d)-%d ;; cM: %.3f-(%.3f-%.3f)-%.3f\n",
						test_sample_i,
						win_start_var_i,
						win_mid_start_var_i,
						win_mid_end_var_i,
						win_end_var_i,
						cur_chr_testing_haplocoded_tag_geno_regs->at(win_start_var_i)->start,
						cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_start_var_i)->start,
						cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_end_var_i)->start,
						cur_chr_testing_haplocoded_tag_geno_regs->at(win_end_var_i)->start,
						cur_chr_testing_haplocoded_tag_geno_regs->at(win_start_var_i)->dbl_score,
						cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_start_var_i)->dbl_score,
						cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_end_var_i)->dbl_score,
						cur_chr_testing_haplocoded_tag_geno_regs->at(win_end_var_i)->dbl_score);

				/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				// Count the number of targets in this window.				
				int n_target_vars_in_cur_win = 0;
				for (int i_ref_var = 0; i_ref_var < cur_chrom_focus_ref_target_var_regs->size(); i_ref_var++)
				{
					double target_cM = cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->dbl_score;
					double left_win_cM = cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_start_var_i)->dbl_score;
					double right_win_cM = cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_end_var_i)->dbl_score;

					// Make sure this is an untyped variant and it is covered well by tags.
					if (cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->score == VAR_UNTYPED_TARGET_IMPUTED ||
						cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->score == VAR_TYPED_TAG ||
						target_cM < left_win_cM ||
						target_cM > right_win_cM)
					{

					}
					else
					{
						n_target_vars_in_cur_win++;
					}

					if (cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->start > cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_end_var_i)->end)
					{
						break;
					}
				} // i_ref_var loop.

				if (n_target_vars_in_cur_win == 0)
				{
					if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
					{
						fprintf(stderr, "Skipping window %d-%d\n",
							cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_start_var_i)->start,
							cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_end_var_i)->start);
					}

					continue;
				}
				else
				{
					if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
					{
						fprintf(stderr, "%d target variants in window\n", n_target_vars_in_cur_win);
					}
				}
				/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

				start_cM = cur_chr_testing_haplocoded_tag_geno_regs->at(win_start_var_i)->dbl_score;
				end_cM = cur_chr_testing_haplocoded_tag_geno_regs->at(win_end_var_i)->dbl_score;

				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				// Set the current window's tag variant regions.
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				vector<t_annot_region*>* cur_win_var_regs = new vector<t_annot_region*>();

				// Subselect windows if requested.
				if (max_n_tag_vars_per_window == 0 ||
					(win_end_var_i - win_start_var_i + 1) < max_n_tag_vars_per_window)
				{
					cur_win_var_regs->insert(cur_win_var_regs->end(),
						cur_chr_testing_haplocoded_tag_geno_regs->begin() + win_start_var_i,
						cur_chr_testing_haplocoded_tag_geno_regs->begin() + win_end_var_i + 1);
				}
				else
				{
					cur_win_var_regs->push_back(cur_chr_testing_haplocoded_tag_geno_regs->at(win_start_var_i));
					cur_win_var_regs->push_back(cur_chr_testing_haplocoded_tag_geno_regs->at(win_end_var_i));

					double n_vars_2_fit = max_n_tag_vars_per_window - 2;
					double n_vars_in_window = win_end_var_i - win_start_var_i + 1 - 2;
					double frac_i_update = n_vars_in_window / n_vars_2_fit;
					double frac_i = win_start_var_i + 1;

					if (frac_i_update < 1)
					{
						fprintf(stderr, "Sanity check failed: frac_i_update < 1: %lf\n", frac_i_update);
						exit(0);
					}

					int last_var_i = -1;
					while (1)
					{
						int var_i = (int)(floor(frac_i));

						if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
						{
							fprintf(stderr, "Selecting var %d (%.4f) in window @ (%d-%d) [%d windows already]\n", var_i, frac_i, win_start_var_i, win_end_var_i,
								cur_win_var_regs->size());
						}

						if (var_i >= win_end_var_i - 1)
						{
							break;
						}

						if (last_var_i == var_i)
						{
							fprintf(stderr, "Sanity check failed while subsampling windows: %d is repeated.\n");
							exit(0);
						}

						cur_win_var_regs->push_back(cur_chr_testing_haplocoded_tag_geno_regs->at(var_i));
						last_var_i = var_i;

						frac_i += frac_i_update;
					} // frac_i loop.

					if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
					{
						fprintf(stderr, "Subsampled %d/%d::%d tag variants in the window.\n",
							win_end_var_i - win_start_var_i + 1,
							cur_win_var_regs->size(), max_n_tag_vars_per_window);
					}

					sort(cur_win_var_regs->begin(), cur_win_var_regs->end(), sort_regions);
				} // var subsampling check.			

				// Make sure we inserted the correct elements.
				if (!t_string::compare_strings(cur_win_var_regs->at(0)->name, cur_chr_testing_haplocoded_tag_geno_regs->at(win_start_var_i)->name) ||
					!t_string::compare_strings(cur_win_var_regs->back()->name, cur_chr_testing_haplocoded_tag_geno_regs->at(win_end_var_i)->name))
				{
					fprintf(stderr, "Start and end did not copy as expected.\n");
					exit(0);
				}

				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
				{
					fprintf(stderr, "Computing ML scores for sample %d\n", test_sample_i);
				}

				// Generate the blocks and the reduced states: These must be generated before ends are added.
				vector<t_var_block*>* haplo_var_blocks = generate_reduced_state_blocks_constant_size_blocks(cur_win_var_regs,
					ref_sample_ids,
					1,
					l_blocks);

				fprintf(stderr, "Window (%d/%d/%d tags)::%d targets: %s:%d-%d (%d-%d) spanning %.4f cM between %.4f-%.4f cMs [%d haplo-blocks of %d variants]\n",
					cur_win_var_regs->size(), max_n_tag_vars_per_window, win_end_var_i - win_start_var_i + 1,
					n_target_vars_in_cur_win,
					cur_win_var_regs->at(0)->chrom, cur_win_var_regs->at(0)->start, cur_win_var_regs->back()->start,
					win_start_var_i, win_end_var_i,
					fabs(start_cM - end_cM),
					start_cM, end_cM,
					haplo_var_blocks->size(), l_blocks);

				if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
				{
					fprintf(stderr, "%d haplotype blocks:\n", haplo_var_blocks->size());
					for (int i_block = 0; i_block < haplo_var_blocks->size(); i_block++)
					{
						fprintf(stderr, "%d-%d: %d haplotypes\n",
							haplo_var_blocks->at(i_block)->start_var_i, haplo_var_blocks->at(i_block)->end_var_i,
							haplo_var_blocks->at(i_block)->haplogroup_haplo_indices->size());
					} // i_block loop.
				}

				// Set the start and end state regions.
				t_annot_region* start_state_reg = NULL;
				cur_win_var_regs->insert(cur_win_var_regs->begin(), start_state_reg);

				// Set the number of variants for indexing purposes here.
				int n_vars = cur_win_var_regs->size() - 1;

				// Add the end state region.
				t_annot_region* end_state_reg = NULL;
				cur_win_var_regs->push_back(end_state_reg);

				// Allocate and initialize the forward/backward arrays.
				double*** ML_scores_per_hap = new double**[n_vars + 2];
				int*** ML_prev_state_per_hap = new int**[n_vars + 2];
				for (int var_i = 0; var_i <= n_vars + 1; var_i++)
				{
					ML_scores_per_hap[var_i] = new double*[2];
					ML_prev_state_per_hap[var_i] = new int*[2];

					for (int test_hap = 0; test_hap < 2; test_hap++)
					{
						ML_scores_per_hap[var_i][test_hap] = new double[n_states + 2];
						memset(ML_scores_per_hap[var_i][test_hap], 0, sizeof(double) * (n_states + 1));

						ML_prev_state_per_hap[var_i][test_hap] = new int[n_states + 2];
						memset(ML_prev_state_per_hap[var_i][test_hap], 0, sizeof(int) * (n_states + 1));
					}
				} // i loop.

				// Initialize the state probabilities for both haplotypes.
				for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
				{
					ML_scores_per_hap[0][0][state_i] = XCONVERT_LIN((double)1.0 / n_ref_haplotypes);
					ML_scores_per_hap[0][1][state_i] = XCONVERT_LIN((double)1.0 / n_ref_haplotypes);
				} // state_i loop.

				// Process each haplotype of the test sample separately.
				for (int test_hap_i = 0; test_hap_i < 2; test_hap_i++)
				{
					// Process each block sequentially.
					for (int block_i = 0; block_i < haplo_var_blocks->size(); block_i++)
					{
						// Process the current block.
						int cur_block_start_var_i = haplo_var_blocks->at(block_i)->start_var_i;
						int cur_block_end_var_i = haplo_var_blocks->at(block_i)->end_var_i;

						int n_cur_block_haplogrps = haplo_var_blocks->at(block_i)->haplogroup_haplo_indices->size();
						vector<vector<int>*>* per_hplgrp_hap_i = haplo_var_blocks->at(block_i)->haplogroup_haplo_indices;

						if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
						{
							fprintf(stderr, "Processing test_hap_i:%d;; Block:%d-%d with %d haplogroups.\n",
								test_hap_i,
								cur_block_start_var_i, cur_block_end_var_i,
								n_cur_block_haplogrps);
						}

						////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
						// 1. Select the representative (top-scoring) sequence per haplogroup.
						vector<int>* per_hplgrp_top_scoring_hap_i = new vector<int>();
						for (int i_grp = 0; i_grp < n_cur_block_haplogrps; i_grp++)
						{
							int cur_grp_top_scoring_hap_i = -1;
							double cur_grp_top_scoring_hap_score = -1;
							for (int i_hap_i = 0;
								i_hap_i < per_hplgrp_hap_i->at(i_grp)->size();
								i_hap_i++)
							{
								int ref_hap_i = per_hplgrp_hap_i->at(i_grp)->at(i_hap_i);
								if (cur_grp_top_scoring_hap_score < ML_scores_per_hap[cur_block_start_var_i - 1][test_hap_i][ref_hap_i])
								{
									cur_grp_top_scoring_hap_i = ref_hap_i;
									cur_grp_top_scoring_hap_score = ML_scores_per_hap[cur_block_start_var_i - 1][test_hap_i][ref_hap_i];
								}
							} // i_sample_i loop.

							per_hplgrp_top_scoring_hap_i->push_back(cur_grp_top_scoring_hap_i);

							if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
							{
								fprintf(stderr, "Block %d::Haplogroup %d; max_scoring_hap: %d (%.4f)\n",
									block_i,
									i_grp,
									cur_grp_top_scoring_hap_i,
									cur_grp_top_scoring_hap_score);
							}
						} // i_grp loop.

						////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
						// 2. Run viterbi on only the top scoring haplotypes (states) and update their scores.
						for (int var_i = cur_block_start_var_i;
							var_i <= cur_block_end_var_i;
							var_i++)
						{
							if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
							{
								if (var_i % 10 == 0)
								{
									fprintf(stderr, "Viterbi: sample_i: %d/%d: var_i: %d         \r", test_sample_i, testing_sample_ids->size(), var_i);
								}
							}

							void** cur_tag_var_info = NULL;
							char* test_sample_geno = NULL;
							char* ref_sample_geno = NULL;

							if (var_i <= n_vars)
							{
								void** cur_tag_var_info = (void**)(cur_win_var_regs->at(var_i)->data);
								test_sample_geno = (char*)(cur_tag_var_info[0]);
								ref_sample_geno = (char*)(cur_tag_var_info[1]);
							}

							/////////////////////////////////////////////////////
							// Pre-compute the transition probabilities:
							double prev_var_cM = cur_win_var_regs->at(1)->dbl_score;
							if (var_i > 1)
							{
								prev_var_cM = cur_win_var_regs->at(var_i - 1)->dbl_score;
							}

							double cur_var_cM = cur_win_var_regs->at(n_vars)->dbl_score;
							if (var_i <= n_vars)
							{
								cur_var_cM = cur_win_var_regs->at(var_i)->dbl_score;
							}

							double r_m = fabs(cur_var_cM - prev_var_cM);
							double rho_m = 4 * N_e * r_m;

							double tau_m = 1 - exp(-1 * rho_m / n_ref_haplotypes);

							double other_prob = tau_m / n_ref_haplotypes;
							double self_prob = (1 - tau_m) + (tau_m / n_ref_haplotypes);

							//fprintf(stderr, "Self prob: %.5f, other_prob: %.5f (%.4f, %.4f)\n", self_prob, other_prob, prev_var_cM, cur_var_cM);
							/////////////////////////////////////////////////////

							double max_ML_score_per_cur_var = ZERO;
							double min_ML_score_per_cur_var = pow(10, 20);
							for (int hapgrp_i = 0; hapgrp_i < n_cur_block_haplogrps; hapgrp_i++)
							{
								for (int i_state_i = 0; i_state_i < per_hplgrp_hap_i->at(hapgrp_i)->size(); i_state_i++)
								{
									int state_i = per_hplgrp_hap_i->at(hapgrp_i)->at(i_state_i);

									// Recurse over the previous states.
									ML_scores_per_hap[var_i][test_hap_i][state_i] = ZERO;

									//////////////////////////////////////////////////////////////////////////////////////////////////////////
									//////////////////////////////////////////////////////////////////////////////////////////////////////////
									// First, do a self transition of the current state.
									double trans_prob = self_prob;
									double emit_prob = ZERO;

									// There is no conditional below.
									int ref_sample_i = state_i / 2;
									int ref_hap_i = state_i % 2;
									int cur_ref_allele = get_allele_per_haplotype(ref_sample_geno[ref_sample_i], ref_hap_i);
									int cur_test_allele = get_allele_per_haplotype(test_sample_geno[test_sample_i], test_hap_i);

									double emit_prob_per_ref_allele[2];
									emit_prob_per_ref_allele[cur_ref_allele] = 1 - allele_err;
									emit_prob_per_ref_allele[1 - cur_ref_allele] = allele_err;
									emit_prob = emit_prob_per_ref_allele[cur_test_allele];

									double trans_emit_prob = XMUL(XCONVERT_LIN(trans_prob), XCONVERT_LIN(emit_prob));

									// Add the scaler for this position.
									trans_emit_prob = XMUL(trans_emit_prob, XCONVERT_LIN(global_lin_scaler));

									ML_scores_per_hap[var_i][test_hap_i][state_i] = MAX(ML_scores_per_hap[var_i][test_hap_i][state_i],
										XMUL(trans_emit_prob, ML_scores_per_hap[var_i - 1][test_hap_i][state_i]));

									if (XCOMPARE(ML_scores_per_hap[var_i][test_hap_i][state_i], XMUL(trans_emit_prob, ML_scores_per_hap[var_i - 1][test_hap_i][state_i])))
									{
										ML_prev_state_per_hap[var_i][test_hap_i][state_i] = state_i;
									}

									//////////////////////////////////////////////////////////////////////////////////////////////////////////
									//////////////////////////////////////////////////////////////////////////////////////////////////////////
									// Now do a transition to all the top scoring haplotypes.
									for (int i_prev_state_i = 0; i_prev_state_i < n_cur_block_haplogrps; i_prev_state_i++)
									{
										int prev_state_i = per_hplgrp_top_scoring_hap_i->at(i_prev_state_i);

										// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
										// Set the transition probabilities.
										double trans_prob = 0;
										if (state_i == prev_state_i)
										{
											trans_prob = self_prob;
										}
										else
										{
											trans_prob = other_prob;
										}

										// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
										// Set the emission prob: If the ref allele is matching to test, set to 0.99999, otherwise set to the error prob.
										double emit_prob = ZERO;

										// There is no conditional below.
										int ref_sample_i = state_i / 2;
										int ref_hap_i = state_i % 2;
										int cur_ref_allele = get_allele_per_haplotype(ref_sample_geno[ref_sample_i], ref_hap_i);
										int cur_test_allele = get_allele_per_haplotype(test_sample_geno[test_sample_i], test_hap_i);

										double emit_prob_per_ref_allele[2];
										emit_prob_per_ref_allele[cur_ref_allele] = 1 - allele_err;
										emit_prob_per_ref_allele[1 - cur_ref_allele] = allele_err;
										emit_prob = emit_prob_per_ref_allele[cur_test_allele];
										// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 

										double trans_emit_prob = XMUL(XCONVERT_LIN(trans_prob), XCONVERT_LIN(emit_prob));

										// Add the scaler for this position.
										trans_emit_prob = XMUL(trans_emit_prob, XCONVERT_LIN(global_lin_scaler));

										ML_scores_per_hap[var_i][test_hap_i][state_i] = MAX(ML_scores_per_hap[var_i][test_hap_i][state_i],
											XMUL(trans_emit_prob, ML_scores_per_hap[var_i - 1][test_hap_i][prev_state_i]));

										if (XCOMPARE(ML_scores_per_hap[var_i][test_hap_i][state_i], XMUL(trans_emit_prob, ML_scores_per_hap[var_i - 1][test_hap_i][prev_state_i])))
										{
											ML_prev_state_per_hap[var_i][test_hap_i][state_i] = prev_state_i;
										}

										if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
										{
											fprintf(stderr, "ML[%d][%d]: %.5f: P_emit=%.5f; P_trans(%d->%d): %.5f\n", var_i, state_i, ML_scores_per_hap[var_i][test_hap_i][state_i],
												emit_prob, prev_state_i, state_i, trans_prob);
										}
									} // prev_state_i loop.

									// Save.
									if (f_viterbi_arrays != NULL)
									{
										fprintf(f_viterbi_arrays, "%d\t%d\t%d\t%d\t%d\t%.5f\n", win_start_var_i, test_sample_i, test_hap_i, var_i, state_i, XCONVERT_2_LOG(ML_scores_per_hap[var_i][test_hap_i][state_i]));
									}

									// Update the maximum and the minimum values. We need the maximum for scaling scores below.
									if (min_ML_score_per_cur_var > ML_scores_per_hap[var_i][test_hap_i][state_i])
									{
										min_ML_score_per_cur_var = ML_scores_per_hap[var_i][test_hap_i][state_i];
									}

									if (max_ML_score_per_cur_var < ML_scores_per_hap[var_i][test_hap_i][state_i])
									{
										max_ML_score_per_cur_var = ML_scores_per_hap[var_i][test_hap_i][state_i];
									}
								} // i_state_i loop.
							} // hapgrp_i loop.

							// Keep track of the highest among haplogroups.
							if (math_mode == MATH_MODE_LIN_SPACE)
							{
								double largest_max_log10_value = 50;

								double log10_max_ML_score = log(max_ML_score_per_cur_var) / log(10);
								if (log10_max_ML_score > largest_max_log10_value ||
									log10_max_ML_score < -largest_max_log10_value)
								{
									double rescaling_factor = pow(10, largest_max_log10_value - log10_max_ML_score);

									fprintf(stderr, "Top score out of range @ %d: %.5f; Re-scaling maximum to 10^(%.5f) with rescaling value %.5f\n",
										var_i, log(max_ML_score_per_cur_var) / log(10), largest_max_log10_value, log(rescaling_factor) / log(10));

									for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
									{
										ML_scores_per_hap[var_i][test_hap_i][state_i] = XMUL(ML_scores_per_hap[var_i][test_hap_i][state_i], rescaling_factor);
									} // state_i loop.
								} // overflow/underflow check for the top path's score.
							}
						} // var_i loop.
						////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

						delete per_hplgrp_top_scoring_hap_i;
					} // block_i loop.
				} // test_hap_i loop.

				  // Free haplogroup and block memories.
				for (int i_block = 0; i_block < haplo_var_blocks->size(); i_block++)
				{
					for (int i_grp = 0; i_grp < haplo_var_blocks->at(i_block)->haplogroup_haplo_indices->size(); i_grp++)
					{
						delete haplo_var_blocks->at(i_block)->haplogroup_haplo_indices->at(i_grp);
					}
				} // i_block loop.

				// Compute the total log forward and backward probabilities.
				double per_hap_total_log_ML_prob[2];
				per_hap_total_log_ML_prob[0] = ZERO;
				per_hap_total_log_ML_prob[1] = ZERO;

				int per_hap_ML_state_at_end[2];
				per_hap_ML_state_at_end[0] = -1;
				per_hap_ML_state_at_end[1] = -1;

				for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
				{
					per_hap_total_log_ML_prob[0] = MAX(per_hap_total_log_ML_prob[0], ML_scores_per_hap[n_vars][0][state_i]);
					if (per_hap_total_log_ML_prob[0] == ML_scores_per_hap[n_vars][0][state_i])
					{
						per_hap_ML_state_at_end[0] = state_i;
					}

					per_hap_total_log_ML_prob[1] = MAX(per_hap_total_log_ML_prob[1], ML_scores_per_hap[n_vars][1][state_i]);
					if (per_hap_total_log_ML_prob[1] == ML_scores_per_hap[n_vars][1][state_i])
					{
						per_hap_ML_state_at_end[1] = state_i;
					}
				} // state_i loop.

				if (per_hap_total_log_ML_prob[0] == 0 ||
					per_hap_total_log_ML_prob[1] == 0)
				{
					fprintf(stderr, "Final score is 0!\n");
					exit(0);
				}

				fprintf(stderr, "Sample %d: Haplotype 0 ML probabilities: ML=%.5f\n", test_sample_i, log(per_hap_total_log_ML_prob[0]));
				fprintf(stderr, "Sample %d: Haplotype 1 ML probabilities: ML=%.5f\n", test_sample_i, log(per_hap_total_log_ML_prob[1]));

				// Trace the states back to identify the optimal path.
				fprintf(stderr, "Tracing back the Viterbi path for each haplotype.\n");
				int** per_hap_per_variant_viterbi_haplotype = new int*[2];
				for (int hap_i = 0; hap_i < 2; hap_i++)
				{
					int cur_state = per_hap_ML_state_at_end[hap_i];
					int* per_variant_viterbi_haplotype = new int[n_vars + 2];
					for (int var_i = n_vars; var_i >= 1; var_i--)
					{
						// Set the optimal state for the current index.
						per_variant_viterbi_haplotype[var_i] = cur_state;

						// Get the new state for the previous variant.
						int new_state = ML_prev_state_per_hap[var_i][hap_i][cur_state];

						// Reset the current state for the previous variant.
						cur_state = new_state; // This is the state for the previous value..
					} // var_i loop.

					per_hap_per_variant_viterbi_haplotype[hap_i] = per_variant_viterbi_haplotype;
				} // hap_i loop.

				// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
				// Start computing the posterior probabilities for the untyped markers.
				// Now go over the untyped genotypes and compute the posterior probability of all the untyped markers.				
				int n_untyped = 0;

				if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
				{
					fprintf(stderr, "Computing the probabilities of the %d untyped markers in %d-%d.\n",
						cur_chrom_focus_ref_target_var_regs->size(),
						n_target_vars_in_cur_win,
						win_mid_start_var_i, win_mid_end_var_i);
				}

				//// +1 is necessary to match the 1-based coordinates.
				//int rel_win_mid_start_i = win_mid_start_var_i - win_start_var_i + 1;
				//int rel_win_mid_end_i = win_mid_end_var_i - win_start_var_i + 1;

				//// Relative coordinates must match exactly to the absolute coordinates.
				//if (cur_win_var_regs->at(rel_win_mid_end_i)->dbl_score != cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_end_var_i)->dbl_score ||
				//	cur_win_var_regs->at(rel_win_mid_start_i)->dbl_score != cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_start_var_i)->dbl_score)
				//{
				//	fprintf(stderr, "Sanuty check failed: Start: %lf, %lf; End: %lf, %lf\n",
				//		cur_win_var_regs->at(rel_win_mid_start_i)->dbl_score, cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_start_var_i)->dbl_score,
				//		cur_win_var_regs->at(rel_win_mid_end_i)->dbl_score, cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_end_var_i)->dbl_score);
				//	exit(0);
				//}

				if (win_mid_start_var_i > win_end_var_i ||
					win_mid_end_var_i > win_end_var_i)
				{
					fprintf(stderr, "Sanity check failed: %d-(%d-%d)-%d\n", win_start_var_i, win_mid_start_var_i, win_mid_end_var_i, win_end_var_i);
					exit(0);
				}

				// Go over all the reference target variants in the focus region.
				for (int i_ref_var = 0; i_ref_var < cur_chrom_focus_ref_target_var_regs->size(); i_ref_var++)
				{
					double target_cM = cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->dbl_score;
					/*double left_mid_cM = cur_win_var_regs->at(rel_win_mid_start_i)->dbl_score;
					double right_mid_cM = cur_win_var_regs->at(rel_win_mid_end_i)->dbl_score;*/

					double left_mid_cM = cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_start_var_i)->dbl_score;
					double right_mid_cM = cur_chr_testing_haplocoded_tag_geno_regs->at(win_mid_end_var_i)->dbl_score;

					// Make sure this is an untyped variant and it is covered well by tags.
					if (cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->score == VAR_UNTYPED_TARGET_IMPUTED ||
						cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->score == VAR_TYPED_TAG ||
						target_cM < left_mid_cM ||
						target_cM > right_mid_cM)
					{
						if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
						{
							fprintf(stderr, "Skipping %s:%d @ cM: %.4f (%.4f-%.4f) typed: %d (%.4f, %.4f)\n",
								cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->chrom, cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->start,
								target_cM, left_mid_cM, right_mid_cM,
								cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->score,
								fabs(right_mid_cM - target_cM),
								fabs(left_mid_cM - target_cM));
						}
						continue;
					}

					// This is the untyped variant's information.
					void** cur_untyped_var_ref_var_info = (void**)(cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->data);
					char* untyped_var_ref_geno = (char*)(cur_untyped_var_ref_var_info[0]);
					char* known_var_ref_geno = (char*)(cur_untyped_var_ref_var_info[1]);
					double*** imputed_allele_prob_per_sample_per_hap = (double***)(cur_untyped_var_ref_var_info[2]);

					// If there is probability assigned to this variant for this sample, do not process it.
					if (imputed_allele_prob_per_sample_per_hap[test_sample_i][0][0] != ZERO ||
						imputed_allele_prob_per_sample_per_hap[test_sample_i][0][1] != ZERO)
					{
						continue;
					}

					n_untyped++;

					//fprintf(stderr, "@%d. untyped marker         \r", n_untyped);

					// Find the closest tag variants to this target:
					//for (int var_i = rel_win_mid_start_i; var_i < rel_win_mid_end_i; var_i++)
					for (int var_i = 1; var_i < n_vars; var_i++)
					{
						// Make sure the tag variants span the middle region.
						if (cur_win_var_regs->at(var_i + 1)->dbl_score < left_mid_cM ||
							cur_win_var_regs->at(var_i)->dbl_score > right_mid_cM)
						{
							continue;
						}

						// Find the untyped target variant in between the current tags.
						if (cur_win_var_regs->at(var_i)->start < cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->start &&
							cur_win_var_regs->at(var_i + 1)->start > cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->start)
						{
							if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
							{
								fprintf(stderr, "Found %d-%d-%d\n",
									cur_win_var_regs->at(var_i)->start,
									cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->start,
									cur_win_var_regs->at(var_i + 1)->start);
							}

							void** cur_tag_var_info = (void**)(cur_win_var_regs->at(var_i)->data);
							char* test_sample_geno = (char*)(cur_tag_var_info[0]);
							char* ref_sample_geno = (char*)(cur_tag_var_info[1]);

							// Get the self transition probability.
							double cur_cM = cur_win_var_regs->at(var_i)->dbl_score;
							double next_cM = cur_cM;
							if (var_i <= n_vars - 1)
							{
								next_cM = cur_win_var_regs->at(var_i + 1)->dbl_score;
							}

							double target_cM = cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->dbl_score;

							double left_weight = .5;
							double right_weight = .5;

							if (cur_cM != next_cM)
							{
								left_weight = fabs(target_cM - cur_cM) / (fabs(cur_cM - next_cM));
								right_weight = fabs(target_cM - next_cM) / (fabs(cur_cM - next_cM));
							}

							/*
							double r_m = fabs(cur_cM - next_cM);
							double rho_m = 4 * N_e * r_m;

							double tau_m = 1 - exp(-1 * rho_m / n_ref_haplotypes);

							double other_prob = tau_m / n_ref_haplotypes;
							double self_prob = (1 - tau_m) + (tau_m / n_ref_haplotypes);
							*/
							//////////////////////////////////////////////////////////////////////////////////////////////////

							for (int hap_i = 0; hap_i < 2; hap_i++)
							{
								int cur_hap_allele = -1;

								int var_ML_hap_i = per_hap_per_variant_viterbi_haplotype[hap_i][var_i];

								// For this untyped state in the refernce, we would like to sum all the paths that pass through here.
								int ref_sample_i = var_ML_hap_i / 2;
								int ref_hap_i = var_ML_hap_i % 2;

								int cur_ML_allele = get_allele_per_haplotype(untyped_var_ref_geno[ref_sample_i], ref_hap_i);

								// Get the haplotype next variant.
								int next_var_ML_hap_i = per_hap_per_variant_viterbi_haplotype[hap_i][var_i + 1];
								ref_sample_i = next_var_ML_hap_i / 2;
								ref_hap_i = next_var_ML_hap_i % 2;

								int next_ML_allele = get_allele_per_haplotype(untyped_var_ref_geno[ref_sample_i], ref_hap_i);

								// Set the allelic probabilities based on weights, i.e., genetic distance ratio to left and right closest variants.
								double per_allele_prob[2];
								per_allele_prob[0] = 0;
								per_allele_prob[1] = 0;
								//per_allele_prob[cur_ML_allele]++;
								per_allele_prob[cur_ML_allele] += left_weight;
								per_allele_prob[next_ML_allele] += right_weight;

								int known_geno = -1;
								if (known_var_ref_geno)
								{
									known_geno = (int)(known_var_ref_geno[test_sample_i]);
								}

								if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
								{
									fprintf(stderr, "Var %d on sample %d: %s:%d (%s): Hap%d: [Viterbi Allele Probs: %.3f/%.3f (Tot. weight: %.3f)]; (known_geno=%d)\n",
										i_ref_var,
										test_sample_i,
										cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->chrom,
										cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->start,
										cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->name,
										hap_i,
										per_allele_prob[0], per_allele_prob[1], left_weight + right_weight,
										known_geno);
								}

								imputed_allele_prob_per_sample_per_hap[test_sample_i][hap_i][0] = per_allele_prob[0];
								imputed_allele_prob_per_sample_per_hap[test_sample_i][hap_i][1] = per_allele_prob[1];
							} // hap_i loop.

							// Set the inputated variant type.
							if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
							{
								fprintf(stderr, "Setting target ref_var_i=%d (%d @ %.4f)\n",
									i_ref_var,
									cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->start,
									cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->dbl_score);
							}

							// Set the inputated variant type.
							cur_chrom_focus_ref_target_var_regs->at(i_ref_var)->score = VAR_UNTYPED_TARGET_IMPUTED;

							break; // Breaks out of position search.
						} // middle check.
					} // i_var loop.
				} // i_ref_var loop.

				for (int hap_i = 0; hap_i < 2; hap_i++)
				{
					delete[] per_hap_per_variant_viterbi_haplotype[hap_i];
				}
				delete[] per_hap_per_variant_viterbi_haplotype;

				// Allocate and initialize the ML arrays.
				for (int var_i = 0; var_i <= n_vars + 1; var_i++)
				{
					for (int test_hap = 0; test_hap < 2; test_hap++)
					{
						delete[] ML_scores_per_hap[var_i][test_hap];
						delete[] ML_prev_state_per_hap[var_i][test_hap];
					} // test_hap loop.

					delete[] ML_scores_per_hap[var_i];
					delete[] ML_prev_state_per_hap[var_i];
				} // i loop.

				delete[] ML_scores_per_hap;
				delete[] ML_prev_state_per_hap;

				// Copy the current windows; this is necessary to make sure we move forward.
				if (prev_win_var_regs != NULL)
				{
					delete prev_win_var_regs;
				}

				prev_win_var_regs = cur_win_var_regs;
			} // win_mid_start_var_i loop.
		} // test_sample_i loop.
	} // i_chr loop.

	if (f_viterbi_arrays != NULL)
	{
		close_f(f_viterbi_arrays, arrays_op_fp);
	}

	// Save the imputed genotype probabilities.
	FILE* f_geno_probs = open_f(geno_probs_op_fp, "w");
	for (int i_ref_var = 0; i_ref_var < reference_haplocoded_tag_target_geno_regs->size(); i_ref_var++)
	{
		// Make sure this is an untyped variant and it is covered well by tags.
		if (reference_haplocoded_tag_target_geno_regs->at(i_ref_var)->score != VAR_UNTYPED_TARGET_IMPUTED)
		{
			continue;
		}

		// This is the untyped variant's information.
		void** cur_untyped_var_ref_var_info = (void**)(reference_haplocoded_tag_target_geno_regs->at(i_ref_var)->data);
		double*** imputed_allele_prob_per_sample_per_hap = (double***)(cur_untyped_var_ref_var_info[2]);

		fprintf(f_geno_probs, "%s\t%d\t%d\t%s", reference_haplocoded_tag_target_geno_regs->at(i_ref_var)->chrom,
			translate_coord(reference_haplocoded_tag_target_geno_regs->at(i_ref_var)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(reference_haplocoded_tag_target_geno_regs->at(i_ref_var)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			reference_haplocoded_tag_target_geno_regs->at(i_ref_var)->name);

		for (int i_s = 0; i_s < testing_sample_ids->size(); i_s++)
		{
			fprintf(f_geno_probs, "\t%.3f;%.3f;%.3f;%.3f",
				imputed_allele_prob_per_sample_per_hap[i_s][0][0],
				imputed_allele_prob_per_sample_per_hap[i_s][0][1],
				imputed_allele_prob_per_sample_per_hap[i_s][1][0],
				imputed_allele_prob_per_sample_per_hap[i_s][1][1]);
		} // i_s loop.

		fprintf(f_geno_probs, "\n");
	} // i_ref_var loop.
	close_f(f_geno_probs, geno_probs_op_fp);

	fprintf(stderr, "Done!\n");
} // run_Imputation_Viterbi_State_Reduction_Sliding_Windows function

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
	char* geno_probs_op_fp)
{
	fprintf(stderr, "Running Viterbi with ad-hoc pruning of Viterbi entries @ %.5f of max value.\n", frac_2_max_2_prune);

	double l_target_2_center_pred_buffer_in_cM = 0.1;

	double allele_err = pow(10, -4);

	double(*XSUM)(double num1, double num2) = NULL;
	double(*XMUL)(double num1, double num2) = NULL;
	double(*XDIV)(double num1, double num2) = NULL;
	double(*XSUB)(double num1, double num2) = NULL;
	bool(*XCOMPARE)(double num1, double num2) = NULL;
	double(*XCONVERT_LIN)(double num1) = NULL;
	double ZERO = 0.0;

	if (math_mode == MATH_MODE_LOG_SPACE)
	{
		fprintf(stderr, "Math mode log.\n");
		XSUM = xlog_sum;
		XMUL = xlog_mul;
		XDIV = xlog_div;
		XSUB = xlog_sub;
		XCOMPARE = xlog_comp;
		XCONVERT_LIN = xlog;
		ZERO = xlog(0);
	}
	else if (math_mode == MATH_MODE_LIN_SPACE)
	{
		fprintf(stderr, "Math mode linear.\n");
		XSUM = lin_sum;
		XMUL = lin_mul;
		XDIV = lin_div;
		XSUB = lin_sub;
		XCOMPARE = lin_compare;
		XCONVERT_LIN = get_self;
		ZERO = 0.0;
	}

	double global_lin_scaler = exp(global_scaler_in_log);
	if (math_mode == MATH_MODE_LOG_SPACE)
	{
		global_lin_scaler = xlog(1.0);
	}

	fprintf(stderr, "Running Viterbi-ML sliding window with tag2tag distance of %.3f cMs with effective pop. size N_e=%d (Scaler: %.3f)\n", tag_2_tag_distance_cM, (int)N_e, global_scaler_in_log);

	vector<t_annot_region*>* reference_haplocoded_tag_target_geno_regs = load_variant_signal_regions_wrapper(reference_tag_target_haplocoded_genotypes_fp, ref_sample_ids_list_fp);
	vector<char*>* ref_sample_ids = buffer_file(ref_sample_ids_list_fp);
	int n_ref_haplotypes = ref_sample_ids->size() * 2;
	fprintf(stderr, "Loaded %d variant regions for %d reference samples, %d haplotypes.\n",
		reference_haplocoded_tag_target_geno_regs->size(),
		ref_sample_ids->size(),
		n_ref_haplotypes);

	// Reset the scores to 0 to indicates all untyped.
	for (int i_reg = 0; i_reg < reference_haplocoded_tag_target_geno_regs->size(); i_reg++)
	{
		reference_haplocoded_tag_target_geno_regs->at(i_reg)->score = VAR_UNTYPED_TARGET;
	} // i_reg loop.

	vector<t_annot_region*>* testing_haplocoded_tag_geno_regs = load_variant_signal_regions_wrapper(testing_tag_haplocoded_genotypes_fp, testing_sample_ids_list_fp);
	vector<char*>* testing_sample_ids = buffer_file(testing_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d testing variant regions for %d testing samples.\n", testing_haplocoded_tag_geno_regs->size(), testing_sample_ids->size());

	fprintf(stderr, "Assigning testing tag regions with the reference variants.\n");
	vector<t_annot_region*>* intersects = intersect_annot_regions(testing_haplocoded_tag_geno_regs, reference_haplocoded_tag_target_geno_regs, true);
	fprintf(stderr, "Processing %d intersections.\n", intersects->size());
	for (int i_int = 0; i_int < intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* testing_reg = int_info->src_reg;
		t_annot_region* ref_reg = int_info->dest_reg;

		// Set this reference variant as "typed". 
		ref_reg->score = VAR_TYPED_TAG;

		void** testing_reg_info = (void**)(testing_reg->data);
		void** ref_reg_info = (void**)(ref_reg->data);
		testing_reg_info[1] = ref_reg_info[0];

		delete int_info;
	} // i_int loop.
	delete_annot_regions(intersects);

	vector<t_annot_region*>* focus_regs = NULL;
	if (check_file(target_focus_reg_BED_fp))
	{
		focus_regs = load_BED(target_focus_reg_BED_fp);
		fprintf(stderr, "Found %d focus regions.\n", focus_regs->size());
	}

	// Assign the haplocoded regions.
	vector<t_annot_region*>* known_haplocoded_tag_target_geno_regs = NULL;
	if (check_file(known_tag_target_haplocoded_genotypes_fp))
	{
		fprintf(stderr, "Loading known genotypes from %s\n", known_tag_target_haplocoded_genotypes_fp);

		known_haplocoded_tag_target_geno_regs = load_variant_signal_regions_wrapper(known_tag_target_haplocoded_genotypes_fp, testing_sample_ids_list_fp);
		fprintf(stderr, "Loaded %d known tag+target variants for %d samples.\n", known_haplocoded_tag_target_geno_regs->size(), testing_sample_ids->size());

		fprintf(stderr, "Assigning the known genotype regions to reference genotypes.\n");
		intersects = intersect_annot_regions(reference_haplocoded_tag_target_geno_regs, known_haplocoded_tag_target_geno_regs, true);
		fprintf(stderr, "Processing %d known genotype reference intersects.\n", intersects->size());
		for (int i_int = 0; i_int < intersects->size(); i_int++)
		{
			t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
			t_annot_region* ref_reg = int_info->src_reg;
			t_annot_region* known_reg = int_info->dest_reg;

			void** ref_reg_info = (void**)(ref_reg->data);
			void** known_reg_info = (void**)(known_reg->data);
			ref_reg_info[1] = known_reg_info[0];

			delete int_info;
		} // i_int loop.

		delete_annot_regions(intersects);
	} // known target check.

	  // Set the imputed genotypes probabilities for untyped variants.
	for (int i_var = 0; i_var < reference_haplocoded_tag_target_geno_regs->size(); i_var++)
	{
		void** cur_ref_var_info = (void**)(reference_haplocoded_tag_target_geno_regs->at(i_var)->data);
		void** new_info = new void*[5];
		new_info[0] = cur_ref_var_info[0];
		new_info[1] = cur_ref_var_info[1];

		double*** imputed_allele_prob_per_sample_per_hap = new double**[testing_sample_ids->size() + 2];
		for (int i_s = 0; i_s < testing_sample_ids->size(); i_s++)
		{
			imputed_allele_prob_per_sample_per_hap[i_s] = new double*[2];

			imputed_allele_prob_per_sample_per_hap[i_s][0] = new double[2];
			imputed_allele_prob_per_sample_per_hap[i_s][0][0] = ZERO;
			imputed_allele_prob_per_sample_per_hap[i_s][0][1] = ZERO;

			imputed_allele_prob_per_sample_per_hap[i_s][1] = new double[2];
			imputed_allele_prob_per_sample_per_hap[i_s][1][0] = ZERO;
			imputed_allele_prob_per_sample_per_hap[i_s][1][1] = ZERO;
		} // i_s loop.
		new_info[2] = imputed_allele_prob_per_sample_per_hap;

		reference_haplocoded_tag_target_geno_regs->at(i_var)->data = new_info;
	} // i_var loop.

	int n_states = n_ref_haplotypes;
	int n_symbols = 2;

	t_restr_annot_region_list* restr_testing_haplocoded_tag_geno_regs = restructure_annot_regions(testing_haplocoded_tag_geno_regs);
	t_restr_annot_region_list* restr_reference_haplocoded_tag_target_geno_regs = restructure_annot_regions(reference_haplocoded_tag_target_geno_regs, restr_testing_haplocoded_tag_geno_regs->chr_ids);

	for (int i_chr = 0; i_chr < restr_testing_haplocoded_tag_geno_regs->chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Running Viterbi on %s\n", restr_testing_haplocoded_tag_geno_regs->chr_ids->at(i_chr));

		char cur_chr_recombination_rate_fp[1000];
		sprintf(cur_chr_recombination_rate_fp, "%s/%s.map", recombination_rate_dir, restr_testing_haplocoded_tag_geno_regs->chr_ids->at(i_chr));
		vector<t_annot_region*>* cur_chrom_recomb_regs = load_recombination_rates(cur_chr_recombination_rate_fp);
		if (cur_chrom_recomb_regs == NULL)
		{
			fprintf(stderr, "Could not load recombination rates from %s\n", cur_chr_recombination_rate_fp);
			exit(0);
		}

		sort_set_sorting_info(cur_chrom_recomb_regs, sort_regions);

		// Add one element to the beginning of the regions for the current chromosome.
		vector<t_annot_region*>* cur_chr_testing_haplocoded_tag_geno_regs = restr_testing_haplocoded_tag_geno_regs->regions_per_chrom[i_chr];
		vector<t_annot_region*>* cur_chrom_ref_tag_target_var_regs = restr_reference_haplocoded_tag_target_geno_regs->regions_per_chrom[i_chr];

		// Assign the recomb rates.
		fprintf(stderr, "Setting recombination rates using %d recombination mapped regions.\n", cur_chrom_recomb_regs->size());
		for (int i_reg = 0; i_reg < cur_chr_testing_haplocoded_tag_geno_regs->size(); i_reg++)
		{
			double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant_optimized(cur_chr_testing_haplocoded_tag_geno_regs->at(i_reg), cur_chrom_recomb_regs);
			cur_chr_testing_haplocoded_tag_geno_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
		} // i_reg loop.

		for (int i_reg = 0; i_reg < cur_chrom_ref_tag_target_var_regs->size(); i_reg++)
		{
			double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant_optimized(cur_chrom_ref_tag_target_var_regs->at(i_reg), cur_chrom_recomb_regs);
			cur_chrom_ref_tag_target_var_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
		} // i_reg loop.

		  // Get the focus regions on this chromosome.
		int focus_start = cur_chrom_ref_tag_target_var_regs->at(0)->start;
		int focus_end = cur_chrom_ref_tag_target_var_regs->back()->start;
		if (focus_regs != NULL)
		{
			vector<t_annot_region*>* cur_chr_focus_regs = get_regions_per_chromosome(focus_regs, restr_testing_haplocoded_tag_geno_regs->chr_ids->at(i_chr));
			sort(cur_chr_focus_regs->begin(), cur_chr_focus_regs->end(), sort_regions);
			focus_start = cur_chr_focus_regs->at(0)->start;
			focus_end = cur_chr_focus_regs->back()->end;
			delete cur_chr_focus_regs;
		}

		fprintf(stderr, "Focusing on variants @ %d-%d\n", focus_start, focus_end);

		vector<t_annot_region*>* prev_win_var_regs = NULL;
		for (int test_sample_i = 0; test_sample_i < testing_sample_ids->size(); test_sample_i++)
		{
			prev_win_var_regs = NULL;

			for (int win_start_var_i = 0; win_start_var_i < cur_chr_testing_haplocoded_tag_geno_regs->size(); win_start_var_i++)
			{
				double start_cM = cur_chr_testing_haplocoded_tag_geno_regs->at(win_start_var_i)->dbl_score;
				int win_end_var_i = win_start_var_i;
				while (win_end_var_i < cur_chr_testing_haplocoded_tag_geno_regs->size() &&
					fabs(start_cM - cur_chr_testing_haplocoded_tag_geno_regs->at(win_end_var_i)->dbl_score) < tag_2_tag_distance_cM)
				{
					win_end_var_i++;
				} // win_end_var_iloop.

				  // Fix end if passed over the end.
				if (win_end_var_i >= cur_chr_testing_haplocoded_tag_geno_regs->size())
				{
					win_end_var_i = cur_chr_testing_haplocoded_tag_geno_regs->size() - 1;
				}

				// Make sure we push the window forward at the end.
				if (prev_win_var_regs != NULL)
				{
					int n_prev_vars = prev_win_var_regs->size() - 2;
					double prev_win_mid_cM = (prev_win_var_regs->at(1)->dbl_score + prev_win_var_regs->at(n_prev_vars)->dbl_score) / 2;

					double cur_win_mid_cM = (cur_chr_testing_haplocoded_tag_geno_regs->at(win_start_var_i)->dbl_score + cur_chr_testing_haplocoded_tag_geno_regs->at(win_end_var_i)->dbl_score) / 2;

					// If we are to the very left of the mid point, move forward.
					if (cur_win_mid_cM < (prev_win_mid_cM + l_target_2_center_pred_buffer_in_cM / 2))
					{
						continue;
					}

					//int n_prev_vars = prev_win_var_regs->size() - 2;
					//if (prev_win_var_regs->at(n_prev_vars)->end == cur_chr_testing_haplocoded_tag_geno_regs->at(win_end_var_i)->end)
					//{
					//	fprintf(stderr, "Skipping %d-%d\n", win_start_var_i, win_end_var_i);
					//	continue;
					//}					
				}

				// Make sure we are overlapping with the focus region at the borders, at least.
				if (cur_chr_testing_haplocoded_tag_geno_regs->at(win_start_var_i)->start > (focus_end + tag_2_tag_distance_cM / 2) ||
					cur_chr_testing_haplocoded_tag_geno_regs->at(win_end_var_i)->end < (focus_start - tag_2_tag_distance_cM / 2))
				{
					continue;
				}

				double end_cM = cur_chr_testing_haplocoded_tag_geno_regs->at(win_end_var_i)->dbl_score;

				//fprintf(stderr, "Computing fore-ward scores for sample %d; hap: %d\n", test_sample_i, test_hap_i);
				vector<t_annot_region*>* cur_win_var_regs = new vector<t_annot_region*>();
				cur_win_var_regs->insert(cur_win_var_regs->end(),
					cur_chr_testing_haplocoded_tag_geno_regs->begin() + win_start_var_i,
					cur_chr_testing_haplocoded_tag_geno_regs->begin() + win_end_var_i + 1);

				// Make sure we inserted the correct elements.
				if (!t_string::compare_strings(cur_win_var_regs->at(0)->name, cur_chr_testing_haplocoded_tag_geno_regs->at(win_start_var_i)->name) ||
					!t_string::compare_strings(cur_win_var_regs->back()->name, cur_chr_testing_haplocoded_tag_geno_regs->at(win_end_var_i)->name))
				{
					fprintf(stderr, "Start and end did not copy as expected.\n");
					exit(0);
				}

				fprintf(stderr, "Window (%d variants): %s:%d-%d (%d-%d) spanning %.4f cM between %.4f-%.4f cMs\n",
					win_end_var_i - win_start_var_i + 1,
					cur_win_var_regs->at(0)->chrom, cur_win_var_regs->at(0)->start, cur_win_var_regs->back()->start,
					win_start_var_i, win_end_var_i,
					fabs(start_cM - end_cM),
					start_cM, end_cM);

				// Set the start and end state regions.
				t_annot_region* start_state_reg = NULL;
				cur_win_var_regs->insert(cur_win_var_regs->begin(), start_state_reg);

				// Set the number of variants for indexing purposes here.
				int n_vars = cur_win_var_regs->size() - 1;

				// Add the end state region.
				t_annot_region* end_state_reg = NULL;
				cur_win_var_regs->push_back(end_state_reg);

				// Allocate and initialize the forward/backward arrays.
				double*** ML_scores_per_hap = new double**[n_vars + 2];
				int*** ML_prev_state_per_hap = new int**[n_vars + 2];
				for (int var_i = 0; var_i <= n_vars + 1; var_i++)
				{
					ML_scores_per_hap[var_i] = new double*[2];
					ML_prev_state_per_hap[var_i] = new int*[2];

					for (int test_hap = 0; test_hap < 2; test_hap++)
					{
						ML_scores_per_hap[var_i][test_hap] = new double[n_states + 2];
						memset(ML_scores_per_hap[var_i][test_hap], 0, sizeof(double) * (n_states + 1));

						ML_prev_state_per_hap[var_i][test_hap] = new int[n_states + 2];
						memset(ML_prev_state_per_hap[var_i][test_hap], 0, sizeof(int) * (n_states + 1));
					}
				} // i loop.

				bool PERFORM_FULL_COMPUTE = true;

				if (PERFORM_FULL_COMPUTE ||
					prev_win_var_regs == NULL)
				{
					// Initialize the state probabilities for both haplotypes.
					for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
					{
						ML_scores_per_hap[0][0][state_i] = XCONVERT_LIN((double)1.0 / n_ref_haplotypes);
						ML_scores_per_hap[0][1][state_i] = XCONVERT_LIN((double)1.0 / n_ref_haplotypes);
					} // state_i loop.

					  // Set the currently tracked states.
					vector<int>** cur_states_per_haps = new vector<int>*[2];
					for (int hap_i = 0; hap_i < 2; hap_i++)
					{
						cur_states_per_haps[hap_i] = new vector<int>();
						for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
						{
							cur_states_per_haps[hap_i]->push_back(state_i);
						} // state_i loop.
					} // hap_i loop.

					// Start recursing over the variants; following is the main loop that compues the viterbi scores.
					for (int var_i = 1; var_i <= n_vars + 1; var_i++)
					{
						if (var_i % 10 == 0)
						{
							fprintf(stderr, "Forward: sample_i: %d/%d: var_i: %d         \r", test_sample_i, testing_sample_ids->size(), var_i);
						}

						void** cur_tag_var_info = NULL;
						char* test_sample_geno = NULL;
						char* ref_sample_geno = NULL;

						if (var_i <= n_vars)
						{
							void** cur_tag_var_info = (void**)(cur_win_var_regs->at(var_i)->data);
							test_sample_geno = (char*)(cur_tag_var_info[0]);
							ref_sample_geno = (char*)(cur_tag_var_info[1]);
						}

						/////////////////////////////////////////////////////
						// Pre-compute the transition probabilities:
						double prev_var_cM = cur_win_var_regs->at(1)->dbl_score;
						if (var_i > 1)
						{
							prev_var_cM = cur_win_var_regs->at(var_i - 1)->dbl_score;
						}

						double cur_var_cM = cur_win_var_regs->at(n_vars)->dbl_score;
						if (var_i <= n_vars)
						{
							cur_var_cM = cur_win_var_regs->at(var_i)->dbl_score;
						}

						double r_m = fabs(cur_var_cM - prev_var_cM);
						double rho_m = 4 * N_e * r_m;

						double tau_m = 1 - exp(-1 * rho_m / n_ref_haplotypes);

						double other_prob = tau_m / n_ref_haplotypes;
						double self_prob = (1 - tau_m) + (tau_m / n_ref_haplotypes);

						//fprintf(stderr, "Self prob: %.5f, other_prob: %.5f (%.4f, %.4f)\n", self_prob, other_prob, prev_var_cM, cur_var_cM);
						/////////////////////////////////////////////////////

						// Loop over all the states.
						double per_hap_max_scores[2];
						per_hap_max_scores[0] = ZERO;
						per_hap_max_scores[1] = ZERO;
						double per_hap_min_scores[2];
						per_hap_min_scores[0] = (double)(1 << 30);
						per_hap_min_scores[1] = (double)(1 << 30);
						for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
						{
							// Recurse over the previous states.
							ML_scores_per_hap[var_i][1][state_i] = ZERO;
							ML_scores_per_hap[var_i][0][state_i] = ZERO;

							//for (int prev_state_i = 0; prev_state_i < n_ref_haplotypes; prev_state_i++)
							for (int hap_i = 0; hap_i < 2; hap_i++)
							{
								for (int i_state_i = 0; i_state_i < cur_states_per_haps[hap_i]->size(); i_state_i++)
								{
									int prev_state_i = cur_states_per_haps[hap_i]->at(i_state_i);

									// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
									// Set the transition probabilities.
									double trans_prob = 0;
									if (state_i == prev_state_i)
									{
										trans_prob = self_prob;
									}
									else
									{
										trans_prob = other_prob;
									}

									// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
									// Set the emission prob: If the ref allele is matching to test, set to 0.99999, otherwise set to the error prob.
									double emit_prob = 0;

									if (var_i == n_vars + 1)
									{
										emit_prob = 1.0;
									}
									else
									{
										// There is no conditional below.
										int ref_sample_i = state_i / 2;
										int ref_hap_i = state_i % 2;
										int cur_ref_allele = get_allele_per_haplotype(ref_sample_geno[ref_sample_i], ref_hap_i);
										int cur_test_allele = get_allele_per_haplotype(test_sample_geno[test_sample_i], hap_i);

										double emit_prob_per_ref_allele[2];
										emit_prob_per_ref_allele[cur_ref_allele] = 1 - allele_err;
										emit_prob_per_ref_allele[1 - cur_ref_allele] = allele_err;
										emit_prob = emit_prob_per_ref_allele[cur_test_allele];
									}
									// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 

									double trans_emit_prob = XMUL(XCONVERT_LIN(trans_prob), XCONVERT_LIN(emit_prob));

									// Add the scaler for this position.
									trans_emit_prob = XMUL(trans_emit_prob, XCONVERT_LIN(global_lin_scaler));

									ML_scores_per_hap[var_i][hap_i][state_i] = MAX(ML_scores_per_hap[var_i][hap_i][state_i], XMUL(trans_emit_prob, ML_scores_per_hap[var_i - 1][hap_i][prev_state_i]));

									if (XCOMPARE(ML_scores_per_hap[var_i][hap_i][state_i], XMUL(trans_emit_prob, ML_scores_per_hap[var_i - 1][hap_i][prev_state_i])))
									{
										ML_prev_state_per_hap[var_i][hap_i][state_i] = prev_state_i;
									}

									if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
									{
										fprintf(stderr, "ML[%d][%d]: %.5f: P_emit=%.5f; P_trans(%d->%d): %.5f\n", var_i, state_i, ML_scores_per_hap[var_i][hap_i][state_i],
											emit_prob, prev_state_i, state_i, trans_prob);
									}
								} // prev_state_i loop.
							} // hap_i loop.

							if (math_mode == MATH_MODE_LIN_SPACE)
							{
								if (log(ML_scores_per_hap[var_i][0][state_i]) / log(10) > 50 ||
									log(ML_scores_per_hap[var_i][1][state_i]) / log(10) > 50)
								{
									fprintf(stderr, "Out of range @ %d: %.5f, %.5f\n", var_i, ML_scores_per_hap[var_i][0][state_i],
										ML_scores_per_hap[var_i][1][state_i]);
									exit(0);
								}
							}

							per_hap_max_scores[0] = MAX(per_hap_max_scores[0], ML_scores_per_hap[var_i][0][state_i]);
							per_hap_max_scores[1] = MAX(per_hap_max_scores[1], ML_scores_per_hap[var_i][1][state_i]);

							per_hap_min_scores[0] = MIN(per_hap_min_scores[0], ML_scores_per_hap[var_i][0][state_i]);
							per_hap_min_scores[1] = MIN(per_hap_min_scores[1], ML_scores_per_hap[var_i][1][state_i]);
						} // state_i loop.

						if (per_hap_min_scores[0] > 0 &&
							per_hap_min_scores[1] > 0)
						{
							fprintf(stderr, "var_i %d: max_2_min[0]=%.5f ;; max_2_min[1]=%.5f\n",
								var_i,
								XDIV(per_hap_max_scores[0], per_hap_min_scores[0]),
								XDIV(per_hap_max_scores[1], per_hap_min_scores[1]));
						}

						// Update the list of states to recurse on.
						cur_states_per_haps[0]->clear();
						cur_states_per_haps[1]->clear();

						// Select the "good" set of haplotypes.
						for (int hap_i = 0; hap_i < 2; hap_i++)
						{
							for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
							{
								// If this value is in the list, use it.
								if (ML_scores_per_hap[var_i][hap_i][state_i] > 0 &&
									XDIV(per_hap_max_scores[hap_i], ML_scores_per_hap[var_i][hap_i][state_i]) < frac_2_max_2_prune)
								{
									cur_states_per_haps[hap_i]->push_back(state_i);
								}
							} // state_i loop.

							fprintf(stderr, "Variant %d: Keeping %d/%d haplotype states\n", var_i, cur_states_per_haps[hap_i]->size(), n_ref_haplotypes);
						} // hap_i loop.
					} // var_i loop.

					for (int hap_i = 0; hap_i < 2; hap_i++)
					{
						delete cur_states_per_haps[hap_i];
					} // hap_i loop.
					delete[] cur_states_per_haps;
				} //  full compute check.

				  // Compute the total log forward and backward probabilities.
				double per_hap_total_log_ML_prob[2];
				per_hap_total_log_ML_prob[0] = ZERO;
				per_hap_total_log_ML_prob[1] = ZERO;

				int per_hap_ML_state_at_end[2];
				per_hap_ML_state_at_end[0] = -1;
				per_hap_ML_state_at_end[1] = -1;

				for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
				{
					per_hap_total_log_ML_prob[0] = MAX(per_hap_total_log_ML_prob[0], ML_scores_per_hap[n_vars + 1][0][state_i]);
					if (per_hap_total_log_ML_prob[0] == ML_scores_per_hap[n_vars + 1][0][state_i])
					{
						per_hap_ML_state_at_end[0] = state_i;
					}

					per_hap_total_log_ML_prob[1] = MAX(per_hap_total_log_ML_prob[1], ML_scores_per_hap[n_vars + 1][1][state_i]);
					if (per_hap_total_log_ML_prob[1] == ML_scores_per_hap[n_vars + 1][1][state_i])
					{
						per_hap_ML_state_at_end[1] = state_i;
					}
				} // state_i loop.

				if (per_hap_total_log_ML_prob[0] == 0 ||
					per_hap_total_log_ML_prob[1] == 0)
				{
					fprintf(stderr, "Final score is 0!\n");
					exit(0);
				}

				fprintf(stderr, "Haplotype 0 total probabilities: ML=%.5f\n", log(per_hap_total_log_ML_prob[0]));
				fprintf(stderr, "Haplotype 1 total probabilities: ML=%.5f\n", log(per_hap_total_log_ML_prob[1]));

				// Trace the states back to identify the optimal path.
				fprintf(stderr, "Tracing back the Viterbi path for each haplotype.\n");
				int** per_hap_per_variant_viterbi_haplotype = new int*[2];
				for (int hap_i = 0; hap_i < 2; hap_i++)
				{
					int cur_state = per_hap_ML_state_at_end[hap_i];
					int* per_variant_viterbi_haplotype = new int[n_vars + 2];
					for (int var_i = n_vars; var_i >= 1; var_i--)
					{
						// Set the optimal state for the current index.
						per_variant_viterbi_haplotype[var_i] = cur_state;

						// Get the new state for the previous variant.
						int new_state = ML_prev_state_per_hap[var_i][hap_i][cur_state];

						// Reset the current state for the previous variant.
						cur_state = new_state; // This is the state for the previous value..
					} // var_i loop.

					per_hap_per_variant_viterbi_haplotype[hap_i] = per_variant_viterbi_haplotype;
				} // hap_i loop.

				// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
				// Start computing the posterior probabilities for the untyped markers.
				// Now go over the untyped genotypes and compute the posterior probability of all the untyped markers.				
				int n_untyped = 0;
				double min_tar2tag_cM = MAX(0, tag_2_tag_distance_cM / 2 - 0.1);
				fprintf(stderr, "Computing the probabilities of the untyped markers out of %d markers (min dist: %.4f).\n",
					cur_chrom_ref_tag_target_var_regs->size(),
					min_tar2tag_cM);

				for (int i_ref_var = 0; i_ref_var < cur_chrom_ref_tag_target_var_regs->size(); i_ref_var++)
				{
					double target_cM = cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->dbl_score;
					double left_cM = cur_win_var_regs->at(1)->dbl_score;
					double right_cM = cur_win_var_regs->at(n_vars)->dbl_score;

					// Make sure this is an untyped variant and it is covered well by tags.
					if (cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->score == VAR_TYPED_TAG ||
						fabs(right_cM - target_cM) < min_tar2tag_cM ||
						fabs(left_cM - target_cM) < min_tar2tag_cM ||
						cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->start > focus_end ||
						cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->start < focus_start)
					{
						if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
						{
							fprintf(stderr, "Skipping %s:%d @ cM: %.4f (%.4f-%.4f) typed: %d (%.4f, %.4f)\n",
								cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->chrom, cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->start,
								target_cM, left_cM, right_cM,
								cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->score,
								fabs(right_cM - target_cM),
								fabs(left_cM - target_cM));
						}
						continue;
					}

					// This is the untyped variant's information.
					void** cur_untyped_var_ref_var_info = (void**)(cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->data);
					char* untyped_var_ref_geno = (char*)(cur_untyped_var_ref_var_info[0]);
					char* known_var_ref_geno = (char*)(cur_untyped_var_ref_var_info[1]);
					double*** imputed_allele_prob_per_sample_per_hap = (double***)(cur_untyped_var_ref_var_info[2]);

					// If there is probability assigned to this variant for this sample, do not process it.
					if (imputed_allele_prob_per_sample_per_hap[test_sample_i][0][0] != ZERO ||
						imputed_allele_prob_per_sample_per_hap[test_sample_i][0][1] != ZERO)
					{
						continue;
					}					

					n_untyped++;

					//fprintf(stderr, "@%d. untyped marker         \r", n_untyped);

					// Find the closest tag variants to this target:
					for (int var_i = 1; var_i < n_vars; var_i++)
					{
						// Find the untyped target variant in between the current tags.
						if (cur_win_var_regs->at(var_i)->start < cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->start &&
							cur_win_var_regs->at(var_i + 1)->start > cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->start)
						{
							if (__DUMP_SHMMR_GIMP_STATO_REDUCTO_MESSAGES__)
							{
								fprintf(stderr, "Found %d-%d-%d\n",
									cur_win_var_regs->at(var_i)->start,
									cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->start,
									cur_win_var_regs->at(var_i + 1)->start);
							}

							void** cur_tag_var_info = (void**)(cur_win_var_regs->at(var_i)->data);
							char* test_sample_geno = (char*)(cur_tag_var_info[0]);
							char* ref_sample_geno = (char*)(cur_tag_var_info[1]);

							// Get the self transition probability.
							double cur_cM = cur_win_var_regs->at(var_i)->dbl_score;
							double next_cM = cur_cM;
							if (var_i <= n_vars - 1)
							{
								next_cM = cur_win_var_regs->at(var_i + 1)->dbl_score;
							}

							double target_cM = cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->dbl_score;

							double left_weight = .5;
							double right_weight = .5;

							if (cur_cM != next_cM)
							{
								left_weight = fabs(target_cM - cur_cM) / (fabs(cur_cM - next_cM));
								right_weight = fabs(target_cM - next_cM) / (fabs(cur_cM - next_cM));
							}

							/*
							double r_m = fabs(cur_cM - next_cM);
							double rho_m = 4 * N_e * r_m;

							double tau_m = 1 - exp(-1 * rho_m / n_ref_haplotypes);

							double other_prob = tau_m / n_ref_haplotypes;
							double self_prob = (1 - tau_m) + (tau_m / n_ref_haplotypes);
							*/
							//////////////////////////////////////////////////////////////////////////////////////////////////

							for (int hap_i = 0; hap_i < 2; hap_i++)
							{
								int cur_hap_allele = -1;

								int var_ML_hap_i = per_hap_per_variant_viterbi_haplotype[hap_i][var_i];

								// For this untyped state in the refernce, we would like to sum all the paths that pass through here.
								int ref_sample_i = var_ML_hap_i / 2;
								int ref_hap_i = var_ML_hap_i % 2;

								int cur_ML_allele = get_allele_per_haplotype(untyped_var_ref_geno[ref_sample_i], ref_hap_i);

								// Get the haplotype next variant.
								int next_var_ML_hap_i = per_hap_per_variant_viterbi_haplotype[hap_i][var_i + 1];
								ref_sample_i = next_var_ML_hap_i / 2;
								ref_hap_i = next_var_ML_hap_i % 2;

								int next_ML_allele = get_allele_per_haplotype(untyped_var_ref_geno[ref_sample_i], ref_hap_i);

								// Set the allelic probabilities based on weights, i.e., genetic distance ratio to left and right closest variants.
								double per_allele_prob[2];
								per_allele_prob[0] = 0;
								per_allele_prob[1] = 0;
								//per_allele_prob[cur_ML_allele]++;
								per_allele_prob[cur_ML_allele] += left_weight;
								per_allele_prob[next_ML_allele] += right_weight;

								int known_geno = -1;
								if (known_var_ref_geno)
								{
									known_geno = (int)(known_var_ref_geno[test_sample_i]);
								}

								fprintf(stderr, "Var %d on sample %d: %s:%d (%s): Hap%d: [Viterbi Allele Probs: %.3f/%.3f (Tot. weight: %.3f)]; (known_geno=%d)\n",
									i_ref_var,
									test_sample_i,
									cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->chrom,
									cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->start,
									cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->name,
									hap_i,
									per_allele_prob[0], per_allele_prob[1], left_weight + right_weight,
									known_geno);

								imputed_allele_prob_per_sample_per_hap[test_sample_i][hap_i][0] = per_allele_prob[0];
								imputed_allele_prob_per_sample_per_hap[test_sample_i][hap_i][1] = per_allele_prob[1];
							} // hap_i loop.

							  // Set the inputated variant type.
							cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->score = VAR_UNTYPED_TARGET_IMPUTED;

							break; // Breaks out of position search.
						} // middle check.
					} // i_var loop.
				} // i_ref_var loop.

				for (int hap_i = 0; hap_i < 2; hap_i++)
				{
					delete[] per_hap_per_variant_viterbi_haplotype[hap_i];
				}
				delete[] per_hap_per_variant_viterbi_haplotype;

				// Allocate and initialize the forward/backward arrays.
				for (int var_i = 0; var_i <= n_vars + 1; var_i++)
				{
					for (int test_hap = 0; test_hap < 2; test_hap++)
					{
						delete[] ML_scores_per_hap[var_i][test_hap];
						delete[] ML_prev_state_per_hap[var_i][test_hap];
					} // test_hap loop.

					delete[] ML_scores_per_hap[var_i];
					delete[] ML_prev_state_per_hap[var_i];
				} // i loop.

				delete[] ML_scores_per_hap;
				delete[] ML_prev_state_per_hap;

				// Copy the current windows; this is necessary to make sure we move forward.
				if (prev_win_var_regs != NULL)
				{
					delete prev_win_var_regs;
				}

				prev_win_var_regs = cur_win_var_regs;
			} // win_start_i loop.
		} // test_sample_i loop.
	} // i_chr loop.

	  // Save the imputed genotype probabilities.
	FILE* f_geno_probs = open_f(geno_probs_op_fp, "w");
	for (int i_ref_var = 0; i_ref_var < reference_haplocoded_tag_target_geno_regs->size(); i_ref_var++)
	{
		// Make sure this is an untyped variant and it is covered well by tags.
		if (reference_haplocoded_tag_target_geno_regs->at(i_ref_var)->score != VAR_UNTYPED_TARGET_IMPUTED)
		{
			continue;
		}

		// This is the untyped variant's information.
		void** cur_untyped_var_ref_var_info = (void**)(reference_haplocoded_tag_target_geno_regs->at(i_ref_var)->data);
		double*** imputed_allele_prob_per_sample_per_hap = (double***)(cur_untyped_var_ref_var_info[2]);

		fprintf(f_geno_probs, "%s\t%d\t%d\t%s", reference_haplocoded_tag_target_geno_regs->at(i_ref_var)->chrom,
			translate_coord(reference_haplocoded_tag_target_geno_regs->at(i_ref_var)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(reference_haplocoded_tag_target_geno_regs->at(i_ref_var)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			reference_haplocoded_tag_target_geno_regs->at(i_ref_var)->name);

		for (int i_s = 0; i_s < testing_sample_ids->size(); i_s++)
		{
			fprintf(f_geno_probs, "\t%.3f;%.3f;%.3f;%.3f",
				imputed_allele_prob_per_sample_per_hap[i_s][0][0],
				imputed_allele_prob_per_sample_per_hap[i_s][0][1],
				imputed_allele_prob_per_sample_per_hap[i_s][1][0],
				imputed_allele_prob_per_sample_per_hap[i_s][1][1]);
		} // i_s loop.

		fprintf(f_geno_probs, "\n");
	} // i_ref_var loop.
	close_f(f_geno_probs, geno_probs_op_fp);

	fprintf(stderr, "Done!\n");
}


