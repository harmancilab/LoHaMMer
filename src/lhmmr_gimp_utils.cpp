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

bool __DUMP_SHMMR_GIMP_MESSAGES__ = false;

bool sort_viterbi_array_lines(double* array1_entry, double* array2_entry)
{
	for (int i = 0; i < 6; i++)
	{
		if (array1_entry[i] < array2_entry[i])
		{
			return(true);
		}
		else if (array1_entry[i] > array2_entry[i])
		{
			return(false);
		}
	} // i loop.

	return(false);
}

void compare_Viterbi_arrays(char* array1_fp, char* array2_fp)
{
	fprintf(stderr, "Loading Viterbi arrays.\n");
	vector<char*>* array1_lines = buffer_file(array1_fp);
	vector<char*>* array2_lines = buffer_file(array2_fp);
	fprintf(stderr, "Loaded %d, %d lines.\n", array1_lines->size(), array2_lines->size());
	if (array1_lines->size() != array2_lines->size())
	{
		fprintf(stderr, "# of lines are different, cannot compare.\n");
		exit(0);
	}

	vector<double*>* array1_entries = new vector<double*>();
	for (int i_l = 0; i_l < array1_lines->size(); i_l++)
	{
		double* cur_arr = new double[6];
		if(sscanf(array1_lines->at(i_l), "%lf %lf %lf %lf %lf %lf",
			&(cur_arr[0]),
			&(cur_arr[1]),
			&(cur_arr[2]),
			&(cur_arr[3]),
			&(cur_arr[4]),
			&(cur_arr[5])) != 6)
		{
			fprintf(stderr, "Could not parse: %s\n", array2_lines->at(i_l));
			exit(0);
		}

		array1_entries->push_back(cur_arr);
	} // i_l loop.
	fprintf(stderr, "Parsed %d entries in array1\n", array1_entries->size());

	vector<double*>* array2_entries = new vector<double*>();
	for (int i_l = 0; i_l < array1_lines->size(); i_l++)
	{
		double* cur_arr = new double[6];
		if (sscanf(array2_lines->at(i_l), "%lf %lf %lf %lf %lf %lf",
			&(cur_arr[0]),
			&(cur_arr[1]),
			&(cur_arr[2]),
			&(cur_arr[3]),
			&(cur_arr[4]),
			&(cur_arr[5])) != 6)
		{
			fprintf(stderr, "Could not parse: %s\n", array2_lines->at(i_l));
			exit(0);
		}

		array2_entries->push_back(cur_arr);
	} // i_l loop.
	fprintf(stderr, "Parsed %d entries in array2\n", array2_entries->size());

	sort(array1_entries->begin(), array1_entries->end(), sort_viterbi_array_lines);
	sort(array2_entries->begin(), array2_entries->end(), sort_viterbi_array_lines);

	fprintf(stderr, "Sorted entries, comparing.\n");

	for (int i_ent = 0; i_ent < array1_entries->size(); i_ent++)
	{
		if (i_ent % 100000 == 0)
		{
			fprintf(stderr, "@ %d. entry           \n", i_ent);
		}

		if (fabs(array1_entries->at(i_ent)[5] - array2_entries->at(i_ent)[5]) > 0.05)
		{
			fprintf(stderr, "%lf <-> %lf\n",
				array1_entries->at(i_ent)[5],
				array2_entries->at(i_ent)[5]);

			exit(0);
		}
	} // i_ent loop.
}

void save_foreback_matrix(double** foreback_matrix, int n_vars, int n_haplotypes, char* op_fp)
{
	FILE* f_op = open_f(op_fp, "wb");
	for (int i_var = 0; i_var < n_vars; i_var++)
	{
		fwrite(foreback_matrix[i_var], sizeof(double), n_haplotypes, f_op);
	} // i_var loop.

	close_f(f_op, op_fp);
}

void run_Imputation_ForeBack_Full_States(char* reference_tag_target_haplocoded_genotypes_fp, char* ref_sample_ids_list_fp,
							char* testing_tag_haplocoded_genotypes_fp, char* testing_sample_ids_list_fp,
							char* known_tag_target_haplocoded_genotypes_fp,
							char* recombination_rate_dir,
							double min_tar2tag_cM,
							char* fore_back_output_fp)
{
	fprintf(stderr, "Running fore-back with minimum tar2tag distance of %.3f cMs\n", min_tar2tag_cM);

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
		reference_haplocoded_tag_target_geno_regs->at(i_reg)->score = 0;
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
		ref_reg->score = 1;

		void** testing_reg_info = (void**)(testing_reg->data);
		void** ref_reg_info = (void**)(ref_reg->data);
		testing_reg_info[1] = ref_reg_info[0];

		delete int_info;
	} // i_int loop.
	delete_annot_regions(intersects);

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
	int n_states = n_ref_haplotypes;
	int n_symbols = 2;

	double N_e = pow(10, 6);
	double allele_err = pow(10, -4);

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

		// Add one element to the beginning of the regions for the current chromosome.
		vector<t_annot_region*>* cur_chr_testing_haplocoded_tag_geno_regs = restr_testing_haplocoded_tag_geno_regs->regions_per_chrom[i_chr];
		vector<t_annot_region*>* cur_chrom_ref_tag_target_var_regs = restr_reference_haplocoded_tag_target_geno_regs->regions_per_chrom[i_chr];

		// Assign the recomb rates.
		fprintf(stderr, "Setting recombination rates using %d recombination mapped regions.\n", cur_chrom_recomb_regs->size());
		for (int i_reg = 0; i_reg < cur_chr_testing_haplocoded_tag_geno_regs->size(); i_reg++)
		{
			double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant(cur_chr_testing_haplocoded_tag_geno_regs->at(i_reg), cur_chrom_recomb_regs);
			cur_chr_testing_haplocoded_tag_geno_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
		} // i_reg loop.

		for (int i_reg = 0; i_reg < cur_chrom_ref_tag_target_var_regs->size(); i_reg++)
		{
			double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant(cur_chrom_ref_tag_target_var_regs->at(i_reg), cur_chrom_recomb_regs);
			cur_chrom_ref_tag_target_var_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
		} // i_reg loop.		

		// Set the start and end state regions.
		t_annot_region* start_state_reg = get_empty_region();
		start_state_reg->chrom = t_string::copy_me_str(restr_testing_haplocoded_tag_geno_regs->chr_ids->at(i_chr));
		start_state_reg->start = 1;
		start_state_reg->end = 1;
		start_state_reg->strand = '+';
		start_state_reg->dbl_score = cur_chr_testing_haplocoded_tag_geno_regs->at(0)->dbl_score;
		cur_chr_testing_haplocoded_tag_geno_regs->insert(cur_chr_testing_haplocoded_tag_geno_regs->begin(), start_state_reg);

		// Set the number of variants for indexing purposes here.
		int n_vars = cur_chr_testing_haplocoded_tag_geno_regs->size() - 1;

		// Add the end state region.
		t_annot_region* end_state_reg = get_empty_region();
		end_state_reg->chrom = t_string::copy_me_str(restr_testing_haplocoded_tag_geno_regs->chr_ids->at(i_chr));
		end_state_reg->start = 1;
		end_state_reg->end = 1;
		end_state_reg->strand = '+';
		end_state_reg->dbl_score = cur_chr_testing_haplocoded_tag_geno_regs->back()->dbl_score;
		cur_chr_testing_haplocoded_tag_geno_regs->push_back(end_state_reg);

		for (int test_sample_i = 0; test_sample_i < testing_sample_ids->size(); test_sample_i++)
		{
			fprintf(stderr, "Computing forward probabilities for sample %d\n", test_sample_i);

			for (int test_hap_i = 0; test_hap_i < 2; test_hap_i++)
			{
				fprintf(stderr, "Computing fore-ward scores for sample %d; hap: %d\n", test_sample_i, test_hap_i);

				// Allocate and initialize the forward array.
				double** fore_scores = new double*[n_vars + 2];
				for (int var_i = 0; var_i <= n_vars + 1; var_i++)
				{
					fore_scores[var_i] = new double[n_states + 2];
					memset(fore_scores[var_i], 0, sizeof(double) * (n_states + 1));
				} // i loop.

				// Initialize the state probabilities.
				for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
				{
					fore_scores[0][state_i] = xlog((double)1.0 / n_ref_haplotypes);
				} // state_i loop.

				// Start recursing over the variants.
				for (int var_i = 1; var_i <= n_vars+1; var_i++)
				{
					fprintf(stderr, "Forward: sample_i: %d/%d [%d]: var_i: %d         \r", test_sample_i, testing_sample_ids->size(), test_hap_i, var_i);

					void** cur_tag_var_info = (void**)(cur_chr_testing_haplocoded_tag_geno_regs->at(var_i)->data);

					char* test_sample_geno = NULL;
					char* ref_sample_geno = NULL;
					if (var_i <= n_vars)
					{
						test_sample_geno = (char*)(cur_tag_var_info[0]);
						ref_sample_geno = (char*)(cur_tag_var_info[1]);
					}

					/////////////////////////////////////////////////////
					// Pre-compute the transition probabilities:
					double prev_var_cM = 0;
					if (var_i > 1)
					{
						prev_var_cM = cur_chr_testing_haplocoded_tag_geno_regs->at(var_i - 1)->dbl_score;
					}

					double cur_var_cM = cur_chr_testing_haplocoded_tag_geno_regs->at(var_i)->dbl_score;

					double r_m = fabs(cur_var_cM - prev_var_cM);
					double rho_m = 4 * N_e * r_m;

					double tau_m = 1 - exp(-1 * rho_m / n_ref_haplotypes);

					double other_prob = tau_m / n_ref_haplotypes;
					double self_prob = (1 - tau_m) + (tau_m / n_ref_haplotypes);

					//fprintf(stderr, "Self prob: %.5f, other_prob: %.5f (%.4f, %.4f)\n", self_prob, other_prob, prev_var_cM, cur_var_cM);
					/////////////////////////////////////////////////////

					// Loop over all the states.
					for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
					{
						// Recurse over the previous states.
						fore_scores[var_i][state_i] = xlog(0);
						for (int prev_state_i = 0; prev_state_i < n_ref_haplotypes; prev_state_i++)
						{
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
								int ref_sample_i = state_i / 2;
								int ref_hap_i = state_i % 2;
								int cur_ref_allele = get_allele_per_haplotype(ref_sample_geno[ref_sample_i], ref_hap_i);
								int cur_test_allele = get_allele_per_haplotype(test_sample_geno[test_sample_i], test_hap_i);
								if (cur_ref_allele == cur_test_allele)
								{
									emit_prob = 1 - allele_err;
								}
								else
								{
									emit_prob = allele_err;
								}
							}
							// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 

							double trans_emit_prob = xlog_mul(xlog(trans_prob), xlog(emit_prob));

							fore_scores[var_i][state_i] = xlog_sum(fore_scores[var_i][state_i], xlog_mul(trans_emit_prob, fore_scores[var_i - 1][prev_state_i]));

							if (__DUMP_SHMMR_GIMP_MESSAGES__)
							{
								//fprintf(stderr, "fore[%d][%d]: %.5f: P_emit(%d, %d)=%.5f; P_trans(%d->%d): %.5f\n", var_i, state_i, fore_scores[var_i][state_i], cur_ref_allele, cur_test_allele, emit_prob, prev_state_i, state_i, trans_prob);
							}
						} // prev_state_i loop.
					} // state_i loop.
				} // var_i loop.

				// Move the backward scores.
				fprintf(stderr, "Computing back-ward scores for sample %d; hap: %d\n", test_sample_i, test_hap_i);

				// Compute the backward probabilities.
				// Allocate and initialize the backward array.
				double** back_scores = new double*[n_vars + 2];
				for (int var_i = 0; var_i <= n_vars + 1; var_i++)
				{
					back_scores[var_i] = new double[n_states + 2];
					memset(back_scores[var_i], 0, sizeof(double) * (n_states + 2));
				} // i loop.

				for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
				{
					back_scores[n_vars+1][state_i] = xlog((double)1.0 / n_ref_haplotypes);
				} // state_i loop.

				// Compute backward probabilities: Our backward is a little different to be able to perform below marginalization.
				for (int var_i = n_vars; var_i >= 0; var_i--)
				{
					fprintf(stderr, "Backward: sample_i: %d/%d [%d]: var_i: %d         \r", test_sample_i, testing_sample_ids->size(), test_hap_i, var_i);

					// Do not set the genotypes for the after-last index.
					char* test_sample_geno = NULL;
					char* ref_sample_geno = NULL;
					if (var_i > 0)
					{
						void** cur_tag_var_info = (void**)(cur_chr_testing_haplocoded_tag_geno_regs->at(var_i)->data);
						test_sample_geno = (char*)(cur_tag_var_info[0]);
						ref_sample_geno = (char*)(cur_tag_var_info[1]);
					}

					// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
					// Set the transition probabilities.
					double cur_cM = cur_chr_testing_haplocoded_tag_geno_regs->at(var_i)->dbl_score;
					double next_cM = cur_cM;
					if (var_i <= n_vars-1)
					{
						next_cM = cur_chr_testing_haplocoded_tag_geno_regs->at(var_i + 1)->dbl_score;
					}

					double r_m = fabs(cur_cM - next_cM);
					double rho_m = 4 * N_e * r_m;

					double tau_m = 1 - exp(-1 * rho_m / n_ref_haplotypes);

					double other_prob = tau_m / n_ref_haplotypes;
					double self_prob = (1 - tau_m) + (tau_m / n_ref_haplotypes);
					
					// Go over all the reference haplotypes.
					for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
					{
						back_scores[var_i][state_i] = xlog(0);
						for (int next_state_i = 0; next_state_i < n_ref_haplotypes; next_state_i++)
						{
							// Same state ---, different state ---.
							double trans_prob = 0;
							if (state_i == next_state_i)
							{
								trans_prob = self_prob;
							}
							else
							{
								trans_prob = other_prob;
							}

							// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
							// Set the emission prob: If the ref allele is matching to test, set to 0.99999, otherwise set to the error prob.
							// This is the emission of the current allele.
							double emit_prob = 0; // 0.999 for error, 
							if (var_i == 0)
							{
								emit_prob = 1.0;
							}
							else
							{
								int ref_sample_i = state_i / 2;
								int ref_hap_i = state_i % 2;
								int ref_allele = get_allele_per_haplotype(ref_sample_geno[ref_sample_i], ref_hap_i);
								int test_allele = get_allele_per_haplotype(test_sample_geno[test_sample_i], test_hap_i);

								if (ref_allele == test_allele)
								{
									emit_prob = (double)(1.0) - allele_err;
								}
								else
								{
									emit_prob = allele_err;
								}
							}
							// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 

							double trans_emit_prob = xlog_mul(xlog(trans_prob), xlog(emit_prob));

							// Update the backward score.
							back_scores[var_i][state_i] = xlog_sum(back_scores[var_i][state_i], xlog_mul(trans_emit_prob, back_scores[var_i + 1][next_state_i]));
						} // prev_state_i loop.
					} // state_i loop.
				} // var_i loop.

				// Compute the total log forward and backward probabilities.
				double total_log_fore_prob = xlog(0.0);
				double total_log_back_prob = xlog(0.0);
				for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
				{
					total_log_fore_prob = xlog_sum(total_log_fore_prob, fore_scores[n_vars+1][state_i]);
					total_log_back_prob = xlog_sum(total_log_back_prob, back_scores[0][state_i]);
				} // state_i loop.
				fprintf(stderr, "Total probabilities: fore=%.5f ;; back=%.5f\n", total_log_fore_prob, total_log_back_prob);

				// Now go over the untyped genotypes and compute the posterior probability of all the untyped markers.
				fprintf(stderr, "Computing the probabilities of the untyped markers out of %d markers.\n", cur_chrom_ref_tag_target_var_regs->size());
				int n_untyped = 0;
				for (int i_ref_var = 0; i_ref_var < cur_chrom_ref_tag_target_var_regs->size(); i_ref_var++)
				{
					double target_cM = cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->dbl_score;
					double left_cM = cur_chr_testing_haplocoded_tag_geno_regs->at(0)->dbl_score;
					double right_cM = cur_chr_testing_haplocoded_tag_geno_regs->back()->dbl_score;

					// Make sure this is an untyped variant and it is covered well by tags.
					if (cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->score == 1 ||
						fabs(right_cM - target_cM) < min_tar2tag_cM || 
						fabs(left_cM - target_cM) < min_tar2tag_cM)
					{
						continue;
					}

					n_untyped++;

					fprintf(stderr, "@%d. untyped marker         \r", n_untyped);

					// Find the closest tag variants to this target: This is the fore-back coordinates.
					for (int var_i = 1; var_i <= n_vars; var_i++)
					{
						//if (!t_string::compare_strings("rs141571109_G_A_0.0109824", cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->name))
						//{
						//	continue;
						//}

						if (cur_chr_testing_haplocoded_tag_geno_regs->at(var_i)->start < cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->start &&
							cur_chr_testing_haplocoded_tag_geno_regs->at(var_i+1)->start > cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->start)
						{
							fprintf(stderr, "Found %d-%d-%d\n", 
									cur_chr_testing_haplocoded_tag_geno_regs->at(var_i)->start, 
									cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->start,
									cur_chr_testing_haplocoded_tag_geno_regs->at(var_i + 1)->start);			

							void** cur_tag_var_info = (void**)(testing_haplocoded_tag_geno_regs->at(var_i)->data);
							char* test_sample_geno = (char*)(cur_tag_var_info[0]);
							char* ref_sample_geno = (char*)(cur_tag_var_info[1]);

							// This is the untyped variant's information.
							void** cur_untyped_var_ref_var_info = (void**)(cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->data);
							char* untyped_var_ref_geno = (char*)(cur_untyped_var_ref_var_info[0]);
							char* known_var_ref_geno = (char*)(cur_untyped_var_ref_var_info[1]);

							// Get the self transition probability.
							double cur_cM = cur_chr_testing_haplocoded_tag_geno_regs->at(var_i)->dbl_score;
							double next_cM = cur_cM;
							if (var_i <= n_vars - 1)
							{
								next_cM = cur_chr_testing_haplocoded_tag_geno_regs->at(var_i + 1)->dbl_score;
							}
							double r_m = fabs(cur_cM - next_cM);
							double rho_m = 4 * N_e * r_m;

							double tau_m = 1 - exp(-1 * rho_m / n_ref_haplotypes);

							double other_prob = tau_m / n_ref_haplotypes;
							double self_prob = (1 - tau_m) + (tau_m / n_ref_haplotypes);
							//////////////////////////////////////////////////////////////////////////////////////////////////

							double per_allele_probs[2];
							per_allele_probs[0] = xlog(0);
							per_allele_probs[1] = xlog(0);
							for (int untyped_state_i = 0; untyped_state_i < n_ref_haplotypes; untyped_state_i++)
							{
								// For this untyped state in the refernce, we would like to sum all the paths that pass through here.
								int ref_sample_i = untyped_state_i / 2;
								int ref_hap_i = untyped_state_i % 2;
								int cur_untyped_ref_allele = get_allele_per_haplotype(untyped_var_ref_geno[ref_sample_i], ref_hap_i);

								double cur_untyped_state_fore_prob = xlog_mul(self_prob, fore_scores[var_i][untyped_state_i]);

								double cur_untyped_state_back_prob = back_scores[var_i + 1][untyped_state_i];

								per_allele_probs[cur_untyped_ref_allele] = xlog_sum(per_allele_probs[cur_untyped_ref_allele], xlog_mul(cur_untyped_state_fore_prob, cur_untyped_state_back_prob));
							} // untyped_state_i loop.

							  // We went through all the states, normalize the probabilities.
							per_allele_probs[0] = xlog_div(per_allele_probs[0], total_log_fore_prob);
							per_allele_probs[1] = xlog_div(per_allele_probs[1], total_log_fore_prob);

							////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
							// FOLLOWING DOES NOT WORK.
							//double per_allele_probs[2];
							//per_allele_probs[0] = xlog(0);
							//per_allele_probs[1] = xlog(0);
							//for (int untyped_state_i = 0; untyped_state_i < n_ref_haplotypes; untyped_state_i++)
							//{
							//	// For this untyped state in the refernce, we would like to sum all the paths that pass through here.
							//	int ref_sample_i = untyped_state_i / 2;
							//	int ref_hap_i = untyped_state_i % 2;
							//	int cur_untyped_ref_allele = get_allele_per_haplotype(untyped_var_ref_geno[ref_sample_i], ref_hap_i);

							//	double cur_untyped_state_fore_prob = xlog(0);
							//	for (int prev_state_i = 0; prev_state_i < n_ref_haplotypes; prev_state_i++)
							//	{
							//		double fore_trans_prob = xlog(0);
							//		// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
							//		double r_m = fabs(cur_chr_testing_haplocoded_tag_geno_regs->at(var_i)->dbl_score - 
							//							cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->dbl_score);
							//		double rho_m = 4 * N_e * r_m;

							//		double tau_m = 1 - exp(-1 * rho_m / n_ref_haplotypes);

							//		double other_prob = tau_m / n_ref_haplotypes;
							//		double self_prob = (1 - tau_m) + (tau_m / n_ref_haplotypes);

							//		if (prev_state_i == untyped_state_i)
							//		{
							//			fore_trans_prob = self_prob;
							//		}
							//		else
							//		{
							//			fore_trans_prob = other_prob;
							//		}

							//		cur_untyped_state_fore_prob = xlog_sum(cur_untyped_state_fore_prob, xlog_mul(fore_trans_prob, fore_scores[var_i][prev_state_i]));
							//	} // prev_state_i loop.

							//	double cur_untyped_state_back_prob = xlog(0);
							//	for(int next_state_i = 0; next_state_i < n_ref_haplotypes; next_state_i++)
							//	{
							//		// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
							//		double back_trans_prob = xlog(0);

							//		// Next tag and current reference distance.
							//		double r_m = fabs(cur_chr_testing_haplocoded_tag_geno_regs->at(var_i+1)->dbl_score - 
							//							cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->dbl_score);
							//		double rho_m = 4 * N_e * r_m;

							//		double tau_m = 1 - exp(-1 * rho_m / n_ref_haplotypes);

							//		double other_prob = tau_m / n_ref_haplotypes;
							//		double self_prob = (1 - tau_m) + (tau_m / n_ref_haplotypes);

							//		if (next_state_i == untyped_state_i)
							//		{
							//			back_trans_prob = self_prob;
							//		}
							//		else
							//		{
							//			back_trans_prob = other_prob;
							//		}

							//		// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 

							//		cur_untyped_state_back_prob = xlog_sum(cur_untyped_state_back_prob, xlog_mul(back_trans_prob, back_scores[var_i+1][next_state_i]));
							//	} // next_state_i loop.

							//	per_allele_probs[cur_untyped_ref_allele] = xlog_sum(per_allele_probs[cur_untyped_ref_allele], xlog_mul(cur_untyped_state_fore_prob, cur_untyped_state_back_prob));
							//} // untyped_state_i loop.

							//// We went through all the states, normalize the probabilities.
							//per_allele_probs[0] = xlog_div(per_allele_probs[0], total_log_fore_prob);
							//per_allele_probs[1] = xlog_div(per_allele_probs[1], total_log_fore_prob);
							////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

							int known_geno = -1;
							if (known_var_ref_geno)
							{
								known_geno = (int)(known_var_ref_geno[test_sample_i]);
							}

							fprintf(stderr, "Var %d: %s:%d (%s) (hap: %d): %.4f, %.4f (known_geno=%d)\n", 
									i_ref_var,
									cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->chrom,
									cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->start,
									cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->name, 
									test_hap_i, 
									per_allele_probs[0], per_allele_probs[1], known_geno);
							break; // Breaks out of position search.
 						} // middle check.
					} // i_var loop.
				} // i_ref_var loop.
			} // test_hap_i loop.
		} // test_sample_i loop.
	} // i_chr loop.

	fprintf(stderr, "Done!\n");
}

void run_Imputation_ForeBack_Full_States_Sliding_Windows_Math_Mode(char* reference_tag_target_haplocoded_genotypes_fp, char* ref_sample_ids_list_fp,
	char* testing_tag_haplocoded_genotypes_fp, char* testing_sample_ids_list_fp,
	char* known_tag_target_haplocoded_genotypes_fp,
	char* recombination_rate_dir,
	double tag_2_tag_distance_cM,
	double N_e,
	char* target_focus_reg_BED_fp,
	bool math_mode,
	double global_scaler_in_log,
	char* geno_probs_op_fp)
{
	double l_target_2_center_pred_buffer_in_cM = 0.1;

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

	double global_lin_scaler = exp(global_scaler_in_log);
	if (math_mode == MATH_MODE_LOG_SPACE)
	{
		global_lin_scaler = xlog(1.0);
	}

	fprintf(stderr, "Running fore-back sliding window with tag2tag distance of %.3f cMs with effective pop. size N_e=%d (Math_Mode: %d; Scaler: %.3f)\n", tag_2_tag_distance_cM, (int)N_e, math_mode, global_scaler_in_log);

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

		// Add one element to the beginning of the regions for the current chromosome.
		vector<t_annot_region*>* cur_chr_testing_haplocoded_tag_geno_regs = restr_testing_haplocoded_tag_geno_regs->regions_per_chrom[i_chr];
		vector<t_annot_region*>* cur_chrom_ref_tag_target_var_regs = restr_reference_haplocoded_tag_target_geno_regs->regions_per_chrom[i_chr];

		// Assign the recomb rates.
		fprintf(stderr, "Setting recombination rates using %d recombination mapped regions.\n", cur_chrom_recomb_regs->size());
		for (int i_reg = 0; i_reg < cur_chr_testing_haplocoded_tag_geno_regs->size(); i_reg++)
		{
			double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant(cur_chr_testing_haplocoded_tag_geno_regs->at(i_reg), cur_chrom_recomb_regs);
			cur_chr_testing_haplocoded_tag_geno_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
		} // i_reg loop.

		for (int i_reg = 0; i_reg < cur_chrom_ref_tag_target_var_regs->size(); i_reg++)
		{
			double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant(cur_chrom_ref_tag_target_var_regs->at(i_reg), cur_chrom_recomb_regs);
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

				////////////////////////////////////////////////////////////////////////////////////
				// BELOW DOES NOT WORK YET.
				// Copy the previous window if they exist.
				bool PERFORM_FULL_COMPUTE = true;

				if (PERFORM_FULL_COMPUTE ||
					prev_win_var_regs == NULL)
				{
					// Initialize the state probabilities for both haplotypes.
					for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
					{
						fore_scores_per_hap[0][0][state_i] = XCONVERT_LIN((double)1.0 / n_ref_haplotypes);
						fore_scores_per_hap[0][1][state_i] = XCONVERT_LIN((double)1.0 / n_ref_haplotypes);
					} // state_i loop.

					  // Start recursing over the variants.
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
						for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
						{
							// Recurse over the previous states.
							fore_scores_per_hap[var_i][1][state_i] = ZERO;
							fore_scores_per_hap[var_i][0][state_i] = ZERO;
							for (int prev_state_i = 0; prev_state_i < n_ref_haplotypes; prev_state_i++)
							{
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
								double emit_prob0 = ZERO;
								double emit_prob1 = ZERO;

								if (var_i == n_vars + 1)
								{
									emit_prob0 = XCONVERT_LIN(1.0);
									emit_prob1 = XCONVERT_LIN(1.0);
								}
								else
								{
									int ref_sample_i = state_i / 2;
									int ref_hap_i = state_i % 2;
									int cur_ref_allele = get_allele_per_haplotype(ref_sample_geno[ref_sample_i], ref_hap_i);
									int cur_test_allele0 = get_allele_per_haplotype(test_sample_geno[test_sample_i], 0);
									int cur_test_allele1 = get_allele_per_haplotype(test_sample_geno[test_sample_i], 1);

									double emit_prob_per_ref_allele[2];
									emit_prob_per_ref_allele[cur_ref_allele] = 1 - allele_err;
									emit_prob_per_ref_allele[1 - cur_ref_allele] = allele_err;
									emit_prob0 = emit_prob_per_ref_allele[cur_test_allele0];
									emit_prob1 = emit_prob_per_ref_allele[cur_test_allele1];
								}
								// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 

								double trans_emit_prob0 = XMUL(XCONVERT_LIN(trans_prob), XCONVERT_LIN(emit_prob0));
								double trans_emit_prob1 = XMUL(XCONVERT_LIN(trans_prob), XCONVERT_LIN(emit_prob1));

								// Apply the scaler.
								trans_emit_prob0 = XMUL(trans_emit_prob0, XCONVERT_LIN(global_lin_scaler));
								trans_emit_prob1 = XMUL(trans_emit_prob1, XCONVERT_LIN(global_lin_scaler));

								fore_scores_per_hap[var_i][0][state_i] = XSUM(fore_scores_per_hap[var_i][0][state_i], 
																				XMUL(trans_emit_prob0, fore_scores_per_hap[var_i - 1][0][prev_state_i]));

								fore_scores_per_hap[var_i][1][state_i] = XSUM(fore_scores_per_hap[var_i][1][state_i], 
																				XMUL(trans_emit_prob1, fore_scores_per_hap[var_i - 1][1][prev_state_i]));

								if (__DUMP_SHMMR_GIMP_MESSAGES__)
								{
									//fprintf(stderr, "fore[%d][%d]: %.5f: P_emit(%d, %d)=%.5f; P_trans(%d->%d): %.5f\n", var_i, state_i, fore_scores[var_i][state_i], cur_ref_allele, cur_test_allele, emit_prob, prev_state_i, state_i, trans_prob);
								}
							} // prev_state_i loop.
						} // state_i loop.
					} // var_i loop.

					  // Move the backward scores.
					fprintf(stderr, "Computing back-ward scores for sample %d\n", test_sample_i);

					for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
					{
						back_scores_per_hap[n_vars + 1][0][state_i] = XCONVERT_LIN((double)1.0 / n_ref_haplotypes);
						back_scores_per_hap[n_vars + 1][1][state_i] = XCONVERT_LIN((double)1.0 / n_ref_haplotypes);
					} // state_i loop.

					  // Compute backward probabilities: Our backward is a little different to be able to perform below marginalization.
					for (int var_i = n_vars; var_i >= 0; var_i--)
					{
						if (var_i % 10 == 0)
						{
							fprintf(stderr, "Backward: sample_i: %d/%d: var_i: %d         \r", test_sample_i, testing_sample_ids->size(), var_i);
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

						// Go over all the reference haplotypes.
						for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
						{
							back_scores_per_hap[var_i][0][state_i] = ZERO;
							back_scores_per_hap[var_i][1][state_i] = ZERO;
							for (int next_state_i = 0; next_state_i < n_ref_haplotypes; next_state_i++)
							{
								// Same state ---, different state ---.
								double trans_prob = 0;
								if (state_i == next_state_i)
								{
									trans_prob = self_prob;
								}
								else
								{
									trans_prob = other_prob;
								}

								// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
								// Set the emission prob: If the ref allele is matching to test, set to 0.99999, otherwise set to the error prob.
								// This is the emission of the current allele.
								double emit_prob0 = ZERO; // 0.999 for error, 
								double emit_prob1 = ZERO; // 0.999 for error, 
								if (var_i == 0)
								{
									emit_prob0 = XCONVERT_LIN(1.0);
									emit_prob1 = XCONVERT_LIN(1.0);
								}
								else
								{
									int ref_sample_i = state_i / 2;
									int ref_hap_i = state_i % 2;
									int cur_ref_allele = get_allele_per_haplotype(ref_sample_geno[ref_sample_i], ref_hap_i);
									//int test_allele = get_allele_per_haplotype(test_sample_geno[test_sample_i], test_hap_i);

									int cur_test_allele0 = get_allele_per_haplotype(test_sample_geno[test_sample_i], 0);
									int cur_test_allele1 = get_allele_per_haplotype(test_sample_geno[test_sample_i], 1);

									double emit_prob_per_ref_allele[2];
									emit_prob_per_ref_allele[cur_ref_allele] = 1 - allele_err;
									emit_prob_per_ref_allele[1 - cur_ref_allele] = allele_err;
									emit_prob0 = emit_prob_per_ref_allele[cur_test_allele0];
									emit_prob1 = emit_prob_per_ref_allele[cur_test_allele1];
								}
								// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 

								double trans_emit_prob0 = XMUL(XCONVERT_LIN(trans_prob), XCONVERT_LIN(emit_prob0));
								double trans_emit_prob1 = XMUL(XCONVERT_LIN(trans_prob), XCONVERT_LIN(emit_prob1));

								// Apply the scaler.
								trans_emit_prob0 = XMUL(trans_emit_prob0, XCONVERT_LIN(global_lin_scaler));
								trans_emit_prob1 = XMUL(trans_emit_prob1, XCONVERT_LIN(global_lin_scaler));

								// Update the backward score.
								back_scores_per_hap[var_i][0][state_i] = XSUM(back_scores_per_hap[var_i][0][state_i], 
																			XMUL(trans_emit_prob0, back_scores_per_hap[var_i + 1][0][next_state_i]));

								back_scores_per_hap[var_i][1][state_i] = XSUM(back_scores_per_hap[var_i][1][state_i], 
																			XMUL(trans_emit_prob1, back_scores_per_hap[var_i + 1][1][next_state_i]));
							} // prev_state_i loop.
						} // state_i loop.
					} // var_i loop.
				} //  full compute check.

				  // Compute the total log forward and backward probabilities.
				double per_hap_total_log_fore_prob[2];
				per_hap_total_log_fore_prob[0] = XCONVERT_LIN(0.0);
				per_hap_total_log_fore_prob[1] = XCONVERT_LIN(0.0);
				double per_hap_total_log_back_prob[2];
				per_hap_total_log_back_prob[0] = XCONVERT_LIN(0.0);
				per_hap_total_log_back_prob[1] = XCONVERT_LIN(0.0);

				for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
				{
					per_hap_total_log_fore_prob[0] = XSUM(per_hap_total_log_fore_prob[0], fore_scores_per_hap[n_vars + 1][0][state_i]);
					per_hap_total_log_back_prob[0] = XSUM(per_hap_total_log_back_prob[0], back_scores_per_hap[0][0][state_i]);

					per_hap_total_log_fore_prob[1] = XSUM(per_hap_total_log_fore_prob[1], fore_scores_per_hap[n_vars + 1][1][state_i]);
					per_hap_total_log_back_prob[1] = XSUM(per_hap_total_log_back_prob[1], back_scores_per_hap[0][1][state_i]);
				} // state_i loop.

				fprintf(stderr, "Haplotype 0 total probabilities: fore=%.5f ;; back=%.5f\n", XCONVERT_2_LOG(per_hap_total_log_fore_prob[0]), XCONVERT_2_LOG(per_hap_total_log_back_prob[0]));
				fprintf(stderr, "Haplotype 1 total probabilities: fore=%.5f ;; back=%.5f\n", XCONVERT_2_LOG(per_hap_total_log_fore_prob[1]), XCONVERT_2_LOG(per_hap_total_log_back_prob[1]));

				// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
				// Start computing the posterior probabilities for the untyped markers.
				// Now go over the untyped genotypes and compute the posterior probability of all the untyped markers.				
				int n_untyped = 0;
				double min_tar2tag_cM = tag_2_tag_distance_cM / 2 - l_target_2_center_pred_buffer_in_cM;
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
						if (__DUMP_SHMMR_GIMP_MESSAGES__)
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

					fprintf(stderr, "@%d. untyped marker         \r", n_untyped);

					// Find the closest tag variants to this target: This is the fore-back coordinates.
					for (int var_i = 1; var_i < n_vars; var_i++)
					{
						// Find the untyped target variant in between the current tags.
						if (cur_win_var_regs->at(var_i)->start < cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->start &&
							cur_win_var_regs->at(var_i + 1)->start > cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->start)
						{
							fprintf(stderr, "Found %d-%d-%d\n",
								cur_win_var_regs->at(var_i)->start,
								cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->start,
								cur_win_var_regs->at(var_i + 1)->start);

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
								for (int untyped_state_i = 0; untyped_state_i < n_ref_haplotypes; untyped_state_i++)
								{
									int untyped_state_j = untyped_state_i;
									//for (int untyped_state_j = 0; untyped_state_j < n_ref_haplotypes; untyped_state_j++)
									{
										// For this untyped state in the refernce, we would like to sum all the paths that pass through here.
										int ref_sample_i = untyped_state_i / 2;
										int ref_hap_i = untyped_state_i % 2;
										int cur_untyped_ref_allele = get_allele_per_haplotype(untyped_var_ref_geno[ref_sample_i], ref_hap_i);

										double cur_untyped_state_fore_prob = XMUL(other_prob, fore_scores_per_hap[var_i][hap_i][untyped_state_i]);;
										if (untyped_state_i == untyped_state_j)
										{
											cur_untyped_state_fore_prob = XMUL(self_prob, fore_scores_per_hap[var_i][hap_i][untyped_state_i]);
										}

										double cur_untyped_state_back_prob = back_scores_per_hap[var_i + 1][hap_i][untyped_state_j];

										per_allele_probs[cur_untyped_ref_allele] = XSUM(per_allele_probs[cur_untyped_ref_allele], XMUL(cur_untyped_state_fore_prob, cur_untyped_state_back_prob));
									} // untyped_state_j loop/block.
								} // untyped_state_i loop.

								  // We went through all the states, normalize the probabilities.
								per_allele_probs[0] = XDIV(per_allele_probs[0], per_hap_total_log_fore_prob[hap_i]);
								per_allele_probs[1] = XDIV(per_allele_probs[1], per_hap_total_log_fore_prob[hap_i]);

								int known_geno = -1;
								if (known_var_ref_geno)
								{
									known_geno = (int)(known_var_ref_geno[test_sample_i]);
								}

								fprintf(stderr, "Var %d on sample %d: %s:%d (%s): Hap%d: [0:%.4f, 1:%.4f]; (known_geno=%d)\n",
									i_ref_var,
									test_sample_i,
									cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->chrom,
									cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->start,
									cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->name,
									hap_i,
									XCONVERT_2_LOG(per_allele_probs[0]), 
									XCONVERT_2_LOG(per_allele_probs[1]), 
									known_geno);

								imputed_allele_prob_per_sample_per_hap[test_sample_i][hap_i][0] = per_allele_probs[0];
								imputed_allele_prob_per_sample_per_hap[test_sample_i][hap_i][1] = per_allele_probs[1];
							} // hap_i loop.

							  // Set the inputated variant type.
							cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->score = VAR_UNTYPED_TARGET_IMPUTED;

							break; // Breaks out of position search.
						} // middle check.
					} // i_var loop.
				} // i_ref_var loop.

				  // Copy the current windows.
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

void get_n_local_window_tag_variants_per_target(char* testing_tag_regions_BED_fp,
											char* recombination_rate_dir,
											double tag_2_tag_distance_cM,
											char* op_fp)
{
	fprintf(stderr, "Getting local tag variant stats with tag2tag distance of %.3f cMs with effective pop. size.\n", tag_2_tag_distance_cM);

	vector<t_annot_region*>* testing_haplocoded_tag_geno_regs = load_BED(testing_tag_regions_BED_fp);
	fprintf(stderr, "Loaded %d testing variant regions.\n", testing_haplocoded_tag_geno_regs->size());

	t_restr_annot_region_list* restr_testing_haplocoded_tag_geno_regs = restructure_annot_regions(testing_haplocoded_tag_geno_regs);

	FILE* f_op = open_f(op_fp, "w");
	for (int i_chr = 0; i_chr < restr_testing_haplocoded_tag_geno_regs->chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Processing regions on %s\n", restr_testing_haplocoded_tag_geno_regs->chr_ids->at(i_chr));

		char cur_chr_recombination_rate_fp[1000];
		sprintf(cur_chr_recombination_rate_fp, "%s/%s.map", recombination_rate_dir, restr_testing_haplocoded_tag_geno_regs->chr_ids->at(i_chr));
		vector<t_annot_region*>* cur_chrom_recomb_regs = load_recombination_rates(cur_chr_recombination_rate_fp);
		if (cur_chrom_recomb_regs == NULL)
		{
			fprintf(stderr, "Could not load recombination rates from %s\n", cur_chr_recombination_rate_fp);
			continue;
		}

		sort_set_sorting_info(cur_chrom_recomb_regs, sort_regions);

		// Add one element to the beginning of the regions for the current chromosome.
		vector<t_annot_region*>* cur_chr_testing_haplocoded_tag_geno_regs = restr_testing_haplocoded_tag_geno_regs->regions_per_chrom[i_chr];
		//vector<t_annot_region*>* cur_chrom_ref_tag_target_var_regs = restr_reference_haplocoded_tag_target_geno_regs->regions_per_chrom[i_chr];

		// Assign the recomb rates.
		fprintf(stderr, "Setting recombination rates using %d recombination mapped regions.\n", cur_chrom_recomb_regs->size());
		for (int i_reg = 0; i_reg < cur_chr_testing_haplocoded_tag_geno_regs->size(); i_reg++)
		{
			double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant_optimized(cur_chr_testing_haplocoded_tag_geno_regs->at(i_reg), cur_chrom_recomb_regs);			
			cur_chr_testing_haplocoded_tag_geno_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
		} // i_reg loop.

		//for (int i_reg = 0; i_reg < cur_chrom_ref_tag_target_var_regs->size(); i_reg++)
		//{
		//	double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant(cur_chrom_ref_tag_target_var_regs->at(i_reg), cur_chrom_recomb_regs);
		//	cur_chrom_ref_tag_target_var_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
		//} // i_reg loop.

		for (int tag_var_i = 0; tag_var_i < cur_chr_testing_haplocoded_tag_geno_regs->size(); tag_var_i++)
		{
			int end_tag_var_i = tag_var_i;
			while (end_tag_var_i < cur_chr_testing_haplocoded_tag_geno_regs->size() &&
				fabs(cur_chr_testing_haplocoded_tag_geno_regs->at(end_tag_var_i)->dbl_score - cur_chr_testing_haplocoded_tag_geno_regs->at(tag_var_i)->dbl_score) < tag_2_tag_distance_cM)
			{
				end_tag_var_i++;
			}

			// Make sure the end does not pass over the end.
			end_tag_var_i = MIN(end_tag_var_i, cur_chr_testing_haplocoded_tag_geno_regs->size()-1);

			fprintf(f_op, "%s\t%d\t%d\t%d\t%.4f\t%.4f\n", 
				cur_chr_testing_haplocoded_tag_geno_regs->at(tag_var_i)->chrom, cur_chr_testing_haplocoded_tag_geno_regs->at(tag_var_i)->start, 
				tag_var_i, end_tag_var_i, 
				cur_chr_testing_haplocoded_tag_geno_regs->at(tag_var_i)->dbl_score,
				cur_chr_testing_haplocoded_tag_geno_regs->at(end_tag_var_i)->dbl_score);
		} // tag_var_i loop.
	} // i_chr loop.
	fclose(f_op);
	fprintf(stderr, "Done!\n");
}

void run_Imputation_ForeBack_Full_States_Sliding_Windows(char* reference_tag_target_haplocoded_genotypes_fp, char* ref_sample_ids_list_fp,
	char* testing_tag_haplocoded_genotypes_fp, char* testing_sample_ids_list_fp,
	char* known_tag_target_haplocoded_genotypes_fp,
	char* recombination_rate_dir,
	double tag_2_tag_distance_cM,
	double N_e,
	char* target_focus_reg_BED_fp,
	char* geno_probs_op_fp)
{
	double l_target_2_center_pred_buffer_in_cM = 0.1;

	//double N_e = pow(10, 4);
	double allele_err = pow(10, -4);

	fprintf(stderr, "Running fore-back sliding window with tag2tag distance of %.3f cMs with effective pop. size N_e=%d\n", tag_2_tag_distance_cM, (int)N_e);

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
			imputed_allele_prob_per_sample_per_hap[i_s][0][0] = xlog(0.0);
			imputed_allele_prob_per_sample_per_hap[i_s][0][1] = xlog(0.0);

			imputed_allele_prob_per_sample_per_hap[i_s][1] = new double[2];
			imputed_allele_prob_per_sample_per_hap[i_s][1][0] = xlog(0.0);
			imputed_allele_prob_per_sample_per_hap[i_s][1][1] = xlog(0.0);
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

		// Add one element to the beginning of the regions for the current chromosome.
		vector<t_annot_region*>* cur_chr_testing_haplocoded_tag_geno_regs = restr_testing_haplocoded_tag_geno_regs->regions_per_chrom[i_chr];
		vector<t_annot_region*>* cur_chrom_ref_tag_target_var_regs = restr_reference_haplocoded_tag_target_geno_regs->regions_per_chrom[i_chr];

		// Assign the recomb rates.
		fprintf(stderr, "Setting recombination rates using %d recombination mapped regions.\n", cur_chrom_recomb_regs->size());
		for (int i_reg = 0; i_reg < cur_chr_testing_haplocoded_tag_geno_regs->size(); i_reg++)
		{
			double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant(cur_chr_testing_haplocoded_tag_geno_regs->at(i_reg), cur_chrom_recomb_regs);
			cur_chr_testing_haplocoded_tag_geno_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
		} // i_reg loop.

		for (int i_reg = 0; i_reg < cur_chrom_ref_tag_target_var_regs->size(); i_reg++)
		{
			double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant(cur_chrom_ref_tag_target_var_regs->at(i_reg), cur_chrom_recomb_regs);
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

			for(int win_start_var_i = 0; win_start_var_i < cur_chr_testing_haplocoded_tag_geno_regs->size(); win_start_var_i++)
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

				////////////////////////////////////////////////////////////////////////////////////
				// BELOW DOES NOT WORK YET.
				// Copy the previous window if they exist.
				bool PERFORM_FULL_COMPUTE = true;

				if (PERFORM_FULL_COMPUTE ||
					prev_win_var_regs == NULL)
				{
					// Initialize the state probabilities for both haplotypes.
					for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
					{
						fore_scores_per_hap[0][0][state_i] = xlog((double)1.0 / n_ref_haplotypes);
						fore_scores_per_hap[0][1][state_i] = xlog((double)1.0 / n_ref_haplotypes);
					} // state_i loop.

					// Start recursing over the variants.
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
						for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
						{
							// Recurse over the previous states.
							fore_scores_per_hap[var_i][1][state_i] = xlog(0);
							fore_scores_per_hap[var_i][0][state_i] = xlog(0);
							for (int prev_state_i = 0; prev_state_i < n_ref_haplotypes; prev_state_i++)
							{
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
								double emit_prob0 = 0;
								double emit_prob1 = 0;

								if (var_i == n_vars + 1)
								{
									emit_prob0 = 1.0;
									emit_prob1 = 1.0;
								}
								else
								{
									int ref_sample_i = state_i / 2;
									int ref_hap_i = state_i % 2;
									int cur_ref_allele = get_allele_per_haplotype(ref_sample_geno[ref_sample_i], ref_hap_i);
									int cur_test_allele0 = get_allele_per_haplotype(test_sample_geno[test_sample_i], 0);
									int cur_test_allele1 = get_allele_per_haplotype(test_sample_geno[test_sample_i], 1);

									double emit_prob_per_ref_allele[2];
									emit_prob_per_ref_allele[cur_ref_allele] = 1 - allele_err;
									emit_prob_per_ref_allele[1 - cur_ref_allele] = allele_err;
									emit_prob0 = emit_prob_per_ref_allele[cur_test_allele0];
									emit_prob1 = emit_prob_per_ref_allele[cur_test_allele1];
								}
								// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 

								double trans_emit_prob0 = xlog_mul(xlog(trans_prob), xlog(emit_prob0));
								double trans_emit_prob1 = xlog_mul(xlog(trans_prob), xlog(emit_prob1));

								fore_scores_per_hap[var_i][0][state_i] = xlog_sum(fore_scores_per_hap[var_i][0][state_i], xlog_mul(trans_emit_prob0, fore_scores_per_hap[var_i - 1][0][prev_state_i]));
								fore_scores_per_hap[var_i][1][state_i] = xlog_sum(fore_scores_per_hap[var_i][1][state_i], xlog_mul(trans_emit_prob1, fore_scores_per_hap[var_i - 1][1][prev_state_i]));

								if (__DUMP_SHMMR_GIMP_MESSAGES__)
								{
									//fprintf(stderr, "fore[%d][%d]: %.5f: P_emit(%d, %d)=%.5f; P_trans(%d->%d): %.5f\n", var_i, state_i, fore_scores[var_i][state_i], cur_ref_allele, cur_test_allele, emit_prob, prev_state_i, state_i, trans_prob);
								}
							} // prev_state_i loop.
						} // state_i loop.
					} // var_i loop.

					// Move the backward scores.
					fprintf(stderr, "Computing back-ward scores for sample %d\n", test_sample_i);			
				
					for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
					{
						back_scores_per_hap[n_vars + 1][0][state_i] = xlog((double)1.0 / n_ref_haplotypes);
						back_scores_per_hap[n_vars + 1][1][state_i] = xlog((double)1.0 / n_ref_haplotypes);
					} // state_i loop.

					// Compute backward probabilities: Our backward is a little different to be able to perform below marginalization.
					for (int var_i = n_vars; var_i >= 0; var_i--)
					{
						if (var_i % 10 == 0)
						{
							fprintf(stderr, "Backward: sample_i: %d/%d: var_i: %d         \r", test_sample_i, testing_sample_ids->size(), var_i);
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
						if(var_i > 0)
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

						// Go over all the reference haplotypes.
						for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
						{
							back_scores_per_hap[var_i][0][state_i] = xlog(0);
							back_scores_per_hap[var_i][1][state_i] = xlog(0);
							for (int next_state_i = 0; next_state_i < n_ref_haplotypes; next_state_i++)
							{
								// Same state ---, different state ---.
								double trans_prob = 0;
								if (state_i == next_state_i)
								{
									trans_prob = self_prob;
								}
								else
								{
									trans_prob = other_prob;
								}

								// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
								// Set the emission prob: If the ref allele is matching to test, set to 0.99999, otherwise set to the error prob.
								// This is the emission of the current allele.
								double emit_prob0 = 0; // 0.999 for error, 
								double emit_prob1 = 0; // 0.999 for error, 
								if (var_i == 0)
								{
									emit_prob0 = 1.0;
									emit_prob1 = 1.0;
								}
								else
								{
									int ref_sample_i = state_i / 2;
									int ref_hap_i = state_i % 2;
									int cur_ref_allele = get_allele_per_haplotype(ref_sample_geno[ref_sample_i], ref_hap_i);
									//int test_allele = get_allele_per_haplotype(test_sample_geno[test_sample_i], test_hap_i);

									int cur_test_allele0 = get_allele_per_haplotype(test_sample_geno[test_sample_i], 0);
									int cur_test_allele1 = get_allele_per_haplotype(test_sample_geno[test_sample_i], 1);

									double emit_prob_per_ref_allele[2];
									emit_prob_per_ref_allele[cur_ref_allele] = 1 - allele_err;
									emit_prob_per_ref_allele[1 - cur_ref_allele] = allele_err;
									emit_prob0 = emit_prob_per_ref_allele[cur_test_allele0];
									emit_prob1 = emit_prob_per_ref_allele[cur_test_allele1];
								}
								// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 

								double trans_emit_prob0 = xlog_mul(xlog(trans_prob), xlog(emit_prob0));
								double trans_emit_prob1 = xlog_mul(xlog(trans_prob), xlog(emit_prob1));

								// Update the backward score.
								back_scores_per_hap[var_i][0][state_i] = xlog_sum(back_scores_per_hap[var_i][0][state_i], xlog_mul(trans_emit_prob0, back_scores_per_hap[var_i + 1][0][next_state_i]));
								back_scores_per_hap[var_i][1][state_i] = xlog_sum(back_scores_per_hap[var_i][1][state_i], xlog_mul(trans_emit_prob1, back_scores_per_hap[var_i + 1][1][next_state_i]));
							} // prev_state_i loop.
						} // state_i loop.
					} // var_i loop.
				} //  full compute check.

				// Compute the total log forward and backward probabilities.
				double per_hap_total_log_fore_prob[2];
				per_hap_total_log_fore_prob[0] = xlog(0.0);
				per_hap_total_log_fore_prob[1] = xlog(0.0);
				double per_hap_total_log_back_prob[2];
				per_hap_total_log_back_prob[0] = xlog(0.0);
				per_hap_total_log_back_prob[1] = xlog(0.0);

				for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
				{
					per_hap_total_log_fore_prob[0] = xlog_sum(per_hap_total_log_fore_prob[0], fore_scores_per_hap[n_vars + 1][0][state_i]);
					per_hap_total_log_back_prob[0] = xlog_sum(per_hap_total_log_back_prob[0], back_scores_per_hap[0][0][state_i]);

					per_hap_total_log_fore_prob[1] = xlog_sum(per_hap_total_log_fore_prob[1], fore_scores_per_hap[n_vars + 1][1][state_i]);
					per_hap_total_log_back_prob[1] = xlog_sum(per_hap_total_log_back_prob[1], back_scores_per_hap[0][1][state_i]);
				} // state_i loop.

				fprintf(stderr, "Haplotype 0 total probabilities: fore=%.5f ;; back=%.5f\n", per_hap_total_log_fore_prob[0], per_hap_total_log_back_prob[0]);
				fprintf(stderr, "Haplotype 1 total probabilities: fore=%.5f ;; back=%.5f\n", per_hap_total_log_fore_prob[1], per_hap_total_log_back_prob[1]);

				// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
				// Start computing the posterior probabilities for the untyped markers.
				// Now go over the untyped genotypes and compute the posterior probability of all the untyped markers.				
				int n_untyped = 0;
				double min_tar2tag_cM = tag_2_tag_distance_cM / 2 - l_target_2_center_pred_buffer_in_cM;
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
						if (__DUMP_SHMMR_GIMP_MESSAGES__)
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
					if (imputed_allele_prob_per_sample_per_hap[test_sample_i][0][0] != xlog(0.0) ||
						imputed_allele_prob_per_sample_per_hap[test_sample_i][0][1] != xlog(0.0))
					{
						continue;
					}

					n_untyped++;

					fprintf(stderr, "@%d. untyped marker         \r", n_untyped);

					// Find the closest tag variants to this target: This is the fore-back coordinates.
					for (int var_i = 1; var_i < n_vars; var_i++)
					{
						// Find the untyped target variant in between the current tags.
						if (cur_win_var_regs->at(var_i)->start < cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->start &&
							cur_win_var_regs->at(var_i + 1)->start > cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->start)
						{
							fprintf(stderr, "Found %d-%d-%d\n",
								cur_win_var_regs->at(var_i)->start,
								cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->start,
								cur_win_var_regs->at(var_i + 1)->start);

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
							double r_m = fabs(cur_cM - next_cM);
							double rho_m = 4 * N_e * r_m;

							double tau_m = 1 - exp(-1 * rho_m / n_ref_haplotypes);

							double other_prob = tau_m / n_ref_haplotypes;
							double self_prob = (1 - tau_m) + (tau_m / n_ref_haplotypes);
							//////////////////////////////////////////////////////////////////////////////////////////////////

							for (int hap_i = 0; hap_i < 2; hap_i++)
							{
								double per_allele_probs[2];
								per_allele_probs[0] = xlog(0);
								per_allele_probs[1] = xlog(0);
								for (int untyped_state_i = 0; untyped_state_i < n_ref_haplotypes; untyped_state_i++)
								{
									//for (int untyped_state_j = 0; untyped_state_j < n_ref_haplotypes; untyped_state_j++)
									int untyped_state_j = untyped_state_i;
									{
										// For this untyped state in the refernce, we would like to sum all the paths that pass through here.
										int ref_sample_i = untyped_state_i / 2;
										int ref_hap_i = untyped_state_i % 2;
										int cur_untyped_ref_allele = get_allele_per_haplotype(untyped_var_ref_geno[ref_sample_i], ref_hap_i);

										double cur_untyped_state_fore_prob = xlog_mul(other_prob, fore_scores_per_hap[var_i][hap_i][untyped_state_i]);;
										if (untyped_state_i == untyped_state_j)
										{
											cur_untyped_state_fore_prob = xlog_mul(self_prob, fore_scores_per_hap[var_i][hap_i][untyped_state_i]);
										}

										double cur_untyped_state_back_prob = back_scores_per_hap[var_i + 1][hap_i][untyped_state_j];

										per_allele_probs[cur_untyped_ref_allele] = xlog_sum(per_allele_probs[cur_untyped_ref_allele], xlog_mul(cur_untyped_state_fore_prob, cur_untyped_state_back_prob));
									} // untyped_state_j loop/block.
								} // untyped_state_i loop.

								// We went through all the states, normalize the probabilities.
								per_allele_probs[0] = xlog_div(per_allele_probs[0], per_hap_total_log_fore_prob[hap_i]);
								per_allele_probs[1] = xlog_div(per_allele_probs[1], per_hap_total_log_fore_prob[hap_i]);

								int known_geno = -1;
								if (known_var_ref_geno)
								{
									known_geno = (int)(known_var_ref_geno[test_sample_i]);
								}

								fprintf(stderr, "Var %d on sample %d: %s:%d (%s): Hap%d: [0:%.4f, 1:%.4f]; (known_geno=%d)\n",
									i_ref_var,
									test_sample_i,
									cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->chrom,
									cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->start,
									cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->name,
									hap_i,
									per_allele_probs[0], per_allele_probs[1], known_geno);

								imputed_allele_prob_per_sample_per_hap[test_sample_i][hap_i][0] = per_allele_probs[0];
								imputed_allele_prob_per_sample_per_hap[test_sample_i][hap_i][1] = per_allele_probs[1];
							} // hap_i loop.

							// Set the inputated variant type.
							cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->score = VAR_UNTYPED_TARGET_IMPUTED;

							break; // Breaks out of position search.
						} // middle check.
					} // i_var loop.
				} // i_ref_var loop.

				// Copy the current windows.
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

double get_self(double a)
{
	return(a);
}

void run_Imputation_Viterbi_Full_States_Sliding_Windows(char* reference_tag_target_haplocoded_genotypes_fp, char* ref_sample_ids_list_fp,
	char* testing_tag_haplocoded_genotypes_fp, char* testing_sample_ids_list_fp,
	char* known_tag_target_haplocoded_genotypes_fp,
	char* recombination_rate_dir,
	double tag_2_tag_distance_cM,
	double N_e,
	char* target_focus_reg_BED_fp,
	bool math_mode,
	double global_scaler_in_log,
	char* geno_probs_op_fp)
{
	double l_target_2_center_pred_buffer_in_cM = 0.1;

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

		// Add one element to the beginning of the regions for the current chromosome.
		vector<t_annot_region*>* cur_chr_testing_haplocoded_tag_geno_regs = restr_testing_haplocoded_tag_geno_regs->regions_per_chrom[i_chr];
		vector<t_annot_region*>* cur_chrom_ref_tag_target_var_regs = restr_reference_haplocoded_tag_target_geno_regs->regions_per_chrom[i_chr];

		// Assign the recomb rates.
		fprintf(stderr, "Setting recombination rates using %d recombination mapped regions.\n", cur_chrom_recomb_regs->size());
		for (int i_reg = 0; i_reg < cur_chr_testing_haplocoded_tag_geno_regs->size(); i_reg++)
		{
			double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant(cur_chr_testing_haplocoded_tag_geno_regs->at(i_reg), cur_chrom_recomb_regs);
			cur_chr_testing_haplocoded_tag_geno_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
		} // i_reg loop.

		for (int i_reg = 0; i_reg < cur_chrom_ref_tag_target_var_regs->size(); i_reg++)
		{
			double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant(cur_chrom_ref_tag_target_var_regs->at(i_reg), cur_chrom_recomb_regs);
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

					// Start recursing over the variants.
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
						for (int state_i = 0; state_i < n_ref_haplotypes; state_i++)
						{
							// Recurse over the previous states.
							ML_scores_per_hap[var_i][1][state_i] = ZERO;
							ML_scores_per_hap[var_i][0][state_i] = ZERO;
							for (int prev_state_i = 0; prev_state_i < n_ref_haplotypes; prev_state_i++)
							{
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
								double emit_prob0 = 0;
								double emit_prob1 = 0;

								if (var_i == n_vars + 1)
								{
									emit_prob0 = 1.0;
									emit_prob1 = 1.0;
								}
								else
								{
									// There is no conditional below.
									int ref_sample_i = state_i / 2;
									int ref_hap_i = state_i % 2;
									int cur_ref_allele = get_allele_per_haplotype(ref_sample_geno[ref_sample_i], ref_hap_i);
									int cur_test_allele0 = get_allele_per_haplotype(test_sample_geno[test_sample_i], 0);
									int cur_test_allele1 = get_allele_per_haplotype(test_sample_geno[test_sample_i], 1);

									double emit_prob_per_ref_allele[2];
									emit_prob_per_ref_allele[cur_ref_allele] = 1 - allele_err;
									emit_prob_per_ref_allele[1 - cur_ref_allele] = allele_err;
									emit_prob0 = emit_prob_per_ref_allele[cur_test_allele0];
									emit_prob1 = emit_prob_per_ref_allele[cur_test_allele1];
								}
								// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 

								double trans_emit_prob0 = XMUL(XCONVERT_LIN(trans_prob), XCONVERT_LIN(emit_prob0));
								double trans_emit_prob1 = XMUL(XCONVERT_LIN(trans_prob), XCONVERT_LIN(emit_prob1));

								// Add the scaler for this position.
								trans_emit_prob0 = XMUL(trans_emit_prob0, XCONVERT_LIN(global_lin_scaler));
								trans_emit_prob1 = XMUL(trans_emit_prob1, XCONVERT_LIN(global_lin_scaler));
								
								ML_scores_per_hap[var_i][0][state_i] = MAX(ML_scores_per_hap[var_i][0][state_i], XMUL(trans_emit_prob0, ML_scores_per_hap[var_i - 1][0][prev_state_i]));
								ML_scores_per_hap[var_i][1][state_i] = MAX(ML_scores_per_hap[var_i][1][state_i], XMUL(trans_emit_prob1, ML_scores_per_hap[var_i - 1][1][prev_state_i]));


								if (XCOMPARE(ML_scores_per_hap[var_i][0][state_i], XMUL(trans_emit_prob0, ML_scores_per_hap[var_i - 1][0][prev_state_i])))
								{
									ML_prev_state_per_hap[var_i][0][state_i] = prev_state_i;
								}

								if (XCOMPARE(ML_scores_per_hap[var_i][1][state_i], XMUL(trans_emit_prob1, ML_scores_per_hap[var_i - 1][1][prev_state_i])))
								{
									ML_prev_state_per_hap[var_i][1][state_i] = prev_state_i;
								}

								if (__DUMP_SHMMR_GIMP_MESSAGES__)
								{
									fprintf(stderr, "ML[%d][%d]: %.5f;;%.5f: P_emit=%.5f, %.5f; P_trans(%d->%d): %.5f\n", 
											var_i, state_i, ML_scores_per_hap[var_i][0][state_i], ML_scores_per_hap[var_i][1][state_i],
											emit_prob0, emit_prob1, prev_state_i, state_i, trans_prob);
								}
							} // prev_state_i loop.

							// Save.
							if (var_i <= n_vars &&
								f_viterbi_arrays != NULL)
							{
								fprintf(f_viterbi_arrays, "%d\t%d\t0\t%d\t%d\t%.5f\n", win_start_var_i, test_sample_i, var_i, state_i, XCONVERT_2_LOG(ML_scores_per_hap[var_i][0][state_i]));
								fprintf(f_viterbi_arrays, "%d\t%d\t1\t%d\t%d\t%.5f\n", win_start_var_i, test_sample_i, var_i, state_i, XCONVERT_2_LOG(ML_scores_per_hap[var_i][1][state_i]));
							}

							if (math_mode == MATH_MODE_LIN_SPACE)
							{
								if (log(ML_scores_per_hap[var_i][0][state_i]) / log(10) > 50 ||
									log(ML_scores_per_hap[var_i][0][state_i]) / log(10) < -50 ||
									log(ML_scores_per_hap[var_i][1][state_i]) / log(10) > 50 ||
									log(ML_scores_per_hap[var_i][1][state_i]) / log(10) < -50)
								{
									fprintf(stderr, "Out of range @ %d: %.5f\n", var_i, ML_scores_per_hap[var_i][0][state_i],
										ML_scores_per_hap[var_i][1][state_i]);
									exit(0);
								}
							}
						} // state_i loop.
					} // var_i loop.
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
				double min_tar2tag_cM = MAX(0, tag_2_tag_distance_cM / 2 - l_target_2_center_pred_buffer_in_cM);
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
						if (__DUMP_SHMMR_GIMP_MESSAGES__)
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

					fprintf(stderr, "@%d. untyped marker         \r", n_untyped);

					// Find the closest tag variants to this target:
					for (int var_i = 1; var_i < n_vars; var_i++)
					{
						// Find the untyped target variant in between the current tags.
						if (cur_win_var_regs->at(var_i)->start < cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->start &&
							cur_win_var_regs->at(var_i + 1)->start > cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->start)
						{
							fprintf(stderr, "Found %d-%d-%d\n",
								cur_win_var_regs->at(var_i)->start,
								cur_chrom_ref_tag_target_var_regs->at(i_ref_var)->start,
								cur_win_var_regs->at(var_i + 1)->start);

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

				  // Copy the current windows.
				if (prev_win_var_regs != NULL)
				{
					delete prev_win_var_regs;
				}

				prev_win_var_regs = cur_win_var_regs;
			} // win_start_i loop.
		} // test_sample_i loop.
	} // i_chr loop.

	close_f(f_viterbi_arrays, arrays_op_fp);

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

void get_PR_stats_per_3entry_genotype_probs(char* imputed_genotypes_fp, char* imputed_sample_ids_list_fp, char* known_genotypes_fp, char* known_sample_ids_list_fp)
{
	vector<t_annot_region*>* known_genotype_regs = load_variant_signal_regions_wrapper(known_genotypes_fp, known_sample_ids_list_fp);
	vector<char*>* known_sample_ids = buffer_file(known_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d known genotype regions for %d samples.\n", known_genotype_regs->size(), known_sample_ids->size());

	vector<char*>* imputed_sample_ids = buffer_file(imputed_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d imputed sample ids\n", imputed_sample_ids->size());

	vector<t_annot_region*>* imputed_genotype_regs = load_BED_with_line_information(imputed_genotypes_fp);
	fprintf(stderr, "Loaded %d imputed variants.\n", imputed_genotype_regs->size());
	for (int i_reg = 0; i_reg < imputed_genotype_regs->size(); i_reg++)
	{
		imputed_genotype_regs->at(i_reg)->score = 0;
	} // i_reg loop.

	if (imputed_genotype_regs->size() == 0)
	{
		fprintf(stderr, "Found no imputed variants.\n");
		return;
	}

	// Set the mapping between sample id's.
	vector<int>* imp_2_known_sample_i = new vector<int>();
	int n_matched_samples = 0;
	for (int imp_i = 0; imp_i < imputed_sample_ids->size(); imp_i++)
	{
		int cur_imp_known_i = t_string::get_i_str(known_sample_ids, imputed_sample_ids->at(imp_i));
		if (cur_imp_known_i < known_sample_ids->size())
		{
			n_matched_samples++;
		}
		imp_2_known_sample_i->push_back(cur_imp_known_i);
	} // imp_i loop.
	fprintf(stderr, "Matched %d samples.\n", n_matched_samples);

	fprintf(stderr, "Intersecting imputed and known variants.\n");

	// Intersect and process.
	double** known_imp_sample_geno = new double*[2];
	known_imp_sample_geno[0] = new double[n_matched_samples + 2];
	known_imp_sample_geno[1] = new double[n_matched_samples + 2];
	vector<t_annot_region*>* intersects = intersect_annot_regions(imputed_genotype_regs, known_genotype_regs, true);

	enum {
		N_ALL_CORRECT, N_ALL_IMPUTED, N_ALL_TOTAL,
		N_NR_CORRECT, N_NR_IMPUTED, N_NR_TOTAL, N_PR_ENTRIES
	};

	int n_bins = 100;
	double MIN_SCORE = 0.0;
	double MAX_SCORE = 1.0;
	double** all_PR_info = new double*[n_bins];
	for (int i_bin = 0; i_bin < n_bins; i_bin++)
	{
		all_PR_info[i_bin] = new double[N_PR_ENTRIES];
		all_PR_info[i_bin][N_ALL_CORRECT] = 0; // This is the correct
		all_PR_info[i_bin][N_ALL_IMPUTED] = 0; // This is the total predicted @ this probability threshold.
		all_PR_info[i_bin][N_ALL_TOTAL] = 0; // This is the total known

		all_PR_info[i_bin][N_NR_CORRECT] = 0; // This is the total known
		all_PR_info[i_bin][N_NR_IMPUTED] = 0; // This is the total imputed
		all_PR_info[i_bin][N_NR_TOTAL] = 0; // This is the total known
	} // i_bin loop.
	fprintf(stderr, "Setup %d PR-bins\n", n_bins);

	double l_bin = (MAX_SCORE - MIN_SCORE) / n_bins;

	fprintf(stderr, "Found %d intersections\n", intersects->size());
	int n_processed_imputed_targets = 0;
	for (int i_int = 0; i_int < intersects->size(); i_int++)
	{
		if (i_int % 1000 == 0)
		{
			fprintf(stderr, "@ %d. intersect         \r", i_int);
		}

		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* imp_reg = int_info->src_reg;
		t_annot_region* known_reg = int_info->dest_reg;

		void** known_reg_info = (void**)(known_reg->data);
		char* known_reg_geno = (char*)(known_reg_info[0]);

		if (t_string::compare_strings(imp_reg->name, known_reg->name) &&
			imp_reg->score == 0)
		{
			imp_reg->score = 1;
			n_processed_imputed_targets++;

			double total_n_non_refs = 0;
			double n_matching_non_refs = 0;
			double n_matching_all = 0;
			double n_all = 0;

			int geno_i = 0;
			char* imp_var_line = (char*)(imp_reg->data);
			t_string_tokens* toks = t_string::tokenize_by_chars(imp_var_line, "\t");
			for (int imp_i = 0; imp_i < imputed_sample_ids->size(); imp_i++)
			{
				int tok_i = imp_i + 4;
				t_string_tokens* prob_toks = toks->at(tok_i)->tokenize_by_chars(";");
				if (prob_toks->size() != 3)
				{
					fprintf(stderr, "Could not parse the probability tokens @ %s on %s\n", toks->at(tok_i)->str(), imp_var_line);
					exit(0);
				}

				double imp_geno_probs[3];
				imp_geno_probs[0] = atof(prob_toks->at(0)->str());
				imp_geno_probs[1] = atof(prob_toks->at(1)->str());
				imp_geno_probs[2] = atof(prob_toks->at(2)->str());
				t_string::clean_tokens(prob_toks);
				
				// Get the MAP genotype.
				double map_geno_prob = -1;
				int map_geno = -1;
				for (int geno = 0; geno < 3; geno++)
				{
					if (map_geno_prob < imp_geno_probs[geno])
					{
						map_geno_prob = imp_geno_probs[geno];
						map_geno = geno;
					}
				} // geno loop.

				int cur_prob_i_bin = (int)((map_geno_prob - MIN_SCORE) / l_bin);
				if (cur_prob_i_bin >= n_bins)
				{
					cur_prob_i_bin = n_bins - 1;
				}

				if (__DUMP_SHMMR_GIMP_MESSAGES__)
				{
					fprintf(stderr, "%s: Geno probs: %.4f, %.4f, %.4f;; MAP_geno: %d (%.4f);; i_bin: %d\n",
						toks->at(tok_i)->str(),
						imp_geno_probs[0], imp_geno_probs[1], imp_geno_probs[2],
						map_geno, map_geno_prob,
						cur_prob_i_bin);
				}

				if (imp_2_known_sample_i->at(imp_i) < known_sample_ids->size())
				{
					int known_geno = (int)(known_reg_geno[imp_2_known_sample_i->at(imp_i)]);

					// Update non-ref PR stats.
					if (known_geno > 0)
					{
						if (known_geno == map_geno)
						{
							for (int i_bin = 0; i_bin <= cur_prob_i_bin; i_bin++)
							{
								all_PR_info[i_bin][N_NR_CORRECT]++;
							} // i_bin loop.
						}

						for (int i_bin = 0; i_bin <= cur_prob_i_bin; i_bin++)
						{
							all_PR_info[i_bin][N_NR_IMPUTED]++;
						} // i_bin loop.

						for (int i_bin = 0; i_bin < n_bins; i_bin++)
						{
							all_PR_info[i_bin][N_NR_TOTAL]++;
						} // i_bin loop.
					} // non-ref check.

					  // Update all PR stats.
					if (known_geno == map_geno)
					{
						for (int i_bin = 0; i_bin <= cur_prob_i_bin; i_bin++)
						{
							all_PR_info[i_bin][N_ALL_CORRECT]++;
						} // i_bin loop.
					}

					for (int i_bin = 0; i_bin <= cur_prob_i_bin; i_bin++)
					{
						all_PR_info[i_bin][N_ALL_IMPUTED]++;
					} // i_bin loop.

					for (int i_bin = 0; i_bin < n_bins; i_bin++)
					{
						all_PR_info[i_bin][N_ALL_TOTAL]++;
					} // i_bin loop.					
				}
			} // imp_i loop.

			t_string::clean_tokens(toks);
		} // overlapping region name comparison.
	} // i_int loop.

	  // Write the PR statistics.
	FILE* f_op = open_f("PR_stats.txt", "w");
	fprintf(f_op, "SCORE_THRESHOLD\tALL_CORRECT\tALL_IMPUTED\tALL_TOTAL\tNR_CORRECT\tNR_IMPUTED\tNR_TOTAL\n");
	for (int i_bin = 0; i_bin < n_bins; i_bin++)
	{
		fprintf(f_op, "%.3f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
			MIN_SCORE + i_bin * l_bin,
			all_PR_info[i_bin][N_ALL_CORRECT], all_PR_info[i_bin][N_ALL_IMPUTED], all_PR_info[i_bin][N_ALL_TOTAL],
			all_PR_info[i_bin][N_NR_CORRECT], all_PR_info[i_bin][N_NR_IMPUTED], all_PR_info[i_bin][N_NR_TOTAL]);
	} // i_bin loop.
	fclose(f_op);
}

void get_PR_stats_per_GIMP_4entry_allelic_probs(char* imputed_genotypes_fp, char* imputed_sample_ids_list_fp, char* known_genotypes_fp, char* known_sample_ids_list_fp,
	bool flag_imputed_probs_are_linear)
{
	fprintf(stderr, "Computing PR stats for linear=%d probs.\n", flag_imputed_probs_are_linear);

	vector<t_annot_region*>* known_genotype_regs = load_variant_signal_regions_wrapper(known_genotypes_fp, known_sample_ids_list_fp);
	vector<char*>* known_sample_ids = buffer_file(known_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d known genotype regions for %d samples.\n", known_genotype_regs->size(), known_sample_ids->size());

	vector<char*>* imputed_sample_ids = buffer_file(imputed_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d imputed sample ids\n", imputed_sample_ids->size());

	vector<t_annot_region*>* imputed_genotype_regs = load_BED_with_line_information(imputed_genotypes_fp);
	fprintf(stderr, "Loaded %d imputed variants.\n", imputed_genotype_regs->size());
	for (int i_reg = 0; i_reg < imputed_genotype_regs->size(); i_reg++)
	{
		imputed_genotype_regs->at(i_reg)->score = 0;
	} // i_reg loop.

	if (imputed_genotype_regs->size() == 0)
	{
		fprintf(stderr, "Found no imputed variants.\n");
		return;
	}

	// Set the mapping between sample id's.
	vector<int>* imp_2_known_sample_i = new vector<int>();
	int n_matched_samples = 0;
	for (int imp_i = 0; imp_i < imputed_sample_ids->size(); imp_i++)
	{
		int cur_imp_known_i = t_string::get_i_str(known_sample_ids, imputed_sample_ids->at(imp_i));
		if (cur_imp_known_i < known_sample_ids->size())
		{
			n_matched_samples++;
		}
		imp_2_known_sample_i->push_back(cur_imp_known_i);
	} // imp_i loop.
	fprintf(stderr, "Matched %d samples.\n", n_matched_samples);

	fprintf(stderr, "Intersecting imputed and known variants.\n");

	// Intersect and process.
	double** known_imp_sample_geno = new double*[2];
	known_imp_sample_geno[0] = new double[n_matched_samples + 2];
	known_imp_sample_geno[1] = new double[n_matched_samples + 2];
	vector<t_annot_region*>* intersects = intersect_annot_regions(imputed_genotype_regs, known_genotype_regs, true);

	enum { N_ALL_CORRECT , N_ALL_IMPUTED, N_ALL_TOTAL,
			N_NR_CORRECT, N_NR_IMPUTED, N_NR_TOTAL, N_PR_ENTRIES};

	int n_bins = 100;	
	double MIN_SCORE = 0.0;
	double MAX_SCORE = 1.0;
	double** all_PR_info = new double*[n_bins];
	for (int i_bin = 0; i_bin < n_bins; i_bin++)
	{
		all_PR_info[i_bin] = new double[N_PR_ENTRIES];
		all_PR_info[i_bin][N_ALL_CORRECT] = 0; // This is the correct
		all_PR_info[i_bin][N_ALL_IMPUTED] = 0; // This is the total predicted @ this probability threshold.
		all_PR_info[i_bin][N_ALL_TOTAL] = 0; // This is the total known

		all_PR_info[i_bin][N_NR_CORRECT] = 0; // This is the total known
		all_PR_info[i_bin][N_NR_IMPUTED] = 0; // This is the total imputed
		all_PR_info[i_bin][N_NR_TOTAL] = 0; // This is the total known
	} // i_bin loop.
	fprintf(stderr, "Setup %d PR-bins\n", n_bins);

	double l_bin = (MAX_SCORE - MIN_SCORE) / n_bins;

	fprintf(stderr, "Found %d intersections\n", intersects->size());
	int n_processed_imputed_targets = 0;
	for (int i_int = 0; i_int < intersects->size(); i_int++)
	{
		if (i_int % 1000 == 0)
		{
			fprintf(stderr, "@ %d. intersect         \r", i_int);
		}

		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* imp_reg = int_info->src_reg;
		t_annot_region* known_reg = int_info->dest_reg;

		void** known_reg_info = (void**)(known_reg->data);
		char* known_reg_geno = (char*)(known_reg_info[0]);

		if (t_string::compare_strings(imp_reg->name, known_reg->name) &&
			imp_reg->score == 0)
		{
			imp_reg->score = 1;
			n_processed_imputed_targets++;

			double total_n_non_refs = 0;
			double n_matching_non_refs = 0;
			double n_matching_all = 0;
			double n_all = 0;

			int geno_i = 0;
			char* imp_var_line = (char*)(imp_reg->data);
			t_string_tokens* toks = t_string::tokenize_by_chars(imp_var_line, "\t");
			for (int imp_i = 0; imp_i < imputed_sample_ids->size(); imp_i++)
			{
				int tok_i = imp_i + 4;
				t_string_tokens* prob_toks = toks->at(tok_i)->tokenize_by_chars(";");
				if (prob_toks->size() != 4)
				{
					fprintf(stderr, "Could not parse the probability tokens @ %s on %s\n", toks->at(tok_i)->str(), imp_var_line);
					exit(0);
				}

				double hap_00_score = atof(prob_toks->at(0)->str());
				double hap_01_score = atof(prob_toks->at(1)->str());
				double hap_10_score = atof(prob_toks->at(2)->str());
				double hap_11_score = atof(prob_toks->at(3)->str());
				t_string::clean_tokens(prob_toks);

				double hap0_AA_prob = exp(hap_01_score) / (exp(hap_00_score) + exp(hap_01_score));
				double hap1_AA_prob = exp(hap_11_score) / (exp(hap_10_score) + exp(hap_11_score));
				if (flag_imputed_probs_are_linear)
				{
					hap0_AA_prob = hap_01_score / (hap_00_score + hap_01_score);
					hap1_AA_prob = hap_11_score / (hap_10_score + hap_11_score);
				}

				double imp_geno_probs[3];
				imp_geno_probs[0] = (1 - hap0_AA_prob) * (1 - hap1_AA_prob);
				imp_geno_probs[1] = hap0_AA_prob * (1 - hap1_AA_prob) + hap1_AA_prob * (1 - hap0_AA_prob);
				imp_geno_probs[2] = hap0_AA_prob * hap1_AA_prob;

				// Get the MAP genotype.
				double map_geno_prob = -1;
				int map_geno = -1;
				for (int geno = 0; geno < 3; geno++)
				{
					if (map_geno_prob < imp_geno_probs[geno])
					{
						map_geno_prob = imp_geno_probs[geno];
						map_geno = geno;
					}
				} // geno loop.

				int cur_prob_i_bin = (int)((map_geno_prob - MIN_SCORE) / l_bin);
				if (cur_prob_i_bin >= n_bins)
				{
					cur_prob_i_bin = n_bins - 1;
				}
 				
				if (__DUMP_SHMMR_GIMP_MESSAGES__)
				{
					fprintf(stderr, "%s: hap0_AA: %lf / hap1_AA: %lf; Geno probs: %.4f, %.4f, %.4f;; MAP_geno: %d (%.4f);; i_bin: %d\n",
						toks->at(tok_i)->str(),
						hap0_AA_prob, hap1_AA_prob,
						imp_geno_probs[0], imp_geno_probs[1], imp_geno_probs[2],
						map_geno, map_geno_prob,
						cur_prob_i_bin);
				}

				if (imp_2_known_sample_i->at(imp_i) < known_sample_ids->size())
				{
					int known_geno = (int)(known_reg_geno[imp_2_known_sample_i->at(imp_i)]);

					// Update non-ref PR stats.
					if (known_geno > 0)
					{
						if (known_geno == map_geno)
						{
							for (int i_bin = 0; i_bin <= cur_prob_i_bin; i_bin++)
							{
								all_PR_info[i_bin][N_NR_CORRECT]++;
							} // i_bin loop.
						}

						for (int i_bin = 0; i_bin <= cur_prob_i_bin; i_bin++)
						{
							all_PR_info[i_bin][N_NR_IMPUTED]++;
						} // i_bin loop.

						for (int i_bin = 0; i_bin < n_bins; i_bin++)
						{
							all_PR_info[i_bin][N_NR_TOTAL]++;
						} // i_bin loop.
					} // non-ref check.

					// Update all PR stats.
					if (known_geno == map_geno)
					{
						for (int i_bin = 0; i_bin <= cur_prob_i_bin; i_bin++)
						{
							all_PR_info[i_bin][N_ALL_CORRECT]++;
						} // i_bin loop.
					}
					
					for (int i_bin = 0; i_bin <= cur_prob_i_bin; i_bin++)
					{
						all_PR_info[i_bin][N_ALL_IMPUTED]++;
					} // i_bin loop.

					for (int i_bin = 0; i_bin < n_bins; i_bin++)
					{
						all_PR_info[i_bin][N_ALL_TOTAL]++;
					} // i_bin loop.					
				}
			} // imp_i loop.

			t_string::clean_tokens(toks);
		} // overlapping region name comparison.
	} // i_int loop.
	
	// Write the PR statistics.
	FILE* f_op = open_f("PR_stats.txt", "w");
	fprintf(f_op, "SCORE_THRESHOLD\tALL_CORRECT\tALL_IMPUTED\tALL_TOTAL\tNR_CORRECT\tNR_IMPUTED\tNR_TOTAL\n");
	for (int i_bin = 0; i_bin < n_bins; i_bin++)
	{
		fprintf(f_op, "%.3f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
				MIN_SCORE + i_bin * l_bin,
				all_PR_info[i_bin][N_ALL_CORRECT], all_PR_info[i_bin][N_ALL_IMPUTED], all_PR_info[i_bin][N_ALL_TOTAL],
				all_PR_info[i_bin][N_NR_CORRECT], all_PR_info[i_bin][N_NR_IMPUTED], all_PR_info[i_bin][N_NR_TOTAL]);
	} // i_bin loop.
	fclose(f_op);
}

void get_R2_per_imputed_genotypes(char* imputed_genotypes_fp, char* imputed_sample_ids_list_fp,
	char* known_genotypes_fp, char* known_sample_ids_list_fp)
{
	vector<t_annot_region*>* known_genotype_regs = load_variant_signal_regions_wrapper(known_genotypes_fp, known_sample_ids_list_fp);
	vector<char*>* known_sample_ids = buffer_file(known_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d known genotype regions for %d samples.\n", known_genotype_regs->size(), known_sample_ids->size());

	vector<char*>* imputed_sample_ids = buffer_file(imputed_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d imputed sample ids\n", imputed_sample_ids->size());

	vector<t_annot_region*>* imputed_genotype_regs = load_variant_signal_regions_wrapper(imputed_genotypes_fp, imputed_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d imputed variants.\n", imputed_genotype_regs->size());
	for (int i_reg = 0; i_reg < imputed_genotype_regs->size(); i_reg++)
	{
		imputed_genotype_regs->at(i_reg)->score = 0;
	} // i_reg loop.

	  // Set the mapping between sample id's.
	vector<int>* imp_2_known_sample_i = new vector<int>();
	int n_matched_samples = 0;
	for (int imp_i = 0; imp_i < imputed_sample_ids->size(); imp_i++)
	{
		int cur_imp_known_i = t_string::get_i_str(known_sample_ids, imputed_sample_ids->at(imp_i));
		if (cur_imp_known_i < known_sample_ids->size())
		{
			n_matched_samples++;
		}
		imp_2_known_sample_i->push_back(cur_imp_known_i);
	} // imp_i loop.
	fprintf(stderr, "Matched %d samples.\n", n_matched_samples);

	fprintf(stderr, "Intersecting known regions with imputed regions.\n");
	FILE* f_op = open_f("R2_stats.txt", "w");
	// Intersect and process.
	double** known_imp_sample_geno = new double*[2];
	known_imp_sample_geno[0] = new double[n_matched_samples + 2];
	known_imp_sample_geno[1] = new double[n_matched_samples + 2];
	vector<t_annot_region*>* intersects = intersect_annot_regions(imputed_genotype_regs, known_genotype_regs, true);
	fprintf(stderr, "Found %d intersections\n", intersects->size());
	int n_processed_imputed_targets = 0;
	for (int i_int = 0; i_int < intersects->size(); i_int++)
	{
		if (i_int % 1000 == 0)
		{
			fprintf(stderr, "@ %d. intersect         \r", i_int);
		}

		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* imp_reg = int_info->src_reg;
		t_annot_region* known_reg = int_info->dest_reg;

		void** known_reg_info = (void**)(known_reg->data);
		char* known_reg_geno = (char*)(known_reg_info[0]);

		void** imp_reg_info = (void**)(imp_reg->data);
		char* imp_reg_geno = (char*)(imp_reg_info[0]);

		if (t_string::compare_strings(imp_reg->name, known_reg->name) &&
			imp_reg->score == 0)
		{
			imp_reg->score = 1;
			n_processed_imputed_targets++;
			
			double total_n_non_refs = 0;
			double n_matching_non_refs = 0;
			double n_matching_all = 0;
			double n_all = 0;

			int geno_i = 0;
			for (int imp_i = 0; imp_i < imputed_sample_ids->size(); imp_i++)
			{
				if (imp_2_known_sample_i->at(imp_i) < known_sample_ids->size())
				{
					double known_geno = (double)(known_reg_geno[imp_2_known_sample_i->at(imp_i)]);
					double imp_geno = (double)(imp_reg_geno[imp_i]);
					known_imp_sample_geno[0][geno_i] = known_geno;
					known_imp_sample_geno[1][geno_i] = imp_geno;
					geno_i++;

					// Update non-ref concordance.
					if (known_geno > 0)
					{
						if (known_geno == imp_geno)
						{
							n_matching_non_refs++;
						}

						total_n_non_refs++;
					}

					// Update all concordance.
					if (known_geno == imp_geno)
					{
						n_matching_all++;
					}
					n_all++;
				} // matching check.
			} // imp_i loop.			

			double cur_geno_corr = 0;
			get_correlation(known_imp_sample_geno[0], known_imp_sample_geno[1], geno_i, cur_geno_corr);

			fprintf(f_op, "%s\t%d\t%d\t%s\t%.4f\t%.0f\t%.0f\t%.0f\t%.0f\n", imp_reg->chrom,
					translate_coord(imp_reg->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
					translate_coord(imp_reg->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
					imp_reg->name, cur_geno_corr * cur_geno_corr,
					n_matching_all, n_all, 
					n_matching_non_refs, total_n_non_refs);
		} // overlapping region name comparison.
	} // i_int loop.
	fclose(f_op);

	fprintf(stderr, "\nDone.\n");
}

void get_R2_per_GIMP_4entry_allelic_probs(char* imputed_genotypes_fp, char* imputed_sample_ids_list_fp, 
											char* known_genotypes_fp, char* known_sample_ids_list_fp)
{
	vector<t_annot_region*>* known_genotype_regs = load_variant_signal_regions_wrapper(known_genotypes_fp, known_sample_ids_list_fp);
	vector<char*>* known_sample_ids = buffer_file(known_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d known genotype regions for %d samples.\n", known_genotype_regs->size(), known_sample_ids->size());

	vector<char*>* imputed_sample_ids = buffer_file(imputed_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d imputed sample ids\n", imputed_sample_ids->size());

	vector<t_annot_region*>* imputed_genotype_regs = load_BED_with_line_information(imputed_genotypes_fp);
	fprintf(stderr, "Loaded %d imputed variants.\n", imputed_genotype_regs->size());
	for (int i_reg = 0; i_reg < imputed_genotype_regs->size(); i_reg++)
	{
		imputed_genotype_regs->at(i_reg)->score = 0;
	} // i_reg loop.

	// Set the mapping between sample id's.
	vector<int>* imp_2_known_sample_i = new vector<int>();
	int n_matched_samples = 0;
	for (int imp_i = 0; imp_i < imputed_sample_ids->size(); imp_i++)
	{
		int cur_imp_known_i = t_string::get_i_str(known_sample_ids, imputed_sample_ids->at(imp_i));
		if (cur_imp_known_i < known_sample_ids->size())
		{
			n_matched_samples++;
		}
		imp_2_known_sample_i->push_back(cur_imp_known_i);
	} // imp_i loop.
	fprintf(stderr, "Matched %d samples.\n", n_matched_samples);

	fprintf(stderr, "Intersecting imputed and known variants.\n");
	FILE* f_op = open_f("R2_stats.txt", "w");
	// Intersect and process.
	double** known_imp_sample_geno = new double*[2];
	known_imp_sample_geno[0] = new double[n_matched_samples + 2];
	known_imp_sample_geno[1] = new double[n_matched_samples + 2];
	vector<t_annot_region*>* intersects = intersect_annot_regions(imputed_genotype_regs, known_genotype_regs, true);

	fprintf(stderr, "Found %d intersections\n", intersects->size());
	int n_processed_imputed_targets = 0;
	for (int i_int = 0; i_int < intersects->size(); i_int++)
	{
		if (i_int % 1000 == 0)
		{
			fprintf(stderr, "@ %d. intersect         \r", i_int);
		}

		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* imp_reg = int_info->src_reg;
		t_annot_region* known_reg = int_info->dest_reg;

		void** known_reg_info = (void**)(known_reg->data);
		char* known_reg_geno = (char*)(known_reg_info[0]);

		if (t_string::compare_strings(imp_reg->name, known_reg->name) &&
			imp_reg->score == 0)
		{
			imp_reg->score = 1;
			n_processed_imputed_targets++;

			double total_n_non_refs = 0;
			double n_matching_non_refs = 0;
			double n_matching_all = 0;
			double n_all = 0;

			int geno_i = 0;
			char* imp_var_line = (char*)(imp_reg->data);
			t_string_tokens* toks = t_string::tokenize_by_chars(imp_var_line, "\t");
			for (int imp_i = 0; imp_i < imputed_sample_ids->size(); imp_i++)
			{
				int tok_i = imp_i + 4;
				t_string_tokens* prob_toks = toks->at(tok_i)->tokenize_by_chars(";");
				if (prob_toks->size() != 4)
				{
					fprintf(stderr, "Could not parse the probability tokens @ %s on %s\n", toks->at(tok_i)->str(), imp_var_line);
					exit(0);
				}

				double hap_00_prob = atof(prob_toks->at(0)->str());
				double hap_01_prob = atof(prob_toks->at(1)->str());
				double hap_10_prob = atof(prob_toks->at(2)->str());
				double hap_11_prob = atof(prob_toks->at(3)->str());
				t_string::clean_tokens(prob_toks);

				int imp_geno = 0;
				if (hap_01_prob > hap_00_prob)
				{
					imp_geno++;
				}

				if (hap_11_prob > hap_10_prob)
				{
					imp_geno++;
				}

				if (imp_2_known_sample_i->at(imp_i) < known_sample_ids->size())
				{
					int known_geno = (int)(known_reg_geno[imp_2_known_sample_i->at(imp_i)]);
					known_imp_sample_geno[0][geno_i] = known_geno;
					known_imp_sample_geno[1][geno_i] = imp_geno;
					geno_i++;

					// Update non-ref concordance.
					if (known_geno > 0)
					{
						if (known_geno == imp_geno)
						{
							n_matching_non_refs++;
						}

						total_n_non_refs++;
					}

					// Update all concordance.
					if (known_geno == imp_geno)
					{
						n_matching_all++;
					}
					n_all++;
				}
			} // imp_i loop.

			t_string::clean_tokens(toks);

			double cur_geno_corr = 0;
			get_correlation(known_imp_sample_geno[0], known_imp_sample_geno[1], geno_i, cur_geno_corr);

			fprintf(f_op, "%s\t%d\t%d\t%s\t%.4f\t%.0f\t%.0f\t%.0f\t%.0f\n", imp_reg->chrom,
				translate_coord(imp_reg->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
				translate_coord(imp_reg->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
				imp_reg->name, cur_geno_corr * cur_geno_corr,
				n_matching_all, n_all,
				n_matching_non_refs, total_n_non_refs);
		} // overlapping region name comparison.
	} // i_int loop.
	fclose(f_op);
}