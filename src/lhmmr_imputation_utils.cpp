#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "lhmmr_variation_tools.h"
#include "lhmmr_genomics_coords.h"
#include "lhmmr_annot_region_tools.h"
#include "lhmmr_genome_sequence_tools.h"
#include "lhmmr_signal_track_tools.h"
#include "lhmmr_gff_utils.h"
#include "lhmmr_file_utils.h"
#include "lhmmr_histogram.h"
#include "lhmmr_ansi_string.h"
#include "lhmmr_nucleotide.h"
#include <string.h>
#include <ctype.h>
#include <algorithm>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "lhmmr_nomenclature.h"
#include "lhmmr_ansi_thread.h"
#include "lhmmr_rng.h"
#include "lhmmr_seed_manager.h"
#include <ctype.h>
#include <math.h>
#include <time.h>

#include "lhmmr_imputation_utils.h"
#include "lhmmr_variation_tools.h"

#ifdef __unix__
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_linalg.h>
#endif 

#define MIN(x,y) ((x) < (y)?(x):(y))
#define MAX(x,y) ((x) > (y)?(x):(y))

bool __DUMP_INPUTATION_UTILS_MSGS__ = false;

struct t_LMSE_imputation_thread_info
{
	int which;
	int outof;

	double target_normalizer;
	double tag_normalizer;

	//vector<t_annot_region*>* testing_tag_genotype_signal_regs;
	vector<t_annot_region*>* tag_genotype_signal_regs;
	//vector<t_annot_region*>* testing_target_genotype_signal_regs;
	vector<t_annot_region*>* target_genotype_signal_regs;

	t_restr_annot_region_list* restr_target_genotype_signal_regs;
	t_restr_annot_region_list* restr_tag_genotype_signal_regs;

	int testing_flag;
	int start_pos;
	int end_pos;
	char* op_dir;

	vector<char*>* tag_sample_ids;
	vector<char*>* target_sample_ids;

	vector<char*>* testing_tag_sample_ids;
	vector<char*>* testing_target_sample_ids;

	int n_tags_vars_per_side;
};

bool compare_haplotypes(double* haplo1, double* haplo2)
{
	int i_var = 0;
	while (haplo1[i_var] != -1)
	{
		if (haplo1[i_var] != haplo2[i_var])
		{
			return(false);
		}
		i_var++;
	}

	return(true);
}

bool sort_haplotypes(double* haplo1, double* haplo2)
{
	int i_var = 0;
	while (haplo1[i_var] >= 0)
	{
		if (haplo1[i_var] < haplo2[i_var])
		{
			return(true);
		}
		if (haplo1[i_var] > haplo2[i_var])
		{
			return(false);
		}

		i_var++;
	} // i_var loop.

	return(false);
}

void get_signal_level_vs_genotype_level_R2(char* imputed_genotype_signal_regs_BED_fp, char* known_genotype_regs_BED_fp, char* sample_ids_list_fp)
{
	vector<t_annot_region*>* known_geno_regs = load_variant_signal_regions_wrapper(known_genotype_regs_BED_fp, sample_ids_list_fp);
	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	fprintf(stderr, "Loaded %d known genotype regions for %d samples.\n", known_geno_regs->size(), sample_ids->size());

	int n_loaded_samples = 0;
	vector<t_annot_region*>* imp_target_geno_sig_regs = load_signal_regs_BED(imputed_genotype_signal_regs_BED_fp, n_loaded_samples);

	fprintf(stderr, "Loaded %d known genotype regions for %d samples.\n", imp_target_geno_sig_regs->size(), n_loaded_samples, sample_ids->size());

	fprintf(stderr, "Assigning signals to variant regions.\n");
	vector<t_annot_region*>* target_intersects = intersect_regions_per_names(known_geno_regs, imp_target_geno_sig_regs, true);
	for (int i_int = 0; i_int < target_intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(target_intersects->at(i_int)->data);
		t_annot_region* known_geno_reg = int_info->src_reg;
		t_annot_region* imp_sig_reg = int_info->dest_reg;

		void** known_geno_reg_info = (void**)(known_geno_reg->data);
		double* test_reg_info = (double*)(imp_sig_reg->data);
		known_geno_reg_info[1] = test_reg_info;
	} // i_int loop.

	vector<double>* sig_corrs = new vector<double>();
	vector<double>* geno_corrs = new vector<double>();
	for (int i_reg = 0; i_reg < known_geno_regs->size(); i_reg++)
	{
		void** cur_reg_info = (void**)(known_geno_regs->at(i_reg)->data);
		char* cur_known_geno = (char*)(cur_reg_info[0]);
		double* cur_imp_sig = (double*)(cur_reg_info[1]);

		if (cur_imp_sig == NULL)
		{
			continue;
		}

		double* cur_known_geno_dbl = new double[sample_ids->size() + 1];
		double* cur_impsig2geno_dbl = new double[sample_ids->size() + 1];
		for (int i_s = 0; i_s < sample_ids->size(); i_s++)
		{
			cur_known_geno_dbl[i_s] = cur_known_geno[i_s];

			double imp_sig = (cur_imp_sig[i_s]*2 / 100);
			if (imp_sig < 0.66)
			{
				imp_sig = 0;
			}
			else if (imp_sig < 1.32)
			{
				imp_sig = 1;
			}
			else
			{
				imp_sig = 2;
			}

			cur_impsig2geno_dbl[i_s] = imp_sig;
		} // i_s loop.

		double sig_corr = 0;
		double geno_corr = 0;
		get_correlation(cur_imp_sig, cur_known_geno_dbl, sample_ids->size(), sig_corr);
		get_correlation(cur_impsig2geno_dbl, cur_known_geno_dbl, sample_ids->size(), geno_corr);

		sig_corrs->push_back(sig_corr);
		geno_corrs->push_back(geno_corr);
	} // i_reg loop.

	double sig_mean, sig_std_dev;
	get_stats(sig_corrs, sig_mean, sig_std_dev);

	double geno_mean, geno_std_dev;
	get_stats(geno_corrs, geno_mean, geno_std_dev);

	fprintf(stderr, "Sig: %.4f (%.4f) ;; Geno: %.4f (%.4f)\n", sig_mean, sig_std_dev, geno_mean, geno_std_dev);
}


bool sort_haplotype_info(void** haplo1_info, void** haplo2_info)
{
	double* haplo1 = (double*)(haplo1_info[0]);
	double* haplo2 = (double*)(haplo2_info[0]);

	int i_var = 0;
	while (haplo1[i_var] >= 0)
	{
		if (haplo1[i_var] < haplo2[i_var])
		{
			return(true);
		}
		if (haplo1[i_var] > haplo2[i_var])
		{
			return(false);
		}

		i_var++;
	} // i_var loop.

	return(false);
}

bool compare_haplotype_info(void** haplo1_info, void** haplo2_info)
{
	double* haplo1 = (double*)(haplo1_info[0]);
	double* haplo2 = (double*)(haplo2_info[0]);

	int i_var = 0;
	while (haplo1[i_var] != -1)
	{
		if (haplo1[i_var] != haplo2[i_var])
		{
			return(false);
		}
		i_var++;
	}

	return(true);
}

void get_unique_haplotype_indices(vector<double*>* haplotype_alleles, vector<vector<int>*>* per_uniq_haplo_indices)
{
	vector<void**>* haplo_info = new vector<void**>();
	for (int haplo_i = 0; haplo_i < haplotype_alleles->size(); haplo_i++)
	{
		void** cur_haplo_info = new void*[2];
		cur_haplo_info[0] = haplotype_alleles->at(haplo_i);
		int* haplo_i_ptr = new int[1];
		haplo_i_ptr[0] = haplo_i;
		cur_haplo_info[1] = haplo_i_ptr;

		haplo_info->push_back(cur_haplo_info);
	} // haplo_i loop.

	sort(haplo_info->begin(), haplo_info->end(), sort_haplotype_info);

	if (__DUMP_INPUTATION_UTILS_MSGS__)
	{
		for (int hap_i = 0; hap_i < haplo_info->size(); hap_i++)
		{
			double* cur_hap = (double*)(((void**)(haplo_info->at(hap_i)))[0]);

			fprintf(stderr, "Hap %d: ", hap_i);
			int var_i = 0;
			while (1)
			{
				if (cur_hap[var_i] == -1)
				{
					break;
				}
				else
				{
					fprintf(stderr, "%d", (int)(cur_hap[var_i]));
				}
				var_i++;
			}
			fprintf(stderr, "\n");
		} // hap_i loop.
	}

	int i_hap = 0;
	while (i_hap < haplo_info->size())
	{
		int j_hap = i_hap;

		int cur_hap_cnt = 0;
		vector<int>* cur_unique_haplo_indices = new vector<int>();
		while (j_hap < haplo_info->size() &&
				compare_haplotype_info(haplo_info->at(i_hap), haplo_info->at(j_hap)))
		{
			void** cur_haplo_info = (void**)(haplo_info->at(j_hap));
			int* haplo_i_ptr = (int*)(cur_haplo_info[1]);
			cur_unique_haplo_indices->push_back(haplo_i_ptr[0]);
			j_hap++;
		} // j_hap loop.

		per_uniq_haplo_indices->push_back(cur_unique_haplo_indices);

		i_hap = j_hap;
	} // i_hap loop.

	//vector<void**>* haplo_info = new vector<void**>();
	for (int haplo_i = 0; haplo_i < haplo_info->size(); haplo_i++)
	{
		void** cur_hap_info = haplo_info->at(haplo_i);
		delete[] ((int*)(cur_hap_info[1]));

		delete[] cur_hap_info;
	} // haplo_i loop.
	delete haplo_info;
}


void count_unique_haplotypes(vector<double*>* lowMAF_containing_haplotypes, vector<int>* n_cnt_per_uniq_haplotypes)
{
	sort(lowMAF_containing_haplotypes->begin(), lowMAF_containing_haplotypes->end(), sort_haplotypes);

	if (__DUMP_INPUTATION_UTILS_MSGS__)
	{
		for (int hap_i = 0; hap_i < lowMAF_containing_haplotypes->size(); hap_i++)
		{
			fprintf(stderr, "Hap %d: ", hap_i);
			int var_i = 0;
			while (1)
			{
				if (lowMAF_containing_haplotypes->at(hap_i)[var_i] == -1)
				{
					break;
				}
				else
				{
					fprintf(stderr, "%d", (int)(lowMAF_containing_haplotypes->at(hap_i)[var_i]));
				}
				var_i++;
			}
			fprintf(stderr, "\n");
		} // hap_i loop.
	}

	int i_hap = 0;
	while (i_hap < lowMAF_containing_haplotypes->size())
	{
		int j_hap = i_hap;

		int cur_hap_cnt = 0;
		while (j_hap < lowMAF_containing_haplotypes->size() &&
			compare_haplotypes(lowMAF_containing_haplotypes->at(i_hap), lowMAF_containing_haplotypes->at(j_hap)))
		{
			cur_hap_cnt++;
			j_hap++;
		} // j_hap loop.

		n_cnt_per_uniq_haplotypes->push_back(cur_hap_cnt);

		i_hap = j_hap;
	} // i_hap loop.
}

vector<double*>* get_unique_haplotypes(vector<double*>* lowMAF_containing_haplotypes, vector<int>* n_cnt_per_uniq_haplotypes)
{
	vector<double*>* unique_lowMAF_containing_haplotypes = new vector<double*>();
	sort(lowMAF_containing_haplotypes->begin(), lowMAF_containing_haplotypes->end(), sort_haplotypes);

	int i_hap = 0;
	while(i_hap < lowMAF_containing_haplotypes->size())
	{
		int j_hap = i_hap;		

		int cur_hap_cnt = 0;
		while (j_hap < lowMAF_containing_haplotypes->size() &&
				compare_haplotypes(lowMAF_containing_haplotypes->at(i_hap), lowMAF_containing_haplotypes->at(j_hap)))
		{
			cur_hap_cnt++;
			j_hap++;
		} // j_hap loop.

		unique_lowMAF_containing_haplotypes->push_back(lowMAF_containing_haplotypes->at(i_hap));
		n_cnt_per_uniq_haplotypes->push_back(cur_hap_cnt);

		i_hap = j_hap;
	} // i_hap loop.

	return(unique_lowMAF_containing_haplotypes);
}

void is_haplo_coded(vector<t_annot_region*>* geno_regs, int n_samples, bool& found_haplocoded_genotype)
{
	for (int i_reg = 0; i_reg < geno_regs->size(); i_reg++)
	{
		void** old_info = (void**)(geno_regs->at(i_reg)->data);
		char* genos = (char*)(old_info[0]);

		for (int i_s = 0; i_s < n_samples; i_s++)
		{
			if (genos[i_s] == 3)
			{
				found_haplocoded_genotype = true;
				return;
			}
		}
	} // i_reg loop.

	found_haplocoded_genotype = false;
	return;
}

vector<t_annot_region*>* load_recombination_rates(char* recombination_rate_fp)
{
	if (!check_file(recombination_rate_fp))
	{
		return NULL;
	}

	vector<t_annot_region*>* recomb_regs = new vector<t_annot_region*>();

	FILE* f_recombination_rate = open_f(recombination_rate_fp, "r");
	char* header_line = getline(f_recombination_rate);
	fprintf(stderr, "Read header: %s\n", header_line);
	while (1)
	{
		char* cur_line = getline(f_recombination_rate);

		if (cur_line == NULL)
		{
			break;
		}

		// 16052618 9.7078062447 0.0123386217370137
		int cur_posn = 0;
		double cur_cM = 0;
		if (sscanf(cur_line, "%d %*s %lf", &cur_posn, &cur_cM) != 2)
		{
			fprintf(stderr, "Could not parse line: %s\n", cur_line);
			exit(0);
		}

		t_annot_region* cur_recomb_reg = get_empty_region();
		cur_recomb_reg->chrom = t_string::copy_me_str("XX");
		cur_recomb_reg->start = cur_posn;
		cur_recomb_reg->end = cur_posn;
		cur_recomb_reg->dbl_score = cur_cM;
		cur_recomb_reg->strand = '+';

		recomb_regs->push_back(cur_recomb_reg);

		delete[] cur_line;
	} // file reading loop.

	fprintf(stderr, "Loaded %d recombination regions.\n", recomb_regs->size());
	close_f(f_recombination_rate, recombination_rate_fp);

	return(recomb_regs);
}

double get_cumulative_recomb_rate_per_variant_optimized(t_annot_region* var_reg, vector<t_annot_region*>* sorted_recomb_regs)
{
	// If the variant is to the left of the first recomb variant, return the first recomb rate.
	if (var_reg->end < sorted_recomb_regs->at(0)->start)
	{
		return(sorted_recomb_regs->at(0)->dbl_score);
	}

	// If the region is beyond the end of all the regions, return the last region.
	if (sorted_recomb_regs->back()->end < var_reg->start)
	{
		return(sorted_recomb_regs->back()->dbl_score);
	}

	int closest_reg_i = locate_posn_region_per_region_starts(var_reg->start, sorted_recomb_regs, 0, sorted_recomb_regs->size());
	
	int i_reg = closest_reg_i;
	while (i_reg > 0 && 
		i_reg < sorted_recomb_regs->size())
	{
		if (sorted_recomb_regs->at(i_reg)->end < var_reg->start)
		{
			break;
		}

		i_reg--;
	} // i_reg loop.

	double interpolated_cM = -1;
	while (i_reg >= 0 &&
			(i_reg+1) < sorted_recomb_regs->size())
	{
		if (sorted_recomb_regs->at(i_reg)->start > var_reg->end)
		{
			break;
		}

		if (sorted_recomb_regs->at(i_reg)->start <= var_reg->start &&
			sorted_recomb_regs->at(i_reg + 1)->start >= var_reg->start)
		{
			int l_recomb_reg = sorted_recomb_regs->at(i_reg + 1)->start - sorted_recomb_regs->at(i_reg)->start;
			double delta_cM = sorted_recomb_regs->at(i_reg + 1)->dbl_score - sorted_recomb_regs->at(i_reg)->dbl_score;

			double cM_slope = delta_cM / l_recomb_reg;

			double linear_interpolated_cM = sorted_recomb_regs->at(i_reg)->dbl_score + (var_reg->start - sorted_recomb_regs->at(i_reg)->start) * cM_slope;

			//// Get the distances and scale.
			//double linear_interpolated_cM = recomb_regs->at(i_reg)->dbl_score * (double)(fabs(recomb_regs->at(i_reg + 1)->start - var_reg->start)) / l_recomb_reg +
			//								recomb_regs->at(i_reg+1)->dbl_score * (double)(fabs(recomb_regs->at(i_reg)->start - var_reg->start)) / l_recomb_reg;

			bool cur_reg_is_not_side = (sorted_recomb_regs->at(i_reg)->start != var_reg->start &&
										sorted_recomb_regs->at(i_reg + 1)->start != var_reg->start);

			if (cur_reg_is_not_side &&
				(linear_interpolated_cM < MIN(sorted_recomb_regs->at(i_reg + 1)->dbl_score, sorted_recomb_regs->at(i_reg)->dbl_score) ||
					linear_interpolated_cM > MAX(sorted_recomb_regs->at(i_reg + 1)->dbl_score, sorted_recomb_regs->at(i_reg)->dbl_score)))
			{
				fprintf(stderr, "Sanity check failed for interpolation.\n");
				exit(0);
			}

			interpolated_cM = linear_interpolated_cM;
			break;
		} // overlap check.

		i_reg++;
	} // i_reg loop.

	// Following is brute-force search.
	bool CHECK_RECOMB_VAL = false;
	if (CHECK_RECOMB_VAL)
	{
		for (int i_reg = 0; i_reg < sorted_recomb_regs->size() - 1; i_reg++)
		{
			if (sorted_recomb_regs->at(i_reg)->start <= var_reg->start &&
				sorted_recomb_regs->at(i_reg + 1)->start >= var_reg->start)
			{
				int l_recomb_reg = sorted_recomb_regs->at(i_reg + 1)->start - sorted_recomb_regs->at(i_reg)->start;
				double delta_cM = sorted_recomb_regs->at(i_reg + 1)->dbl_score - sorted_recomb_regs->at(i_reg)->dbl_score;

				double cM_slope = delta_cM / l_recomb_reg;

				double linear_interpolated_cM = sorted_recomb_regs->at(i_reg)->dbl_score + (var_reg->start - sorted_recomb_regs->at(i_reg)->start) * cM_slope;

				//// Get the distances and scale.
				//double linear_interpolated_cM = recomb_regs->at(i_reg)->dbl_score * (double)(fabs(recomb_regs->at(i_reg + 1)->start - var_reg->start)) / l_recomb_reg +
				//								recomb_regs->at(i_reg+1)->dbl_score * (double)(fabs(recomb_regs->at(i_reg)->start - var_reg->start)) / l_recomb_reg;

				if (linear_interpolated_cM < MIN(sorted_recomb_regs->at(i_reg + 1)->dbl_score, sorted_recomb_regs->at(i_reg)->dbl_score) ||
					linear_interpolated_cM > MAX(sorted_recomb_regs->at(i_reg + 1)->dbl_score, sorted_recomb_regs->at(i_reg)->dbl_score))
				{
					fprintf(stderr, "Sanity check failed for interpolation.\n");
					exit(0);
				}

				if (interpolated_cM != linear_interpolated_cM)
				{
					fprintf(stderr, "Recomb rate is assigned incorrectly for %s:%d: %.5f, %.5f\n",
						var_reg->chrom, var_reg->start,
						interpolated_cM, linear_interpolated_cM);
					exit(0);
				}
				break;
			} // overlap check
		} // i_reg loop.
	} // check for assignment.

	return(interpolated_cM);
}

double get_cumulative_recomb_rate_per_variant(t_annot_region* var_reg, vector<t_annot_region*>* recomb_regs)
{
	for (int i_reg = 0; i_reg < recomb_regs->size()-1; i_reg++)
	{
		if (recomb_regs->at(i_reg)->start <= var_reg->start &&
			recomb_regs->at(i_reg+1)->start >= var_reg->start)
		{
			int l_recomb_reg = recomb_regs->at(i_reg + 1)->start - recomb_regs->at(i_reg)->start;
			double delta_cM = recomb_regs->at(i_reg+1)->dbl_score - recomb_regs->at(i_reg)->dbl_score;

			double cM_slope = delta_cM / l_recomb_reg;

			double linear_interpolated_cM = recomb_regs->at(i_reg)->dbl_score + (var_reg->start - recomb_regs->at(i_reg)->start) * cM_slope;

			//// Get the distances and scale.
			//double linear_interpolated_cM = recomb_regs->at(i_reg)->dbl_score * (double)(fabs(recomb_regs->at(i_reg + 1)->start - var_reg->start)) / l_recomb_reg +
			//								recomb_regs->at(i_reg+1)->dbl_score * (double)(fabs(recomb_regs->at(i_reg)->start - var_reg->start)) / l_recomb_reg;

			bool cur_reg_is_not_side = (recomb_regs->at(i_reg)->start != var_reg->start &&
										recomb_regs->at(i_reg + 1)->start != var_reg->start);

			if (cur_reg_is_not_side &&
				(linear_interpolated_cM < MIN(recomb_regs->at(i_reg + 1)->dbl_score, recomb_regs->at(i_reg)->dbl_score) ||
				linear_interpolated_cM > MAX(recomb_regs->at(i_reg + 1)->dbl_score, recomb_regs->at(i_reg)->dbl_score)))
			{
				fprintf(stderr, "Sanity check failed for interpolation: %lf, %lf, %lf\n",
					recomb_regs->at(i_reg + 1)->dbl_score, linear_interpolated_cM, recomb_regs->at(i_reg)->dbl_score);
				exit(0);
			}

			return(linear_interpolated_cM);
		}
	}

	return(recomb_regs->back()->dbl_score);
}

void get_child_geno_scores_per_parent_geno(char par1_geno, char par2_geno, double* per_geno_scores)
{
	// Reset the scores.
	per_geno_scores[0] = 0;
	per_geno_scores[1] = 0;
	per_geno_scores[2] = 0;

	// 00
	if (par1_geno == 0 && par2_geno == 0)
	{
		per_geno_scores[0] = 1.0;
	}

	//01
	if (par1_geno == 0 && par2_geno == 1)
	{
		per_geno_scores[0] = 0.5;
		per_geno_scores[1] = 0.5;
	}

	//02
	if (par1_geno == 0 && par2_geno == 2)
	{
		per_geno_scores[1] = 1.0;
	}

	//10
	if (par1_geno == 1 && par2_geno == 0)
	{
		per_geno_scores[0] = 0.5;
		per_geno_scores[1] = 0.5;
	}

	//11
	if (par1_geno == 1 && par2_geno == 1)
	{
		per_geno_scores[0] = 0.25;
		per_geno_scores[1] = .5;
		per_geno_scores[2] = .25;
	}

	//12
	if (par1_geno == 1 && par2_geno == 2)
	{
		per_geno_scores[1] = 0.5;
		per_geno_scores[2] = 0.5;
	}

	//20
	if (par1_geno == 2 && par2_geno == 0)
	{
		per_geno_scores[1] = 1.0;
	}
	//21
	if (par1_geno == 2 && par2_geno == 1)
	{
		per_geno_scores[1] = 0.5;
		per_geno_scores[2] = 0.5;
	}

	//22
	if (par1_geno == 2 && par2_geno == 2)
	{
		per_geno_scores[2] = 1;
	}

	double tot_geno_scores = per_geno_scores[0] + per_geno_scores[1] + per_geno_scores[2];
	if (tot_geno_scores > 1.01 || 
		tot_geno_scores < 0.99)
	{
		fprintf(stderr, "Sanity check failed: Total prior genotype probabilities sum up to %lf\n", tot_geno_scores);
	}
	
	return;
}

void update_parents_2_child_imputed_genotype_probabilities(char* trio_list_fp, char* genocoded_genotype_matrix_fp, char* sample_ids_list_fp,
															char* imputed_geno_signals_matrix_fp, char* imputed_signal_sample_ids_list_fp)
{
	fprintf(stderr, "Updating the imputed genotype probabilities.\n");
	vector<t_annot_region*>* genocoded_geno_regs = load_variant_signal_regions_wrapper(genocoded_genotype_matrix_fp, sample_ids_list_fp);
	vector<char*>* genotype_sample_ids = buffer_file(sample_ids_list_fp);
	fprintf(stderr, "Loaded %d regions with %d samples.\n", genocoded_geno_regs->size(), genotype_sample_ids->size());

	// 
	int n_loaded_samples = 0;
	vector<t_annot_region*>* imp_target_geno_sig_regs = load_signal_regs_BED(imputed_geno_signals_matrix_fp, n_loaded_samples);
	vector<char*>* imputed_sample_ids = buffer_file(imputed_signal_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d samples from %s\n", imputed_sample_ids->size(), imputed_signal_sample_ids_list_fp);

	if (n_loaded_samples != imputed_sample_ids->size())
	{
		fprintf(stderr, "Could not match the number of loaded samples in imputed signal matrix with the sample list.\n");
		exit(0);
	}

	for (int i_reg = 0; i_reg < imp_target_geno_sig_regs->size(); i_reg++)
	{
		void** cur_reg_info = new void*[10];
		cur_reg_info[0] = imp_target_geno_sig_regs->at(i_reg)->data;
		cur_reg_info[1] = NULL;
		cur_reg_info[2] = NULL;

		imp_target_geno_sig_regs->at(i_reg)->data = cur_reg_info;
	} // i_reg loop.

	// Assign the imputed regions to the genocoded regions.
	fprintf(stderr, "Assigning the imputed signals to known genotypes.\n");
	vector<t_annot_region*>* intersects = intersect_annot_regions(imp_target_geno_sig_regs, genocoded_geno_regs, true);
	fprintf(stderr, "Processing %d intersections.\n", intersects->size());
	for (int i_int = 0; i_int < intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* sig_reg = int_info->src_reg;
		t_annot_region* geno_reg = int_info->dest_reg;

		if (t_string::compare_strings(sig_reg->name, geno_reg->name))
		{
			void** cur_sig_reg_sigs = (void**)(sig_reg->data);
			void** cur_geno_reg_info = (void**)(geno_reg->data);
			cur_sig_reg_sigs[1] = cur_geno_reg_info[0]; // Assign the char* here.
		}
	} // i_int loop.

	fprintf(stderr, "Loaded %d imputed target regions for %d samples.\n", imp_target_geno_sig_regs->size(), imputed_sample_ids->size());

	vector<char*>* trio_lines = buffer_file(trio_list_fp);
	fprintf(stderr, "Processing %d trios.\n", trio_lines->size());

	for (int i_trio = 0; i_trio < trio_lines->size(); i_trio++)
	{
		fprintf(stderr, "Processing trio: %s\n", trio_lines->at(i_trio));
		char parent1_id[100];
		char parent2_id[100];
		char child_id[100];
		if (sscanf(trio_lines->at(i_trio), "%s %s %s", parent1_id, parent2_id, child_id) != 3)
		{
			fprintf(stderr, "Could not parse trio sample id's.\n");
			exit(0);
		}

		int par1_sample_i = t_string::get_i_str(genotype_sample_ids, parent1_id);
		int par2_sample_i = t_string::get_i_str(genotype_sample_ids, parent2_id);
		int child_known_geno_sample_i = t_string::get_i_str(genotype_sample_ids, child_id);
		int child_imp_signal_sample_i = t_string::get_i_str(imputed_sample_ids, child_id);

		if (par1_sample_i == genotype_sample_ids->size() ||
			par2_sample_i == genotype_sample_ids->size() ||
			child_known_geno_sample_i == genotype_sample_ids->size() ||
			child_imp_signal_sample_i == imputed_sample_ids->size())
		{
			fprintf(stderr, "Could not locate one of the members of the trio: %s (%s:%d, %s:%d, %s:%d,%d)\n", trio_lines->at(i_trio),
					parent1_id, par1_sample_i,
					parent2_id, par2_sample_i,
					child_id, child_known_geno_sample_i, child_imp_signal_sample_i);
			exit(0);
		}

		// Start comparing the known genotypes of parents and accumulate to the child:
		double** ML_confusion_matrix = new double*[3];
		double** MAP_confusion_matrix = new double*[3];
		for (int geno = 0; geno < 3; geno++)
		{
			ML_confusion_matrix[geno] = new double[3];
			MAP_confusion_matrix[geno] = new double[3];
			memset(ML_confusion_matrix[geno], 0, sizeof(double) * 3);
			memset(MAP_confusion_matrix[geno], 0, sizeof(double) * 3);
		} // geno loop.
		
		double* known_geno_sig = new double[imp_target_geno_sig_regs->size() + 2];
		double* MAP_geno_sig = new double[imp_target_geno_sig_regs->size() + 2];
		double* ML_geno_sig = new double[imp_target_geno_sig_regs->size() + 2];
		for (int i_var = 0; i_var < imp_target_geno_sig_regs->size(); i_var++)
		{
			void** cur_var_info = (void**)(imp_target_geno_sig_regs->at(i_var)->data);
			double* imp_target_signal_val = (double*)(cur_var_info[0]);
			double child_imp_signal_val = imp_target_signal_val[child_imp_signal_sample_i];

			char* geno_signals = (char*)(cur_var_info[1]);
			char par1_geno = geno_signals[par1_sample_i];
			char par2_geno = geno_signals[par2_sample_i];
			char child_known_geno = geno_signals[child_known_geno_sample_i];


			known_geno_sig[i_var] = child_known_geno;

			double per_geno_scores[10];
			get_child_geno_scores_per_parent_geno(par1_geno, par2_geno, per_geno_scores);
			
			// For the current sample, get the child's prior genotype probabilities.
			double tot_weighted_probs = 0;
			double ML_geno_prob = -1;
			double ML_geno = -1;
			double per_geno_LH[3];
			for (int geno = 0; geno < 3; geno++)
			{
				double cur_LH = exp(-1 * fabs(2 * (child_imp_signal_val / 100) - geno));
				per_geno_LH[geno] = cur_LH;

				// Update the numerator of posterior.
				tot_weighted_probs += cur_LH * per_geno_scores[geno];

				if (cur_LH > ML_geno_prob)
				{
					ML_geno_prob = cur_LH;
					ML_geno = geno;
				}
			} // geno loop.

			// Integrate the posterior.

			// Compute the probabilities.
			double MAP_geno = -1;
			double MAP_geno_prob = -1;
			double post_probs[3];
			for (int geno = 0; geno < 3; geno++)
			{
				double cur_LH = per_geno_LH[geno];
				post_probs[geno] = per_geno_scores[geno] * cur_LH / tot_weighted_probs;

				if (MAP_geno_prob < post_probs[geno])
				{
					MAP_geno_prob = post_probs[geno];
					MAP_geno = geno;
				}
			} // geno loop.

			MAP_geno_sig[i_var] = MAP_geno;
			ML_geno_sig[i_var] = ML_geno;

			  // Update the confusion matrix.
			MAP_confusion_matrix[(int)MAP_geno][child_known_geno]++;
			ML_confusion_matrix[(int)ML_geno][child_known_geno]++;
		} // i_var loop.

		// Dump the confusion matrix.
		fprintf(stderr, "Finished processing trio: %s, %s, %s\n", parent1_id, parent2_id, child_id);
		fprintf(stderr, "MAP Confusion Matrix:\n");
		double tot_MAP_corr = 0;
		double tot_MAP = 0;
		for (int pred_geno = 0; pred_geno < 3; pred_geno++)
		{
			fprintf(stderr, "Imp %d:", pred_geno);
			for (int known_geno = 0; known_geno < 3; known_geno++)
			{
				fprintf(stderr, "\t%.0f", MAP_confusion_matrix[pred_geno][known_geno]);

				if (known_geno == pred_geno)
				{
					tot_MAP_corr += MAP_confusion_matrix[pred_geno][known_geno];
				}

				tot_MAP += MAP_confusion_matrix[pred_geno][known_geno];
			} // known_geno loop.

			fprintf(stderr, "\n");
		} // pred_geno loop.

		fprintf(stderr, "ML Confusion Matrix:\n");
		double tot_ML_corr = 0;
		double tot_ML = 0;
		for (int pred_geno = 0; pred_geno < 3; pred_geno++)
		{
			fprintf(stderr, "Imp %d:", pred_geno);
			for (int known_geno = 0; known_geno < 3; known_geno++)
			{
				fprintf(stderr, "\t%.0f", ML_confusion_matrix[pred_geno][known_geno]);

				if (known_geno == pred_geno)
				{
					tot_ML_corr += ML_confusion_matrix[pred_geno][known_geno];
				}

				tot_ML += ML_confusion_matrix[pred_geno][known_geno];
			} // known_geno loop.

			fprintf(stderr, "\n");
		} // pred_geno loop.

		double MAP_corr = 0;
		get_correlation(known_geno_sig, MAP_geno_sig, imp_target_geno_sig_regs->size(), MAP_corr);
		double ML_corr = 0;
		get_correlation(known_geno_sig, ML_geno_sig, imp_target_geno_sig_regs->size(), ML_corr);

		fprintf(stderr, "R^2: MAP: %.4f;; ML: %.4f (MAP=%lf, ML=%lf)\n", MAP_corr * MAP_corr, ML_corr * ML_corr,
				tot_MAP_corr / tot_MAP, tot_ML_corr / tot_ML);
	} // i_trio loop.
}

static void* resampling_thread_callback(void* thread_info_ptr)
{
	void** thread_ptrs_list = (void**)(thread_info_ptr);

	t_restr_annot_region_list* restr_geno_regs = (t_restr_annot_region_list*)(thread_ptrs_list[0]);

	char* recombination_rate_dir = (char*)(thread_ptrs_list[1]);

	int* thread_i_ptr = (int*)(thread_ptrs_list[2]);
	int thread_i = thread_i_ptr[0];

	int* n_threads_ptr = (int*)(thread_ptrs_list[3]);
	int n_threads = n_threads_ptr[0];

	double* N_e_ptr = (double*)(thread_ptrs_list[4]);
	double N_e = N_e_ptr[0];

	double* allele_error_ptr = (double*)(thread_ptrs_list[5]);
	double allele_error = allele_error_ptr[0];

	int* n_original_haps_ptr = (int*)(thread_ptrs_list[6]);
	int n_original_haps = n_original_haps_ptr[0];

	int* n_resampled_sample_size_ptr = (int*)(thread_ptrs_list[7]);
	int n_resampled_sample_size = n_resampled_sample_size_ptr[0];

	// Load the recombination rates.
	int cur_thread_seed = time(NULL) + thread_i;
	t_rng* rng = new t_rng(cur_thread_seed);

	fprintf(stderr, "Setting thread %d/%d: N_e: %d, error_prob: %.5f, n_ref_haps: %d, n_resampled_size: %d;; random_seed: %d\n",
			thread_i,
			n_threads, (int)N_e, allele_error, n_original_haps, n_resampled_sample_size, cur_thread_seed);

	fprintf(stderr, "First 5 random numbers: ");
	for (int r_i = 0; r_i < 5; r_i++)
	{
		fprintf(stderr, "%.4f, ", rng->random_double_ran3());
	} // r_i loop.
	fprintf(stderr, "\n");

	for (int i_chr = 0; i_chr < restr_geno_regs->chr_ids->size(); i_chr++)
	{
		//fprintf(stderr, "Re-sampling variants on %s\n", restr_geno_regs->chr_ids->at(i_chr));
		vector<t_annot_region*>* cur_chrom_var_regs = restr_geno_regs->regions_per_chrom[i_chr];

		/*char cur_chr_recombination_rate_fp[1000];
		sprintf(cur_chr_recombination_rate_fp, "%s/%s.map", recombination_rate_dir, restr_geno_regs->chr_ids->at(i_chr));
		vector<t_annot_region*>* cur_chrom_recomb_regs = load_recombination_rates(cur_chr_recombination_rate_fp);
		if (cur_chrom_recomb_regs == NULL)
		{
			fprintf(stderr, "Could not load recombination rates from %s\n", cur_chr_recombination_rate_fp);
			exit(0);
		}*/

		//// Assign the recomb rates.
		//fprintf(stderr, "Setting recombination rates.\n");
		//for (int i_reg = 0; i_reg < cur_chrom_var_regs->size(); i_reg++)
		//{
		//	double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant(cur_chrom_var_regs->at(i_reg), cur_chrom_recomb_regs);
		//	cur_chrom_var_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
		//} // i_reg loop.

		//for (int i_reg = 0; i_reg < cur_chrom_var_regs->size(); i_reg++)
		//{
		//	void** cur_reg_info = (void**)(cur_chrom_var_regs->at(i_reg)->data);
		//	char* cur_reg_resampled_geno = new char[n_resampled_sample_size + 2];
		//	memset(cur_reg_resampled_geno, 0, sizeof(char) * (n_resampled_sample_size + 2));
		//	cur_reg_info[1] = cur_reg_resampled_geno;
		//} // i_reg loop.

		for (int i_s = 0; i_s < n_resampled_sample_size; i_s++)
		{
			if (i_s % n_threads != thread_i)
			{
				continue;
			}

			// Sample the two haplotypes.
			for (int i_hap = 0; i_hap < 2; i_hap++)
			{
				int cur_haplo_state = -1;
				int n_recombs = 0;
				for (int i_reg = 0; i_reg < cur_chrom_var_regs->size(); i_reg++)
				{
					if (i_reg == 0)
					{
						cur_haplo_state = MIN((n_original_haps - 1), floor(rng->random_double_ran3() * n_original_haps));
					}
					else
					{
						// Get the recombination rate between these positions.
						double r_m = fabs(cur_chrom_var_regs->at(i_reg)->dbl_score - cur_chrom_var_regs->at(i_reg - 1)->dbl_score);
						double rho_m = 4 * N_e * r_m;

						double tau_m = 1 - exp(-1 * rho_m / n_original_haps);

						double other_prob = tau_m / n_original_haps;
						double self_prob = (1 - tau_m) + (tau_m / n_original_haps);

						double rand_cumul_val = rng->random_double_ran3();

						if (__DUMP_INPUTATION_UTILS_MSGS__)
						{
							fprintf(stderr, "Sample %d: Var %d (%s:%d); Haplo_state: %d; dcM:%.4f (%.5f - %.5f), other_prob=%.3f, self_prob=%.3f\n",
								i_s, i_reg,
								cur_chrom_var_regs->at(i_reg)->chrom, cur_chrom_var_regs->at(i_reg)->start,
								cur_haplo_state, r_m,
								cur_chrom_var_regs->at(i_reg)->dbl_score, cur_chrom_var_regs->at(i_reg - 1)->dbl_score,
								other_prob, self_prob);
						}

						double cur_cumul_val = 0;
						for (int hap_state_i = 0; hap_state_i < n_original_haps; hap_state_i++)
						{
							// Update the cumulative.
							if (hap_state_i == cur_haplo_state)
							{
								cur_cumul_val += self_prob;
							}
							else
							{
								cur_cumul_val += other_prob;
							}

							// Check the cumulative.
							if (cur_cumul_val > rand_cumul_val)
							{
								if (cur_haplo_state != hap_state_i)
								{
									n_recombs++;
								}

								cur_haplo_state = hap_state_i;
								break;
							}
						} // hap_state_i loop.
					} // region check for initing the haplo state.

					  // Copy the allele.
					int sample_i = (cur_haplo_state - (cur_haplo_state % 2)) / 2;
					int sampled_hap_i = (cur_haplo_state % 2);
					void** cur_reg_info = (void**)(cur_chrom_var_regs->at(i_reg)->data);
					char* cur_reg_geno_sig = (char*)(cur_reg_info[0]);
					char* cur_reg_sampled_geno = (char*)(cur_reg_info[1]);
					int cur_allele = get_allele_per_haplotype(cur_reg_geno_sig[sample_i], sampled_hap_i);

					// Copy the allele to the haplotype 
					char cur_val = (char)(cur_reg_sampled_geno[i_s]);
					cur_reg_sampled_geno[i_s] = cur_val | (cur_allele << i_hap);
				} // i_reg loop.

				if ((i_s - thread_i) % 5000 == 0)
				{
					fprintf(stderr, "Thread %d: Re-sampled sample %d (%d): %d recombinations.\n", thread_i, i_s, i_hap, n_recombs);
				}
			} // i_hap loop.
		} // i_s loop.
	} // i_chr loop.
}

void resample_phased_haplotypes_per_recombination_rates_multithreaded(char* haplocoded_genotype_matrix_fp,
	char* sample_ids_list_fp,
	char* recombination_rate_dir,
	double upsampling_rate,
	double N_e,
	double allele_error_prob,
	int n_threads,
	char* op_fp)
{
	fprintf(stderr, "%d-thread Re-sampling recombination-aware genotypes using N_e=%d and allelic error=%.6f\n", n_threads, (int)N_e, allele_error_prob);
	vector<t_annot_region*>* haplocoded_geno_regs = load_variant_signal_regions_wrapper(haplocoded_genotype_matrix_fp, sample_ids_list_fp);

	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	fprintf(stderr, "Loaded %d haplocoded variants for %d individuals.\n", haplocoded_geno_regs->size(), sample_ids->size());
	double n_original_haps = 2*sample_ids->size();

	t_restr_annot_region_list* restr_geno_regs = restructure_annot_regions(haplocoded_geno_regs);

	int n_resampled_sample_size = (int)(upsampling_rate * sample_ids->size());
	fprintf(stderr, "Re-Sampling %d samples.\n", n_resampled_sample_size);

	for (int i_chr = 0; i_chr < restr_geno_regs->chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Re-sampling variants on %s\n", restr_geno_regs->chr_ids->at(i_chr));
		vector<t_annot_region*>* cur_chrom_var_regs = restr_geno_regs->regions_per_chrom[i_chr];

		char cur_chr_recombination_rate_fp[1000];
		sprintf(cur_chr_recombination_rate_fp, "%s/%s.map", recombination_rate_dir, restr_geno_regs->chr_ids->at(i_chr));
		vector<t_annot_region*>* cur_chrom_recomb_regs = load_recombination_rates(cur_chr_recombination_rate_fp);
		if (cur_chrom_recomb_regs == NULL)
		{
			fprintf(stderr, "Could not load recombination rates from %s\n", cur_chr_recombination_rate_fp);
			exit(0);
		}

		// Assign the recomb rates.
		fprintf(stderr, "Setting recombination rates.\n");
		for (int i_reg = 0; i_reg < cur_chrom_var_regs->size(); i_reg++)
		{
			double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant(cur_chrom_var_regs->at(i_reg), cur_chrom_recomb_regs);
			cur_chrom_var_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
		} // i_reg loop.

		for (int i_reg = 0; i_reg < cur_chrom_var_regs->size(); i_reg++)
		{
			void** cur_reg_info = (void**)(cur_chrom_var_regs->at(i_reg)->data);
			char* cur_reg_resampled_geno = new char[n_resampled_sample_size + 2];
			memset(cur_reg_resampled_geno, 0, sizeof(char) * (n_resampled_sample_size + 2));
			cur_reg_info[1] = cur_reg_resampled_geno;
		} // i_reg loop.
	} // Assign recombination rates

	vector<t_ansi_thread*>* resampling_threads = new vector<t_ansi_thread*>();
	for (int thread_i = 0; thread_i < n_threads; thread_i++)
	{
		fprintf(stderr, "Starting %d. re-sampling thread..\n", thread_i);

		void** thread_ptrs_list = new void*[10];

		thread_ptrs_list[0] = restr_geno_regs;

		thread_ptrs_list[1] = recombination_rate_dir;

		int* thread_i_ptr = new int[2];
		thread_i_ptr[0] = thread_i;
		thread_ptrs_list[2] = thread_i_ptr;

		int* n_threads_ptr = new int[2];
		n_threads_ptr[0] = n_threads;
		thread_ptrs_list[3] = n_threads_ptr;

		double* N_e_ptr = new double[2];
		N_e_ptr[0] = N_e;
		thread_ptrs_list[4] = N_e_ptr;

		double* allele_error_ptr = new double[2];
		allele_error_ptr[0] = allele_error_prob;
		thread_ptrs_list[5] = allele_error_ptr;

		int* n_original_haps_ptr = new int[2];
		n_original_haps_ptr[0] = n_original_haps;
		thread_ptrs_list[6] = n_original_haps_ptr;

		int* n_resampled_sample_size_ptr = new int[2];
		n_resampled_sample_size_ptr[0] = n_resampled_sample_size;
		thread_ptrs_list[7] = n_resampled_sample_size_ptr;

		t_ansi_thread* cur_thread = new t_ansi_thread(resampling_thread_callback, thread_ptrs_list);
		cur_thread->run_thread();

		resampling_threads->push_back(cur_thread);
	} // thread_i loop.

	fprintf(stderr, "Stared %d/%d threads; waiting.\n", resampling_threads->size(), n_threads);

	for (int thread_i = 0; thread_i < resampling_threads->size(); thread_i++)
	{
		resampling_threads->at(thread_i)->wait_thread();
		fprintf(stderr, "%d. thread finished.\n", thread_i);
	} // thread_i waiting loop.

	fprintf(stderr, "Saving re-sampled genotypes.\n");

	// Save the genotypes.
	FILE* f_op = open_f(op_fp, "w");
	for (int i_reg = 0; i_reg < haplocoded_geno_regs->size(); i_reg++)
	{
		void** cur_reg_info = (void**)(haplocoded_geno_regs->at(i_reg)->data);
		char* cur_resampled_geno_sig = (char*)(cur_reg_info[1]);

		fprintf(f_op, "%s\t%d\t%d\t%s", haplocoded_geno_regs->at(i_reg)->chrom,
			translate_coord(haplocoded_geno_regs->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			translate_coord(haplocoded_geno_regs->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			haplocoded_geno_regs->at(i_reg)->name);

		for (int i_s = 0; i_s < n_resampled_sample_size; i_s++)
		{
			fprintf(f_op, "\t%d", (int)(cur_resampled_geno_sig[i_s]));
		} // i_s loop.

		fprintf(f_op, "\n");
	} // i_reg loop.
	fclose(f_op);

	FILE* f_samples = open_f("sample_ids.list", "w");
	for (int i_s = 0; i_s < n_resampled_sample_size; i_s++)
	{
		fprintf(f_samples, "sample_%d\n", i_s);
	} // i_s loop.
	fclose(f_samples);
}

void resample_phased_haplotypes_per_recombination_rates(char* haplocoded_genotype_matrix_fp, 
														char* sample_ids_list_fp, 
														char* recombination_rate_dir, 
														double upsampling_rate,
														double N_e,
														double allele_error_prob,
														char* op_fp)
{
	//double N_e = 10 ^ 6;
	//double allele_error_prob = pow(10, -4);

	fprintf(stderr, "Re-sampling recombination-aware genotypes using N_e=%d and allelic error=%.6f\n", (int)N_e, allele_error_prob);
	vector<t_annot_region*>* haplocoded_geno_regs = load_variant_signal_regions_wrapper(haplocoded_genotype_matrix_fp, sample_ids_list_fp);

	vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
	fprintf(stderr, "Loaded %d haplocoded variants for %d individuals.\n", haplocoded_geno_regs->size(), sample_ids->size());
	double n_original_haps = 2*sample_ids->size();

	t_restr_annot_region_list* restr_geno_regs = restructure_annot_regions(haplocoded_geno_regs);

	int n_resampled_sample_size = (int)(upsampling_rate * sample_ids->size());
	fprintf(stderr, "Re-Sampling %d samples.\n", n_resampled_sample_size);

	// Load the recombination rates.
	t_rng* rng = new t_rng(t_seed_manager::seed_me());
	for(int i_chr = 0; i_chr < restr_geno_regs->chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Re-sampling variants on %s\n", restr_geno_regs->chr_ids->at(i_chr));
		vector<t_annot_region*>* cur_chrom_var_regs = restr_geno_regs->regions_per_chrom[i_chr];	

		char cur_chr_recombination_rate_fp[1000];
		sprintf(cur_chr_recombination_rate_fp, "%s/%s.map", recombination_rate_dir, restr_geno_regs->chr_ids->at(i_chr));
		vector<t_annot_region*>* cur_chrom_recomb_regs = load_recombination_rates(cur_chr_recombination_rate_fp);
		if (cur_chrom_recomb_regs == NULL)
		{
			fprintf(stderr, "Could not load recombination rates from %s\n", cur_chr_recombination_rate_fp);
			exit(0);
		}
		
		// Assign the recomb rates.
		fprintf(stderr, "Setting recombination rates.\n");
		for (int i_reg = 0; i_reg < cur_chrom_var_regs->size(); i_reg++)
		{
			double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant(cur_chrom_var_regs->at(i_reg), cur_chrom_recomb_regs);
			cur_chrom_var_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
		} // i_reg loop.

		for (int i_reg = 0; i_reg < cur_chrom_var_regs->size(); i_reg++)
		{
			void** cur_reg_info = (void**)(cur_chrom_var_regs->at(i_reg)->data);
			char* cur_reg_resampled_geno = new char[n_resampled_sample_size + 2];
			memset(cur_reg_resampled_geno, 0, sizeof(char) * (n_resampled_sample_size + 2));
			cur_reg_info[1] = cur_reg_resampled_geno;
		} // i_reg loop.

		for (int i_s = 0; i_s < n_resampled_sample_size; i_s++)
		{
			// Sample the two haplotypes.
			for (int i_hap = 0; i_hap < 2; i_hap++)
			{
				fprintf(stderr, "Re-sampling sample %d (%d)\n", i_s, i_hap);				

				int cur_haplo_state = -1;
				int n_recombs = 0;
				for (int i_reg = 0; i_reg < cur_chrom_var_regs->size(); i_reg++)
				{
					if (i_reg == 0)
					{
						cur_haplo_state = MIN((n_original_haps - 1), floor(rng->random_double_ran3() * n_original_haps));
					}
					else
					{
						// Get the recombination rate between these positions.
						double r_m = fabs(cur_chrom_var_regs->at(i_reg)->dbl_score - cur_chrom_var_regs->at(i_reg -1 )->dbl_score);
						double rho_m = 4 * N_e * r_m;

						double tau_m = 1 - exp(-1 * rho_m / n_original_haps);

						double other_prob = tau_m / n_original_haps;
						double self_prob = (1 - tau_m) + (tau_m / n_original_haps);

						double rand_cumul_val = rng->random_double_ran3();

						if (__DUMP_INPUTATION_UTILS_MSGS__)
						{
							fprintf(stderr, "Sample %d: Var %d (%s:%d); Haplo_state: %d; dcM:%.4f (%.5f - %.5f), other_prob=%.3f, self_prob=%.3f\n",
								i_s, i_reg,
								cur_chrom_var_regs->at(i_reg)->chrom, cur_chrom_var_regs->at(i_reg)->start,
								cur_haplo_state, r_m,
								cur_chrom_var_regs->at(i_reg)->dbl_score, cur_chrom_var_regs->at(i_reg - 1)->dbl_score,
								other_prob, self_prob);
						}

						double cur_cumul_val = 0;
						for (int hap_state_i = 0; hap_state_i < n_original_haps; hap_state_i++)
						{
							// Update the cumulative.
							if (hap_state_i == cur_haplo_state)
							{
								cur_cumul_val += self_prob;
							}
							else
							{
								cur_cumul_val += other_prob;
							}
							
							// Check the cumulative.
							if (cur_cumul_val > rand_cumul_val)
							{
								if (cur_haplo_state != hap_state_i)
								{
									n_recombs++;
								}

								cur_haplo_state = hap_state_i;
								break;
							}
						} // hap_state_i loop.
					} // region check for initing the haplo state.

					// Copy the allele.
					int sample_i = (cur_haplo_state - (cur_haplo_state % 2)) / 2;
					int sampled_hap_i = (cur_haplo_state % 2);
					void** cur_reg_info = (void**)(cur_chrom_var_regs->at(i_reg)->data);
					char* cur_reg_geno_sig = (char*)(cur_reg_info[0]);
					char* cur_reg_sampled_geno = (char*)(cur_reg_info[1]);
					int cur_allele = get_allele_per_haplotype(cur_reg_geno_sig[sample_i], sampled_hap_i);

					// Copy the allele to the haplotype 
					char cur_val = (char)(cur_reg_sampled_geno[i_s]);
					cur_reg_sampled_geno[i_s] = cur_val | (cur_allele << i_hap);
				} // i_reg loop.

				fprintf(stderr, "Re-sampling sample %d (%d): %d recombinations.\n", i_s, i_hap, n_recombs);
			} // i_hap loop.
 		} // i_s loop.
	} // i_chr loop.

	fprintf(stderr, "Saving re-sampled genotypes.\n");

	// Save the genotypes.
	FILE* f_op = open_f(op_fp, "w");
	for (int i_reg = 0; i_reg < haplocoded_geno_regs->size(); i_reg++)
	{
		void** cur_reg_info = (void**)(haplocoded_geno_regs->at(i_reg)->data);
		char* cur_resampled_geno_sig = (char*)(cur_reg_info[1]);

		fprintf(f_op, "%s\t%d\t%d\t%s", haplocoded_geno_regs->at(i_reg)->chrom, 
				translate_coord(haplocoded_geno_regs->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base), 
				translate_coord(haplocoded_geno_regs->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
				haplocoded_geno_regs->at(i_reg)->name);

		for (int i_s = 0; i_s < n_resampled_sample_size; i_s++)
		{
			fprintf(f_op, "\t%d", (int)(cur_resampled_geno_sig[i_s]));
		} // i_s loop.

		fprintf(f_op, "\n");
	} // i_reg loop.
	fclose(f_op);

	FILE* f_samples = open_f("sample_ids.list", "w");
	for (int i_s = 0; i_s < n_resampled_sample_size; i_s++)
	{
		fprintf(f_samples, "sample_%d\n", i_s);
	} // i_s loop.
	fclose(f_samples);
}

/*
void assign_posterior_haplotype_recombination_scores_per_imputation_signal(char* train_tag_haplocoded_genotype_matrix_fp,
	char* train_target_haplocoded_genotype_matrix_fp,
	char* train_sample_ids_list_fp,
	char* imp_tag_genocoded_genotype_matrix_fp,
	char* imp_target_geno_signal_matrix_fp,
	char* imp_sample_ids_list_fp,
	char* test_target_genocoded_genotype_matrix_fp,
	char* test_sample_ids_list_fp,
	char* recombination_map_dir,
	int n_vic_vars,
	double allelic_error_prob,
	char* probs_5col_fp)
{
	fprintf(stderr, "Assigning posterior haplotype recombination scores using %d vicinity variants using maps @ %s, allelic error probability of %.4f.\n", n_vic_vars, recombination_map_dir, allelic_error_prob);

	// Load the training data.
	fprintf(stderr, "Loading training tag genotypes.\n");
	vector<t_annot_region*>* train_tag_geno_regs = load_variant_signal_regions_wrapper(train_tag_haplocoded_genotype_matrix_fp, train_sample_ids_list_fp);

	for (int i_reg = 0; i_reg < train_tag_geno_regs->size(); i_reg++)
	{
		void** old_info = (void**)(train_tag_geno_regs->at(i_reg)->data);
		void** new_info = new void*[10];
		new_info[0] = old_info[0];
		train_tag_geno_regs->at(i_reg)->data = new_info;

		train_tag_geno_regs->at(i_reg)->score = 0;
	} // i_reg loop.
	fprintf(stderr, "Loaded %d training tag variants.\n", train_tag_geno_regs->size());

	// Load the training target genotypes, these represent the haplotype information.
	fprintf(stderr, "Loading training target genotypes.\n");
	vector<t_annot_region*>* train_target_geno_regs = load_variant_signal_regions_wrapper(train_target_haplocoded_genotype_matrix_fp, train_sample_ids_list_fp);
	for (int i_reg = 0; i_reg < train_target_geno_regs->size(); i_reg++)
	{
		void** old_info = (void**)(train_target_geno_regs->at(i_reg)->data);
		void** new_info = new void*[10];
		new_info[0] = old_info[0];

		train_target_geno_regs->at(i_reg)->score = 1;
	} // i_reg loop.
	fprintf(stderr, "Loaded %d training target variants.\n", train_target_geno_regs->size());

	// Load the test data: This is the known tag genotypes: Must be the same data as the imputed data.
	fprintf(stderr, "Loading testing (known) tag genotypes.\n");
	vector<t_annot_region*>* test_tag_geno_regs = load_variant_signal_regions_wrapper(imp_tag_genocoded_genotype_matrix_fp, test_sample_ids_list_fp);
	for (int i_reg = 0; i_reg < test_tag_geno_regs->size(); i_reg++)
	{
		test_tag_geno_regs->at(i_reg)->score = 0;
	} // i_reg loop.
	fprintf(stderr, "Loaded %d testing (known) tag variants.\n", test_tag_geno_regs->size());

	fprintf(stderr, "Loading testing (known) target genotypes.\n");
	vector<t_annot_region*>* test_target_geno_regs = load_variant_signal_regions_wrapper(test_target_genocoded_genotype_matrix_fp, test_sample_ids_list_fp);
	for (int i_reg = 0; i_reg < test_target_geno_regs->size(); i_reg++)
	{
		test_target_geno_regs->at(i_reg)->score = 1;
	} // i_reg loop.
	fprintf(stderr, "Loaded %d testing (known) target variants.\n", test_target_geno_regs->size());

	// Load the imputed data.
	fprintf(stderr, "Loading imputed tag variants.\n");
	vector<t_annot_region*>* imp_tag_geno_regs = load_variant_signal_regions_wrapper(imp_tag_genocoded_genotype_matrix_fp, imp_sample_ids_list_fp);
	for (int i_reg = 0; i_reg < imp_tag_geno_regs->size(); i_reg++)
	{
		imp_tag_geno_regs->at(i_reg)->score = 0;
	} // i_reg loop.
	fprintf(stderr, "Loaded %d imputed tag variants.\n", imp_tag_geno_regs->size());

	fprintf(stderr, "Loading imputed target variant (soft) signal regions.\n");
	int n_loaded_samples = 0;
	vector<t_annot_region*>* imp_target_geno_sig_regs = load_signal_regs_BED(imp_target_geno_signal_matrix_fp, n_loaded_samples);
	for (int i_reg = 0; i_reg < imp_target_geno_sig_regs->size(); i_reg++)
	{
		imp_target_geno_sig_regs->at(i_reg)->score = 1;
	} // i_reg loop.
	fprintf(stderr, "Loaded %d imputed target variant (soft) signal regions.\n", imp_target_geno_sig_regs->size());

	// Assign the test data to train data.
	// 0 : Training samples.
	// 1 : Testing samples: Known sample data for accuracy stats.
	// 2 : Imputed samples.
	vector<t_annot_region*>* target_intersects = intersect_regions_per_names(train_target_geno_regs, test_target_geno_regs, true);
	fprintf(stderr, "Processing %d intersects between target genotype regions of training and test (known) data.\n", target_intersects->size());
	for (int i_int = 0; i_int < target_intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(target_intersects->at(i_int)->data);
		t_annot_region* train_target = int_info->src_reg;
		t_annot_region* test_target = int_info->dest_reg;

		void** train_reg_info = (void**)(train_target->data);
		void** test_reg_info = (void**)(test_target->data);
		train_reg_info[1] = test_reg_info[0];
	} // i_int loop.

	vector<t_annot_region*>* tag_intersects = intersect_regions_per_names(train_tag_geno_regs, test_tag_geno_regs, true);
	fprintf(stderr, "Processing %d intersects between tag genotype regions of training and test (known) data.\n", tag_intersects->size());
	for (int i_int = 0; i_int < tag_intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(tag_intersects->at(i_int)->data);
		t_annot_region* train_tag = int_info->src_reg;
		t_annot_region* test_tag = int_info->dest_reg;

		void** train_reg_info = (void**)(train_tag->data);
		void** test_reg_info = (void**)(test_tag->data);
		train_reg_info[1] = test_reg_info[0];
	} // i_int loop.

	  // Assign the genotype signal data to train data.
	vector<t_annot_region*>* target_signal_intersects = intersect_regions_per_names(train_target_geno_regs, imp_target_geno_sig_regs, true);
	fprintf(stderr, "Processing %d intersects between target genotype regions of training and imputed signal data.\n", target_signal_intersects->size());
	for (int i_int = 0; i_int < target_signal_intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(target_signal_intersects->at(i_int)->data);
		t_annot_region* train_target = int_info->src_reg;
		t_annot_region* imp_signal_target = int_info->dest_reg;

		void** train_reg_info = (void**)(train_target->data);
		train_reg_info[2] = imp_signal_target->data;
	} // i_int loop.

	vector<t_annot_region*>* imp_tag_intersects = intersect_regions_per_names(train_tag_geno_regs, imp_tag_geno_regs, true);
	fprintf(stderr, "Processing %d intersects between tag genotype regions of training and imputed data.\n", imp_tag_intersects->size());
	for (int i_int = 0; i_int < imp_tag_intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(imp_tag_intersects->at(i_int)->data);
		t_annot_region* train_tag = int_info->src_reg;
		t_annot_region* test_tag = int_info->dest_reg;

		void** train_reg_info = (void**)(train_tag->data);
		void** test_reg_info = (void**)(test_tag->data);
		train_reg_info[2] = test_reg_info[0];
	} // i_int loop.

	vector<char*>* train_sample_ids = buffer_file(train_sample_ids_list_fp);
	vector<char*>* test_sample_ids = buffer_file(test_sample_ids_list_fp);
	vector<char*>* imp_sample_ids = buffer_file(imp_sample_ids_list_fp);

	if (test_sample_ids->size() != imp_sample_ids->size())
	{
		for (int i_s = 0; i_s < test_sample_ids->size(); i_s++)
		{
			if (!t_string::compare_strings(test_sample_ids->at(i_s), imp_sample_ids->at(i_s)))
			{
				fprintf(stderr, "The sample id's of the testing (known) samples are not matching imputed matrix sample id's.\n");
				exit(0);
			}
		} // i_s loop.
	}

	fprintf(stderr, "-----------------------------------\n");
	fprintf(stderr, "Attempting to verify coding of genotypes:\n");
	bool found_haplocoded_genotype = false;
	is_haplo_coded(train_tag_geno_regs, train_sample_ids->size(), found_haplocoded_genotype);
	if (!found_haplocoded_genotype)
	{
		fprintf(stderr, "***Training tag genotypes may not be haplocoded.\n ***");
	}
	is_haplo_coded(train_target_geno_regs, train_sample_ids->size(), found_haplocoded_genotype);
	if (!found_haplocoded_genotype)
	{
		fprintf(stderr, "***Training target genotypes may not be haplocoded.\n ***");
	}

	is_haplo_coded(test_tag_geno_regs, test_sample_ids->size(), found_haplocoded_genotype);
	if (found_haplocoded_genotype)
	{
		fprintf(stderr, "Testing tag genotypes are haplocoded.\n");
		exit(0);
	}

	is_haplo_coded(test_target_geno_regs, test_sample_ids->size(), found_haplocoded_genotype);
	if (found_haplocoded_genotype)
	{
		fprintf(stderr, "Testing target genotypes are haplocoded.\n");
		exit(0);
	}

	is_haplo_coded(imp_tag_geno_regs, imp_sample_ids->size(), found_haplocoded_genotype);
	if (found_haplocoded_genotype)
	{
		fprintf(stderr, "Imputed tag genotypes are haplocoded.\n");
		exit(0);
	}

	fprintf(stderr, "-----------------------------------\n");

	t_restr_annot_region_list* restr_train_target_var_regs = restructure_annot_regions(train_target_geno_regs);
	t_restr_annot_region_list* restr_train_tag_var_regs = restructure_annot_regions(train_tag_geno_regs);

	// All the non-ref variants.
	double n_nr_tot = 0;
	double n_tot = 0;

	double ML_n_nr_corr = 0;
	double ML_n_nr_tot = 0;

	double ML_n_corr = 0;
	double ML_n_tot = 0;

	double MAP_n_nr_corr = 0;
	double MAP_n_nr_tot = 0;

	double MAP_n_corr = 0;
	double MAP_n_tot = 0;

	double** ML_confusion_matrix = new double*[3];
	double** MAP_confusion_matrix = new double*[3];
	for (int geno = 0; geno < 3; geno++)
	{
		ML_confusion_matrix[geno] = new double[3];
		MAP_confusion_matrix[geno] = new double[3];
		memset(ML_confusion_matrix[geno], 0, sizeof(double) * 3);
		memset(MAP_confusion_matrix[geno], 0, sizeof(double) * 3);
	} // geno loop.

	FILE* f_probs_5col = open_f(probs_5col_fp, "w");

	FILE* f_per_variant_stats = open_f("per_variant_stats.txt", "w");
	for (int i_chr = 0; i_chr < restr_train_target_var_regs->chr_ids->size(); i_chr++)
	{
		vector<t_annot_region*>* cur_chr_target_var_regs = restr_train_target_var_regs->regions_per_chrom[i_chr];

		char cur_chr_recombination_rate_fp[1000];
		sprintf(cur_chr_recombination_rate_fp, "%s/%s.map", recombination_map_dir, restr_train_target_var_regs->chr_ids->at(i_chr));
		vector<t_annot_region*>* cur_chrom_recomb_regs = load_recombination_rates(cur_chr_recombination_rate_fp);
		if (cur_chrom_recomb_regs == NULL)
		{
			fprintf(stderr, "Could not load recombination rates from %s\n", cur_chr_recombination_rate_fp);
			exit(0);
		}

		int tag_var_chr_i = t_string::get_i_str(restr_train_tag_var_regs->chr_ids, restr_train_target_var_regs->chr_ids->at(i_chr));
		vector<t_annot_region*>* cur_chr_tag_var_regs = restr_train_tag_var_regs->regions_per_chrom[tag_var_chr_i];

		// Assign the recomb rates.
		fprintf(stderr, "Setting recombination rates using %d recombination mapped regions.\n", cur_chrom_recomb_regs->size());
		for (int i_reg = 0; i_reg < cur_chr_tag_var_regs->size(); i_reg++)
		{
			double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant(cur_chr_tag_var_regs->at(i_reg), cur_chrom_recomb_regs);
			cur_chr_tag_var_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
		} // i_reg loop.

		for (int i_reg = 0; i_reg < cur_chr_target_var_regs->size(); i_reg++)
		{
			double cur_reg_recomb_rate = get_cumulative_recomb_rate_per_variant(cur_chr_target_var_regs->at(i_reg), cur_chrom_recomb_regs);
			cur_chr_target_var_regs->at(i_reg)->dbl_score = cur_reg_recomb_rate;
		} // i_reg loop.

		for (int target_var_i = 0; target_var_i < cur_chr_target_var_regs->size(); target_var_i++)
		{
			if (target_var_i > 0 &&
				target_var_i % 1 == 0)
			{
				fprintf(stderr, "Processing %d. variant:\n\
ALL: ML sens: %.4f, ML PPV: %.4f\n\
ALL: MAP sens: %.4f, MAP PPV: %.4f\n\
NR: ML sens: %.4f, ML PPV: %.4f\n\
NR: MAP sens: %.4f, MAP PPV: %.4f\n",
target_var_i,
ML_n_corr / n_tot,
ML_n_corr / ML_n_tot,
MAP_n_corr / n_tot,
MAP_n_corr / MAP_n_tot,
ML_n_nr_corr / n_nr_tot,
ML_n_nr_corr / ML_n_nr_tot,
MAP_n_nr_corr / n_nr_tot,
MAP_n_nr_corr / MAP_n_nr_tot);
			}

			t_string_tokens* toks = t_string::tokenize_by_chars(cur_chr_target_var_regs->at(target_var_i)->name, "_");
			if (atof(toks->at(3)->str()) > 0.05)
			{
				fprintf(stderr, "Skipping: %s\n", cur_chr_target_var_regs->at(target_var_i)->name);
				continue;
			}
			t_string::clean_tokens(toks);

			// Find the tag variants around the target variant.
			int tag_j_var = 0;
			for (int j_var = 0;
				j_var < cur_chr_tag_var_regs->size();
				j_var++)
			{
				if (cur_chr_tag_var_regs->at(j_var)->start > cur_chr_target_var_regs->at(target_var_i)->start)
				{
					tag_j_var = j_var;
					break;
				}
			} // j_var loop.

			  // Make sure we have the right number of tag variants.
			if (tag_j_var < n_vic_vars ||
				(tag_j_var + n_vic_vars) >= cur_chr_tag_var_regs->size())
			{
				continue;
			}

			int n_total_haplotypes = train_sample_ids->size() * 2;

			double** per_haplo_per_var_alleles = new double*[n_total_haplotypes + 2];
			for (int hap_i = 0; hap_i < n_total_haplotypes; hap_i++)
			{
				per_haplo_per_var_alleles[hap_i] = new double[2 * n_vic_vars + 5];

				// Reset all entries to -1.
				for (int i_var = 0;
					i_var < (2 * n_vic_vars + 5);
					i_var++)
				{
					per_haplo_per_var_alleles[hap_i][i_var] = -1;
				} // i_var loop.				
			} // hap_i loop.

			vector<t_annot_region*>* pooled_tag_target_vars = new vector<t_annot_region*>();
			pooled_tag_target_vars->push_back(cur_chr_target_var_regs->at(target_var_i));
			for (int j_var = tag_j_var - n_vic_vars;
				j_var < tag_j_var + n_vic_vars;
				j_var++)
			{
				pooled_tag_target_vars->push_back(cur_chr_tag_var_regs->at(j_var));
			} // j_var loop.
			sort(pooled_tag_target_vars->begin(), pooled_tag_target_vars->end(), sort_regions);

			if (!t_string::compare_strings(pooled_tag_target_vars->at(n_vic_vars)->name, cur_chr_target_var_regs->at(target_var_i)->name))
			{
				fprintf(stderr, "Could not match target variant names: %s, %s\n\n", pooled_tag_target_vars->at(n_vic_vars)->name, cur_chr_target_var_regs->at(target_var_i)->name);
				exit(0);
			}

			// Setup the tag matrix for the current vicinity.
			for (int j_var = 0;
				j_var < pooled_tag_target_vars->size();
				j_var++)
			{
				void** cur_var_info = (void**)(pooled_tag_target_vars->at(j_var)->data);
				char* cur_train_sample_haplocoded_geno = (char*)(cur_var_info[0]);
				for (int i_s = 0; i_s < train_sample_ids->size(); i_s++)
				{
					int rel_var_i = j_var;
					per_haplo_per_var_alleles[i_s * 2][rel_var_i] = get_allele_per_haplotype(cur_train_sample_haplocoded_geno[i_s], 0);
					per_haplo_per_var_alleles[i_s * 2 + 1][rel_var_i] = get_allele_per_haplotype(cur_train_sample_haplocoded_geno[i_s], 1);
				} // i_s loop.
			} // j_var loop.

			// Find the haplotypes with the rare alleles at the target.
			vector<double*>* all_haplotypes = new vector<double*>();
			for (int i_hap = 0; i_hap < n_total_haplotypes; i_hap++)
			{
				all_haplotypes->push_back(per_haplo_per_var_alleles[i_hap]);
			} // i_hap loop.			

			// Get the unique haplotypes.
			vector<int>* n_cnt_per_uniq_haplotypes = new vector<int>();
			vector<double*>* all_unique_haplotypes = get_unique_haplotypes(all_haplotypes, n_cnt_per_uniq_haplotypes);

			// Sanity check on the counts from the unique haplotypes.
			int n_total_haps = 0;
			for (int i_uniq_hap = 0; i_uniq_hap < all_unique_haplotypes->size(); i_uniq_hap++)
			{
				n_total_haps += n_cnt_per_uniq_haplotypes->at(i_uniq_hap);
			} // i_uniq_hap loop.

			if (n_total_haps != all_haplotypes->size())
			{
				fprintf(stderr, "Could not match the total haplotypes per uniq counts to the number of all haplotypes: %d, %d\n",
					n_total_haps, all_haplotypes->size());

				exit(0);
			}

			// These are the accuracy stats.
			double cur_var_n_nr_tot = 0;
			double cur_var_ML_n_nr_corr = 0;
			double cur_var_ML_n_nr_tot = 0;

			double cur_var_MAP_n_nr_corr = 0;
			double cur_var_MAP_n_nr_tot = 0;

			double cur_var_n_tot = 0;
			double cur_var_ML_n_corr = 0;
			double cur_var_ML_n_tot = 0;

			double cur_var_MAP_n_corr = 0;
			double cur_var_MAP_n_tot = 0;

			// Following is a very expensive loop.
			for (int i_s = 0; i_s < test_sample_ids->size(); i_s++)
			{
				double per_geno_scores[3];
				per_geno_scores[0] = 0;
				per_geno_scores[1] = 0;
				per_geno_scores[2] = 0;

				double per_geno_cnts[3];
				per_geno_cnts[0] = 1;
				per_geno_cnts[1] = 1;
				per_geno_cnts[2] = 1;

				// Do "crossings" of all haplotypes and find the high scoring ones with the Li-Stephens HM probability.
				for (int recomb_posn = 0; recomb_posn < pooled_tag_target_vars->size(); recomb_posn++)
				{
					// Hapi and hapj recombine at recomb_posn'th position.
					for (int hap_i = 0; hap_i < all_unique_haplotypes->size(); hap_i++)
					{
						for (int hap_j = 0; hap_j < all_unique_haplotypes->size(); hap_j++)
						{
							// Get the score for this 

							// Sum the current haplotypes.
							double allelic_dist = 0;
							double tot_n_alleles = 1;
							int n_vars_processed = 0;
							for (int var_i = 0; var_i < pooled_tag_target_vars->size(); var_i++)
							{
								// Compute distances using tags only.
								if (pooled_tag_target_vars->at(var_i)->score == 0)
								{
									void** cur_var_info = (void**)(pooled_tag_target_vars->at(var_i)->data);
									char* cur_imp_var_geno = (char*)(cur_var_info[2]);

									int cur_tot_alleles = (all_unique_haplotypes->at(hap_i))[var_i] + (all_unique_haplotypes->at(hap_j))[var_i];

									if (cur_tot_alleles > 2)
									{
										fprintf(stderr, "Sanity check failed; the total number of alleles is greater than 2.\n");
										exit(0);
									}

									// This is the distance between the vicinity for the current haplotypes.
									allelic_dist += fabs(cur_imp_var_geno[i_s] - cur_tot_alleles);
									tot_n_alleles += cur_imp_var_geno[i_s];

									n_vars_processed++;
								} // tag check.
							} // var_i loop.

							if (n_vars_processed != pooled_tag_target_vars->size() - 1)
							{
								fprintf(stderr, "Sanity check failed: Could not match the number of tags to expected: %d, %d\n",
									n_vars_processed, pooled_tag_target_vars->size() - 1);
								exit(0);
							}

							// Compute the multiplicity of this haplotype combination.
							int cur_ij_multiplicity = n_cnt_per_uniq_haplotypes->at(hap_i) * n_cnt_per_uniq_haplotypes->at(hap_j);

							if (hap_i == hap_j)
							{
								cur_ij_multiplicity = (n_cnt_per_uniq_haplotypes->at(hap_i) - 1) * n_cnt_per_uniq_haplotypes->at(hap_j) / 2 + n_cnt_per_uniq_haplotypes->at(hap_j);
							}

							// Update the target genotype.
							int target_geno_per_hapij = (all_unique_haplotypes->at(hap_i))[n_vic_vars] + (all_unique_haplotypes->at(hap_j))[n_vic_vars];
							per_geno_scores[target_geno_per_hapij] += cur_ij_multiplicity * exp(-1 * exponential_allelic_prior_weight * (allelic_dist / tot_n_alleles));
							per_geno_cnts[target_geno_per_hapij] += cur_ij_multiplicity;
						} // hap_j loop.
					} // hap_i loop.
				} // recomb_posn loop.

				  // Normalize the scores.
				per_geno_scores[0] /= per_geno_cnts[0];
				per_geno_scores[1] /= per_geno_cnts[1];
				per_geno_scores[2] /= per_geno_cnts[2];

				// 
				double total_n_pairwise_haplo_combinations = per_geno_cnts[0] + per_geno_cnts[1] + per_geno_cnts[2];
				double n_expected_combinations = all_haplotypes->size() * (all_haplotypes->size() - 1) / 2 + all_haplotypes->size();

				if (total_n_pairwise_haplo_combinations != n_expected_combinations + 3)
				{
					fprintf(stderr, "Expected # of haplotype combinations do not hold: %.0f vs %.0f\n", total_n_pairwise_haplo_combinations, n_expected_combinations);
				}

				//per_geno_scores[0] = ((double)(rand() % 1000)) / 1000;
				//per_geno_scores[1] = ((double)(rand() % 1000)) / 1000;
				//per_geno_scores[2] = ((double)(rand() % 1000)) / 1000;

				// Get the imputed target signal and known target signal.
				void** cur_target_info = (void**)(cur_chr_target_var_regs->at(target_var_i)->data);

				char* test_target_signal_val = (char*)(cur_target_info[1]);
				double* imp_target_signal_val = (double*)(cur_target_info[2]);

				int testing_sample_geno = test_target_signal_val[i_s];
				if (testing_sample_geno > 2)
				{
					fprintf(stderr, "Sanity check failed: The known genotype is haplocoded.\n");
					exit(0);
				}

				// Update the signal.
				double tot_weighted_probs = 0;
				double ML_geno_prob = -1;
				double ML_geno = -1;
				double per_geno_LH[3];
				for (int geno = 0; geno < 3; geno++)
				{
					double cur_LH = exp(-1 * fabs(2 * (imp_target_signal_val[i_s] / 100) - geno));
					tot_weighted_probs += per_geno_scores[geno] * cur_LH;
					per_geno_LH[geno] = cur_LH;

					if (cur_LH > ML_geno_prob)
					{
						ML_geno_prob = cur_LH;
						ML_geno = geno;
					}
				} // geno loop.

				  // Compute the probabilities.
				double MAP_geno = -1;
				double MAP_geno_prob = -1;
				double post_probs[3];
				for (int geno = 0; geno < 3; geno++)
				{
					double cur_LH = per_geno_LH[geno];
					post_probs[geno] = per_geno_scores[geno] * cur_LH / tot_weighted_probs;

					if (MAP_geno_prob < post_probs[geno])
					{
						MAP_geno_prob = post_probs[geno];
						MAP_geno = geno;
					}
				} // geno loop.

				if (MAP_geno == -1 || ML_geno == -1)
				{
					fprintf(stderr, "Sanity check failed: Could not set MAP or ML genotype.\n");
					exit(0);
				}

				// Update the confusion matrix.
				MAP_confusion_matrix[(int)MAP_geno][testing_sample_geno]++;
				ML_confusion_matrix[(int)ML_geno][testing_sample_geno]++;

				fprintf(f_probs_5col, "%d,%d,%.4f,%.4f,%.4f\n",
					i_s, cur_chr_target_var_regs->at(target_var_i)->start,
					post_probs[0],
					post_probs[1],
					post_probs[2]);

				if (__DUMP_INPUTATION_UTILS_MSGS__)
				{
					if (testing_sample_geno > 0)
					{
						fprintf(stderr, "i_s [[%d]]: MAP_geno: %.0f/%d (%.4f)\n", i_s, MAP_geno, testing_sample_geno, MAP_geno_prob);
					}
				}

				//////////////////////////////////////////////////////////////////////////////////////////
				// Update all statistics.
				if (MAP_geno == testing_sample_geno)
				{
					MAP_n_corr++;
					cur_var_MAP_n_corr++;
				}

				cur_var_MAP_n_tot++;
				MAP_n_tot++;

				if (ML_geno == testing_sample_geno)
				{
					ML_n_corr++;
					cur_var_ML_n_corr++;
				}

				cur_var_ML_n_tot++;
				ML_n_tot++;

				n_tot++;
				cur_var_n_tot++;

				//////////////////////////////////////////////////////////////////////////////////////////
				// Update non-ref statistics.
				if (MAP_geno > 0)
				{
					if (MAP_geno == testing_sample_geno)
					{
						MAP_n_nr_corr++;
						cur_var_MAP_n_nr_corr++;
					}

					MAP_n_nr_tot++;
					cur_var_MAP_n_nr_tot++;
				}

				if (ML_geno > 0)
				{
					if (ML_geno == testing_sample_geno)
					{
						ML_n_nr_corr++;
						cur_var_ML_n_nr_corr++;
					}

					ML_n_nr_tot++;
					cur_var_ML_n_nr_tot++;
				}

				if (testing_sample_geno > 0)
				{
					n_nr_tot++;
					cur_var_n_nr_tot++;
				}
				//////////////////////////////////////////////////////////////////////////////////////////
			} // i_s loop.

			fprintf(stderr, "%s\t%d\t%s\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%d\n",
				cur_chr_target_var_regs->at(target_var_i)->chrom, cur_chr_target_var_regs->at(target_var_i)->start, cur_chr_target_var_regs->at(target_var_i)->name,
				cur_var_ML_n_corr, cur_var_ML_n_tot, cur_var_MAP_n_corr, cur_var_MAP_n_tot, cur_var_n_tot,
				cur_var_ML_n_nr_corr, cur_var_ML_n_nr_tot, cur_var_MAP_n_nr_corr, cur_var_MAP_n_nr_tot, cur_var_n_nr_tot,
				all_unique_haplotypes->size());

			fprintf(f_per_variant_stats, "%s\t%d\t%s\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%d\n",
				cur_chr_target_var_regs->at(target_var_i)->chrom, cur_chr_target_var_regs->at(target_var_i)->start, cur_chr_target_var_regs->at(target_var_i)->name,
				cur_var_ML_n_corr, cur_var_ML_n_tot, cur_var_MAP_n_corr, cur_var_MAP_n_tot, cur_var_n_tot,
				cur_var_ML_n_nr_corr, cur_var_ML_n_nr_tot, cur_var_MAP_n_nr_corr, cur_var_MAP_n_nr_tot, cur_var_n_nr_tot,
				all_unique_haplotypes->size());

			// Dump the confusion matrix.
			fprintf(stderr, "MAP Confusion Matrix:\n");
			for (int pred_geno = 0; pred_geno < 3; pred_geno++)
			{
				fprintf(stderr, "Imp %d:", pred_geno);
				for (int known_geno = 0; known_geno < 3; known_geno++)
				{
					fprintf(stderr, "\t%.0f", MAP_confusion_matrix[pred_geno][known_geno]);
				} // known_geno loop.

				fprintf(stderr, "\n");
			} // pred_geno loop.

			fprintf(stderr, "ML Confusion Matrix:\n");
			for (int pred_geno = 0; pred_geno < 3; pred_geno++)
			{
				fprintf(stderr, "Imp %d:", pred_geno);
				for (int known_geno = 0; known_geno < 3; known_geno++)
				{
					fprintf(stderr, "\t%.0f", ML_confusion_matrix[pred_geno][known_geno]);
				} // known_geno loop.

				fprintf(stderr, "\n");
			} // pred_geno loop.

			  //delete lowMAF_containing_haplotypes;
			  //delete lowMAF_missing_haplotypes;
			delete all_haplotypes;
			//delete unique_lowMAF_containing_haplotypes;
			//delete unique_lowMAF_missing_haplotypes;

			for (int hap_i = 0; hap_i < n_total_haplotypes; hap_i++)
			{
				delete[] per_haplo_per_var_alleles[hap_i];
			} // hap_i loop.
			delete[] per_haplo_per_var_alleles;

			//// Set the per variant testing genotypes.
			//for (int i_s = 0; i_s < test_sample_ids->size(); i_s++)
			//{
			//	delete[] per_testing_sample_per_var_testing_geno[i_s];
			//} // i_s loop
			//delete[] per_testing_sample_per_var_testing_geno;

			//for (int i_s = 0; i_s < imp_sample_ids->size(); i_s++)
			//{
			//	delete[] per_imp_sample_per_var_testing_geno[i_s];
			//}
			//delete[] per_imp_sample_per_var_testing_geno;
		} // i_var loop.
	} // i_chr loop.

	close_f(f_probs_5col, probs_5col_fp);
	fclose(f_per_variant_stats);

	FILE* f_MAP_confusion_matrix = open_f("MAP_confusion_matrices.txt", "w");
	for (int MAP_geno = 0; MAP_geno < 3; MAP_geno++)
	{
		for (int known_geno = 0; known_geno < 3; known_geno++)
		{
			fprintf(f_MAP_confusion_matrix, "%d\t%d\t%.0f\n", MAP_geno, known_geno, MAP_confusion_matrix[MAP_geno][known_geno]);
		} // known geno loop.
	} // map geno loop. 
	fclose(f_MAP_confusion_matrix);

	FILE* f_ML_confusion_matrix = open_f("ML_confusion_matrices.txt", "w");
	for (int ML_geno = 0; ML_geno < 3; ML_geno++)
	{
		for (int known_geno = 0; known_geno < 3; known_geno++)
		{
			fprintf(f_ML_confusion_matrix, "%d\t%d\t%.0f\n", ML_geno, known_geno, ML_confusion_matrix[ML_geno][known_geno]);
		} // known geno loop.
	} // map geno loop. 
	fclose(f_ML_confusion_matrix);

	//// Re-save the imputed variants.
	//char output_genotype_matrix_bed_fp[1000];
	//sprintf(output_genotype_matrix_bed_fp, "consistency_calibrated_genotypes.txt");
	//dump_geno_sig_regs_plain(imp_target_geno_regs, imp_sample_ids, false, output_genotype_matrix_bed_fp);

	//// Re-save the imputed variants.
	//char output_genotype_signals_bed_fp[1000];
	//sprintf(output_genotype_signals_bed_fp, "consistency_calibrated_genotypes.sig");
	//FILE* f_calib_geno_sigs = open_f(output_genotype_signals_bed_fp, "w");
	//for (int i_reg = 0; i_reg < imp_target_geno_sig_regs->size(); i_reg++)
	//{
	//	double* cur_reg_signals = (double*)(imp_target_geno_sig_regs->at(i_reg)->data);
	//	fprintf(f_calib_geno_sigs, "%s\t%d\t%d\t%s", imp_target_geno_sig_regs->at(i_reg)->chrom,
	//		translate_coord(imp_target_geno_sig_regs->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
	//		translate_coord(imp_target_geno_sig_regs->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
	//		imp_target_geno_sig_regs->at(i_reg)->name);

	//	for (int i_s = 0; i_s < imp_sample_ids->size(); i_s++)
	//	{
	//		fprintf(f_calib_geno_sigs, "\t%.0f", cur_reg_signals[i_s]);
	//	}

	//	fprintf(f_calib_geno_sigs, "\n");
	//} // i_reg loop.
	//close_f(f_calib_geno_sigs, output_genotype_signals_bed_fp);
}
*/

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
	char* probs_5col_fp)
{
	fprintf(stderr, "Assigning posterior haplotype consistency scores using %d vicinity variants and %.4f allelic weight.\n", n_vic_vars, exponential_allelic_prior_weight);

	// Load the training data.
	fprintf(stderr, "Loading training tag genotypes.\n");
	vector<t_annot_region*>* train_tag_geno_regs = load_variant_signal_regions_wrapper(train_tag_haplocoded_genotype_matrix_fp, train_sample_ids_list_fp);

	for (int i_reg = 0; i_reg < train_tag_geno_regs->size(); i_reg++)
	{
		void** old_info = (void**)(train_tag_geno_regs->at(i_reg)->data);
		void** new_info = new void*[10];
		new_info[0] = old_info[0];
		train_tag_geno_regs->at(i_reg)->data = new_info;

		train_tag_geno_regs->at(i_reg)->score = 0;
	} // i_reg loop.
	fprintf(stderr, "Loaded %d training tag variants.\n", train_tag_geno_regs->size());

	// Load the training target genotypes, these represent the haplotype information.
	fprintf(stderr, "Loading training target genotypes.\n");
	vector<t_annot_region*>* train_target_geno_regs = load_variant_signal_regions_wrapper(train_target_haplocoded_genotype_matrix_fp, train_sample_ids_list_fp);
	for (int i_reg = 0; i_reg < train_target_geno_regs->size(); i_reg++)
	{
		void** old_info = (void**)(train_target_geno_regs->at(i_reg)->data);
		void** new_info = new void*[10];
		new_info[0] = old_info[0];

		train_target_geno_regs->at(i_reg)->score = 1;
	} // i_reg loop.
	fprintf(stderr, "Loaded %d training target variants.\n", train_target_geno_regs->size());

	// Load the test data: This is the known tag genotypes: Must be the same data as the imputed data.
	fprintf(stderr, "Loading testing (known) tag genotypes.\n");
	vector<t_annot_region*>* test_tag_geno_regs = load_variant_signal_regions_wrapper(imp_tag_genocoded_genotype_matrix_fp, test_sample_ids_list_fp);
	for (int i_reg = 0; i_reg < test_tag_geno_regs->size(); i_reg++)
	{
		test_tag_geno_regs->at(i_reg)->score = 0;
	} // i_reg loop.
	fprintf(stderr, "Loaded %d testing (known) tag variants.\n", test_tag_geno_regs->size());

	fprintf(stderr, "Loading testing (known) target genotypes.\n");
	vector<t_annot_region*>* test_target_geno_regs = load_variant_signal_regions_wrapper(test_target_genocoded_genotype_matrix_fp, test_sample_ids_list_fp);
	for (int i_reg = 0; i_reg < test_target_geno_regs->size(); i_reg++)
	{
		test_target_geno_regs->at(i_reg)->score = 1;
	} // i_reg loop.
	fprintf(stderr, "Loaded %d testing (known) target variants.\n", test_target_geno_regs->size());

	// Load the imputed data.
	fprintf(stderr, "Loading imputed tag variants.\n");
	vector<t_annot_region*>* imp_tag_geno_regs = load_variant_signal_regions_wrapper(imp_tag_genocoded_genotype_matrix_fp, imp_sample_ids_list_fp);
	for (int i_reg = 0; i_reg < imp_tag_geno_regs->size(); i_reg++)
	{
		imp_tag_geno_regs->at(i_reg)->score = 0;
	} // i_reg loop.
	fprintf(stderr, "Loaded %d imputed tag variants.\n", imp_tag_geno_regs->size());

	fprintf(stderr, "Loading imputed target variant (soft) signal regions.\n");
	int n_loaded_samples = 0;
	vector<t_annot_region*>* imp_target_geno_sig_regs = load_signal_regs_BED(imp_target_geno_signal_matrix_fp, n_loaded_samples);
	for (int i_reg = 0; i_reg < imp_target_geno_sig_regs->size(); i_reg++)
	{
		imp_target_geno_sig_regs->at(i_reg)->score = 1;
	} // i_reg loop.
	fprintf(stderr, "Loaded %d imputed target variant (soft) signal regions.\n", imp_target_geno_sig_regs->size());

	// Assign the test data to train data.
	// 0 : Training samples.
	// 1 : Testing samples: Known sample data for accuracy stats.
	// 2 : Imputed samples.
	vector<t_annot_region*>* target_intersects = intersect_regions_per_names(train_target_geno_regs, test_target_geno_regs, true);
	fprintf(stderr, "Processing %d intersects between target genotype regions of training and test (known) data.\n", target_intersects->size());
	for (int i_int = 0; i_int < target_intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(target_intersects->at(i_int)->data);
		t_annot_region* train_target = int_info->src_reg;
		t_annot_region* test_target = int_info->dest_reg;

		void** train_reg_info = (void**)(train_target->data);
		void** test_reg_info = (void**)(test_target->data);
		train_reg_info[1] = test_reg_info[0];
	} // i_int loop.

	vector<t_annot_region*>* tag_intersects = intersect_regions_per_names(train_tag_geno_regs, test_tag_geno_regs, true);
	fprintf(stderr, "Processing %d intersects between tag genotype regions of training and test (known) data.\n", tag_intersects->size());
	for (int i_int = 0; i_int < tag_intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(tag_intersects->at(i_int)->data);
		t_annot_region* train_tag = int_info->src_reg;
		t_annot_region* test_tag = int_info->dest_reg;

		void** train_reg_info = (void**)(train_tag->data);
		void** test_reg_info = (void**)(test_tag->data);
		train_reg_info[1] = test_reg_info[0];
	} // i_int loop.

	// Assign the genotype signal data to train data.
	vector<t_annot_region*>* target_signal_intersects = intersect_regions_per_names(train_target_geno_regs, imp_target_geno_sig_regs, true);
	fprintf(stderr, "Processing %d intersects between target genotype regions of training and imputed signal data.\n", target_signal_intersects->size());
	for (int i_int = 0; i_int < target_signal_intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(target_signal_intersects->at(i_int)->data);
		t_annot_region* train_target = int_info->src_reg;
		t_annot_region* imp_signal_target = int_info->dest_reg;

		void** train_reg_info = (void**)(train_target->data);
		train_reg_info[2] = imp_signal_target->data;
	} // i_int loop.

	vector<t_annot_region*>* imp_tag_intersects = intersect_regions_per_names(train_tag_geno_regs, imp_tag_geno_regs, true);
	fprintf(stderr, "Processing %d intersects between tag genotype regions of training and imputed data.\n", imp_tag_intersects->size());
	for (int i_int = 0; i_int < imp_tag_intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(imp_tag_intersects->at(i_int)->data);
		t_annot_region* train_tag = int_info->src_reg;
		t_annot_region* test_tag = int_info->dest_reg;

		void** train_reg_info = (void**)(train_tag->data);
		void** test_reg_info = (void**)(test_tag->data);
		train_reg_info[2] = test_reg_info[0];
	} // i_int loop.

	vector<char*>* train_sample_ids = buffer_file(train_sample_ids_list_fp);
	vector<char*>* test_sample_ids = buffer_file(test_sample_ids_list_fp);
	vector<char*>* imp_sample_ids = buffer_file(imp_sample_ids_list_fp);

	if (test_sample_ids->size() != imp_sample_ids->size())
	{
		for (int i_s = 0; i_s < test_sample_ids->size(); i_s++)
		{
			if (!t_string::compare_strings(test_sample_ids->at(i_s), imp_sample_ids->at(i_s)))
			{
				fprintf(stderr, "The sample id's of the testing (known) samples are not matching imputed matrix sample id's.\n");
				exit(0);
			}
		} // i_s loop.
	}

	fprintf(stderr, "-----------------------------------\n");
	fprintf(stderr, "Attempting to verify coding of genotypes:\n");
	bool found_haplocoded_genotype = false;
	is_haplo_coded(train_tag_geno_regs, train_sample_ids->size(), found_haplocoded_genotype);
	if (!found_haplocoded_genotype)
	{
		fprintf(stderr, "***Training tag genotypes may not be haplocoded.\n ***");
	}
	is_haplo_coded(train_target_geno_regs, train_sample_ids->size(), found_haplocoded_genotype);
	if (!found_haplocoded_genotype)
	{
		fprintf(stderr, "***Training target genotypes may not be haplocoded.\n ***");
	}

	is_haplo_coded(test_tag_geno_regs, test_sample_ids->size(), found_haplocoded_genotype);
	if (found_haplocoded_genotype)
	{
		fprintf(stderr, "Testing tag genotypes are haplocoded.\n");
		exit(0);
	}

	is_haplo_coded(test_target_geno_regs, test_sample_ids->size(), found_haplocoded_genotype);
	if (found_haplocoded_genotype)
	{
		fprintf(stderr, "Testing target genotypes are haplocoded.\n");
		exit(0);
	}

	is_haplo_coded(imp_tag_geno_regs, imp_sample_ids->size(), found_haplocoded_genotype);
	if (found_haplocoded_genotype)
	{
		fprintf(stderr, "Imputed tag genotypes are haplocoded.\n");
		exit(0);
	}

	fprintf(stderr, "-----------------------------------\n");

	t_restr_annot_region_list* restr_train_target_var_regs = restructure_annot_regions(train_target_geno_regs);
	t_restr_annot_region_list* restr_train_tag_var_regs = restructure_annot_regions(train_tag_geno_regs);

	// All the non-ref variants.
	double n_nr_tot = 0;
	double n_tot = 0;

	double ML_n_nr_corr = 0;
	double ML_n_nr_tot = 0;

	double ML_n_corr = 0;
	double ML_n_tot = 0;
	
	double MAP_n_nr_corr = 0;
	double MAP_n_nr_tot = 0;

	double MAP_n_corr = 0;
	double MAP_n_tot = 0;
	
	double** ML_confusion_matrix = new double*[3];
	double** MAP_confusion_matrix = new double*[3];
	for (int geno = 0; geno < 3; geno++)
	{
		ML_confusion_matrix[geno] = new double[3];
		MAP_confusion_matrix[geno] = new double[3];
		memset(ML_confusion_matrix[geno], 0, sizeof(double) * 3);
		memset(MAP_confusion_matrix[geno], 0, sizeof(double) * 3);
	} // geno loop.

	FILE* f_probs_5col = open_f(probs_5col_fp, "w");

	FILE* f_per_variant_stats = open_f("per_variant_stats.txt", "w");
	for (int i_chr = 0; i_chr < restr_train_target_var_regs->chr_ids->size(); i_chr++)
	{
		vector<t_annot_region*>* cur_chr_target_var_regs = restr_train_target_var_regs->regions_per_chrom[i_chr];

		int tag_var_chr_i = t_string::get_i_str(restr_train_tag_var_regs->chr_ids, restr_train_target_var_regs->chr_ids->at(i_chr));
		vector<t_annot_region*>* cur_chr_tag_var_regs = restr_train_tag_var_regs->regions_per_chrom[tag_var_chr_i];
		for (int target_var_i = 0; target_var_i < cur_chr_target_var_regs->size(); target_var_i++)
		{
			if (target_var_i > 0 &&
				target_var_i % 1 == 0)
			{
				fprintf(stderr, "Processing %d. variant:\n\
ALL: ML sens: %.4f, ML PPV: %.4f\n\
ALL: MAP sens: %.4f, MAP PPV: %.4f\n\
NR: ML sens: %.4f, ML PPV: %.4f\n\
NR: MAP sens: %.4f, MAP PPV: %.4f\n",
						target_var_i,
						ML_n_corr / n_tot,
						ML_n_corr / ML_n_tot,
						MAP_n_corr / n_tot,
						MAP_n_corr / MAP_n_tot,
						ML_n_nr_corr / n_nr_tot,
						ML_n_nr_corr / ML_n_nr_tot,
						MAP_n_nr_corr / n_nr_tot,
						MAP_n_nr_corr / MAP_n_nr_tot);
			}

			t_string_tokens* toks = t_string::tokenize_by_chars(cur_chr_target_var_regs->at(target_var_i)->name, "_");
			if (atof(toks->at(3)->str()) > 0.05)
			{
				fprintf(stderr, "Skipping: %s\n", cur_chr_target_var_regs->at(target_var_i)->name);
				continue;
			}
			t_string::clean_tokens(toks);

			// Find the tag variants around the target variant.
			int tag_j_var = 0;
			for (int j_var = 0;
				j_var < cur_chr_tag_var_regs->size();
				j_var++)
			{
				if (cur_chr_tag_var_regs->at(j_var)->start > cur_chr_target_var_regs->at(target_var_i)->start)
				{
					tag_j_var = j_var;
					break;
				}
			} // j_var loop.

			// Make sure we have the right number of tag variants.
			if (tag_j_var < n_vic_vars || 
				(tag_j_var + n_vic_vars) >= cur_chr_tag_var_regs->size())
			{
				continue;
			}

			int n_total_haplotypes = train_sample_ids->size() * 2;

			double** per_haplo_per_var_alleles = new double*[n_total_haplotypes + 2];
			for (int hap_i = 0; hap_i < n_total_haplotypes; hap_i++)
			{
				per_haplo_per_var_alleles[hap_i] = new double[2 * n_vic_vars + 5];

				// Reset all entries to -1.
				for (int i_var = 0;
					i_var < (2 * n_vic_vars + 5);
					i_var++)
				{
					per_haplo_per_var_alleles[hap_i][i_var] = -1;
				} // i_var loop.				
			} // hap_i loop.

			vector<t_annot_region*>* pooled_tag_target_vars = new vector<t_annot_region*>();
			pooled_tag_target_vars->push_back(cur_chr_target_var_regs->at(target_var_i));
			for (int j_var = tag_j_var - n_vic_vars;
				j_var < tag_j_var + n_vic_vars;
				j_var++)
			{
				pooled_tag_target_vars->push_back(cur_chr_tag_var_regs->at(j_var));
			} // j_var loop.
			sort(pooled_tag_target_vars->begin(), pooled_tag_target_vars->end(), sort_regions);

			if (!t_string::compare_strings(pooled_tag_target_vars->at(n_vic_vars)->name, cur_chr_target_var_regs->at(target_var_i)->name))
			{
				fprintf(stderr, "Could not match target variant names: %s, %s\n\n", pooled_tag_target_vars->at(n_vic_vars)->name, cur_chr_target_var_regs->at(target_var_i)->name);
				exit(0);
			}

			// Setup the tag matrix for the current vicinity.
			for (int j_var = 0;
				j_var < pooled_tag_target_vars->size();
				j_var++)
			{
				void** cur_var_info = (void**)(pooled_tag_target_vars->at(j_var)->data);
				char* cur_train_sample_haplocoded_geno = (char*)(cur_var_info[0]);
				for (int i_s = 0; i_s < train_sample_ids->size(); i_s++)
				{
					int rel_var_i = j_var;
					per_haplo_per_var_alleles[i_s * 2][rel_var_i] = get_allele_per_haplotype(cur_train_sample_haplocoded_geno[i_s], 0);
					per_haplo_per_var_alleles[i_s * 2 + 1][rel_var_i] = get_allele_per_haplotype(cur_train_sample_haplocoded_geno[i_s], 1);
				} // i_s loop.
			} // j_var loop.

			// Find the haplotypes with the rare alleles at the target.
			vector<double*>* all_haplotypes = new vector<double*>();
			for (int i_hap = 0; i_hap < n_total_haplotypes; i_hap++)
			{
				all_haplotypes->push_back(per_haplo_per_var_alleles[i_hap]);
			} // i_hap loop.			

			// Get the unique haplotypes.
			vector<int>* n_cnt_per_uniq_haplotypes = new vector<int>();
			vector<double*>* all_unique_haplotypes = get_unique_haplotypes(all_haplotypes, n_cnt_per_uniq_haplotypes);

			// Sanity check on the counts from the unique haplotypes.
			int n_total_haps = 0;
			for (int i_uniq_hap = 0; i_uniq_hap < all_unique_haplotypes->size(); i_uniq_hap++)
			{
				n_total_haps += n_cnt_per_uniq_haplotypes->at(i_uniq_hap);
			} // i_uniq_hap loop.

			if (n_total_haps != all_haplotypes->size())
			{
				fprintf(stderr, "Could not match the total haplotypes per uniq counts to the number of all haplotypes: %d, %d\n",
						n_total_haps, all_haplotypes->size());

				exit(0);
			}

			// These are the accuracy stats.
			double cur_var_n_nr_tot = 0;
			double cur_var_ML_n_nr_corr = 0;
			double cur_var_ML_n_nr_tot = 0;

			double cur_var_MAP_n_nr_corr = 0;
			double cur_var_MAP_n_nr_tot = 0;
	
			double cur_var_n_tot = 0;
			double cur_var_ML_n_corr = 0;
			double cur_var_ML_n_tot = 0;

			double cur_var_MAP_n_corr = 0;
			double cur_var_MAP_n_tot = 0;

			// Following is a very expensive loop.
			for (int i_s = 0; i_s < test_sample_ids->size(); i_s++)
			{
				double per_geno_scores[3];
				per_geno_scores[0] = 0;
				per_geno_scores[1] = 0;
				per_geno_scores[2] = 0;

				double per_geno_cnts[3];
				per_geno_cnts[0] = 1;
				per_geno_cnts[1] = 1;
				per_geno_cnts[2] = 1;

				for (int hap_i = 0; hap_i < all_unique_haplotypes->size(); hap_i++)
				{
					for (int hap_j = hap_i; hap_j < all_unique_haplotypes->size(); hap_j++)
					{
						// Sum the current haplotypes.
						double allelic_dist = 0;
						double tot_n_alleles = 1;
						int n_vars_processed = 0;
						for (int var_i = 0; var_i < pooled_tag_target_vars->size(); var_i++)
						{
							// Compute distances using tags only.
							if (pooled_tag_target_vars->at(var_i)->score == 0)
							{
								void** cur_var_info = (void**)(pooled_tag_target_vars->at(var_i)->data);
								char* cur_imp_var_geno = (char*)(cur_var_info[2]);

								int cur_tot_alleles = (all_unique_haplotypes->at(hap_i))[var_i] + (all_unique_haplotypes->at(hap_j))[var_i];

								if (cur_tot_alleles > 2)
								{
									fprintf(stderr, "Sanity check failed; the total number of alleles is greater than 2.\n");
									exit(0);
								}

								// This is the distance between the vicinity for the current haplotypes.
								allelic_dist += fabs(cur_imp_var_geno[i_s] - cur_tot_alleles);
								tot_n_alleles += cur_imp_var_geno[i_s];

								n_vars_processed++;
							} // tag check.
						} // var_i loop.

						if (n_vars_processed != pooled_tag_target_vars->size() - 1)
						{
							fprintf(stderr, "Sanity check failed: Could not match the number of tags to expected: %d, %d\n",
								n_vars_processed, pooled_tag_target_vars->size() - 1);
							exit(0);
						}

						// Compute the multiplicity of this haplotype combination.
						int cur_ij_multiplicity = n_cnt_per_uniq_haplotypes->at(hap_i) * n_cnt_per_uniq_haplotypes->at(hap_j);

						if (hap_i == hap_j)
						{
							cur_ij_multiplicity = (n_cnt_per_uniq_haplotypes->at(hap_i) - 1) * n_cnt_per_uniq_haplotypes->at(hap_j) / 2 + n_cnt_per_uniq_haplotypes->at(hap_j);
						}

						// Update the target genotype.
						int target_geno_per_hapij = (all_unique_haplotypes->at(hap_i))[n_vic_vars] + (all_unique_haplotypes->at(hap_j))[n_vic_vars];
						per_geno_scores[target_geno_per_hapij] += cur_ij_multiplicity * exp(-1 * exponential_allelic_prior_weight * (allelic_dist / tot_n_alleles));
						per_geno_cnts[target_geno_per_hapij] += cur_ij_multiplicity;
					} // hap_j loop.
				} // hap_i loop.

				// Normalize the scores.
				per_geno_scores[0] /= per_geno_cnts[0];
				per_geno_scores[1] /= per_geno_cnts[1];
				per_geno_scores[2] /= per_geno_cnts[2];

				// 
				double total_n_pairwise_haplo_combinations = per_geno_cnts[0] + per_geno_cnts[1] + per_geno_cnts[2];
				double n_expected_combinations = all_haplotypes->size() * (all_haplotypes->size() - 1) / 2 + all_haplotypes->size();

				if (total_n_pairwise_haplo_combinations != n_expected_combinations + 3)
				{
					fprintf(stderr, "Expected # of haplotype combinations do not hold: %.0f vs %.0f\n", total_n_pairwise_haplo_combinations, n_expected_combinations);
				}

				//per_geno_scores[0] = ((double)(rand() % 1000)) / 1000;
				//per_geno_scores[1] = ((double)(rand() % 1000)) / 1000;
				//per_geno_scores[2] = ((double)(rand() % 1000)) / 1000;

				// Get the imputed target signal and known target signal.
				void** cur_target_info = (void**)(cur_chr_target_var_regs->at(target_var_i)->data);

				char* test_target_signal_val = (char*)(cur_target_info[1]);
				double* imp_target_signal_val = (double*)(cur_target_info[2]);

				int testing_sample_geno = test_target_signal_val[i_s];
				if (testing_sample_geno > 2)
				{
					fprintf(stderr, "Sanity check failed: The known genotype is haplocoded.\n");
					exit(0);
				}

				// Update the signal.
				double tot_weighted_probs = 0;
				double ML_geno_prob = -1;
				double ML_geno = -1;
				double per_geno_LH[3];
				for(int geno = 0; geno < 3; geno++)
				{
					double cur_LH = exp(-1 * fabs(2 * (imp_target_signal_val[i_s] / 100) - geno));
					tot_weighted_probs += per_geno_scores[geno] * cur_LH;
					per_geno_LH[geno] = cur_LH;

					if (cur_LH > ML_geno_prob)
					{
						ML_geno_prob = cur_LH;
						ML_geno = geno;
					}
				} // geno loop.

				// Compute the probabilities.
				double MAP_geno = -1;
				double MAP_geno_prob = -1;
				double post_probs[3];
				for (int geno = 0; geno < 3; geno++)
				{
					double cur_LH = per_geno_LH[geno];
					post_probs[geno] = per_geno_scores[geno] * cur_LH / tot_weighted_probs;

					if (MAP_geno_prob < post_probs[geno])
					{
						MAP_geno_prob = post_probs[geno];
						MAP_geno = geno;
					}
				} // geno loop.

				if (MAP_geno == -1 || ML_geno == -1)
				{
					fprintf(stderr, "Sanity check failed: Could not set MAP or ML genotype.\n");
					exit(0);
				}

				// Update the confusion matrix.
				MAP_confusion_matrix[(int)MAP_geno][testing_sample_geno]++;
				ML_confusion_matrix[(int)ML_geno][testing_sample_geno]++;

				fprintf(f_probs_5col, "%d,%d,%.4f,%.4f,%.4f\n",
						i_s, cur_chr_target_var_regs->at(target_var_i)->start,
						post_probs[0],
						post_probs[1],
						post_probs[2]);

				if (__DUMP_INPUTATION_UTILS_MSGS__)
				{
					if (testing_sample_geno > 0)
					{
						fprintf(stderr, "i_s [[%d]]: MAP_geno: %.0f/%d (%.4f)\n", i_s, MAP_geno, testing_sample_geno, MAP_geno_prob);
					}
				}

				//////////////////////////////////////////////////////////////////////////////////////////
				// Update all statistics.
				if (MAP_geno == testing_sample_geno)
				{
					MAP_n_corr++;
					cur_var_MAP_n_corr++;
				}

				cur_var_MAP_n_tot++;
				MAP_n_tot++;

				if (ML_geno == testing_sample_geno)
				{
					ML_n_corr++;
					cur_var_ML_n_corr++;
				}

				cur_var_ML_n_tot++;
				ML_n_tot++;

				n_tot++;
				cur_var_n_tot++;

				//////////////////////////////////////////////////////////////////////////////////////////
				// Update non-ref statistics.
				if (MAP_geno > 0)
				{
					if (MAP_geno == testing_sample_geno)
					{
						MAP_n_nr_corr++;
						cur_var_MAP_n_nr_corr++;
					}

					MAP_n_nr_tot++;
					cur_var_MAP_n_nr_tot++;
				}

				if (ML_geno > 0)
				{
					if (ML_geno == testing_sample_geno)
					{
						ML_n_nr_corr++;
						cur_var_ML_n_nr_corr++;
					}

					ML_n_nr_tot++;
					cur_var_ML_n_nr_tot++;
				}

				if (testing_sample_geno > 0)
				{
					n_nr_tot++;
					cur_var_n_nr_tot++;
				}
				//////////////////////////////////////////////////////////////////////////////////////////
			} // i_s loop.

			fprintf(stderr, "%s\t%d\t%s\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%d\n",
				cur_chr_target_var_regs->at(target_var_i)->chrom, cur_chr_target_var_regs->at(target_var_i)->start, cur_chr_target_var_regs->at(target_var_i)->name,
				cur_var_ML_n_corr, cur_var_ML_n_tot, cur_var_MAP_n_corr, cur_var_MAP_n_tot, cur_var_n_tot,
				cur_var_ML_n_nr_corr, cur_var_ML_n_nr_tot, cur_var_MAP_n_nr_corr, cur_var_MAP_n_nr_tot, cur_var_n_nr_tot,
				all_unique_haplotypes->size());

			fprintf(f_per_variant_stats, "%s\t%d\t%s\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%d\n",
				cur_chr_target_var_regs->at(target_var_i)->chrom, cur_chr_target_var_regs->at(target_var_i)->start, cur_chr_target_var_regs->at(target_var_i)->name,
				cur_var_ML_n_corr, cur_var_ML_n_tot, cur_var_MAP_n_corr, cur_var_MAP_n_tot, cur_var_n_tot,
				cur_var_ML_n_nr_corr, cur_var_ML_n_nr_tot, cur_var_MAP_n_nr_corr, cur_var_MAP_n_nr_tot, cur_var_n_nr_tot, 
				all_unique_haplotypes->size());

			// Dump the confusion matrix.
			fprintf(stderr, "MAP Confusion Matrix:\n");
			for (int pred_geno = 0; pred_geno < 3; pred_geno++)
			{
				fprintf(stderr, "Imp %d:", pred_geno);
				for (int known_geno = 0; known_geno < 3; known_geno++)
				{
					fprintf(stderr, "\t%.0f", MAP_confusion_matrix[pred_geno][known_geno]);
				} // known_geno loop.

				fprintf(stderr, "\n");
			} // pred_geno loop.

			fprintf(stderr, "ML Confusion Matrix:\n");
			for (int pred_geno = 0; pred_geno < 3; pred_geno++)
			{
				fprintf(stderr, "Imp %d:", pred_geno);
				for (int known_geno = 0; known_geno < 3; known_geno++)
				{
					fprintf(stderr, "\t%.0f", ML_confusion_matrix[pred_geno][known_geno]);
				} // known_geno loop.

				fprintf(stderr, "\n");
			} // pred_geno loop.

			//delete lowMAF_containing_haplotypes;
			//delete lowMAF_missing_haplotypes;
			delete all_haplotypes;
			//delete unique_lowMAF_containing_haplotypes;
			//delete unique_lowMAF_missing_haplotypes;

			for (int hap_i = 0; hap_i < n_total_haplotypes; hap_i++)
			{
				delete[] per_haplo_per_var_alleles[hap_i];
			} // hap_i loop.
			delete[] per_haplo_per_var_alleles;

			//// Set the per variant testing genotypes.
			//for (int i_s = 0; i_s < test_sample_ids->size(); i_s++)
			//{
			//	delete[] per_testing_sample_per_var_testing_geno[i_s];
			//} // i_s loop
			//delete[] per_testing_sample_per_var_testing_geno;

			//for (int i_s = 0; i_s < imp_sample_ids->size(); i_s++)
			//{
			//	delete[] per_imp_sample_per_var_testing_geno[i_s];
			//}
			//delete[] per_imp_sample_per_var_testing_geno;
		} // i_var loop.
	} // i_chr loop.

	close_f(f_probs_5col, probs_5col_fp);
	fclose(f_per_variant_stats);

	FILE* f_MAP_confusion_matrix = open_f("MAP_confusion_matrices.txt", "w");
	for (int MAP_geno = 0; MAP_geno < 3; MAP_geno++)
	{
		for (int known_geno = 0; known_geno < 3; known_geno++)
		{
			fprintf(f_MAP_confusion_matrix, "%d\t%d\t%.0f\n", MAP_geno, known_geno, MAP_confusion_matrix[MAP_geno][known_geno]);
		} // known geno loop.
	} // map geno loop. 
	fclose(f_MAP_confusion_matrix);

	FILE* f_ML_confusion_matrix = open_f("ML_confusion_matrices.txt", "w");
	for (int ML_geno = 0; ML_geno < 3; ML_geno++)
	{
		for (int known_geno = 0; known_geno < 3; known_geno++)
		{
			fprintf(f_ML_confusion_matrix, "%d\t%d\t%.0f\n", ML_geno, known_geno, ML_confusion_matrix[ML_geno][known_geno]);
		} // known geno loop.
	} // map geno loop. 
	fclose(f_ML_confusion_matrix);

	//// Re-save the imputed variants.
	//char output_genotype_matrix_bed_fp[1000];
	//sprintf(output_genotype_matrix_bed_fp, "consistency_calibrated_genotypes.txt");
	//dump_geno_sig_regs_plain(imp_target_geno_regs, imp_sample_ids, false, output_genotype_matrix_bed_fp);

	//// Re-save the imputed variants.
	//char output_genotype_signals_bed_fp[1000];
	//sprintf(output_genotype_signals_bed_fp, "consistency_calibrated_genotypes.sig");
	//FILE* f_calib_geno_sigs = open_f(output_genotype_signals_bed_fp, "w");
	//for (int i_reg = 0; i_reg < imp_target_geno_sig_regs->size(); i_reg++)
	//{
	//	double* cur_reg_signals = (double*)(imp_target_geno_sig_regs->at(i_reg)->data);
	//	fprintf(f_calib_geno_sigs, "%s\t%d\t%d\t%s", imp_target_geno_sig_regs->at(i_reg)->chrom,
	//		translate_coord(imp_target_geno_sig_regs->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
	//		translate_coord(imp_target_geno_sig_regs->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
	//		imp_target_geno_sig_regs->at(i_reg)->name);

	//	for (int i_s = 0; i_s < imp_sample_ids->size(); i_s++)
	//	{
	//		fprintf(f_calib_geno_sigs, "\t%.0f", cur_reg_signals[i_s]);
	//	}

	//	fprintf(f_calib_geno_sigs, "\n");
	//} // i_reg loop.
	//close_f(f_calib_geno_sigs, output_genotype_signals_bed_fp);
}

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
	int n_vic_vars)
{
	// Load the training data.
	fprintf(stderr, "Loading tag genotypes.\n");
	vector<t_annot_region*>* train_tag_geno_regs = load_variant_signal_regions_wrapper(train_tag_haplocoded_genotype_matrix_fp, train_sample_ids_list_fp);
	for (int i_reg = 0; i_reg < train_tag_geno_regs->size(); i_reg++)
	{
		void** old_info = (void**)(train_tag_geno_regs->at(i_reg)->data);
		void** new_info = new void*[10];
		new_info[0] = old_info[0];
		train_tag_geno_regs->at(i_reg)->data = new_info;

		train_tag_geno_regs->at(i_reg)->score = 0;
	} // i_reg loop.

	fprintf(stderr, "Loaded %d tag variants.\n", train_tag_geno_regs->size());

	fprintf(stderr, "Loading tag genotypes.\n");
	vector<t_annot_region*>* train_target_geno_regs = load_variant_signal_regions_wrapper(train_target_haplocoded_genotype_matrix_fp, train_sample_ids_list_fp);
	for (int i_reg = 0; i_reg < train_target_geno_regs->size(); i_reg++)
	{
		void** old_info = (void**)(train_target_geno_regs->at(i_reg)->data);
		void** new_info = new void*[10];
		new_info[0] = old_info[0];

		train_target_geno_regs->at(i_reg)->score = 1;
	} // i_reg loop.
	fprintf(stderr, "Loaded %d target variants.\n", train_target_geno_regs->size());

	// Load the test data.
	fprintf(stderr, "Loading testing tag genotypes.\n");
	vector<t_annot_region*>* test_tag_geno_regs = load_variant_signal_regions_wrapper(test_tag_genocoded_genotype_matrix_fp, test_sample_ids_list_fp);
	for (int i_reg = 0; i_reg < test_tag_geno_regs->size(); i_reg++)
	{
		test_tag_geno_regs->at(i_reg)->score = 0;
	} // i_reg loop.

	fprintf(stderr, "Loaded %d tag variants.\n", test_tag_geno_regs->size());

	fprintf(stderr, "Loading tag genotypes.\n");
	vector<t_annot_region*>* test_target_geno_regs = load_variant_signal_regions_wrapper(test_target_genocoded_genotype_matrix_fp, test_sample_ids_list_fp);
	for (int i_reg = 0; i_reg < test_target_geno_regs->size(); i_reg++)
	{
		test_target_geno_regs->at(i_reg)->score = 1;
	} // i_reg loop.
	fprintf(stderr, "Loaded %d target variants.\n", test_target_geno_regs->size());

	// Load the imputed data.
	vector<t_annot_region*>* imp_tag_geno_regs = load_variant_signal_regions_wrapper(imp_tag_genocoded_genotype_matrix_fp, imp_sample_ids_list_fp);
	for (int i_reg = 0; i_reg < imp_tag_geno_regs->size(); i_reg++)
	{
		imp_tag_geno_regs->at(i_reg)->score = 0;
	} // i_reg loop.

	fprintf(stderr, "Loaded %d imputed tag variants.\n", imp_tag_geno_regs->size());

	fprintf(stderr, "Loading imputed target genotypes.\n");
	vector<t_annot_region*>* imp_target_geno_regs = load_variant_signal_regions_wrapper(imp_target_genocoded_genotype_matrix_fp, imp_sample_ids_list_fp);
	for (int i_reg = 0; i_reg < imp_target_geno_regs->size(); i_reg++)
	{
		imp_target_geno_regs->at(i_reg)->score = 1;
	} // i_reg loop.
	fprintf(stderr, "Loaded %d target variants.\n", imp_target_geno_regs->size());

	//vector<t_annot_region*>* imp_target_geno_sig_regs = load_variant_signal_regions_wrapper(imp_target_geno_signal_matrix_fp, imp_sample_ids_list_fp);
	int n_loaded_samples = 0;
	vector<t_annot_region*>* imp_target_geno_sig_regs = load_signal_regs_BED(imp_target_geno_signal_matrix_fp, n_loaded_samples);
	for (int i_reg = 0; i_reg < imp_target_geno_regs->size(); i_reg++)
	{
		imp_target_geno_regs->at(i_reg)->score = 2;
	} // i_reg loop.
	fprintf(stderr, "Loaded %d target variant signal regions.\n", imp_target_geno_sig_regs->size());

	// Assign the test data to train data.
	vector<t_annot_region*>* target_intersects = intersect_regions_per_names(train_target_geno_regs, test_target_geno_regs, true);
	for (int i_int = 0; i_int < target_intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(target_intersects->at(i_int)->data);
		t_annot_region* train_target = int_info->src_reg;
		t_annot_region* test_target = int_info->dest_reg;

		void** train_reg_info = (void**)(train_target->data);
		void** test_reg_info = (void**)(test_target->data);
		train_reg_info[1] = test_reg_info[0];
	} // i_int loop.

	vector<t_annot_region*>* tag_intersects = intersect_regions_per_names(train_tag_geno_regs, test_tag_geno_regs, true);
	fprintf(stderr, "Processing %d intersects between tag genotype regions of training and imputed data.\n", tag_intersects->size());
	for (int i_int = 0; i_int < tag_intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(tag_intersects->at(i_int)->data);
		t_annot_region* train_tag = int_info->src_reg;
		t_annot_region* test_tag = int_info->dest_reg;

		void** train_reg_info = (void**)(train_tag->data);
		void** test_reg_info = (void**)(test_tag->data);
		train_reg_info[1] = test_reg_info[0];
	} // i_int loop.

	  // Assign the imp data to train data.
	vector<t_annot_region*>* imp_target_intersects = intersect_regions_per_names(train_target_geno_regs, imp_target_geno_regs, true);
	fprintf(stderr, "Processing %d intersects between target genotype regions of training and imputed data.\n", imp_target_intersects->size());
	for (int i_int = 0; i_int < imp_target_intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(imp_target_intersects->at(i_int)->data);
		t_annot_region* train_target = int_info->src_reg;
		t_annot_region* test_target = int_info->dest_reg;

		void** train_reg_info = (void**)(train_target->data);
		void** test_reg_info = (void**)(test_target->data);
		train_reg_info[2] = test_reg_info[0];
	} // i_int loop.

	  // Assign the genotype signal data to train data.
	vector<t_annot_region*>* target_signal_intersects = intersect_regions_per_names(train_target_geno_regs, imp_target_geno_sig_regs, true);
	fprintf(stderr, "Processing %d intersects between target genotype regions of training and imputed signal data.\n", target_signal_intersects->size());
	for (int i_int = 0; i_int < target_signal_intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(target_signal_intersects->at(i_int)->data);
		t_annot_region* train_target = int_info->src_reg;
		t_annot_region* imp_signal_target = int_info->dest_reg;

		void** train_reg_info = (void**)(train_target->data);
		//void** test_reg_info = (void**)(imp_signal_target->data);
		train_reg_info[3] = imp_signal_target->data;
	} // i_int loop.

	vector<t_annot_region*>* imp_tag_intersects = intersect_regions_per_names(train_tag_geno_regs, imp_tag_geno_regs, true);
	fprintf(stderr, "Processing %d intersects between tag genotype regions of training and imputed data.\n", imp_tag_intersects->size());
	for (int i_int = 0; i_int < imp_tag_intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(imp_tag_intersects->at(i_int)->data);
		t_annot_region* train_tag = int_info->src_reg;
		t_annot_region* test_tag = int_info->dest_reg;

		void** train_reg_info = (void**)(train_tag->data);
		void** test_reg_info = (void**)(test_tag->data);
		train_reg_info[2] = test_reg_info[0];
	} // i_int loop.

	vector<char*>* train_sample_ids = buffer_file(train_sample_ids_list_fp);
	vector<char*>* test_sample_ids = buffer_file(test_sample_ids_list_fp);
	vector<char*>* imp_sample_ids = buffer_file(imp_sample_ids_list_fp);

	t_restr_annot_region_list* restr_train_target_var_regs = restructure_annot_regions(train_target_geno_regs);
	t_restr_annot_region_list* restr_train_tag_var_regs = restructure_annot_regions(train_tag_geno_regs);

	double prev_n_nr_corr = 0;
	double n_false_negs_nr = 0;
	double n_nr_tot = 0;
	double calib_n_nr_corr = 0;
	double random_n_nr_corr = 0;
	for (int i_chr = 0; i_chr < restr_train_target_var_regs->chr_ids->size(); i_chr++)
	{
		vector<t_annot_region*>* cur_chr_target_var_regs = restr_train_target_var_regs->regions_per_chrom[i_chr];

		int tag_var_chr_i = t_string::get_i_str(restr_train_tag_var_regs->chr_ids, restr_train_target_var_regs->chr_ids->at(i_chr));
		vector<t_annot_region*>* cur_chr_tag_var_regs = restr_train_tag_var_regs->regions_per_chrom[tag_var_chr_i];
		for (int i_var = 0; i_var < cur_chr_target_var_regs->size(); i_var++)
		{
			if (i_var > 1 &&
				i_var % 1000 == 0)
			{
				fprintf(stderr, "Processing %d. variant: Orig: %.4f, Calib: %.4f, FN: %.4f, Rand: %.4f\n", i_var,
					prev_n_nr_corr / n_nr_tot,
					calib_n_nr_corr / n_nr_tot,
					n_false_negs_nr / n_nr_tot,
					random_n_nr_corr / n_nr_tot);
			}

			// Find the tag variants around the target variant.
			int tag_j_var = 0;
			for (int j_var = 0;
				j_var < cur_chr_tag_var_regs->size();
				j_var++)
			{
				if (cur_chr_tag_var_regs->at(j_var)->start > cur_chr_target_var_regs->at(i_var)->start)
				{
					tag_j_var = j_var;
					break;
				}
			} // j_var loop.

			  // Make sure we have the right number of tag variants.
			if (tag_j_var < n_vic_vars ||
				(tag_j_var + n_vic_vars) >= cur_chr_tag_var_regs->size())
			{
				continue;
			}

			int n_total_haplotypes = train_sample_ids->size() * 2;

			double** per_haplo_per_var_alleles = new double*[n_total_haplotypes + 2];
			for (int hap_i = 0; hap_i < n_total_haplotypes; hap_i++)
			{
				per_haplo_per_var_alleles[hap_i] = new double[2 * n_vic_vars + 5];

				// Reset all entries to -1.
				for (int i_var = 0;
					i_var < (2 * n_vic_vars + 5);
					i_var++)
				{
					per_haplo_per_var_alleles[hap_i][i_var] = -1;
				} // i_var loop.				
			} // hap_i loop.

			  // Set the per variant testing genotypes.
			double** per_testing_sample_per_var_testing_geno = new double*[test_sample_ids->size() + 1];
			for (int i_s = 0; i_s < test_sample_ids->size(); i_s++)
			{
				per_testing_sample_per_var_testing_geno[i_s] = new double[2 * n_vic_vars + 5];
				for (int i_var = 0;
					i_var < (2 * n_vic_vars + 5);
					i_var++)
				{
					per_testing_sample_per_var_testing_geno[i_s][i_var] = -1;
				} // i_var loop.				
			} // i_s loop

			double** per_imp_sample_per_var_testing_geno = new double*[imp_sample_ids->size() + 1];
			for (int i_s = 0; i_s < imp_sample_ids->size(); i_s++)
			{
				per_imp_sample_per_var_testing_geno[i_s] = new double[2 * n_vic_vars + 5];
				for (int i_var = 0;
					i_var < (2 * n_vic_vars + 5);
					i_var++)
				{
					per_imp_sample_per_var_testing_geno[i_s][i_var] = -1;
				} // i_var loop.				
			}

			vector<t_annot_region*>* pooled_tag_target_vars = new vector<t_annot_region*>();
			pooled_tag_target_vars->push_back(cur_chr_target_var_regs->at(i_var));
			for (int j_var = tag_j_var - n_vic_vars;
				j_var < tag_j_var + n_vic_vars;
				j_var++)
			{
				pooled_tag_target_vars->push_back(cur_chr_tag_var_regs->at(j_var));
			} // j_var loop.
			sort(pooled_tag_target_vars->begin(), pooled_tag_target_vars->end(), sort_regions);

			if (!t_string::compare_strings(pooled_tag_target_vars->at(n_vic_vars)->name, cur_chr_target_var_regs->at(i_var)->name))
			{
				fprintf(stderr, "Could not match target variant names: %s, %s\n\n", pooled_tag_target_vars->at(n_vic_vars)->name, cur_chr_target_var_regs->at(i_var)->name);
				exit(0);
			}

			// Setup the matrix.
			for (int j_var = 0;
				j_var < pooled_tag_target_vars->size();
				j_var++)
			{
				void** cur_var_info = (void**)(pooled_tag_target_vars->at(j_var)->data);
				char* cur_var_geno = (char*)(cur_var_info[0]);
				char* cur_test_var_geno = (char*)(cur_var_info[1]);
				char* cur_imp_var_geno = (char*)(cur_var_info[2]);

				for (int i_s = 0; i_s < train_sample_ids->size(); i_s++)
				{
					int rel_var_i = j_var;
					per_haplo_per_var_alleles[i_s * 2][rel_var_i] = get_allele_per_haplotype(cur_var_geno[i_s], 0);
					per_haplo_per_var_alleles[i_s * 2 + 1][rel_var_i] = get_allele_per_haplotype(cur_var_geno[i_s], 1);
				} // i_s loop.

				  // Copy the testing sample signal.
				for (int i_s = 0; i_s < test_sample_ids->size(); i_s++)
				{
					int rel_var_i = j_var;
					per_testing_sample_per_var_testing_geno[i_s][rel_var_i] = cur_test_var_geno[i_s];
				} // i_s loop.

				for (int i_s = 0; i_s < imp_sample_ids->size(); i_s++)
				{
					int rel_var_i = j_var;
					per_imp_sample_per_var_testing_geno[i_s][rel_var_i] = cur_imp_var_geno[i_s];
				} // i_s loop.
			} // j_var loop.

			  // Find the haplotypes with the rare alleles at the target.
			vector<double*>* lowMAF_containing_haplotypes = new vector<double*>();
			vector<double*>* lowMAF_missing_haplotypes = new vector<double*>();
			for (int i_hap = 0; i_hap < n_total_haplotypes; i_hap++)
			{
				double cur_target_allele = per_haplo_per_var_alleles[i_hap][n_vic_vars];
				if (cur_target_allele == 1)
				{
					lowMAF_containing_haplotypes->push_back(per_haplo_per_var_alleles[i_hap]);
				}
				else
				{
					lowMAF_missing_haplotypes->push_back(per_haplo_per_var_alleles[i_hap]);
				}
			} // i_hap loop.			

			  // Get the unique haplotypes.
			vector<int>* n_cnt_per_uniq_haplotypes = new vector<int>();
			vector<double*>* unique_lowMAF_containing_haplotypes = get_unique_haplotypes(lowMAF_containing_haplotypes, n_cnt_per_uniq_haplotypes);
			vector<double*>* unique_lowMAF_missing_haplotypes = get_unique_haplotypes(lowMAF_missing_haplotypes, n_cnt_per_uniq_haplotypes);

			//fprintf(stderr, "%s: %d/%d unique haplotypes are carrying the rare allele\n",
			//	cur_chr_target_var_regs->at(i_var)->name,
			//	lowMAF_containing_haplotypes->size(),
			//	unique_lowMAF_containing_haplotypes->size());

			//for (int i_hap = 0; i_hap < unique_lowMAF_containing_haplotypes->size(); i_hap++)
			//{
			//	fprintf(stderr, "%4d (%3d):", i_hap, n_cnt_per_uniq_haplotypes->at(i_hap));
			//	for (int var_i = 0; var_i < pooled_tag_target_vars->size(); var_i++)
			//	{
			//		if (var_i == n_vic_vars)
			//		{
			//			fprintf(stderr, " [[%.0f]]", (unique_lowMAF_containing_haplotypes->at(i_hap))[var_i]);
			//		}
			//		else
			//		{
			//			fprintf(stderr, " %.0f", (unique_lowMAF_containing_haplotypes->at(i_hap))[var_i]);
			//		}

			//	} // var_i loop.

			//	fprintf(stderr, "\n");
			//} // i_hap loop.

			// Go over the imputed haplotypes and evaluate the scores.
			for (int i_s = 0; i_s < test_sample_ids->size(); i_s++)
			{
				// Does this sample have non-ref genotype?
				double cur_target_allele = per_imp_sample_per_var_testing_geno[i_s][n_vic_vars];
				if (cur_target_allele > 0)
				{
					int testing_sample_geno = get_genotype_per_haplocoded_genotype(per_testing_sample_per_var_testing_geno[i_s][n_vic_vars]);
					double max_consistency = 0;
					for (int i_hap = 0; i_hap < unique_lowMAF_containing_haplotypes->size(); i_hap++)
					{
						// Go over all the variants and correlate.
						int cur_hap_consistency = 0;
						int cur_hap_tot_non_refs = 0;
						for (int var_i = 0; var_i < pooled_tag_target_vars->size(); var_i++)
						{
							if ((unique_lowMAF_containing_haplotypes->at(i_hap))[var_i] == 1)
							{
								cur_hap_tot_non_refs++;

								// Does this imputed genotype comply with the haplotype from training data?
								if (per_imp_sample_per_var_testing_geno[i_s][var_i] > 0)
								{
									cur_hap_consistency++;
								}
							}
						} // var_i loop.

						  //fprintf(stderr, "%5d [1/%d]: hap %d: %d/%d\n", i_s, testing_sample_geno, i_hap, cur_hap_consistency, cur_hap_tot_non_refs);

						max_consistency = MAX(max_consistency, (double)cur_hap_consistency / (double)cur_hap_tot_non_refs);
					} // hap_i loop

					  //fprintf(stderr, "%5d [1/%d]: %.3f\n", i_s, testing_sample_geno, max_consistency);

					// Get the imputed value.
					void** cur_target_info = (void**)(cur_chr_target_var_regs->at(i_var)->data);
					char* target_val = (char*)(cur_target_info[2]);
					double* imp_target_signal_val = (double*)(cur_target_info[3]);

					// Must count the number of variants that were correct but turned to ref.
					char orig_geno = target_val[i_s];

					n_nr_tot++;

					if (testing_sample_geno == orig_geno)
					{
						prev_n_nr_corr++;
					}

					char rand_val = orig_geno;
					if ((rand() % 5000) < 1000)
					{
						rand_val = 0;
					}

					if (testing_sample_geno == rand_val)
					{
						random_n_nr_corr++;
					}

					// Update value based on consistency.
					if (max_consistency != 1.0)
					{
						// Set it to 0.
						target_val[i_s] = 0;
						imp_target_signal_val[i_s] = 0;

						if (testing_sample_geno > 0)
						{
							n_false_negs_nr++;
						}
					}

					if (testing_sample_geno == target_val[i_s])
					{
						calib_n_nr_corr++;
					}
				} // non-ref check for the target.
			} // i_s loop.

			  /////////////// ///////////// ///////////// ///////////// ///////////// ///////////// ///////////// ///////////// 
			  /////////////// Process the unique haplotype missing genotypes.
			  //fprintf(stderr, "%s: %d/%d unique haplotypes are carrying the rare allele\n",
			  //	cur_chr_target_var_regs->at(i_var)->name,
			  //	unique_lowMAF_missing_haplotypes->size(),
			  //	unique_lowMAF_missing_haplotypes->size());

			  //for (int i_hap = 0; i_hap < unique_lowMAF_missing_haplotypes->size(); i_hap++)
			  //{
			  //	fprintf(stderr, "%4d (%3d):", i_hap, n_cnt_per_uniq_haplotypes->at(i_hap));
			  //	for (int var_i = 0; var_i < pooled_tag_target_vars->size(); var_i++)
			  //	{
			  //		if (var_i == n_vic_vars)
			  //		{
			  //			fprintf(stderr, " [[%.0f]]", (unique_lowMAF_missing_haplotypes->at(i_hap))[var_i]);
			  //		}
			  //		else
			  //		{
			  //			fprintf(stderr, " %.0f", (unique_lowMAF_missing_haplotypes->at(i_hap))[var_i]);
			  //		}

			  //	} // var_i loop.

			  //	fprintf(stderr, "\n");
			  //} // i_hap loop.

			  //fprintf(stderr, "%s: %d/%d unique missing haplotypes are carrying the rare allele\n",
			  //	cur_chr_target_var_regs->at(i_var)->name,
			  //	lowMAF_missing_haplotypes->size(),
			  //	unique_lowMAF_missing_haplotypes->size());

			  //getc(stdin);

			delete lowMAF_containing_haplotypes;
			delete lowMAF_missing_haplotypes;
			delete unique_lowMAF_containing_haplotypes;
			delete unique_lowMAF_missing_haplotypes;

			for (int hap_i = 0; hap_i < n_total_haplotypes; hap_i++)
			{
				delete[] per_haplo_per_var_alleles[hap_i];
			} // hap_i loop.
			delete[] per_haplo_per_var_alleles;

			// Set the per variant testing genotypes.
			for (int i_s = 0; i_s < test_sample_ids->size(); i_s++)
			{
				delete[] per_testing_sample_per_var_testing_geno[i_s];
			} // i_s loop
			delete[] per_testing_sample_per_var_testing_geno;

			for (int i_s = 0; i_s < imp_sample_ids->size(); i_s++)
			{
				delete[] per_imp_sample_per_var_testing_geno[i_s];
			}
			delete[] per_imp_sample_per_var_testing_geno;
		} // i_var loop.
	} // i_chr loop.

	  // Re-save the imputed variants.
	char output_genotype_matrix_bed_fp[1000];
	sprintf(output_genotype_matrix_bed_fp, "consistency_calibrated_genotypes.txt");
	dump_geno_sig_regs_plain(imp_target_geno_regs, imp_sample_ids, false, output_genotype_matrix_bed_fp);

	// Re-save the imputed variants.
	char output_genotype_signals_bed_fp[1000];
	sprintf(output_genotype_signals_bed_fp, "consistency_calibrated_genotypes.sig");
	FILE* f_calib_geno_sigs = open_f(output_genotype_signals_bed_fp, "w");
	for (int i_reg = 0; i_reg < imp_target_geno_sig_regs->size(); i_reg++)
	{
		double* cur_reg_signals = (double*)(imp_target_geno_sig_regs->at(i_reg)->data);
		fprintf(f_calib_geno_sigs, "%s\t%d\t%d\t%s", imp_target_geno_sig_regs->at(i_reg)->chrom, 
				translate_coord(imp_target_geno_sig_regs->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base), 
				translate_coord(imp_target_geno_sig_regs->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base), 
				imp_target_geno_sig_regs->at(i_reg)->name);

		for (int i_s = 0; i_s < imp_sample_ids->size(); i_s++)
		{
			fprintf(f_calib_geno_sigs, "\t%.0f", cur_reg_signals[i_s]);
		}

		fprintf(f_calib_geno_sigs, "\n");
	} // i_reg loop.
	close_f(f_calib_geno_sigs, output_genotype_signals_bed_fp);
}

void assign_haplotype_consistency_scores_per_low_MAF_imputations(char* train_tag_haplocoded_genotype_matrix_fp,
	char* train_target_haplocoded_genotype_matrix_fp,
	char* train_sample_ids_list_fp,	
	char* imp_tag_genocoded_genotype_matrix_fp,
	char* imp_target_genocoded_genotype_matrix_fp,
	char* imp_sample_ids_list_fp,
	char* test_tag_genocoded_genotype_matrix_fp, 
	char* test_target_genocoded_genotype_matrix_fp,
	char* test_sample_ids_list_fp,
	int n_vic_vars)
{
	// Load the training data.
	fprintf(stderr, "Loading tag genotypes.\n");
	vector<t_annot_region*>* train_tag_geno_regs = load_variant_signal_regions_wrapper(train_tag_haplocoded_genotype_matrix_fp, train_sample_ids_list_fp);
	for (int i_reg = 0; i_reg < train_tag_geno_regs->size(); i_reg++)
	{
		void** old_info = (void**)(train_tag_geno_regs->at(i_reg)->data);
		void** new_info = new void*[10];
		new_info[0] = old_info[0];
		train_tag_geno_regs->at(i_reg)->data = new_info;

		train_tag_geno_regs->at(i_reg)->score = 0;
	} // i_reg loop.

	fprintf(stderr, "Loaded %d tag variants.\n", train_tag_geno_regs->size());

	fprintf(stderr, "Loading tag genotypes.\n");
	vector<t_annot_region*>* train_target_geno_regs = load_variant_signal_regions_wrapper(train_target_haplocoded_genotype_matrix_fp, train_sample_ids_list_fp);
	for (int i_reg = 0; i_reg < train_target_geno_regs->size(); i_reg++)
	{
		void** old_info = (void**)(train_target_geno_regs->at(i_reg)->data);
		void** new_info = new void*[10];
		new_info[0] = old_info[0];

		train_target_geno_regs->at(i_reg)->score = 1;
	} // i_reg loop.
	fprintf(stderr, "Loaded %d target variants.\n", train_target_geno_regs->size());

	// Load the test data.
	fprintf(stderr, "Loading testing tag genotypes.\n");
	vector<t_annot_region*>* test_tag_geno_regs = load_variant_signal_regions_wrapper(test_tag_genocoded_genotype_matrix_fp, test_sample_ids_list_fp);
	for (int i_reg = 0; i_reg < test_tag_geno_regs->size(); i_reg++)
	{
		test_tag_geno_regs->at(i_reg)->score = 0;
	} // i_reg loop.

	fprintf(stderr, "Loaded %d tag variants.\n", test_tag_geno_regs->size());

	fprintf(stderr, "Loading tag genotypes.\n");
	vector<t_annot_region*>* test_target_geno_regs = load_variant_signal_regions_wrapper(test_target_genocoded_genotype_matrix_fp, test_sample_ids_list_fp);
	for (int i_reg = 0; i_reg < test_target_geno_regs->size(); i_reg++)
	{
		test_target_geno_regs->at(i_reg)->score = 1;
	} // i_reg loop.
	fprintf(stderr, "Loaded %d target variants.\n", test_target_geno_regs->size());

	// Load the imputed data.
	vector<t_annot_region*>* imp_tag_geno_regs = load_variant_signal_regions_wrapper(imp_tag_genocoded_genotype_matrix_fp, imp_sample_ids_list_fp);
	for (int i_reg = 0; i_reg < imp_tag_geno_regs->size(); i_reg++)
	{
		imp_tag_geno_regs->at(i_reg)->score = 0;
	} // i_reg loop.

	fprintf(stderr, "Loaded %d tag variants.\n", imp_tag_geno_regs->size());

	fprintf(stderr, "Loading tag genotypes.\n");
	vector<t_annot_region*>* imp_target_geno_regs = load_variant_signal_regions_wrapper(imp_target_genocoded_genotype_matrix_fp, imp_sample_ids_list_fp);
	for (int i_reg = 0; i_reg < imp_target_geno_regs->size(); i_reg++)
	{
		imp_target_geno_regs->at(i_reg)->score = 1;
	} // i_reg loop.
	fprintf(stderr, "Loaded %d target variants.\n", imp_target_geno_regs->size());

	// Assign the test data to train data.
	vector<t_annot_region*>* target_intersects = intersect_regions_per_names(train_target_geno_regs, test_target_geno_regs, true);
	for (int i_int = 0; i_int < target_intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(target_intersects->at(i_int)->data);
		t_annot_region* train_target = int_info->src_reg;
		t_annot_region* test_target = int_info->dest_reg;

		void** train_reg_info = (void**)(train_target->data);
		void** test_reg_info = (void**)(test_target->data);
		train_reg_info[1] = test_reg_info[0];
	} // i_int loop.

	vector<t_annot_region*>* tag_intersects = intersect_regions_per_names(train_tag_geno_regs, test_tag_geno_regs, true);
	for (int i_int = 0; i_int < tag_intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(tag_intersects->at(i_int)->data);
		t_annot_region* train_tag = int_info->src_reg;
		t_annot_region* test_tag = int_info->dest_reg;

		void** train_reg_info = (void**)(train_tag->data);
		void** test_reg_info = (void**)(test_tag->data);
		train_reg_info[1] = test_reg_info[0];
	} // i_int loop.

	// Assign the imp data to train data.
	vector<t_annot_region*>* imp_target_intersects = intersect_regions_per_names(train_target_geno_regs, imp_target_geno_regs, true);
	for (int i_int = 0; i_int < imp_target_intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(imp_target_intersects->at(i_int)->data);
		t_annot_region* train_target = int_info->src_reg;
		t_annot_region* test_target = int_info->dest_reg;

		void** train_reg_info = (void**)(train_target->data);
		void** test_reg_info = (void**)(test_target->data);
		train_reg_info[2] = test_reg_info[0];
	} // i_int loop.

	vector<t_annot_region*>* imp_tag_intersects = intersect_regions_per_names(train_tag_geno_regs, imp_tag_geno_regs, true);
	for (int i_int = 0; i_int < imp_tag_intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(imp_tag_intersects->at(i_int)->data);
		t_annot_region* train_tag = int_info->src_reg;
		t_annot_region* test_tag = int_info->dest_reg;

		void** train_reg_info = (void**)(train_tag->data);
		void** test_reg_info = (void**)(test_tag->data);
		train_reg_info[2] = test_reg_info[0];
	} // i_int loop.

	vector<char*>* train_sample_ids = buffer_file(train_sample_ids_list_fp);
	vector<char*>* test_sample_ids = buffer_file(test_sample_ids_list_fp);
	vector<char*>* imp_sample_ids = buffer_file(imp_sample_ids_list_fp);

	t_restr_annot_region_list* restr_train_target_var_regs = restructure_annot_regions(train_target_geno_regs);
	t_restr_annot_region_list* restr_train_tag_var_regs = restructure_annot_regions(train_tag_geno_regs);

	double prev_n_nr_corr = 0;
	double n_false_negs_nr = 0;
	double n_nr_tot = 0;
	double calib_n_nr_corr = 0;
	double random_n_nr_corr = 0;
	for (int i_chr = 0; i_chr < restr_train_target_var_regs->chr_ids->size(); i_chr++)
	{
		vector<t_annot_region*>* cur_chr_target_var_regs = restr_train_target_var_regs->regions_per_chrom[i_chr];

		int tag_var_chr_i = t_string::get_i_str(restr_train_tag_var_regs->chr_ids, restr_train_target_var_regs->chr_ids->at(i_chr));
		vector<t_annot_region*>* cur_chr_tag_var_regs = restr_train_tag_var_regs->regions_per_chrom[tag_var_chr_i];
		for (int i_var = 0; i_var < cur_chr_target_var_regs->size(); i_var++)
		{
			if (i_var > 1 &&
				i_var % 1000 == 0)
			{
				fprintf(stderr, "Processing %d. variant: Orig: %.4f, Calib: %.4f, FN: %.4f, Rand: %.4f\n", i_var,
					prev_n_nr_corr / n_nr_tot,
					calib_n_nr_corr/ n_nr_tot,
					n_false_negs_nr / n_nr_tot,
					random_n_nr_corr / n_nr_tot);
			}

			// Find the tag variants around the target variant.
			int tag_j_var = 0;
			for (int j_var = 0;
				j_var < cur_chr_tag_var_regs->size();
				j_var++)
			{
				if (cur_chr_tag_var_regs->at(j_var)->start > cur_chr_target_var_regs->at(i_var)->start)
				{
					tag_j_var = j_var;
					break;
				}
			} // j_var loop.

			  // Make sure we have the right number of tag variants.
			if (tag_j_var < n_vic_vars ||
				(tag_j_var + n_vic_vars) >= cur_chr_tag_var_regs->size())
			{
				continue;
			}

			int n_total_haplotypes = train_sample_ids->size() * 2;

			double** per_haplo_per_var_alleles = new double*[n_total_haplotypes + 2];			
			for (int hap_i = 0; hap_i < n_total_haplotypes; hap_i++)
			{
				per_haplo_per_var_alleles[hap_i] = new double[2 * n_vic_vars + 5];

				// Reset all entries to -1.
				for (int i_var = 0;
					i_var < (2 * n_vic_vars + 5);
					i_var++)
				{
					per_haplo_per_var_alleles[hap_i][i_var] = -1;
				} // i_var loop.				
			} // hap_i loop.

			// Set the per variant testing genotypes.
			double** per_testing_sample_per_var_testing_geno = new double*[test_sample_ids->size() + 1];
			for (int i_s = 0; i_s < test_sample_ids->size(); i_s++)
			{
				per_testing_sample_per_var_testing_geno[i_s] = new double[2 * n_vic_vars + 5];
				for (int i_var = 0;
					i_var < (2 * n_vic_vars + 5);
					i_var++)
				{
					per_testing_sample_per_var_testing_geno[i_s][i_var] = -1;
				} // i_var loop.				
			} // i_s loop

			double** per_imp_sample_per_var_testing_geno = new double*[imp_sample_ids->size() + 1];
			for (int i_s = 0; i_s < imp_sample_ids->size(); i_s++)
			{
				per_imp_sample_per_var_testing_geno[i_s] = new double[2 * n_vic_vars + 5];
				for (int i_var = 0;
					i_var < (2 * n_vic_vars + 5);
					i_var++)
				{
					per_imp_sample_per_var_testing_geno[i_s][i_var] = -1;
				} // i_var loop.				
			}

			vector<t_annot_region*>* pooled_tag_target_vars = new vector<t_annot_region*>();
			pooled_tag_target_vars->push_back(cur_chr_target_var_regs->at(i_var));
			for (int j_var = tag_j_var - n_vic_vars;
				j_var < tag_j_var + n_vic_vars;
				j_var++)
			{
				pooled_tag_target_vars->push_back(cur_chr_tag_var_regs->at(j_var));
			} // j_var loop.
			sort(pooled_tag_target_vars->begin(), pooled_tag_target_vars->end(), sort_regions);

			if (!t_string::compare_strings(pooled_tag_target_vars->at(n_vic_vars)->name, cur_chr_target_var_regs->at(i_var)->name))
			{
				fprintf(stderr, "Could not match target variant names: %s, %s\n\n", pooled_tag_target_vars->at(n_vic_vars)->name, cur_chr_target_var_regs->at(i_var)->name);
				exit(0);
			}

			// Setup the matrix.
			for (int j_var = 0;
				j_var < pooled_tag_target_vars->size();
				j_var++)
			{
				void** cur_var_info = (void**)(pooled_tag_target_vars->at(j_var)->data);
				char* cur_var_geno = (char*)(cur_var_info[0]);
				char* cur_test_var_geno = (char*)(cur_var_info[1]);
				char* cur_imp_var_geno = (char*)(cur_var_info[2]);

				for (int i_s = 0; i_s < train_sample_ids->size(); i_s++)
				{
					int rel_var_i = j_var;
					per_haplo_per_var_alleles[i_s * 2][rel_var_i] = get_allele_per_haplotype(cur_var_geno[i_s], 0);
					per_haplo_per_var_alleles[i_s * 2 + 1][rel_var_i] = get_allele_per_haplotype(cur_var_geno[i_s], 1);
				} // i_s loop.

				// Copy the testing sample signal.
				for (int i_s = 0; i_s < test_sample_ids->size(); i_s++)
				{
					int rel_var_i = j_var;
					per_testing_sample_per_var_testing_geno[i_s][rel_var_i] = cur_test_var_geno[i_s];
				} // i_s loop.

				for (int i_s = 0; i_s < imp_sample_ids->size(); i_s++)
				{
					int rel_var_i = j_var;
					per_imp_sample_per_var_testing_geno[i_s][rel_var_i] = cur_imp_var_geno[i_s];
				} // i_s loop.
			} // j_var loop.

			// Find the haplotypes with the rare alleles at the target.
			vector<double*>* lowMAF_containing_haplotypes = new vector<double*>();
			vector<double*>* lowMAF_missing_haplotypes = new vector<double*>();
			for (int i_hap = 0; i_hap < n_total_haplotypes; i_hap++)
			{
				double cur_target_allele = per_haplo_per_var_alleles[i_hap][n_vic_vars];
				if (cur_target_allele == 1)
				{
					lowMAF_containing_haplotypes->push_back(per_haplo_per_var_alleles[i_hap]);
				}
				else
				{
					lowMAF_missing_haplotypes->push_back(per_haplo_per_var_alleles[i_hap]);
				}
			} // i_hap loop.			

			// Get the unique haplotypes.
			vector<int>* n_cnt_per_uniq_haplotypes = new vector<int>();
			vector<double*>* unique_lowMAF_containing_haplotypes = get_unique_haplotypes(lowMAF_containing_haplotypes, n_cnt_per_uniq_haplotypes);
			vector<double*>* unique_lowMAF_missing_haplotypes = get_unique_haplotypes(lowMAF_missing_haplotypes, n_cnt_per_uniq_haplotypes);

			//fprintf(stderr, "%s: %d/%d unique haplotypes are carrying the rare allele\n",
			//	cur_chr_target_var_regs->at(i_var)->name,
			//	lowMAF_containing_haplotypes->size(),
			//	unique_lowMAF_containing_haplotypes->size());

			//for (int i_hap = 0; i_hap < unique_lowMAF_containing_haplotypes->size(); i_hap++)
			//{
			//	fprintf(stderr, "%4d (%3d):", i_hap, n_cnt_per_uniq_haplotypes->at(i_hap));
			//	for (int var_i = 0; var_i < pooled_tag_target_vars->size(); var_i++)
			//	{
			//		if (var_i == n_vic_vars)
			//		{
			//			fprintf(stderr, " [[%.0f]]", (unique_lowMAF_containing_haplotypes->at(i_hap))[var_i]);
			//		}
			//		else
			//		{
			//			fprintf(stderr, " %.0f", (unique_lowMAF_containing_haplotypes->at(i_hap))[var_i]);
			//		}

			//	} // var_i loop.

			//	fprintf(stderr, "\n");
			//} // i_hap loop.

			// Go over the imputed haplotypes and evaluate the scores.
			for (int i_s = 0; i_s < test_sample_ids->size(); i_s++)
			{
				// Does this sample have 
				double cur_target_allele = per_imp_sample_per_var_testing_geno[i_s][n_vic_vars];
				if (cur_target_allele > 0)
				{
					int testing_sample_geno = get_genotype_per_haplocoded_genotype(per_testing_sample_per_var_testing_geno[i_s][n_vic_vars]);
					double max_consistency = 0;
					for (int i_hap = 0; i_hap < unique_lowMAF_containing_haplotypes->size(); i_hap++)
					{
						// Go over all the variants and correlate.
						int cur_hap_consistency = 0;
						int cur_hap_tot_non_refs = 0;
						for (int var_i = 0; var_i < pooled_tag_target_vars->size(); var_i++)
						{
							if((unique_lowMAF_containing_haplotypes->at(i_hap))[var_i] == 1)
							{
								cur_hap_tot_non_refs++;

								// Does this imputed genotype comply with the haplotype from training data?
								if (per_imp_sample_per_var_testing_geno[i_s][var_i] > 0)
								{
									cur_hap_consistency++;
								}
							}
						} // var_i loop.

						//fprintf(stderr, "%5d [1/%d]: hap %d: %d/%d\n", i_s, testing_sample_geno, i_hap, cur_hap_consistency, cur_hap_tot_non_refs);

						max_consistency = MAX(max_consistency, (double)cur_hap_consistency / (double)cur_hap_tot_non_refs);
					} // hap_i loop

					//fprintf(stderr, "%5d [1/%d]: %.3f\n", i_s, testing_sample_geno, max_consistency);

					// Get the imputed value.
					void** cur_target_info = (void**)(cur_chr_target_var_regs->at(i_var)->data);
					char* target_val = (char*)(cur_target_info[2]);

					// Must count the number of variants that were correct but turned to ref.

					char orig_geno = target_val[i_s];

					n_nr_tot++;

					if(testing_sample_geno == orig_geno)
					{
						prev_n_nr_corr++;
					}

					char rand_val = orig_geno;
					if ((rand() % 5000) < 1000)
					{
						rand_val = 0;
					}

					if (testing_sample_geno == rand_val)
					{
						random_n_nr_corr++;
					}

					// Update value based on consistency.
					if (max_consistency != 1.0)
					{
						// Set it to 0.
						target_val[i_s] = 0;

						if (testing_sample_geno > 0)
						{
							n_false_negs_nr++;
						}
					}

					if (testing_sample_geno == target_val[i_s])
					{
						calib_n_nr_corr++;
					}
				} // non-ref check for the target.
			} // i_s loop.

			/////////////// ///////////// ///////////// ///////////// ///////////// ///////////// ///////////// ///////////// 
			/////////////// Process the unique haplotype missing genotypes.
			//fprintf(stderr, "%s: %d/%d unique haplotypes are carrying the rare allele\n",
			//	cur_chr_target_var_regs->at(i_var)->name,
			//	unique_lowMAF_missing_haplotypes->size(),
			//	unique_lowMAF_missing_haplotypes->size());

			//for (int i_hap = 0; i_hap < unique_lowMAF_missing_haplotypes->size(); i_hap++)
			//{
			//	fprintf(stderr, "%4d (%3d):", i_hap, n_cnt_per_uniq_haplotypes->at(i_hap));
			//	for (int var_i = 0; var_i < pooled_tag_target_vars->size(); var_i++)
			//	{
			//		if (var_i == n_vic_vars)
			//		{
			//			fprintf(stderr, " [[%.0f]]", (unique_lowMAF_missing_haplotypes->at(i_hap))[var_i]);
			//		}
			//		else
			//		{
			//			fprintf(stderr, " %.0f", (unique_lowMAF_missing_haplotypes->at(i_hap))[var_i]);
			//		}

			//	} // var_i loop.

			//	fprintf(stderr, "\n");
			//} // i_hap loop.

			//fprintf(stderr, "%s: %d/%d unique missing haplotypes are carrying the rare allele\n",
			//	cur_chr_target_var_regs->at(i_var)->name,
			//	lowMAF_missing_haplotypes->size(),
			//	unique_lowMAF_missing_haplotypes->size());

			//getc(stdin);

			delete lowMAF_containing_haplotypes;
			delete lowMAF_missing_haplotypes;
			delete unique_lowMAF_containing_haplotypes;
			delete unique_lowMAF_missing_haplotypes;

			for (int hap_i = 0; hap_i < n_total_haplotypes; hap_i++)
			{
				delete [] per_haplo_per_var_alleles[hap_i];
			} // hap_i loop.
			delete[] per_haplo_per_var_alleles;

			  // Set the per variant testing genotypes.
			for (int i_s = 0; i_s < test_sample_ids->size(); i_s++)
			{
				delete [] per_testing_sample_per_var_testing_geno[i_s];
			} // i_s loop
			delete[] per_testing_sample_per_var_testing_geno;

			for (int i_s = 0; i_s < imp_sample_ids->size(); i_s++)
			{
				delete [] per_imp_sample_per_var_testing_geno[i_s];
			}
			delete[] per_imp_sample_per_var_testing_geno;
		} // i_var loop.
	} // i_chr loop.

	// Re-save the imputed variants.
	char output_genotype_matrix_bed_fp[1000];
	sprintf(output_genotype_matrix_bed_fp, "consistency_calibrated_genotypes.txt");
	dump_geno_sig_regs_plain(imp_target_geno_regs, imp_sample_ids, false, output_genotype_matrix_bed_fp);
}

void analyze_low_MAF_haplotype_statistics_per_target_variants (char* train_tag_haplocoded_genotype_matrix_fp,
																char* train_target_haplocoded_genotype_matrix_fp,
																char* train_sample_ids_list_fp,
																int n_vic_vars)
{
	fprintf(stderr, "Loading tag genotypes.\n");
	vector<t_annot_region*>* train_tag_geno_regs = load_variant_signal_regions_wrapper(train_tag_haplocoded_genotype_matrix_fp, train_sample_ids_list_fp);
	for (int i_reg = 0; i_reg < train_tag_geno_regs->size(); i_reg++)
	{
		train_tag_geno_regs->at(i_reg)->score = 0;
	} // i_reg loop.
	fprintf(stderr, "Loaded %d tag variants.\n", train_tag_geno_regs->size());

	fprintf(stderr, "Loading tag genotypes.\n");
	vector<t_annot_region*>* train_target_geno_regs = load_variant_signal_regions_wrapper(train_target_haplocoded_genotype_matrix_fp, train_sample_ids_list_fp);
	for (int i_reg = 0; i_reg < train_target_geno_regs->size(); i_reg++)
	{
		train_target_geno_regs->at(i_reg)->score = 1;
	} // i_reg loop.
	fprintf(stderr, "Loaded %d target variants.\n", train_target_geno_regs->size());

	vector<char*>* train_sample_ids = buffer_file(train_sample_ids_list_fp);

	//// Put the target and tag genotypes and 
	//vector<t_annot_region*>* train_target_tag_geno_regs = new vector<t_annot_region*>();
	//train_target_tag_geno_regs->insert(train_target_tag_geno_regs->end(), train_target_geno_regs->begin(), train_target_geno_regs->end());
	//train_target_tag_geno_regs->insert(train_target_tag_geno_regs->end(), train_tag_geno_regs->begin(), train_tag_geno_regs->end());

	t_restr_annot_region_list* restr_train_target_var_regs = restructure_annot_regions(train_target_geno_regs);
	t_restr_annot_region_list* restr_train_tag_var_regs = restructure_annot_regions(train_tag_geno_regs);

	for (int i_chr = 0; i_chr < restr_train_target_var_regs->chr_ids->size(); i_chr++)
	{
		vector<t_annot_region*>* cur_chr_target_var_regs = restr_train_target_var_regs->regions_per_chrom[i_chr];

		int tag_var_chr_i = t_string::get_i_str(restr_train_tag_var_regs->chr_ids, restr_train_target_var_regs->chr_ids->at(i_chr));
		vector<t_annot_region*>* cur_chr_tag_var_regs = restr_train_tag_var_regs->regions_per_chrom[tag_var_chr_i];
		for (int i_var = 0; i_var < cur_chr_target_var_regs->size(); i_var++)
		{
			// Find the tag variants around the target variant.
			int tag_j_var = 0;
			for (int j_var = 0;
				j_var < cur_chr_tag_var_regs->size();
				j_var++)
			{
				if (cur_chr_tag_var_regs->at(j_var)->start > cur_chr_target_var_regs->at(i_var)->start)
				{
					tag_j_var = j_var;
					break;
				}
			} // j_var loop.

			// Make sure we have the right number of tag variants.
			if (tag_j_var < n_vic_vars ||
				(tag_j_var + n_vic_vars) >= cur_chr_tag_var_regs->size())
			{
				continue;
			}

			int n_total_haplotypes = train_sample_ids->size() * 2;

			double** per_haplo_per_var_alleles = new double*[n_total_haplotypes + 2];
			for (int hap_i = 0; hap_i < n_total_haplotypes; hap_i++)
			{
				per_haplo_per_var_alleles[hap_i] = new double[2 * n_vic_vars + 5];

				// Reset all entries to -1.
				for(int i_var = 0; 
					i_var < (2 * n_vic_vars + 5); 
					i_var++)
				{ 
					per_haplo_per_var_alleles[hap_i][i_var] = -1;
				} // i_var loop.				
			} // hap_i loop.

			vector<t_annot_region*>* pooled_tag_target_vars = new vector<t_annot_region*>();
			pooled_tag_target_vars->push_back(cur_chr_target_var_regs->at(i_var));
			for (int j_var = tag_j_var - n_vic_vars; 
				j_var < tag_j_var + n_vic_vars;
				j_var++)
			{
				pooled_tag_target_vars->push_back(cur_chr_tag_var_regs->at(j_var));
			} // j_var loop.
			sort(pooled_tag_target_vars->begin(), pooled_tag_target_vars->end(), sort_regions);

			if (!t_string::compare_strings(pooled_tag_target_vars->at(n_vic_vars)->name, cur_chr_target_var_regs->at(i_var)->name))
			{
				fprintf(stderr, "Could not match target variant names: %s, %s\n\n", pooled_tag_target_vars->at(n_vic_vars)->name, cur_chr_target_var_regs->at(i_var)->name);
				exit(0);
			}

			// Setup the matrix.
			for (int j_var = 0;
				j_var < pooled_tag_target_vars->size();
				j_var++)
			{
				void** cur_var_info = (void**)(pooled_tag_target_vars->at(j_var)->data);
				char* cur_var_geno = (char*)(cur_var_info[0]);

				for (int i_s = 0; i_s < train_sample_ids->size(); i_s++)
				{
					int rel_var_i = j_var;
					per_haplo_per_var_alleles[i_s * 2][rel_var_i] = get_allele_per_haplotype(cur_var_geno[i_s], 0);
					per_haplo_per_var_alleles[i_s * 2 + 1][rel_var_i] = get_allele_per_haplotype(cur_var_geno[i_s], 1);
				} // i_s loop.
			} // j_var loop.

			// Find the haplotypes with the rare alleles at the target.
			vector<double*>* lowMAF_containing_haplotypes = new vector<double*>();
			vector<double*>* lowMAF_missing_haplotypes = new vector<double*>();
			for (int i_hap = 0; i_hap < n_total_haplotypes; i_hap++)
			{
				double cur_target_allele = per_haplo_per_var_alleles[i_hap][n_vic_vars];
				if (cur_target_allele == 1)
				{
					lowMAF_containing_haplotypes->push_back(per_haplo_per_var_alleles[i_hap]);
				}
				else
				{
					lowMAF_missing_haplotypes->push_back(per_haplo_per_var_alleles[i_hap]);
				}
			} // i_hap loop.			

			// Get the unique haplotypes.
			vector<int>* n_cnt_per_uniq_haplotypes = new vector<int>();
			vector<double*>* unique_lowMAF_containing_haplotypes = get_unique_haplotypes(lowMAF_containing_haplotypes, n_cnt_per_uniq_haplotypes);
			vector<double*>* unique_lowMAF_missing_haplotypes = get_unique_haplotypes(lowMAF_missing_haplotypes, n_cnt_per_uniq_haplotypes);

			fprintf(stderr, "%s: %d/%d unique haplotypes are carrying the rare allele\n",
					cur_chr_target_var_regs->at(i_var)->name,
					lowMAF_containing_haplotypes->size(),
					unique_lowMAF_containing_haplotypes->size());

			for (int i_hap = 0; i_hap < unique_lowMAF_containing_haplotypes->size(); i_hap++)
			{
				fprintf(stderr, "%4d (%3d):", i_hap, n_cnt_per_uniq_haplotypes->at(i_hap));
				for (int var_i = 0; var_i < pooled_tag_target_vars->size(); var_i++)
				{
					if (var_i == n_vic_vars)
					{
						fprintf(stderr, " [[%.0f]]", (unique_lowMAF_containing_haplotypes->at(i_hap))[var_i]);
					}
					else
					{
						fprintf(stderr, " %.0f", (unique_lowMAF_containing_haplotypes->at(i_hap))[var_i]);
					}

				} // var_i loop.

				fprintf(stderr, "\n");
			} // i_hap loop.

			fprintf(stderr, "%s: %d/%d unique haplotypes are carrying the rare allele\n",
				cur_chr_target_var_regs->at(i_var)->name,
				lowMAF_containing_haplotypes->size(),
				unique_lowMAF_containing_haplotypes->size());

			for (int i_hap = 0; i_hap < unique_lowMAF_missing_haplotypes->size(); i_hap++)
			{
				fprintf(stderr, "%4d (%3d):", i_hap, n_cnt_per_uniq_haplotypes->at(i_hap));
				for (int var_i = 0; var_i < pooled_tag_target_vars->size(); var_i++)
				{
					if (var_i == n_vic_vars)
					{
						fprintf(stderr, " [[%.0f]]", (unique_lowMAF_missing_haplotypes->at(i_hap))[var_i]);
					}
					else
					{
						fprintf(stderr, " %.0f", (unique_lowMAF_missing_haplotypes->at(i_hap))[var_i]);
					}

				} // var_i loop.

				fprintf(stderr, "\n");
			} // i_hap loop.

			fprintf(stderr, "%s: %d/%d unique missing haplotypes are carrying the rare allele\n",
				cur_chr_target_var_regs->at(i_var)->name,
				lowMAF_missing_haplotypes->size(),
				unique_lowMAF_missing_haplotypes->size());

			getc(stdin);
		} // i_var loop.
	} // i_chr loop.
}

void generate_trio_genotyping_imputability_stats(char* genotype_signal_matrix_fp, 
												char* sample_ids_fp, 
												char* trio_fmc_sample_ids_fp)
{
	fprintf(stderr, "Loading trio genotype signals.\n");
	vector<char*>* sample_ids = buffer_file(sample_ids_fp);
	fprintf(stderr, "Loading genotype matrix with %d testing samples.\n", sample_ids->size());
	vector<t_annot_region*>* trio_genotype_signal_regs = load_variant_genotype_signal_regions(genotype_signal_matrix_fp, sample_ids);
	fprintf(stderr, "Loaded %d tag genotype signal regions with %d testing samples.\n", trio_genotype_signal_regs->size(), sample_ids->size());

	vector<char*>* trio_fmc_list = buffer_file(trio_fmc_sample_ids_fp);

	t_restr_annot_region_list* restr_trio_genotype_signal_regs = restructure_annot_regions(trio_genotype_signal_regs);

	FILE* f_fam_allele_stats = open_f("family_allelic_stats.txt", "w");
	for (int i_fam = 0; i_fam < trio_fmc_list->size(); i_fam++)
	{
		int child_sample_i = -1;
		int pat_sample_i = -1;
		int mat_sample_i = -1;
		char pat_sample_id[100];
		char mat_sample_id[100];
		char child_sample_id[100];
		if (sscanf(trio_fmc_list->at(i_fam), "%s %s %s", child_sample_id, mat_sample_id, pat_sample_id) != 3)
		{
			fprintf(stderr, "Could not parse the child-pat-mat sample id's: %s\n", trio_fmc_list->at(i_fam));
			continue;
		}
		else
		{			
			child_sample_i = t_string::get_i_str(sample_ids, child_sample_id);
			pat_sample_i = t_string::get_i_str(sample_ids, pat_sample_id);
			mat_sample_i = t_string::get_i_str(sample_ids, mat_sample_id);

			if (child_sample_i == sample_ids->size() ||
				mat_sample_i == sample_ids->size() ||
				pat_sample_i == sample_ids->size())
			{
				fprintf(stderr, "Could not find the child-pat-mat sample id's: %s\n", trio_fmc_list->at(i_fam));
				continue;
			}
		}

		// If we are here, we have identified the samples.
		fprintf(stderr, "Processing %s, %s, %s\n", child_sample_id, mat_sample_id, pat_sample_id);
		int n_4_copies = 0;
		int n_3_copies = 0;
		int n_2_copies = 0;
		int n_hom_hom = 0;
		int n_het_het = 0;
		int n_all_vars = 0;
		int n_consistent = 0;

		for (int i_chr = 0; i_chr < restr_trio_genotype_signal_regs->chr_ids->size(); i_chr++)
		{
			vector<t_annot_region*>* genotype_sig_regs = restr_trio_genotype_signal_regs->regions_per_chrom[i_chr];

			for (int i_reg = 0; i_reg < genotype_sig_regs->size(); i_reg++)
			{
				n_all_vars++;

				void** cur_reg_info = (void**)(genotype_sig_regs->at(i_reg)->data);
				char* cur_reg_geno_sigs = (char*)(cur_reg_info[0]);

				// 4 copies check.
				if ((cur_reg_geno_sigs[mat_sample_i] + cur_reg_geno_sigs[pat_sample_i]) == 4 ||
					(cur_reg_geno_sigs[mat_sample_i] + cur_reg_geno_sigs[pat_sample_i]) == 0)
				{
					if ((cur_reg_geno_sigs[mat_sample_i] + cur_reg_geno_sigs[pat_sample_i]) == 4 &&
						cur_reg_geno_sigs[child_sample_i] != 2)
					{
						fprintf(stderr, "Sanity check failed: Child is not 2.\n");
					}
					else if ((cur_reg_geno_sigs[mat_sample_i] + cur_reg_geno_sigs[pat_sample_i]) == 0 &&
						cur_reg_geno_sigs[child_sample_i] != 0)
					{
						fprintf(stderr, "Sanity check failed: Child is not 0.\n");
					}
					else
					{
						n_consistent++;
						n_4_copies++;
					}
				}

				// 3 copies check.
				if ((cur_reg_geno_sigs[mat_sample_i] + cur_reg_geno_sigs[pat_sample_i]) == 3 || 
					(cur_reg_geno_sigs[mat_sample_i] + cur_reg_geno_sigs[pat_sample_i]) == 1)
				{
					if ((cur_reg_geno_sigs[mat_sample_i] + cur_reg_geno_sigs[pat_sample_i]) == 3 &&
						cur_reg_geno_sigs[child_sample_i] == 0)
					{
						fprintf(stderr, "Sanity check failed: Child is 0 @ 2-1 parent genotype.\n");
					}
					else if ((cur_reg_geno_sigs[mat_sample_i] + cur_reg_geno_sigs[pat_sample_i]) == 1 &&
						cur_reg_geno_sigs[child_sample_i] == 2)
					{
						fprintf(stderr, "Sanity check failed: Child is 2 @ 0-1 parent genotype.\n");
					}
					else
					{
						n_consistent++;
						n_3_copies++;
					}
				}

				if ((cur_reg_geno_sigs[mat_sample_i] + cur_reg_geno_sigs[pat_sample_i]) == 2)
				{
					if (cur_reg_geno_sigs[mat_sample_i] == 0 || cur_reg_geno_sigs[mat_sample_i] == 2)
					{
						if (cur_reg_geno_sigs[child_sample_i] != 1)
						{
							fprintf(stderr, "Sanity check failed: Child is not 1 @ hom-hom.\n");
						}
						else
						{
							n_2_copies++;
							n_consistent++;
							n_hom_hom++;
						}
						
					}
					// 2 copies check: Hardest to decode.
					else if (cur_reg_geno_sigs[mat_sample_i] == 1 && cur_reg_geno_sigs[pat_sample_i] == 1)
					{
						n_consistent++;
						n_2_copies++;
						n_het_het++;
					}
				}
			} // i_reg loop.
		} // i_chr loop.

		fprintf(f_fam_allele_stats, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", trio_fmc_list->at(i_fam), n_hom_hom, n_het_het, n_2_copies, n_3_copies, n_4_copies, n_consistent, n_all_vars);
	} // i_fam loop.

	fclose(f_fam_allele_stats);
}

void save_train_test_matrices_per_variant(char* roi_target_vars_BED_fp,
										int n_half_vicinity_tags,
										char* training_tag_geno_sig_regs_fp,
										char* training_target_geno_sig_regs_fp,
										char* training_sample_ids_list_fp,
										char* testing_tag_geno_sig_regs_fp,
										char* testing_target_geno_sig_regs_fp,
										char* testing_sample_ids_list_fp,
										char* op_prefix)
{
	vector<t_annot_region*>* roi_target_var_regs = load_BED(roi_target_vars_BED_fp);
	fprintf(stderr, "Loaded %d ROI target variants regions.\n", roi_target_var_regs->size());

	// Load the testing data.
	fprintf(stderr, "Loading testing datasets.\n");
	vector<char*>* testing_sample_ids = buffer_file(testing_sample_ids_list_fp);
	fprintf(stderr, "Loading tag genotype matrix with %d testing samples.\n", testing_sample_ids->size());
	vector<t_annot_region*>* all_testing_tag_genotype_signal_regs = load_variant_genotype_signal_regions(testing_tag_geno_sig_regs_fp, testing_sample_ids);
	fprintf(stderr, "Loaded %d tag genotype signal regions with %d testing samples.\n", all_testing_tag_genotype_signal_regs->size(), testing_sample_ids->size());
	t_restr_annot_region_list* restr_testing_tag_genotype_signal_regs = restructure_annot_regions(all_testing_tag_genotype_signal_regs);

	fprintf(stderr, "Loading target genotype matrix with %d testing samples.\n", testing_sample_ids->size());
	vector<t_annot_region*>* all_testing_target_genotype_signal_regs = load_variant_genotype_signal_regions(testing_target_geno_sig_regs_fp, testing_sample_ids);
	fprintf(stderr, "Loaded %d target genotype signal regions with %d testing samples.\n", all_testing_target_genotype_signal_regs->size(), testing_sample_ids->size());
	t_restr_annot_region_list* restr_testing_target_genotype_signal_regs = restructure_annot_regions(all_testing_target_genotype_signal_regs);

	// Load the training data.
	fprintf(stderr, "Loading training datasets.\n");
	vector<char*>* training_sample_ids = buffer_file(training_sample_ids_list_fp);

	fprintf(stderr, "Loading tag genotype matrix with %d training samples.\n", training_sample_ids->size());
	vector<t_annot_region*>* all_training_tag_genotype_signal_regs = load_variant_genotype_signal_regions(training_tag_geno_sig_regs_fp, training_sample_ids);
	fprintf(stderr, "Loaded %d tag genotype signal regions with %d training samples.\n", all_training_tag_genotype_signal_regs->size(), training_sample_ids->size());
	t_restr_annot_region_list* restr_training_tag_genotype_signal_regs = restructure_annot_regions(all_training_tag_genotype_signal_regs);

	fprintf(stderr, "Loading target genotype matrix with %d training samples.\n", training_sample_ids->size());
	vector<t_annot_region*>* all_training_target_genotype_signal_regs = load_variant_genotype_signal_regions(training_target_geno_sig_regs_fp, training_sample_ids);
	fprintf(stderr, "Loaded %d target genotype signal regions with %d training samples.\n", all_training_target_genotype_signal_regs->size(), training_sample_ids->size());
	t_restr_annot_region_list* restr_training_target_genotype_signal_regs = restructure_annot_regions(all_training_target_genotype_signal_regs);

	for (int i_roi = 0; i_roi < roi_target_var_regs->size(); i_roi++)
	{
		fprintf(stderr, "Saving matrices for %s:%d (%s)\n", roi_target_var_regs->at(i_roi)->chrom, roi_target_var_regs->at(i_roi)->start, roi_target_var_regs->at(i_roi)->name);

		int i_chr = t_string::get_i_str(restr_training_target_genotype_signal_regs->chr_ids, roi_target_var_regs->at(i_roi)->chrom);

		if (i_chr != t_string::get_i_str(restr_testing_target_genotype_signal_regs->chr_ids, roi_target_var_regs->at(i_roi)->chrom) ||
			i_chr != t_string::get_i_str(restr_training_tag_genotype_signal_regs->chr_ids, roi_target_var_regs->at(i_roi)->chrom) ||
			i_chr != t_string::get_i_str(restr_testing_tag_genotype_signal_regs->chr_ids, roi_target_var_regs->at(i_roi)->chrom))
		{
			fprintf(stderr, "Sanity check failed @ chromosome checks for %s\n", roi_target_var_regs->at(i_roi)->chrom);
			exit(0);
		}

		vector<t_annot_region*>* training_target_genotype_signal_regs = restr_training_target_genotype_signal_regs->regions_per_chrom[i_chr];
		vector<t_annot_region*>* testing_target_genotype_signal_regs = restr_testing_target_genotype_signal_regs->regions_per_chrom[i_chr];
		vector<t_annot_region*>* training_tag_genotype_signal_regs = restr_training_tag_genotype_signal_regs->regions_per_chrom[i_chr];
		vector<t_annot_region*>* testing_tag_genotype_signal_regs = restr_testing_tag_genotype_signal_regs->regions_per_chrom[i_chr];

		// Find the training data for the current VOI.
		for (int i_var = 0; i_var < training_target_genotype_signal_regs->size(); i_var++)
		{
			if (training_target_genotype_signal_regs->at(i_var)->start == roi_target_var_regs->at(i_roi)->start)
			{
				fprintf(stderr, "Found training target region.\n");

				// Write.
				char target_training_geno_fp[1000];
				sprintf(target_training_geno_fp, "%s_%s_%d_training_target.txt", 
						op_prefix, 
						roi_target_var_regs->at(i_roi)->chrom, roi_target_var_regs->at(i_roi)->start);

				void** cur_var_info = (void**)(training_target_genotype_signal_regs->at(i_var)->data);
				char* cur_var_geno_signal = (char*)(cur_var_info[0]);

				FILE* f_target_training = open_f(target_training_geno_fp, "w");
				fprintf(f_target_training, "%s", roi_target_var_regs->at(i_roi)->name);
				for (int i_s = 0; i_s < training_sample_ids->size(); i_s++)
				{
					fprintf(f_target_training, "\t%d", (int)(cur_var_geno_signal[i_s]));
				} // i_s loop.
				fclose(f_target_training);
				break;
			} // coordinate comparison.
		} // i_train_target loop.

		// Search for the testing target region.
		for (int i_var = 0; i_var < testing_target_genotype_signal_regs->size(); i_var++)
		{
			if (testing_target_genotype_signal_regs->at(i_var)->start == roi_target_var_regs->at(i_roi)->start)
			{
				fprintf(stderr, "Found testing target region.\n");

				// Write.
				char target_testing_geno_fp[1000];
				sprintf(target_testing_geno_fp, "%s_%s_%d_testing_target.txt",
						op_prefix,
						roi_target_var_regs->at(i_roi)->chrom, roi_target_var_regs->at(i_roi)->start);

				void** cur_var_info = (void**)(testing_target_genotype_signal_regs->at(i_var)->data);
				char* cur_var_geno_signal = (char*)(cur_var_info[0]);

				FILE* f_target_testing = open_f(target_testing_geno_fp, "w");
				fprintf(f_target_testing, "%s", roi_target_var_regs->at(i_roi)->name);
				for (int i_s = 0; i_s < testing_sample_ids->size(); i_s++)
				{
					fprintf(f_target_testing, "\t%d", (int)(cur_var_geno_signal[i_s]));
				} // i_s loop.
				fclose(f_target_testing);
				break;
			} // coordinate comparison.
		} // i_var loop.

		// Search for the training tags.
		vector<t_annot_region*>* cur_var_training_tag_regs = new vector<t_annot_region*>();
		for (int i_var = 0; i_var < training_tag_genotype_signal_regs->size()-1; i_var++)
		{
			if (training_tag_genotype_signal_regs->at(i_var)->start < roi_target_var_regs->at(i_roi)->start &&
				training_tag_genotype_signal_regs->at(i_var+1)->start > roi_target_var_regs->at(i_roi)->start)
			{
				fprintf(stderr, "Found training tag region.\n");

				int min_j_var = MAX(0, i_var - n_half_vicinity_tags + 1);
				int max_j_var = MIN(training_tag_genotype_signal_regs->size()-1, i_var + n_half_vicinity_tags);

				// Write.
				char training_tag_geno_fp[1000];
				sprintf(training_tag_geno_fp, "%s_%s_%d_training_tags.txt",
						op_prefix,
						roi_target_var_regs->at(i_roi)->chrom, roi_target_var_regs->at(i_roi)->start);

				FILE* f_training_tag = open_f(training_tag_geno_fp, "w");
				for (int j_var = min_j_var; j_var <= max_j_var; j_var++)
				{
					cur_var_training_tag_regs->push_back(training_tag_genotype_signal_regs->at(j_var));

					void** cur_var_info = (void**)(training_tag_genotype_signal_regs->at(j_var)->data);
					char* cur_var_geno_signal = (char*)(cur_var_info[0]);

					fprintf(f_training_tag, "%d", training_tag_genotype_signal_regs->at(j_var)->start);
					for (int i_s = 0; i_s < training_sample_ids->size(); i_s++)
					{
						fprintf(f_training_tag, "\t%d", (int)(cur_var_geno_signal[i_s]));
					} // i_s loop.

					fprintf(f_training_tag, "\n");
				} // j_var loop.
				
				fclose(f_training_tag);
				break;
			} // coordinate comparison.
		} // i_var loop.

		// Search for the testing tags.
		vector<t_annot_region*>* cur_var_testing_tag_regs = new vector<t_annot_region*>();
		for (int i_var = 0; i_var < testing_tag_genotype_signal_regs->size() - 1; i_var++)
		{
			if (testing_tag_genotype_signal_regs->at(i_var)->start < roi_target_var_regs->at(i_roi)->start &&
				testing_tag_genotype_signal_regs->at(i_var + 1)->start > roi_target_var_regs->at(i_roi)->start)
			{
				fprintf(stderr, "Found testing tag region.\n");

				int min_j_var = MAX(0, i_var - n_half_vicinity_tags + 1);
				int max_j_var = MIN(testing_tag_genotype_signal_regs->size() - 1, i_var + n_half_vicinity_tags);

				// Write.
				char testing_tag_geno_fp[1000];
				sprintf(testing_tag_geno_fp, "%s_%s_%d_testing_tags.txt",
						op_prefix,
						roi_target_var_regs->at(i_roi)->chrom, roi_target_var_regs->at(i_roi)->start);

				FILE* f_testing_tag = open_f(testing_tag_geno_fp, "w");
				for (int j_var = min_j_var; j_var <= max_j_var; j_var++)
				{
					cur_var_testing_tag_regs->push_back(testing_tag_genotype_signal_regs->at(j_var));

					void** cur_var_info = (void**)(testing_tag_genotype_signal_regs->at(j_var)->data);
					char* cur_var_geno_signal = (char*)(cur_var_info[0]);

					fprintf(f_testing_tag, "%d", testing_tag_genotype_signal_regs->at(j_var)->start);
					for (int i_s = 0; i_s < testing_sample_ids->size(); i_s++)
					{
						fprintf(f_testing_tag, "\t%d", (int)(cur_var_geno_signal[i_s]));
					} // i_s loop.

					fprintf(f_testing_tag, "\n");
				} // j_var loop.
				
				fclose(f_testing_tag);
				break;
			} // coordinate comparison.
		} // i_var loop.

		// Compare the tag variants for training and testing sets.
		if (cur_var_testing_tag_regs->size() == 0 ||
			cur_var_testing_tag_regs->size() != cur_var_training_tag_regs->size())
		{
			fprintf(stderr, "Sanity check failed @ # of tags: %d, %d\n", cur_var_testing_tag_regs->size(), cur_var_training_tag_regs->size());
			exit(0);
		}
		else
		{
			for (int i_tag = 0; i_tag < cur_var_testing_tag_regs->size(); i_tag++)
			{
				if (cur_var_testing_tag_regs->at(i_tag)->start != cur_var_training_tag_regs->at(i_tag)->start)
				{
					fprintf(stderr, "%d. tag is different: %s:%d vs %s:%d\n", i_tag,
							cur_var_testing_tag_regs->at(i_tag)->name, cur_var_testing_tag_regs->at(i_tag)->start,
							cur_var_training_tag_regs->at(i_tag)->name, cur_var_training_tag_regs->at(i_tag)->start);

					exit(0);
				}
			} // i_tag loop.

			fprintf(stderr, "Found matching %d tags.\n", cur_var_training_tag_regs->size());
		}

		delete cur_var_testing_tag_regs;
		delete cur_var_training_tag_regs;
	} // i_roi loop.
}

static void* LMSE_imputation_thread_callback(void* thread_info_ptr)
{
	t_LMSE_imputation_thread_info* impute_thread_info = (t_LMSE_imputation_thread_info*)(thread_info_ptr);

	// Following are the thread info.
	double target_normalizer = impute_thread_info->target_normalizer;
	double tag_normalizer = impute_thread_info->tag_normalizer;

	//vector<t_annot_region*>* testing_tag_genotype_signal_regs = impute_thread_info->testing_tag_genotype_signal_regs;
	vector<t_annot_region*>* tag_genotype_signal_regs = impute_thread_info->tag_genotype_signal_regs;
	//vector<t_annot_region*>* testing_target_genotype_signal_regs = impute_thread_info->testing_target_genotype_signal_regs;
	vector<t_annot_region*>* target_genotype_signal_regs = impute_thread_info->target_genotype_signal_regs;

	t_restr_annot_region_list* restr_target_genotype_signal_regs = impute_thread_info->restr_target_genotype_signal_regs;
	t_restr_annot_region_list* restr_tag_genotype_signal_regs = impute_thread_info->restr_tag_genotype_signal_regs;

	int testing_flag = impute_thread_info->testing_flag;
	int start_pos = impute_thread_info->start_pos;
	int end_pos = impute_thread_info->end_pos;
	char* op_dir = impute_thread_info->op_dir;

	vector<char*> * tag_sample_ids = impute_thread_info->tag_sample_ids;
	vector<char*>* target_sample_ids = impute_thread_info->target_sample_ids;

	vector<char*>* testing_tag_sample_ids = impute_thread_info->testing_tag_sample_ids;
	vector<char*>* testing_target_sample_ids = impute_thread_info->testing_target_sample_ids;

	int n_tags_vars_per_side = impute_thread_info->n_tags_vars_per_side;
	int which = impute_thread_info->which;
	int outof = impute_thread_info->outof;

	fprintf(stderr, "Started Imputation thread %d/%d\n", which, outof);

	// This is the number of dimensions.
	int n_params_w_intercept = 2 * n_tags_vars_per_side + 1;

	// This is the number of data points in training.
	int training_sample_size = target_sample_ids->size();

	// Allocate the fitting data.
	double xi, yi, ei, chisq;
	gsl_matrix *X, *cov;
	gsl_vector *y, *w, *c;
	gsl_vector *testing_x;

	// SVD associated matrices/vectors.
	gsl_matrix* fullX = gsl_matrix_alloc(training_sample_size, n_params_w_intercept);
	gsl_matrix * svdV = gsl_matrix_alloc(n_params_w_intercept, n_params_w_intercept);
	gsl_vector * svdS = gsl_vector_alloc(n_params_w_intercept);
	gsl_vector * svdWork = gsl_vector_alloc(n_params_w_intercept);

	// Tag usage flag at the variant coordinates; this is necessary to save file and write per variant messages.
	int* tag_usage_flag = new int[n_params_w_intercept + 1];
	memset(tag_usage_flag, 0, sizeof(int) * (n_params_w_intercept + 1));

	char per_var_stats_fp[1000];
	sprintf(per_var_stats_fp, "per_var_stats_thread_%d.txt", which);
	FILE* f_per_var_stats = open_f(per_var_stats_fp, "w");

	// Process all the target variants.
	for (int target_chr_i = 0;
		target_chr_i < restr_target_genotype_signal_regs->chr_ids->size();
		target_chr_i++)
	{
		fprintf(stderr, "Building the models on chromosome %s (%d/%d)\n", restr_target_genotype_signal_regs->chr_ids->at(target_chr_i), which, outof);

		// Get the chromosome index for testring variant regions.
		int tag_chr_i = t_string::get_i_str(restr_tag_genotype_signal_regs->chr_ids, restr_target_genotype_signal_regs->chr_ids->at(target_chr_i));

		vector<t_annot_region*>* target_var_regs = restr_target_genotype_signal_regs->regions_per_chrom[target_chr_i];
		vector<t_annot_region*>* tag_var_regs = restr_tag_genotype_signal_regs->regions_per_chrom[tag_chr_i];

		// Following sets up n data points. We set this up over the samples.
		for (int target_var_i = 0; target_var_i < target_var_regs->size(); target_var_i++)
		{
			if (start_pos < 0 ||
				end_pos < 0 ||
				target_var_regs->at(target_var_i)->start >= start_pos && target_var_regs->at(target_var_i)->start <= end_pos)
			{

			}
			else
			{
				continue;
			}

			if (target_var_i % outof != which)
			{
				continue;
			}

			int closest_tag_var_left_i = -1;
			for (int tag_var_i = 0; (tag_var_i + 1)< tag_var_regs->size(); tag_var_i++)
			{
				if (tag_var_regs->at(tag_var_i)->start < target_var_regs->at(target_var_i)->start &&
					tag_var_regs->at(tag_var_i + 1)->start > target_var_regs->at(target_var_i)->start)
				{
					closest_tag_var_left_i = tag_var_i;
					break;
				}
			} // tag_var_i loop.

			if (closest_tag_var_left_i == -1)
			{
				continue;
			}

			if (__DUMP_INPUTATION_UTILS_MSGS__)
			{
				fprintf(stderr, "Closest tag to target %d: (%d, %d)\n", target_var_regs->at(target_var_i)->start,
					tag_var_regs->at(closest_tag_var_left_i)->start, tag_var_regs->at(closest_tag_var_left_i + 1)->start);
			}

			// Add the tag variant to the right.
			vector<t_annot_region*>* tag_var_regs_per_cur_target_var = new vector<t_annot_region*>();
			for (int tag_var_i = MAX(0, closest_tag_var_left_i - n_tags_vars_per_side + 1);
				tag_var_i <= closest_tag_var_left_i;
				tag_var_i++)
			{
				tag_var_regs_per_cur_target_var->push_back(tag_var_regs->at(tag_var_i));
			} // tag_var_i loop.

			for (int tag_var_i = closest_tag_var_left_i + 1;
				tag_var_i <= (closest_tag_var_left_i + n_tags_vars_per_side) && tag_var_i < tag_var_regs->size();
				tag_var_i++)
			{
				tag_var_regs_per_cur_target_var->push_back(tag_var_regs->at(tag_var_i));
			} // tag_var_i loop.

			  // Check the number of tag SNVs around the target SNVs, make sure they satisfy the requested number.
			if (tag_var_regs_per_cur_target_var->size() != 2 * n_tags_vars_per_side)
			{
				if (__DUMP_INPUTATION_UTILS_MSGS__)
				{
					fprintf(stderr, "Could not find the correct number of tag variants for %s:%d (%d), will skip.\n",
						target_var_regs->at(target_var_i)->chrom, target_var_regs->at(target_var_i)->start,
						tag_var_regs_per_cur_target_var->size());
				}
			}
			else
			{
				// Set the tag usage flag.
				memset(tag_usage_flag, 0, sizeof(int) * n_params_w_intercept);
				int n_tags_2_use = 0;
				for (int par_i = 0; par_i < n_params_w_intercept; par_i++)
				{
					tag_usage_flag[par_i] = 1;
				} // par_i loop.

				  // Fill fullX.
				for (int sample_i = 0; sample_i < target_sample_ids->size(); sample_i++)
				{
					// Add the intercept; 1.0
					gsl_matrix_set(fullX, sample_i, 0, 1.0);

					// Jump over the intercept, set the tag genotype matrix entry.
					for (int tag_dim_i = 1; tag_dim_i <= 2 * n_tags_vars_per_side; tag_dim_i++)
					{
						void** reg_data = (void**)(tag_var_regs_per_cur_target_var->at(tag_dim_i - 1)->data);
						char* cur_tag_genotypes = (char*)(reg_data[0]);

						// Set the current entry for the full X matrix.
						double tagg_ij = cur_tag_genotypes[sample_i];

						// Normalize: Use tag normalization factor.
						tagg_ij /= tag_normalizer;

						gsl_matrix_set(fullX, sample_i, tag_dim_i, tagg_ij);
					} // tag_dim_i loop.
				} // sample_i loop.			

				  // https://stackoverflow.com/questions/12304963/using-eigenvalues-to-test-for-singularity-identifying-collinear-columns
				  // Generate the SVD of X.
				gsl_linalg_SV_decomp(fullX, svdV, svdS, svdWork);

				// Look for the low singular values that can create problems.
				double top_singular_value = gsl_vector_get(svdS, 0);
				for (int dim_i = 0; dim_i < n_params_w_intercept; dim_i++)
				{
					// Evaluate the current singular value.
					double cur_singular_val = gsl_vector_get(svdS, dim_i);
					//fprintf(stderr, "%d: %d: %lf\n", target_var_regs->at(target_var_i)->start, i, cur_singular_val);

					double sv_tolerance_fraction = 0.001;
					double min_component_strength = 0.01;

					// If the singular value is 0, identify which variants have components on the eigenvectors with 0 singular:
					if (cur_singular_val < top_singular_value * sv_tolerance_fraction &&
						cur_singular_val < 0.001)
					{
						if (__DUMP_INPUTATION_UTILS_MSGS__)
						{
							fprintf(stderr, "Negligible singular value: %lf (dim: %d)\n", cur_singular_val, dim_i);
						}

						// Check the entries that are on this eigenvector.
						vector<int>* colinear_var_par_i = new vector<int>();
						for (int par_i = 1; par_i < n_params_w_intercept; par_i++)
						{
							double cur_comp = gsl_matrix_get(svdV, par_i, dim_i);

							if (fabs(cur_comp) > min_component_strength)
							{
								if (__DUMP_INPUTATION_UTILS_MSGS__)
								{
									fprintf(stderr, "%d has component on 0 eigenvalue: %lf\n", tag_var_regs_per_cur_target_var->at(par_i - 1)->start, cur_comp);
								}

								colinear_var_par_i->push_back(par_i - 1);

								tag_usage_flag[par_i - 1] = 0;
							}
						} // par_i loop.

						  //// Following performs variant selection to optimize the selections.
						  //int n_colinear_vars = colinear_var_par_i->size();
						  //int tag_i_2_set_unused = n_colinear_vars - 1;
						  //while (tag_i_2_set_unused >= 0)
						  //{
						  //	if (tag_usage_flag[colinear_var_par_i->at(tag_i_2_set_unused)] == 0)
						  //	{
						  //		tag_i_2_set_unused--;
						  //	}
						  //	else
						  //	{
						  //		tag_usage_flag[colinear_var_par_i->at(tag_i_2_set_unused)] = 0;
						  //		fprintf(stderr, "%d colinear tags; setting %d to be unused: %lf\n",
						  //				n_colinear_vars,
						  //				tag_var_regs_per_cur_target_var->at(colinear_var_par_i->at(tag_i_2_set_unused))->start);
						  //		break;
						  //	}
						  //}
						delete colinear_var_par_i;
					}
				} // i loop.

				n_tags_2_use = 0;
				for (int par_i = 0; par_i < n_params_w_intercept - 1; par_i++)
				{
					if (tag_usage_flag[par_i] == 1)
					{
						n_tags_2_use++;
					}
				} // par_i loop.

				if (__DUMP_INPUTATION_UTILS_MSGS__)
				{
					fprintf(stderr, "%s:%d (%d tags)\n", target_var_regs->at(target_var_i)->chrom, target_var_regs->at(target_var_i)->start, n_tags_2_use);
				}

				int n_tags_2_use_w_intercept = n_tags_2_use + 1;

				// Allocate the structures.
				X = gsl_matrix_alloc(training_sample_size, n_tags_2_use_w_intercept);
				testing_x = gsl_vector_alloc(n_tags_2_use_w_intercept);

				gsl_matrix * svdX = gsl_matrix_alloc(training_sample_size, n_tags_2_use_w_intercept);
				gsl_matrix * V = gsl_matrix_alloc(n_tags_2_use_w_intercept, n_tags_2_use_w_intercept);
				gsl_vector * S = gsl_vector_alloc(n_tags_2_use_w_intercept);
				gsl_vector * svd_work = gsl_vector_alloc(n_tags_2_use_w_intercept);

				// Target SNP genotypes vector.
				y = gsl_vector_alloc(training_sample_size);

				// Weights of each data point.
				w = gsl_vector_alloc(training_sample_size);

				// This is the coefficients matrix.
				c = gsl_vector_alloc(n_tags_2_use_w_intercept);
				cov = gsl_matrix_alloc(n_tags_2_use_w_intercept, n_tags_2_use_w_intercept);

				// Set the matrices and vectors for regression.
				for (int sample_i = 0; sample_i < target_sample_ids->size(); sample_i++)
				{
					// Add the intercept; 1.0
					gsl_matrix_set(X, sample_i, 0, 1.0);

					// Jump over the intercept, set the tag genotype matrix entry.
					int matrix_dim_i = 1;
					for (int tag_dim_i = 1; tag_dim_i <= 2 * n_tags_vars_per_side; tag_dim_i++)
					{
						if (tag_usage_flag[tag_dim_i - 1] == 1)
						{
							void** reg_data = (void**)(tag_var_regs_per_cur_target_var->at(tag_dim_i - 1)->data);
							char* cur_tag_genotypes = (char*)(reg_data[0]);

							// Set the current X matrix value.
							double tagg_ij = cur_tag_genotypes[sample_i];

							// Normalize with tag normalizer.
							tagg_ij /= tag_normalizer;

							//gsl_matrix_set(X, sample_i, tag_dim_i, tagg_ij);
							gsl_matrix_set(X, sample_i, matrix_dim_i, tagg_ij);
							matrix_dim_i++;
						}
					} // tag_dim_i loop.

					if (matrix_dim_i != n_tags_2_use_w_intercept)
					{
						fprintf(stderr, "Could not match the matrix dimension to the tags 2 use: %d, %d\n", matrix_dim_i, n_tags_2_use_w_intercept);
						exit(0);
					}

					// Set the target genotype vector entry.
					void** reg_data = (void**)(target_var_regs->at(target_var_i)->data);
					char* target_var_genotypes = (char*)(reg_data[0]);
					double targetg_ij = target_var_genotypes[sample_i];

					// Normalize with target normalizer.
					targetg_ij /= target_normalizer;

					gsl_vector_set(y, sample_i, targetg_ij);
				} // sample_i loop.				

				  // At this point, the target genotype array and the tag genotype matrix are setup.
				  // Now do the fitting.
				  // Do fitting.
				{
					// Do the fit.
					gsl_multifit_linear_workspace* work = gsl_multifit_linear_alloc(training_sample_size, n_tags_2_use_w_intercept);
					//gsl_multifit_linear_workspace* work = gsl_multifit_linear_alloc(training_sample_size, n_params_w_intercept);
					//gsl_multifit_linear_tsvd(X, y, c, cov, &chisq, &sv_rank, work); // Does not link!
					gsl_multifit_linear(X, y, c, cov, &chisq, work);
					gsl_multifit_linear_free(work);
				} // fitting block.

				  // Save the current parameters.
				  //#define C(i) (gsl_vector_get(c,(i)))
#define COV(i,j) (gsl_matrix_get(cov,(i),(j)))

			
				if (__DUMP_INPUTATION_UTILS_MSGS__)
				{
					fprintf(stderr, "Saving Parameters for %s:%d (%d-%d): %d tags used.          \r",
						target_var_regs->at(target_var_i)->chrom,
						target_var_regs->at(target_var_i)->start,
						tag_var_regs->at(closest_tag_var_left_i)->start, tag_var_regs->at(closest_tag_var_left_i + 1)->start, n_tags_2_use_w_intercept);

					int used_par_i = 0;
					for (int par_i = 0; par_i < n_params_w_intercept; par_i++)
					{
						if (par_i == 0)
						{
							fprintf(stderr, "Intercept\t%.17f\n", gsl_vector_get(c, used_par_i));
							used_par_i++;
						}
						else if (tag_usage_flag[par_i - 1] == 1)
						{
							fprintf(stderr, "%d\t%.17f\n", tag_var_regs_per_cur_target_var->at(par_i - 1)->start, gsl_vector_get(c, used_par_i));
							used_par_i++;
						}
						else
						{
							fprintf(stderr, "%d\tNA\n", tag_var_regs_per_cur_target_var->at(par_i - 1)->start);
						}
					} // par_i loop.
				}

				// Save the parameters file.
				char op_fp[1000];
				sprintf(op_fp, "%s/%d.params", op_dir, target_var_regs->at(target_var_i)->start);
				FILE* f_op = open_f(op_fp, "w");
				int used_par_i = 0;
				for (int par_i = 0; par_i < n_params_w_intercept; par_i++)
				{
					if (par_i == 0)
					{
						fprintf(f_op, "%.17f\n", gsl_vector_get(c, used_par_i));
						used_par_i++;
					}
					else if (tag_usage_flag[par_i - 1] == 1)
					{
						fprintf(f_op, "%.17f\n", gsl_vector_get(c, used_par_i));
						used_par_i++;
					}
					else
					{
						fprintf(f_op, "NA\n");
					}
				} // par_i loop.

				for (int tag_i = 0; tag_i < tag_var_regs_per_cur_target_var->size(); tag_i++)
				{
					fprintf(f_op, "%d\n", tag_var_regs_per_cur_target_var->at(tag_i)->start);
				} // par_i loop.

				fprintf(f_op, "%d\n", target_var_regs->at(target_var_i)->start);
				fclose(f_op);

				// Test if available.
				if (testing_flag)
				{
					vector<double>* per_testing_sample_errors = new vector<double>();
					vector<double>* per_testing_sample_non_ref_errors = new vector<double>();
					for (int testing_sample_i = 0; testing_sample_i < testing_tag_sample_ids->size(); testing_sample_i++)
					{
						// Reset the testing vector x: First element of x is 2.
						gsl_vector_set(testing_x, 0, 1.0);

						int matrix_dim_i = 1;
						for (int tag_dim_i = 1; tag_dim_i <= 2 * n_tags_vars_per_side; tag_dim_i++)
						{
							if (tag_usage_flag[tag_dim_i - 1] == 1)
							{
								void** cur_tag_reg_data = (void**)(tag_var_regs_per_cur_target_var->at(tag_dim_i - 1)->data);

								char* cur_testing_tag_genotypes = (char*)(cur_tag_reg_data[1]);
								if (cur_testing_tag_genotypes == NULL)
								{
									fprintf(stderr, "Testing Target region is null, exiting.\n");
									exit(0);
								}

								double cur_testing_sample_geno = (double)(cur_testing_tag_genotypes[testing_sample_i]);

								// Normalize with tag normalizer.
								cur_testing_sample_geno /= tag_normalizer;

								gsl_vector_set(testing_x, matrix_dim_i, cur_testing_sample_geno);
								matrix_dim_i++;
							} // tag usage flag check.
						} // par_i loop.

						if (matrix_dim_i != n_tags_2_use_w_intercept)
						{
							fprintf(stderr, "Could not match the matrix dimensions: %d, %d\n",
								matrix_dim_i, n_tags_2_use_w_intercept);
							exit(0);
						}

						double geno_est = 0;
						double geno_err_est = 0;
						gsl_multifit_linear_est(testing_x, c, cov, &geno_est, &geno_err_est);

						// Estimate the real error.
						void** target_reg_data = (void**)(target_var_regs->at(target_var_i)->data);
						char* testing_target_var_genotypes = (char*)(target_reg_data[1]);
						if (testing_target_var_genotypes == NULL)
						{
							fprintf(stderr, "Testing Target region is null, exiting.\n");
							exit(0);
						}

						// Set the predicted genotype.
						char* pred_genoypes = (char*)(target_reg_data[2]);
						pred_genoypes[testing_sample_i] = (char)(100 * geno_est);

						double cur_testing_sample_target_geno = (double)(testing_target_var_genotypes[testing_sample_i]);

						// Normalize.
						cur_testing_sample_target_geno /= target_normalizer;

						double real_err = fabs(cur_testing_sample_target_geno - geno_est);

						if (__DUMP_INPUTATION_UTILS_MSGS__)
						{
							fprintf(stderr, "%s:%d-%d (%s): Error: %lf (%lf)\n",
								target_var_regs->at(target_var_i)->chrom,
								target_var_regs->at(target_var_i)->start,
								target_var_regs->at(target_var_i)->end,
								target_var_regs->at(target_var_i)->name,
								real_err,
								geno_err_est);
						}

						per_testing_sample_errors->push_back(real_err);

						if (__DUMP_INPUTATION_UTILS_MSGS__)
						{
							fprintf(stderr, "%s:%d: Sample %d: %.4f\n", target_var_regs->at(target_var_i)->chrom, target_var_regs->at(target_var_i)->start, testing_sample_i, real_err);
						}

						if (cur_testing_sample_target_geno != 0)
						{
							per_testing_sample_non_ref_errors->push_back(real_err);
						}
					} // testing_sample_i loop.

					  // Get the error statistics.
					double mean_err = 0;
					double std_dev_err = 0;
					get_stats(per_testing_sample_errors, mean_err, std_dev_err);

					double mean_non_ref_err = 0;
					double stddev_non_ref_err = 0;
					get_stats(per_testing_sample_non_ref_errors, mean_non_ref_err, stddev_non_ref_err);

					if (__DUMP_INPUTATION_UTILS_MSGS__)
					{
						fprintf(stderr, "%s:%d-%d (%s): Avg Error: %lf (%d); Avg Non-ref Error: %lf (%d)          \n",
							target_var_regs->at(target_var_i)->chrom,
							target_var_regs->at(target_var_i)->start,
							target_var_regs->at(target_var_i)->end,
							target_var_regs->at(target_var_i)->name,
							mean_err,
							per_testing_sample_errors->size(),
							mean_non_ref_err,
							per_testing_sample_non_ref_errors->size());
					}

					fprintf(f_per_var_stats, "%s\t%d\t%d\t%s\t%lf\t%d\t%lf\t%d\t%d\n",
						target_var_regs->at(target_var_i)->chrom,
						target_var_regs->at(target_var_i)->start,
						target_var_regs->at(target_var_i)->end,
						target_var_regs->at(target_var_i)->name,
						mean_err,
						per_testing_sample_errors->size(),
						mean_non_ref_err,
						per_testing_sample_non_ref_errors->size(),
						n_tags_2_use_w_intercept);

					delete per_testing_sample_non_ref_errors;
					delete per_testing_sample_errors;
				} // testing check.

				  // Free memory.
				gsl_matrix_free(X);
				gsl_vector_free(testing_x);
				gsl_matrix_free(svdX);
				gsl_matrix_free(V);
				gsl_vector_free(S);
				gsl_vector_free(svd_work);
				gsl_vector_free(y);
				gsl_vector_free(w);
				gsl_vector_free(c);
				gsl_matrix_free(cov);

				//printf("# covariance matrix:\n");
				//printf("[ %+.5e, %+.5e, %+.5e  \n",
				//	COV(0, 0), COV(0, 1), COV(0, 2));
				//printf("  %+.5e, %+.5e, %+.5e  \n",
				//	COV(1, 0), COV(1, 1), COV(1, 2));
				//printf("  %+.5e, %+.5e, %+.5e ]\n",
				//	COV(2, 0), COV(2, 1), COV(2, 2));
				//printf("# chisq = %g\n", chisq);
			} // tag snp count check.

			delete tag_var_regs_per_cur_target_var;
		} // target_var_i loop.
	} // target_chr_i loop.

	fclose(f_per_var_stats);
}

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

	char* op_dir)
{
	fprintf(stderr, "Building LMSE models using matrices %s, %s with sample id's in %s, %s for SNVs in range [%d-%d] using %d threads.\n",
		tag_genotype_matrix_fp, target_genotype_matrix_fp,
		tag_sample_ids_list_fp, target_sample_ids_list_fp,
		start_pos, end_pos, n_threads);

	vector<char*>* tag_sample_ids = buffer_file(tag_sample_ids_list_fp);
	fprintf(stderr, "Loading tag genotype matrix with %d training samples.\n", tag_sample_ids->size());
	vector<t_annot_region*>* tag_genotype_signal_regs = load_variant_genotype_signal_regions(tag_genotype_matrix_fp, tag_sample_ids);
	fprintf(stderr, "Loaded %d tag genotype signal regions with %d training samples.\n", tag_genotype_signal_regs->size(), tag_sample_ids->size());
	t_restr_annot_region_list* restr_tag_genotype_signal_regs = restructure_annot_regions(tag_genotype_signal_regs);

	vector<char*>* target_sample_ids = buffer_file(target_sample_ids_list_fp);
	fprintf(stderr, "Loading target genotype matrix with %d training samples.\n", target_sample_ids->size());
	vector<t_annot_region*>* target_genotype_signal_regs = load_variant_genotype_signal_regions(target_genotype_matrix_fp, target_sample_ids);
	fprintf(stderr, "Loaded %d target genotype signal regions with %d training samples.\n", target_genotype_signal_regs->size(), target_sample_ids->size());
	t_restr_annot_region_list* restr_target_genotype_signal_regs = restructure_annot_regions(target_genotype_signal_regs);

	if (target_sample_ids->size() != tag_sample_ids->size())
	{
		fprintf(stderr, "Target and tag sample id's are not the same size: %d, %d\n", target_sample_ids->size(), tag_sample_ids->size());
		exit(0);
	}

	// Make sure the sample id's are matching.
	for (int i_s = 0; i_s < target_sample_ids->size(); i_s++)
	{
		if (!t_string::compare_strings(target_sample_ids->at(i_s), tag_sample_ids->at(i_s)))
		{
			fprintf(stderr, "Target and tag sample id's are not matching: %s, %s; Re-order and run again.\n", target_sample_ids->at(i_s), tag_sample_ids->at(i_s));
			exit(0);
		}
	} // i_s loop.

	  // Decide on haplocoded signal; tag and target normalizers.
	double tag_normalizer = -1;
	double target_normalizer = -1;
	for (int var_i = 0; var_i < tag_genotype_signal_regs->size(); var_i++)
	{
		// Get the maximum.
		void** cur_var_sig_info = (void**)(tag_genotype_signal_regs->at(var_i)->data);
		char* cur_var_sig = (char*)(cur_var_sig_info[0]);

		for (int i_s = 0; i_s < target_sample_ids->size(); i_s++)
		{
			tag_normalizer = MAX(cur_var_sig[i_s], tag_normalizer);
		} // i_s loop.
	} // var_i loop.

	for (int var_i = 0; var_i < target_genotype_signal_regs->size(); var_i++)
	{
		// Get the maximum.
		void** cur_var_sig_info = (void**)(target_genotype_signal_regs->at(var_i)->data);
		char* cur_var_sig = (char*)(cur_var_sig_info[0]);

		for (int i_s = 0; i_s < target_sample_ids->size(); i_s++)
		{
			target_normalizer = MAX(cur_var_sig[i_s], target_normalizer);
		} // i_s loop.
	} // var_i loop.

	bool testing_flag = false;
	vector<char*>* testing_tag_sample_ids = NULL;
	vector<char*>* testing_target_sample_ids = NULL;
	if (check_file(testing_tag_genotype_matrix_fp) &&
		check_file(testing_target_genotype_matrix_fp))
	{
		fprintf(stderr, "Testing data exists, loading and setting it.\n");
		testing_flag = true;

		fprintf(stderr, "Loading and assigning testing target data.\n");
		testing_tag_sample_ids = buffer_file(testing_tag_sample_ids_list_fp);
		fprintf(stderr, "Loading testing tag genotype matrix with %d testing samples.\n", testing_tag_sample_ids->size());
		vector<t_annot_region*>* testing_tag_genotype_signal_regs = load_variant_genotype_signal_regions(testing_tag_genotype_matrix_fp, testing_tag_sample_ids);
		fprintf(stderr, "Loaded %d testing tag genotype signal regions with %d testing samples.\n", testing_tag_genotype_signal_regs->size(), testing_tag_sample_ids->size());

		vector<t_annot_region*>* tag_intersects = intersect_annot_regions(tag_genotype_signal_regs, testing_tag_genotype_signal_regs, true);
		fprintf(stderr, "Processing %d tag intersects.\n", tag_intersects->size());
		for (int i_int = 0; i_int < tag_intersects->size(); i_int++)
		{
			t_intersect_info* int_info = (t_intersect_info*)(tag_intersects->at(i_int)->data);
			t_annot_region* tag_reg = int_info->src_reg;
			t_annot_region* testing_tag_reg = int_info->dest_reg;

			if (t_string::compare_strings(tag_reg->name, testing_tag_reg->name))
			{
				void** cur_tag_reg_info = (void**)(tag_reg->data);
				void** cur_testing_tag_reg_info = (void**)(testing_tag_reg->data);
				cur_tag_reg_info[1] = cur_testing_tag_reg_info[0];
			}

			delete int_info;
		} // i_int loop.
		delete_annot_regions(tag_intersects);

		fprintf(stderr, "Testing assignment of tag regions to all training tag regions.\n");
		for (int i_reg = 0; i_reg < tag_genotype_signal_regs->size(); i_reg++)
		{
			void** cur_tag_reg_info = (void**)(tag_genotype_signal_regs->at(i_reg)->data);

			if (cur_tag_reg_info[1] == NULL)
			{
				fprintf(stderr, "Could not assign a testing tag genotype data for %s:%d-%d (%s)\n",
					tag_genotype_signal_regs->at(i_reg)->chrom,
					tag_genotype_signal_regs->at(i_reg)->start,
					tag_genotype_signal_regs->at(i_reg)->end,
					tag_genotype_signal_regs->at(i_reg)->name);

				exit(0);
			}
		} // i_reg loop.

		fprintf(stderr, "Loading and assigning testing target data.\n");
		testing_target_sample_ids = buffer_file(testing_target_sample_ids_list_fp);
		fprintf(stderr, "Loading testing target genotype matrix with %d testing samples.\n", testing_target_sample_ids->size());
		vector<t_annot_region*>* testing_target_genotype_signal_regs = load_variant_genotype_signal_regions(testing_target_genotype_matrix_fp, testing_target_sample_ids);
		fprintf(stderr, "Loaded %d testing target genotype signal regions with %d testing samples.\n", testing_target_genotype_signal_regs->size(), testing_target_sample_ids->size());

		vector<t_annot_region*>* target_intersects = intersect_annot_regions(target_genotype_signal_regs, testing_target_genotype_signal_regs, true);
		fprintf(stderr, "Processing %d target intersects.\n", target_intersects->size());
		for (int i_int = 0; i_int < target_intersects->size(); i_int++)
		{
			t_intersect_info* int_info = (t_intersect_info*)(target_intersects->at(i_int)->data);
			t_annot_region* target_reg = int_info->src_reg;
			t_annot_region* testing_target_reg = int_info->dest_reg;

			if (t_string::compare_strings(target_reg->name, testing_target_reg->name))
			{
				void** cur_target_reg_info = (void**)(target_reg->data);
				void** cur_testing_target_reg_info = (void**)(testing_target_reg->data);
				cur_target_reg_info[1] = cur_testing_target_reg_info[0];

				void** updated_info = new void*[5];
				updated_info[0] = cur_target_reg_info[0];
				updated_info[1] = cur_target_reg_info[1];
				updated_info[2] = new char[testing_target_sample_ids->size() + 2];
				memset(updated_info[2], 0, sizeof(char)  * (testing_target_sample_ids->size() + 2));

				target_reg->data = updated_info;
			}

			delete int_info;
		} // i_int loop.
		delete_annot_regions(target_intersects);

		fprintf(stderr, "Testing assignment of target regions to all training target regions.\n");
		for (int i_reg = 0; i_reg < target_genotype_signal_regs->size(); i_reg++)
		{
			void** cur_target_reg_info = (void**)(target_genotype_signal_regs->at(i_reg)->data);

			if (cur_target_reg_info[1] == NULL)
			{
				fprintf(stderr, "Could not assign a testing target genotype data for %s:%d-%d (%s)\n",
					target_genotype_signal_regs->at(i_reg)->chrom,
					target_genotype_signal_regs->at(i_reg)->start,
					target_genotype_signal_regs->at(i_reg)->end,
					target_genotype_signal_regs->at(i_reg)->name);

				exit(0);
			}
		} // i_reg loop.

		if (testing_target_sample_ids->size() != testing_tag_sample_ids->size())
		{
			fprintf(stderr, "Testing target and tag sample id's are not the same size: %d, %d\n", testing_target_sample_ids->size(), testing_tag_sample_ids->size());
			exit(0);
		}

		// Make sure the sample id's are matching.
		for (int i_s = 0; i_s < testing_target_sample_ids->size(); i_s++)
		{
			if (!t_string::compare_strings(testing_target_sample_ids->at(i_s), testing_tag_sample_ids->at(i_s)))
			{
				fprintf(stderr, "Testing target and tag sample id's are not matching: %s, %s; Re-order and run again.\n", testing_target_sample_ids->at(i_s), testing_tag_sample_ids->at(i_s));
				exit(0);
			}
		} // i_s loop.

		  // Go over the testing data too for updating normalizers.
		for (int var_i = 0; var_i < testing_tag_genotype_signal_regs->size(); var_i++)
		{
			// Get the maximum.
			void** cur_var_sig_info = (void**)(testing_tag_genotype_signal_regs->at(var_i)->data);
			char* cur_var_sig = (char*)(cur_var_sig_info[0]);

			for (int i_s = 0; i_s < testing_tag_sample_ids->size(); i_s++)
			{
				tag_normalizer = MAX(cur_var_sig[i_s], tag_normalizer);
			} // i_s loop.
		} // var_i loop.

		for (int var_i = 0; var_i < testing_target_genotype_signal_regs->size(); var_i++)
		{
			// Get the maximum.
			void** cur_var_sig_info = (void**)(testing_target_genotype_signal_regs->at(var_i)->data);
			char* cur_var_sig = (char*)(cur_var_sig_info[0]);

			for (int i_s = 0; i_s < testing_target_sample_ids->size(); i_s++)
			{
				target_normalizer = MAX(cur_var_sig[i_s], target_normalizer);
			} // i_s loop.
		} // var_i loop.
	} // Testing data existence check.
	else
	{
		fprintf(stderr, "Not going to perform testing statistics.\n");
	}

	if (tag_normalizer < 2 || tag_normalizer > 3 ||
		target_normalizer < 2 || target_normalizer > 3)
	{
		fprintf(stderr, "tag or target normalizer is illegal or very improbably unlucky: %.1f, %.1f\n", tag_normalizer, target_normalizer);
		exit(0);
	}

	fprintf(stderr, "----------------------------------\n");
	if (tag_normalizer == 2)
	{
		fprintf(stderr, "Tag genotypes are unphased.\n");
	}
	else if (tag_normalizer == 3)
	{
		fprintf(stderr, "Tag genotypes are phased.\n");
	}

	if (target_normalizer == 2)
	{
		fprintf(stderr, "Target genotypes are unphased.\n");
	}
	else if (target_normalizer == 3)
	{
		fprintf(stderr, "Target genotypes are phased.\n");
	}
	fprintf(stderr, "----------------------------------\n");

	// This is the number of dimensions.
	int n_params_w_intercept = 2 * n_tags_vars_per_side + 1;

	// This is the number of data points in training.
	int training_sample_size = target_sample_ids->size();

	vector<t_ansi_thread*>* impute_threads = new vector<t_ansi_thread*>();
	for (int i_thread = 0; i_thread < n_threads; i_thread++)
	{
		t_LMSE_imputation_thread_info* cur_thread_info = new t_LMSE_imputation_thread_info();
		cur_thread_info->which = i_thread;
		cur_thread_info->outof = n_threads;

		cur_thread_info->target_normalizer = target_normalizer;
		cur_thread_info->tag_normalizer = tag_normalizer;

		//cur_thread_info->testing_tag_genotype_signal_regs = testing_tag_genotype_signal_regs;
		cur_thread_info->tag_genotype_signal_regs = tag_genotype_signal_regs;
		//cur_thread_info->testing_target_genotype_signal_regs = testing_target_genotype_signal_regs;
		cur_thread_info->target_genotype_signal_regs = target_genotype_signal_regs;

		cur_thread_info->restr_target_genotype_signal_regs = restr_target_genotype_signal_regs;
		cur_thread_info->restr_tag_genotype_signal_regs = restr_tag_genotype_signal_regs;

		cur_thread_info->testing_flag = testing_flag;
		cur_thread_info->start_pos = start_pos;
		cur_thread_info->end_pos = end_pos;
		cur_thread_info->op_dir = op_dir;

		cur_thread_info->tag_sample_ids = tag_sample_ids;
		cur_thread_info->target_sample_ids = target_sample_ids;

		cur_thread_info->testing_tag_sample_ids = testing_tag_sample_ids;
		cur_thread_info->testing_target_sample_ids = testing_target_sample_ids;

		cur_thread_info->n_tags_vars_per_side = n_tags_vars_per_side;

		t_ansi_thread* new_thread = new t_ansi_thread(LMSE_imputation_thread_callback, cur_thread_info);

		new_thread->run_thread();
		impute_threads->push_back(new_thread);
	} // thread_i loop.

	// Wait for the threads.
	for (int i_thread = 0; i_thread < impute_threads->size(); i_thread++)
	{
		impute_threads->at(i_thread)->wait_thread();
	} // i_thread loop.

	// Save results.
	char geno_op_fp[1000];
	sprintf(geno_op_fp, "%s/LMSE_imputed_genotypes.sig", op_dir);

	// Following sets up n data points. We set this up over the samples.
	vector<t_annot_region*>* target_regs_2_dump = new vector<t_annot_region*>();
	for (int target_var_i = 0; target_var_i < target_genotype_signal_regs->size(); target_var_i++)
	{
		void** target_reg_info = (void**)(target_genotype_signal_regs->at(target_var_i)->data);

		if (target_reg_info[1] == NULL || target_reg_info[2] == NULL)
		{

		}
		else
		{
			t_annot_region* cur_reg = duplicate_region(target_genotype_signal_regs->at(target_var_i));
			void** cur_reg_info = new void*[5];
			cur_reg_info[0] = (char*)(target_reg_info[2]);
			cur_reg->data = cur_reg_info;

			target_regs_2_dump->push_back(cur_reg);
		}
	} // target_var_i loop.

	fprintf(stderr, "Saving %d regions to %s\n", target_regs_2_dump->size(), geno_op_fp);
	dump_geno_sig_regs_plain(target_regs_2_dump, testing_tag_sample_ids, false, geno_op_fp);
}


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

	char* op_dir)
{
	fprintf(stderr, "Building LMSE models using matrices %s, %s with sample id's in %s, %s for SNVs in range [%d-%d]\n",
			tag_genotype_matrix_fp, target_genotype_matrix_fp,
			tag_sample_ids_list_fp, target_sample_ids_list_fp,
			start_pos, end_pos);

	vector<char*>* tag_sample_ids = buffer_file(tag_sample_ids_list_fp);
	fprintf(stderr, "Loading tag genotype matrix with %d testing samples.\n", tag_sample_ids->size());
	vector<t_annot_region*>* tag_genotype_signal_regs = load_variant_genotype_signal_regions(tag_genotype_matrix_fp, tag_sample_ids);
	fprintf(stderr, "Loaded %d tag genotype signal regions with %d testing samples.\n", tag_genotype_signal_regs->size(), tag_sample_ids->size());
	t_restr_annot_region_list* restr_tag_genotype_signal_regs = restructure_annot_regions(tag_genotype_signal_regs);

	vector<char*>* target_sample_ids = buffer_file(target_sample_ids_list_fp);
	fprintf(stderr, "Loading target genotype matrix with %d testing samples.\n", target_sample_ids->size());
	vector<t_annot_region*>* target_genotype_signal_regs = load_variant_genotype_signal_regions(target_genotype_matrix_fp, target_sample_ids);
	fprintf(stderr, "Loaded %d target genotype signal regions with %d testing samples.\n", target_genotype_signal_regs->size(), target_sample_ids->size());
	t_restr_annot_region_list* restr_target_genotype_signal_regs = restructure_annot_regions(target_genotype_signal_regs);

	if (target_sample_ids->size() != tag_sample_ids->size())
	{
		fprintf(stderr, "Target and tag sample id's are not the same size: %d, %d\n", target_sample_ids->size(), tag_sample_ids->size());
		exit(0);
	}

	// Make sure the sample id's are matching.
	for (int i_s = 0; i_s < target_sample_ids->size(); i_s++)
	{
		if (!t_string::compare_strings(target_sample_ids->at(i_s), tag_sample_ids->at(i_s)))
		{
			fprintf(stderr, "Target and tag sample id's are not matching: %s, %s; Re-order and run again.\n", target_sample_ids->at(i_s), tag_sample_ids->at(i_s));
			exit(0);
		}
	} // i_s loop.

	// Decide on haplocoded signal; tag and target normalizers.
	double tag_normalizer = -1;
	double target_normalizer = -1;
	for (int var_i = 0; var_i < tag_genotype_signal_regs->size(); var_i++)
	{
		// Get the maximum.
		void** cur_var_sig_info = (void**)(tag_genotype_signal_regs->at(var_i)->data);
		char* cur_var_sig = (char*)(cur_var_sig_info[0]);

		for (int i_s = 0; i_s < target_sample_ids->size(); i_s++)
		{
			tag_normalizer = MAX(cur_var_sig[i_s], tag_normalizer);
		} // i_s loop.
	} // var_i loop.

	for (int var_i = 0; var_i < target_genotype_signal_regs->size(); var_i++)
	{
		// Get the maximum.
		void** cur_var_sig_info = (void**)(target_genotype_signal_regs->at(var_i)->data);
		char* cur_var_sig = (char*)(cur_var_sig_info[0]);

		for (int i_s = 0; i_s < target_sample_ids->size(); i_s++)
		{
			target_normalizer = MAX(cur_var_sig[i_s], target_normalizer);
		} // i_s loop.
	} // var_i loop.

	bool testing_flag = false;
	vector<char*>* testing_tag_sample_ids = NULL;
	vector<char*>* testing_target_sample_ids = NULL;
	if (check_file(testing_tag_genotype_matrix_fp) &&
		check_file(testing_target_genotype_matrix_fp))
	{
		fprintf(stderr, "Testing data exists, loading and setting it.\n");
		testing_flag = true;

		fprintf(stderr, "Loading and assigning testing target data.\n");
		testing_tag_sample_ids = buffer_file(testing_tag_sample_ids_list_fp);
		fprintf(stderr, "Loading testing tag genotype matrix with %d testing samples.\n", testing_tag_sample_ids->size());
		vector<t_annot_region*>* testing_tag_genotype_signal_regs = load_variant_genotype_signal_regions(testing_tag_genotype_matrix_fp, testing_tag_sample_ids);
		fprintf(stderr, "Loaded %d testing tag genotype signal regions with %d testing samples.\n", testing_tag_genotype_signal_regs->size(), testing_tag_sample_ids->size());

		vector<t_annot_region*>* tag_intersects = intersect_annot_regions(tag_genotype_signal_regs, testing_tag_genotype_signal_regs, true);
		fprintf(stderr, "Processing %d tag intersects.\n", tag_intersects->size());
		for (int i_int = 0; i_int < tag_intersects->size(); i_int++)
		{
			t_intersect_info* int_info = (t_intersect_info*)(tag_intersects->at(i_int)->data);
			t_annot_region* tag_reg = int_info->src_reg;
			t_annot_region* testing_tag_reg = int_info->dest_reg;

			if (t_string::compare_strings(tag_reg->name, testing_tag_reg->name))
			{
				void** cur_tag_reg_info = (void**)(tag_reg->data);
				void** cur_testing_tag_reg_info = (void**)(testing_tag_reg->data);
				cur_tag_reg_info[1] = cur_testing_tag_reg_info[0];
			}

			delete int_info;
		} // i_int loop.
		delete_annot_regions(tag_intersects);

		fprintf(stderr, "Testing assignment of tag regions to all training tag regions.\n");
		for (int i_reg = 0; i_reg < tag_genotype_signal_regs->size(); i_reg++)
		{
			void** cur_tag_reg_info = (void**)(tag_genotype_signal_regs->at(i_reg)->data);

			if (cur_tag_reg_info[1] == NULL)
			{
				fprintf(stderr, "Could not assign a testing tag genotype data for %s:%d-%d (%s)\n",
					tag_genotype_signal_regs->at(i_reg)->chrom,
					tag_genotype_signal_regs->at(i_reg)->start,
					tag_genotype_signal_regs->at(i_reg)->end,
					tag_genotype_signal_regs->at(i_reg)->name);

				exit(0);
			}
		} // i_reg loop.

		fprintf(stderr, "Loading and assigning testing target data.\n");
		testing_target_sample_ids = buffer_file(testing_target_sample_ids_list_fp);
		fprintf(stderr, "Loading testing target genotype matrix with %d testing samples.\n", testing_target_sample_ids->size());
		vector<t_annot_region*>* testing_target_genotype_signal_regs = load_variant_genotype_signal_regions(testing_target_genotype_matrix_fp, testing_target_sample_ids);
		fprintf(stderr, "Loaded %d testing target genotype signal regions with %d testing samples.\n", testing_target_genotype_signal_regs->size(), testing_target_sample_ids->size());

		vector<t_annot_region*>* target_intersects = intersect_annot_regions(target_genotype_signal_regs, testing_target_genotype_signal_regs, true);
		fprintf(stderr, "Processing %d target intersects.\n", target_intersects->size());
		for (int i_int = 0; i_int < target_intersects->size(); i_int++)
		{
			t_intersect_info* int_info = (t_intersect_info*)(target_intersects->at(i_int)->data);
			t_annot_region* target_reg = int_info->src_reg;
			t_annot_region* testing_target_reg = int_info->dest_reg;

			if (t_string::compare_strings(target_reg->name, testing_target_reg->name))
			{
				void** cur_target_reg_info = (void**)(target_reg->data);
				void** cur_testing_target_reg_info = (void**)(testing_target_reg->data);
				cur_target_reg_info[1] = cur_testing_target_reg_info[0];

				void** updated_info = new void*[5];
				updated_info[0] = cur_target_reg_info[0];
				updated_info[1] = cur_target_reg_info[1];
				updated_info[2] = new char[testing_target_sample_ids->size() + 2];
				memset(updated_info[2], 0, sizeof(char)  * (testing_target_sample_ids->size() + 2));

				target_reg->data = updated_info;
			}

			delete int_info;
		} // i_int loop.
		delete_annot_regions(target_intersects);

		fprintf(stderr, "Testing assignment of target regions to all training target regions.\n");
		for (int i_reg = 0; i_reg < target_genotype_signal_regs->size(); i_reg++)
		{
			void** cur_target_reg_info = (void**)(target_genotype_signal_regs->at(i_reg)->data);

			if (cur_target_reg_info[1] == NULL)
			{
				fprintf(stderr, "Could not assign a testing target genotype data for %s:%d-%d (%s)\n",
					target_genotype_signal_regs->at(i_reg)->chrom,
					target_genotype_signal_regs->at(i_reg)->start,
					target_genotype_signal_regs->at(i_reg)->end,
					target_genotype_signal_regs->at(i_reg)->name);

				exit(0);
			}
		} // i_reg loop.

		if (testing_target_sample_ids->size() != testing_tag_sample_ids->size())
		{
			fprintf(stderr, "Testing target and tag sample id's are not the same size: %d, %d\n", testing_target_sample_ids->size(), testing_tag_sample_ids->size());
			exit(0);
		}

		// Make sure the sample id's are matching.
		for (int i_s = 0; i_s < testing_target_sample_ids->size(); i_s++)
		{
			if (!t_string::compare_strings(testing_target_sample_ids->at(i_s), testing_tag_sample_ids->at(i_s)))
			{
				fprintf(stderr, "Testing target and tag sample id's are not matching: %s, %s; Re-order and run again.\n", testing_target_sample_ids->at(i_s), testing_tag_sample_ids->at(i_s));
				exit(0);
			}
		} // i_s loop.

		// Go over the testing data too for updating normalizers.
		for (int var_i = 0; var_i < testing_tag_genotype_signal_regs->size(); var_i++)
		{
			// Get the maximum.
			void** cur_var_sig_info = (void**)(testing_tag_genotype_signal_regs->at(var_i)->data);
			char* cur_var_sig = (char*)(cur_var_sig_info[0]);

			for (int i_s = 0; i_s < testing_tag_sample_ids->size(); i_s++)
			{
				tag_normalizer = MAX(cur_var_sig[i_s], tag_normalizer);
			} // i_s loop.
		} // var_i loop.

		for (int var_i = 0; var_i < testing_target_genotype_signal_regs->size(); var_i++)
		{
			// Get the maximum.
			void** cur_var_sig_info = (void**)(testing_target_genotype_signal_regs->at(var_i)->data);
			char* cur_var_sig = (char*)(cur_var_sig_info[0]);

			for (int i_s = 0; i_s < testing_target_sample_ids->size(); i_s++)
			{
				target_normalizer = MAX(cur_var_sig[i_s], target_normalizer);
			} // i_s loop.
		} // var_i loop.
	} // Testing data existence check.

	if (tag_normalizer < 2 || tag_normalizer > 3 ||
		target_normalizer < 2 || target_normalizer > 3)
	{
		fprintf(stderr, "tag or target normalizer is illegal or very improbably unlucky: %.1f, %.1f\n", tag_normalizer, target_normalizer);
		exit(0);
	}
	
	fprintf(stderr, "----------------------------------\n");
	if(tag_normalizer == 2)
	{
		fprintf(stderr, "Tag genotypes are unphased.\n");
	}
	else if (tag_normalizer == 3)
	{
		fprintf(stderr, "Tag genotypes are phased.\n");
	}

	if (target_normalizer == 2)
	{
		fprintf(stderr, "Target genotypes are unphased.\n");
	}
	else if (target_normalizer == 3)
	{
		fprintf(stderr, "Target genotypes are phased.\n");
	}
	fprintf(stderr, "----------------------------------\n");

	// This is the number of dimensions.
	int n_params_w_intercept = 2 * n_tags_vars_per_side + 1;

	// This is the number of data points in training.
	int training_sample_size = target_sample_ids->size();

	// Allocate the fitting data.
	double xi, yi, ei, chisq;
	gsl_matrix *X, *cov;
	gsl_vector *y, *w, *c;
	gsl_vector *testing_x;

	// SVD associated matrices/vectors.
	gsl_matrix* fullX = gsl_matrix_alloc(training_sample_size, n_params_w_intercept);
	gsl_matrix * svdV = gsl_matrix_alloc(n_params_w_intercept, n_params_w_intercept);
	gsl_vector * svdS = gsl_vector_alloc(n_params_w_intercept);
	gsl_vector * svdWork = gsl_vector_alloc(n_params_w_intercept);

	// Tag usage flag at the variant coordinates; this is necessary to save file and write per variant messages.
	int* tag_usage_flag = new int[n_params_w_intercept + 1];
	memset(tag_usage_flag, 0, sizeof(int) * (n_params_w_intercept + 1));

	FILE* f_per_var_stats = open_f("per_var_stats.txt", "w");

	// Process all the target variants.
	for (int target_chr_i = 0;
		target_chr_i < restr_target_genotype_signal_regs->chr_ids->size();
		target_chr_i++)
	{
		fprintf(stderr, "Building the models on chromosome %s\n", restr_target_genotype_signal_regs->chr_ids->at(target_chr_i));

		// Get the chromosome index for testring variant regions.
		int tag_chr_i = t_string::get_i_str(restr_tag_genotype_signal_regs->chr_ids, restr_target_genotype_signal_regs->chr_ids->at(target_chr_i));

		vector<t_annot_region*>* target_var_regs = restr_target_genotype_signal_regs->regions_per_chrom[target_chr_i];
		vector<t_annot_region*>* tag_var_regs = restr_tag_genotype_signal_regs->regions_per_chrom[tag_chr_i];

		// Following sets up n data points. We set this up over the samples.
		for (int target_var_i = 0; target_var_i < target_var_regs->size(); target_var_i++)
		{
			if (start_pos < 0 ||
				end_pos < 0 ||
				target_var_regs->at(target_var_i)->start >= start_pos && target_var_regs->at(target_var_i)->start <= end_pos)
			{

			}
			else
			{
				continue;
			}

			int closest_tag_var_left_i = -1;
			for (int tag_var_i = 0; (tag_var_i + 1)< tag_var_regs->size(); tag_var_i++)
			{
				if (tag_var_regs->at(tag_var_i)->start < target_var_regs->at(target_var_i)->start &&
					tag_var_regs->at(tag_var_i + 1)->start > target_var_regs->at(target_var_i)->start)
				{
					closest_tag_var_left_i = tag_var_i;
					break;
				}
			} // tag_var_i loop.

			if (closest_tag_var_left_i == -1)
			{
				continue;
			}

			if (__DUMP_INPUTATION_UTILS_MSGS__)
			{
				fprintf(stderr, "Closest tag to target %d: (%d, %d)\n", target_var_regs->at(target_var_i)->start,
					tag_var_regs->at(closest_tag_var_left_i)->start, tag_var_regs->at(closest_tag_var_left_i + 1)->start);
			}

			// Add the tag variant to the right.
			vector<t_annot_region*>* tag_var_regs_per_cur_target_var = new vector<t_annot_region*>();
			for (int tag_var_i = MAX(0, closest_tag_var_left_i - n_tags_vars_per_side + 1);
				tag_var_i <= closest_tag_var_left_i;
				tag_var_i++)
			{
				tag_var_regs_per_cur_target_var->push_back(tag_var_regs->at(tag_var_i));
			} // tag_var_i loop.

			for (int tag_var_i = closest_tag_var_left_i + 1;
				tag_var_i <= (closest_tag_var_left_i + n_tags_vars_per_side) && tag_var_i < tag_var_regs->size();
				tag_var_i++)
			{
				tag_var_regs_per_cur_target_var->push_back(tag_var_regs->at(tag_var_i));
			} // tag_var_i loop.

			  // Check the number of tag SNVs around the target SNVs, make sure they satisfy the requested number.
			if (tag_var_regs_per_cur_target_var->size() != 2 * n_tags_vars_per_side)
			{
				fprintf(stderr, "Could not find the correct number of tag variants for %s:%d (%d), will skip.\n",
					target_var_regs->at(target_var_i)->chrom, target_var_regs->at(target_var_i)->start,
					tag_var_regs_per_cur_target_var->size());
			}
			else
			{
				// Set the tag usage flag.
				memset(tag_usage_flag, 0, sizeof(int) * n_params_w_intercept);
				int n_tags_2_use = 0;
				for (int par_i = 0; par_i < n_params_w_intercept; par_i++)
				{
					tag_usage_flag[par_i] = 1;
				} // par_i loop.

				// Fill fullX.
				for (int sample_i = 0; sample_i < target_sample_ids->size(); sample_i++)
				{
					// Add the intercept; 1.0
					gsl_matrix_set(fullX, sample_i, 0, 1.0);

					// Jump over the intercept, set the tag genotype matrix entry.
					for (int tag_dim_i = 1; tag_dim_i <= 2 * n_tags_vars_per_side; tag_dim_i++)
					{
						void** reg_data = (void**)(tag_var_regs_per_cur_target_var->at(tag_dim_i - 1)->data);
						char* cur_tag_genotypes = (char*)(reg_data[0]);

						// Set the current entry for the full X matrix.
						double tagg_ij = cur_tag_genotypes[sample_i];

						// Normalize: Use tag normalization factor.
						tagg_ij /= tag_normalizer;

						gsl_matrix_set(fullX, sample_i, tag_dim_i, tagg_ij);
					} // tag_dim_i loop.
				} // sample_i loop.			

				// https://stackoverflow.com/questions/12304963/using-eigenvalues-to-test-for-singularity-identifying-collinear-columns
				// Generate the SVD of X.
				gsl_linalg_SV_decomp(fullX, svdV, svdS, svdWork);

				// Look for the low singular values that can create problems.
				double top_singular_value = gsl_vector_get(svdS, 0);
				for (int dim_i = 0; dim_i < n_params_w_intercept; dim_i++)
				{
					// Evaluate the current singular value.
					double cur_singular_val = gsl_vector_get(svdS, dim_i);
					//fprintf(stderr, "%d: %d: %lf\n", target_var_regs->at(target_var_i)->start, i, cur_singular_val);

					double sv_tolerance_fraction = 0.001;
					double min_component_strength = 0.01;

					// If the singular value is 0, identify which variants have components on the eigenvectors with 0 singular:
					if (cur_singular_val < top_singular_value * sv_tolerance_fraction &&
						cur_singular_val < 0.001)
					{
						fprintf(stderr, "Negligible singular value: %lf (dim: %d)\n", cur_singular_val, dim_i);

						// Check the entries that are on this eigenvector.
						vector<int>* colinear_var_par_i = new vector<int>();
						for (int par_i = 1; par_i < n_params_w_intercept; par_i++)
						{
							double cur_comp = gsl_matrix_get(svdV, par_i, dim_i);
							
							if (fabs(cur_comp) > min_component_strength)
							{
								fprintf(stderr, "%d has component on 0 eigenvalue: %lf\n", tag_var_regs_per_cur_target_var->at(par_i-1)->start, cur_comp);
								colinear_var_par_i->push_back(par_i-1);

								tag_usage_flag[par_i - 1] = 0;
							}
						} // par_i loop.

						//// Following performs variant selection to optimize the selections.
						//int n_colinear_vars = colinear_var_par_i->size();
						//int tag_i_2_set_unused = n_colinear_vars - 1;
						//while (tag_i_2_set_unused >= 0)
						//{
						//	if (tag_usage_flag[colinear_var_par_i->at(tag_i_2_set_unused)] == 0)
						//	{
						//		tag_i_2_set_unused--;
						//	}
						//	else
						//	{
						//		tag_usage_flag[colinear_var_par_i->at(tag_i_2_set_unused)] = 0;
						//		fprintf(stderr, "%d colinear tags; setting %d to be unused: %lf\n",
						//				n_colinear_vars,
						//				tag_var_regs_per_cur_target_var->at(colinear_var_par_i->at(tag_i_2_set_unused))->start);
						//		break;
						//	}
						//}
						delete colinear_var_par_i;
					}
				} // i loop.

				n_tags_2_use = 0;
				for (int par_i = 0; par_i < n_params_w_intercept-1; par_i++)
				{
					if (tag_usage_flag[par_i] == 1)
					{
						n_tags_2_use++;
					}
				} // par_i loop.

				fprintf(stderr, "%s:%d (%d tags)\n", target_var_regs->at(target_var_i)->chrom, target_var_regs->at(target_var_i)->start, n_tags_2_use);

				int n_tags_2_use_w_intercept = n_tags_2_use + 1;

				// Allocate the structures.
				X = gsl_matrix_alloc(training_sample_size, n_tags_2_use_w_intercept);
				testing_x = gsl_vector_alloc(n_tags_2_use_w_intercept);

				gsl_matrix * svdX = gsl_matrix_alloc(training_sample_size, n_tags_2_use_w_intercept);
				gsl_matrix * V = gsl_matrix_alloc(n_tags_2_use_w_intercept, n_tags_2_use_w_intercept);
				gsl_vector * S = gsl_vector_alloc(n_tags_2_use_w_intercept);
				gsl_vector * svd_work = gsl_vector_alloc(n_tags_2_use_w_intercept);

				// Target SNP genotypes vector.
				y = gsl_vector_alloc(training_sample_size);

				// Weights of each data point.
				w = gsl_vector_alloc(training_sample_size);

				// This is the coefficients matrix.
				c = gsl_vector_alloc(n_tags_2_use_w_intercept);
				cov = gsl_matrix_alloc(n_tags_2_use_w_intercept, n_tags_2_use_w_intercept);

				// Set the matrices and vectors for regression.
				for (int sample_i = 0; sample_i < target_sample_ids->size(); sample_i++)
				{
					// Add the intercept; 1.0
					gsl_matrix_set(X, sample_i, 0, 1.0);

					// Jump over the intercept, set the tag genotype matrix entry.
					int matrix_dim_i = 1;
					for (int tag_dim_i = 1; tag_dim_i <= 2 * n_tags_vars_per_side; tag_dim_i++)
					{
						if (tag_usage_flag[tag_dim_i - 1] == 1)
						{
							void** reg_data = (void**)(tag_var_regs_per_cur_target_var->at(tag_dim_i - 1)->data);
							char* cur_tag_genotypes = (char*)(reg_data[0]);

							// Set the current X matrix value.
							double tagg_ij = cur_tag_genotypes[sample_i];

							// Normalize with tag normalizer.
							tagg_ij /= tag_normalizer;

							//gsl_matrix_set(X, sample_i, tag_dim_i, tagg_ij);
							gsl_matrix_set(X, sample_i, matrix_dim_i, tagg_ij);
							matrix_dim_i++;
						}
					} // tag_dim_i loop.

					if (matrix_dim_i != n_tags_2_use_w_intercept)
					{
						fprintf(stderr, "Could not match the matrix dimension to the tags 2 use: %d, %d\n", matrix_dim_i, n_tags_2_use_w_intercept);
						exit(0);
					}

					// Set the target genotype vector entry.
					void** reg_data = (void**)(target_var_regs->at(target_var_i)->data);
					char* target_var_genotypes = (char*)(reg_data[0]);
					double targetg_ij = target_var_genotypes[sample_i];

					// Normalize with target normalizer.
					targetg_ij /= target_normalizer;

					gsl_vector_set(y, sample_i, targetg_ij);
				} // sample_i loop.				

				  // At this point, the target genotype array and the tag genotype matrix are setup.
				  // Now do the fitting.
				  // Do fitting.
				{
					// Do the fit.
					gsl_multifit_linear_workspace* work = gsl_multifit_linear_alloc(training_sample_size, n_tags_2_use_w_intercept);
					//gsl_multifit_linear_workspace* work = gsl_multifit_linear_alloc(training_sample_size, n_params_w_intercept);
					//gsl_multifit_linear_tsvd(X, y, c, cov, &chisq, &sv_rank, work); // Does not link!
					gsl_multifit_linear(X, y, c, cov, &chisq, work);
					gsl_multifit_linear_free(work);
				} // fitting block.

				  // Save the current parameters.
				  //#define C(i) (gsl_vector_get(c,(i)))
#define COV(i,j) (gsl_matrix_get(cov,(i),(j)))

				fprintf(stderr, "Saving Parameters for %s:%d (%d-%d): %d tags used.          \r",
					target_var_regs->at(target_var_i)->chrom,
					target_var_regs->at(target_var_i)->start,
					tag_var_regs->at(closest_tag_var_left_i)->start, tag_var_regs->at(closest_tag_var_left_i + 1)->start, n_tags_2_use_w_intercept);

				if (__DUMP_INPUTATION_UTILS_MSGS__)
				{
					int used_par_i = 0;
					for (int par_i = 0; par_i < n_params_w_intercept; par_i++)
					{
						if (par_i == 0)
						{
							fprintf(stderr, "Intercept\t%.17f\n", gsl_vector_get(c, used_par_i));
							used_par_i++;
						}
						else if (tag_usage_flag[par_i - 1] == 1)
						{
							fprintf(stderr, "%d\t%.17f\n", tag_var_regs_per_cur_target_var->at(par_i - 1)->start, gsl_vector_get(c, used_par_i));
							used_par_i++;
						}
						else
						{
							fprintf(stderr, "%d\tNA\n", tag_var_regs_per_cur_target_var->at(par_i - 1)->start);
						}
					} // par_i loop.
				}

				// Save the parameters file.
				char op_fp[1000];
				sprintf(op_fp, "%s/%d.params", op_dir, target_var_regs->at(target_var_i)->start);
				FILE* f_op = open_f(op_fp, "w");
				int used_par_i = 0;
				for (int par_i = 0; par_i < n_params_w_intercept; par_i++)
				{
					if (par_i == 0)
					{
						fprintf(f_op, "%.17f\n", gsl_vector_get(c, used_par_i));
						used_par_i++;
					}
					else if (tag_usage_flag[par_i - 1] == 1)
					{
						fprintf(f_op, "%.17f\n", gsl_vector_get(c, used_par_i));
						used_par_i++;
					}
					else
					{
						fprintf(f_op, "NA\n");
					}
				} // par_i loop.

				for (int tag_i = 0; tag_i < tag_var_regs_per_cur_target_var->size(); tag_i++)
				{
					fprintf(f_op, "%d\n", tag_var_regs_per_cur_target_var->at(tag_i)->start);
				} // par_i loop.

				fprintf(f_op, "%d\n", target_var_regs->at(target_var_i)->start);
				fclose(f_op);

				// Test if available.
				if (testing_flag)
				{
					vector<double>* per_testing_sample_errors = new vector<double>();
					vector<double>* per_testing_sample_non_ref_errors = new vector<double>();
					for (int testing_sample_i = 0; testing_sample_i < testing_tag_sample_ids->size(); testing_sample_i++)
					{
						// Reset the testing vector x: First element of x is 2.
						gsl_vector_set(testing_x, 0, 1.0);

						int matrix_dim_i = 1;
						for (int tag_dim_i = 1; tag_dim_i <= 2 * n_tags_vars_per_side; tag_dim_i++)
						{
							if (tag_usage_flag[tag_dim_i - 1] == 1)
							{
								void** cur_tag_reg_data = (void**)(tag_var_regs_per_cur_target_var->at(tag_dim_i - 1)->data);

								char* cur_testing_tag_genotypes = (char*)(cur_tag_reg_data[1]);
								if (cur_testing_tag_genotypes == NULL)
								{
									fprintf(stderr, "Testing Target region is null, exiting.\n");
									exit(0);
								}

								double cur_testing_sample_geno = (double)(cur_testing_tag_genotypes[testing_sample_i]);

								// Normalize with tag normalizer.
								cur_testing_sample_geno /= tag_normalizer;

								gsl_vector_set(testing_x, matrix_dim_i, cur_testing_sample_geno);
								matrix_dim_i++;
							} // tag usage flag check.
						} // par_i loop.

						if (matrix_dim_i != n_tags_2_use_w_intercept)
						{
							fprintf(stderr, "Could not match the matrix dimensions: %d, %d\n", 
								matrix_dim_i, n_tags_2_use_w_intercept);
							exit(0);
						}

						double geno_est = 0;
						double geno_err_est = 0;
						gsl_multifit_linear_est(testing_x, c, cov, &geno_est, &geno_err_est);

						// Estimate the real error.
						void** target_reg_data = (void**)(target_var_regs->at(target_var_i)->data);
						char* testing_target_var_genotypes = (char*)(target_reg_data[1]);
						if (testing_target_var_genotypes == NULL)
						{
							fprintf(stderr, "Testing Target region is null, exiting.\n");
							exit(0);
						}

						// Set the predicted genotype.
						char* pred_genoypes = (char*)(target_reg_data[2]);
						pred_genoypes[testing_sample_i] = (char)(100 * geno_est);

						double cur_testing_sample_target_geno = (double)(testing_target_var_genotypes[testing_sample_i]);

						// Normalize.
						cur_testing_sample_target_geno /= target_normalizer;

						double real_err = fabs(cur_testing_sample_target_geno - geno_est);

						if (__DUMP_INPUTATION_UTILS_MSGS__)
						{
							fprintf(stderr, "%s:%d-%d (%s): Error: %lf (%lf)\n",
								target_var_regs->at(target_var_i)->chrom,
								target_var_regs->at(target_var_i)->start,
								target_var_regs->at(target_var_i)->end,
								target_var_regs->at(target_var_i)->name,
								real_err,
								geno_err_est);
						}

						per_testing_sample_errors->push_back(real_err);

						if (__DUMP_INPUTATION_UTILS_MSGS__)
						{
							fprintf(stderr, "%s:%d: Sample %d: %.4f\n", target_var_regs->at(target_var_i)->chrom, target_var_regs->at(target_var_i)->start, testing_sample_i, real_err);
						}

						if (cur_testing_sample_target_geno != 0)
						{
							per_testing_sample_non_ref_errors->push_back(real_err);
						}
					} // testing_sample_i loop.

					  // Get the error statistics.
					double mean_err = 0;
					double std_dev_err = 0;
					get_stats(per_testing_sample_errors, mean_err, std_dev_err);

					double mean_non_ref_err = 0;
					double stddev_non_ref_err = 0;
					get_stats(per_testing_sample_non_ref_errors, mean_non_ref_err, stddev_non_ref_err);

					fprintf(stderr, "%s:%d-%d (%s): Avg Error: %lf (%d); Avg Non-ref Error: %lf (%d)          \n",
						target_var_regs->at(target_var_i)->chrom,
						target_var_regs->at(target_var_i)->start,
						target_var_regs->at(target_var_i)->end,
						target_var_regs->at(target_var_i)->name,
						mean_err,
						per_testing_sample_errors->size(),
						mean_non_ref_err, 
						per_testing_sample_non_ref_errors->size());

					fprintf(f_per_var_stats, "%s\t%d\t%d\t%s\t%lf\t%d\t%lf\t%d\t%d\n",
						target_var_regs->at(target_var_i)->chrom,
						target_var_regs->at(target_var_i)->start,
						target_var_regs->at(target_var_i)->end,
						target_var_regs->at(target_var_i)->name,
						mean_err,
						per_testing_sample_errors->size(),
						mean_non_ref_err,
						per_testing_sample_non_ref_errors->size(),
						n_tags_2_use_w_intercept);

					delete per_testing_sample_non_ref_errors;
					delete per_testing_sample_errors;
				} // testing check.

				  // Free memory.
				gsl_matrix_free(X);
				gsl_vector_free(y);
				gsl_vector_free(w);
				gsl_vector_free(c);
				gsl_matrix_free(cov);

				//printf("# covariance matrix:\n");
				//printf("[ %+.5e, %+.5e, %+.5e  \n",
				//	COV(0, 0), COV(0, 1), COV(0, 2));
				//printf("  %+.5e, %+.5e, %+.5e  \n",
				//	COV(1, 0), COV(1, 1), COV(1, 2));
				//printf("  %+.5e, %+.5e, %+.5e ]\n",
				//	COV(2, 0), COV(2, 1), COV(2, 2));
				//printf("# chisq = %g\n", chisq);
			} // tag snp count check.
		} // target_var_i loop.
	} // target_chr_i loop.

	fclose(f_per_var_stats);

	// Save results.
	char geno_op_fp[1000];
	sprintf(geno_op_fp, "%s/LMSE_imputed_genotypes.sig", op_dir);

	// Following sets up n data points. We set this up over the samples.
	vector<t_annot_region*>* target_regs_2_dump = new vector<t_annot_region*>();
	for (int target_var_i = 0; target_var_i < target_genotype_signal_regs->size(); target_var_i++)
	{
		void** target_reg_info = (void**)(target_genotype_signal_regs->at(target_var_i)->data);

		if (target_reg_info[1] == NULL || target_reg_info[2] == NULL)
		{

		}
		else
		{
			t_annot_region* cur_reg = duplicate_region(target_genotype_signal_regs->at(target_var_i));
			void** cur_reg_info = new void*[5];
			cur_reg_info[0] = (char*)(target_reg_info[2]);
			cur_reg->data = cur_reg_info;

			target_regs_2_dump->push_back(cur_reg);
		}
	} // target_var_i loop.

	fprintf(stderr, "Saving %d regions to %s\n", target_regs_2_dump->size(), geno_op_fp);
	dump_geno_sig_regs_plain(target_regs_2_dump, testing_tag_sample_ids, false, geno_op_fp);
}

void extract_genotype_probability_info_per_IMPUTE2_probabilities_output(char* IMPUTE2_output_fp, 
	t_IMPUTE2_op_col_info* IMPUTE2_op_col_info,
	char* target_geno_sig_regs_fp,
	char* target_sample_ids_list_fp,
	char* chrom, 
	char* op_prefix)
{
	vector<char*>* target_sample_ids = buffer_file(target_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d samples for the target regions.\n", target_sample_ids->size());
	vector<t_annot_region*>* target_geno_sig_regs = load_variant_genotype_signal_regions(target_geno_sig_regs_fp, target_sample_ids);
	fprintf(stderr, "Loaded %d target genotype signal regions.\n", target_geno_sig_regs->size());

	vector<char*>* genotyped_sample_ids = target_sample_ids;

	FILE* f_IMPUTE2_output = open_f(IMPUTE2_output_fp, "r");

	// Start reading each variant; write to the output file.
	vector<t_annot_region*>* impute2_genotyped_regs = new vector<t_annot_region*>();
	while (1)
	{
		char* cur_line = getline(f_IMPUTE2_output);
		if (cur_line == NULL)
		{
			break;
		}

		t_string_tokens* cur_line_toks = t_string::tokenize_by_chars(cur_line, " ");

		if (cur_line_toks->size() != IMPUTE2_op_col_info->geno_prob_starting_0_based_col_i + 3 * genotyped_sample_ids->size())
		{
			fprintf(stderr, "%d columns in %s but %d expected columns (%d, %d)\n",
				cur_line_toks->size(),
				IMPUTE2_op_col_info->geno_prob_starting_0_based_col_i + 3 * genotyped_sample_ids->size(),
				IMPUTE2_op_col_info->geno_prob_starting_0_based_col_i, genotyped_sample_ids->size());

			exit(0);
		}

		// Get the variant position.
		char* cur_var_name = cur_line_toks->at(IMPUTE2_op_col_info->name_col_i)->str();
		int cur_var_posn = atoi(cur_line_toks->at(IMPUTE2_op_col_info->posn_col_i)->str());
		t_annot_region* reg = get_empty_region();
		reg->chrom = t_string::copy_me_str(chrom);
		reg->start = cur_var_posn;
		reg->end = cur_var_posn;
		reg->strand = '+';
		reg->name = t_string::copy_me_str(cur_var_name);
		void** cur_reg_dat = new void*[2];
		reg->data = cur_reg_dat;
		impute2_genotyped_regs->push_back(reg);

		int i_tok = IMPUTE2_op_col_info->geno_prob_starting_0_based_col_i;

		double** cur_var_IMPUTE2_genotype_probs = new double*[genotyped_sample_ids->size() + 2];
		cur_reg_dat[0] = cur_var_IMPUTE2_genotype_probs;

		// Go over all the samples.
		for (int sample_i = 0; sample_i < genotyped_sample_ids->size(); sample_i++)
		{
			double* cur_sample_geno_probs = new double[3];
			int max_prob_geno = 0;
			double max_prob_geno_prob = -1;
			for (int geno = 0; geno < 3; geno++)
			{
				// Set the probability for this genotype.
				cur_sample_geno_probs[geno] = atof(cur_line_toks->at(i_tok)->str());
				//fprintf(stderr, "sample_i: %d; geno: %d, i_tok: %d\n", sample_i, geno, i_tok);

				i_tok++;
			} // geno loop.

			cur_var_IMPUTE2_genotype_probs[sample_i] = cur_sample_geno_probs;
		} // i_tok loop.

		if (i_tok != IMPUTE2_op_col_info->geno_prob_starting_0_based_col_i + 3 * genotyped_sample_ids->size())
		{
			fprintf(stderr, "Could not use all the entries: %d, %d\n", i_tok, IMPUTE2_op_col_info->geno_prob_starting_0_based_col_i + 3 * genotyped_sample_ids->size());
			exit(0);
		}

		t_string::clean_tokens(cur_line_toks);
	} // file reading loop.

	// Close files.
	close_f(f_IMPUTE2_output, IMPUTE2_output_fp);

	fprintf(stderr, "Loaded %d regions from IMPUTE2 output, matching to the target regions and saving.\n", impute2_genotyped_regs->size());

	// Load the target variant regions.
	vector<t_annot_region*>* intersects = intersect_annot_regions(target_geno_sig_regs, impute2_genotyped_regs, true);
	fprintf(stderr, "Processing %d intersects.\n", intersects->size());

	for (int i_int = 0; i_int < intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);

		t_annot_region* cur_target_reg = int_info->src_reg;
		t_annot_region* cur_impute2_geno_reg = int_info->dest_reg;

		void** target_reg_info = (void**)(cur_target_reg->data);
		target_reg_info[1] = cur_impute2_geno_reg;
	} // i_int loop.

	fprintf(stderr, "Saving known genotypes and genotype probabilities.\n");
	char known_geno_op_fp[1000];
	sprintf(known_geno_op_fp, "%s.known.gz", op_prefix);
	char geno_prob_op_fp[1000];
	sprintf(geno_prob_op_fp, "%s.probs.gz", op_prefix);

	FILE* f_known_geno = open_f(known_geno_op_fp, "w");
	FILE* f_prob_info = open_f(geno_prob_op_fp, "w");
	fprintf(stderr, "Saving probability info and known genotypes.\n");
	for (int i_reg = 0; i_reg < target_geno_sig_regs->size(); i_reg++)
	{
		if (i_reg % 1000 == 0)
		{
			fprintf(stderr, "Saving %d. region's info.           \r", i_reg);
		}

		void** target_reg_info = (void**)(target_geno_sig_regs->at(i_reg)->data);
		if (target_reg_info[1] != NULL)
		{
			char* cur_target_reg_geno_sigs = (char*)(target_reg_info[0]);
			t_annot_region* impute2_prob_info_reg = (t_annot_region*)(target_reg_info[1]);
			void** impute2_reg_info = (void**)(impute2_prob_info_reg->data);
			double** cur_target_reg_prob_info = (double**)(impute2_reg_info[0]);

			// Write.
			for (int i_s = 0; i_s < genotyped_sample_ids->size(); i_s++)
			{
				//fprintf(f_known_geno, "%s\t%s\t%d\n", vcf_sample_ids->at(i_s), target_geno_sig_regs->at(i_reg)->name, (int)(cur_target_reg_geno_sigs[i_s]));
				fprintf(f_known_geno, "%d\t%s\t%d\n", i_s, target_geno_sig_regs->at(i_reg)->name, (int)(cur_target_reg_geno_sigs[i_s]));

				//fprintf(f_prob_info, "%s\t%s\t%s\n", vcf_sample_ids->at(i_s), vcf_prob_info_reg->name, cur_target_reg_prob_info[i_s]);
				fprintf(f_prob_info, "%d\t%s\t%.3f,%.3f,%.3f\n", i_s, impute2_prob_info_reg->name, cur_target_reg_prob_info[i_s][0], cur_target_reg_prob_info[i_s][1], cur_target_reg_prob_info[i_s][2]);
			} // i_s loop.
		}
	} // i_reg loop.

	  // Close files.
	close_f(f_known_geno, known_geno_op_fp);
	close_f(f_prob_info, geno_prob_op_fp);
}

// Parse the output and write.
void extract_genotype_probability_info_per_GT_entry_id_VCF(char* vcf_fp,							// This is the VCF file from which the genotypes are read.
													char* vcf_sample_ids_list_fp,		// This is the sample ids list file path.
													char* target_geno_sig_regs_fp, 
													char* target_sample_ids_list_fp,
													char* GT_entry_id,
													char* op_prefix)						// Output file path.
{
	vector<char*>* target_sample_ids = buffer_file(target_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d samples for the target regions.\n", target_sample_ids->size());
	vector<t_annot_region*>* target_geno_sig_regs = load_variant_genotype_signal_regions(target_geno_sig_regs_fp, target_sample_ids);
	fprintf(stderr, "Loaded %d target genotype signal regions.\n", target_geno_sig_regs->size());

	if (!check_file(vcf_sample_ids_list_fp))
	{
		fprintf(stderr, "Could not find sample id's list @ %s\n", vcf_sample_ids_list_fp);
		exit(0);
	}

	vector<char*>* vcf_sample_ids = buffer_file(vcf_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d sample ids.\n", vcf_sample_ids->size());

	// The sample id's must match.
	if (vcf_sample_ids->size() != target_sample_ids->size())
	{
		fprintf(stderr, "The samples ids do not match: %d, %d\n", vcf_sample_ids->size(), target_sample_ids->size());
		exit(0);
	}
	else
	{
		for (int i_s = 0; i_s < target_sample_ids->size(); i_s++)
		{
			if (!t_string::compare_strings(target_sample_ids->at(i_s), vcf_sample_ids->at(i_s)))
			{
				fprintf(stderr, "Could not match the sample id's.\n");
				exit(0);
			}
		} // i_s loop.
	}

	FILE* f_vcf = open_f(vcf_fp, "r");

	char* buff = new char[100 * 1000];

	vector<t_annot_region*>* vcf_geno_prob_regs = new vector<t_annot_region*>();
	while (1)
	{
		char* cur_line = getline(f_vcf);
		if (cur_line == NULL)
		{
			break;
		}

		if (cur_line[0] == '#')
		{
			delete[] cur_line;
			continue;
		}

		if (vcf_geno_prob_regs->size() % 1000 == 0)
		{
			fprintf(stderr, "Processing %d. VCF region.           \r", vcf_geno_prob_regs->size());
		}

		// #CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT
		char chrom[100];
		char posn_str[100];
		char id[100];
		char ref[100];
		char alt[100];
		char qual[100];
		char filter[100];
		char info[10000];
		char format[100];

		int i_cur_char = 0;
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(chrom, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(posn_str, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(id, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(ref, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(alt, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(qual, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(filter, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(info, buff);
		t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);
		strcpy(format, buff);

		// Get the genotype entry: GT
		int format_char_i = 0;
		int format_tok_i = 0;
		int custom_entry_i = -1;
		int l_format = t_string::string_length(format);
		char* cur_tok = new char[l_format + 2];
		while (t_string::get_next_token(format, cur_tok, l_format, ":", format_char_i))
		{
			if (t_string::compare_strings(cur_tok, GT_entry_id))
			{
				custom_entry_i = format_tok_i;
				break;
			}

			format_tok_i++;
		} // format string parsing loop.

		if (__DUMP_INPUTATION_UTILS_MSGS__)
		{
			fprintf(stderr, "Found %s entry @ %d\n", GT_entry_id, custom_entry_i);
		}

		if (custom_entry_i == -1)
		{
			fprintf(stderr, "Could not find %s entry in %s, skipping\n",
					GT_entry_id,
					cur_line);
			delete[] cur_line;
			continue;
		}

		bool correctly_parsed_all_genotypes = true;
		char** per_sample_entries = new char*[vcf_sample_ids->size() + 2];
		for (int i_s = 0; i_s < vcf_sample_ids->size(); i_s++)
		{
			t_string::get_next_token(cur_line, buff, 100 * 1000, "\t", i_cur_char);

			// Get the entry.
			int cur_tok_i = 0;
			int buff_char_i = 0;
			char tok_buff[1000];
			while (t_string::get_next_token(buff, tok_buff, 1000, ":", buff_char_i))
			{
				if (cur_tok_i == custom_entry_i)
				{
					break;
				}

				cur_tok_i++;
			} // genotype entry parsing loop.

			//fprintf(f_op, "%s\t%s\t%s\n", vcf_sample_ids->at(i_s), id, tok_buff);
			per_sample_entries[i_s] = t_string::copy_me_str(tok_buff);
			//fprintf(stderr, "%s: %s (%d)\n", vcf_sample_ids->at(i_s), buff, cur_var_alt_alle_cnt_sig[i_s]);
			//getc(stdin);
		} // i_s loop.

		t_annot_region* new_reg = get_empty_region();
		new_reg->chrom = t_string::copy_me_str(chrom);
		new_reg->start = translate_coord(atoi(posn_str), VCF_COORDS::start_base, CODEBASE_COORDS::start_base);
		int l_ref_allele = t_string::string_length(ref);
		new_reg->end = new_reg->start + l_ref_allele - 1;
		new_reg->name = t_string::copy_me_str(id);
		new_reg->strand = '+';

		new_reg->data = per_sample_entries;

		// Add this region to the other regions.
		vcf_geno_prob_regs->push_back(new_reg);

		if (cur_line[i_cur_char] != 0)
		{
			fprintf(stderr, "Could not finish the whole line.\n");
			exit(0);
		}

		delete[] cur_line;
	} // vcf file reading loop.

	fprintf(stderr, "Loaded %d regions from %s, matching to the target regions and saving.\n", vcf_geno_prob_regs->size(), vcf_fp);

	// Load the target variant regions.
	vector<t_annot_region*>* intersects = intersect_annot_regions(target_geno_sig_regs, vcf_geno_prob_regs, true);
	fprintf(stderr, "Processing %d intersects.\n", intersects->size());

	for (int i_int = 0; i_int < intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);

		t_annot_region* cur_target_reg = int_info->src_reg;
		t_annot_region* cur_vcf_geno_reg = int_info->dest_reg;

		//char** vcf_reg_prob_info = (char**)(cur_vcf_geno_reg->data);
		void** target_reg_info = (void**)(cur_target_reg->data);
		//target_reg_info[1] = vcf_reg_prob_info;
		target_reg_info[1] = cur_vcf_geno_reg;
	} // i_int loop.

	fprintf(stderr, "Saving known genotypes and genotype probabilities.\n");
	char known_geno_op_fp[1000];
	sprintf(known_geno_op_fp, "%s.known.gz", op_prefix);
	char geno_prob_op_fp[1000];
	sprintf(geno_prob_op_fp, "%s.probs.gz", op_prefix);

	FILE* f_known_geno = open_f(known_geno_op_fp, "w");
	FILE* f_prob_info = open_f(geno_prob_op_fp, "w");
	fprintf(stderr, "Saving probability info and known genotypes.\n");
	for (int i_reg = 0; i_reg < target_geno_sig_regs->size(); i_reg++)
	{
		if (i_reg % 1000 == 0)
		{
			fprintf(stderr, "Saving %d. region's info.           \r", i_reg);
		}

		void** target_reg_info = (void**)(target_geno_sig_regs->at(i_reg)->data);
		if (target_reg_info[1] != NULL)
		{
			char* cur_target_reg_geno_sigs = (char*)(target_reg_info[0]);
			t_annot_region* vcf_prob_info_reg = (t_annot_region*)(target_reg_info[1]);
			char** cur_target_reg_prob_info = (char**)(vcf_prob_info_reg->data);

			// Write.
			for (int i_s = 0; i_s < vcf_sample_ids->size(); i_s++)
			{
				//fprintf(f_known_geno, "%s\t%s\t%d\n", vcf_sample_ids->at(i_s), target_geno_sig_regs->at(i_reg)->name, (int)(cur_target_reg_geno_sigs[i_s]));
				fprintf(f_known_geno, "%d\t%s\t%d\n", i_s, target_geno_sig_regs->at(i_reg)->name, (int)(cur_target_reg_geno_sigs[i_s]));

				//fprintf(f_prob_info, "%s\t%s\t%s\n", vcf_sample_ids->at(i_s), vcf_prob_info_reg->name, cur_target_reg_prob_info[i_s]);			
				fprintf(f_prob_info, "%d\t%s\t%s\n", i_s, vcf_prob_info_reg->name, cur_target_reg_prob_info[i_s]);
			} // i_s loop.
		}
	} // i_reg loop.

	// Close files.
	close_f(f_known_geno, known_geno_op_fp);
	close_f(f_prob_info, geno_prob_op_fp);
}

// http://linux.math.tifr.res.in/manuals/html/gsl-ref-html/gsl-ref_35.html
void train_vicinity_based_LMSE_imputation_model_ALL_FEATURES(int n_tags_vars_per_side,
												char* tag_genotype_matrix_fp,
												char* tag_sample_ids_list_fp,
												char* target_genotype_matrix_fp,
												char* target_sample_ids_list_fp,

												// The testing data.
												char* testing_tag_genotype_matrix_fp,
												char* testing_tag_sample_ids_list_fp,
												char* testing_target_genotype_matrix_fp,
												char* testing_target_sample_ids_list_fp,

												char* op_dir)
{
	fprintf(stderr, "Building LMSE models using matrices %s, %s with sample id's in %s, %s\n", 
			tag_genotype_matrix_fp, target_genotype_matrix_fp, 
			tag_sample_ids_list_fp, target_sample_ids_list_fp);

	vector<char*>* tag_sample_ids = buffer_file(tag_sample_ids_list_fp);
	fprintf(stderr, "Loading tag genotype matrix with %d testing samples.\n", tag_sample_ids->size());
	vector<t_annot_region*>* tag_genotype_signal_regs = load_variant_genotype_signal_regions(tag_genotype_matrix_fp, tag_sample_ids);
	fprintf(stderr, "Loaded %d tag genotype signal regions with %d testing samples.\n", tag_genotype_signal_regs->size(), tag_sample_ids->size());
	t_restr_annot_region_list* restr_tag_genotype_signal_regs = restructure_annot_regions(tag_genotype_signal_regs);

	vector<char*>* target_sample_ids = buffer_file(target_sample_ids_list_fp);
	fprintf(stderr, "Loading target genotype matrix with %d testing samples.\n", target_sample_ids->size());
	vector<t_annot_region*>* target_genotype_signal_regs = load_variant_genotype_signal_regions(target_genotype_matrix_fp, target_sample_ids);
	fprintf(stderr, "Loaded %d target genotype signal regions with %d testing samples.\n", target_genotype_signal_regs->size(), target_sample_ids->size());
	t_restr_annot_region_list* restr_target_genotype_signal_regs = restructure_annot_regions(target_genotype_signal_regs);

	if (target_sample_ids->size() != tag_sample_ids->size())
	{
		fprintf(stderr, "Target and tag sample id's are not the same size: %d, %d\n", target_sample_ids->size(), tag_sample_ids->size());
		exit(0);
	}

	// Make sure the sample id's are matching.
	for (int i_s = 0; i_s < target_sample_ids->size(); i_s++)
	{
		if (!t_string::compare_strings(target_sample_ids->at(i_s), tag_sample_ids->at(i_s)))
		{
			fprintf(stderr, "Target and tag sample id's are not matching: %s, %s; Re-order and run again.\n", target_sample_ids->at(i_s), tag_sample_ids->at(i_s));
			exit(0);
		}
	} // i_s loop.

	bool testing_flag = false;
	vector<char*>* testing_tag_sample_ids = NULL;
	vector<char*>* testing_target_sample_ids = NULL;
	if (check_file(testing_tag_genotype_matrix_fp) &&
		check_file(testing_target_genotype_matrix_fp))
	{
		fprintf(stderr, "Testing data exists, loading and setting it.\n");
		testing_flag = true;

		fprintf(stderr, "Loading and assigning testing target data.\n");
		testing_tag_sample_ids = buffer_file(testing_tag_sample_ids_list_fp);
		fprintf(stderr, "Loading testing tag genotype matrix with %d testing samples.\n", testing_tag_sample_ids->size());
		vector<t_annot_region*>* testing_tag_genotype_signal_regs = load_variant_genotype_signal_regions(testing_tag_genotype_matrix_fp, testing_tag_sample_ids);
		fprintf(stderr, "Loaded %d testing tag genotype signal regions with %d testing samples.\n", testing_tag_genotype_signal_regs->size(), testing_tag_sample_ids->size());

		vector<t_annot_region*>* tag_intersects = intersect_annot_regions(tag_genotype_signal_regs, testing_tag_genotype_signal_regs, true);
		fprintf(stderr, "Processing %d tag intersects.\n", tag_intersects->size());
		for (int i_int = 0; i_int < tag_intersects->size(); i_int++)
		{
			t_intersect_info* int_info = (t_intersect_info*)(tag_intersects->at(i_int)->data);
			t_annot_region* tag_reg = int_info->src_reg;
			t_annot_region* testing_tag_reg = int_info->dest_reg;

			if (t_string::compare_strings(tag_reg->name, testing_tag_reg->name))
			{
				void** cur_tag_reg_info = (void**)(tag_reg->data);
				void** cur_testing_tag_reg_info = (void**)(testing_tag_reg->data);
				cur_tag_reg_info[1] = cur_testing_tag_reg_info[0];
			}

			delete int_info;
		} // i_int loop.
		delete_annot_regions(tag_intersects);

		fprintf(stderr, "Testing assignment of tag regions to all training tag regions.\n");
		for (int i_reg = 0; i_reg < tag_genotype_signal_regs->size(); i_reg++)
		{
			void** cur_tag_reg_info = (void**)(tag_genotype_signal_regs->at(i_reg)->data);

			if (cur_tag_reg_info[1] == NULL)
			{
				fprintf(stderr, "Could not assign a testing tag genotype data for %s:%d-%d (%s)\n", 
					tag_genotype_signal_regs->at(i_reg)->chrom, 
					tag_genotype_signal_regs->at(i_reg)->start, 
					tag_genotype_signal_regs->at(i_reg)->end,
					tag_genotype_signal_regs->at(i_reg)->name);

				exit(0);
			}
		} // i_reg loop.

		fprintf(stderr, "Loading and assigning testing target data.\n");
		testing_target_sample_ids = buffer_file(testing_target_sample_ids_list_fp);
		fprintf(stderr, "Loading testing target genotype matrix with %d testing samples.\n", testing_target_sample_ids->size());
		vector<t_annot_region*>* testing_target_genotype_signal_regs = load_variant_genotype_signal_regions(testing_target_genotype_matrix_fp, testing_target_sample_ids);
		fprintf(stderr, "Loaded %d testing target genotype signal regions with %d testing samples.\n", testing_target_genotype_signal_regs->size(), testing_target_sample_ids->size());

		vector<t_annot_region*>* target_intersects = intersect_annot_regions(target_genotype_signal_regs, testing_target_genotype_signal_regs, true);
		fprintf(stderr, "Processing %d target intersects.\n", target_intersects->size());
		for (int i_int = 0; i_int < target_intersects->size(); i_int++)
		{
			t_intersect_info* int_info = (t_intersect_info*)(target_intersects->at(i_int)->data);
			t_annot_region* target_reg = int_info->src_reg;
			t_annot_region* testing_target_reg = int_info->dest_reg;

			if (t_string::compare_strings(target_reg->name, testing_target_reg->name))
			{
				void** cur_target_reg_info = (void**)(target_reg->data);
				void** cur_testing_target_reg_info = (void**)(testing_target_reg->data);
				cur_target_reg_info[1] = cur_testing_target_reg_info[0];
			}

			delete int_info;
		} // i_int loop.
		delete_annot_regions(target_intersects);

		fprintf(stderr, "Testing assignment of target regions to all training target regions.\n");
		for (int i_reg = 0; i_reg < target_genotype_signal_regs->size(); i_reg++)
		{
			void** cur_target_reg_info = (void**)(target_genotype_signal_regs->at(i_reg)->data);

			if (cur_target_reg_info[1] == NULL)
			{
				fprintf(stderr, "Could not assign a testing target genotype data for %s:%d-%d (%s)\n",
					target_genotype_signal_regs->at(i_reg)->chrom,
					target_genotype_signal_regs->at(i_reg)->start,
					target_genotype_signal_regs->at(i_reg)->end,
					target_genotype_signal_regs->at(i_reg)->name);

				exit(0);
			}
		} // i_reg loop.

		if (testing_target_sample_ids->size() != testing_tag_sample_ids->size())
		{
			fprintf(stderr, "Testing target and tag sample id's are not the same size: %d, %d\n", testing_target_sample_ids->size(), testing_tag_sample_ids->size());
			exit(0);
		}

		// Make sure the sample id's are matching.
		for (int i_s = 0; i_s < testing_target_sample_ids->size(); i_s++)
		{
			if (!t_string::compare_strings(testing_target_sample_ids->at(i_s), testing_tag_sample_ids->at(i_s)))
			{
				fprintf(stderr, "Testing target and tag sample id's are not matching: %s, %s; Re-order and run again.\n", testing_target_sample_ids->at(i_s), testing_tag_sample_ids->at(i_s));
				exit(0);
			}
		} // i_s loop.

	} // Testing data existence check.

	// Allocate the fitting data.
	//int i, n;
	double xi, yi, ei, chisq;
	gsl_matrix *X, *cov;
	gsl_vector *y, *w, *c;
	gsl_vector *testing_x;

	// This is the number of data points in training.
	//n = atoi(argv[1]);
	//int n = training_sample_ids->size();
	int training_sample_size = target_sample_ids->size();

	// This is the number of dimensions.
	int n_params_w_intercept = 2 * n_tags_vars_per_side + 1;

	// Tag SNP genotypes
	//X = gsl_matrix_alloc(n, 3);
	X = gsl_matrix_alloc(training_sample_size, n_params_w_intercept);
	testing_x = gsl_vector_alloc(n_params_w_intercept);

	// Target SNP genotypes vector.
	y = gsl_vector_alloc(training_sample_size);
	
	// Weights of each data point.
	w = gsl_vector_alloc(training_sample_size);

	// This is the coefficients matrix.
	c = gsl_vector_alloc(n_params_w_intercept);
	cov = gsl_matrix_alloc(n_params_w_intercept, n_params_w_intercept);

	// Process all the target variants.
	for (int target_chr_i = 0;
		target_chr_i < restr_target_genotype_signal_regs->chr_ids->size();
		target_chr_i++)
	{
		fprintf(stderr, "Building the models on chromosome %s\n", restr_target_genotype_signal_regs->chr_ids->at(target_chr_i));

		// Get the chromosome index for testring variant regions.
		int tag_chr_i = t_string::get_i_str(restr_tag_genotype_signal_regs->chr_ids, restr_target_genotype_signal_regs->chr_ids->at(target_chr_i));

		vector<t_annot_region*>* target_var_regs = restr_target_genotype_signal_regs->regions_per_chrom[target_chr_i];
		vector<t_annot_region*>* tag_var_regs = restr_tag_genotype_signal_regs->regions_per_chrom[tag_chr_i];

		// Following sets up n data points. We set this up over the samples.
		for (int target_var_i = 0; target_var_i < target_var_regs->size(); target_var_i++)
		{
			int closest_tag_var_left_i = -1;
			for (int tag_var_i = 0; (tag_var_i+1)< tag_var_regs->size(); tag_var_i++)
			{
				if (tag_var_regs->at(tag_var_i)->start < target_var_regs->at(target_var_i)->start &&
					tag_var_regs->at(tag_var_i+1)->start > target_var_regs->at(target_var_i)->start)
				{
					closest_tag_var_left_i = tag_var_i;
					break;
				}
			} // tag_var_i loop.

			if (closest_tag_var_left_i == -1)
			{
				continue;
			}

			if (__DUMP_INPUTATION_UTILS_MSGS__)
			{
				fprintf(stderr, "Closest tag to target %d: (%d, %d)\n", target_var_regs->at(target_var_i)->start,
						tag_var_regs->at(closest_tag_var_left_i)->start, tag_var_regs->at(closest_tag_var_left_i + 1)->start);
			}

			// Add the tag variant to the right.
			vector<t_annot_region*>* tag_var_regs_per_cur_target_var = new vector<t_annot_region*>();
			for (int tag_var_i = MAX(0, closest_tag_var_left_i - n_tags_vars_per_side + 1);
				tag_var_i <= closest_tag_var_left_i;
				tag_var_i++)
			{
				tag_var_regs_per_cur_target_var->push_back(tag_var_regs->at(tag_var_i));
			} // tag_var_i loop.
					
			for (int tag_var_i = closest_tag_var_left_i + 1;
				tag_var_i <= (closest_tag_var_left_i + n_tags_vars_per_side) && tag_var_i < tag_var_regs->size();
				tag_var_i++)
			{
				tag_var_regs_per_cur_target_var->push_back(tag_var_regs->at(tag_var_i));
			} // tag_var_i loop.

			// Check the number of tag SNVs around the target SNVs, make sure they satisfy the requested number.
			if (tag_var_regs_per_cur_target_var->size() != 2 * n_tags_vars_per_side)
			{
				fprintf(stderr, "Could not find the correct number of tag variants for %s:%d (%d), will skip.\n",
						target_var_regs->at(target_var_i)->chrom, target_var_regs->at(target_var_i)->start,
						tag_var_regs_per_cur_target_var->size());
			}
			else
			{
				// Check the tag variant genotype correlation among samples.
				double cur_n_per_geno[3];	
				double* dim_i_genotypes = new double[testing_tag_sample_ids->size() + 2];
				double* dim_j_genotypes = new double[testing_tag_sample_ids->size() + 2];
				int* tag_usage_flag = new int[2 * n_tags_vars_per_side + 3];
				memset(tag_usage_flag, 0, sizeof(int) * (2 * n_tags_vars_per_side + 3));
				for (int tag_dim_i = 0; tag_dim_i < 2 * n_tags_vars_per_side; tag_dim_i++)
				{
					memset(cur_n_per_geno, 0, sizeof(double) * 3);

					for (int sample_i = 0;
						sample_i < tag_sample_ids->size();
						sample_i++)
					{
						void** reg_data = (void**)(tag_var_regs_per_cur_target_var->at(tag_dim_i)->data);
						char* cur_tag_genotypes = (char*)(reg_data[0]);

						dim_i_genotypes[sample_i] = cur_tag_genotypes[sample_i];

						// Set the 
						cur_n_per_geno[cur_tag_genotypes[sample_i]]++;
					} // sample_i loop.

					  // Analyze the distribution.
					if (cur_n_per_geno[0] == tag_sample_ids->size() ||
						cur_n_per_geno[1] == tag_sample_ids->size() ||
						cur_n_per_geno[2] == tag_sample_ids->size())
					{
						fprintf(stderr, "Dimension %d provides no info.\n",
							tag_var_regs_per_cur_target_var->at(tag_dim_i)->chrom, tag_var_regs_per_cur_target_var->at(tag_dim_i)->start);
					}

					double min_distance_dim_i = 1000000;
					for (int tag_dim_j = 0; tag_dim_j < tag_dim_i; tag_dim_j++)
					{
						for (int sample_i = 0;
							sample_i < tag_sample_ids->size();
							sample_i++)
						{
							void** reg_data = (void**)(tag_var_regs_per_cur_target_var->at(tag_dim_j)->data);
							char* cur_tag_genotypes = (char*)(reg_data[0]);

							dim_j_genotypes[sample_i] = cur_tag_genotypes[sample_i];
						} // sample_i loop.

						// Compute the distance.
						double cur_dist = 0;
						for (int sample_i = 0; sample_i < tag_sample_ids->size(); sample_i++)
						{
							cur_dist += (dim_i_genotypes[sample_i] != dim_j_genotypes[sample_i]);
						} // samepl_i loop.

						// Update the min distance.
						min_distance_dim_i = MIN(min_distance_dim_i, cur_dist);
					} // tag_dim_j;

					// Do greedy selection.
					if (min_distance_dim_i > 0)
					{
						tag_usage_flag[tag_dim_i] = 1;
					}
				} // tag_dim_i;

				delete[] dim_i_genotypes;
				delete[] dim_j_genotypes;

				// Set the matrices and vectors for regression.
				for (int sample_i = 0; sample_i < target_sample_ids->size(); sample_i++)
				{
					// Add the intercept; 1.0
					gsl_matrix_set(X, sample_i, 0, 2.0);

					// Jump over the intercept, set the tag genotype matrix entry.
					for (int tag_dim_i = 1; tag_dim_i <= 2 * n_tags_vars_per_side; tag_dim_i++)
					{
						void** reg_data = (void**)(tag_var_regs_per_cur_target_var->at(tag_dim_i - 1)->data);
						char* cur_tag_genotypes = (char*)(reg_data[0]);

						// Set the 
						double tagg_ij = cur_tag_genotypes[sample_i];
						gsl_matrix_set(X, sample_i, tag_dim_i, tagg_ij);
					} // tag_dim_i loop.

					// Set the target genotype vector entry.
					void** reg_data = (void**)(target_var_regs->at(target_var_i)->data);
					char* target_var_genotypes = (char*)(reg_data[0]);
					double targetg_ij = target_var_genotypes[sample_i];
					gsl_vector_set(y, sample_i, targetg_ij);	
				} // sample_i loop.

				// At this point, the target genotype array and the tag genotype matrix are setup.
				// Now do the fitting.
				// Do fitting.
				{
					size_t sv_rank;

					gsl_multifit_linear_workspace* work = gsl_multifit_linear_alloc(training_sample_size, n_params_w_intercept);
					//gsl_multifit_wlinear(X, w, y, c, cov, &chisq, work);
					//gsl_multifit_linear_tsvd(X, y, c, cov, &chisq, &sv_rank, work);
					gsl_multifit_linear(X, y, c, cov, &chisq, work);
					gsl_multifit_linear_free(work);
				} // fitting block.

				// Save the current parameters.
//#define C(i) (gsl_vector_get(c,(i)))
#define COV(i,j) (gsl_matrix_get(cov,(i),(j)))
				
				fprintf(stderr, "Saving Parameters for %s:%d (%d-%d)          \r", 
						target_var_regs->at(target_var_i)->chrom, 
						target_var_regs->at(target_var_i)->start,
						tag_var_regs->at(closest_tag_var_left_i)->start, tag_var_regs->at(closest_tag_var_left_i + 1)->start);

				if (__DUMP_INPUTATION_UTILS_MSGS__)
				{
					for (int par_i = 0; par_i < n_params_w_intercept; par_i++)
					{
						if (par_i == 0)
						{
							fprintf(stderr, "Intercept\t%.17f\n", gsl_vector_get(c, par_i));
						}
						else
						{
							fprintf(stderr, "%d\t%.17f\n", tag_var_regs_per_cur_target_var->at(par_i - 1)->start, gsl_vector_get(c, par_i));
						}

						//printf("# best fit: Y = %g + %g X + %g X^2\n",
						//	C(0), C(1), C(2));
					} // par_i loop.
				}

				// Save the parameters file.
				char op_fp[1000];
				sprintf(op_fp, "%s/%d.params", op_dir, target_var_regs->at(target_var_i)->start);
				FILE* f_op = open_f(op_fp, "w");
				for (int par_i = 0; par_i < n_params_w_intercept; par_i++)
				{
					fprintf(f_op, "%.17f\n", gsl_vector_get(c, par_i));
				} // par_i loop.

				for (int tag_i = 0; tag_i < tag_var_regs_per_cur_target_var->size(); tag_i++)
				{
					fprintf(f_op, "%d\n", tag_var_regs_per_cur_target_var->at(tag_i)->start);
				} // par_i loop.

				fprintf(f_op, "%d\n", target_var_regs->at(target_var_i)->start);
				fclose(f_op);

				// Test if available.
				if (testing_flag)
				{
					vector<double>* per_testing_sample_errors = new vector<double>();
					for (int testing_sample_i = 0; testing_sample_i < testing_tag_sample_ids->size(); testing_sample_i++)
					{
						// Reset the testing vector x: First element of x is 2.
						gsl_vector_set(testing_x, 0, 2);

						for (int reg_i = 0; reg_i < tag_var_regs_per_cur_target_var->size(); reg_i++)
						{
							void** cur_tag_reg_data = (void**)(tag_var_regs_per_cur_target_var->at(reg_i)->data);
							char* cur_testing_tag_genotypes = (char*)(cur_tag_reg_data[1]);
							if (cur_testing_tag_genotypes == NULL)
							{
								fprintf(stderr, "Testing Target region is null, exiting.\n");
								exit(0);
							}

							double cur_testing_sample_geno = (double)(cur_testing_tag_genotypes[testing_sample_i]);
							gsl_vector_set(testing_x, reg_i+1, cur_testing_sample_geno);	
						} // par_i loop.

						double geno_est = 0;
						double geno_err_est = 0;
						gsl_multifit_linear_est(testing_x, c, cov, &geno_est, &geno_err_est);

						// Estimate the real error.
						void** target_reg_data = (void**)(target_var_regs->at(target_var_i)->data);
						char* testing_target_var_genotypes = (char*)(target_reg_data[1]);
						if (testing_target_var_genotypes == NULL)
						{
							fprintf(stderr, "Testing Target region is null, exiting.\n");
							exit(0);
						}

						double cur_testing_sample_target_geno = (double)(testing_target_var_genotypes[testing_sample_i]);
						double real_err = fabs(cur_testing_sample_target_geno - geno_est);

						if (__DUMP_INPUTATION_UTILS_MSGS__)
						{
							fprintf(stderr, "%s:%d-%d (%s): Error: %lf (%lf)\n",
								target_var_regs->at(target_var_i)->chrom,
								target_var_regs->at(target_var_i)->start,
								target_var_regs->at(target_var_i)->end,
								target_var_regs->at(target_var_i)->name,
								real_err,
								geno_err_est);
						}

						per_testing_sample_errors->push_back(real_err);
					} // testing_sample_i loop.

					// Get the error statistics.
					double mean_err = 0;
					double std_dev_err = 0;
					get_stats(per_testing_sample_errors, mean_err, std_dev_err);

					fprintf(stderr, "%s:%d-%d (%s): Error: %lf (%lf)\n",
						target_var_regs->at(target_var_i)->chrom,
						target_var_regs->at(target_var_i)->start,
						target_var_regs->at(target_var_i)->end,
						target_var_regs->at(target_var_i)->name,
						mean_err,
						std_dev_err);
				} // testing check.

				//printf("# covariance matrix:\n");
				//printf("[ %+.5e, %+.5e, %+.5e  \n",
				//	COV(0, 0), COV(0, 1), COV(0, 2));
				//printf("  %+.5e, %+.5e, %+.5e  \n",
				//	COV(1, 0), COV(1, 1), COV(1, 2));
				//printf("  %+.5e, %+.5e, %+.5e ]\n",
				//	COV(2, 0), COV(2, 1), COV(2, 2));
				//printf("# chisq = %g\n", chisq);
			} // tag snp count check.
		} // target_var_i loop.
	} // target_chr_i loop.

	// Free memory.
	gsl_matrix_free(X);
	gsl_vector_free(y);
	gsl_vector_free(w);
	gsl_vector_free(c);
	gsl_matrix_free(cov);
	
//#endif 
}

char** get_match_by_identifiers()
{
	char** match_by_identifiers = new char*[N_MATCH_BY_IDS + 2];
	match_by_identifiers[MATCH_BY_NAME] = t_string::copy_me_str("byName");
	match_by_identifiers[MATCH_BY_START_POSN] = t_string::copy_me_str("byStartPosn");
	return(match_by_identifiers);
}

int get_PR_stats_bin_index_per_score(double MAP_geno, 
										double min_score_threshold, 
										double max_score_threshold, 
										double score_delta)
{
	if (MAP_geno < min_score_threshold ||
		MAP_geno > max_score_threshold)
	{
		return(-1);
	}

	int bin_i = (int)(floor((MAP_geno - min_score_threshold) / score_delta));

	return(bin_i);
}

void compute_imputation_stats_per_5col_genotype_probs_w_PR_curve_info(char* imputed_chr_id,
	char* imputed_5col_genotype_probs_fp,
	char* imputed_5col_genotype_probs_sample_ids_list_fp,
	char* known_genotype_regs_fp,
	char* known_genotype_regs_sample_ids_list_fp,
	char* match_by_str,
	double min_score_threshold,
	double max_score_threshold,
	double score_delta,
	char* op_stats_fp)
{
	fprintf(stderr, "Computing imputation statistics with PR information using [%.2f:%.3f:%.2f]\n", 
			min_score_threshold, score_delta, max_score_threshold);

	char** match_by_identifiers = get_match_by_identifiers();
	fprintf(stderr, "Matching by %s\n", match_by_str);

	vector<char*>* imputed_genotype_sample_ids = buffer_file(imputed_5col_genotype_probs_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d sample id's for the imputed genotypes.\n", imputed_genotype_sample_ids->size());

	vector<char*>* known_genotype_sample_ids = buffer_file(known_genotype_regs_sample_ids_list_fp);
	fprintf(stderr, "Loading known genotypes for %d samples from %s\n", known_genotype_sample_ids->size(), known_genotype_regs_fp);
	vector<t_annot_region*>* known_genotype_regs = load_variant_genotype_signal_regions(known_genotype_regs_fp, known_genotype_sample_ids);
	fprintf(stderr, "Loaded %d known genotype regions.\n", known_genotype_regs->size());

	fprintf(stderr, "Mapping imputed sample id's.\n");
	vector<int>* per_imputed_geno_sample_known_geno_sample_id = new vector<int>();
	int n_matched_imp_sample_ids = 0;
	for (int imp_sample_i = 0; imp_sample_i < imputed_genotype_sample_ids->size(); imp_sample_i++)
	{
		int cur_imp_sample_i_in_known_sample_i = t_string::get_i_str(known_genotype_sample_ids, imputed_genotype_sample_ids->at(imp_sample_i));
		per_imputed_geno_sample_known_geno_sample_id->push_back(cur_imp_sample_i_in_known_sample_i);

		if (cur_imp_sample_i_in_known_sample_i < known_genotype_sample_ids->size())
		{
			n_matched_imp_sample_ids++;
		}
	} // imp_sample_i loop.
	fprintf(stderr, "Mapped %d sample id's.\n", n_matched_imp_sample_ids);

	// Re-organize the known variant regions.
	t_restr_annot_region_list* restr_known_genotype_regs = restructure_annot_regions(known_genotype_regs);

	int i_imputed_chr = t_string::get_i_str(restr_known_genotype_regs->chr_ids, imputed_chr_id);
	fprintf(stderr, "Found the imputed chromosome id @ %d.\n", i_imputed_chr);

	FILE* f_5col_geno_probs = open_f(imputed_5col_genotype_probs_fp, "r");
	if (f_5col_geno_probs == NULL)
	{
		fprintf(stderr, "Could not open %s\n", imputed_5col_genotype_probs_fp);
		exit(0);
	}

	vector<char*>* sorted_region_names = new vector<char*>();
	if (t_string::compare_strings(match_by_str, match_by_identifiers[MATCH_BY_NAME]))
	{
		sort(restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->begin(), restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->end(), sort_genes_regions_per_name);

		for (int i_reg = 0; i_reg < restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->size(); i_reg++)
		{
			sorted_region_names->push_back(restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(i_reg)->name);
		} // i_reg loop.
	} // match by name check.

	  // Set the accuracy statistics for known genotypes.
	enum { GENO_CNT_MATCH, GENO_CNT_NON_REF, GENO_CNT_NON_REF_MATCH, GENO_CNT_N_KNOWN_IMPUTED_SAMPLE_MATCH, GENO_CNT_N_ENTRIES };

	// Allocate the PR curve information.
	int n_PR_bins = get_PR_stats_bin_index_per_score(max_score_threshold, min_score_threshold, max_score_threshold, score_delta) + 1;
	for (int i_reg = 0; i_reg < known_genotype_regs->size(); i_reg++)
	{
		int* cur_reg_acc_stats = new int[GENO_CNT_N_ENTRIES + 1];
		memset(cur_reg_acc_stats, 0, sizeof(int) * GENO_CNT_N_ENTRIES);
		void** geno_reg_info = (void**)(known_genotype_regs->at(i_reg)->data);

		// Replace the genomic info.
		void** ext_reg_info = new void*[10];
		known_genotype_regs->at(i_reg)->data = ext_reg_info;
		ext_reg_info[0] = geno_reg_info[0];

		ext_reg_info[1] = cur_reg_acc_stats;

		// All variant non-ref PR stats.
		double** cur_reg_non_ref_PR_stats = new double*[3];
		cur_reg_non_ref_PR_stats[0] = new double[n_PR_bins + 2];
		memset(cur_reg_non_ref_PR_stats[0], 0, sizeof(double) * n_PR_bins);
		cur_reg_non_ref_PR_stats[1] = new double[n_PR_bins + 2];
		memset(cur_reg_non_ref_PR_stats[1], 0, sizeof(double) * n_PR_bins);
		cur_reg_non_ref_PR_stats[2] = new double[n_PR_bins + 2];
		memset(cur_reg_non_ref_PR_stats[2], 0, sizeof(double) * n_PR_bins);
		ext_reg_info[2] = cur_reg_non_ref_PR_stats;

		// All variant PR stats.
		double** cur_reg_PR_stats = new double*[3];
		cur_reg_PR_stats[0] = new double[n_PR_bins + 2];
		memset(cur_reg_PR_stats[0], 0, sizeof(double) * n_PR_bins);
		cur_reg_PR_stats[1] = new double[n_PR_bins + 2];
		memset(cur_reg_PR_stats[1], 0, sizeof(double) * n_PR_bins);
		cur_reg_PR_stats[2] = new double[n_PR_bins + 2];
		memset(cur_reg_PR_stats[2], 0, sizeof(double) * n_PR_bins);
		ext_reg_info[3] = cur_reg_PR_stats;
	} // i_reg loop.

	int n_processed_vars = 0;
	double cumul_total_non_ref = 0;
	double cumul_matching_non_ref = 0;

	// Allocate the confusion matrix.
	double** confusion_matrix = new double*[3];	
	for (int geno = 0; geno < 3; geno++)
	{
		confusion_matrix[geno] = new double[3];
		memset(confusion_matrix[geno], 0, sizeof(double) * 3);
	} // geno loop.

	double*** per_bin_confusion_matrices = new double**[n_PR_bins];
	for (int i_bin = 0; i_bin < n_PR_bins; i_bin++)
	{
		per_bin_confusion_matrices[i_bin] = new double*[3];
		for (int geno = 0; geno < 3; geno++)
		{
			per_bin_confusion_matrices[i_bin][geno] = new double[3];
			memset(per_bin_confusion_matrices[i_bin][geno], 0, sizeof(double) * 3);
		} // geno loop.
	} // i_bin loop.

	while (1)
	{
		char* cur_imputed_var_line = getline(f_5col_geno_probs);
		if (cur_imputed_var_line == NULL)
		{
			break;
		}

		n_processed_vars++;

		if (n_processed_vars % 100000 == 0)
		{
			fprintf(stderr, "Parsing %d. imputed variant: %s  (%.3f)                      \r",
				n_processed_vars,
				cur_imputed_var_line,
				cumul_matching_non_ref / cumul_total_non_ref);
		}

		// Replace the spaces.
		t_string::replace_avoid_list(cur_imputed_var_line, " ", 0);

		// Parse the position, etc.
		/*
		Subject ID, target SNP, 0, 1, 2
		0, 17084716, -0.0763668, 0.0396188, -0.0321174
		*/
		t_string_tokens* toks = t_string::tokenize_by_chars(cur_imputed_var_line, ",");
		if (toks->size() != 5)
		{
			fprintf(stderr, "Could not parse 5 entries from: %s\n", cur_imputed_var_line);
			exit(0);
		}

		int imp_sample_i = atoi(toks->at(0)->str());
		char MAP_geno = -1;
		double MAP_geno_prob = -1000 * 1000 * 100;

		for (char geno = 0; geno < 3; geno++)
		{
			double cur_geno_prob = atof(toks->at(2 + geno)->str());
			if (cur_geno_prob > MAP_geno_prob)
			{
				MAP_geno_prob = cur_geno_prob;
				MAP_geno = geno;
			}
		} // geno loop.

		int match_reg_i = 0;
		bool found_match = false;
		if (t_string::compare_strings_ci(match_by_str, match_by_identifiers[MATCH_BY_START_POSN]))
		{
			int posn = atoi(toks->at(1)->str());

			if (__DUMP_INPUTATION_UTILS_MSGS__)
			{
				fprintf(stderr, "Parsed: %d: %s (MAP geno: %d -- %.3f)\n", posn, cur_imputed_var_line, MAP_geno, MAP_geno_prob);
			}

			// Identify the matching known genotype region for the current imputed variant.
			match_reg_i = locate_posn_region_per_region_starts(posn, restr_known_genotype_regs->regions_per_chrom[i_imputed_chr], 0, restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->size());
			while (match_reg_i > 0 &&
				(restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(match_reg_i)->start == posn ||
					restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(match_reg_i)->start > posn))
			{
				match_reg_i--;
			} // match_reg_i loop.

			while (match_reg_i < restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->size() &&
				(restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(match_reg_i)->start < posn ||
					restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(match_reg_i)->start == posn))
			{
				if (restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(match_reg_i)->start == posn)
				{
					found_match = true;
					break;
				}

				match_reg_i++;
			} // match_reg_i loop.
		}
		else if (t_string::compare_strings_ci(match_by_str, match_by_identifiers[MATCH_BY_NAME]))
		{
			char* cur_var_name = toks->at(1)->str();

			// Identify the matching known genotype region for the current imputed variant.
			match_reg_i = t_string::fast_search_string_per_prefix(cur_var_name, sorted_region_names, 0, sorted_region_names->size());
			while (match_reg_i > 0 &&
				(t_string::compare_strings(restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(match_reg_i)->name, cur_var_name) ||
					t_string::sort_strings_per_prefix(cur_var_name, restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(match_reg_i)->name)))
			{
				match_reg_i--;
			} // match_reg_i loop.

			while (match_reg_i < restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->size() &&
				(t_string::sort_strings_per_prefix(restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(match_reg_i)->name, cur_var_name) ||
					t_string::compare_strings(restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(match_reg_i)->name, cur_var_name)))
			{
				if (t_string::compare_strings(restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(match_reg_i)->name, cur_var_name))
				{
					found_match = true;
					break;
				} // position matching check.

				match_reg_i++;
			} // match_reg_i loop.

			if (!found_match)
			{
				fprintf(stderr, "Could not match %s\n", cur_var_name);
			}
		} // matchBy check.

		  // Match find check.
		if (found_match)
		{
			void** cur_known_var_reg_info = (void**)(restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(match_reg_i)->data);
			char* cur_known_var_genotypes = (char*)(cur_known_var_reg_info[0]);
			int* cur_known_var_acc_stats = (int*)(cur_known_var_reg_info[1]);
			double** cur_known_var_non_ref_PR_stats = (double**)(cur_known_var_reg_info[2]);
			double** cur_known_var_PR_stats = (double**)(cur_known_var_reg_info[3]);

			if (__DUMP_INPUTATION_UTILS_MSGS__)
			{
				fprintf(stderr, "Found the known genotype: %s (MAP geno: %d -- %.3f; Known genotype: %d)\n",
					cur_imputed_var_line, MAP_geno, MAP_geno_prob,
					cur_known_var_genotypes[per_imputed_geno_sample_known_geno_sample_id->at(imp_sample_i)]);
			}

			cur_known_var_acc_stats[GENO_CNT_N_KNOWN_IMPUTED_SAMPLE_MATCH]++;

			// Non-ref update.
			if (cur_known_var_genotypes[per_imputed_geno_sample_known_geno_sample_id->at(imp_sample_i)] != 0)
			{
				cumul_total_non_ref++;
				cur_known_var_acc_stats[GENO_CNT_NON_REF]++;
			}

			if (cur_known_var_genotypes[per_imputed_geno_sample_known_geno_sample_id->at(imp_sample_i)] == MAP_geno)
			{
				// Update the all match statistics.
				cur_known_var_acc_stats[GENO_CNT_MATCH]++;

				// Update the non-ref statistics.
				if (cur_known_var_genotypes[per_imputed_geno_sample_known_geno_sample_id->at(imp_sample_i)] != 0)
				{
					cumul_matching_non_ref++;
					cur_known_var_acc_stats[GENO_CNT_NON_REF_MATCH]++;
				}
			} // genotype matching check.

			// Update the confusion matrix of the whole data.
			int cur_known_geno = (int)(cur_known_var_genotypes[per_imputed_geno_sample_known_geno_sample_id->at(imp_sample_i)]);
			confusion_matrix[(int)MAP_geno][cur_known_geno]++;

			bool UPDATE_ALL_PR_STATS = true;
			if (UPDATE_ALL_PR_STATS)
			{
				// Update the all PR stats.
				int cur_PR_bin_i = get_PR_stats_bin_index_per_score(MAP_geno_prob, min_score_threshold, max_score_threshold, score_delta);
				if (cur_PR_bin_i == -1)
				{
					fprintf(stderr, "Illegal bin index: %s\n", cur_imputed_var_line);
					exit(0);
				}

				if (__DUMP_INPUTATION_UTILS_MSGS__)
					fprintf(stderr, "%s: Score %.4f; PR_bin_i: %d\n", cur_imputed_var_line, MAP_geno_prob, cur_PR_bin_i);

				for (int bin_i = 0; bin_i < cur_PR_bin_i; bin_i++)
				{
					per_bin_confusion_matrices[bin_i][(int)MAP_geno][cur_known_geno]++;
				}

				// Update the base: All the non-ref values. This is constant among all the score thresholds.
				for (int bin_i = 0; bin_i < n_PR_bins; bin_i++)
				{
					cur_known_var_PR_stats[2][bin_i]++;
				} // bin_i loop.

				// Go over the thresholds below this.
				for (int bin_i = 0; bin_i < cur_PR_bin_i; bin_i++)
				{
					// This is the denominator for the thresholded cases.
					cur_known_var_PR_stats[1][bin_i]++;

					if (cur_known_var_genotypes[per_imputed_geno_sample_known_geno_sample_id->at(imp_sample_i)] == MAP_geno)
					{
						cur_known_var_PR_stats[0][bin_i]++;
					}
				} // bin_i loop.
			} // All PR statistics check.

			// Update non-ref PR stats.
			bool UPDATE_NON_REF_PR_STATS = true;
			if (UPDATE_NON_REF_PR_STATS)
			{
				if (cur_known_var_genotypes[per_imputed_geno_sample_known_geno_sample_id->at(imp_sample_i)] != 0)
				{
					int cur_PR_bin_i = get_PR_stats_bin_index_per_score(MAP_geno_prob, min_score_threshold, max_score_threshold, score_delta);
					if (cur_PR_bin_i == -1)
					{
						fprintf(stderr, "Illegal bin index: %s\n", cur_imputed_var_line);
						exit(0);
					}

					if (__DUMP_INPUTATION_UTILS_MSGS__)
						fprintf(stderr, "%s: Score %.4f; PR_bin_i: %d\n", cur_imputed_var_line, MAP_geno_prob, cur_PR_bin_i);

					// Update the base: All the non-ref values. This is constant among all the score thresholds.
					for (int bin_i = 0; bin_i < n_PR_bins; bin_i++)
					{
						cur_known_var_non_ref_PR_stats[2][bin_i]++;
					} // bin_i loop.

					// Go over the thresholds below this.
					for (int bin_i = 0; bin_i < cur_PR_bin_i; bin_i++)
					{
						// This is the denominator for the thresholded cases.
						cur_known_var_non_ref_PR_stats[1][bin_i]++;

						if (cur_known_var_genotypes[per_imputed_geno_sample_known_geno_sample_id->at(imp_sample_i)] == MAP_geno)
						{
							cur_known_var_non_ref_PR_stats[0][bin_i]++;
						}
					} // bin_i loop.
				} // non-ref check.
			} // update non-ref stats check.
		} // match search check.

		// FRee memory.
		t_string::clean_tokens(toks);
		delete[] cur_imputed_var_line;
	} // imputed genotype reading loop.

	// Save file.
	FILE* f_op_stats = open_f(op_stats_fp, "w");

	//char non_ref_PR_stats_op_fp[1000];
	//sprintf(non_ref_PR_stats_op_fp, "%s_non_ref_PR_stats.txt", op_stats_fp);
	//FILE* f_non_ref_PR_stats_op = open_f(non_ref_PR_stats_op_fp, "w");

	double** micro_agg_non_ref_PR_stats = new double*[3];
	micro_agg_non_ref_PR_stats[0] = new double[n_PR_bins + 2];
	memset(micro_agg_non_ref_PR_stats[0], 0, sizeof(double) * n_PR_bins);
	micro_agg_non_ref_PR_stats[1] = new double[n_PR_bins + 2];
	memset(micro_agg_non_ref_PR_stats[1], 0, sizeof(double) * n_PR_bins);
	micro_agg_non_ref_PR_stats[2] = new double[n_PR_bins + 2];
	memset(micro_agg_non_ref_PR_stats[2], 0, sizeof(double) * n_PR_bins);

	double** micro_agg_PR_stats = new double*[3];
	micro_agg_PR_stats[0] = new double[n_PR_bins + 2];
	memset(micro_agg_PR_stats[0], 0, sizeof(double) * n_PR_bins);
	micro_agg_PR_stats[1] = new double[n_PR_bins + 2];
	memset(micro_agg_PR_stats[1], 0, sizeof(double) * n_PR_bins);
	micro_agg_PR_stats[2] = new double[n_PR_bins + 2];
	memset(micro_agg_PR_stats[2], 0, sizeof(double) * n_PR_bins);

	for (int i_reg = 0; i_reg < restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->size(); i_reg++)
	{
		void** cur_known_var_reg_info = (void**)(restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(i_reg)->data);
		int* cur_known_var_acc_stats = (int*)(cur_known_var_reg_info[1]);
		double** cur_known_var_non_ref_PR_stats = (double**)(cur_known_var_reg_info[2]);
		double** cur_known_var_PR_stats = (double**)(cur_known_var_reg_info[3]);

		if (cur_known_var_acc_stats[GENO_CNT_N_KNOWN_IMPUTED_SAMPLE_MATCH] > 0)
		{
			// Save the accuracy statistics.
			fprintf(f_op_stats, "%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\n",
				restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(i_reg)->chrom,
				translate_coord(restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
				translate_coord(restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
				restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(i_reg)->name,
				cur_known_var_acc_stats[GENO_CNT_MATCH],
				cur_known_var_acc_stats[GENO_CNT_N_KNOWN_IMPUTED_SAMPLE_MATCH],
				cur_known_var_acc_stats[GENO_CNT_NON_REF], cur_known_var_acc_stats[GENO_CNT_NON_REF_MATCH]);

			//// Save the PR statistics.
			//fprintf(f_PR_stats_op, "%s\t%d\t%d\t%s",
			//		restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(i_reg)->chrom,
			//		translate_coord(restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
			//		translate_coord(restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
			//		restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(i_reg)->name);

			// Update the pooled PR curve statistics.
			for (int i_bin = 0; i_bin < n_PR_bins; i_bin++)
			{
				//fprintf(f_PR_stats_op, "\t%.2f %.2f %.2f",
				//		cur_known_var_PR_stats[0][i_bin], 
				//		cur_known_var_PR_stats[1][i_bin], 
				//		cur_known_var_PR_stats[2][i_bin]);

				micro_agg_PR_stats[0][i_bin] += cur_known_var_PR_stats[0][i_bin];
				micro_agg_PR_stats[1][i_bin] += cur_known_var_PR_stats[1][i_bin];
				micro_agg_PR_stats[2][i_bin] += cur_known_var_PR_stats[2][i_bin];

				micro_agg_non_ref_PR_stats[0][i_bin] += cur_known_var_non_ref_PR_stats[0][i_bin];
				micro_agg_non_ref_PR_stats[1][i_bin] += cur_known_var_non_ref_PR_stats[1][i_bin];
				micro_agg_non_ref_PR_stats[2][i_bin] += cur_known_var_non_ref_PR_stats[2][i_bin];
			} // i_bin loop.

			//fprintf(f_PR_stats_op, "\n");
		}
	} // i_reg loop.
	close_f(f_op_stats, op_stats_fp);
	//close_f(f_PR_stats_op, PR_stats_op_fp);

	// Macro-accuracy aggregation per threshold.
	fprintf(stderr, "Macro-aggregating sensitivity and PPV for changing thresholds.\n");
	double* macro_agg_per_bin_PR_sens = new double[n_PR_bins + 2];
	memset(macro_agg_per_bin_PR_sens, 0, sizeof(double) * n_PR_bins);
	double* macro_agg_per_bin_PR_PPV = new double[n_PR_bins + 2];
	memset(macro_agg_per_bin_PR_PPV, 0, sizeof(double) * n_PR_bins);

	double* macro_agg_per_bin_non_ref_PR_sens = new double[n_PR_bins + 2];
	memset(macro_agg_per_bin_non_ref_PR_sens, 0, sizeof(double) * n_PR_bins);
	double* macro_agg_per_bin_non_ref_PR_PPV = new double[n_PR_bins + 2];
	memset(macro_agg_per_bin_non_ref_PR_PPV, 0, sizeof(double) * n_PR_bins);

	// Macro aggregation of the variants and the PR curve?
	char macro_agg_sens_PPV_stats_fp[1000];
	sprintf(macro_agg_sens_PPV_stats_fp, "%s_macro_agg_all_sens_PPV.txt", op_stats_fp);
	FILE* f_macro_agg_sens_PPV_stats = open_f(macro_agg_sens_PPV_stats_fp, "w");

	char macro_agg_non_ref_sens_PPV_stats_fp[1000];
	sprintf(macro_agg_non_ref_sens_PPV_stats_fp, "%s_macro_agg_non_ref_sens_PPV.txt", op_stats_fp);
	FILE* f_macro_agg_non_ref_sens_PPV_stats = open_f(macro_agg_non_ref_sens_PPV_stats_fp, "w");
	for (int i_bin = 0; i_bin < n_PR_bins; i_bin++)
	{
		int n_agg_all_vars = 0;
		int n_agg_non_ref_vars = 0;

		for (int i_reg = 0; i_reg < restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->size(); i_reg++)
		{
			void** cur_known_var_reg_info = (void**)(restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(i_reg)->data);
			int* cur_known_var_acc_stats = (int*)(cur_known_var_reg_info[1]);
			double** cur_known_var_non_ref_PR_stats = (double**)(cur_known_var_reg_info[2]);
			double** cur_known_var_PR_stats = (double**)(cur_known_var_reg_info[3]);

			if (cur_known_var_acc_stats[GENO_CNT_N_KNOWN_IMPUTED_SAMPLE_MATCH] > 0)
			{
				if (cur_known_var_PR_stats[2][i_bin] > 0 && cur_known_var_PR_stats[1][i_bin] > 0)
				{
					double cur_var_all_sens = cur_known_var_PR_stats[0][i_bin] / cur_known_var_PR_stats[2][i_bin];
					double cur_var_all_PPV = cur_known_var_PR_stats[0][i_bin] / cur_known_var_PR_stats[1][i_bin];
					macro_agg_per_bin_PR_sens[i_bin] += cur_var_all_sens;
					macro_agg_per_bin_PR_PPV[i_bin] += cur_var_all_PPV;
					n_agg_all_vars++;
				}

				if (cur_known_var_non_ref_PR_stats[2][i_bin] > 0 && cur_known_var_non_ref_PR_stats[1][i_bin] > 0)
				{
					double cur_var_non_ref_sens = cur_known_var_non_ref_PR_stats[0][i_bin] / cur_known_var_non_ref_PR_stats[2][i_bin];
					double cur_var_non_ref_PPV = cur_known_var_non_ref_PR_stats[0][i_bin] / cur_known_var_non_ref_PR_stats[1][i_bin];
					macro_agg_per_bin_non_ref_PR_sens[i_bin] += cur_var_non_ref_sens;
					macro_agg_per_bin_non_ref_PR_PPV[i_bin] += cur_var_non_ref_PPV;
					n_agg_non_ref_vars++;
				}
			} // variant existence check.
		} // i_reg loop.

		// Take the mean.
		if (n_agg_all_vars > 0)
		{
			macro_agg_per_bin_PR_sens[i_bin] /= n_agg_all_vars;
			macro_agg_per_bin_PR_PPV[i_bin] /= n_agg_all_vars;
		}

		if (n_agg_non_ref_vars > 0)
		{
			macro_agg_per_bin_non_ref_PR_sens[i_bin] /= n_agg_non_ref_vars;
			macro_agg_per_bin_non_ref_PR_PPV[i_bin] /= n_agg_non_ref_vars;
		}

		fprintf(f_macro_agg_sens_PPV_stats, "%lf\t%lf\n", macro_agg_per_bin_PR_sens[i_bin], macro_agg_per_bin_PR_PPV[i_bin]);
		fprintf(f_macro_agg_non_ref_sens_PPV_stats, "%lf\t%lf\n", macro_agg_per_bin_non_ref_PR_sens[i_bin], macro_agg_per_bin_non_ref_PR_PPV[i_bin]);
	} // l_bin loop.
	fclose(f_macro_agg_sens_PPV_stats);
	fclose(f_macro_agg_non_ref_sens_PPV_stats);

	fprintf(stderr, "\nDone. Saving accuracy stat files.\n");

	// Save the aggregate PR stats.
	char micro_agg_PR_stats_op_fp[1000];
	sprintf(micro_agg_PR_stats_op_fp, "%s_micro_agg_PR_stats.txt", op_stats_fp);
	FILE* f_micro_agg_PR_stats_op = open_f(micro_agg_PR_stats_op_fp, "w");

	char micro_agg_non_ref_PR_stats_op_fp[1000];
	sprintf(micro_agg_non_ref_PR_stats_op_fp, "%s_micro_agg_non_ref_PR_stats.txt", op_stats_fp);
	FILE* f_micro_agg_non_ref_PR_stats_op = open_f(micro_agg_non_ref_PR_stats_op_fp, "w");

	// Save the aggregated PR statistics.
	for (int i_bin = 0; i_bin < n_PR_bins; i_bin++)
	{
		fprintf(f_micro_agg_PR_stats_op, "%.2f %.2f %.2f\n",
				micro_agg_PR_stats[0][i_bin], 
				micro_agg_PR_stats[1][i_bin],
				micro_agg_PR_stats[2][i_bin]);

		fprintf(f_micro_agg_non_ref_PR_stats_op, "%.2f %.2f %.2f\n",
				micro_agg_non_ref_PR_stats[0][i_bin],
				micro_agg_non_ref_PR_stats[1][i_bin],
				micro_agg_non_ref_PR_stats[2][i_bin]);
	} // i_bin loop.

	close_f(f_micro_agg_PR_stats_op, micro_agg_PR_stats_op_fp);
	close_f(f_micro_agg_non_ref_PR_stats_op, micro_agg_non_ref_PR_stats_op_fp);

	// Save the confusion matrix.
	char confusion_matrix_fp[1000];
	sprintf(confusion_matrix_fp, "%s_confusion_matrix.txt", op_stats_fp);
	FILE* f_confusion_matrix = open_f(confusion_matrix_fp, "w");
	for (int geno_i = 0; geno_i < 3; geno_i++)
	{
		for (int geno_j = 0; geno_j < 3; geno_j++)
		{
			fprintf(f_confusion_matrix, "%d\t%d\t%.1f\n", geno_i, geno_j, confusion_matrix[geno_i][geno_j]);
		} // geno_j loop.
	} // geno_i loop.
	close_f(f_confusion_matrix, confusion_matrix_fp);

	// Save the per bin counted confusion matrices.
	char per_bin_confusion_matrix_fp[1000];
	sprintf(per_bin_confusion_matrix_fp, "%s_per_bin_confusion_matrix.txt", op_stats_fp);
	FILE* f_per_bin_confusion_matrix = open_f(per_bin_confusion_matrix_fp, "w");
	for (int i_bin = 0; i_bin < n_PR_bins; i_bin++)
	{
		for (int geno_i = 0; geno_i < 3; geno_i++)
		{
			for (int geno_j = 0; geno_j < 3; geno_j++)
			{
				fprintf(f_per_bin_confusion_matrix, "%d\t%d\t%d\t%.1f\n", i_bin, geno_i, geno_j, per_bin_confusion_matrices[i_bin][geno_i][geno_j]);
			} // geno_j loop.
		} // geno_i loop.
	} // i_bin loop.
	fclose(f_per_bin_confusion_matrix);
}

void compute_imputation_stats_per_5col_genotype_probs(char* imputed_chr_id,
														char* imputed_5col_genotype_probs_fp,
														char* imputed_5col_genotype_probs_sample_ids_list_fp,
														char* known_genotype_regs_fp,
														char* known_genotype_regs_sample_ids_list_fp,
														char* match_by_str,
														char* op_stats_fp)
{
	char** match_by_identifiers = get_match_by_identifiers();
	fprintf(stderr, "Matching by %s\n", match_by_str);

	vector<char*>* imputed_genotype_sample_ids = buffer_file(imputed_5col_genotype_probs_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d sample id's for the imputed genotypes.\n", imputed_genotype_sample_ids->size());

	vector<char*>* known_genotype_sample_ids = buffer_file(known_genotype_regs_sample_ids_list_fp);
	fprintf(stderr, "Loading known genotypes for %d samples from %s\n", known_genotype_sample_ids->size(), known_genotype_regs_fp);
	vector<t_annot_region*>* known_genotype_regs = load_variant_genotype_signal_regions(known_genotype_regs_fp, known_genotype_sample_ids);
	fprintf(stderr, "Loaded %d known genotype regions.\n", known_genotype_regs->size());

	fprintf(stderr, "Mapping imputed sample id's.\n");
	vector<int>* per_imputed_geno_sample_known_geno_sample_id = new vector<int>();
	int n_matched_imp_sample_ids = 0;
	for (int imp_sample_i = 0; imp_sample_i < imputed_genotype_sample_ids->size(); imp_sample_i++)
	{
		int cur_imp_sample_i_in_known_sample_i = t_string::get_i_str(known_genotype_sample_ids, imputed_genotype_sample_ids->at(imp_sample_i));
		per_imputed_geno_sample_known_geno_sample_id->push_back(cur_imp_sample_i_in_known_sample_i);

		if (cur_imp_sample_i_in_known_sample_i < known_genotype_sample_ids->size())
		{
			n_matched_imp_sample_ids++;
		}
	} // imp_sample_i loop.
	fprintf(stderr, "Mapped %d sample id's.\n", n_matched_imp_sample_ids);

	// Re-organize the known variant regions.
	t_restr_annot_region_list* restr_known_genotype_regs = restructure_annot_regions(known_genotype_regs);

	int i_imputed_chr = t_string::get_i_str(restr_known_genotype_regs->chr_ids, imputed_chr_id);
	fprintf(stderr, "Found the imputed chromosome id @ %d.\n", i_imputed_chr);

	FILE* f_5col_geno_probs = open_f(imputed_5col_genotype_probs_fp, "r");
	if (f_5col_geno_probs == NULL)
	{
		fprintf(stderr, "Could not open %s\n", imputed_5col_genotype_probs_fp);
		exit(0);
	}

	vector<char*>* sorted_region_names = new vector<char*>();
	if (t_string::compare_strings(match_by_str, match_by_identifiers[MATCH_BY_NAME]))
	{
		sort(restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->begin(), restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->end(), sort_genes_regions_per_name);

		for (int i_reg = 0; i_reg < restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->size(); i_reg++)
		{
			sorted_region_names->push_back(restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(i_reg)->name);
		} // i_reg loop.
	} // match by name check.

	// Set the accuracy statistics for known genotypes.
	enum{GENO_CNT_MATCH, GENO_CNT_NON_REF, GENO_CNT_NON_REF_MATCH, GENO_CNT_N_KNOWN_IMPUTED_SAMPLE_MATCH, GENO_CNT_N_ENTRIES};
	for (int i_reg = 0; i_reg < known_genotype_regs->size(); i_reg++)
	{
		int* cur_reg_acc_stats = new int[GENO_CNT_N_ENTRIES + 1];
		memset(cur_reg_acc_stats, 0, sizeof(int) * GENO_CNT_N_ENTRIES);
		void** geno_reg_info = (void**)(known_genotype_regs->at(i_reg)->data);
		geno_reg_info[1] = cur_reg_acc_stats;
	} // i_reg loop.

	//// Allocate the PR curve information.
	//double min_threshold = -10;
	//double max_threshold = 10;
	//double delta = 0.001;
	//double*  


	int n_processed_vars = 0;
	double cumul_total_non_ref = 0;
	double cumul_matching_non_ref = 0;
	while (1)
	{
		char* cur_imputed_var_line = getline(f_5col_geno_probs);
		if (cur_imputed_var_line == NULL)
		{
			break;
		}

		n_processed_vars++;

		if (n_processed_vars % 100000 == 0)
		{
			fprintf(stderr, "Parsing %d. imputed variant: %s  (%.3f)                      \r", 
					n_processed_vars, 
					cur_imputed_var_line, 
					cumul_matching_non_ref / cumul_total_non_ref);
		}

		// Replace the spaces.
		t_string::replace_avoid_list(cur_imputed_var_line, " ", 0);

		// Parse the position, etc.
		/*
		Subject ID, target SNP, 0, 1, 2
		0, 17084716, -0.0763668, 0.0396188, -0.0321174
		*/
		t_string_tokens* toks = t_string::tokenize_by_chars(cur_imputed_var_line, ",");
		if (toks->size() != 5)
		{
			fprintf(stderr, "Could not parse 5 entries from: %s\n", cur_imputed_var_line);
			exit(0);
		}

		int imp_sample_i = atoi(toks->at(0)->str());	
		char MAP_geno = -1;
		double MAP_geno_prob = -1;
		for (char geno = 0; geno < 3; geno++)
		{
			double cur_geno_prob = atof(toks->at(2 + geno)->str());
			if (cur_geno_prob > MAP_geno_prob)
			{
				MAP_geno_prob = cur_geno_prob;
				MAP_geno = geno;
			}
		} // geno loop.

		int match_reg_i = 0;
		bool found_match = false;
		if (t_string::compare_strings_ci(match_by_str, match_by_identifiers[MATCH_BY_START_POSN]))
		{
			int posn = atoi(toks->at(1)->str());

			if (__DUMP_INPUTATION_UTILS_MSGS__)
			{
				fprintf(stderr, "Parsed: %d: %s (MAP geno: %d -- %.3f)\n", posn, cur_imputed_var_line, MAP_geno, MAP_geno_prob);
			}

			// Identify the matching known genotype region for the current imputed variant.
			match_reg_i = locate_posn_region_per_region_starts(posn, restr_known_genotype_regs->regions_per_chrom[i_imputed_chr], 0, restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->size());
			while (match_reg_i > 0 &&
				(restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(match_reg_i)->start == posn ||
					restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(match_reg_i)->start > posn))
			{
				match_reg_i--;
			} // match_reg_i loop.

			while (match_reg_i < restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->size() &&
				(restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(match_reg_i)->start < posn ||
					restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(match_reg_i)->start == posn))
			{
				if (restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(match_reg_i)->start == posn)
				{
					found_match = true;
					break;
				}

				match_reg_i++;
			} // match_reg_i loop.
		}
		else if (t_string::compare_strings_ci(match_by_str, match_by_identifiers[MATCH_BY_NAME]))
		{
			char* cur_var_name = toks->at(1)->str();

			// Identify the matching known genotype region for the current imputed variant.
			match_reg_i = t_string::fast_search_string_per_prefix(cur_var_name, sorted_region_names, 0, sorted_region_names->size());
			while (match_reg_i > 0 &&
					(t_string::compare_strings(restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(match_reg_i)->name, cur_var_name) ||
					t_string::sort_strings_per_prefix(cur_var_name, restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(match_reg_i)->name)))
			{
				match_reg_i--;
			} // match_reg_i loop.

			while (match_reg_i < restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->size() &&
					(t_string::sort_strings_per_prefix(restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(match_reg_i)->name, cur_var_name) ||
					t_string::compare_strings(restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(match_reg_i)->name, cur_var_name)))
			{
				if (t_string::compare_strings(restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(match_reg_i)->name, cur_var_name))
				{
					found_match = true;
					break;
				} // position matching check.

				match_reg_i++;
			} // match_reg_i loop.

			if (!found_match)
			{
				fprintf(stderr, "Could not match %s\n", cur_var_name);
			}
		} // matchBy check.

		// Match find check.
		if (found_match)
		{
			void** cur_known_var_reg_info = (void**)(restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(match_reg_i)->data);
			char* cur_known_var_genotypes = (char*)(cur_known_var_reg_info[0]);
			int* cur_known_var_acc_stats = (int*)(cur_known_var_reg_info[1]);

			if (__DUMP_INPUTATION_UTILS_MSGS__)
			{
				fprintf(stderr, "Found the known genotype: %s (MAP geno: %d -- %.3f; Known genotype: %d)\n",
					cur_imputed_var_line, MAP_geno, MAP_geno_prob,
					cur_known_var_genotypes[per_imputed_geno_sample_known_geno_sample_id->at(imp_sample_i)]);
			}

			cur_known_var_acc_stats[GENO_CNT_N_KNOWN_IMPUTED_SAMPLE_MATCH]++;

			// Non-ref update.
			if (cur_known_var_genotypes[per_imputed_geno_sample_known_geno_sample_id->at(imp_sample_i)] != 0)
			{
				cumul_total_non_ref++;
				cur_known_var_acc_stats[GENO_CNT_NON_REF]++;
			}

			if (cur_known_var_genotypes[per_imputed_geno_sample_known_geno_sample_id->at(imp_sample_i)] == MAP_geno)
			{
				// Update the all match statistics.
				cur_known_var_acc_stats[GENO_CNT_MATCH]++;

				// Update the non-ref statistics.
				if (cur_known_var_genotypes[per_imputed_geno_sample_known_geno_sample_id->at(imp_sample_i)] != 0)
				{
					cumul_matching_non_ref++;
					cur_known_var_acc_stats[GENO_CNT_NON_REF_MATCH]++;
				}
			} // genotype matching check.
		}

		// FRee memory.
		t_string::clean_tokens(toks);
		delete[] cur_imputed_var_line;
	} // imputed genotype reading loop.

	// Save file.
	FILE* f_op_stats = open_f(op_stats_fp, "w");
	for (int i_reg = 0; i_reg < restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->size(); i_reg++)
	{
		void** cur_known_var_reg_info = (void**)(restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(i_reg)->data);
		int* cur_known_var_acc_stats = (int*)(cur_known_var_reg_info[1]);

		if (cur_known_var_acc_stats[GENO_CNT_N_KNOWN_IMPUTED_SAMPLE_MATCH] > 0)
		{
			fprintf(f_op_stats, "%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\n",
				restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(i_reg)->chrom,
				translate_coord(restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(i_reg)->start, CODEBASE_COORDS::start_base, BED_COORDS::start_base),
				translate_coord(restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(i_reg)->end, CODEBASE_COORDS::end_base, BED_COORDS::end_base),
				restr_known_genotype_regs->regions_per_chrom[i_imputed_chr]->at(i_reg)->name,
				cur_known_var_acc_stats[GENO_CNT_MATCH],
				cur_known_var_acc_stats[GENO_CNT_N_KNOWN_IMPUTED_SAMPLE_MATCH],
				cur_known_var_acc_stats[GENO_CNT_NON_REF], cur_known_var_acc_stats[GENO_CNT_NON_REF_MATCH]);
		}
	} // i_reg loop.
	close_f(f_op_stats, op_stats_fp);
}

void compute_imputation_stats_per_IMPUTED_genotypes(char* imputed_genotyped_regs_fp, char* imputed_genotyped_regs_sample_ids_list_fp,
	char* known_genotype_regs_fp, char* known_genotype_regs_sample_ids_list_fp)
{
	vector<char*>* imputed_genotype_sample_ids = buffer_file(imputed_genotyped_regs_sample_ids_list_fp);
	fprintf(stderr, "Loading imputed genotypes for %d samples from %s\n", imputed_genotype_sample_ids->size(), imputed_genotyped_regs_fp);
	vector<t_annot_region*>* imputed_genotype_regs = load_variant_genotype_signal_regions(imputed_genotyped_regs_fp, imputed_genotype_sample_ids);
	fprintf(stderr, "Loaded %d imputed genotype regions.\n", imputed_genotype_regs->size());

	vector<char*>* known_genotype_sample_ids = buffer_file(known_genotype_regs_sample_ids_list_fp);
	fprintf(stderr, "Loading known genotypes for %d samples from %s\n", known_genotype_sample_ids->size(), known_genotype_regs_fp);
	vector<t_annot_region*>* known_genotype_regs = load_variant_genotype_signal_regions(known_genotype_regs_fp, known_genotype_sample_ids);
	fprintf(stderr, "Loaded %d known genotype regions.\n", known_genotype_regs->size());

	fprintf(stderr, "Mapping imputed sample id's.\n");
	vector<int>* per_imputed_geno_sample_known_geno_sample_id = new vector<int>();
	int n_matched_imp_sample_ids = 0;
	for (int imp_sample_i = 0; imp_sample_i < imputed_genotype_sample_ids->size(); imp_sample_i++)
	{
		int cur_imp_sample_i_in_known_sample_i = t_string::get_i_str(known_genotype_sample_ids, imputed_genotype_sample_ids->at(imp_sample_i));
		per_imputed_geno_sample_known_geno_sample_id->push_back(cur_imp_sample_i_in_known_sample_i);

		if (cur_imp_sample_i_in_known_sample_i < known_genotype_sample_ids->size())
		{
			n_matched_imp_sample_ids++;
		}
	} // imp_sample_i loop.
	fprintf(stderr, "Mapped %d sample id's.\n", n_matched_imp_sample_ids);

	fprintf(stderr, "Intersecting known and imputed genotypes.\n");
	vector<t_annot_region*>* intersects = intersect_annot_regions(imputed_genotype_regs, known_genotype_regs, true);
	fprintf(stderr, "Processing %d intersects.\n", intersects->size());

	for (int i_int = 0; i_int < intersects->size(); i_int++)
	{
		t_intersect_info* cur_int_info = (t_intersect_info*)(intersects->at(i_int)->data);

		if (cur_int_info->src_reg->start == cur_int_info->dest_reg->start &&
			cur_int_info->src_reg->end == cur_int_info->dest_reg->end)
		{
			t_annot_region* imputed_reg = cur_int_info->src_reg;
			t_annot_region* known_reg = cur_int_info->dest_reg;
			void** imputed_reg_dat = (void**)(imputed_reg->data);
			imputed_reg_dat[1] = known_reg;
		}
	} // i_int loop.

	  // Analyze the intersects.
	char* mapped_imputed_geno_sig = new char[imputed_genotype_sample_ids->size() + 10];
	memset(mapped_imputed_geno_sig, 0, sizeof(char) * imputed_genotype_sample_ids->size());
	char* mapped_known_geno_sig = new char[imputed_genotype_sample_ids->size() + 10];
	memset(mapped_known_geno_sig, 0, sizeof(char) * imputed_genotype_sample_ids->size());
	FILE* f_op_stats = open_f("imputation_stats.txt", "w");
	for (int i_reg = 0; i_reg < imputed_genotype_regs->size(); i_reg++)
	{
		void** imputed_reg_dat = (void**)(imputed_genotype_regs->at(i_reg)->data);
		char* imputed_reg_geno_sig = (char*)(imputed_reg_dat[0]);
		if (imputed_reg_dat[1] != NULL)
		{
			t_annot_region* cur_imp_reg_known_reg = (t_annot_region*)(imputed_reg_dat[1]);
			void** known_reg_dat = (void**)(cur_imp_reg_known_reg->data);
			char* known_reg_geno_sig = (char*)(known_reg_dat[0]);

			int mapped_sample_i = 0;
			int n_total_matches = 0;
			int n_non_ref_matches = 0;
			int n_non_ref_known_geno = 0;
			int n_non_ref_imputed_geno = 0;
			for (int imp_sample_i = 0; imp_sample_i < imputed_genotype_sample_ids->size(); imp_sample_i++)
			{
				if (per_imputed_geno_sample_known_geno_sample_id->at(imp_sample_i) < known_genotype_sample_ids->size())
				{
					mapped_imputed_geno_sig[mapped_sample_i] = imputed_reg_geno_sig[imp_sample_i];
					mapped_known_geno_sig[mapped_sample_i] = known_reg_geno_sig[per_imputed_geno_sample_known_geno_sample_id->at(imp_sample_i)];

					// We should not have any haplocoded genotypes at this point.
					if (mapped_known_geno_sig[mapped_sample_i] == 3 ||
						mapped_imputed_geno_sig[mapped_sample_i] == 3)
					{
						fprintf(stderr, "Invalid genotype!\n");
						exit(0);
					}

					if (mapped_imputed_geno_sig[mapped_sample_i] > 0)
					{
						n_non_ref_imputed_geno++;
					}

					// Check total consistency.
					if (mapped_known_geno_sig[mapped_sample_i] == mapped_imputed_geno_sig[mapped_sample_i])
					{
						n_total_matches++;
					}

					// Check non-ref consistency.
					if (mapped_known_geno_sig[mapped_sample_i] > 0)
					{
						n_non_ref_known_geno++;

						if (mapped_known_geno_sig[mapped_sample_i] == mapped_imputed_geno_sig[mapped_sample_i])
						{
							n_non_ref_matches++;
						}
					}

					mapped_sample_i++;
				}
			} // imp_sample_i loop.

			fprintf(f_op_stats, "%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\n",
					imputed_genotype_regs->at(i_reg)->chrom, 
					imputed_genotype_regs->at(i_reg)->start, 
					imputed_genotype_regs->at(i_reg)->name,
					n_total_matches, mapped_sample_i,
					n_non_ref_matches, n_non_ref_known_geno, n_non_ref_imputed_geno);
		}
	} // i_reg loop.
	fclose(f_op_stats);
}

/*
imputed_testing_genotype_signal_fp
imputed_training_genotype_signal_fp
known_training_genotype_BED_fp
*/
void extract_genotyped_values_per_continuous_genotypes_signal_w_adaptive_thresholds_per_training_signal(char* imputed_testing_geno_signal_fp,
	char* imputed_testing_geno_signal_sample_ids_list_fp,
	char* imputed_training_geno_signal_fp,
	char* imputed_training_geno_signal_sample_ids_list_fp,
	char* known_training_geno_BED_fp,
	char* known_training_geno_sample_ids_list_fp,
	bool is_signal_haplocoded,
	bool save_haplocoded,
	char* output_genotype_matrix_bed_fp,
	char* prob_5col_op_fp)
{
	fprintf(stderr, "Extracting genotypes and probabilities from imputed signals:\n%s (%s)\nusing imputed training genotype signals:\n%s (%s)\nby using known training genotype signals as reference: \n%s (%s)\n", 
			imputed_testing_geno_signal_fp, imputed_testing_geno_signal_sample_ids_list_fp,
			imputed_training_geno_signal_fp, imputed_training_geno_signal_sample_ids_list_fp,
			known_training_geno_BED_fp, known_training_geno_sample_ids_list_fp);

	// Load the imputed genotype signal for the testing data.
	vector<char*>* imputed_testing_geno_signal_sample_ids = buffer_file(imputed_testing_geno_signal_sample_ids_list_fp);
	fprintf(stderr, "%d individuals in the imputed testing genotype signal regions.\n", imputed_testing_geno_signal_sample_ids->size());

	vector<t_annot_region*>* imputed_testing_geno_signal_regs = load_variant_genotype_signal_regions(imputed_testing_geno_signal_fp, imputed_testing_geno_signal_sample_ids);
	fprintf(stderr, "Loaded %d imputed testing genotype signal regions.\n", imputed_testing_geno_signal_regs->size());

	// Load the imputed genotype signal for the training data.
	vector<char*>* imputed_training_geno_signal_sample_ids = buffer_file(imputed_training_geno_signal_sample_ids_list_fp);
	fprintf(stderr, "%d individuals in the imputed training genotype signal regions.\n", imputed_training_geno_signal_sample_ids->size());

	vector<t_annot_region*>* imputed_training_geno_signal_regs = load_variant_genotype_signal_regions(imputed_training_geno_signal_fp, imputed_training_geno_signal_sample_ids);
	fprintf(stderr, "Loaded %d imputed training genotype signal regions.\n", imputed_training_geno_signal_regs->size());

	// Load the known genotypes from the training data.
	vector<char*>* known_training_geno_regs_sample_ids = buffer_file(known_training_geno_sample_ids_list_fp);
	fprintf(stderr, "%d individuals in the training genotype regions.\n", known_training_geno_regs_sample_ids->size());

	vector<t_annot_region*>* known_training_geno_regs = load_variant_genotype_signal_regions(known_training_geno_BED_fp, known_training_geno_regs_sample_ids);
	fprintf(stderr, "Loaded %d known training genotype regions.\n", known_training_geno_regs->size());

	if (known_training_geno_regs_sample_ids->size() != imputed_training_geno_signal_sample_ids->size())
	{
		fprintf(stderr, "Could not match the training sample ids: %d, %d\n", 
				known_training_geno_regs_sample_ids->size(), imputed_training_geno_signal_sample_ids->size());
		exit(0);
	}
	else
	{
		for (int i_s = 0; i_s < known_training_geno_regs_sample_ids->size(); i_s++)
		{
			if (!t_string::compare_strings(known_training_geno_regs_sample_ids->at(i_s), imputed_training_geno_signal_sample_ids->at(i_s)))
			{
				fprintf(stderr, "Training genotype and signal samples are not matching.\n");
				exit(0);
			}
		} // i_s loop.
	}

	fprintf(stderr, "Intersecting testing and training signal variants.\n");
	vector<t_annot_region*>* intersects = intersect_annot_regions(imputed_testing_geno_signal_regs, known_training_geno_regs, true);
	fprintf(stderr, "Processing %d intersects\n", intersects->size());
	int n_intersects = 0;
	for (int i_int = 0; i_int < intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* imputed_test_reg = int_info->src_reg;
		t_annot_region* known_geno_reg = int_info->dest_reg;

		if (t_string::compare_strings(imputed_test_reg->name, known_geno_reg->name))
		{
			void** test_info = (void**)(imputed_test_reg->data);
			if (test_info[1] != NULL)
			{
				fprintf(stderr, "%s already has a training region assigned.\n", imputed_test_reg->name);
				exit(0);
			}

			test_info[1] = known_geno_reg;
			
			n_intersects++;
		}
	} // i_int loop.
	fprintf(stderr, "Assigned %d training genotype regions to %d imputed testing regions.\n", n_intersects, imputed_testing_geno_signal_regs->size());

	fprintf(stderr, "Intersecting training genotype regions with training signal regions.\n");
	intersects = intersect_annot_regions(known_training_geno_regs, imputed_training_geno_signal_regs, true);
	n_intersects = 0;
	for (int i_int = 0; i_int < intersects->size(); i_int++)
	{
		t_intersect_info* int_info = (t_intersect_info*)(intersects->at(i_int)->data);
		t_annot_region* known_geno_reg = int_info->src_reg;
		t_annot_region* imputed_sig_reg = int_info->dest_reg;

		if (t_string::compare_strings(known_geno_reg->name, imputed_sig_reg->name))
		{
			void** geno_info = (void**)(known_geno_reg->data);
			if (geno_info[1] != NULL)
			{
				fprintf(stderr, "%s already has a training region assigned.\n", known_geno_reg->name);
				exit(0);
			}

			geno_info[1] = imputed_sig_reg;

			n_intersects++;
		}
	} // i_int loop.
	fprintf(stderr, "Assigned %d training genotype regions to %d training signal regions.\n", n_intersects, known_training_geno_regs->size());

	int max_geno = 2;
	if (is_signal_haplocoded)
	{
		fprintf(stderr, "Input signal is haplocoded.\n");
		max_geno = 3;
	}

	if (save_haplocoded)
	{
		fprintf(stderr, "Saving genotyped values as haplocoded values.\n");
	}
	else
	{
		fprintf(stderr, "Saving genotyped values as genocoded values.\n");
	}

	fprintf(stderr, "Rounding the genotype signals.\n");
	FILE* f_prob_5col_op = open_f(prob_5col_op_fp, "w");
	for (int i_reg = 0; i_reg < imputed_testing_geno_signal_regs->size(); i_reg++)
	{
		void** cur_reg_dat = (void**)(imputed_testing_geno_signal_regs->at(i_reg)->data);
		signed char* cur_reg_imputed_geno_sigs = (signed char*)(cur_reg_dat[0]);
		t_annot_region* known_train_geno_reg = (t_annot_region*)(cur_reg_dat[1]);

		void** train_geno_reg_info = (void**)(known_train_geno_reg->data);
		char* known_train_geno = (char*)(train_geno_reg_info[0]);
		t_annot_region* imputed_train_geno_sig_reg = (t_annot_region*)(train_geno_reg_info[1]);
		
		void** train_geno_sig_reg_info = (void**)(imputed_train_geno_sig_reg->data);
		signed char* imputed_train_geno_sigs = (signed char*)(train_geno_sig_reg_info[0]);

		// Find Center-of-mass for 0 and 1 then 1 and 2?
		vector<double>** per_geno_model_sigs = new vector<double>*[max_geno + 2];
		for (int geno = 0; geno <= max_geno; geno++)
		{
			per_geno_model_sigs[geno] = new vector<double>();
		} // geno loop.

		// Set the per genotype values and the maximum predicted genotype signal.
		double max_geno_sig = -1000;
		for (int i_s = 0; i_s < known_training_geno_regs_sample_ids->size(); i_s++)
		{
			per_geno_model_sigs[known_train_geno[i_s]]->push_back((double)(imputed_train_geno_sigs[i_s]));
			if (max_geno_sig < (double)(imputed_train_geno_sigs[i_s]))
			{
				max_geno_sig = (double)(imputed_train_geno_sigs[i_s]);
			}
		} // i_s loop.

		double* per_geno_means = new double[max_geno + 2];
		double* per_geno_std_devs = new double[max_geno + 2];
		for (int geno = 0; geno <= max_geno; geno++)
		{
			double cur_geno_mean, cur_geno_std_dev;
			get_stats(per_geno_model_sigs[geno], cur_geno_mean, cur_geno_std_dev);

			if (per_geno_model_sigs[geno]->size() > 1)
			{
				per_geno_means[geno] = cur_geno_mean;
				per_geno_std_devs[geno] = cur_geno_std_dev;
			}
			else
			{
				per_geno_means[geno] = ((double)geno / (double)max_geno) * 100;
				per_geno_std_devs[geno] = 0;
			}

			if (__DUMP_INPUTATION_UTILS_MSGS__)
			{
				fprintf(stderr, "Stats @ %s:%d (%s)::Geno: %d: %.2f (%.2f)\n",
					imputed_testing_geno_signal_regs->at(i_reg)->chrom, imputed_testing_geno_signal_regs->at(i_reg)->start, imputed_testing_geno_signal_regs->at(i_reg)->name,
					geno, per_geno_means[geno], per_geno_std_devs[geno]);
			}
		} // geno loop.

		// Set the geno-2-geno thresholds (maximum signal for each genotype); following does not set a threshold for the last genotype.
		double* per_geno_left2_right_cumul_thresh = new double[max_geno + 2];
		for (int geno = 0; geno < max_geno; geno++)
		{
			per_geno_left2_right_cumul_thresh[geno] = (per_geno_means[geno] + per_geno_means[geno + 1]) / 2;

			// Compute the false negatives for this genotype:
			int n_cur_geno_tot = 0;
			int n_cur_geno_false = 0;
			for (int i_s = 0; i_s < known_training_geno_regs_sample_ids->size(); i_s++)
			{
				if (known_train_geno[i_s] == geno)
				{
					n_cur_geno_tot++;

					bool false_impute_check = (imputed_train_geno_sigs[i_s] > per_geno_left2_right_cumul_thresh[geno]) ||
						(geno > 0 && imputed_train_geno_sigs[i_s] < per_geno_left2_right_cumul_thresh[geno - 1]);

					if (false_impute_check)
					{
						n_cur_geno_false++;
					}
				}
			} // i_s loop.

			if (__DUMP_INPUTATION_UTILS_MSGS__)
			{
				fprintf(stderr, "FN-rate of Geno %d\t%d\t%d\t%.3f\n", imputed_testing_geno_signal_regs->at(i_reg)->name, geno, n_cur_geno_false, n_cur_geno_tot, ((double)n_cur_geno_false) / n_cur_geno_tot);
			}

			// Reset the genotype 0 threshold.
			if(geno == 0)
			{
				double reset_geno0_thresh = -1;
				//for (int thresh = per_geno_means[geno]; thresh <= per_geno_means[geno + 1]; thresh++)
				for (int thresh = per_geno_means[geno]; thresh <= (geno + 1) * 50; thresh++)
				{
					n_cur_geno_tot = 0;
					n_cur_geno_false = 0;
					double n_cur_geno1_tot = 0;
					double n_cur_geno1_marked_geno0 = 0;
					double prev_fpr = 1.0;
					for (int i_s = 0; i_s < known_training_geno_regs_sample_ids->size(); i_s++)
					{
						if (known_train_geno[i_s] == geno)
						{
							n_cur_geno_tot++;

							bool false_impute_check = (imputed_train_geno_sigs[i_s] > thresh) ||
														(geno > 0 && imputed_train_geno_sigs[i_s] < per_geno_left2_right_cumul_thresh[geno - 1]);

							if (false_impute_check)
							{
								n_cur_geno_false++;
							}
						}
						else if (known_train_geno[i_s] == 1)
						{
							n_cur_geno1_tot++;

							bool false_impute_check = (imputed_train_geno_sigs[i_s] <= thresh);
							if (false_impute_check)
							{
								n_cur_geno1_marked_geno0++;
							}
						}
					} // i_s loop.

					if (reset_geno0_thresh == -1 &&
						((double)n_cur_geno_false) / n_cur_geno_tot < 0.01)
						//((double)n_cur_geno1_marked_geno0) / n_cur_geno1_tot < 0.05)
						//((double)n_cur_geno1_marked_geno0) / n_cur_geno1_tot <= prev_fpr)
					{
						if (__DUMP_INPUTATION_UTILS_MSGS__)
						{
							fprintf(stderr, "Reset FN optimized threshold @ %d\n", reset_geno0_thresh);
						}
						reset_geno0_thresh = thresh;
						prev_fpr = ((double)n_cur_geno1_marked_geno0) / n_cur_geno1_tot;
					}

					if (__DUMP_INPUTATION_UTILS_MSGS__)
					{
						fprintf(stderr, "%s: FN-Rate of Geno %d @ %d \t%d\t%d\t%.3f\t%d\t%d\t%.3f\n", 
							imputed_testing_geno_signal_regs->at(i_reg)->name, geno, thresh, n_cur_geno_false, n_cur_geno_tot, ((double)n_cur_geno_false) / n_cur_geno_tot,
							(int)n_cur_geno1_marked_geno0, (int)n_cur_geno1_tot, n_cur_geno1_marked_geno0 / n_cur_geno1_tot);
					}
				} // tresh loop.
				
				if (reset_geno0_thresh == -1)
				{
					reset_geno0_thresh = per_geno_left2_right_cumul_thresh[geno];
				}

				// Reset the threshold.
				if (per_geno_left2_right_cumul_thresh[geno] < reset_geno0_thresh)
				{
					per_geno_left2_right_cumul_thresh[geno] = reset_geno0_thresh;
				}
			} // false classification compute check.

			if (__DUMP_INPUTATION_UTILS_MSGS__)
			{
				fprintf(stderr, "Threshold @ %s:%d (%s):: %d [%.3f]\n",
					imputed_testing_geno_signal_regs->at(i_reg)->chrom, imputed_testing_geno_signal_regs->at(i_reg)->start, imputed_testing_geno_signal_regs->at(i_reg)->name,
					geno, per_geno_left2_right_cumul_thresh[geno]);
			}
		} // geno loop.
		per_geno_left2_right_cumul_thresh[1] = 75;
		per_geno_left2_right_cumul_thresh[max_geno] = 1000*1000;
		
		// Go over all the imputed testing region scores and identify the genotypes.
		for (int i_s = 0; i_s < imputed_testing_geno_signal_sample_ids->size(); i_s++)
		{
			// Multiply the probabilities: Note that the probabilities are stored after multiplcation by 100.
			//char cur_rounded_geno = (char)(round((max_geno * (double)cur_reg_geno_sigs[i_s]) / 100));
			
			char cur_rounded_geno = -1;
			for (int geno = 0; geno <= max_geno; geno++)
			{
				if (cur_reg_imputed_geno_sigs[i_s] <= per_geno_left2_right_cumul_thresh[geno])
				{
					cur_rounded_geno = geno;
					break;
				}

			} // geno loop.

			double per_geno_score[10];
			for (int geno = 0; geno <= max_geno; geno++)
			{
				per_geno_score[geno] = exp(-1 * fabs(per_geno_means[geno] - (double)(cur_reg_imputed_geno_sigs[i_s])));
			} // geno loop.

			fprintf(f_prob_5col_op, "%d,%d,%.4f,%.4f,%.4f\n", 
				i_s, imputed_testing_geno_signal_regs->at(i_reg)->start,
				per_geno_score[0], 
				per_geno_score[1], 
				per_geno_score[2]);

			if (__DUMP_INPUTATION_UTILS_MSGS__)
			{
				fprintf(stderr, "Sample %d: Genotype signal %d => %d\n", 
						i_s, 
						(int)(cur_reg_imputed_geno_sigs[i_s]),
						cur_rounded_geno);
			}

			if (cur_rounded_geno == -1)
			{
				fprintf(stderr, "Sanity check failed (%d), could not assign a genotype.\n", (int)(cur_reg_imputed_geno_sigs[i_s]));
				for (int geno = 0; geno <= max_geno; geno++)
				{
					fprintf(stderr, "Threshold %d: %.3f\n", geno, per_geno_left2_right_cumul_thresh[geno]);
				} // geno loop.
				exit(0);
			}

			// Convert to genocoded if haplocoded input is given.
			if (!save_haplocoded)
			{
				if (is_signal_haplocoded)
				{
					cur_rounded_geno = get_genotype_per_haplocoded_genotype((char)(cur_rounded_geno));
				}
			}

			cur_reg_imputed_geno_sigs[i_s] = cur_rounded_geno;
		} // i_s loop.

		delete[] per_geno_left2_right_cumul_thresh;
	} // i_reg loop.

	// Close the 5-column genotype probability file.
	close_f(f_prob_5col_op, prob_5col_op_fp);

	// Save the plain file.
	fprintf(stderr, "Saving genotyped signals.\n");

	dump_geno_sig_regs_plain(imputed_testing_geno_signal_regs, imputed_testing_geno_signal_sample_ids, false, output_genotype_matrix_bed_fp);
}

// This function converts the genotype regression model output to discrete genotypes.
void extract_genotyped_values_per_continuous_genotypes_signal_w_preset_thresholds(char* baseline_model_output_file_path,
																						char* input_genotype_matrix_sample_ids_list_fp,
																						bool is_signal_haplocoded,
																						bool save_haplocoded,
																						char* output_genotype_matrix_bed_fp)
{
	vector<char*>* genotyped_sample_ids = buffer_file(input_genotype_matrix_sample_ids_list_fp);
	fprintf(stderr, "%d individuals in the baseline model regressed genotyped data.\n", genotyped_sample_ids->size());

	vector<t_annot_region*>* baseline_genotyped_signal_regs = load_variant_genotype_signal_regions(baseline_model_output_file_path, genotyped_sample_ids);
	fprintf(stderr, "Loaded %d regions.\n", baseline_genotyped_signal_regs->size());

	double max_geno = 2;
	if (is_signal_haplocoded)
	{
		fprintf(stderr, "Input signal is haplocoded.\n");
		max_geno = 3;
	}

	if (save_haplocoded)
	{
		fprintf(stderr, "Saving genotyped values as haplocoded values.\n");
	}
	else
	{
		fprintf(stderr, "Saving genotyped values as genocoded values.\n");
	}

	fprintf(stderr, "Rounding the genotype signals.\n");
	for (int i_reg = 0; i_reg < baseline_genotyped_signal_regs->size(); i_reg++)
	{
		void** cur_reg_dat = (void**)(baseline_genotyped_signal_regs->at(i_reg)->data);
		char* cur_reg_geno_sigs = (char*)(cur_reg_dat[0]);

		for (int i_s = 0; i_s < genotyped_sample_ids->size(); i_s++)
		{
			// Multiply the probabilities: Note that the probabilities are stored after multiplcation by 100.
			char cur_rounded_geno = (char)(round((max_geno * (double)cur_reg_geno_sigs[i_s]) / 100));

			if (cur_rounded_geno > max_geno)
			{
				cur_rounded_geno = max_geno;
			}
			
			if (cur_rounded_geno < 0)
			{
				cur_rounded_geno = 0;
			}

			// Convert to genocoded if haplocoded input is given.
			if (!save_haplocoded)
			{
				if (is_signal_haplocoded)
				{
					cur_rounded_geno = get_genotype_per_haplocoded_genotype((char)(cur_rounded_geno));
				}
			}
			cur_reg_geno_sigs[i_s] = cur_rounded_geno;
		} // i_s loop.
	} // i_reg loop.

	// Save the plain file.
	fprintf(stderr, "Saving genotyped signals.\n");
	
	dump_geno_sig_regs_plain(baseline_genotyped_signal_regs, genotyped_sample_ids, false, output_genotype_matrix_bed_fp);
}

void extract_genotyped_values_per_IMPUTE2_probabilities_output(char* IMPUTE2_output_fp, char* input_genotype_matrix_sample_ids_list_fp, char* output_genotype_matrix_bed_fp,
	double min_prob_2_genotype,
	t_IMPUTE2_op_col_info* IMPUTE2_op_col_info,
	char* chrom)
{
	vector<char*>* genotyped_sample_ids = buffer_file(input_genotype_matrix_sample_ids_list_fp);
	fprintf(stderr, "%d individuals in the IMPUTE2 genotyped data.\n", genotyped_sample_ids->size());

	FILE* f_IMPUTE2_output = open_f(IMPUTE2_output_fp, "r");

	// Start reading each variant; write to the output file.
	vector<t_annot_region*>* impute2_genotyped_regs = new vector<t_annot_region*>();
	fprintf(stderr, "Selecting max prob genotypes and saving to %s..\n", output_genotype_matrix_bed_fp);
	while (1)
	{
		char* cur_line = getline(f_IMPUTE2_output);
		if (cur_line == NULL)
		{
			break;
		}

		t_string_tokens* cur_line_toks = t_string::tokenize_by_chars(cur_line, " ");

		if (cur_line_toks->size() != IMPUTE2_op_col_info->geno_prob_starting_0_based_col_i + 3 * genotyped_sample_ids->size())
		{
			fprintf(stderr, "%d columns in %s but %d expected columns (%d, %d)\n",
				cur_line_toks->size(),
				IMPUTE2_op_col_info->geno_prob_starting_0_based_col_i + 3 * genotyped_sample_ids->size(),
				IMPUTE2_op_col_info->geno_prob_starting_0_based_col_i, genotyped_sample_ids->size());

			exit(0);
		}

		// Get the variant position.
		char* cur_var_name = cur_line_toks->at(IMPUTE2_op_col_info->name_col_i)->str();
		int cur_var_posn = atoi(cur_line_toks->at(IMPUTE2_op_col_info->posn_col_i)->str());
		t_annot_region* reg = get_empty_region();
		reg->chrom = t_string::copy_me_str(chrom);
		reg->start = cur_var_posn;
		reg->end = cur_var_posn;
		reg->strand = '+';
		reg->name = t_string::copy_me_str(cur_var_name);
		void** cur_reg_dat = new void*[2];
		reg->data = cur_reg_dat;
		impute2_genotyped_regs->push_back(reg);

		int i_tok = IMPUTE2_op_col_info->geno_prob_starting_0_based_col_i;

		char* cur_var_IMPUTE2_genotypes = new char[genotyped_sample_ids->size() + 2];
		cur_reg_dat[0] = cur_var_IMPUTE2_genotypes;

		// Go over all the samples.
		for (int sample_i = 0; sample_i < genotyped_sample_ids->size(); sample_i++)
		{
			int max_prob_geno = 0;
			double max_prob_geno_prob = -1;
			for (int geno = 0; geno < 3; geno++)
			{
				if (atof(cur_line_toks->at(i_tok)->str()) > max_prob_geno_prob)
				{
					// Sanity check on the current token before parsing it.
					if (i_tok >= cur_line_toks->size())
					{
						fprintf(stderr, "Sanity check failed: Could not parse the genotypes before line ended.\n");
						exit(0);
					}

					max_prob_geno_prob = atof(cur_line_toks->at(i_tok)->str());
					max_prob_geno = geno;
				}

				//fprintf(stderr, "sample_i: %d; geno: %d, i_tok: %d\n", sample_i, geno, i_tok);

				i_tok++;
			} // geno loop.

			cur_var_IMPUTE2_genotypes[sample_i] = (char)(max_prob_geno);
		} // i_tok loop.

		if (i_tok != IMPUTE2_op_col_info->geno_prob_starting_0_based_col_i + 3 * genotyped_sample_ids->size())
		{
			fprintf(stderr, "Could not use all the entries: %d, %d\n", i_tok, IMPUTE2_op_col_info->geno_prob_starting_0_based_col_i + 3 * genotyped_sample_ids->size());
			exit(0);
		}

		t_string::clean_tokens(cur_line_toks);
	} // file reading loop.

	  // Close files.
	close_f(f_IMPUTE2_output, IMPUTE2_output_fp);

	// Save the plain file.
	dump_geno_sig_regs_plain(impute2_genotyped_regs, genotyped_sample_ids, false, output_genotype_matrix_bed_fp);
}

//mach - d sample.dat - p sample.ped - h hapmap.haplos - s hapmap.snps --rounds 500 --states 200 --weighted --geno
void extract_MACH_data_input_files_per_haplocoded_genotype_matrix(char* reference_haplotype_genotype_matrix_matbed_fp, char* reference_haplotype_genotype_matrix_sample_ids_list_fp,
	char* input_genotype_matrix_matbed_fp, char* input_genotype_matrix_sample_ids_list_fp,
	char* d_option_fp,
	char* p_option_fp,
	char* h_option_fp,
	char* s_option_fp)
{
	fprintf(stderr, "Saving the MACH files using the reference haplotype data @ %s (%s) and input sample genotype data @ %s (%s). Outputs:\n\
	-d %s\n\
	-p %s\n\
	-h %s\n\
	-s %s\n",
		reference_haplotype_genotype_matrix_matbed_fp, reference_haplotype_genotype_matrix_sample_ids_list_fp,
		input_genotype_matrix_matbed_fp, input_genotype_matrix_sample_ids_list_fp,
		d_option_fp,
		p_option_fp,
		h_option_fp,
		s_option_fp);

	vector<char*>* reference_haplo_sample_ids = buffer_file(reference_haplotype_genotype_matrix_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d reference haplotype sample id's.\n", reference_haplo_sample_ids->size());
	vector<t_annot_region*>* reference_haplo_geno_sig_regs = load_variant_genotype_signal_regions(reference_haplotype_genotype_matrix_matbed_fp, reference_haplo_sample_ids);
	fprintf(stderr, "Loaded %d reference haplotype variant regions.\n", (int)(reference_haplo_geno_sig_regs->size()));

	// Sort the regions first.
	t_restr_annot_region_list* restr_ref_panel_var_regs = restructure_annot_regions(reference_haplo_geno_sig_regs);

	// Saving the legend file.
	/*
	-s option (hapmap.snps):
	SNP1
	SNP3
	SNP7
	SNP8
	SNP9
	SNP11
	SNP13
	SNP14
	SNP18
	*/
	FILE* f_s_option = open_f(s_option_fp, "w");

	/*
	-h option (hapmap.haplos):
	H_0001->H_0001 HAPLO1 2221111222222212222122212222222221222222212222222122212222111222222122212212122211121211222222122222121122222112212121211112211111121212222222121221112122212221222212122211221221
	H_0001->H_0001 HAPLO2 2222222122121212221212121121212112211222221221112222222121222122211221121222211122222222222222222122121122222112212121211112211111121212222222121221122122212221222212122211221221
	H_0002->H_0002 HAPLO1 2221111222222212222122212222222221222222212222222122212222111222222122212212122211121211222122222111222221111221122212122221122222222122211111212112221222212221222111122211222122
	H_0002->H_0002 HAPLO2 1121212211222221212121222222222221222112222122222122112221222122211221121222211222222222221212212211222221111221122212122221122221121212222222121221122121122112111222211122121221
	H_0003->H_0003 HAPLO1 2221111222222212222122212222222221222222212222222122212222111222222122212212122211121211222222122211222221111221122212122221122222222122211111212112221222212221222212122211221221
	H_0003->H_0003 HAPLO2 2222222122211212221212121121212122212222121211122221221121222112211211121222211222222222222222121111222221111221122212122221122222222122211111212112221221122112111222211122122212
	H_0004->H_0004 HAPLO1 2222222122211212221212121121212122212221121211122221221121222112211211121222211222222222222222121122121122222122211121211122222221121211222222121221122122212221222212122211222212
	*/
	FILE* f_h_option = open_f(h_option_fp, "w");

	// Write the alleles for each sample.
	for (int i_s = 0; i_s < reference_haplo_sample_ids->size(); i_s++)
	{
		fprintf(stderr, "Processing reference panel sample %s (%d/%d)                       \r", reference_haplo_sample_ids->at(i_s), i_s, reference_haplo_sample_ids->size());

		for (int hap_i = 0; hap_i < 2; hap_i++)
		{
			for (int i_chr = 0; i_chr < restr_ref_panel_var_regs->chr_ids->size(); i_chr++)
			{				
				vector<t_annot_region*>* cur_chr_ref_panel_var_regs = restr_ref_panel_var_regs->regions_per_chrom[i_chr];

				fprintf(f_h_option, "%s->%s\tHAPLO%d\t", reference_haplo_sample_ids->at(i_s), reference_haplo_sample_ids->at(i_s), hap_i + 1);

				for (int i_reg = 0; i_reg < cur_chr_ref_panel_var_regs->size(); i_reg++)
				{
					// Write the legend entry.			
					t_string_tokens* toks = t_string::tokenize_by_chars(cur_chr_ref_panel_var_regs->at(i_reg)->name, "_");
					char* cur_var_name = toks->at(0)->str();
					char* cur_var_ref_str = toks->at(1)->str();
					char* cur_var_alt_str = toks->at(2)->str();

					// Save the current var name; if this is the first sample and first haplotype.
					if (i_s == 0 &&
						hap_i == 0)
					{
						fprintf(f_s_option, "%s\n", cur_var_name);
					}

					// Save the reference haplotype file (-h).
					void** cur_reg_info = (void**)(cur_chr_ref_panel_var_regs->at(i_reg)->data);
					char* per_sample_haplocoded_geno = (char*)(cur_reg_info[0]);

					int cur_hap = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], hap_i);

					fprintf(f_h_option, "%d", cur_hap+1);

					t_string::clean_tokens(toks);
				} // i_reg loop.
			} // chr_i loop.
			fprintf(f_h_option, "\n");
		} // hap_i loop.
	} // i_s loop.
	close_f(f_s_option, s_option_fp);
	close_f(f_h_option, h_option_fp);

	/*
	-d option (sample.dat):
	M SNP1
	M SNP7
	M SNP8
	M SNP11
	M SNP20
	M SNP23
	M SNP27
	M SNP33
	M SNP34
	*/

	/*
	-p option (sample.ped):
	S_0001  S_0001  0       0       2       2/1     2/2     2/1     2/1     2/2     2/2     2/2     1/1     2/2     2/2     2/2     2/2     2/2     2/2     2/2     1/2     1/1     2/2     2/2     2/2     2/2     2/2     2/2     1/2     1/2 2/1      2/2     1/2     1/2     2/2     2/1     2/1     2/2     2/2     1/2     2/1     2/2     2/2     2/2     2/2     2/2     2/2     2/2     1/1     2/2     2/2
	S_0002  S_0002  0       0       2       2/2     2/2     2/2     2/2     1/1     2/2     2/2     1/1     2/2     2/2     2/2     2/2     1/1     2/2     1/1     2/2     2/2     2/2     2/2     2/2     2/2     2/2     2/2     2/2     1/1 1/1      2/2     2/2     2/2     2/2     1/1     1/1     2/2     2/2     2/2     2/2     2/2     2/2     2/2     2/2     2/1     2/1     1/2     2/2     2/1     2/2
	S_0003  S_0003  0       0       2       2/2     2/2     2/1     2/1     2/2     2/2     2/2     1/2     2/1     2/2     2/2     2/2     2/2     2/2     2/2     1/2     1/2     2/2     2/2     2/1     2/2     2/2     2/2     1/2     1/2 1/1      2/2     2/2     2/2     2/2     1/1     1/1     2/2     2/2     2/2     2/2     2/2     1/2     1/2     1/2     2/2     2/2     2/2     1/1     2/2     2/2
	*/
	fprintf(stderr, "Processing study genotype information.\n");
	vector<char*>* input_geno_sample_ids = buffer_file(input_genotype_matrix_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d input genotype sample id's.\n", input_geno_sample_ids->size());
	vector<t_annot_region*>* input_geno_sig_regs = load_variant_genotype_signal_regions(input_genotype_matrix_matbed_fp, input_geno_sample_ids);
	fprintf(stderr, "Loaded %d input genotype variant regions.\n", (int)(input_geno_sig_regs->size()));

	t_restr_annot_region_list* restr_input_geno_sig_regs = restructure_annot_regions(input_geno_sig_regs);

	FILE* f_d_option = open_f(d_option_fp, "w");
	FILE* f_p_option = open_f(p_option_fp, "w");
	for (int i_s = 0; i_s < input_geno_sample_ids->size(); i_s++)
	{
		fprintf(stderr, "Processing input sample %s (%d/%d)                \r", input_geno_sample_ids->at(i_s), i_s, input_geno_sample_ids->size());

		// Write the alleles for each sample.
		//S_0001  S_0001  0       0       2
		fprintf(f_p_option, "%s\t%s\t0\t0\t2", input_geno_sample_ids->at(i_s), input_geno_sample_ids->at(i_s));

		for (int i_chr = 0; i_chr < restr_input_geno_sig_regs->chr_ids->size(); i_chr++)
		{
			vector<t_annot_region*>* cur_chr_input_geno_sig_regs = restr_input_geno_sig_regs->regions_per_chrom[i_chr];
			for (int i_reg = 0; i_reg < cur_chr_input_geno_sig_regs->size(); i_reg++)
			{
				// Write the legend entry.
				t_string_tokens* toks = t_string::tokenize_by_chars(cur_chr_input_geno_sig_regs->at(i_reg)->name, "_");
				char* cur_var_name = toks->at(0)->str();
				char* cur_var_ref_str = toks->at(1)->str();
				char* cur_var_alt_str = toks->at(2)->str();

				// Write the strand.
				if (i_s == 0)
				{
					fprintf(f_d_option, "M\t%s\n", cur_var_name);
				}

				// Save the reference haplotype file (-h).
				void** cur_reg_info = (void**)(cur_chr_input_geno_sig_regs->at(i_reg)->data);
				char* per_sample_haplocoded_geno = (char*)(cur_reg_info[0]);

				int cur_hap0 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 0);
				int cur_hap1 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 1);

				fprintf(f_p_option, "\t%d/%d", cur_hap0+1, cur_hap1+1);
			} // i_reg loop.	
		} // i_chr loop.

		fprintf(f_p_option, "\n");
	} // i_s loop. 

	close_f(f_p_option, p_option_fp);
	close_f(f_d_option, d_option_fp);
}

// ref and gt options must have the same number of variants.
void extract_BEAGLE_data_input_files_per_haplocoded_genotype_matrix(char* reference_haplotype_genotype_matrix_matbed_fp, char* reference_haplotype_genotype_matrix_sample_ids_list_fp,
	char* input_genotype_matrix_matbed_fp, char* input_genotype_matrix_sample_ids_list_fp,
	char* ref_option_fp,
	char* gt_option_fp,
	bool save_phased_gt_option)
{
	fprintf(stderr, "Saving the BEAGLE files using the reference haplotype data @ %s (%s) and input sample genotype data @ %s (%s). Outputs:\n\
	-ref %s\n\
	-gt %s\n",
		reference_haplotype_genotype_matrix_matbed_fp, reference_haplotype_genotype_matrix_sample_ids_list_fp,
		input_genotype_matrix_matbed_fp, input_genotype_matrix_sample_ids_list_fp,
		ref_option_fp, gt_option_fp);

	if (save_phased_gt_option)
	{
		fprintf(stderr, "**Saving phased gt option.**\n");
	}

	vector<char*>* reference_haplo_sample_ids = buffer_file(reference_haplotype_genotype_matrix_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d reference haplotype sample id's.\n", reference_haplo_sample_ids->size());
	vector<t_annot_region*>* reference_haplo_geno_sig_regs = load_variant_genotype_signal_regions(reference_haplotype_genotype_matrix_matbed_fp, reference_haplo_sample_ids);
	fprintf(stderr, "Loaded %d reference haplotype variant regions.\n", (int)(reference_haplo_geno_sig_regs->size()));

	// Sort the regions first.
	t_restr_annot_region_list* restr_ref_panel_var_regs = restructure_annot_regions(reference_haplo_geno_sig_regs);

	// Saving the legend file.
	FILE* f_ref_option = open_f(ref_option_fp, "w");
	fprintf(f_ref_option, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
	for (int i_s = 0; i_s < reference_haplo_sample_ids->size(); i_s++)
	{
		fprintf(f_ref_option, "\t%s", reference_haplo_sample_ids->at(i_s));
	} // i_s loop.
	fprintf(f_ref_option, "\n");

	for (int i_chr = 0; i_chr < restr_ref_panel_var_regs->chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Processing reference panel variants on %s\n", restr_ref_panel_var_regs->chr_ids->at(i_chr));

		vector<t_annot_region*>* cur_chr_ref_panel_var_regs = restr_ref_panel_var_regs->regions_per_chrom[i_chr];
		for (int i_reg = 0; i_reg < cur_chr_ref_panel_var_regs->size(); i_reg++)
		{
			// Write the legend entry.			
			t_string_tokens* toks = t_string::tokenize_by_chars(cur_chr_ref_panel_var_regs->at(i_reg)->name, "_");
			char* cur_var_name = toks->at(0)->str();
			char* cur_var_ref_str = toks->at(1)->str();
			char* cur_var_alt_str = toks->at(2)->str();

			// Save the reference haplotype file (-h).
			void** cur_reg_info = (void**)(cur_chr_ref_panel_var_regs->at(i_reg)->data);
			char* per_sample_haplocoded_geno = (char*)(cur_reg_info[0]);

			// Write the alleles for each sample.
			// 22      20000086        rs138720731     T       C       100     PASS    . GT
			fprintf(f_ref_option, "%s\t%d\t%s\t%s\t%s\t100\tPASS\t.\tGT",
				cur_chr_ref_panel_var_regs->at(i_reg)->chrom,
				cur_chr_ref_panel_var_regs->at(i_reg)->start,
				cur_var_name,
				cur_var_ref_str, cur_var_alt_str);

			// Write the alleles for each sample.
			for (int i_s = 0; i_s < reference_haplo_sample_ids->size(); i_s++)
			{
				int cur_hap0 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 0);
				int cur_hap1 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 1);

				// Setup the phased genotype string.
				char geno_str[10];
				strcpy(geno_str, "0|0");
				if(cur_hap0 == 1)
				{
					geno_str[0] = '1';
				}
				
				if (cur_hap1 == 1)
				{
					geno_str[2] = '1';
				}

				fprintf(f_ref_option, "\t%s", geno_str);
			} // i_s loop.

			fprintf(f_ref_option, "\n");
			t_string::clean_tokens(toks);
		} // i_reg loop.
	} // i_chr loop.
	close_f(f_ref_option, ref_option_fp);

	fprintf(stderr, "Processing study genotype information.\n");
	vector<char*>* input_geno_sample_ids = buffer_file(input_genotype_matrix_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d input genotype sample id's.\n", input_geno_sample_ids->size());
	vector<t_annot_region*>* input_geno_sig_regs = load_variant_genotype_signal_regions(input_genotype_matrix_matbed_fp, input_geno_sample_ids);
	fprintf(stderr, "Loaded %d input genotype variant regions.\n", (int)(input_geno_sig_regs->size()));

	// Reset the scores to 0 to indicate these are not assigned a ref region, yet.
	for (int i_reg = 0; i_reg < input_geno_sig_regs->size(); i_reg++)
	{
		input_geno_sig_regs->at(i_reg)->score = 0;
	} // i_reg loop.

	FILE* f_gt_option = open_f(gt_option_fp, "w");

	// Write the header.
	fprintf(f_gt_option, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
	for (int i_s = 0; i_s < input_geno_sample_ids->size(); i_s++)
	{
		fprintf(f_gt_option, "\t%s", input_geno_sample_ids->at(i_s));
	} // i_s loop.
	fprintf(f_gt_option, "\n");

	// Restructure the genotype signal regions.
	t_restr_annot_region_list* restr_input_geno_sig_regs = restructure_annot_regions(input_geno_sig_regs);

	// Write the genotypes.
	for (int i_chr = 0; i_chr < restr_input_geno_sig_regs->chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Processing input panel variants on %s\n", restr_input_geno_sig_regs->chr_ids->at(i_chr));

		vector<t_annot_region*>* cur_chr_input_geno_sig_regs = restr_input_geno_sig_regs->regions_per_chrom[i_chr];
		for (int i_reg = 0; i_reg < cur_chr_input_geno_sig_regs->size(); i_reg++)
		{
			t_string_tokens* toks = t_string::tokenize_by_chars(cur_chr_input_geno_sig_regs->at(i_reg)->name, "_");
			char* cur_var_name = toks->at(0)->str();
			char* cur_var_ref_str = toks->at(1)->str();
			char* cur_var_alt_str = toks->at(2)->str();

			// Apparently, BEAGLE does not want this.
			// Write the alleles for each sample.
			// 22      20000086        rs138720731     T       C       100     PASS    . GT
			fprintf(f_gt_option, "%s\t%d\t%s\t%s\t%s\t100\tPASS\t.\tGT",
				cur_chr_input_geno_sig_regs->at(i_reg)->chrom,
				cur_chr_input_geno_sig_regs->at(i_reg)->start,
				cur_var_name,
				cur_var_ref_str, cur_var_alt_str);

			// Save the reference haplotype file (-h).
			void** cur_reg_info = (void**)(cur_chr_input_geno_sig_regs->at(i_reg)->data);
			char* per_sample_haplocoded_geno = (char*)(cur_reg_info[0]);

			for (int i_s = 0; i_s < input_geno_sample_ids->size(); i_s++)
			{
				int cur_hap0 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 0);
				int cur_hap1 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 1);

				if(save_phased_gt_option)
				{ 
					fprintf(f_gt_option, "\t%d|%d", cur_hap0, cur_hap1);
				}
				else
				{
					int geno = cur_hap0 + cur_hap1;

					if (geno == 0)
					{
						fprintf(f_gt_option, "\t0/0");
					}
					else if (geno == 1)
					{
						fprintf(f_gt_option, "\t0/1");
					}
					else if (geno == 2)
					{
						fprintf(f_gt_option, "\t1/1");
					}
					else
					{
						fprintf(stderr, "Could not parse the genotype: %s:%d::%d (%d, %d)\n",
							cur_chr_input_geno_sig_regs->at(i_reg)->chrom, cur_chr_input_geno_sig_regs->at(i_reg)->start,
							(int)(per_sample_haplocoded_geno[i_s]), cur_hap0, cur_hap1);

						exit(0);
					}
				}

			} // i_s loop.

			fprintf(f_gt_option, "\n");

			t_string::clean_tokens(toks);
		} // i_reg loop.
	} // i_chr loop. 

	close_f(f_gt_option, gt_option_fp);
}

void extract_IMPUTE2_data_input_files_per_haplocoded_genotype_matrix(char* reference_haplotype_genotype_matrix_matbed_fp, char* reference_haplotype_genotype_matrix_sample_ids_list_fp,
	char* input_genotype_matrix_matbed_fp, char* input_genotype_matrix_sample_ids_list_fp,
	char* h_option_fp,
	char* l_option_fp,
	char* g_option_fp,
	char* strand_g_option_fp)
{
	fprintf(stderr, "Saving the IMPUTE2 files using the reference haplotype data @ %s (%s) and input sample genotype data @ %s (%s). Outputs:\n\
	-h %s\n\
	-l %s\n\
	-g %s\n\
	-strand_g %s\n",
		reference_haplotype_genotype_matrix_matbed_fp, reference_haplotype_genotype_matrix_sample_ids_list_fp,
		input_genotype_matrix_matbed_fp, input_genotype_matrix_sample_ids_list_fp,
		h_option_fp,
		l_option_fp,
		g_option_fp,
		strand_g_option_fp);

	/*
	./impute2 \
	-m ./Example/example.chr22.map :: Recombination map.
	-h ./Example/example.chr22.1kG.haps :: List of haplotypes: one haplotype per row: Generate this for the individuals only: Matches the legend (.legend) file.
	-l ./Example/example.chr22.1kG.legend :: rsID position a0 a1 :: rs4821114 20300810 G C
	-g ./Example/example.chr22.study.gens :: Known genotypes with positions: Must be sorted with respect to position.
	-strand_g ./Example/example.chr22.study.strand :: Strand of each variant in '-g' option.
	-int 20.4e6 20.5e6 :: The start and end positions.
	-Ne 20000 :: Do not change
	-o ./Example/example.chr22.one.phased.impute2 :: Output file.
	*/
	vector<char*>* reference_haplo_sample_ids = buffer_file(reference_haplotype_genotype_matrix_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d reference haplotype sample id's.\n", reference_haplo_sample_ids->size());
	vector<t_annot_region*>* reference_haplo_geno_sig_regs = load_variant_genotype_signal_regions(reference_haplotype_genotype_matrix_matbed_fp, reference_haplo_sample_ids);
	fprintf(stderr, "Loaded %d reference haplotype variant regions.\n", (int)(reference_haplo_geno_sig_regs->size()));

	// Sort the regions first.
	t_restr_annot_region_list* restr_ref_panel_var_regs = restructure_annot_regions(reference_haplo_geno_sig_regs);

	// Saving the legend file.
	FILE* f_l_option = open_f(l_option_fp, "w");
	FILE* f_h_option = open_f(h_option_fp, "w");

	fprintf(f_l_option, "rsID position a0 a1\n");
	for (int i_chr = 0; i_chr < restr_ref_panel_var_regs->chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Processing reference panel variants on %s\n", restr_ref_panel_var_regs->chr_ids->at(i_chr));

		vector<t_annot_region*>* cur_chr_ref_panel_var_regs = restr_ref_panel_var_regs->regions_per_chrom[i_chr];
		for (int i_reg = 0; i_reg < cur_chr_ref_panel_var_regs->size(); i_reg++)
		{
			// Write the legend entry.			
			t_string_tokens* toks = t_string::tokenize_by_chars(cur_chr_ref_panel_var_regs->at(i_reg)->name, "_");
			char* cur_var_name = toks->at(0)->str();
			char* cur_var_ref_str = toks->at(1)->str();
			char* cur_var_alt_str = toks->at(2)->str();

			fprintf(f_l_option, "%s %d %s %s\n",
					cur_var_name,
					cur_chr_ref_panel_var_regs->at(i_reg)->start,
					cur_var_ref_str,
					cur_var_alt_str);

			// Save the reference haplotype file (-h).
			void** cur_reg_info = (void**)(cur_chr_ref_panel_var_regs->at(i_reg)->data);
			char* per_sample_haplocoded_geno = (char*)(cur_reg_info[0]);

			// Write the alleles for each sample.
			for (int i_s = 0; i_s < reference_haplo_sample_ids->size(); i_s++)
			{
				int cur_hap0 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 0);
				int cur_hap1 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 1);

				if (i_s > 0)
				{
					fprintf(f_h_option, " ");
				}

				fprintf(f_h_option, "%d %d", cur_hap0, cur_hap1);
			} // i_s loop.

			fprintf(f_h_option, "\n");

			t_string::clean_tokens(toks);
		} // i_reg loop.
	} // i_chr loop.
	close_f(f_l_option, l_option_fp);
	close_f(f_h_option, h_option_fp);

	fprintf(stderr, "Processing study genotype information.\n");
	vector<char*>* input_geno_sample_ids = buffer_file(input_genotype_matrix_sample_ids_list_fp);
	fprintf(stderr, "Loaded %d input genotype sample id's.\n", input_geno_sample_ids->size());
	vector<t_annot_region*>* input_geno_sig_regs = load_variant_genotype_signal_regions(input_genotype_matrix_matbed_fp, input_geno_sample_ids);
	fprintf(stderr, "Loaded %d input genotype variant regions.\n", (int)(input_geno_sig_regs->size()));

	t_restr_annot_region_list* restr_input_geno_sig_regs = restructure_annot_regions(input_geno_sig_regs);

	FILE* f_g_option = open_f(g_option_fp, "w");
	FILE* f_strand_g_option = open_f(strand_g_option_fp, "w");
	for (int i_chr = 0; i_chr < restr_input_geno_sig_regs->chr_ids->size(); i_chr++)
	{
		fprintf(stderr, "Processing input panel variants on %s\n", restr_input_geno_sig_regs->chr_ids->at(i_chr));

		vector<t_annot_region*>* cur_chr_input_geno_sig_regs = restr_input_geno_sig_regs->regions_per_chrom[i_chr];
		for (int i_reg = 0; i_reg < cur_chr_input_geno_sig_regs->size(); i_reg++)
		{
			// Write the legend entry.
			t_string_tokens* toks = t_string::tokenize_by_chars(cur_chr_input_geno_sig_regs->at(i_reg)->name, "_");
			char* cur_var_name = toks->at(0)->str();
			char* cur_var_ref_str = toks->at(1)->str();
			char* cur_var_alt_str = toks->at(2)->str();

			// Write the strand.
			fprintf(f_strand_g_option, "%d +\n", cur_chr_input_geno_sig_regs->at(i_reg)->start);

			// Save the reference haplotype file (-h).
			void** cur_reg_info = (void**)(cur_chr_input_geno_sig_regs->at(i_reg)->data);
			char* per_sample_haplocoded_geno = (char*)(cur_reg_info[0]);

			// Write the alleles for each sample.
			fprintf(f_g_option, "%s %s %d %s %s",
				cur_var_name, cur_var_name,
				cur_chr_input_geno_sig_regs->at(i_reg)->start,
				cur_var_ref_str, cur_var_alt_str);

			for (int i_s = 0; i_s < input_geno_sample_ids->size(); i_s++)
			{
				int cur_hap0 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 0);
				int cur_hap1 = get_allele_per_haplotype(per_sample_haplocoded_geno[i_s], 1);

				int geno = cur_hap0 + cur_hap1;

				if (geno == 0)
				{
					fprintf(f_g_option, " 1 0 0");
				}
				else if (geno == 1)
				{
					fprintf(f_g_option, " 0 1 0");
				}
				else if (geno == 2)
				{
					fprintf(f_g_option, " 0 0 1");
				}
				else
				{
					fprintf(stderr, "Could not parse the genotype: %s:%d::%d (%d, %d)\n",
						cur_chr_input_geno_sig_regs->at(i_reg)->chrom, cur_chr_input_geno_sig_regs->at(i_reg)->start,
						(int)(per_sample_haplocoded_geno[i_s]), cur_hap0, cur_hap1);

					exit(0);
				}

			} // i_s loop.

			fprintf(f_g_option, "\n");
		} // i_reg loop.
	} // i_chr loop. 
	close_f(f_g_option, g_option_fp);
	close_f(f_strand_g_option, strand_g_option_fp);
}
