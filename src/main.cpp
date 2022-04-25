#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lhmmr_file_utils.h"
#include "lhmmr_gimp_utils.h"
#include "lhmmr_ansi_cli.h"
#include "lhmmr_variation_tools.h"
#include "lhmmr_annot_region_tools.h"
#include "lhmmr_ansi_string.h"
#include <string.h>
#include <time.h>

int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "USAGE: %s [options] [arguments]\n\
Options:\n\
Preprocessing Options:\n\
	-extract_genotype_signals_per_VCF\n\
	-extract_genotype_signals_per_region_list\n\
	-extract_genotype_signals_per_subsample_list\n\
	-binarize_genotype_signals_BED\n\
Imputation HMM:\n\
	-run_GIMP_sliding_window_Viterbi_State_Reduction\n\
	-run_GIMP_sliding_window_ForeBack_State_Reduction\n\
Accuracy Evaluation:\n\
	-get_R2_per_GIMP_4entry_allelic_probs\n\
	-get_R2_per_imputed_genotypes\n\
	-get_R2_per_imputed_genotypes_signal_level\n\
	-get_PR_stats_per_GIMP_4entry_allelic_probs\n\
	-get_PR_stats_per_3entry_genotype_probs\n", argv[0]);

		exit(0);
	}

	clock_t start_c = clock();

	if (strcmp(argv[1], "-get_R2_per_imputed_genotypes_signal_level") == 0)
	{
		if (argc < 6)
		{
			fprintf(stderr, "USAGE: %s %s --imputed_dir [Imputed genotypes directory] --known_dir [Known genotypes directory]\n", argv[0], argv[1]);
			exit(0);
		}

		t_ansi_cli* cli = new t_ansi_cli(argc, argv, "--");
		bool res = false;

		char* imputed_dir = t_string::copy_me_str(cli->get_value_by_option("--imputed_dir", res));
		if (!res)
		{
			fprintf(stderr, "Need to specify the output directory with \"--imputed_dir\"\n");
			exit(0);
		}

		char* known_dir = t_string::copy_me_str(cli->get_value_by_option("--known_dir", res));
		if (!res)
		{
			fprintf(stderr, "Need to specify the model parameters directory with \"--known_dir\"\n");
			exit(0);
		}

		char imp_samples_list_fp[1000];
		sprintf(imp_samples_list_fp, "%s/sample_ids.list", imputed_dir);
		if (!check_file(imp_samples_list_fp))
		{
			fprintf(stderr, "Could not find imputed sample id's @ %s\n", imp_samples_list_fp);
			exit(0);
		}

		char known_samples_list_fp[1000];
		sprintf(known_samples_list_fp, "%s/sample_ids.list", known_dir);
		if (!check_file(known_samples_list_fp))
		{
			fprintf(stderr, "Could not find known sample id's @ %s\n", known_samples_list_fp);
			exit(0);
		}

		char chr_ids_list_fp[1000];
		sprintf(chr_ids_list_fp, "%s/chr_ids.list", imputed_dir);

		if (!check_file(chr_ids_list_fp))
		{
			fprintf(stderr, "Could not fine the chromosome id's @ %s\n", chr_ids_list_fp);
			exit(0);
		}

		vector<char*>* chr_ids = buffer_file(chr_ids_list_fp);
		fprintf(stderr, "Loaded %d chromosomes.\n", chr_ids->size());

		// Process the chromosomes in order.		
		for (int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
		{
			char cur_chr_imp_geno_fp[1000];
			sprintf(cur_chr_imp_geno_fp, "%s/%s.imp", imputed_dir, chr_ids->at(i_chr));

			char cur_chr_known_geno_fp[1000];
			sprintf(cur_chr_known_geno_fp, "%s/%s.matbed", known_dir, chr_ids->at(i_chr));

			char stats_fp[1000];
			sprintf(stats_fp, "R2_stats_%s.txt", chr_ids->at(i_chr));
			get_R2_per_imputed_genotypes_signal_level(cur_chr_imp_geno_fp, imp_samples_list_fp,
				cur_chr_known_geno_fp, known_samples_list_fp, stats_fp);
		} // i_chr loop.
	} // -get_accuracy_statistics option.
	if (strcmp(argv[1], "-generate_reduced_state_haplotype_blocks") == 0)
	{
		if (argc != 5)
		{
			fprintf(stderr, "%s %s [Reference haplocoded genotypes file path] [Reference genotypes sample ids path] [# of variants per block]\n", argv[0], argv[1]);
			exit(0);
		}

		char* reference_haplocoded_tag_target_geno_regs_fp = argv[2];
		char* ref_sample_ids_list_fp = argv[3];
		int n_vars_per_block = atoi(argv[4]);

		generate_reduced_state_blocks_constant_size_blocks(reference_haplocoded_tag_target_geno_regs_fp,
			ref_sample_ids_list_fp, 0, n_vars_per_block);
	} // --generate_reduced_state_haplotype_blocks option.
	else if (strcmp(argv[1], "-get_PR_stats_per_GIMP_4entry_allelic_probs") == 0)
	{
		if (argc != 7)
		{
			fprintf(stderr, "%s %s [Imputed genotypes matrix file path] [Imputed genotypes sample ids path] \
[Known genotypes matrix bed file path] [Known genotypes samples path] [Imputed probs are linear? (>0 for linear)]\n", argv[0], argv[1]);
			exit(0);
		}

		char* imputed_genotypes_fp = argv[2];
		char* imputed_sample_ids_list_fp = argv[3];
		char* known_genotypes_fp = argv[4];
		char* known_sample_ids_list_fp = argv[5];
		bool imputed_probs_are_linear = (atoi(argv[6]) > 0);

		get_PR_stats_per_GIMP_4entry_allelic_probs(imputed_genotypes_fp, imputed_sample_ids_list_fp, known_genotypes_fp, known_sample_ids_list_fp, imputed_probs_are_linear);
	}
	else if (strcmp(argv[1], "-get_PR_stats_per_3entry_genotype_probs") == 0)
	{
		if (argc != 6)
		{
			fprintf(stderr, "%s %s [Imputed genotypes matrix file path] [Imputed genotypes sample ids path] \
[Known genotypes matrix bed file path] [Known genotypes samples path]\n", argv[0], argv[1]);
			exit(0);
		}

		char* imputed_genotypes_fp = argv[2];
		char* imputed_sample_ids_list_fp = argv[3];
		char* known_genotypes_fp = argv[4];
		char* known_sample_ids_list_fp = argv[5];

		get_PR_stats_per_3entry_genotype_probs(imputed_genotypes_fp, imputed_sample_ids_list_fp, known_genotypes_fp, known_sample_ids_list_fp);
	} // -get_PR_stats_per_3entry_genotype_probs option.
	else if (strcmp(argv[1], "-get_R2_per_imputed_genotypes") == 0)
	{
		if (argc != 6)
		{
			fprintf(stderr, "%s %s [Imputed genotypes matrix file path] [Imputed genotypes sample ids path] \
[Known genotypes matrix bed file path] [Known genotypes samples path]\n", argv[0], argv[1]);
			exit(0);
		}

		char* imputed_genotypes_fp = argv[2];
		char* imputed_sample_ids_list_fp = argv[3];
		char* known_genotypes_fp = argv[4];
		char* known_sample_ids_list_fp = argv[5];

		get_R2_per_imputed_genotypes(imputed_genotypes_fp, imputed_sample_ids_list_fp, known_genotypes_fp, known_sample_ids_list_fp);
	}
	else if (strcmp(argv[1], "-get_R2_per_GIMP_4entry_allelic_probs") == 0)
	{
		if (argc != 6)
		{
			fprintf(stderr, "%s %s [LoHaMMer/SHiMMer generated 4entry probabilities path] [Imputed genotypes sample ids path] \
[Known genotypes matrix bed file path] [Known genotypes samples path]\n", argv[0], argv[1]);
			exit(0);
		}

		char* imputed_genotypes_fp = argv[2];
		char* imputed_sample_ids_list_fp = argv[3];
		char* known_genotypes_fp = argv[4];
		char* known_sample_ids_list_fp = argv[5];

		get_R2_per_GIMP_4entry_allelic_probs(imputed_genotypes_fp, imputed_sample_ids_list_fp, known_genotypes_fp, known_sample_ids_list_fp);
	}
	else if (strcmp(argv[1], "-run_GIMP_sliding_window_Viterbi_State_Reduction") == 0)
	{
		if (argc != 18)
		{
			fprintf(stderr, "%s %s [Phased Reference genotypes path] [Reference sample ids list path] \
[Testing Tag genotypes path] [Testing sample ids list path] \
[Known haplocoded genotypes path] [Minimum tag-target genetic distance (cM)] \
[Genetic distance maps directory] [Target variant focus regions BED path] \
[Math mode : 0: Log, 1 : Linear] [Global scaler in log domain (Linear mode only)] \
[# variants per block] \
[Max # tag variants in window] \
[Effective pop. size (N_e)] \
[Posterior Mode : 0 (Single path) / 1 (Weighted PW)] \
[Center - 2 - untyped target maximum buffer length] \
[Genotype probability output path]\n", argv[0], argv[1]);
			exit(0);
		}

		char* reference_tag_target_haplocoded_genotypes_fp = argv[2];
		char* ref_sample_ids_list_fp = argv[3];
		char* testing_tag_haplocoded_genotypes_fp = argv[4];
		char* testing_sample_ids_list_fp = argv[5];
		char* known_tag_target_haplocoded_genotypes_fp = argv[6];
		double tag_2_tag_distance_cM = atof(argv[7]);
		char* recombination_rate_dir = argv[8];
		char* target_focus_reg_BED_fp = argv[9];
		int math_mode = atoi(argv[10]);
		double global_scaler_in_log = atof(argv[11]);
		int l_blocks = atoi(argv[12]);
		int max_n_tag_vars_per_window = atoi(argv[13]);
		double N_e = atof(argv[14]);
		int posterior_mode = atoi(argv[15]);
		double l_target_2_center_pred_buffer_in_cM = atof(argv[16]);
		char* geno_probs_op_fp = argv[17];

		run_Imputation_Viterbi_State_Reduction_Sliding_Windows(reference_tag_target_haplocoded_genotypes_fp, ref_sample_ids_list_fp,
			testing_tag_haplocoded_genotypes_fp, testing_sample_ids_list_fp,
			known_tag_target_haplocoded_genotypes_fp,
			recombination_rate_dir,
			tag_2_tag_distance_cM,
			N_e,
			target_focus_reg_BED_fp,
			math_mode,
			global_scaler_in_log,
			l_blocks,
			max_n_tag_vars_per_window,
			posterior_mode,
			l_target_2_center_pred_buffer_in_cM,
			geno_probs_op_fp);
	}
	else if (strcmp(argv[1], "-run_GIMP_sliding_window_ForeBack_State_Reduction") == 0)
	{
		if (argc != 18)
		{
			fprintf(stderr, "%s %s [Phased Reference genotypes path] [Reference sample ids list path] \
[Testing Tag genotypes path] [Testing sample ids list path] \
[Known haplocoded genotypes path] [Minimum tag-target genetic distance (cM)] \
[Genetic distance maps directory] [Target variant focus regions BED path] \
[Math mode : 0 : Log, 1 : Linear] [Global scaler in log domain(Linear mode only)] \
[# variants per block] \
[Max # tag variants in window] \
[Effective pop. size (N_e)] \
[Posterior Mode: 0 (Single path) / 1 (Weighted PW)] \
[Center-2-untyped target maximum buffer length] \
[Genotype probability output path]\n", argv[0], argv[1]);
			exit(0);
		}

		char* reference_tag_target_haplocoded_genotypes_fp = argv[2];
		char* ref_sample_ids_list_fp = argv[3];
		char* testing_tag_haplocoded_genotypes_fp = argv[4];
		char* testing_sample_ids_list_fp = argv[5];
		char* known_tag_target_haplocoded_genotypes_fp = argv[6];
		double tag_2_tag_distance_cM = atof(argv[7]);
		char* recombination_rate_dir = argv[8];
		char* target_focus_reg_BED_fp = argv[9];
		int math_mode = atoi(argv[10]);
		double global_scaler_in_log = atof(argv[11]);
		int l_blocks = atoi(argv[12]);
		int max_n_tag_vars_per_window = atoi(argv[13]);
		double N_e = atof(argv[14]);
		int posterior_mode = atoi(argv[15]);
		double l_target_2_center_pred_buffer_in_cM = atof(argv[16]);
		char* geno_probs_op_fp = argv[17];

		run_Imputation_ForeBack_State_Reduction_Sliding_Windows_Math_Mode(reference_tag_target_haplocoded_genotypes_fp, ref_sample_ids_list_fp,
			testing_tag_haplocoded_genotypes_fp, testing_sample_ids_list_fp,
			known_tag_target_haplocoded_genotypes_fp,
			recombination_rate_dir,
			tag_2_tag_distance_cM,
			N_e,
			target_focus_reg_BED_fp,
			math_mode,
			global_scaler_in_log,
			l_blocks,
			max_n_tag_vars_per_window,
			posterior_mode,
			l_target_2_center_pred_buffer_in_cM,
			geno_probs_op_fp);
	} // -run_GIMP_sliding_window_ForeBack_State_Reduction option.
	else if (strcmp(argv[1], "-extract_genotype_signals_per_VCF") == 0)
	{
		if (argc != 11)
		{
			fprintf(stderr, "USAGE: %s -extract_genotype_signals_per_VCF [VCF file path] [VCF sample ids list file path] \
[Variant regions BED file path] [chromosome id to process] \
[Binary sequence directory (Necessary for ref matching)] [Match reference allele? (0/1)] \
[Match region names? (0/1)] \
[Haplotype specific encoding (0/1)] \
[Output file path]\n", argv[0]);
			exit(0);
		}

		char* vcf_fp = argv[2];
		char* vcf_sample_ids_list_fp = argv[3];
		char* var_regions_BED_fp = argv[4];
		char* chr_id_2_process = argv[5];
		char* bin_seq_dir = argv[6];
		bool match_ref_alleles_flag = (argv[7][0] == '1');
		bool match_region_names_flag = (argv[8][0] == '1');
		bool haplotype_specific_encoding = (argv[9][0] == '1');
		char* op_fp = argv[10];

		extract_genotype_signals_per_VCF(vcf_fp,
			vcf_sample_ids_list_fp,
			var_regions_BED_fp,
			chr_id_2_process,
			bin_seq_dir,
			match_ref_alleles_flag,
			match_region_names_flag,
			haplotype_specific_encoding,
			op_fp);
	} // -extract_genotype_signals_per_1kg_VCF option.
	else if (strcmp(argv[1], "-extract_genotype_signals_per_region_list") == 0)
	{
		if (argc != 6)
		{
			fprintf(stderr, "USAGE: %s -extract_genotype_signals_per_region_list [Genotype signals BED file path] [sample ids list file path] [BED file with regions] [Output file path]\n", argv[0]);
			exit(0);
		}

		char* geno_sig_regs_BED_fp = argv[2];
		char* sample_ids_list_fp = argv[3];
		char* regions_BED_fp = argv[4];
		char* op_fp = argv[5];

		extract_genotype_signals_per_region_list(geno_sig_regs_BED_fp, sample_ids_list_fp, regions_BED_fp, op_fp);
	} // -extract_genotype_signals_per_region_list option.
	else if (strcmp(argv[1], "-extract_genotype_signals_per_subsample_list") == 0)
	{
		if (argc != 6)
		{
			fprintf(stderr, "USAGE: %s -extract_genotype_signals_per_subsample_list [Genotype signals BED file path] [sample ids list file path] [Subsample ids list file path] [Output file path]\n", argv[0]);
			exit(0);
		}

		char* geno_sig_regs_BED_fp = argv[2];
		char* sample_ids_list_fp = argv[3];
		char* subsample_ids_list_fp = argv[4];
		char* op_fp = argv[5];

		extract_genotype_signals_per_subsample_list(geno_sig_regs_BED_fp, sample_ids_list_fp, subsample_ids_list_fp, op_fp);
	}
	else if (strcmp(argv[1], "-binarize_genotype_signals_BED") == 0)
	{
		if (argc != 5)
		{
			fprintf(stderr, "USAGE: %s -binarize_genotype_signals_BED [Genotype signals BED file path] [VCF sample ids list file path (Use EpiLeak to extract)] [Output file path]\n", argv[0]);
			exit(0);
		}

		char* variant_geno_sig_regs_BED_fp = argv[2];
		char* vcf_sample_ids_list_fp = argv[3];
		char* op_fp = argv[4];

		vector<char*>* vcf_sample_ids = buffer_file(vcf_sample_ids_list_fp);
		fprintf(stderr, "Loaded %d sample ids.\n", vcf_sample_ids->size());

		binarize_variant_genotype_signal_regions(NULL, variant_geno_sig_regs_BED_fp, vcf_sample_ids, op_fp);
	} // -binarize_genotype_signals_BED option.
	else if (strcmp(argv[1], "-exclude") == 0)
	{
		// Exclude the regions in second region file from the first region file.
		if (argc != 5)
		{
			printf("USAGE: %s -exclude [region1 file path] [region2 file path] [Strand specific (y/n)]\n", argv[0]);
			exit(0);
		}

		char* reg1_fp = argv[2];
		char* reg2_fp = argv[3];
		bool strand_specific = t_string::compare_strings_ci(argv[4], "y") ? (true) : (false);

		// Load the regions: Depends on the file format.
		vector<t_annot_region*>* region_list1 = load_BED_with_line_information(reg1_fp);
		vector<t_annot_region*>* region_list2 = load_BED_with_line_information(reg2_fp);

		vector<t_annot_region*>* excluded_region_list = exclude_annot_regions(region_list1,
			region_list2,
			strand_specific);

		fprintf(stderr, "Found %d excluded regions.\n", excluded_region_list->size());

		printf("Dumping excluded regions.\n");
		FILE* f_exc = fopen("excluded.bed", "w");
		for (int i = 0; i < excluded_region_list->size(); i++)
		{
			t_annot_region* src_region = (t_annot_region*)(excluded_region_list->at(i)->data);
			fprintf(f_exc, "%s\n", (char*)(src_region->data));
		}
		fclose(f_exc);
	}
	else if (strcmp(argv[1], "-convert_haplocoded_2_genocoded") == 0)
	{
		if (argc != 5)
		{
			fprintf(stderr, "USAGE: %s -convert_haplocoded_2_genocoded [Haplocoded genotype matbed file path] \
[Haplocoded genotype matrix sample ids list path] \
[Output genotype signal regions bed file path]\n", argv[0]);
			exit(0);
		}

		char* haplocoded_geno_sig_regs_fp = argv[2];
		char* sample_ids_list_fp = argv[3];
		char* op_matbed_fp = argv[4];

		vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);
		vector<t_annot_region*>* haplocoded_geno_regs = load_variant_signal_regions_wrapper(haplocoded_geno_sig_regs_fp, sample_ids_list_fp);

		fprintf(stderr, "Loaded %d haplocoded genotype regions for %d samples.\n", haplocoded_geno_regs->size(), sample_ids->size());

		for (int i_reg = 0; i_reg < haplocoded_geno_regs->size(); i_reg++)
		{
			void** cur_reg_info = (void**)(haplocoded_geno_regs->at(i_reg)->data);
			char* haplocoded_geno_sigs = (char*)(cur_reg_info[0]);

			for (int i_s = 0; i_s < sample_ids->size(); i_s++)
			{
				int genocoded_geno = get_genotype_per_haplocoded_genotype(haplocoded_geno_sigs[i_s]);
				haplocoded_geno_sigs[i_s] = genocoded_geno;
			} // i_s loop.
		} // i_reg loop.

		  // Save.
		fprintf(stderr, "Saving to %s.\n", op_matbed_fp);
		binarize_variant_genotype_signal_regions(haplocoded_geno_regs, NULL, sample_ids, op_matbed_fp);
	} // -convert_haplocoded_2_genocoded option.
	else if (strcmp(argv[1], "-dump_plain_geno_signal_regions") == 0)
	{
		if (argc != 6)
		{
			fprintf(stderr, "USAGE: %s -dump_plain_geno_signal_regions [Genotype signals matBED file path] [sample ids list file path] [Save regions only? (0/1)] [Output file path]\n", argv[0]);
			exit(0);
		}

		char* geno_sig_regs_BED_fp = argv[2];
		char* sample_ids_list_fp = argv[3];
		bool regs_only_flag = (argv[4][0] == '1');
		char* op_fp = argv[5];

		vector<t_annot_region*>* geno_sig_regs = load_variant_signal_regions_wrapper(geno_sig_regs_BED_fp, sample_ids_list_fp);

		vector<char*>* sample_ids = buffer_file(sample_ids_list_fp);

		dump_geno_sig_regs_plain(geno_sig_regs, sample_ids, regs_only_flag, op_fp);
	} // -dump_plain_geno_signal_regions option.
	else if (strcmp(argv[1], "-convert_genotype_signal_regions_2_VCF") == 0)
	{
		if (argc != 6)
		{
			fprintf(stderr, "USAGE: %s -convert_genotype_signal_regions_2_VCF [Genotype signal regions file path] [Sample ids list file path] [Phased output? (0/1)] [Output VCF file path]\n", argv[0]);
			exit(0);
		}

		char* geno_sig_regs_fp = argv[2];
		char* sample_ids_list_fp = argv[3];
		bool phase_op = (argv[4][0] == '1');
		char* op_vcf_fp = argv[5];

		convert_genotype_signal_regs_2_VCF(geno_sig_regs_fp, sample_ids_list_fp, phase_op, op_vcf_fp);
	} // -convert_genotype_signal_regions_2_VCF option.
	else if (strcmp(argv[1], "-separate_LoHaMMer_jobs") == 0)
	{
		if (argc != 5)
		{
			fprintf(stderr, "%s -separate_LoHaMMer_jobs [cmd list (stdin for command line input)] [# jobs] [job dir. prefix]\n", argv[0]);
			exit(0);
		}

		vector<char*>* cmd_list = buffer_file(argv[2]);
		if (cmd_list == NULL)
		{
			fprintf(stderr, "Could not load the command lines from %s\n", argv[2]);
			exit(0);
		}

		int n_jobs = atoi(argv[3]);
		char* job_dir_prefix = argv[4];		

		char cur_cmd_dir[1000];
		int i_cmd = 0;

		// Process all the commands.
		int n_created_dirs = 0;
		while (i_cmd < cmd_list->size())
		{
			// Go over all the directories and save one line per directory job.
			for (int i_dir = 0;
				i_cmd < cmd_list->size() &&
				i_dir < n_jobs;
				i_dir++)
			{
				sprintf(cur_cmd_dir, "%s_%d", job_dir_prefix, i_dir);

				// Create directory.
				char mkdir_cmd[1000];
				sprintf(mkdir_cmd, "mkdir %s", cur_cmd_dir);
				char cur_dir_job_fp[1000];
				sprintf(cur_dir_job_fp, "%s/job.csh", cur_cmd_dir);

				if (!check_file(cur_dir_job_fp))
				{
					system(mkdir_cmd);
					n_created_dirs++;
					fprintf(stderr, "Created %s\n", cur_cmd_dir);
				}

				// Open the job file.
				FILE* f_cur_dir_job = open_f(cur_dir_job_fp, "a");
				fprintf(f_cur_dir_job, "%s\n", cmd_list->at(i_cmd));
				i_cmd++;
				fclose(f_cur_dir_job);
			} // i_dir loop.
		} // cmd's loop.

		// Write the submission script.
		char submission_script_fp[1000];
		sprintf(submission_script_fp, "%s_submission_script.csh", job_dir_prefix);
		FILE* f_sub_script = open_f(submission_script_fp, "w");
		for (int i_dir = 0; i_dir < n_created_dirs; i_dir++)
		{
			sprintf(cur_cmd_dir, "%s_%d", job_dir_prefix, i_dir);

			if (i_dir > 0)
			{
				fprintf(f_sub_script, "cd ../%s\nchmod 755 job.csh;nohup ./job.csh > op.txt &\n", cur_cmd_dir);
			}
			else
			{
				fprintf(f_sub_script, "cd %s\nchmod 755 job.csh;nohup ./job.csh > op.txt &\n", cur_cmd_dir);
			}
		} // i_dir loop.
		fclose(f_sub_script);

		// Change the permission to executable.
		char chmod_cmd[1000];
		sprintf(chmod_cmd, "chmod 755 %s", submission_script_fp);

		if (system(chmod_cmd) != 0)
		{
			fprintf(stderr, "chmod failed.\n");
		}
	} // -separate_write_PBS_job_scripts_per_cmd_list option.

	FILE* f_beacon = open_f("beacon.log", "a");
	clock_t end_c = clock();
	fprintf(f_beacon, "%s finished (%s) in %d seconds.\n", argv[0], argv[1], (int)((end_c - start_c) / CLOCKS_PER_SEC));
	fprintf(stderr, "%s finished (%s) in %d seconds.\n", argv[0], argv[1], (int)((end_c - start_c) / CLOCKS_PER_SEC));
	fclose(f_beacon);
}
