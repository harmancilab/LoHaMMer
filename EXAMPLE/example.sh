#!/bin/bash

KG_dir=1kG_DATA

# Download the Illumina Duo V3 manifest.
PARSE_ILLUMINADUOV3_MANIFEST=1
if [ ${PARSE_ILLUMINADUOV3_MANIFEST} == 1 ]
then
        wget -c ftp://webdata2:webdata2@ussd-ftp.illumina.com/MyIllumina/c1859532-43c5-4df2-a1b4-c9f0374476d2/Human1M-Duov3_H.csv

        head -n 8 Human1M-Duov3_H.csv | tail -n 1 | awk {'n_ents=split($0, arr, ",");for(i=1;i<=n_ents;i++){print i":"arr[i]}'} > col_ids.list

        echo "Parsing SNV loci from manifest"
        grep -v cnv Human1M-Duov3_H.csv | awk 'BEGIN{FS=","}{if(NR>18 && $1!="CHR" && $10!=0 && $11>1)print $10"\t"$11-1"\t"$11"\t"$1" "$2" "$17"\t.\t"$22}' > Human1M-Duov3_H.csv_SNV_loci.bed

        exit
fi

KG_PARSE=1
if [ ${KG_PARSE} == 1 ]
then
	mkdir $KG_dir

        echo "Downloading data"
        wget -c ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
        gzip -cd ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | awk 'BEGIN{FS="\t"}{if($1=="#CHROM"){for(i=10;i<=NF;i++){print $i};exit}}' > 1kg_sample_ids.list

        # Extract the common SNVs.
        echo "Parsing min MAF regions"
        min_MAF=0.01
        gzip -cd ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | awk -v min_MAF=${min_MAF} {'curAF=0;split($8, arr, ";");for(i=1; i <= length(arr); i++){split(arr[i], arr2, "=");if(arr2[1]=="AF"){curAF=arr2[2];break;}};curMAF=curAF;if(curAF>0.5){curMAF=1-curAF;};if(curMAF>min_MAF && $1==20 && length($4)==1 && length($5)==1)print $1"\t"$2-1"\t"$2"\t"$3"_"$4"_"$5"_"curAF"\t.\t+"'} > snv_region_${min_MAF}.bed

        # Extract the matrix.
        echo "Extracting MAF > ${min_MAF} variants on snv_region_${min_MAF}.bed"
        gzip -cd ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | awk {'if($1==20 && length($4)==1 && length($5)==1)print $0'} | ./bin/LoHaMMer -extract_genotype_signals_per_VCF stdin 1kg_sample_ids.list snv_region_${min_MAF}.bed hg19.list nogenome 0 1 $KG_DIR/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_haplocoded.matbed.gz

	exit
fi


PARSE_TAG_TARGET_SNPs=0
if [ ${PARSE_TAG_TARGET_SNPs} == 1 ]
then
	# Extract the common SNVs.
	min_MAF=0.01

	# Extract the tag and target SNP matrices.
	echo "Extracting the genotypes for the tag SNVs"
	../bin/LoHaMMer -extract_genotype_signals_per_region_list ${KG_dir}/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_haplocoded.matbed.gz ${KG_dir}/1kg_sample_ids.list Human1M-Duov3_H.csv_SNV_loci.bed tag_snvs_haplocoded.matbed.gz	

	echo "Saving target variants"
	../bin/LoHaMMer -dump_plain_geno_signal_regions ${KG_dir}/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_haplocoded.matbed.gz ${KG_dir}/1kg_sample_ids.list 0 temp.txt
	../bin/LoHaMMer -exclude temp.txt Human1M-Duov3_H.csv_SNV_loci.bed n
	rm -f temp.txt
	mv excluded.bed target_snvs.matbed.gz.txt
	awk 'BEGIN{FS="\t"}{n_toks=split($4, arr, "_");curAF=arr[n_toks];curMAF=curAF;if(curMAF>.5){curMAF=1-curMAF;}if(curMAF>0.01)print $0}' target_snvs.matbed.gz.txt > temp.txt
	mv temp.txt target_snvs.matbed.gz.txt

	echo "Binarizing target variants."
	../bin/LoHaMMer -binarize_genotype_signals_BED target_snvs.matbed.gz.txt ${KG_dir}/1kg_sample_ids.list target_snvs_haplocoded.matbed.gz

	exit
fi

N_TRAINING=1000
N_TESTING=200

BUILD_TAG_TARGET_TRAINING_TESTING_DATA=0
if [ ${BUILD_TAG_TARGET_TRAINING_TESTING_DATA} == 1 ]
then
	echo "Separating data"

	GET_NEW_TRAINING_TESTING=1
	if [ ${GET_NEW_TRAINING_TESTING} == 1 ]
	then
		shuf $KG_dir/1kg_sample_ids.list | head -n $N_TRAINING > training_sample_ids.list
		grep -v -f training_sample_ids.list $KG_dir/1kg_sample_ids.list | head -n $N_TESTING > testing_sample_ids.list
	fi

	../bin/LoHaMMer -extract_genotype_signals_per_subsample_list tag_snvs_haplocoded.matbed.gz ${KG_dir}/1kg_sample_ids.list training_sample_ids.list tag_snvs_haplocoded.matbed.gz_training.matbed.gz
	../bin/LoHaMMer -convert_haplocoded_2_genocoded tag_snvs_haplocoded.matbed.gz_training.matbed.gz training_sample_ids.list tag_snvs_genocoded.matbed.gz_training.matbed.gz	
	../bin/LoHaMMer -dump_plain_geno_signal_regions tag_snvs_genocoded.matbed.gz_training.matbed.gz training_sample_ids.list 0 tag_snvs_genocoded.matbed.gz_training.matbed.gz.txt
	sort -n -k2,2 tag_snvs_genocoded.matbed.gz_training.matbed.gz.txt > sorted_tag_snvs_genocoded.matbed.gz_training.matbed.gz.txt

        ../bin/LoHaMMer -extract_genotype_signals_per_subsample_list tag_snvs_haplocoded.matbed.gz ${KG_dir}/1kg_sample_ids.list testing_sample_ids.list tag_snvs_haplocoded.matbed.gz_testing.matbed.gz
        ../bin/LoHaMMer -convert_haplocoded_2_genocoded tag_snvs_haplocoded.matbed.gz_testing.matbed.gz testing_sample_ids.list tag_snvs_genocoded.matbed.gz_testing.matbed.gz
	../bin/LoHaMMer -dump_plain_geno_signal_regions tag_snvs_genocoded.matbed.gz_testing.matbed.gz testing_sample_ids.list 0 tag_snvs_genocoded.matbed.gz_testing.matbed.gz.txt
	sort -n -k2,2 tag_snvs_genocoded.matbed.gz_testing.matbed.gz.txt > sorted_tag_snvs_genocoded.matbed.gz_testing.matbed.gz.txt

        ../bin/LoHaMMer -extract_genotype_signals_per_subsample_list target_snvs_haplocoded.matbed.gz ${KG_dir}/1kg_sample_ids.list training_sample_ids.list target_snvs_haplocoded.matbed.gz_training.matbed.gz
        ../bin/LoHaMMer -convert_haplocoded_2_genocoded target_snvs_haplocoded.matbed.gz_training.matbed.gz training_sample_ids.list target_snvs_genocoded.matbed.gz_training.matbed.gz
	../bin/LoHaMMer -dump_plain_geno_signal_regions target_snvs_genocoded.matbed.gz_training.matbed.gz training_sample_ids.list 0 target_snvs_genocoded.matbed.gz_training.matbed.gz.txt
	sort -n -k2,2 target_snvs_genocoded.matbed.gz_training.matbed.gz.txt > sorted_target_snvs_genocoded.matbed.gz_training.matbed.gz.txt

        ../bin/LoHaMMer -extract_genotype_signals_per_subsample_list target_snvs_haplocoded.matbed.gz ${KG_dir}/1kg_sample_ids.list testing_sample_ids.list target_snvs_haplocoded.matbed.gz_testing.matbed.gz
        ../bin/LoHaMMer -convert_haplocoded_2_genocoded target_snvs_haplocoded.matbed.gz_testing.matbed.gz testing_sample_ids.list target_snvs_genocoded.matbed.gz_testing.matbed.gz
	../bin/LoHaMMer -dump_plain_geno_signal_regions target_snvs_genocoded.matbed.gz_testing.matbed.gz testing_sample_ids.list 0 target_snvs_genocoded.matbed.gz_testing.matbed.gz.txt
	sort -n -k2,2 target_snvs_genocoded.matbed.gz_testing.matbed.gz.txt > sorted_target_snvs_genocoded.matbed.gz_testing.matbed.gz.txt

	../bin/LoHaMMer -dump_plain_geno_signal_regions target_snvs_haplocoded.matbed.gz_training.matbed.gz training_sample_ids.list 0 target_snvs_haplocoded.matbed.gz_training.matbed.gz.txt
	../bin/LoHaMMer -dump_plain_geno_signal_regions tag_snvs_haplocoded.matbed.gz_training.matbed.gz training_sample_ids.list 0 tag_snvs_haplocoded.matbed.gz_training.matbed.gz.txt
	cat tag_snvs_haplocoded.matbed.gz_training.matbed.gz.txt target_snvs_haplocoded.matbed.gz_training.matbed.gz.txt > tag_target_snvs_haplocoded.matbed.gz_training.matbed.gz.txt

	exit
fi

# Run EAGLE; Update the tag data.
RUN_EAGLE=0
if [ $RUN_EAGLE == 1 ]
then
	./run_eagle.sh

	gzip -cd tag_snvs_haplocoded.matbed.gz_testing.matbed.gz.txt_name_fix.txt.vcf_header.vcf.gz_eagle_phase.vcf.gz | awk 'BEGIN{FS="\t";OFS="\t"}{if($1==20){$3=$3"_"$4"_"$5;print $1"\t"$2-1"\t"$2"\t"$3"\t.\t+"}}' > targets.bed

	gzip -cd tag_snvs_haplocoded.matbed.gz_testing.matbed.gz.txt_name_fix.txt.vcf_header.vcf.gz_eagle_phase.vcf.gz | awk 'BEGIN{FS="\t";OFS="\t"}{if($1==20 && length($4)==1 && length($5)==1){$3=$3"_"$4"_"$5;print $0}}' | ../bin/LoHaMMer -extract_genotype_signals_per_VCF stdin testing_sample_ids.list targets.bed hg19.list nogenome 0 1 1 tag_snvs_haplocoded.matbed.gz_testing.matbed.gz_EAGLE_phased.matbed.gz

	../bin/LoHaMMer -dump_plain_geno_signal_regions tag_snvs_haplocoded.matbed.gz_testing.matbed.gz_EAGLE_phased.matbed.gz testing_sample_ids.list 0 tag_snvs_haplocoded.matbed.gz_testing.matbed.gz.txt

	exit
fi

RUN_HMM=1
if [ $RUN_HMM == 1 ]
then
	#snp_pos=`grep "rs139969097_G_C_0.0101837" "tag_target_snvs_haplocoded.matbed.gz_training.matbed.gz.txt" | cut -f2`

	#l_vic=5000000
	#awk -v snp_pos=$snp_pos -v l_vic=$l_vic {'if($2>snp_pos-l_vic && $2<snp_pos+l_vic){print $0}'} tag_target_snvs_haplocoded.matbed.gz_training.matbed.gz.txt > tag_target_ref_geno.txt
        #awk -v snp_pos=$snp_pos -v l_vic=$l_vic {'if($2>snp_pos-l_vic && $2<snp_pos+l_vic){print $0}'} tag_snvs_haplocoded.matbed.gz_testing.matbed.gz.txt > tag_testing_ref_geno.txt

	focus_start=25000000
	focus_end=35000000
	l_win=250000
	starts=`seq $focus_start $l_win $focus_end`

        RUN_FOREBACK_STATE_REDUCTION=1
        if [ $RUN_FOREBACK_STATE_REDUCTION == 1 ]
        then
		rm -f temp_imp.sh
		for cur_start in ${starts[@]}
		do
			cur_end=`awk -v l_win=$l_win -v cur_start=$cur_start 'BEGIN{print cur_start+l_win}'`
			#print win_cMs[1]"\t"N_e[i]"\t"l_center_2_target_buffers[1]"\t"l_block"\t"lin_scaler"\t"math_mode"\t"posterior_mode
			awk -f enumerate_parameters.awk | awk -v workdir=$PWD -v cur_start=$cur_start -v cur_end=$cur_end 'BEGIN{FS="\t"}{cur_win_cM=$1;cur_N_e=$2;cur_l_center_2_target_buffer=$3;l_block=$4;lin_scaler=$5;math_mode=$6;posterior_mode=$7;max_n_vars=$8;print "echo -e \"20\\t"cur_start"\\t"cur_end"\" > focus_reg.bed;/usr/bin/time -o timing_memory.log --append -f %e\"\\t\"%M SHiMMer -run_GIMP_sliding_window_ForeBack_State_Reduction "workdir"/tag_target_snvs_haplocoded.matbed.gz_training.matbed.gz.txt "workdir"/training_sample_ids.list "workdir"/tag_snvs_haplocoded.matbed.gz_testing.matbed.gz.txt "workdir"/testing_sample_ids.list "workdir"/target_snvs_haplocoded.matbed.gz_testing.matbed.gz "cur_win_cM" "workdir"/genetic_maps focus_reg.bed "math_mode" "lin_scaler" "l_block" "max_n_vars" "cur_N_e" "posterior_mode" "cur_l_center_2_target_buffer" gimp_"cur_win_cM"_"cur_N_e"_"cur_l_center_2_target_buffer"_"max_n_vars"_"cur_start".probs"}' >> temp_imp.sh
		done

		cat temp_imp.sh | sort -u > temp
		mv temp temp_imp.sh

		rm -f -r q_foreback_reduc*
		job_tools -separate_write_PBS_job_scripts_per_cmd_list temp_imp.sh q_foreback_reduc 60 -1 -1 -1 -1

		exit
        fi

	RUN_VITERBI_STATE_REDUCTION=1
        if [ $RUN_VITERBI_STATE_REDUCTION == 1 ]
        then
                rm -f temp_imp.sh
		for cur_start in ${starts[@]}
		do
			cur_end=`awk -v l_win=$l_win -v cur_start=$cur_start 'BEGIN{print cur_start+l_win}'`
			#print win_cMs[1]"\t"N_e[i]"\t"l_center_2_target_buffers[1]"\t"l_block"\t"lin_scaler"\t"math_mode"\t"posterior_mode
			awk -f enumerate_parameters.awk | awk -v workdir=$PWD -v cur_start=$cur_start -v cur_end=$cur_end 'BEGIN{FS="\t"}{cur_win_cM=$1;cur_N_e=$2;cur_l_center_2_target_buffer=$3;l_block=$4;lin_scaler=$5;math_mode=$6;posterior_mode=$7;max_n_vars=$8;print "echo -e \"20\\t"cur_start"\\t"cur_end"\" > focus_reg.bed;/usr/bin/time -o timing_memory.log --append -f %e\"\\t\"%M SHiMMer -run_GIMP_sliding_window_Viterbi_State_Reduction "workdir"/tag_target_snvs_haplocoded.matbed.gz_training.matbed.gz.txt "workdir"/training_sample_ids.list "workdir"/tag_snvs_haplocoded.matbed.gz_testing.matbed.gz.txt "workdir"/testing_sample_ids.list "workdir"/target_snvs_haplocoded.matbed.gz_testing.matbed.gz "cur_win_cM" "workdir"/genetic_maps focus_reg.bed "math_mode" "lin_scaler" "l_block" "max_n_vars" "cur_N_e" "posterior_mode" "cur_l_center_2_target_buffer" gimp_"cur_win_cM"_"cur_N_e"_"cur_l_center_2_target_buffer"_"max_n_vars"_"cur_start".probs"}' >> temp_imp.sh
		done

		cat temp_imp.sh | sort -u > temp
		mv temp temp_imp.sh

                rm -f -r q_viterbimp_reduc*
                job_tools -separate_write_PBS_job_scripts_per_cmd_list temp_imp.sh q_viterbimp_reduc 60 -1 -1 -1 -1

		exit
        fi

fi




