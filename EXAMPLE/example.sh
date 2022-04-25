#!/bin/bash

if [[ $# -lt 1 ]]
then
	echo "USAGE: $0 [option]
Options:
	-get_typed_vars
	-get_parse_1kG
	-parse_tag_target
	-build_query_reference_panels
	-eagle_phase_typed_variants
	-run_vicinity_HMM
	-get_R2_stats"

	exit
fi

cmd_option=$1

KG_dir=1kG_DATA

sed -i $'s/\r$//' run_eagle.sh
chmod 755 run_eagle.sh

# Check eagle.
if [[ ! -d "EAGLE/Eagle_v2.4.1" ]]
then
	sed -i $'s/\r$//' EAGLE/get_eagle.sh
	chmod 755 EAGLE/get_eagle.sh
	echo "You need to download eagle2.4.1 under EAGLE/"
	exit
fi

# Extract the common SNVs.
min_MAF=0.05

# Sizes of query and reference panels.
N_TRAINING=200
N_TESTING=10

# Download the Illumina Duo V3 manifest.
PARSE_ILLUMINADUOV3_MANIFEST=1
#if [ ${PARSE_ILLUMINADUOV3_MANIFEST} == 1 ]
if [[ ${cmd_option} == "-get_typed_vars" ]]
then
        wget -c ftp://webdata2:webdata2@ussd-ftp.illumina.com/MyIllumina/c1859532-43c5-4df2-a1b4-c9f0374476d2/Human1M-Duov3_H.csv

        head -n 8 Human1M-Duov3_H.csv | tail -n 1 | awk {'n_ents=split($0, arr, ",");for(i=1;i<=n_ents;i++){print i":"arr[i]}'} > col_ids.list

        echo "Parsing SNV loci from manifest"
        grep -v cnv Human1M-Duov3_H.csv | awk 'BEGIN{FS=","}{if(NR>18 && $1!="CHR" && $10!=0 && $11>1)print $10"\t"$11-1"\t"$11"\t"$1" "$2" "$17"\t.\t"$22}' > Human1M-Duov3_H.csv_SNV_loci.bed

        exit
fi

#KG_PARSE=1
if [[ ${cmd_option} == "-get_parse_1kG" ]]
then
	mkdir $KG_dir

	wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes
	chmod 755 fetchChromSizes
	./fetchChromSizes hg19 > hg19.list

    echo "Downloading data"
    wget -c ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
    gzip -cd ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | awk 'BEGIN{FS="\t"}{if($1=="#CHROM"){for(i=10;i<=NF;i++){print $i};exit}}' > ${KG_dir}/1kg_sample_ids.list

    # Extract the common SNVs.
    echo "Parsing min MAF regions"
    gzip -cd ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | awk -v min_MAF=${min_MAF} {'curAF=0;split($8, arr, ";");for(i=1; i <= length(arr); i++){split(arr[i], arr2, "=");if(arr2[1]=="AF"){curAF=arr2[2];break;}};curMAF=curAF;if(curAF>0.5){curMAF=1-curAF;};if(curMAF>min_MAF && $1==20 && length($4)==1 && length($5)==1)print $1"\t"$2-1"\t"$2"\t"$3"_"$4"_"$5"_"curAF"\t.\t+"'} > snv_region_${min_MAF}.bed

    # Extract the matrix.
    echo "Extracting MAF > ${min_MAF} variants on snv_region_${min_MAF}.bed"
    gzip -cd ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | awk {'if($1==20 && length($4)==1 && length($5)==1)print $0'} | ../bin/LoHaMMer -extract_genotype_signals_per_VCF stdin ${KG_dir}/1kg_sample_ids.list snv_region_${min_MAF}.bed hg19.list nogenome 0 1 1 $KG_dir/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz_haplocoded.matbed.gz

	exit
fi

if [[ ! -d "${KG_dir}" ]]
then
	echo "Could not find 1kG directory @ ${KG_dir}/"
	exit
fi

#PARSE_TAG_TARGET_SNPs=0
if [[ ${cmd_option} == "-parse_tag_target" ]]
then
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

if [[ ${cmd_option} == "-build_query_reference_panels" ]]
then
	echo "Building query/reference panels."

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
if [[ ${cmd_option} == "-eagle_phase_typed_variants" ]]
then
	./run_eagle.sh

	gzip -cd tag_snvs_haplocoded.matbed.gz_testing.matbed.gz.txt_name_fix.txt.vcf_header.vcf.gz_eagle_phase.vcf.gz | awk 'BEGIN{FS="\t";OFS="\t"}{if($1==20){$3=$3"_"$4"_"$5;print $1"\t"$2-1"\t"$2"\t"$3"\t.\t+"}}' > targets.bed

	gzip -cd tag_snvs_haplocoded.matbed.gz_testing.matbed.gz.txt_name_fix.txt.vcf_header.vcf.gz_eagle_phase.vcf.gz | awk 'BEGIN{FS="\t";OFS="\t"}{if($1==20 && length($4)==1 && length($5)==1){$3=$3"_"$4"_"$5;print $0}}' | ../bin/LoHaMMer -extract_genotype_signals_per_VCF stdin testing_sample_ids.list targets.bed hg19.list nogenome 0 1 1 tag_snvs_haplocoded.matbed.gz_testing.matbed.gz_EAGLE_phased.matbed.gz

	../bin/LoHaMMer -dump_plain_geno_signal_regions tag_snvs_haplocoded.matbed.gz_testing.matbed.gz_EAGLE_phased.matbed.gz testing_sample_ids.list 0 tag_snvs_haplocoded.matbed.gz_testing.matbed.gz.txt

	exit
fi

if [[ ${cmd_option} == "-run_vicinity_HMM" ]]
then
	focus_start=25000000
	focus_end=35000000
	l_win=250000
	starts=`seq $focus_start $l_win $focus_end`

	cur_win_cM=0.3
	cur_N_e=10000
	cur_l_center_2_target_buffer=0.05
	l_block=10
	lin_scaler=0.2
	math_mode=1
	posterior_mode=0
	max_n_vars=500

	n_parallel_jobs=60

        RUN_FOREBACK_STATE_REDUCTION=1
        if [ $RUN_FOREBACK_STATE_REDUCTION == 1 ]
        then
                rm -f temp_imp_cmds.sh
                for cur_start in ${starts[@]}
                do
                        echo "Processing ${cur_start}"
                        cur_end=`awk -v l_win=$l_win -v cur_start=$cur_start 'BEGIN{print cur_start+l_win}'`
                        echo "echo -e \"20\\t${focus_start}\\t${focus_end}\" > focus_reg.bed;${PWD}/../bin/LoHaMMer -run_GIMP_sliding_window_ForeBack_State_Reduction ${PWD}/tag_target_snvs_haplocoded.matbed.gz_training.matbed.gz.txt ${PWD}/training_sample_ids.list ${PWD}/tag_snvs_haplocoded.matbed.gz_testing.matbed.gz.txt ${PWD}/testing_sample_ids.list ${PWD}/target_snvs_haplocoded.matbed.gz_testing.matbed.gz ${cur_win_cM} ${PWD}/genetic_maps ./focus_reg.bed ${math_mode} ${lin_scaler} ${l_block} ${max_n_vars} ${cur_N_e} ${posterior_mode} ${cur_l_center_2_target_buffer} gimp_${cur_win_cM}_${cur_N_e}_${cur_l_center_2_target_buffer}_${max_n_vars}_cur_start.probs" >> temp_imp_cmds.sh
                done

                rm -f -r q_foreback_*
                ../bin/LoHaMMer -separate_LoHaMMer_jobs temp_imp_cmds.sh ${n_parallel_jobs} q_foreback

                exit
    fi

	RUN_VITERBI_STATE_REDUCTION=1
	if [ $RUN_VITERBI_STATE_REDUCTION == 1 ]
	then
		rm -f temp_imp.sh
		for cur_start in ${starts[@]}
		do
			cur_end=`awk -v l_win=$l_win -v cur_start=$cur_start 'BEGIN{print cur_start+l_win}'`
			echo "echo -e \"20\\t${focus_start}\\t${focus_end}\" > focus_reg.bed;${PWD}/../bin/LoHaMMer -run_GIMP_sliding_window_Viterbi_State_Reduction ${PWD}/tag_target_snvs_haplocoded.matbed.gz_training.matbed.gz.txt ${PWD}/training_sample_ids.list ${PWD}/tag_snvs_haplocoded.matbed.gz_testing.matbed.gz.txt ${PWD}/testing_sample_ids.list ${PWD}/target_snvs_haplocoded.matbed.gz_testing.matbed.gz ${cur_win_cM} ${PWD}/genetic_maps ./focus_reg.bed ${math_mode} ${lin_scaler} ${l_block} ${max_n_vars} ${cur_N_e} ${posterior_mode} ${cur_l_center_2_target_buffer} gimp_${cur_win_cM}_${cur_N_e}_${cur_l_center_2_target_buffer}_${max_n_vars}_cur_start.probs" >> temp_imp_cmds.sh
		done

		cat temp_imp.sh | sort -u > temp
		mv temp temp_imp.sh

        rm -f -r q_viterbimp_*
        ../bin/LoHaMMer -separate_write_PBS_job_scripts_per_cmd_list temp_imp.sh ${n_parallel_jobs} q_viterbimp

		exit
	fi

fi

if [[ ${cmd_option} == "-get_R2_stats" ]]
then
	echo "Pooling genotype probabilities.."
	find . -maxdepth 2 -name 'gimp_*.probs' | grep q_foreback | xargs -Ifiles cat files | sort -u > gimp.probs

	echo "Computing R2 statistics.."
	rm -f R2_stats.txt
	../bin/LoHaMMer -get_R2_per_GIMP_4entry_allelic_probs gimp.probs testing_sample_ids.list sorted_target_snvs_genocoded.matbed.gz_testing.matbed.gz.txt testing_sample_ids.list
fi


