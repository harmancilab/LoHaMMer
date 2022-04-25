RUN_EAGLE=1
if [ $RUN_EAGLE == 1 ]
then
	EAGLE_BIN=${PWD}/EAGLE/Eagle_v2.4.1/eagle
	EAGLE_MAPS=${PWD}/EAGLE//Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz
	#BCFTOOLS_BIN=/home/aharmanci1/bcftools_installation/installation/bin/bcftools
	BCFTOOLS_BIN=bcftools
	MAIN_DATA_DIR=$PWD

		echo "Concatenating tag and target SNVs in the training data."
		../bin/LoHaMMer -dump_plain_geno_signal_regions ${MAIN_DATA_DIR}/tag_snvs_haplocoded.matbed.gz_training.matbed.gz ${MAIN_DATA_DIR}/training_sample_ids.list 0 tag_snvs_haplocoded.matbed.gz_training.matbed.gz.txt
		awk 'BEGIN{FS="\t";OFS="\t"}{split($4, arr_sub, "_");prefix_name=arr_sub[1];split($4, arrl, "[");split(arrl[2], arrr, "]");split(arrr[1], alleles, "/");new_name=prefix_name"_"alleles[1]"_"alleles[2];$4=new_name;print $0}' tag_snvs_haplocoded.matbed.gz_training.matbed.gz.txt | sort -u -k2,2 > tag_snvs_haplocoded.matbed.gz_training.matbed.gz.txt_name_fix.txt

		../bin/LoHaMMer -dump_plain_geno_signal_regions ${MAIN_DATA_DIR}/target_snvs_haplocoded.matbed.gz_training.matbed.gz ${MAIN_DATA_DIR}/training_sample_ids.list 0 target_snvs_haplocoded.matbed.gz_training.matbed.gz.txt
		cat tag_snvs_haplocoded.matbed.gz_training.matbed.gz.txt_name_fix.txt target_snvs_haplocoded.matbed.gz_training.matbed.gz.txt | sort -u -k2,2 | sort -n -k2,2 > ${MAIN_DATA_DIR}/sorted_tag_target_snvs_haplocoded.matbed.gz_training.matbed.gz.txt

		../bin/LoHaMMer -convert_genotype_signal_regions_2_VCF ${MAIN_DATA_DIR}/sorted_tag_target_snvs_haplocoded.matbed.gz_training.matbed.gz.txt ${MAIN_DATA_DIR}/training_sample_ids.list 1 ${MAIN_DATA_DIR}/sorted_tag_target_snvs_haplocoded.matbed.gz_training.matbed.gz.txt.vcf

		# Add the header.
		echo "##fileformat=VCFv4.1
##fileDate=20150121
##source=TCGA
##contig=<ID=1,length=249250621,assembly=b37>
##contig=<ID=2,length=243199373,assembly=b37>
##contig=<ID=3,length=198022430,assembly=b37>
##contig=<ID=4,length=191154276,assembly=b37>
##contig=<ID=5,length=180915260,assembly=b37>
##contig=<ID=6,length=171115067,assembly=b37>
##contig=<ID=7,length=159138663,assembly=b37>
##contig=<ID=8,length=146364022,assembly=b37>
##contig=<ID=9,length=141213431,assembly=b37>
##contig=<ID=10,length=135534747,assembly=b37>
##contig=<ID=11,length=135006516,assembly=b37>
##contig=<ID=12,length=133851895,assembly=b37>
##contig=<ID=13,length=115169878,assembly=b37>
##contig=<ID=14,length=107349540,assembly=b37>
##contig=<ID=15,length=102531392,assembly=b37>
##contig=<ID=16,length=90354753,assembly=b37>
##contig=<ID=17,length=81195210,assembly=b37>
##contig=<ID=18,length=78077248,assembly=b37>
##contig=<ID=19,length=59128983,assembly=b37>
##contig=<ID=20,length=63025520,assembly=b37>
##contig=<ID=21,length=48129895,assembly=b37>
##contig=<ID=22,length=51304566,assembly=b37>
##contig=<ID=X,length=155270560,assembly=b37>
##contig=<ID=Y,length=59373566,assembly=b37>
##contig=<ID=MT,length=16569,assembly=b37>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" > vcf_header.txt

		cat vcf_header.txt ${MAIN_DATA_DIR}/sorted_tag_target_snvs_haplocoded.matbed.gz_training.matbed.gz.txt.vcf > ${MAIN_DATA_DIR}/sorted_tag_target_snvs_haplocoded.matbed.gz_training.matbed.gz.txt.vcf_header.vcf

		# Eagle only takes bcf files.
		bgzip -f ${MAIN_DATA_DIR}/sorted_tag_target_snvs_haplocoded.matbed.gz_training.matbed.gz.txt.vcf_header.vcf
		rm -f ${MAIN_DATA_DIR}/sorted_tag_target_snvs_haplocoded.matbed.gz_training.matbed.gz.txt.vcf_header.vcf.gz.tbi
		tabix -p vcf ${MAIN_DATA_DIR}/sorted_tag_target_snvs_haplocoded.matbed.gz_training.matbed.gz.txt.vcf_header.vcf.gz
		${BCFTOOLS_BIN} view -Ob -o ${MAIN_DATA_DIR}/sorted_tag_target_snvs_haplocoded.matbed.gz_training.matbed.gz.txt.vcf_header.vcf.gz.bcf ${MAIN_DATA_DIR}/sorted_tag_target_snvs_haplocoded.matbed.gz_training.matbed.gz.txt.vcf_header.vcf.gz
		${BCFTOOLS_BIN} index ${MAIN_DATA_DIR}/sorted_tag_target_snvs_haplocoded.matbed.gz_training.matbed.gz.txt.vcf_header.vcf.gz.bcf

        ../bin/LoHaMMer -dump_plain_geno_signal_regions ${MAIN_DATA_DIR}/tag_snvs_haplocoded.matbed.gz_testing.matbed.gz ${MAIN_DATA_DIR}/testing_sample_ids.list 0 tag_snvs_haplocoded.matbed.gz_testing.matbed.gz.txt
        awk 'BEGIN{FS="\t";OFS="\t"}{split($4, arr_sub, "_");prefix_name=arr_sub[1];split($4, arrl, "[");split(arrl[2], arrr, "]");split(arrr[1], alleles, "/");new_name=prefix_name"_"alleles[1]"_"alleles[2];$4=new_name;print $0}' tag_snvs_haplocoded.matbed.gz_testing.matbed.gz.txt | sort -u -k2,2 > tag_snvs_haplocoded.matbed.gz_testing.matbed.gz.txt_name_fix.txt

		../bin/LoHaMMer -convert_genotype_signal_regions_2_VCF tag_snvs_haplocoded.matbed.gz_testing.matbed.gz.txt_name_fix.txt ${MAIN_DATA_DIR}/testing_sample_ids.list 0 ${MAIN_DATA_DIR}/tag_snvs_haplocoded.matbed.gz_testing.matbed.gz.txt_name_fix.txt.vcf
		cat vcf_header.txt ${MAIN_DATA_DIR}/tag_snvs_haplocoded.matbed.gz_testing.matbed.gz.txt_name_fix.txt.vcf > ${MAIN_DATA_DIR}/tag_snvs_haplocoded.matbed.gz_testing.matbed.gz.txt_name_fix.txt.vcf_header.vcf
		bgzip -f ${MAIN_DATA_DIR}/tag_snvs_haplocoded.matbed.gz_testing.matbed.gz.txt_name_fix.txt.vcf_header.vcf
		rm -f ${MAIN_DATA_DIR}/tag_snvs_haplocoded.matbed.gz_testing.matbed.gz.txt_name_fix.txt.vcf_header.vcf.gz.tbi
		tabix -p vcf ${MAIN_DATA_DIR}/tag_snvs_haplocoded.matbed.gz_testing.matbed.gz.txt_name_fix.txt.vcf_header.vcf.gz

	        # Run eagle
	        ${EAGLE_BIN} \
--vcfRef=${MAIN_DATA_DIR}/sorted_tag_target_snvs_haplocoded.matbed.gz_training.matbed.gz.txt.vcf_header.vcf.gz.bcf \
--vcfTarget=${MAIN_DATA_DIR}/tag_snvs_haplocoded.matbed.gz_testing.matbed.gz.txt_name_fix.txt.vcf_header.vcf.gz \
--geneticMapFile=${EAGLE_MAPS} \
--outPrefix=${MAIN_DATA_DIR}/tag_snvs_haplocoded.matbed.gz_testing.matbed.gz.txt_name_fix.txt.vcf_header.vcf.gz_eagle_phase \
--noImpMissing \
--numThreads=40 2>&1 | tee eagle_phase.op
fi


