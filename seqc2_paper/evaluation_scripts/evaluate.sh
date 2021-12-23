#!/bin/bash

# Purpose: Perform evaluation on SEQC-II data and run varsim compare with fixed query VCF
# Usage:   bash evaluate_seqc2.sh -r <ref-fasta> -t <truthset-vcf> -l <low-confidence vcf> -b <region-bed>
#						    -q <query-vcf> -u <tumor-bam> -p <tumor-purity> -c <coverage>
# Date:    2021

set -e
# Set defaults --- {{{
# ref=GRCh38.d1.vd1.fa (https://api.gdc.cancer.gov/data/254f697d-310d-4d7d-a27b-27fbf767a834)
# truth=HighMed.v1.0.vcf
# low_truth=Low.v1.0.vcf
# }}} ---

# Help message --- {{{
usage () {
	echo -e \\n"Fix evaluation on SEQC-II data"\\n
	echo -e """Usage: $0 [-r ref.fa] [-w work] [-t truth.vcf] [-l lowConf.vcf]
                   [-q in.vcf] [-u tumor.bam] [-b region.bed]
                   [-p NUM] [-c NUM]"""\\n
	echo "Options:"
	echo "  -r FILE        Reference FASTA file"
	echo "  -w PATH        work folder"
	echo "  -t FILE        Truth-set VCF.gz file (HighMed.v1.0.vcf.gz)"
	echo "  -l FILE        VCF.gz file containing low confidence variants (Low.v1.0.vcf.gz)"
	echo "  -q FILE        Query VCF file. This VCF file should contain total of 10 columns i.e., only 1 sample column"
	echo "  -u FILE        Tumor sample BAM file"
	echo "  -b FILE        Region BED file"
	echo "  -p NUM         Tumor purity"
	echo "  -c NUM         Coverage"
	echo -e "  -h             Prints help message"\\n
}
# }}} ---

# getopts --- {{{
while getopts r:w:t:l:q:u:b:p:c:h flags
do
	case "${flags}" in
		r) ref=${OPTARG};;
		w) work=${OPTARG};;
		t) truth=${OPTARG};;
		l) low_truth=${OPTARG};;
		q) vcf=${OPTARG};;
		u) tumor_bam=${OPTARG};;
		b) bed=${OPTARG};;
		p) Purity=${OPTARG};;
		c) coverage=${OPTARG};;
		h) usage; exit 2;;
		?) echo -e \\n"Use -h to see the help documentation."\\n; exit 2;;
	esac
done
# }}} ---


if [[ -z $ref ]]; then echo -e \\n"[ERRR] Required -r arguments missing."; usage; exit 1; fi
if [[ -z $work ]]; then echo -e \\n"[ERRR] Required -w arguments missing."; usage; exit 1; fi
if [[ -z $truth ]]; then echo -e \\n"[ERRR] Required -t arguments missing."; usage; exit 1; fi
if [[ -z $low_truth ]]; then echo -e \\n"[ERRR] Required -l arguments missing."; usage; exit 1; fi
if [[ -z $vcf ]]; then echo -e \\n"[ERRR] Required -q arguments missing."; usage; exit 1; fi
if [[ -z $tumor_bam ]]; then echo -e \\n"[ERRR] Required -u arguments missing."; usage; exit 1; fi
if [[ -z $bed ]]; then echo -e \\n"[ERRR] Required -b arguments missing."; usage; exit 1; fi
if [[ -z $Purity ]]; then echo -e \\n"[ERRR] Required -p arguments missing."; usage; exit 1; fi
if [[ -z $coverage ]]; then echo -e \\n"[ERRR] Required -c arguments missing."; usage; exit 1; fi

# Tracking --- {{{
echo """
Query:            ${vcf}
Tumor bam:        ${tumor_bam}
bed:              $bed
Purity, coverage: ${Purity}% , ${coverage}x
----------------------------------------------------"""
# }}} ---

# Fix evaluation and run som.py --- {{{{
mkdir -p ${work}/tmp/
python varsim/compare_vcf.py \
--reference ${ref} \
--out_dir ${work}/tmp \
--vcfs ${vcf} \
--true_vcf ${truth} \
--regions ${bed} \
--vcfcompare_options "-wig 1" --sv_length 100 --var_types SNP Insertion Complex Deletion TandemDup \
2>&1|tee ${work}/tmp/eval.log


fp=${work}/tmp/augmented_fp.vcf.gz
tp=${work}/tmp/augmented_tp.vcf.gz
info_extracted=${work}/tmp/fp.info_extracted.vcf

python extract_info_bam.py \
--reference ${ref} \
--variants ${fp} \
--bam ${tumor_bam} \
--num_threads 6 \
--output ${info_extracted}


filtered_fp=${work}/tmp/fp.dp3.vcf
if (( ((50*100)>(${coverage}*${Purity}*3)) ))
then
zcat ${fp} > ${filtered_fp}
else
cat <(zcat ${fp}|grep "#") \
<(bedtools window -a ${fp} -b ${info_extracted} \
|awk -v Purity=${Purity} -v coverage=${coverage} '(($2==$12)&&($4==$14)&&($5==$15)&&((($17*50*100)>=(coverage*Purity*3)))){print}'|cut -f 1,2,3,4,5,6,7,8,9,10
) > ${filtered_fp}
fi


filtered_fp_low=${work}/tmp/fp.dp3.low.vcf
python exclude_from_vcf.py \
--input ${filtered_fp} \
--ex ${low_truth} \
--output ${filtered_fp_low}

filtered_vcf=${work}/tmp/dp3.low.vcf
cat <(cat ${vcf}|grep "#") \
<(zcat  ${tp}| cat - ${filtered_fp_low}| grep -v "#"|bedtools window -a ${vcf} -b - | \
 awk '(($2==$12)&&($4==$14)&&($5==$15)){print}'|cut -f 1,2,3,4,5,6,7,8,9,10) \
 > ${filtered_vcf}

mkdir -p ${work}/final_eval/
python varsim/compare_vcf.py \
--reference ${ref} \
--out_dir ${work}/final_eval \
--vcfs ${filtered_vcf} \
--true_vcf ${truth} \
--regions ${bed} \
--vcfcompare_options "-wig 1" --sv_length 100 --var_types SNP Insertion Complex Deletion TandemDup \
2>&1|tee ${work}/final_eval/eval.log

echo "----------------------------------------------------"
# }}} ---
