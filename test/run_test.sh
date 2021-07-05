#!/bin/bash
set -e

test_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
neusomatic_dir="$( dirname ${test_dir} )"

cd ${test_dir}
mkdir -p example
cd example
if [ ! -f Homo_sapiens.GRCh37.75.dna.chromosome.22.fa ]
then
	if [ ! -f Homo_sapiens.GRCh37.75.dna.chromosome.22.fa.gz ]
	then
		wget ftp://ftp.ensembl.org/pub/release-75//fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.chromosome.22.fa.gz
	fi
	gunzip -f Homo_sapiens.GRCh37.75.dna.chromosome.22.fa.gz
fi
if [ ! -f Homo_sapiens.GRCh37.75.dna.chromosome.22.fa.fai ]
then
	samtools faidx Homo_sapiens.GRCh37.75.dna.chromosome.22.fa
fi
rm -rf work_standalone
#Stand-alone NeuSomatic test 
python ${neusomatic_dir}/neusomatic/python/preprocess.py \
	--mode call \
	--reference Homo_sapiens.GRCh37.75.dna.chromosome.22.fa \
	--region_bed ${test_dir}/region.bed \
	--tumor_bam ${test_dir}/tumor.bam \
	--normal_bam ${test_dir}/normal.bam \
	--work work_standalone \
	--scan_maf 0.05 \
	--min_mapq 10 \
	--snp_min_af 0.05 \
	--snp_min_bq 20 \
	--snp_min_ao 10 \
	--ins_min_af 0.05 \
	--del_min_af 0.05 \
	--num_threads 1 \
	--no_seq_complexity \
	--scan_alignments_binary ${neusomatic_dir}/neusomatic/bin/scan_alignments

CUDA_VISIBLE_DEVICES= python ${neusomatic_dir}/neusomatic/python/call.py \
		--candidates_tsv work_standalone/dataset/*/candidates*.tsv \
		--reference Homo_sapiens.GRCh37.75.dna.chromosome.22.fa \
		--out work_standalone \
		--checkpoint ${neusomatic_dir}/neusomatic/models/NeuSomatic_v0.1.0_standalone_Dream3_70purity.pth \
		--num_threads 1 \
		--batch_size 100

python ${neusomatic_dir}/neusomatic/python/postprocess.py \
		--reference Homo_sapiens.GRCh37.75.dna.chromosome.22.fa \
		--tumor_bam ${test_dir}/tumor.bam \
		--pred_vcf work_standalone/pred.vcf \
		--candidates_vcf work_standalone/work_tumor/filtered_candidates.vcf \
		--output_vcf work_standalone/NeuSomatic_standalone.vcf \
		--work work_standalone 


rm -rf work_ensemble
#Ensemble NeuSomatic test 
python ${neusomatic_dir}/neusomatic/python/preprocess.py \
	--mode call \
	--reference Homo_sapiens.GRCh37.75.dna.chromosome.22.fa \
	--region_bed ${test_dir}/region.bed \
	--tumor_bam ${test_dir}/tumor.bam \
	--normal_bam ${test_dir}/normal.bam \
	--work work_ensemble \
	--scan_maf 0.05 \
	--min_mapq 10 \
	--snp_min_af 0.05 \
	--snp_min_bq 20 \
	--snp_min_ao 10 \
	--ins_min_af 0.05 \
	--del_min_af 0.05 \
	--num_threads 1 \
	--ensemble_tsv ${test_dir}/ensemble.tsv \
	--no_seq_complexity \
	--scan_alignments_binary ${neusomatic_dir}/neusomatic/bin/scan_alignments

CUDA_VISIBLE_DEVICES= python ${neusomatic_dir}/neusomatic/python/call.py \
		--candidates_tsv work_ensemble/dataset/*/candidates*.tsv \
		--reference Homo_sapiens.GRCh37.75.dna.chromosome.22.fa \
		--out work_ensemble \
		--checkpoint ${neusomatic_dir}/neusomatic/models/NeuSomatic_v0.1.0_ensemble_Dream3_70purity.pth \
		--num_threads 1 \
		--ensemble \
		--batch_size 100

python ${neusomatic_dir}/neusomatic/python/postprocess.py \
		--reference Homo_sapiens.GRCh37.75.dna.chromosome.22.fa \
		--tumor_bam ${test_dir}/tumor.bam \
		--pred_vcf work_ensemble/pred.vcf \
		--candidates_vcf work_ensemble/work_tumor/filtered_candidates.vcf \
		--ensemble_tsv ${test_dir}/ensemble.tsv \
		--output_vcf work_ensemble/NeuSomatic_ensemble.vcf \
		--work work_ensemble 


cd ..

file1=${test_dir}/example/work_standalone/NeuSomatic_standalone.vcf
file2=${test_dir}/NeuSomatic_standalone.vcf

cmp --silent $file1 $file2 && echo "### NeuSomatic stand-alone: SUCCESS! ###" \
|| echo "### NeuSomatic stand-alone FAILED: Files ${file1} and ${file2} Are Different! ###"


file1=${test_dir}/example/work_ensemble/NeuSomatic_ensemble.vcf
file2=${test_dir}/NeuSomatic_ensemble.vcf

cmp --silent $file1 $file2 && echo "### NeuSomatic ensemble: SUCCESS! ###" \
|| echo "### NeuSomatic ensemble FAILED: Files ${file1} and ${file2} Are Different! ###"
