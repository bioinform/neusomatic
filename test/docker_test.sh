#!/bin/bash
set -e

test_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
cd ${test_dir}
mkdir -p example
cd example

if [ ! -f Homo_sapiens.GRCh37.75.dna.chromosome.22.fa ]
then
	if [ ! -f Homo_sapiens.GRCh37.75.dna.chromosome.22.fa.gz ]
	then
		docker run -v ${test_dir}:/mnt -u $UID --memory 30G  msahraeian/neusomatic:0.1.4 /bin/bash -c \
		"cd /mnt/example/ && wget ftp://ftp.ensembl.org/pub/release-75//fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.chromosome.22.fa.gz"
	fi
	docker run -v ${test_dir}:/mnt -u $UID --memory 30G  msahraeian/neusomatic:0.1.4 /bin/bash -c \
	"cd /mnt/example/ && gunzip -f Homo_sapiens.GRCh37.75.dna.chromosome.22.fa.gz"
	
fi
if [ ! -f Homo_sapiens.GRCh37.75.dna.chromosome.22.fa.fai ]
then
	docker run -v ${test_dir}:/mnt -u $UID --memory 30G  msahraeian/neusomatic:0.1.4 /bin/bash -c \
	"samtools faidx /mnt/example/Homo_sapiens.GRCh37.75.dna.chromosome.22.fa"
fi
rm -rf work_standalone



#Stand-alone NeuSomatic test 
docker run -v ${test_dir}:/mnt -u $UID --memory 30G  msahraeian/neusomatic:0.1.4 /bin/bash -c \
"python /opt/neusomatic/neusomatic/python/preprocess.py \
	--mode call \
	--reference /mnt/example/Homo_sapiens.GRCh37.75.dna.chromosome.22.fa \
	--region_bed /mnt/region.bed \
	--tumor_bam /mnt/tumor.bam \
	--normal_bam /mnt/normal.bam \
	--work /mnt/example/work_standalone \
	--scan_maf 0.05 \
	--min_mapq 10 \
	--snp_min_af 0.05 \
	--snp_min_bq 20 \
	--snp_min_ao 10 \
	--ins_min_af 0.05 \
	--del_min_af 0.05 \
	--num_threads 1 \
	--scan_alignments_binary /opt/neusomatic/neusomatic/bin/scan_alignments"

docker run -v ${test_dir}:/mnt -u $UID --memory 30G --shm-size 8G msahraeian/neusomatic:0.1.4 /bin/bash -c \
"CUDA_VISIBLE_DEVICES= python /opt/neusomatic/neusomatic/python/call.py \
		--candidates_tsv /mnt/example/work_standalone/dataset/*/candidates*.tsv \
		--reference /mnt/example/Homo_sapiens.GRCh37.75.dna.chromosome.22.fa \
		--out /mnt/example/work_standalone \
		--checkpoint /opt/neusomatic/neusomatic/models/NeuSomatic_v0.1.0_standalone_Dream3_70purity.pth \
		--num_threads 1 \
		--batch_size 100"

docker run -v ${test_dir}:/mnt -u $UID --memory 30G  msahraeian/neusomatic:0.1.4 /bin/bash -c \
"python /opt/neusomatic/neusomatic/python/postprocess.py \
		--reference /mnt/example/Homo_sapiens.GRCh37.75.dna.chromosome.22.fa \
		--tumor_bam /mnt/tumor.bam \
		--pred_vcf /mnt/example/work_standalone/pred.vcf \
		--candidates_vcf /mnt/example/work_standalone/work_tumor/filtered_candidates.vcf \
		--output_vcf /mnt/example/work_standalone/NeuSomatic_standalone.vcf \
		--work /mnt/example/work_standalone "


rm -rf /mnt/example/work_ensemble
#Ensemble NeuSomatic test 
docker run -v ${test_dir}:/mnt -u $UID --memory 30G  msahraeian/neusomatic:0.1.4 /bin/bash -c \
"python /opt/neusomatic/neusomatic/python/preprocess.py \
	--mode call \
	--reference /mnt/example/Homo_sapiens.GRCh37.75.dna.chromosome.22.fa \
	--region_bed /mnt/region.bed \
	--tumor_bam /mnt/tumor.bam \
	--normal_bam /mnt/normal.bam \
	--work /mnt/example/work_ensemble \
	--scan_maf 0.05 \
	--min_mapq 10 \
	--snp_min_af 0.05 \
	--snp_min_bq 20 \
	--snp_min_ao 10 \
	--ins_min_af 0.05 \
	--del_min_af 0.05 \
	--num_threads 1 \
	--ensemble_tsv /mnt/ensemble.tsv \
	--scan_alignments_binary /opt/neusomatic/neusomatic/bin/scan_alignments"

docker run -v ${test_dir}:/mnt -u $UID --memory 30G --shm-size 8G msahraeian/neusomatic:0.1.4 /bin/bash -c \
"CUDA_VISIBLE_DEVICES= python /opt/neusomatic/neusomatic/python/call.py \
		--candidates_tsv /mnt/example/work_ensemble/dataset/*/candidates*.tsv \
		--reference /mnt/example/Homo_sapiens.GRCh37.75.dna.chromosome.22.fa \
		--out /mnt/example/work_ensemble \
		--checkpoint /opt/neusomatic/neusomatic/models/NeuSomatic_v0.1.0_ensemble_Dream3_70purity.pth \
		--num_threads 1 \
		--ensemble \
		--batch_size 100"

docker run -v ${test_dir}:/mnt -u $UID --memory 30G  msahraeian/neusomatic:0.1.4 /bin/bash -c \
"python /opt/neusomatic/neusomatic/python/postprocess.py \
		--reference /mnt/example/Homo_sapiens.GRCh37.75.dna.chromosome.22.fa \
		--tumor_bam /mnt/tumor.bam \
		--pred_vcf /mnt/example/work_ensemble/pred.vcf \
		--candidates_vcf /mnt/example/work_ensemble/work_tumor/filtered_candidates.vcf \
		--ensemble_tsv /mnt/ensemble.tsv \
		--output_vcf /mnt/example/work_ensemble/NeuSomatic_ensemble.vcf \
		--work /mnt/example/work_ensemble" 

cd ..

file1=${test_dir}/example/work_standalone/NeuSomatic_standalone.vcf
file2=${test_dir}/NeuSomatic_standalone.vcf

cmp --silent $file1 $file2 && echo "### NeuSomatic stand-alone: SUCCESS! ###" \
|| echo "### NeuSomatic stand-alone FAILED: Files ${file1} and ${file2} Are Different! ###"


file1=${test_dir}/example/work_ensemble/NeuSomatic_ensemble.vcf
file2=${test_dir}/NeuSomatic_ensemble.vcf

cmp --silent $file1 $file2 && echo "### NeuSomatic ensemble: SUCCESS! ###" \
|| echo "### NeuSomatic ensemble FAILED: Files ${file1} and ${file2} Are Different! ###"
