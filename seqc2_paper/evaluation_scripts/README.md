# Evaluation Scripts for the SEQC-II deepLearning for somatic mutation detection project
Here you can find evaluation scripts for the following SEQC-II paper:
Sahraeian SME, et al., Achieving Robust Somatic Mutation Detection with Deep Learning Models Derived from Reference Data Sets of a Cancer Sample, Genome Biology, 2022.

## Pre requisite

You need to install [VarSim](https://github.com/bioinform/varsim)

## Example Usage

```
./evaluate.sh \
-r GRCh38.d1.vd1.fa \
-t data/HighMed.v1.0.vcf.gz \
-l data/Low.v1.0.vcf.gz \
-q somatic.vcf \ # predicted somatic VCF
-u tumor.bam \
-b data/test_region_WGS.bed \
-p 100 \ # Tumor purity
-c 80 \ # Tumor average coverage
-w work
```
