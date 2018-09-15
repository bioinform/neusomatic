## Docker pipelines to prepare other callers output for Ensemble mode

NeuSomatic can be used universally as a stand-alone somatic mutation detection method or with an ensemble of existing methods. In the ensemble mode NeuSomatic currently supports outputs from MuTect2, MuSE, Strelka2, SomaticSniper, VarDict, and VarScan2. For ensemble mode, the ensembled outputs of different somatic callers (as a single `.tsv` file) should be prepared and inputed using `--ensemble_tsv` argument in `preprocess.py`. 

The following steps use docker pipelines to prepare the ensemble `.tsv` file.

This is an adaptation of SomaticSeq's [scripts](https://github.com/bioinform/somaticseq/tree/master/utilities/dockered_pipelines). 

### Requirements
* Have internet connection and docker daemon. Be able to pull and run docker images from Docker Hub.
* **Highly recommended**: Have cluster management system with valid `qsub` command, such as Sun Grid Engine.


### 1. Prepare run scripts
To run the individual somatic callers (MuTect2, MuSE, Strelka2, SomaticSniper, VarDict, and VarScan2), you can use the following command that create 10 (if we set `splits` to 10) equal-size regions in 10 bed files, and parallelize the jobs into 10 regions.
```
prepare_callers_scripts.sh \
--normal-bam      /ABSOLUTE/PATH/TO/normal_sample.bam \
--tumor-bam       /ABSOLUTE/PATH/TO/tumor_sample.bam \
--human-reference /ABSOLUTE/PATH/TO/ref.fa \
--output-dir      /ABSOLUTE/PATH/TO/output \
--dbsnp           /ABSOLUTE/PATH/TO/dbSNP.vcf \
--splits        10 \
--selector /ABSOLUTE/PATH/TO/region.bed \
--mutect2 --somaticsniper --vardict --varscan2 --muse --strelka --wrapper
```
This command will create sub-folders with the following sctructure for the split regions:
```
output
├── genome.bed
├── 1
├── 2
├── .
├── .
├── .
├── i                 		 # i'th folder
│   ├── i.bed         		 # i'th sub-region
│   ├── logs          
│   │    ├── strelka.cmd 	 # Strelka2 run script for region i.bed
│   │    ├── muse.cmd        # MuSE run script for region i.bed
│   │    ├── mutect2.cmd     # MuTect2 run script for region i.bed
│   │    ├── vardict.cmd     # VarDict run script for region i.bed
│   │    └── varscan2.cm     # VarScan2 run script for region i.bed
│   └── Wrapper         
│       └── logs           	
│            └── wrapper.cmd # Wrapper script to combine results
├── .
├── .
├── .
├── 10
└── logs                
     └── somaticsniper.cmd   # SomaticSniper run script for whole region genome.bed
```

### 2. Run all individual callers scripts
You should first run all individual callers `.cmd` run scripts for all regions. For instance with `qsub` command:
```
for tool in strelka muse mutect2 vardict varscan2
do
	for i in {1..10}
	do
		for script in output/${i}/logs/${tool}.cmd
		do
			qsub $script
		done
	done
done
qsub output/logs/somaticsniper.cmd
```
### 3. Combine individual callers outputs
Once all these scripts finished successfully, the respective VCF files for each tool will be available under each region sub-folder. 

Then, you can run wrapper.cmd scripts to combine calls made by each caller:
```
for i in {1..10}
do
	qsub output/${i}/Wrapper/logs/wrapper.cmd
done
```
This will generate "Ensemble.sSNV.tsv" and "Ensemble.sINDEL.tsv" under each region subfolder (e.g. under 'output/{1,2,3,...}/Wrapper/'). 

Now you can combine these files to generate `enemble_ann.tsv` file.
```
cat <(cat output/*/Wrapper/Ensemble.s*.tsv |grep CHROM|head -1) \
    <(cat output/*/Wrapper/Ensemble.s*.tsv |grep -v CHROM) | sed "s/nan/0/g" > ensemble_ann.tsv

```
and provide `enemble_ann.tsv` as `--enemble_ann` argument in `preprocess.py`.



## Options and Parameters
**prepare_callers_scripts.sh** can prepare dockerized somatic mutation calling jobs. The following options are available:
* `--normal-bam`                  /ABSOLUTE/PATH/TO/normal_sample.bam (Required)
* `--tumor-bam`                   /ABSOLUTE/PATH/TO/tumor_sample.bam  (Required)
* `--output-dir`                  /ABSOLUTE/PATH/TO/OUTPUT_DIRECTORY (Required)
* `--human-reference`             /ABSOLUTE/PATH/TO/human_reference.fa (Required)
* `--dbsnp`                       /ABSOLUTE/PATH/TO/dbsnp.vcf (Required for MuSE and LoFreq)
* `--cosmic`                      /ABSOLUTE/PATH/TO/cosmic.vcf (Optional)
* `--selector`                    /ABSOLUTE/PATH/TO/Capture_region.bed (Optional. Will create genome.bed from the .fa.fai file when not specified.)
* `--exclude`                     /ABSOLUTE/PATH/TO/Blacklist_Region.bed (Optional)
* `--min-af`                      (Optional. The minimum VAF cutoff for VarDict. Defulat is 0.05.)
* `--splits`                     N (Optional: evenly split the genome into N BED files. Default = 12).
* `--mutect2`                     (Optional flag to invoke MuTect2)
* `--varscan2`                    (Optional flag to invoke VarScan2)
* `--somaticsniper`               (Optional flag to invoke SomaticSniper)
* `--vardict`                     (Optional flag to invoke VarDict)
* `--muse`                        (Optional flag to invoke MuSE)
* `--strelka`                     (Optional flag to invoke Strelka)
* `--wrapper`                     (Optional flag to invoke wrapper to combine outputs of all callers).
* `--exome`                       (Optional flag for Strelka and MuSE when data is WES)
* `--mutect2-arguments`           (Extra parameters to pass onto Mutect2, e.g., --mutect2-arguments '--initial_tumor_lod 3.0 --log_somatic_prior -5.0 --min_base_quality_score 20')
* `--mutect2-filter-arguments`    (Extra parameters to pass onto FilterMutectCalls)
* `--varscan-arguments`           (Extra parameters to pass onto VarScan2)
* `--varscan-pileup-arguments`    (Extra parameters to pass onto samtools mpileup that creates pileup files for VarScan)
* `--somaticsniper-arguments`     (Extra parameters to pass onto SomaticSniper)
* `--vardict-arguments`           (Extra parameters to pass onto VarDict)
* `--muse-arguments`              (Extra parameters to pass onto MuSE)
* `--strelka-config-arguments`    (Extra parameters to pass onto Strelka's config command)
* `--strelka-run-arguments`       (Extra parameters to pass onto Strekla's run command)
* `--wrapper-arguments`           (Extra parameters to pass onto SomaticSeq.Wrapper.sh)




### NOTES
* Due to the way those run scripts are written, the Sun Grid Engine's standard error log will record the time the task completes (i.e., `Done at 2017/10/30 29:03:02`), and it will only do so when the task is completed with an exit code of 0. It can be a quick way to check if a task is done, by looking at the final line of the standard error log file.
* If you don't provide a region file, each split BED file represents 1/N (For N splits) of the total base pairs in the human genome (obtained from the .fa.fai file, but only including 1, 2, 3, ..., MT, or chr1, chr2, ..., chrM contigs).
* Parallelization (i.e., splitting) is not turned on for SomaticSniper because 1) it's manageable on a single split (even for WGS), and 2) it doesn't support partial processing with BED file, so it may not be worth the time to split the BAM. After SomaticSniper finishes, the result VCF files will be split into each of the regions sub-folders `output/{1,2,3,...}`.
* After specifying the reference fasta (must have extensions of .fa or .fasta), it must also include the .dict and .fa.fai (or .fasta.fai) files in the same directory.
* When specifying `/ABSOLUTE/PATH/TO/dbSNP.vcf`, there also needs to be `dbSNP.vcf.idx`, `dbSNP.vcf.gz`, and `dbSNP.vcf.gz.tbi` present at the same directory because MuSE and LoFreq are expecting them.
* We also have no distribution rights for VarScan2, so our script points to a 3rd-party version. Only run it if you are licensed to do so.
