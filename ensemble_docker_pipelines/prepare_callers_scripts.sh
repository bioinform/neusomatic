#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o: --long output-dir:,tumor-bam:,normal-bam:,tumor-name:,normal-name:,human-reference:,selector:,exclude:,dbsnp:,cosmic:,min-vaf:,splits:,mutect2,varscan2,somaticsniper,vardict,muse,strelka,wrapper,exome,mutect2-arguments:,mutect2-filter-arguments:,varscan-arguments:,varscan-pileup-arguments:,somaticsniper-arguments:,vardict-arguments:,muse-arguments:,strelka-config-arguments:,strelka-run-arguments:,wrapper-arguments:, -n 'prepare_callers_scripts.sh'  -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

echo "$OPTS"
eval set -- "$OPTS"

MYDIR="$( cd "$( dirname "$0" )" && pwd )"


tumor_name='TUMOR'
normal_name='NORMAL'
wrapper_dir='Wrapper'
min_vaf=0.05
splits=12
muse_extra_arguments=''
action='echo'

while true; do
    case "$1" in
    -o | --output-dir )
        case "$2" in
            "") shift 2 ;;
            *)  outdir=$2 ; shift 2 ;;
    esac ;;

    --tumor-bam )
        case "$2" in
            "") shift 2 ;;
            *)  tumor_bam=$2 ; shift 2 ;;
        esac ;;

    --normal-bam )
        case "$2" in
            "") shift 2 ;;
            *)  normal_bam=$2 ; shift 2 ;;
        esac ;;

    --tumor-name )
        case "$2" in
            "") shift 2 ;;
            *)  tumor_name=$2 ; shift 2 ;;
        esac ;;

    --normal-name )
        case "$2" in
            "") shift 2 ;;
            *)  normal_name=$2 ; shift 2 ;;
        esac ;;

    --human-reference )
        case "$2" in
            "") shift 2 ;;
            *)  HUMAN_REFERENCE=$2 ; shift 2 ;;
        esac ;;

    --selector )
        case "$2" in
            "") shift 2 ;;
            *)  SELECTOR=$2 ; shift 2 ;;
        esac ;;

    --exclude )
        case "$2" in
            "") shift 2 ;;
            *)  EXCLUSION=$2 ; shift 2 ;;
        esac ;;

    --dbsnp )
        case "$2" in
            "") shift 2 ;;
            *)  dbsnp=$2 ; shift 2 ;;
        esac ;;

    --cosmic )
        case "$2" in
            "") shift 2 ;;
            *)  cosmic=$2 ; shift 2 ;;
        esac ;;

    --min-vaf )
        case "$2" in
            "") shift 2 ;;
            *)  min_vaf=$2 ; shift 2 ;;
        esac ;;
                
    --splits )
        case "$2" in
            "") shift 2 ;;
            *)  splits=$2 ; shift 2 ;;
        esac ;;

    --mutect2 )
            mutect2=1 ; shift ;;

    --varscan2 )
            varscan2=1 ; shift ;;

    --somaticsniper )
            somaticsniper=1 ; shift ;;

    --vardict )
            vardict=1 ; shift ;;

    --muse )
            muse=1 ; shift ;;

    --strelka )
            strelka=1 ; shift ;;

    --wrapper )
            wrapper=1 ; shift ;;
        
    --mutect2-arguments )
        case "$2" in
            "") shift 2 ;;
            *)  mutect2_arguments=$2 ; shift 2 ;;
        esac ;;

    --mutect2-filter-arguments )
        case "$2" in
            "") shift 2 ;;
            *)  mutect2_filter_arguments=$2 ; shift 2 ;;
        esac ;;

    --varscan-arguments )
        case "$2" in
            "") shift 2 ;;
            *)  varscan_arguments=$2 ; shift 2 ;;
        esac ;;

    --varscan-pileup-arguments )
        case "$2" in
            "") shift 2 ;;
            *)  varscan_pileup_arguments=$2 ; shift 2 ;;
        esac ;;

    --somaticsniper-arguments )
        case "$2" in
            "") shift 2 ;;
            *)  somaticsniper_arguments=$2 ; shift 2 ;;
        esac ;;

    --vardict-arguments )
        case "$2" in
            "") shift 2 ;;
            *)  vardict_arguments=$2 ; shift 2 ;;
        esac ;;

    --muse-arguments )
        case "$2" in
            "") shift 2 ;;
            *)  muse_extra_arguments=$2 ; shift 2 ;;
        esac ;;

    --strelka-config-arguments )
        case "$2" in
            "") shift 2 ;;
            *)  strelka_config_arguments=$2 ; shift 2 ;;
        esac ;;

    --strelka-run-arguments )
        case "$2" in
            "") shift 2 ;;
            *)  strelka_run_arguments=$2 ; shift 2 ;;
        esac ;;

    --wrapper-arguments )
        case "$2" in
            "") shift 2 ;;
            *)  wrapper_arguments=$2 ; shift 2 ;;
        esac ;;

    --exome )
        exome_stat=1 ; shift ;;

    -- ) shift; break ;;
    * ) break ;;
    esac

done


logdir=${outdir}/logs
mkdir -p ${logdir}


if [[ $SELECTOR ]]
then
    docker run --rm -v /:/mnt -u $UID --memory 10G lethalfang/bedtools:2.26.0 bash -c \
    "bedtools sort -i /mnt/$SELECTOR -faidx /mnt/${HUMAN_REFERENCE}.fai  > /mnt/${outdir}/genome.bed"   
else
    cat ${HUMAN_REFERENCE}.fai | awk -F "\t" '{print $1 "\t0\t" $2}' | awk -F "\t" '$1 ~ /^(chr)?[0-9XYMT]+$/' > ${outdir}/genome.bed
fi
    

docker run --rm -v /:/mnt -u $UID lethalfang/somaticseq:2.7.0 \
/opt/somaticseq/utilities/split_Bed_into_equal_regions.py \
-infile /mnt/${outdir}/genome.bed -num $splits -outfiles /mnt/${outdir}/bed

ith_split=1
while [[ $ith_split -le $splits ]]
do

    mkdir -p ${outdir}/${ith_split}
    mv ${outdir}/${ith_split}.bed ${outdir}/${ith_split}

    
    if [[ $somaticsniper -eq 1 ]]
    then
        sniper_input="--sniper ${outdir}/${ith_split}/SomaticSniper.vcf"
    fi
    
    if [[ $mutect2 -eq 1 ]]
    then
    
        mutect2_arguments=''
        input_mutect2_filter_arguments=''
    
        if [[ ${mutect2_arguments} ]]
        then
            input_mutect2_arguments="--extra-arguments ${mutect2_arguments}"
        fi
        
        if [[ ${mutect2_filter_arguments} ]]
        then
            input_mutect2_filter_arguments="--extra-filter-arguments ${mutect2_filter_arguments}"
        fi
    
        $MYDIR/mutation_callers/submit_MuTect2.sh \
        --normal-bam ${normal_bam} \
        --tumor-bam ${tumor_bam} \
        --out-dir ${outdir}/${ith_split} \
        --out-vcf MuTect2.vcf \
        --selector ${outdir}/${ith_split}/${ith_split}.bed \
        --human-reference ${HUMAN_REFERENCE} \
        --dbsnp ${dbsnp} \
        ${input_mutect2_arguments} \
        ${input_mutect2_filter_arguments} \
        --action $action
    
        mutect2_input="--mutect2 ${outdir}/${ith_split}/MuTect2.vcf"
    fi
        

    if [[ $varscan2 -eq 1 ]]
    then
    
        if [[ ${varscan_pileup_arguments} ]]
        then
            input_varscan_pileup_arguments="--extra-pileup-arguments ${varscan_pileup_arguments}"
        fi
        
        if [[ ${varscan_arguments} ]]
        then
            input_varscan_arguments="--extra-arguments ${varscan_arguments}"
        fi
    
        $MYDIR/mutation_callers/submit_VarScan2.sh \
        --normal-bam ${normal_bam} \
        --tumor-bam ${tumor_bam} \
        --out-dir ${outdir}/${ith_split} \
        --out-vcf VarScan2.vcf \
        --selector ${outdir}/${ith_split}/${ith_split}.bed \
        --human-reference ${HUMAN_REFERENCE} \
        ${input_varscan_pileup_arguments} \
        ${input_varscan_arguments} \
        --action $action
    
        varscan_snv_input="--varscan-snv ${outdir}/${ith_split}/VarScan2.snp.vcf"
        varscan_indel_input="--varscan-indel ${outdir}/${ith_split}/VarScan2.indel.vcf"
    fi


        
    if [[ $vardict -eq 1 ]]
    then
    
        if [[ ${vardict_arguments} ]]
        then
            input_vardict_arguments="--extra-arguments ${vardict_arguments}"
        fi
    
        $MYDIR/mutation_callers/submit_VarDictJava.sh \
        --normal-bam ${normal_bam} \
        --tumor-bam ${tumor_bam} \
        --out-dir ${outdir}/${ith_split} \
        --selector ${outdir}/${ith_split}/${ith_split}.bed \
        --out-vcf VarDict.vcf \
        --human-reference ${HUMAN_REFERENCE} \
        --VAF ${min_vaf} \
        ${input_vardict_arguments} \
        --action $action
        
        vardict_input="--vardict ${outdir}/${ith_split}/VarDict.vcf"
    fi
    
    
    if [[ $muse -eq 1 ]]
    then

        if [[ ! ${muse_extra_arguments} ]]
        then
            if [[ $exome_stat ]]
            then
                muse_extra_arguments='-E'
            else
                muse_extra_arguments='-G'
            fi
        fi


        $MYDIR/mutation_callers/submit_MuSE.sh \
        --normal-bam ${normal_bam} \
        --tumor-bam ${tumor_bam} \
        --out-dir ${outdir}/${ith_split} \
        --selector ${outdir}/${ith_split}/${ith_split}.bed \
        --out-vcf MuSE.vcf \
        --human-reference ${HUMAN_REFERENCE} \
        --dbsnp ${dbsnp} \
        --extra-arguments "${muse_extra_arguments}" \
        --action $action
        
        muse_input="--muse ${outdir}/${ith_split}/MuSE.vcf"
    fi
    
    if [[ $strelka -eq 1 ]]
    then
    
        if [[ $exome_stat ]]
        then
            strelka_exome_stat='--exome'
        fi
    
        if [[ ${strelka_config_arguments} ]]
        then
            input_strelka_config_arguments=" --extra-config-arguments ${strelka_config_arguments}"
        fi
        
        if [[ ${strelka_run_arguments} ]]
        then
            input_strelka_run_arguments="--extra-run-arguments ${strelka_run_arguments}"
        fi
    
        $MYDIR/mutation_callers/submit_Strelka.sh \
        --normal-bam ${normal_bam} \
        --tumor-bam ${tumor_bam} \
        --out-dir ${outdir}/${ith_split} \
        --selector ${outdir}/${ith_split}/${ith_split}.bed \
        --out-vcf Strelka.vcf \
        --human-reference ${HUMAN_REFERENCE} \
        ${strelka_exome_stat} \
        ${input_strelka_config_arguments} \
        ${input_strelka_run_arguments} \
        --action $action
        
        strelka_snv_input="--strelka-snv ${outdir}/${ith_split}/Strelka/results/variants/somatic.snvs.vcf.gz"
        strelka_indel_input="--strelka-indel ${outdir}/${ith_split}/Strelka/results/variants/somatic.indels.vcf.gz"        
    fi
    
    
    # Wrapper
    if [[ $wrapper -eq 1 ]]
    then
        if [[ $EXCLUSION ]];        then exclusion_text="--exclude ${EXCLUSION}"                            ; fi

        if [[ ${dbsnp} ]];          then dbsnp_input="--dbsnp ${dbsnp}"                                     ; fi
        if [[ ${cosmic} ]];         then cosmic_input="--cosmic ${cosmic}"                                  ; fi
        
        if [[ ${wrapper_arguments} ]]
        then
            input_wrapper_arguments="--extra-arguments ${wrapper_arguments}"
        fi
        
        $MYDIR/mutation_callers/submit_Wrapper.sh \
        --normal-bam ${normal_bam} \
        --tumor-bam ${tumor_bam} \
        --out-dir ${outdir}/${ith_split}/${wrapper_dir} \
        --human-reference ${HUMAN_REFERENCE} \
        $dbsnp_input \
        $cosmic_input \
        --selector ${outdir}/${ith_split}/${ith_split}.bed \
        $exclusion_text \
        $mutect2_input \
        $varscan_snv_input \
        $varscan_indel_input \
        $sniper_input \
        $vardict_input \
        $muse_input \
        $strelka_snv_input \
        $strelka_indel_input \
        ${input_wrapper_arguments} \
        --action ${action}
    fi
        
    ith_split=$(( $ith_split + 1))

done


# SomaticSniper is very fast, so no need to parallelize
if [[ $somaticsniper -eq 1 ]]
then

    if [[ ${somaticsniper_arguments} ]]
    then
        input_somaticsniper_arguments="--extra-arguments ${somaticsniper_arguments}"
    fi

    $MYDIR/mutation_callers/submit_SomaticSniper.sh \
    --normal-bam ${normal_bam} \
    --tumor-bam ${tumor_bam} \
    --out-dir ${outdir} \
    --out-vcf SomaticSniper.vcf \
    --human-reference ${HUMAN_REFERENCE} \
    --split $splits \
    ${input_somaticsniper_arguments} \
    --action $action
fi
