#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o: --long out-dir:,tumor-bam:,normal-bam:,human-reference:,selector:,exclude:,dbsnp:,cosmic:,MEM:,action:,mutect2:,varscan-snv:,varscan-indel:,sniper:,vardict:,muse:,strelka-snv:,strelka-indel:,extra-arguments: -n 'submit_SomaticSeq.sh'  -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

#echo "$OPTS"
eval set -- "$OPTS"

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

action=echo
MEM=6

while true; do
    case "$1" in

    --out-dir )
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

    --MEM )
        case "$2" in
            "") shift 2 ;;
            *) MEM=$2 ; shift 2 ;;
        esac ;;

    --action )
        case "$2" in
            "") shift 2 ;;
            *)  action=$2 ; shift 2 ;;
        esac ;;

    --mutect2 )
        case "$2" in
            "") shift 2 ;;
            *)  mutect2_vcf=$2 ; shift 2 ;;
        esac ;;

    --varscan-snv )
        case "$2" in
            "") shift 2 ;;
            *)  varscan_snv_vcf=$2 ; shift 2 ;;
        esac ;;

    --varscan-indel )
        case "$2" in
            "") shift 2 ;;
            *)  varscan_indel_vcf=$2 ; shift 2 ;;
        esac ;;

    --sniper )
        case "$2" in
            "") shift 2 ;;
            *)  sniper_vcf=$2 ; shift 2 ;;
        esac ;;

    --vardict )
        case "$2" in
            "") shift 2 ;;
            *)  vardict_vcf=$2 ; shift 2 ;;
        esac ;;

    --muse )
        case "$2" in
            "") shift 2 ;;
             *)  muse_vcf=$2 ; shift 2 ;;
        esac ;;

    --strelka-snv )
        case "$2" in
            "") shift 2 ;;
            *)  strelka_snv_vcf=$2 ; shift 2 ;
        esac ;;

    --strelka-indel )
        case "$2" in
            "") shift 2 ;;
            *)  strelka_indel_vcf=$2 ; shift 2 ;;
        esac ;;

    --extra-arguments )
        case "$2" in
            "") shift 2 ;;
            *) extra_arguments=$2 ; shift 2 ;;
        esac ;;
        
    -- ) shift; break ;;
    * ) break ;;
    esac

done


logdir=${outdir}/logs
mkdir -p ${logdir}

out_script=${outdir}/logs/wrapper.cmd

selector_text=''
if [[ -r $SELECTOR ]]; then
    selector_text="--inclusion-region /mnt/${SELECTOR}"
fi

exclusion_text=''
if [[ -r $EXCLUSION ]]; then
    exclusion_text="--exclusion-region /mnt/${EXCLUSION}"
fi

dbsnp_text=''
if [[ -r $dbsnp ]]; then
    dbsnp_text="--dbsnp /mnt/${dbsnp}"
fi

cosmic_text=''
if [[ -r $cosmic ]]; then
    cosmic_text="--cosmic /mnt/${cosmic}"
fi

# VCF inputs
if [[ $mutect2_vcf ]];       then mutect2_text="--mutect2 /mnt/${mutect2_vcf}";                   fi
if [[ $varscan_snv_vcf ]];   then varscan_snv_text="--varscan-snv /mnt/${varscan_snv_vcf}";       fi
if [[ $varscan_indel_vcf ]]; then varscan_indel_text="--varscan-indel /mnt/${varscan_indel_vcf}"; fi
if [[ $sniper_vcf ]];        then sniper_text="--sniper /mnt/${sniper_vcf}";                      fi
if [[ $vardict_vcf ]];       then vardict_text="--vardict /mnt/${vardict_vcf}";                   fi
if [[ $muse_vcf ]];          then muse_text="--muse /mnt/${muse_vcf}";                            fi
if [[ $strelka_snv_vcf ]];   then strelka_snv_text="--strelka-snv /mnt/${strelka_snv_vcf}";       fi
if [[ $strelka_indel_vcf ]]; then strelka_indel_text="--strelka-indel /mnt/${strelka_indel_vcf}"; fi

echo "#!/bin/bash" > $out_script
echo "" >> $out_script

echo "#$ -o ${logdir}" >> $out_script
echo "#$ -e ${logdir}" >> $out_script
echo "#$ -S /bin/bash" >> $out_script
echo "#$ -l h_vmem=${MEM}G" >> $out_script # ML may require lots of RAM
echo 'set -e' >> $out_script
echo "" >> $out_script

echo 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $out_script
echo "" >> $out_script

echo "docker pull lethalfang/somaticseq:2.7.0" >> $out_script
echo "" >> $out_script

echo "docker run --rm -v /:/mnt -u $UID --memory 24g lethalfang/somaticseq:2.7.0 \\" >> $out_script
echo "/opt/somaticseq/SomaticSeq.Wrapper.sh \\" >> $out_script
echo "--output-dir       /mnt/${outdir} \\" >> $out_script
echo "--genome-reference /mnt/${HUMAN_REFERENCE} \\" >> $out_script
echo "--tumor-bam        /mnt/${tumor_bam} \\" >> $out_script
echo "--normal-bam       /mnt/${normal_bam} \\" >> $out_script
echo "$mutect2_text \\" >> $out_script
echo "$varscan_snv_text \\" >> $out_script
echo "$varscan_indel_text \\" >> $out_script
echo "$sniper_text \\" >> $out_script
echo "$vardict_text \\" >> $out_script
echo "$muse_text \\" >> $out_script
echo "$strelka_snv_text \\" >> $out_script
echo "$strelka_indel_text \\" >> $out_script
echo "$selector_text \\" >> $out_script
echo "$exclusion_text \\" >> $out_script
echo "$cosmic_text \\" >> $out_script
echo "$dbsnp_text \\" >> $out_script
echo "${extra_arguments} \\" >> $out_script
echo "--gatk /opt/GATK/GenomeAnalysisTK.jar" >> $out_script

echo "" >> $out_script
echo 'echo -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $out_script

$action $out_script
