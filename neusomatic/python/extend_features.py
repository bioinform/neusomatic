#!/usr/bin/env python
#-------------------------------------------------------------------------
# extend_features.py
# add extra features for standalone mode
#-------------------------------------------------------------------------
import argparse
import traceback
import logging
import multiprocessing
import os
import gzip

import pysam
import numpy as np

import sequencing_features
import genomic_file_handlers as genome
from read_info_extractor import rescale
from utils import skip_empty, get_chromosomes_order


def extract_features(candidate_record):
    reference, tumor_bam, normal_bam, min_mapq, min_bq, dbsnp, no_seq_complexity, batch = candidate_record
    thread_logger = logging.getLogger(
        "{} ({})".format(extract_features.__name__, multiprocessing.current_process().name))
    try:
        tbam = pysam.AlignmentFile(tumor_bam)
        nbam = pysam.AlignmentFile(normal_bam)
        ref_fa = pysam.FastaFile(reference)
        if dbsnp:
            dbsnp_tb = pysam.TabixFile(dbsnp)

        ext_features = []
        for nei_cluster in batch:
            n_cluster_reads = sequencing_features.ClusterReads(nbam, nei_cluster)
            t_cluster_reads = sequencing_features.ClusterReads(tbam, nei_cluster)
            for var_i, [chrom, pos, ref, alt, if_cosmic, num_cosmic_cases] in enumerate(nei_cluster):
                var_id = "-".join([chrom, str(pos), ref, alt])
                pos = int(pos)
                my_coordinate = [chrom, pos]
                nBamFeatures = n_cluster_reads.get_alignment_features(var_i, ref, alt, min_mapq, min_bq)
                tBamFeatures = t_cluster_reads.get_alignment_features(var_i, ref, alt, min_mapq, min_bq)

                sor = sequencing_features.somaticOddRatio(nBamFeatures.nref, nBamFeatures.nalt, tBamFeatures.nref,
                                                          tBamFeatures.nalt)

                homopolymer_length, site_homopolymer_length = sequencing_features.from_genome_reference(
                    ref_fa, my_coordinate, ref, alt)

                indel_length = len(alt) - len(ref)

                if not no_seq_complexity:
                    seq_span_80bp = ref_fa.fetch(my_coordinate[0], max(
                        0, my_coordinate[1] - 41), my_coordinate[1] + 40)
                    seq_left_80bp = ref_fa.fetch(my_coordinate[0], max(
                        0, my_coordinate[1] - 81), my_coordinate[1])
                    seq_right_80bp = ref_fa.fetch(my_coordinate[0], my_coordinate[
                                                  1], my_coordinate[1] + 81)
                    LC_spanning = sequencing_features.subLC(seq_span_80bp, 20)
                    LC_adjacent = min(sequencing_features.subLC(
                        seq_left_80bp, 20), sequencing_features.subLC(seq_right_80bp, 20))
                    LC_spanning_phred = genome.p2phred(1 - LC_spanning, 40)
                    LC_adjacent_phred = genome.p2phred(1 - LC_adjacent, 40)

                if_dbsnp = 0
                if_common = 0
                if dbsnp:
                    region = "{}:{}-{}".format(chrom, pos, pos + 1)
                    dbsnp_vars = {}
                    for x in dbsnp_tb.fetch(region=region):
                        chrom_, pos_, _, ref_, alts_, _, _, info_ = x.strip().split("\t")[
                            0:8]
                        for alt_ in alts_.split(","):
                            dbsnp_var_id = "-".join([chrom_, pos_, ref_, alt_])
                            dbsnp_vars[
                                dbsnp_var_id] = 1 if "COMMON=1" in info_ else 0
                    if var_id in dbsnp_vars:
                        if_dbsnp = 1
                        if_common = dbsnp_vars[var_id]

                p_scale = None
                CHROM = my_coordinate[0]
                POS = my_coordinate[1]
                REF = ref
                ALT = alt
                if_dbsnp = if_dbsnp
                COMMON = if_common
                if_COSMIC = if_cosmic
                COSMIC_CNT = num_cosmic_cases
                Consistent_Mates = tBamFeatures.consistent_mates
                Inconsistent_Mates = tBamFeatures.inconsistent_mates
                if not no_seq_complexity:
                    Seq_Complexity_Span = LC_spanning_phred
                    Seq_Complexity_Adj = LC_adjacent_phred

                N_DP = nBamFeatures.dp
                nBAM_REF_MQ = '%g' % nBamFeatures.ref_mq
                nBAM_ALT_MQ = '%g' % nBamFeatures.alt_mq
                nBAM_Z_Ranksums_MQ = '%g' % nBamFeatures.z_ranksums_mq
                nBAM_REF_BQ = '%g' % nBamFeatures.ref_bq
                nBAM_ALT_BQ = '%g' % nBamFeatures.alt_bq
                nBAM_Z_Ranksums_BQ = '%g' % nBamFeatures.z_ranksums_bq
                nBAM_REF_NM = '%g' % nBamFeatures.ref_NM
                nBAM_ALT_NM = '%g' % nBamFeatures.alt_NM
                nBAM_NM_Diff = '%g' % nBamFeatures.NM_Diff
                nBAM_REF_Concordant = nBamFeatures.ref_concordant_reads
                nBAM_REF_Discordant = nBamFeatures.ref_discordant_reads
                nBAM_ALT_Concordant = nBamFeatures.alt_concordant_reads
                nBAM_ALT_Discordant = nBamFeatures.alt_discordant_reads
                nBAM_Concordance_FET = rescale(
                    nBamFeatures.concordance_fet, 'fraction', p_scale, 1001)
                N_REF_FOR = nBamFeatures.ref_for
                N_REF_REV = nBamFeatures.ref_rev
                N_ALT_FOR = nBamFeatures.alt_for
                N_ALT_REV = nBamFeatures.alt_rev
                nBAM_StrandBias_FET = rescale(
                    nBamFeatures.strandbias_fet, 'fraction', p_scale, 1001)
                nBAM_Z_Ranksums_EndPos = '%g' % nBamFeatures.z_ranksums_endpos
                nBAM_REF_Clipped_Reads = nBamFeatures.ref_SC_reads
                nBAM_ALT_Clipped_Reads = nBamFeatures.alt_SC_reads
                nBAM_Clipping_FET = rescale(
                    nBamFeatures.clipping_fet, 'fraction', p_scale, 1001)
                nBAM_MQ0 = nBamFeatures.MQ0
                nBAM_Other_Reads = nBamFeatures.noise_read_count
                nBAM_Poor_Reads = nBamFeatures.poor_read_count
                nBAM_REF_InDel_3bp = nBamFeatures.ref_indel_3bp
                nBAM_REF_InDel_2bp = nBamFeatures.ref_indel_2bp
                nBAM_REF_InDel_1bp = nBamFeatures.ref_indel_1bp
                nBAM_ALT_InDel_3bp = nBamFeatures.alt_indel_3bp
                nBAM_ALT_InDel_2bp = nBamFeatures.alt_indel_2bp
                nBAM_ALT_InDel_1bp = nBamFeatures.alt_indel_1bp
                SOR = sor
                MaxHomopolymer_Length = homopolymer_length
                SiteHomopolymer_Length = site_homopolymer_length
                T_DP = tBamFeatures.dp
                tBAM_REF_MQ = '%g' % tBamFeatures.ref_mq
                tBAM_ALT_MQ = '%g' % tBamFeatures.alt_mq
                tBAM_Z_Ranksums_MQ = '%g' % tBamFeatures.z_ranksums_mq
                tBAM_REF_BQ = '%g' % tBamFeatures.ref_bq
                tBAM_ALT_BQ = '%g' % tBamFeatures.alt_bq
                tBAM_Z_Ranksums_BQ = '%g' % tBamFeatures.z_ranksums_bq
                tBAM_REF_NM = '%g' % tBamFeatures.ref_NM
                tBAM_ALT_NM = '%g' % tBamFeatures.alt_NM
                tBAM_NM_Diff = '%g' % tBamFeatures.NM_Diff
                tBAM_REF_Concordant = tBamFeatures.ref_concordant_reads
                tBAM_REF_Discordant = tBamFeatures.ref_discordant_reads
                tBAM_ALT_Concordant = tBamFeatures.alt_concordant_reads
                tBAM_ALT_Discordant = tBamFeatures.alt_discordant_reads
                tBAM_Concordance_FET = rescale(
                    tBamFeatures.concordance_fet, 'fraction', p_scale, 1001)
                T_REF_FOR = tBamFeatures.ref_for
                T_REF_REV = tBamFeatures.ref_rev
                T_ALT_FOR = tBamFeatures.alt_for
                T_ALT_REV = tBamFeatures.alt_rev
                tBAM_StrandBias_FET = rescale(
                    tBamFeatures.strandbias_fet, 'fraction', p_scale, 1001)
                tBAM_Z_Ranksums_EndPos = '%g' % tBamFeatures.z_ranksums_endpos
                tBAM_REF_Clipped_Reads = tBamFeatures.ref_SC_reads
                tBAM_ALT_Clipped_Reads = tBamFeatures.alt_SC_reads
                tBAM_Clipping_FET = rescale(
                    tBamFeatures.clipping_fet, 'fraction', p_scale, 1001)
                tBAM_MQ0 = tBamFeatures.MQ0
                tBAM_Other_Reads = tBamFeatures.noise_read_count
                tBAM_Poor_Reads = tBamFeatures.poor_read_count
                tBAM_REF_InDel_3bp = tBamFeatures.ref_indel_3bp
                tBAM_REF_InDel_2bp = tBamFeatures.ref_indel_2bp
                tBAM_REF_InDel_1bp = tBamFeatures.ref_indel_1bp
                tBAM_ALT_InDel_3bp = tBamFeatures.alt_indel_3bp
                tBAM_ALT_InDel_2bp = tBamFeatures.alt_indel_2bp
                tBAM_ALT_InDel_1bp = tBamFeatures.alt_indel_1bp
                InDel_Length = indel_length

                features = [CHROM, POS, ".", REF, ALT, if_dbsnp, COMMON, if_COSMIC, COSMIC_CNT,
                            Consistent_Mates, Inconsistent_Mates]
                if not no_seq_complexity:
                    features.extend([Seq_Complexity_Span, Seq_Complexity_Adj])
                features.extend([N_DP, nBAM_REF_MQ, nBAM_ALT_MQ, nBAM_Z_Ranksums_MQ,
                                 nBAM_REF_BQ, nBAM_ALT_BQ, nBAM_Z_Ranksums_BQ, nBAM_REF_NM, nBAM_ALT_NM, nBAM_NM_Diff,
                                 nBAM_REF_Concordant, nBAM_REF_Discordant, nBAM_ALT_Concordant, nBAM_ALT_Discordant,
                                 nBAM_Concordance_FET, N_REF_FOR, N_REF_REV, N_ALT_FOR, N_ALT_REV, nBAM_StrandBias_FET,
                                 nBAM_Z_Ranksums_EndPos, nBAM_REF_Clipped_Reads, nBAM_ALT_Clipped_Reads, nBAM_Clipping_FET,
                                 nBAM_MQ0, nBAM_Other_Reads, nBAM_Poor_Reads, nBAM_REF_InDel_3bp, nBAM_REF_InDel_2bp,
                                 nBAM_REF_InDel_1bp, nBAM_ALT_InDel_3bp, nBAM_ALT_InDel_2bp, nBAM_ALT_InDel_1bp, SOR,
                                 MaxHomopolymer_Length, SiteHomopolymer_Length, T_DP, tBAM_REF_MQ, tBAM_ALT_MQ, tBAM_Z_Ranksums_MQ,
                                 tBAM_REF_BQ, tBAM_ALT_BQ, tBAM_Z_Ranksums_BQ, tBAM_REF_NM, tBAM_ALT_NM, tBAM_NM_Diff,
                                 tBAM_REF_Concordant, tBAM_REF_Discordant, tBAM_ALT_Concordant, tBAM_ALT_Discordant,
                                 tBAM_Concordance_FET, T_REF_FOR, T_REF_REV, T_ALT_FOR, T_ALT_REV, tBAM_StrandBias_FET,
                                 tBAM_Z_Ranksums_EndPos, tBAM_REF_Clipped_Reads, tBAM_ALT_Clipped_Reads, tBAM_Clipping_FET,
                                 tBAM_MQ0, tBAM_Other_Reads, tBAM_Poor_Reads, tBAM_REF_InDel_3bp, tBAM_REF_InDel_2bp,
                                 tBAM_REF_InDel_1bp, tBAM_ALT_InDel_3bp, tBAM_ALT_InDel_2bp, tBAM_ALT_InDel_1bp, InDel_Length])

                ext_features.append(features)
        return ext_features

    except Exception as ex:
        thread_logger.error(traceback.format_exc())
        thread_logger.error(ex)
        return None


def extend_features(candidates_vcf,
                    exclude_variants,
                    add_variants,
                    output_tsv,
                    reference, tumor_bam, normal_bam,
                    min_mapq, min_bq,
                    dbsnp, cosmic,
                    no_seq_complexity,
                    window_extend,
                    num_threads):

    logger = logging.getLogger(extend_features.__name__)

    logger.info(
        "----------------------Extend Standalone Features------------------------")

    if not os.path.exists(tumor_bam):
        logger.error("Aborting!")
        raise Exception("No tumor BAM file {}".format(tumor_bam))
    if not os.path.exists(normal_bam):
        logger.error("Aborting!")
        raise Exception("No normal BAM file {}".format(normal_bam))
    if not os.path.exists(tumor_bam + ".bai"):
        logger.error("Aborting!")
        raise Exception(
            "No tumor .bai index file {}".format(tumor_bam + ".bai"))
    if not os.path.exists(normal_bam + ".bai"):
        logger.error("Aborting!")
        raise Exception(
            "No normal .bai index file {}".format(normal_bam + ".bai"))

    if dbsnp:
        if not os.path.exists(dbsnp):
            logger.error("Aborting!")
            raise Exception(
                "No dbSNP file {}".format(dbsnp))

        if dbsnp[-6:] != "vcf.gz":
            logger.error("Aborting!")
            raise Exception(
                "The dbSNP file should be a tabix indexed file with .vcf.gz format")
        if not os.path.exists(dbsnp + ".tbi"):
            logger.error("Aborting!")
            raise Exception(
                "The dbSNP file should be a tabix indexed file with .vcf.gz format. No {}.tbi file exists.".format(dbsnp))

    chrom_order = get_chromosomes_order(reference)
    if cosmic:
        cosmic_vars = {}
        with open(cosmic) as i_f:
            for line in skip_empty(i_f):
                x = line.strip().split("\t")
                chrom, pos, _, ref, alts, _, _, info = x[0:8]
                num_cases = info.split("CNT=")[1].split(
                    ";")[0] if "CNT=" in info else float('nan')
                for alt in alts.split(","):
                    var_id = "-".join([chrom, pos, ref, alt])
                    cosmic_vars[var_id] = num_cases

    exclude_vars = set([])
    if exclude_variants:
        with open(exclude_variants) as i_f:
            for line in skip_empty(i_f):
                if exclude_variants.split(".")[-1] == "tsv" and line[0:5] == "CHROM":
                    continue
                x = line.strip().split("\t")
                chrom, pos, _, ref, alt = x[0:5]
                var_id = "-".join([chrom, pos, ref, alt])
                exclude_vars.add(var_id)

    add_vars = set([])
    if add_variants:
        with open(add_variants) as i_f:
            for line in skip_empty(i_f):
                if add_variants.split(".")[-1] == "tsv" and line[0:5] == "CHROM":
                    continue
                x = line.strip().split("\t")
                chrom, pos, _, ref, alt = x[0:5]
                var_id = "-".join([chrom, pos, ref, alt])
                add_vars.add(var_id)

    all_variants=[]
    with open(candidates_vcf) as i_f:
        for line in skip_empty(i_f):
            chrom, pos, _, ref, alt = line.strip().split("\t")[0:5]
            var_id = "-".join([chrom, pos, ref, alt])
            if exclude_variants:
                if var_id in exclude_vars:
                    continue
            if add_variants:
                if var_id in add_vars:
                    add_vars = add_vars - set([var_id])
            num_cosmic_cases = float('nan')
            if_cosmic = 0
            if cosmic and var_id in cosmic_vars:
                if_cosmic = 1
                num_cosmic_cases = cosmic_vars[var_id]
            all_variants.append([chrom, int(pos), ref, alt, if_cosmic, num_cosmic_cases])

    if add_variants and len(add_vars) > 0:
        for var_id in add_vars - set(exclude_vars):
            v = var_id.split("-")
            pos, ref, alt = v[-3:]
            chrom = "-".join(v[:-3])
            num_cosmic_cases = float('nan')
            if_cosmic = 0
            if cosmic and var_id in cosmic_vars:
                if_cosmic = 1
                num_cosmic_cases = cosmic_vars[var_id]
            all_variants.append([chrom, int(pos), ref, alt, if_cosmic, num_cosmic_cases])

    all_variants = sorted(all_variants,key=lambda x:[chrom_order[x[0]],x[1]])
    n_variants = len(all_variants)
    logger.info("Number of variants: {}".format(n_variants))
    split_len = (n_variants + num_threads - 1) // num_threads
    pool = multiprocessing.Pool(num_threads)
    map_args = []
    nei_cluster = []
    batch = []
    n_batch = 0
    curr_pos = None
    for i, [chrom, pos, ref, alt, if_cosmic, num_cosmic_cases] in enumerate(all_variants):
        if curr_pos is None:
            curr_pos = [chrom, pos]
            nei_cluster = [[chrom, pos, ref, alt, if_cosmic, num_cosmic_cases]]
        else:
            if chrom == curr_pos[0] and abs(curr_pos[1]-pos)<window_extend:
                nei_cluster.append([chrom, pos, ref, alt, if_cosmic, num_cosmic_cases])
            else:
                batch.append(nei_cluster)
                n_batch += len(nei_cluster)
                curr_pos = [chrom, pos]
                nei_cluster = [[chrom, pos, ref, alt, if_cosmic, num_cosmic_cases]]
        if n_batch >= split_len or i == n_variants-1:
            if i == n_variants-1:
                batch.append(nei_cluster)
                curr_pos = None
                nei_cluster = []
            if batch:
                map_args.append((reference, tumor_bam, normal_bam,
                                 min_mapq, min_bq, dbsnp, no_seq_complexity, batch))
            batch = []
    assert(n_variants == sum([len(y) for x in map_args for y in x[-1]]))

    logger.info("Number of batches: {}".format(len(map_args)))
    header = ["CHROM", "POS", "ID", "REF", "ALT", "if_dbsnp", "COMMON", "if_COSMIC", "COSMIC_CNT",
              "Consistent_Mates", "Inconsistent_Mates"]
    if not no_seq_complexity:
        header.extend(["Seq_Complexity_Span", "Seq_Complexity_Adj"])
    header.extend(["N_DP", "nBAM_REF_MQ", "nBAM_ALT_MQ", "nBAM_Z_Ranksums_MQ",
                   "nBAM_REF_BQ", "nBAM_ALT_BQ", "nBAM_Z_Ranksums_BQ", "nBAM_REF_NM", "nBAM_ALT_NM", "nBAM_NM_Diff",
                   "nBAM_REF_Concordant", "nBAM_REF_Discordant", "nBAM_ALT_Concordant", "nBAM_ALT_Discordant",
                   "nBAM_Concordance_FET", "N_REF_FOR", "N_REF_REV", "N_ALT_FOR", "N_ALT_REV", "nBAM_StrandBias_FET",
                   "nBAM_Z_Ranksums_EndPos", "nBAM_REF_Clipped_Reads", "nBAM_ALT_Clipped_Reads", "nBAM_Clipping_FET",
                   "nBAM_MQ0", "nBAM_Other_Reads", "nBAM_Poor_Reads", "nBAM_REF_InDel_3bp", "nBAM_REF_InDel_2bp",
                   "nBAM_REF_InDel_1bp", "nBAM_ALT_InDel_3bp", "nBAM_ALT_InDel_2bp", "nBAM_ALT_InDel_1bp", "SOR",
                   "MaxHomopolymer_Length", "SiteHomopolymer_Length", "T_DP", "tBAM_REF_MQ", "tBAM_ALT_MQ", "tBAM_Z_Ranksums_MQ",
                   "tBAM_REF_BQ", "tBAM_ALT_BQ", "tBAM_Z_Ranksums_BQ", "tBAM_REF_NM", "tBAM_ALT_NM", "tBAM_NM_Diff",
                   "tBAM_REF_Concordant", "tBAM_REF_Discordant", "tBAM_ALT_Concordant", "tBAM_ALT_Discordant",
                   "tBAM_Concordance_FET", "T_REF_FOR", "T_REF_REV", "T_ALT_FOR", "T_ALT_REV", "tBAM_StrandBias_FET",
                   "tBAM_Z_Ranksums_EndPos", "tBAM_REF_Clipped_Reads", "tBAM_ALT_Clipped_Reads", "tBAM_Clipping_FET",
                   "tBAM_MQ0", "tBAM_Other_Reads", "tBAM_Poor_Reads", "tBAM_REF_InDel_3bp", "tBAM_REF_InDel_2bp",
                   "tBAM_REF_InDel_1bp", "tBAM_ALT_InDel_3bp", "tBAM_ALT_InDel_2bp", "tBAM_ALT_InDel_1bp", "InDel_Length"])

    try:
        ext_features = pool.map_async(extract_features, map_args).get()
        pool.close()
        with open(output_tsv, "w") as o_f:
            o_f.write("\t".join(header) + "\n")
            for features in ext_features:
                for w in features:
                    o_f.write(
                        "\t".join(map(lambda x: str(x).replace("nan", "0"), w)) + "\n")
    except Exception as inst:
        logger.error(inst)
        pool.close()
        traceback.print_exc()
        raise Exception

    logger.info("Done Extending Standalone Features.")
    return ext_features


if __name__ == '__main__':
    FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    logging.basicConfig(level=logging.INFO, format=FORMAT)
    logger = logging.getLogger(__name__)

    parser = argparse.ArgumentParser(
        description='extract extra features for standalone mode')
    parser.add_argument('--candidates_vcf', type=str, help='candidates vcf',
                        required=True)
    parser.add_argument('--exclude_variants', type=str, help='variants to exclude',
                        default=None)
    parser.add_argument('--add_variants', type=str, help='variants to add if not exist in vcf. (Lower priority than --exclude_variants)',
                        default=None)
    parser.add_argument('--output_tsv', type=str, help='output features tsv',
                        required=True)
    parser.add_argument('--reference', type=str, help='reference fasta filename',
                        required=True)
    parser.add_argument('--tumor_bam', type=str,
                        help='tumor bam', required=True)
    parser.add_argument('--normal_bam', type=str,
                        help='normal bam', required=True)
    parser.add_argument('--min_mapq', type=int,
                        help='minimum mapping quality', default=1)
    parser.add_argument('--min_bq', type=float,
                        help='minimum base quality', default=5)
    parser.add_argument('--dbsnp', type=str,
                        help='dbSNP vcf (to annotate candidate variants)', default=None)
    parser.add_argument('--cosmic', type=str,
                        help='COSMIC vcf (to annotate candidate variants)', default=None)
    parser.add_argument('--no_seq_complexity',
                        help='Dont compute linguistic sequence complexity features',
                        action="store_true")
    parser.add_argument('--window_extend', type=int,
                        help='window size for extending input features (should be in the order of readlength)',
                         default=1000)
    parser.add_argument('--num_threads', type=int,
                        help='number of threads', default=1)
    args = parser.parse_args()
    logger.info(args)

    try:
        output = extend_features(args.candidates_vcf,
                                 args.exclude_variants,
                                 args.add_variants,
                                 args.output_tsv,
                                 args.reference, args.tumor_bam, args.normal_bam,
                                 args.min_mapq, args.min_bq,
                                 args.dbsnp, args.cosmic,
                                 args.no_seq_complexity,
                                 args.window_extend,
                                 args.num_threads,
                                 )
        if output is None:
            raise Exception("extend_features failed!")
    except Exception as e:
        logger.error(traceback.format_exc())
        logger.error("Aborting!")
        logger.error(
            "extend_features.py failure on arguments: {}".format(args))
        raise e
