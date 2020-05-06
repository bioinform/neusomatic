#!/usr/bin/env python
#-------------------------------------------------------------------------
# preprocess.py
# A wrapper that
# 1- scans the tumor and normal alignments to extract features and raw candidates (call to 'scan_alignments.py')
# 2- filters candidates to met cut-offs (call to 'filter_candidates.py')
# 3- generate datasets for candidates to be used by network (call to 'generate_dataset.py')
#-------------------------------------------------------------------------

import multiprocessing
import argparse
import os
import shutil
import traceback
import logging

import tempfile
import numpy as np

from filter_candidates import filter_candidates
from generate_dataset import generate_dataset, extract_ensemble
from scan_alignments import scan_alignments
from extend_features import extend_features
from utils import concatenate_vcfs, run_bedtools_cmd, bedtools_sort, bedtools_merge, bedtools_intersect, bedtools_slop, get_tmp_file, skip_empty, vcf_2_bed


def process_split_region(tn, work, region, reference, mode, alignment_bam, dbsnp,
                         scan_window_size, scan_maf, min_mapq,
                         filtered_candidates_vcf, min_dp, max_dp,
                         filter_duplicate,
                         good_ao, min_ao, snp_min_af, snp_min_bq, snp_min_ao,
                         ins_min_af, del_min_af, del_merge_min_af,
                         ins_merge_min_af, merge_r,
                         scan_alignments_binary, restart, num_threads, calc_qual, regions=[]):

    logger = logging.getLogger(process_split_region.__name__)
    logger.info("Scan bam.")
    scan_outputs = scan_alignments(work, scan_alignments_binary, alignment_bam,
                                   region, reference, num_threads, scan_window_size, scan_maf,
                                   min_mapq, max_dp, filter_duplicate, restart=restart, split_region_files=regions,
                                   calc_qual=calc_qual)
    if filtered_candidates_vcf:
        logger.info("Filter candidates.")
        if restart or not os.path.exists(filtered_candidates_vcf):
            pool = multiprocessing.Pool(num_threads)
            map_args = []
            for i, (raw_vcf, count_bed, split_region_bed) in enumerate(scan_outputs):
                filtered_vcf = os.path.join(os.path.dirname(
                    os.path.realpath(raw_vcf)), "filtered_candidates.vcf")
                map_args.append((raw_vcf, filtered_vcf, reference, dbsnp, min_dp, max_dp, good_ao,
                                 min_ao, snp_min_af, snp_min_bq, snp_min_ao, ins_min_af, del_min_af, del_merge_min_af,
                                 ins_merge_min_af, merge_r))
            try:
                filtered_candidates_vcfs = pool.map_async(
                    filter_candidates, map_args).get()
                pool.close()
            except Exception as inst:
                logger.error(inst)
                pool.close()
                traceback.print_exc()
                raise Exception

            for o in filtered_candidates_vcfs:
                if o is None:
                    raise Exception("filter_candidates failed!")

            concatenate_vcfs(filtered_candidates_vcfs,
                             filtered_candidates_vcf, check_file_existence=True)

        else:
            filtered_candidates_vcfs = []
            for raw_vcf, _, _ in scan_outputs:
                filtered_vcf = os.path.join(os.path.dirname(
                    os.path.realpath(raw_vcf)), "filtered_candidates.vcf")
                filtered_candidates_vcfs.append(filtered_vcf)
    else:
        filtered_candidates_vcfs = None
    return list(map(lambda x: x[1], scan_outputs)), list(map(lambda x: x[2], scan_outputs)), filtered_candidates_vcfs


def generate_dataset_region(work, truth_vcf, mode, filtered_candidates_vcf, region, tumor_count_bed, normal_count_bed, reference,
                            matrix_width, matrix_base_pad, min_ev_frac_per_col, min_cov, num_threads, ensemble_bed, seq_complexity, 
                            no_feature_recomp_for_ensemble, tsv_batch_size):
    logger = logging.getLogger(generate_dataset_region.__name__)
    generate_dataset(work, truth_vcf, mode, filtered_candidates_vcf, region, tumor_count_bed, normal_count_bed, reference,
                     matrix_width, matrix_base_pad, min_ev_frac_per_col, min_cov, num_threads, None, ensemble_bed,
                     seq_complexity,
                     no_feature_recomp_for_ensemble,
                     tsv_batch_size)


def get_ensemble_region(record):
    reference, ensemble_bed, region, ensemble_bed_region_file, matrix_base_pad = record
    thread_logger = logging.getLogger(
        "{} ({})".format(get_ensemble_region.__name__, multiprocessing.current_process().name))
    try:
        ensemble_bed_region_file_tmp = bedtools_slop(
            region, reference + ".fai", args=" -b {}".format(matrix_base_pad + 3),
            run_logger=thread_logger)
        bedtools_intersect(
            ensemble_bed, ensemble_bed_region_file_tmp, args=" -u -header",
            output_fn=ensemble_bed_region_file, run_logger=thread_logger)
        return ensemble_bed_region_file

    except Exception as ex:
        thread_logger.error(traceback.format_exc())
        thread_logger.error(ex)
        return None


def get_ensemble_beds(work, reference, ensemble_bed, split_regions, matrix_base_pad, num_threads):
    logger = logging.getLogger(get_ensemble_beds.__name__)

    work_ensemble = os.path.join(work, "ensemble_anns")
    if not os.path.exists(work_ensemble):
        os.mkdir(work_ensemble)
    map_args = []
    ensemble_beds = []
    for i, split_region_ in enumerate(split_regions):
        ensemble_bed_region_file = os.path.join(
            work_ensemble, "ensemble_ann_{}.bed".format(i))
        ensemble_beds.append(ensemble_bed_region_file)
        map_args.append((reference, ensemble_bed, split_region_,
                         ensemble_bed_region_file, matrix_base_pad))
    pool = multiprocessing.Pool(num_threads)
    try:
        outputs = pool.map_async(get_ensemble_region, map_args).get()
        pool.close()
    except Exception as inst:
        logger.error(inst)
        pool.close()
        traceback.print_exc()
        raise Exception
    for o in outputs:
        if o is None:
            raise Exception("get_ensemble_region failed!")
    return ensemble_beds


def extract_candidate_split_regions(
        work, filtered_candidates_vcfs, split_regions, ensemble_beds,
        reference, matrix_base_pad, merge_d_for_short_read):
    logger = logging.getLogger(extract_candidate_split_regions.__name__)

    candidates_split_regions = []
    for i, (filtered_vcf, split_region_) in enumerate(zip(filtered_candidates_vcfs,
                                                          split_regions)):
        candidates_region_file = os.path.join(
            work, "candidates_region_{}.bed".format(i))

        is_empty = True
        with open(filtered_vcf) as f_:
            for line in skip_empty(f_):
                is_empty = False
                break
        logger.info([filtered_vcf, is_empty])
        if not is_empty:
            candidates_bed = get_tmp_file()
            vcf_2_bed(filtered_vcf, candidates_bed,
                      len_ref=True, keep_ref_alt=False)

            candidates_bed = bedtools_sort(candidates_bed, run_logger=logger)
            candidates_bed = bedtools_slop(
                candidates_bed, reference + ".fai", args=" -b {}".format(matrix_base_pad + 3),
                run_logger=logger)
            candidates_bed = bedtools_merge(
                candidates_bed, args=" -d {}".format(merge_d_for_short_read), run_logger=logger)
        else:
            candidates_bed = get_tmp_file()
        if ensemble_beds:
            cmd = "cat {} {}".format(
                candidates_bed,
                ensemble_beds[i])
            candidates_bed = run_bedtools_cmd(cmd, run_logger=logger)
            cmd = "cut -f 1,2,3 {}".format(
                candidates_bed)
            candidates_bed = run_bedtools_cmd(cmd, run_logger=logger)
            candidates_bed = bedtools_sort(candidates_bed, run_logger=logger)
            candidates_bed = bedtools_merge(
                candidates_bed, args=" -d {}".format(merge_d_for_short_read), run_logger=logger)

        candidates_bed = bedtools_intersect(
            candidates_bed, split_region_, run_logger=logger)
        bedtools_sort(candidates_bed,
                      output_fn=candidates_region_file, run_logger=logger)

        candidates_split_regions.append(candidates_region_file)
    return candidates_split_regions


def preprocess(work, mode, reference, region_bed, tumor_bam, normal_bam, dbsnp,
               scan_window_size, scan_maf, min_mapq,
               min_dp, max_dp, good_ao, min_ao, snp_min_af, snp_min_bq, snp_min_ao,
               ins_min_af, del_min_af, del_merge_min_af,
               ins_merge_min_af, merge_r, truth_vcf, tsv_batch_size,
               matrix_width, matrix_base_pad, min_ev_frac_per_col,
               ensemble_tsv, long_read, restart, first_do_without_qual,
               filter_duplicate,
               add_extra_features,
               seq_complexity,
               no_feature_recomp_for_ensemble,
               num_threads,
               scan_alignments_binary,):
    logger = logging.getLogger(preprocess.__name__)

    logger.info("----------------------Preprocessing------------------------")
    if restart or not os.path.exists(work):
        os.mkdir(work)

    original_tempdir = tempfile.tempdir
    bed_tempdir = os.path.join(work, "bed_tempdir_preprocess")
    if not os.path.exists(bed_tempdir):
        os.mkdir(bed_tempdir)
    tempfile.tempdir = bed_tempdir

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
        if dbsnp[-6:] != "vcf.gz":
            logger.error("Aborting!")
            raise Exception(
                "The dbSNP file should be a tabix indexed file with .vcf.gz format")
        if not os.path.exists(dbsnp + ".tbi"):
            logger.error("Aborting!")
            raise Exception(
                "The dbSNP file should be a tabix indexed file with .vcf.gz format. No {}.tbi file exists.".format(dbsnp))

    ensemble_bed = None
    if ensemble_tsv:
        ensemble_bed = os.path.join(work, "ensemble.bed")
        logger.info("Extract ensemble info.")
        if restart or not os.path.exists(ensemble_bed):
            extract_ensemble(ensemble_tsv, ensemble_bed, seq_complexity, no_feature_recomp_for_ensemble, False)

    merge_d_for_short_read = 100
    candidates_split_regions = []
    ensemble_beds = []
    if not long_read and first_do_without_qual:
        logger.info("Scan tumor bam (first without quality scores).")
        work_tumor_without_q = os.path.join(work, "work_tumor_without_q")
        if restart or not os.path.exists(work_tumor_without_q):
            os.mkdir(work_tumor_without_q)
        filtered_candidates_vcf_without_q = os.path.join(
            work_tumor_without_q, "filtered_candidates.vcf")

        tumor_outputs_without_q = process_split_region("tumor", work_tumor_without_q, region_bed, reference, mode,
                                                       tumor_bam, dbsnp, scan_window_size, scan_maf, min_mapq,
                                                       filtered_candidates_vcf_without_q, min_dp, max_dp,
                                                       filter_duplicate,
                                                       good_ao, min_ao,
                                                       snp_min_af, -10000, snp_min_ao,
                                                       ins_min_af, del_min_af, del_merge_min_af,
                                                       ins_merge_min_af, merge_r,
                                                       scan_alignments_binary, restart, num_threads,
                                                       calc_qual=False)
        tumor_counts_without_q, split_regions, filtered_candidates_vcfs_without_q = tumor_outputs_without_q

        if ensemble_tsv:
            ensemble_beds = get_ensemble_beds(
                work, reference, ensemble_bed, split_regions, matrix_base_pad, num_threads)
        candidates_split_regions = extract_candidate_split_regions(
            work_tumor_without_q, filtered_candidates_vcfs_without_q, split_regions, ensemble_beds,
            reference, matrix_base_pad, merge_d_for_short_read)
    work_tumor = os.path.join(work, "work_tumor")
    if restart or not os.path.exists(work_tumor):
        os.mkdir(work_tumor)
    filtered_candidates_vcf = os.path.join(
        work_tumor, "filtered_candidates.vcf")

    logger.info("Scan tumor bam (and extracting quality scores).")
    tumor_outputs = process_split_region("tumor", work_tumor, region_bed, reference, mode,
                                         tumor_bam, dbsnp, scan_window_size, scan_maf, min_mapq,
                                         filtered_candidates_vcf, min_dp, max_dp,
                                         filter_duplicate,
                                         good_ao, min_ao,
                                         snp_min_af, snp_min_bq, snp_min_ao,
                                         ins_min_af, del_min_af, del_merge_min_af,
                                         ins_merge_min_af, merge_r,
                                         scan_alignments_binary, restart, num_threads,
                                         calc_qual=True,
                                         regions=candidates_split_regions)
    tumor_counts, split_regions, filtered_candidates_vcfs = tumor_outputs

    if ensemble_tsv and not ensemble_beds:
        ensemble_beds = get_ensemble_beds(
            work, reference, ensemble_bed, split_regions, matrix_base_pad, num_threads)

    if (not long_read):
        candidates_split_regions = extract_candidate_split_regions(
            work_tumor, filtered_candidates_vcfs, split_regions, ensemble_beds,
            reference, matrix_base_pad, merge_d_for_short_read)

    if not candidates_split_regions:
        candidates_split_regions = split_regions
    work_normal = os.path.join(work, "work_normal")
    if restart or not os.path.exists(work_normal):
        os.mkdir(work_normal)
    logger.info("Scan normal bam (and extracting quality scores).")
    normal_counts, _, _ = process_split_region("normal", work_normal, region_bed, reference, mode, normal_bam,
                                               None, scan_window_size, 0.2, min_mapq,
                                               None, min_dp, max_dp,
                                               filter_duplicate,
                                               good_ao, min_ao, snp_min_af, snp_min_bq, snp_min_ao,
                                               ins_min_af, del_min_af, del_merge_min_af,
                                               ins_merge_min_af, merge_r,
                                               scan_alignments_binary, restart, num_threads,
                                               calc_qual=True,
                                               regions=candidates_split_regions)

    work_dataset = os.path.join(work, "dataset")
    if restart or not os.path.exists(work_dataset):
        os.mkdir(work_dataset)
    logger.info("Generate dataset.")
    for i, (tumor_count, normal_count, filtered_vcf, candidates_split_region) in enumerate(zip(tumor_counts, normal_counts, filtered_candidates_vcfs, candidates_split_regions)):
        logger.info("Dataset for region {}".format(candidates_split_region))
        work_dataset_split = os.path.join(work_dataset, "work.{}".format(i))
        if restart or not os.path.exists("{}/done.txt".format(work_dataset_split)):
            if os.path.exists(work_dataset_split):
                shutil.rmtree(work_dataset_split)
            os.mkdir(work_dataset_split)
            ensemble_bed_i = ensemble_beds[i] if ensemble_tsv else None
            if add_extra_features:
                work_tumor_i = os.path.dirname(filtered_vcf)
                extra_features_tsv = os.path.join(
                    work_tumor_i, "extra_features.tsv")
                if not os.path.exists(extra_features_tsv) or restart:
                    extend_features(filtered_vcf,
                                    ensemble_beds[
                                        i] if (ensemble_tsv and no_feature_recomp_for_ensemble) else None,
                                    ensemble_beds[
                                        i] if (ensemble_tsv and not no_feature_recomp_for_ensemble) else None,
                                    extra_features_tsv,
                                    reference, tumor_bam, normal_bam,
                                    min_mapq, snp_min_bq,
                                    dbsnp, None,
                                    seq_complexity,
                                    num_threads)
                extra_features_bed = os.path.join(
                    work_dataset_split, "extra_features.bed")
                if not os.path.exists(extra_features_bed) or restart:
                    extract_ensemble(extra_features_tsv,
                                     extra_features_bed, seq_complexity, True, True)
                if ensemble_tsv:
                    merged_features_bed = os.path.join(
                        work_dataset_split, "merged_features.bed")
                    if not os.path.exists(merged_features_bed) or restart:
                        exclude_ens_variants = []
                        header_line = ""
                        if no_feature_recomp_for_ensemble:
                            with open(merged_features_bed, "w") as o_f, open(ensemble_beds[i]) as i_f_1, open(extra_features_bed) as i_f_2:
                                for line in skip_empty(i_f_1, skip_header=False):
                                    if line.startswith("#"):
                                        if not header_line:
                                            header_line = line
                                            o_f.write(line)
                                        else:
                                            if header_line != line:
                                                logger.error(
                                                    "{}!={}".format(header_line, line))
                                                raise Exception
                                        continue
                                    chrom, pos, _, ref, alt = line.strip().split("\t")[
                                        0:5]
                                    var_id = "-".join([chrom, pos, ref, alt])
                                    exclude_ens_variants.append(var_id)
                                    o_f.write(line)
                                for line in skip_empty(i_f_2, skip_header=False):
                                    if line.startswith("#"):
                                        if header_line != line:
                                            logger.error(
                                                "{}!={}".format(header_line, line))
                                            raise Exception
                                        continue
                                    chrom, pos, _, ref, alt = line.strip().split("\t")[
                                        0:5]
                                    var_id = "-".join([chrom, pos, ref, alt])
                                    if var_id in exclude_ens_variants:
                                        continue
                                    o_f.write(line)
                        else:
                            callers_features = ["if_MuTect", "if_VarScan2", "if_JointSNVMix2", "if_SomaticSniper", "if_VarDict", "MuSE_Tier",
                                                "if_LoFreq", "if_Scalpel", "if_Strelka", "if_TNscope", "Strelka_Score", "Strelka_QSS",
                                                "Strelka_TQSS", "VarScan2_Score", "SNVMix2_Score", "Sniper_Score", "VarDict_Score",
                                                "M2_NLOD", "M2_TLOD", "M2_STR", "M2_ECNT", "MSI", "MSILEN", "SHIFT3"]
                            with open(merged_features_bed, "w") as o_f, open(ensemble_beds[i]) as i_f_1, open(extra_features_bed) as i_f_2:
                                ens_variants_info = {}
                                header_1_found = False
                                header_2_found = False
                                for line in skip_empty(i_f_1, skip_header=False):
                                    if line.startswith("#"):
                                        if not header_line:
                                            header_line = line
                                        else:
                                            if header_line != line:
                                                logger.error(
                                                    "{}!={}".format(header_line, line))
                                                raise Exception
                                        header_ = line.strip().split()[5:]
                                        header_caller = list(filter(
                                            lambda x: x[1] in callers_features, enumerate(header_)))
                                        header_caller_ = list(
                                            map(lambda x: x[1], header_caller))
                                        header_i = list(
                                            map(lambda x: x[0], header_caller))
                                        header_1_found = True
                                        continue
                                    assert header_1_found
                                    fields = line.strip().split("\t")
                                    chrom, pos, _, ref, alt = fields[0:5]
                                    var_id = "-".join([chrom, pos, ref, alt])
                                    ens_variants_info[var_id] = np.array(fields[5:])[
                                        header_i]
                                for line in skip_empty(i_f_2, skip_header=False):
                                    if line.startswith("#"):
                                        if header_line != line:
                                            logger.error(
                                                "{}!={}".format(header_line, line))
                                        if not header_2_found:
                                            header_2 = line.strip().split()[5:]
                                            logger.info(header_2)
                                            order_header = []
                                            for f in header_caller_:
                                                if f not in header_2:
                                                    logger.info("Missing header field {}".format(f))
                                                    raise Exception
                                                order_header.append(header_2.index(f))
                                            o_f.write(line)
                                        header_2_found = True

                                    assert header_2_found
                                    fields = line.strip().split("\t")
                                    chrom, pos, _, ref, alt = fields[0:5]
                                    var_id = "-".join([chrom, pos, ref, alt])
                                    if var_id in ens_variants_info:
                                        fields_ = np.array(fields[5:])
                                        fields_[order_header] = ens_variants_info[
                                            var_id]
                                        fields[5:] = fields_.tolist()
                                    o_f.write(
                                        "\t".join(list(map(str, fields))) + "\n")
                    ensemble_bed_i = merged_features_bed
                else:
                    ensemble_bed_i = extra_features_bed

            generate_dataset_region(work_dataset_split, truth_vcf, mode, filtered_vcf,
                                    candidates_split_region, tumor_count, normal_count, reference,
                                    matrix_width, matrix_base_pad, min_ev_frac_per_col, min_dp, num_threads,
                                    ensemble_bed_i, seq_complexity, no_feature_recomp_for_ensemble, tsv_batch_size)

    shutil.rmtree(bed_tempdir)
    tempfile.tempdir = original_tempdir

    logger.info("Preprocessing is Done.")


if __name__ == '__main__':
    FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    logging.basicConfig(level=logging.INFO, format=FORMAT)
    logger = logging.getLogger(__name__)

    parser = argparse.ArgumentParser(
        description='Preprocess input alignments for train/call')
    parser.add_argument('--mode', type=str, help='train/call mode',
                        choices=["train", "call"], required=True)
    parser.add_argument('--reference', type=str,
                        help='reference fasta filename', required=True)
    parser.add_argument('--region_bed', type=str,
                        help='region bed', required=True)
    parser.add_argument('--tumor_bam', type=str,
                        help='tumor bam', required=True)
    parser.add_argument('--normal_bam', type=str,
                        help='normal bam', required=True)
    parser.add_argument('--work', type=str,
                        help='work directory', required=True)
    parser.add_argument('--dbsnp_to_filter', type=str,
                        help='dbsnp vcf.gz (will be used to filter candidate variants)', default=None)
    parser.add_argument('--scan_window_size', type=int,
                        help='window size to scan the variants', default=2000)
    parser.add_argument('--scan_maf', type=float,
                        help='minimum allele freq for scanning', default=0.01)
    parser.add_argument('--min_mapq', type=int,
                        help='minimum mapping quality', default=1)
    parser.add_argument('--min_dp', type=float, help='min depth', default=5)
    parser.add_argument('--max_dp', type=float,
                        help='max depth', default=100000)
    parser.add_argument('--good_ao', type=float,
                        help='good alternate count (ignores maf)', default=10)
    parser.add_argument('--min_ao', type=float,
                        help='min alternate count', default=1)
    parser.add_argument('--snp_min_af', type=float,
                        help='SNP min allele freq', default=0.05)
    parser.add_argument('--snp_min_bq', type=float,
                        help='SNP min base quality', default=10)
    parser.add_argument('--snp_min_ao', type=float,
                        help='SNP min alternate count for low AF candidates', default=3)
    parser.add_argument('--ins_min_af', type=float,
                        help='INS min allele freq', default=0.05)
    parser.add_argument('--del_min_af', type=float,
                        help='DEL min allele freq', default=0.05)
    parser.add_argument('--del_merge_min_af', type=float,
                        help='min allele freq for merging DELs', default=0)
    parser.add_argument('--ins_merge_min_af', type=float,
                        help='min allele freq for merging INSs', default=0)
    parser.add_argument('--merge_r', type=float,
                        help='merge af ratio to the max af for merging adjacent variants', default=0.5)
    parser.add_argument('--truth_vcf', type=str,
                        help='truth vcf (required for train mode)', default=None)
    parser.add_argument('--tsv_batch_size', type=int,
                        help='output files batch size', default=50000)
    parser.add_argument('--matrix_width', type=int,
                        help='target window width', default=32)
    parser.add_argument('--matrix_base_pad', type=int,
                        help='number of bases to pad around the candidate variant', default=7)
    parser.add_argument('--min_ev_frac_per_col', type=float,
                        help='minimum frac cov per column to keep columm', default=0.06)
    parser.add_argument('--ensemble_tsv', type=str,
                        help='Ensemble annotation tsv file (only for short read)', default=None)
    parser.add_argument('--long_read',
                        help='Enable long_read (high error-rate sequence) indel realignment',
                        action="store_true")
    parser.add_argument('--restart',
                        help='Restart the process. (instead of continuing from where we left)',
                        action="store_true")
    parser.add_argument('--first_do_without_qual',
                        help='Perform initial scan without calculating the quality stats',
                        action="store_true")
    parser.add_argument('--filter_duplicate',
                        help='filter duplicate reads when preparing pileup information',
                        action="store_true")
    parser.add_argument('--add_extra_features',
                        help='add extra input features',
                        action="store_true")
    parser.add_argument('--seq_complexity',
                        help='Compute linguistic sequence complexity features',
                        action="store_true")
    parser.add_argument('--no_feature_recomp_for_ensemble',
                        help='Do not recompute features for ensemble_tsv',
                        action="store_true")
    parser.add_argument('--num_threads', type=int,
                        help='number of threads', default=1)
    parser.add_argument('--scan_alignments_binary', type=str,
                        help='binary for scanning alignment bam', default="../bin/scan_alignments")
    args = parser.parse_args()
    logger.info(args)

    try:
        preprocess(args.work, args.mode, args.reference, args.region_bed, args.tumor_bam, args.normal_bam,
                   args.dbsnp_to_filter,
                   args.scan_window_size, args.scan_maf, args.min_mapq,
                   args.min_dp, args.max_dp, args.good_ao, args.min_ao, args.snp_min_af, args.snp_min_bq, args.snp_min_ao,
                   args.ins_min_af, args.del_min_af, args.del_merge_min_af,
                   args.ins_merge_min_af, args.merge_r,
                   args.truth_vcf, args.tsv_batch_size, args.matrix_width, args.matrix_base_pad, args.min_ev_frac_per_col,
                   args.ensemble_tsv, args.long_read, args.restart, args.first_do_without_qual,
                   args.filter_duplicate,
                   args.add_extra_features,
                   args.seq_complexity,
                   args.no_feature_recomp_for_ensemble,
                   args.num_threads,
                   args.scan_alignments_binary)
    except Exception as e:
        logger.error(traceback.format_exc())
        logger.error("Aborting!")
        logger.error(
            "preprocess.py failure on arguments: {}".format(args))
        raise e
