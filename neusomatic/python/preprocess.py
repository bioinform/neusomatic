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

import pybedtools

from filter_candidates import filter_candidates
from generate_dataset import generate_dataset, extract_ensemble
from scan_alignments import scan_alignments
from utils import concatenate_vcfs


def split_dbsnp(record):
    restart, dbsnp, region_bed, dbsnp_region_vcf = record
    thread_logger = logging.getLogger(
        "{} ({})".format(split_dbsnp.__name__, multiprocessing.current_process().name))
    try:
        if restart or not os.path.exists(dbsnp_region_vcf):
            pybedtools.BedTool(dbsnp).intersect(
                region_bed, u=True).saveas(dbsnp_region_vcf)
        return dbsnp_region_vcf
    except Exception as ex:
        thread_logger.error(traceback.format_exc())
        thread_logger.error(ex)
        return None


def process_split_region(tn, work, region, reference, mode, alignment_bam, dbsnp,
                         scan_window_size, scan_maf, min_mapq,
                         filtered_candidates_vcf, min_dp, max_dp, 
                         filter_duplicate,
                         good_ao, min_ao, snp_min_af, snp_min_bq, snp_min_ao,
                         ins_min_af, del_min_af, del_merge_min_af,
                         ins_merge_min_af, merge_r,
                         scan_alignments_binary, restart, num_threads, calc_qual, regions=[], dbsnp_regions=[]):

    logger = logging.getLogger(process_split_region.__name__)
    logger.info("Scan bam.")
    scan_outputs = scan_alignments(work, scan_alignments_binary, alignment_bam,
                                   region, reference, num_threads, scan_window_size, scan_maf,
                                   min_mapq, max_dp, filter_duplicate, restart=restart, split_region_files=regions,
                                   calc_qual=calc_qual)
    if filtered_candidates_vcf:
        logger.info("Filter candidates.")
        if restart or not os.path.exists(filtered_candidates_vcf):
            if dbsnp and not dbsnp_regions:
                map_args = []
                for raw_vcf, count_bed, split_region_bed in scan_outputs:
                    dbsnp_region_vcf = os.path.join(os.path.dirname(
                        os.path.realpath(raw_vcf)), "dbsnp_region.vcf")
                    map_args.append(
                        (restart, dbsnp, split_region_bed, dbsnp_region_vcf))
                pool = multiprocessing.Pool(num_threads)
                try:
                    dbsnp_regions = pool.map_async(split_dbsnp, map_args).get()
                    pool.close()
                except Exception as inst:
                    logger.error(inst)
                    pool.close()
                    traceback.print_exc()
                    raise Exception

                for o in dbsnp_regions:
                    if o is None:
                        raise Exception("split_dbsnp failed!")

            pool = multiprocessing.Pool(num_threads)
            map_args = []
            for i, (raw_vcf, count_bed, split_region_bed) in enumerate(scan_outputs):
                filtered_vcf = os.path.join(os.path.dirname(
                    os.path.realpath(raw_vcf)), "filtered_candidates.vcf")
                dbsnp_region_vcf = None
                if dbsnp:
                    dbsnp_region_vcf = dbsnp_regions[i]
                map_args.append((raw_vcf, filtered_vcf, reference, dbsnp_region_vcf, min_dp, max_dp, good_ao,
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
            if dbsnp and not dbsnp_regions:
                dbsnp_regions = []
                for raw_vcf, _, _ in scan_outputs:
                    dbsnp_region_vcf = os.path.join(os.path.dirname(
                        os.path.realpath(raw_vcf)), "dbsnp_region.vcf")
                    dbsnp_regions.append(dbsnp_region_vcf)
    else:
        filtered_candidates_vcfs = None
    return list(map(lambda x: x[1], scan_outputs)), list(map(lambda x: x[2], scan_outputs)), filtered_candidates_vcfs, dbsnp_regions


def generate_dataset_region(work, truth_vcf, mode, filtered_candidates_vcf, region, tumor_count_bed, normal_count_bed, reference,
                            matrix_width, matrix_base_pad, min_ev_frac_per_col, min_cov, num_threads, ensemble_bed, tsv_batch_size):
    logger = logging.getLogger(generate_dataset_region.__name__)
    generate_dataset(work, truth_vcf, mode, filtered_candidates_vcf, region, tumor_count_bed, normal_count_bed, reference,
                     matrix_width, matrix_base_pad, min_ev_frac_per_col, min_cov, num_threads, None, ensemble_bed,
                     tsv_batch_size)


def get_ensemble_region(record):
    reference, ensemble_bed, region, ensemble_bed_region_file, matrix_base_pad = record
    thread_logger = logging.getLogger(
        "{} ({})".format(get_ensemble_region.__name__, multiprocessing.current_process().name))
    try:
        pybedtools.BedTool(ensemble_bed).intersect(
            pybedtools.BedTool(region).slop(
                g=reference + ".fai", b=matrix_base_pad + 3), u=True).saveas(ensemble_bed_region_file)
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
        candidates_bed = pybedtools.BedTool(filtered_vcf).each(
            lambda x: pybedtools.Interval(x[0], int(x[1]), int(x[1]) + len(x[3]))).sort().slop(
            g=reference + ".fai", b=matrix_base_pad + 3).merge(d=merge_d_for_short_read)
        if ensemble_beds:
            candidates_bed = candidates_bed.cat(ensemble_beds[i], postmerge=False).sort(
            ).merge(d=merge_d_for_short_read)
        candidates_bed.intersect(split_region_).sort().saveas(
            candidates_region_file)
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
               num_threads,
               scan_alignments_binary,):
    logger = logging.getLogger(preprocess.__name__)

    logger.info("----------------------Preprocessing------------------------")
    if restart or not os.path.exists(work):
        os.mkdir(work)

    original_tempdir = pybedtools.get_tempdir()
    pybedtmp = os.path.join(work, "pybedtmp_preprocess")
    if not os.path.exists(pybedtmp):
        os.mkdir(pybedtmp)
    pybedtools.set_tempdir(pybedtmp)

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

    ensemble_bed = None
    if ensemble_tsv:
        ensemble_bed = os.path.join(work, "ensemble.bed")
        logger.info("Extract ensemble info.")
        if restart or not os.path.exists(ensemble_bed):
            ensemble_bed = extract_ensemble(work, ensemble_tsv)

    merge_d_for_short_read = 100
    candidates_split_regions = []
    dbsnp_regions_q = []
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
                                                       calc_qual=False, dbsnp_regions=[])
        tumor_counts_without_q, split_regions, filtered_candidates_vcfs_without_q, dbsnp_regions_q = tumor_outputs_without_q

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
                                         regions=candidates_split_regions,
                                         dbsnp_regions=dbsnp_regions_q)
    tumor_counts, split_regions, filtered_candidates_vcfs, _ = tumor_outputs

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
    normal_counts, _, _, _ = process_split_region("normal", work_normal, region_bed, reference, mode, normal_bam,
                                                  None, scan_window_size, 0.2, min_mapq,
                                                  None, min_dp, max_dp, 
                                                  filter_duplicate,
                                                  good_ao, min_ao, snp_min_af, snp_min_bq, snp_min_ao,
                                                  ins_min_af, del_min_af, del_merge_min_af,
                                                  ins_merge_min_af, merge_r,
                                                  scan_alignments_binary, restart, num_threads,
                                                  calc_qual=True,
                                                  regions=candidates_split_regions,
                                                  dbsnp_regions=[])

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
            generate_dataset_region(work_dataset_split, truth_vcf, mode, filtered_vcf,
                                    candidates_split_region, tumor_count, normal_count, reference,
                                    matrix_width, matrix_base_pad, min_ev_frac_per_col, min_dp, num_threads,
                                    ensemble_beds[i] if ensemble_tsv else None, tsv_batch_size)

    shutil.rmtree(pybedtmp)
    pybedtools.set_tempdir(original_tempdir)

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
                        help='dbsnp vcf (will be used to filter candidate variants)', default=None)
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
                   args.num_threads,
                   args.scan_alignments_binary)
    except Exception as e:
        logger.error(traceback.format_exc())
        logger.error("Aborting!")
        logger.error(
            "preprocess.py failure on arguments: {}".format(args))
        raise e
