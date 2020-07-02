#!/usr/bin/env python
#-------------------------------------------------------------------------
# scan_alignments.py
# Scan the alignment .bam file, extract A/C/G/T/- counts on augmented alignment,
# as well as different alignment feature matrices such as base quality, mapping
# quality, strandness, clipping, alignment score, ...
# It also outputs a raw .vcf files of potential candidates. This .vcf file should be processed
# by 'filter_candidates.py' before it can be used.
#-------------------------------------------------------------------------

import os
import multiprocessing
import argparse
import glob
import traceback
import logging
import shutil

import pysam
import numpy as np

from utils import concatenate_files, run_shell_command, bedtools_sort, bedtools_merge, get_tmp_file, skip_empty
from split_bed import split_region


def run_scan_alignments(record):
    work, reference, merge_d_for_scan, scan_alignments_binary, split_region_file, \
        input_bam, window_size, maf, min_mapq, max_dp, filter_duplicate, calc_qual = record

    if filter_duplicate:
        filter_duplicate_str = "--filter_duplicate"
    else:
        filter_duplicate_str = ""
    thread_logger = logging.getLogger(
        "{} ({})".format(run_scan_alignments.__name__, multiprocessing.current_process().name))
    try:

        if not os.path.exists(scan_alignments_binary):
            raise IOError("File not found: {}".format(scan_alignments_binary))
        if not os.path.exists(work):
            os.mkdir(work)

        if merge_d_for_scan is not None:
            split_region_file_ = os.path.join(work, "merged_region.bed")
            tmp_ = bedtools_sort(split_region_file, run_logger=thread_logger)
            bedtools_merge(
                tmp_, output_fn=split_region_file_, args=" -d {}".format(merge_d_for_scan), run_logger=thread_logger)
        else:
            split_region_file_ = split_region_file

        if os.path.getsize(split_region_file_) > 0:
            cmd = "{} --ref {} -b {} -L {} --out_vcf_file {}/candidates.vcf --out_count_file {}/count.bed \
                        --window_size {} --min_af {} --min_mapq {} --max_depth {} {}".format(
                scan_alignments_binary, reference, input_bam, split_region_file_,
                work, work, window_size, maf, min_mapq, max_dp * window_size / 100.0, filter_duplicate_str)
            if calc_qual:
                cmd += " --calculate_qual_stat"
            run_shell_command(cmd, stdout=os.path.join(work, "scan.out"),
                              stderr=os.path.join(work, "scan.err"),
                              run_logger=thread_logger)
        else:
            with open(os.path.join(work, "candidates.vcf"), "w") as o_f:
                pass
            with open(os.path.join(work, "count.bed"), "w") as o_f:
                pass

        pysam.tabix_index(os.path.join(work, "count.bed"), preset="bed")
        concatenate_files([split_region_file],
                          os.path.join(work, "region.bed"))
        return os.path.join(work, "candidates.vcf"), os.path.join(work, "count.bed.gz"), os.path.join(work, "region.bed")
    except Exception as ex:
        thread_logger.error(traceback.format_exc())
        thread_logger.error(ex)
        stderr_file = os.path.join(work, "scan.err")
        if os.path.exists(stderr_file) and os.path.getsize(stderr_file):
            thread_logger.error(
                "Please check error log at {}".format(stderr_file))
        return None


def scan_alignments(work, merge_d_for_scan, scan_alignments_binary, input_bam,
                    regions_bed_file, reference, num_splits,
                    num_threads, window_size, maf, min_mapq, max_dp, filter_duplicate, restart=True,
                    split_region_files=[], calc_qual=True):

    logger = logging.getLogger(scan_alignments.__name__)

    logger.info("-------------------Scan Alignment BAM----------------------")

    split_len_ratio = 0.98
    if not split_region_files:
        if regions_bed_file:
            regions_bed = get_tmp_file()
            with open(regions_bed_file) as i_f, open(regions_bed, "w") as o_f:
                for line in skip_empty(i_f):
                    chrom, st, en = line.strip().split()[0:3]
                    o_f.write("\t".join([chrom, st, en]) + "\n")
            regions_bed = bedtools_sort(regions_bed, run_logger=logger)
            regions_bed = bedtools_merge(
                regions_bed, args=" -d 0", run_logger=logger)
        else:
            regions_bed = get_tmp_file()
            with pysam.AlignmentFile(input_bam, "rb") as samfile:
                with open(regions_bed, "w") as tmpfile:
                    for chrom, length in zip(samfile.references, samfile.lengths):
                        tmpfile.write("{}\t{}\t{}\t.\t.\t.\n".format(
                            chrom, 1, length - 1))
        if not os.path.exists(work):
            os.mkdir(work)
        total_len = 0
        with open(regions_bed) as r_f:
            for line in skip_empty(r_f):
                chrom, st, en = line.strip().split("\t")[0:3]
                total_len += int(en) - int(st) + 1
        if not restart:
            split_region_files = glob.glob(os.path.join(work, "region_*.bed"))
            spilt_total_len = 0
            for split_file in split_region_files:
                with open(split_file) as s_f:
                    for line in skip_empty(s_f):
                        chrom, st, en = line.strip().split("\t")[0:3]
                        spilt_total_len += int(en) - int(st)
            if spilt_total_len >= split_len_ratio * total_len:
                split_region_files = sorted(split_region_files,
                                            key=lambda x: int(
                                                os.path.basename(x).split(".bed")[0].split(
                                                    "_")[1]))
        if not split_region_files:
            regions_bed_file = os.path.join(work, "all_regions.bed")
            shutil.move(regions_bed, regions_bed_file)

            if num_splits is not None:
                num_split = num_splits
            else:
                num_split = max(int(np.ceil((total_len // 10000000) //
                                            num_threads) * num_threads), num_threads)
            split_region_files = split_region(work, regions_bed_file, num_split,
                                              min_region=window_size, max_region=1e20)
    else:
        logger.info("split_regions to be used (will ignore region_bed): {}".format(
            " ".join(split_region_files)))

    map_args = []
    all_outputs = [[]] * len(split_region_files)
    not_done = []
    for i, split_region_file in enumerate(split_region_files):
        if restart or not os.path.exists(os.path.join(work, "work.{}".format(i), "region.bed")) \
                or not os.path.exists(os.path.join(work, "work.{}".format(i), "candidates.vcf")) \
                or not os.path.exists(os.path.join(work, "work.{}".format(i), "count.bed.gz")):
            work_ = os.path.join(work, "work.{}".format(i))
            if os.path.exists(work_):
                shutil.rmtree(work_)
            map_args.append((os.path.join(work, "work.{}".format(i)),
                             reference, merge_d_for_scan, scan_alignments_binary, split_region_file,
                             input_bam, window_size, maf, min_mapq, max_dp, filter_duplicate, calc_qual))
            not_done.append(i)
        else:
            all_outputs[i] = [os.path.join(work, "work.{}".format(i), "candidates.vcf"),
                              os.path.join(work, "work.{}".format(
                                  i), "count.bed.gz"),
                              os.path.join(work, "work.{}".format(i), "region.bed")]

    pool = multiprocessing.Pool(num_threads)
    try:
        outputs = pool.map_async(run_scan_alignments, map_args).get()
        pool.close()
    except Exception as inst:
        pool.close()
        logger.error(inst)
        traceback.print_exc()
        raise Exception

    for o in outputs:
        if o is None:
            raise Exception("scan_alignments failed!")

    for i, output in zip(not_done, outputs):
        all_outputs[i] = output
    return all_outputs

if __name__ == '__main__':
    FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    logging.basicConfig(level=logging.INFO, format=FORMAT)
    logger = logging.getLogger(__name__)

    parser = argparse.ArgumentParser(
        description='simple call variants from bam')
    parser.add_argument('--input_bam', type=str,
                        help='input bam', required=True)
    parser.add_argument('--reference', type=str,
                        help='reference fasta filename', required=True)
    parser.add_argument('--work', type=str,
                        help='work directory', required=True)
    parser.add_argument('--regions_bed_file', type=str,
                        help='regions bed file', default="")
    parser.add_argument('--scan_alignments_binary', type=str,
                        help='binary for scanning alignment bam', default="../bin/scan_alignments")
    parser.add_argument('--window_size', type=int, help='window size to scan the variants',
                        default=2000)
    parser.add_argument('--maf', type=float,
                        help='minimum allele freq', default=0.01)
    parser.add_argument('--min_mapq', type=int,
                        help='minimum mapping quality', default=1)
    parser.add_argument('--max_dp', type=float,
                        help='max depth', default=100000)
    parser.add_argument('--filter_duplicate',
                        help='filter duplicate reads when preparing pileup information',
                        action="store_true")
    parser.add_argument('--merge_d_for_scan', type=int,
                        help='-d used to merge regions before scan',
                        default=None)
    parser.add_argument('--num_splits', type=int,
                        help='number of region splits', default=None)
    parser.add_argument('--num_threads', type=int,
                        help='number of threads', default=1)
    args = parser.parse_args()
    logger.info(args)

    try:
        outputs = scan_alignments(args.work, args.merge_d_for_scan, args.scan_alignments_binary, args.input_bam,
                                  args.regions_bed_file, args.reference, args.num_splits,
                                  args.num_threads, args.window_size, args.maf,
                                  args.min_mapq, args.max_dp, args.filter_duplicate)
    except Exception as e:
        logger.error(traceback.format_exc())
        logger.error("Aborting!")
        logger.error(
            "scan_alignments.py failure on arguments: {}".format(args))
        raise e
