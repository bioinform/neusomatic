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

import pybedtools
import pysam
import numpy as np

from utils import concatenate_files, run_shell_command
from split_bed import split_region


FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logFormatter = logging.Formatter(FORMAT)
logger = logging.getLogger(__name__)
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
logger.addHandler(consoleHandler)
logging.getLogger().setLevel(logging.INFO)


def run_scan_alignments((work, reference, scan_alignments_binary, split_region_file,
                         input_bam, window_size, maf, min_mapq, calc_qual, num_threads)):
    if not os.path.exists(work):
        os.mkdir(work)
    if len(pybedtools.BedTool(split_region_file)) > 0:
        cmd = "{} --ref {} -b {} -L {} --out_vcf_file {}/candidates.vcf --out_count_file {}/count.bed \
                    --window_size {} --min_af {} --min_mapq {} --num_thread {}".format(
            scan_alignments_binary, reference, input_bam, split_region_file,
            work, work, window_size, maf, min_mapq, num_threads)
        if calc_qual:
            cmd += " --calculate_qual_stat"
        ret_code = run_shell_command(cmd, stdout=os.path.join(work, "scan.out"),
                                     stderr=os.path.join(work, "scan.err"))
        if ret_code != 0:
            logger.error("scan_alignments failed")
            logger.error("Aborting!")
            raise Exception(
                "scan_alignments failure on command: {}".format(cmd))
    else:
        pybedtools.BedTool([]).saveas(os.path.join(work, "candidates.vcf"))
        pybedtools.BedTool([]).saveas(os.path.join(work, "count.bed"))

    pysam.tabix_index(os.path.join(work, "count.bed"), preset="bed")
    concatenate_files([split_region_file], os.path.join(work, "region.bed"))
    return os.path.join(work, "candidates.vcf"), os.path.join(work, "count.bed.gz"), os.path.join(work, "region.bed")


def scan_alignments(work, scan_alignments_binary, input_bam,
                    regions_bed_file, reference,
                    num_threads, window_size, maf, min_mapq, restart=True,
                    split_region_files=[], calc_qual=True):

    logger.info("-----------------------------------------------------------")
    logger.info("Scan Alignment BAM")
    logger.info("-----------------------------------------------------------")

    if not split_region_files:
        if regions_bed_file:
            regions_bed = pybedtools.BedTool(
                regions_bed_file).sort().merge(d=0)
        else:
            intervals = []
            with pysam.AlignmentFile(input_bam, "rb") as samfile:
                for chrom, length in zip(samfile.references, samfile.lengths):
                    intervals.append(pybedtools.Interval(chrom, 1, length - 1))
            regions_bed = pybedtools.BedTool(intervals)
        if not os.path.exists(work):
            os.mkdir(work)
        total_len = sum(map(lambda x: int(x[2]) - int(x[1]) + 1, regions_bed))

        if not restart:
            split_region_files = glob.glob(os.path.join(work, "region_*.bed"))
            spilt_total_len = sum(map(lambda x: sum(
                [y.length for y in pybedtools.BedTool(x)]), split_region_files))
            if spilt_total_len >= 0.98 * total_len:
                split_region_files = sorted(split_region_files,
                                            key=lambda x: int(
                                                os.path.basename(x).split(".bed")[0].split(
                                                    "_")[1]))
        if not split_region_files:
            regions_bed_file = os.path.join(work, "all_regions.bed")
            regions_bed.saveas(regions_bed_file)

            num_split = max(int(np.ceil((total_len / 10000000) /
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
            map_args.append((os.path.join(work, "work.{}".format(i)),
                             reference, scan_alignments_binary, split_region_file,
                             input_bam, window_size, maf, min_mapq, calc_qual, 1))
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
    for i, output in zip(not_done, outputs):
        all_outputs[i] = output
    return all_outputs

if __name__ == '__main__':
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
    parser.add_argument('--num_threads', type=int,
                        help='number of threads', default=1)
    args = parser.parse_args()
    logger.info(args)

    try:
        outputs = scan_alignments(args.work, args.scan_alignments_binary, args.input_bam,
                                  args.regions_bed_file, args.reference,
                                  args.num_threads, args.window_size, args.maf,
                                  args.min_mapq)
    except:
        traceback.print_exc()
        logger.error("Aborting!")
        raise Exception(
            "scan_alignments.py failure on arguments: {}".format(args))
