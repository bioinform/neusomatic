#!/usr/bin/env python
#-------------------------------------------------------------------------
# resolve_variants.py
# Resolve variants (e.g. exact INDEL sequences) for target variants
# identified by 'extract_postprocess_targets.py'.
#-------------------------------------------------------------------------

import multiprocessing
import argparse
import re
import logging
import traceback

import pysam
import numpy as np

from utils import get_chromosomes_order, read_tsv_file, bedtools_sort, bedtools_merge, get_tmp_file, skip_empty
from defaults import VCF_HEADER

CIGAR_MATCH = 0
CIGAR_INS = 1
CIGAR_DEL = 2
CIGAR_SOFTCLIP = 4
CIGAR_HARDCLIP = 5
CIGAR_EQUAL = 7
CIGAR_DIFF = 8

_CIGAR_OPS = "MIDNSHP=X"
_CIGAR_PATTERN = re.compile(r'([0-9]+)([MIDNSHPX=])')
_CIGAR_OP_DICT = {op: index for index, op in enumerate(_CIGAR_OPS)}
_CIGAR_REFERENCE_OPS = [CIGAR_MATCH, CIGAR_DEL, CIGAR_EQUAL, CIGAR_DIFF]
_CIGAR_READ_ALN_OPS = [CIGAR_MATCH, CIGAR_INS, CIGAR_EQUAL, CIGAR_DIFF]
NUC_to_NUM = {"A": 1, "C": 2, "G": 3, "T": 4, "-": 0, "N": 5}
NUM_to_NUC = {1: "A", 2: "C", 3: "G", 4: "T", 0: "-", 5: "N"}


def extract_del(record):
    logger = logging.getLogger(extract_del.__name__)
    dels = []
    pos = record.pos
    for C, L in record.cigartuples:
        if C in [CIGAR_SOFTCLIP, CIGAR_HARDCLIP, CIGAR_INS]:
            continue
        if C == CIGAR_DEL:
            dels.append([record.reference_name, pos, pos + L])
        pos += L

    return dels


def extract_ins(record):
    logger = logging.getLogger(extract_ins.__name__)
    inss = []
    pos = record.pos
    seq_pos = 0
    for C, L in record.cigartuples:
        if C == CIGAR_SOFTCLIP:
            seq_pos += L
            continue
        elif C == CIGAR_HARDCLIP:
            continue
        if C == CIGAR_INS:
            inss.append([record.reference_name, pos, pos + 1,
                         record.query[seq_pos:seq_pos + L]])
            seq_pos += L
        else:
            if C != CIGAR_DEL:
                seq_pos += L
            pos += L
    return inss


def find_resolved_variants(input_record):
    chrom, start, end, variants, input_bam, filter_duplicate, reference = input_record
    thread_logger = logging.getLogger(
        "{} ({})".format(find_resolved_variants.__name__, multiprocessing.current_process().name))
    try:
        ref = pysam.FastaFile(reference)
        out_variants = []
        start, end = list(map(int, [start, end]))
        region = [chrom, start, end]
        vartypes = list(map(lambda x: x[-1], variants))
        scores = list(map(lambda x: x[5], variants))
        if len(set(vartypes)) > 1:
            out_variants.extend(
                list(map(lambda x: [x[0], int(x[1]), x[3], x[4], x[9].split(":")[0], x[5]], variants)))
        else:
            vartype = vartypes[0]
            score = max(scores)
            if vartype == "DEL":
                dels = []
                with pysam.AlignmentFile(input_bam) as samfile:
                    for record in samfile.fetch(chrom, start, end):
                        if not record.is_duplicate or not filter_duplicate:
                            if record.cigarstring and "D" in record.cigarstring:
                                dels.extend(extract_del(record))
                dels = list(filter(lambda x: (
                    start <= x[1] <= end) or start <= x[2] <= end, dels))
                if dels:
                    del_strs = list(
                        map(lambda x: "---".join(map(str, x[0:3])), dels))
                    uniq_dels = list(set(del_strs))
                    uniq_dels_count = {}
                    for del_ in uniq_dels:
                        uniq_dels_count[del_] = del_strs.count(del_)
                    max_count = max(uniq_dels_count.values())
                    for del_ in uniq_dels:
                        if uniq_dels_count[del_] <= max_count * 0.5:
                            del uniq_dels_count[del_]
                    new_bed = get_tmp_file()
                    with open(new_bed, "w") as f_o:
                        for k in uniq_dels_count.keys():
                            x = k.split("---")
                            f_o.write(
                                "\t".join(map(str, x + [".", "."])) + "\n")
                    new_bed = bedtools_sort(new_bed, run_logger=thread_logger)
                    new_bed = bedtools_merge(
                        new_bed, args=" -c 1 -o count", run_logger=thread_logger)
                    vs = read_tsv_file(new_bed, fields=range(4))
                    vs = list(map(lambda x: [x[0], int(x[1]), ref.fetch(x[0], int(
                        x[1]) - 1, int(x[2])).upper(), ref.fetch(x[0], int(x[1]) - 1, int(x[1])).upper(), "0/1", score], vs))
                    out_variants.extend(vs)
            elif vartype == "INS":
                intervals = []
                inss = []
                with pysam.AlignmentFile(input_bam) as samfile:
                    for record in samfile.fetch(chrom, start, end):
                        if not record.is_duplicate or not filter_duplicate:
                            if record.cigarstring and "I" in record.cigarstring:
                                inss.extend(extract_ins(record))
                inss = list(filter(lambda x: (
                    start <= x[1] <= end) or start <= x[2] <= end, inss))
                if inss:
                    ins_strs = list(
                        map(lambda x: "---".join(map(str, x[0:4])), inss))
                    uniq_inss = list(set(ins_strs))
                    uniq_inss_count = {}
                    for ins_ in uniq_inss:
                        uniq_inss_count[ins_] = ins_strs.count(ins_)
                    max_ins, max_count = sorted(
                        uniq_inss_count.items(), key=lambda x: x[1])[-1]
                    max_pos = int(max_ins.split("---")[1])
                    for ins_ in uniq_inss:
                        if uniq_inss_count[ins_] <= max_count * 0.5 or 0 < abs(int(ins_.split("---")[1]) - max_pos) < 4:
                            del uniq_inss_count[ins_]

                    new_bed = get_tmp_file()
                    with open(new_bed, "w") as f_o:
                        for k in uniq_inss_count.keys():
                            x = k.split("---")
                            f_o.write(
                                "\t".join(map(str, x + [".", "."])) + "\n")
                    new_bed = bedtools_sort(new_bed, run_logger=thread_logger)
                    vs = read_tsv_file(new_bed, fields=range(4))
                    vs = list(map(lambda x: [x[0], int(x[1]), ref.fetch(x[0], int(
                        x[1]) - 1, int(x[1])).upper(), ref.fetch(x[0], int(x[1]) - 1, int(x[1])).upper() + x[3], "0/1", score], vs))
                    out_variants.extend(vs)
        return out_variants
    except Exception as ex:
        thread_logger.error(traceback.format_exc())
        thread_logger.error(ex)
        return None


def resolve_variants(input_bam, resolved_vcf, reference, target_vcf_file,
                     target_bed_file, filter_duplicate, num_threads):
    logger = logging.getLogger(resolve_variants.__name__)

    logger.info("-------Resolve variants (e.g. exact INDEL sequences)-------")

    variants = {}
    with open(target_vcf_file) as tv_f:
        for line in skip_empty(tv_f):
            fields = line.strip().split()
            id_ = int(fields[2])
            if len(fields[4]) < len(fields[3]):
                vartype = "DEL"
            elif len(fields[4]) > len(fields[3]):
                vartype = "INS"
            else:
                vartype = "SNP"
            if id_ not in variants:
                variants[id_] = []
            variants[id_].append(fields + [vartype])
    map_args = []
    with open(target_bed_file) as i_f:
        for line in skip_empty(i_f):
            tb = line.strip().split("\t")
            chrom, start, end, id_ = tb[0:4]
            id_ = int(id_)
            map_args.append([chrom, start, end, variants[id_],
                             input_bam, filter_duplicate, reference])

    pool = multiprocessing.Pool(num_threads)
    try:
        out_variants_list = pool.map_async(
            find_resolved_variants, map_args).get()
        pool.close()
    except Exception as inst:
        logger.error(inst)
        pool.close()
        traceback.print_exc()
        raise Exception

    for o in out_variants_list:
        if o is None:
            raise Exception("resolve_variants failed!")

    out_variants = [x for xs in out_variants_list for x in xs]
    chroms_order = get_chromosomes_order(bam=input_bam)

    out_variants = sorted(out_variants, key=lambda x: [
                          chroms_order[x[0]], int(x[1])])
    with open(resolved_vcf, "w") as o_f:
        o_f.write("{}\n".format(VCF_HEADER))
        o_f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
        for chrom, pos, ref, alt, gt, phred_score in out_variants:
            if ref != alt:
                phred_score = float(phred_score)
                prob = np.round(1 - (10**(-phred_score / 10)), 4)
                o_f.write("\t".join([chrom, str(pos), ".", ref,
                                     alt, "{:.4f}".format(
                                         np.round(phred_score, 4)),
                                     ".", "SCORE={:.4f}".format(prob), "GT", gt]) + "\n")


if __name__ == '__main__':

    FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    logging.basicConfig(level=logging.INFO, format=FORMAT)
    logger = logging.getLogger(__name__)

    parser = argparse.ArgumentParser(
        description='Resolve ambigues variants for high quality reads')
    parser.add_argument('--input_bam', type=str,
                        help='input bam', required=True)
    parser.add_argument('--resolved_vcf', type=str,
                        help='resolved_vcf', required=True)
    parser.add_argument('--target_vcf', type=str,
                        help='resolve target vcf', required=True)
    parser.add_argument('--target_bed', type=str,
                        help='resolve target bed', required=True)
    parser.add_argument('--reference', type=str,
                        help='reference fasta filename', required=True)
    parser.add_argument('--filter_duplicate',
                        help='filter duplicate reads in analysis',
                        action="store_true")
    parser.add_argument('--num_threads', type=int,
                        help='number of threads', default=1)
    args = parser.parse_args()

    try:
        resolve_variants(args.input_bam, args.resolved_vcf,
                         args.reference, args.target_vcf,
                         args.target_bed, args.filter_duplicate,
                         args.num_threads)
    except Exception as e:
        logger.error(traceback.format_exc())
        logger.error("Aborting!")
        logger.error(
            "resolve_variants.py failure on arguments: {}".format(args))
        raise e
