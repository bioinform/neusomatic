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

import pybedtools
import pysam
import numpy as np

from utils import get_chromosomes_order

FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logFormatter = logging.Formatter(FORMAT)
logger = logging.getLogger(__name__)
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
logger.addHandler(consoleHandler)
logging.getLogger().setLevel(logging.INFO)

logger = logging.getLogger(__name__)

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


def find_resolved_variants((chrom, start, end, variant, input_bam, ref)):
    out_variants = []
    start, end = map(int, [start, end])
    region = [chrom, start, end]
    vartypes = map(lambda x: x[-1], variant)
    scores = map(lambda x: x[5], variant)
    if len(set(vartypes)) > 1:
        out_variants.extend(
            map(lambda x: [x[0], int(x[1]), x[3], x[4], x[10], x[5]], variant))
    else:
        vartype = vartypes[0]
        score = scores[0]
        if vartype == "DEL":
            intervals = []
            dels = []
            with pysam.AlignmentFile(input_bam) as samfile:
                for record in samfile.fetch(chrom, start, end):
                    if record.cigarstring and "D" in record.cigarstring:
                        dels.extend(extract_del(record))
            dels = filter(lambda x: (
                start <= x[1] <= end) or start <= x[2] <= end, dels)
            if dels:
                intervals = map(lambda x: pybedtools.Interval(
                    x[0], x[1], x[2]), dels)
                bed = pybedtools.BedTool(intervals)
                del_strs = map(lambda x: "---".join(x[0:3]), bed)
                uniq_dels = list(set(del_strs))
                uniq_dels_count = {}
                for del_ in uniq_dels:
                    uniq_dels_count[del_] = del_strs.count(del_)
                max_count = max(uniq_dels_count.values())
                for del_ in uniq_dels:
                    if uniq_dels_count[del_] <= max_count * 0.5:
                        del uniq_dels_count[del_]
                new_bed = pybedtools.BedTool(map(lambda x: pybedtools.Interval(x[0], int(x[1]), int(x[2])),
                                                 map(lambda x: x.split("---"), uniq_dels_count.keys())))
                new_bed = new_bed.sort().merge(c=[1], o="count")
                out_variants.extend(map(lambda x: [x[0], int(x[1]), ref.fetch(x[0], int(
                    x[1]) - 1, int(x[2])), ref.fetch(x[0], int(x[1]) - 1, int(x[1])), "0/1", score], new_bed))
        elif vartype == "INS":
            intervals = []
            inss = []
            with pysam.AlignmentFile(input_bam) as samfile:
                for record in samfile.fetch(chrom, start, end):
                    if record.cigarstring and "I" in record.cigarstring:
                        inss.extend(extract_ins(record))
            inss = filter(lambda x: (
                start <= x[1] <= end) or start <= x[2] <= end, inss)
            if inss:
                intervals = map(lambda x: pybedtools.Interval(
                    x[0], x[1], x[2], x[3]), inss)
                bed = pybedtools.BedTool(intervals)
                ins_strs = map(lambda x: "---".join(x[0:4]), bed)
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
                new_bed = pybedtools.BedTool(map(lambda x: pybedtools.Interval(x[0], int(x[1]), int(x[2]), x[3]),
                                                 map(lambda x: x.split("---"), uniq_inss_count.keys()))).sort()
                out_variants.extend(map(lambda x: [x[0], int(x[1]), ref.fetch(x[0], int(
                    x[1]) - 1, int(x[1])), ref.fetch(x[0], int(x[1]) - 1, int(x[1])) + x[3], "0/1", score], new_bed))
    return out_variants


def resolve_variants(input_bam, resolved_vcf, reference, target_vcf_file,
                     target_bed_file, num_threads):

    logger.info("-----------------------------------------------------------")
    logger.info("Resolve variants (e.g. exact INDEL sequences)")
    logger.info("-----------------------------------------------------------")

    ref = pysam.FastaFile(reference)
    variants = {}
    with open(target_vcf_file) as tv_f:
        for line in tv_f:
            if line[0] == "#":
                continue
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
    target_bed = pybedtools.BedTool(target_bed_file)
    map_args = []
    for tb in target_bed:
        chrom, start, end, id_ = tb[0:4]
        id_ = int(id_)
        map_args.append([chrom, start, end, variants[id_], input_bam, ref])

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
    out_variants = [x for xs in out_variants_list for x in xs]
    chroms_order = get_chromosomes_order(bam=input_bam)

    out_variants = sorted(out_variants, key=lambda x: [
                          chroms_order[x[0]], int(x[1])])
    with open(resolved_vcf, "w") as o_f:
        o_f.write("##fileformat=VCFv4.2\n")
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
    parser.add_argument('--num_threads', type=int,
                        help='number of threads', default=1)
    args = parser.parse_args()

    resolve_variants(args.input_bam, args.resolved_vcf,
                     args.reference, args.target_vcf,
                     args.target_bed, args.num_threads)
