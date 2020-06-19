#!/usr/bin/env python
#-------------------------------------------------------------------------
# extract_postprocess_targets.py
# Extract variants that need postprocessing (like larege INDELs)
# from the output predictions of 'call.py'.
#-------------------------------------------------------------------------
import argparse
import traceback
import logging
import pysam

from utils import skip_empty
from defaults import VCF_HEADER
from resolve_variants import push_left_var


def extract_postprocess_targets(reference, input_vcf, min_len, max_dist, pad):
    logger = logging.getLogger(extract_postprocess_targets.__name__)

    logger.info("--------------Extract Postprocessing Targets---------------")

    ref_fasta = pysam.FastaFile(reference)

    base_name = ".".join(input_vcf.split(".")[:-1])
    out_vcf = "{}.no_resolve.vcf".format(base_name)
    redo_vcf = "{}.resolve_target.vcf".format(base_name)
    redo_bed = "{}.resolve_target.bed".format(base_name)

    record_sets = []
    record_set = []
    with open(input_vcf) as i_f, open(out_vcf, "w") as o_f, open(redo_vcf, "w") as r_f, open(redo_bed, "w") as r_b:
        r_f.write("{}\n".format(VCF_HEADER))
        r_f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
        for line in skip_empty(i_f):
            if len(line) < 2:
                continue

            chrom, pos, _, ref, alt, _, _, _, _, gt = line.strip().split()
            pos = int(pos)
            record = [chrom, pos, ref, alt, gt, line]
            if not record_set:
                record_set.append(record)
                continue
            chrom_, pos_, ref_, alt_ = push_left_var(
                ref_fasta, chrom, pos, ref, alt)
            if len(list(filter(lambda x: (chrom == x[0] and
                                          (min(abs(x[1] + len(x[2]) - (pos + len(ref))),
                                               abs(x[1] - pos),
                                               abs(min(x[1] + len(x[2]), pos + len(ref)) - max(x[1], pos))) <= max_dist)), record_set))) > 0 or len(
                list(filter(lambda x: (chrom_ == x[0] and
                                       (min(abs(x[1] + len(x[2]) - (pos_ + len(ref_))),
                                            abs(x[1] - pos_),
                                            abs(min(x[1] + len(x[2]), pos_ + len(ref_)) - max(x[1], pos_))) <= max_dist)), record_set))) > 0:
                record_set.append(record)
                continue


            if record_set:
                record_sets.append(record_set)
            record_set = [record]
        record_sets.append(record_set)

        for ii, record_set in enumerate(record_sets):
            if len(record_set) > 1:
                varid_pos = {}
                for chrom, pos, ref, alt, _, _ in record_set:
                    if pos not in varid_pos:
                        varid_pos[pos] = set([])
                    vid = "-".join([ref, alt])
                    varid_pos[pos].add(vid)
                multi_allelic = False
                for vid in varid_pos:
                    if len(varid_pos[vid]) > 1:
                        multi_allelic = True

                if list(filter(lambda x: len(x[2]) != len(x[3]), record_set)) or multi_allelic:
                    for x in record_set:
                        fields = x[-1].strip().split()
                        fields[2] = str(ii)
                        r_f.write("\t".join(fields) + "\n")
                    r_b.write("\t".join(map(str, [record_set[0][0], max(0, min(map(lambda x:x[1], record_set)) - pad),
                                                  max(map(lambda x:x[
                                                      1] + len(x[2]), record_set)) + pad, ii,
                                                  ])) + "\n")
                else:
                    for x in record_set:
                        o_f.write(x[-1])

            elif record_set:
                if abs(len(record_set[0][2]) - len(record_set[0][3])) >= min_len:
                    fields = record_set[0][-1].strip().split()
                    fields[2] = str(ii)
                    r_f.write("\t".join(fields) + "\n")
                    chrom_, pos_, ref_, alt_ = record_set[0][0:4]
                    r_b.write("\t".join(
                        map(str, [chrom_, max(0, pos_ - pad), pos_ + len(ref_) + pad, ii])) + "\n")
                else:
                    o_f.write(record_set[0][-1])


if __name__ == '__main__':

    FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    logging.basicConfig(level=logging.INFO, format=FORMAT)
    logger = logging.getLogger(__name__)

    parser = argparse.ArgumentParser(
        description='infer genotype by ao and ro counts')
    parser.add_argument('--reference', type=str,
                        help='reference fasta filename', required=True)
    parser.add_argument('--input_vcf', type=str,
                        help='input vcf', required=True)
    parser.add_argument('--min_len', type=int,
                        help='minimum INDEL len to resolve', default=4)
    parser.add_argument('--max_dist', type=int,
                        help='max distance to neighboring variant', default=5)
    parser.add_argument(
        '--pad', type=int, help='padding to bed region for extracting reads', default=10)
    args = parser.parse_args()
    logger.info(args)
    try:
        extract_postprocess_targets(
            args.reference, args.input_vcf, args.min_len, args.max_dist, args.pad)
    except Exception as e:
        logger.error(traceback.format_exc())
        logger.error("Aborting!")
        logger.error(
            "extract_postprocess_targets.py failure on arguments: {}".format(args))
        raise e
