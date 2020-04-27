#!/usr/bin/env python
#-------------------------------------------------------------------------
# merge_post_vcfs.py
# Merge resolved variants and other predicted variants and output the final NeuSomatic .vcf
#-------------------------------------------------------------------------
import argparse
import traceback
import fileinput
import logging

import numpy as np

from utils import get_chromosomes_order, skip_empty


def merge_post_vcfs(ref, resolved_vcf, no_resolve_vcf, out_vcf,
                    pass_threshold, lowqual_threshold):

    logger = logging.getLogger(merge_post_vcfs.__name__)

    logger.info("------------------------Merge vcfs-------------------------")

    chroms_order = get_chromosomes_order(reference=ref)

    good_records = []
    # Use fileinput to stream as if the two files were concatenated
    for line in skip_empty(fileinput.input([no_resolve_vcf, resolved_vcf])):
        chrom, pos, _, ref, alt, score, _, info, format_, gt = line.strip().split()
        good_records.append([chrom, pos, ref, alt, gt, score])

    with open(out_vcf, "w") as o_f:
        o_f.write("##fileformat=VCFv4.2\n")
        o_f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
        for record in sorted(good_records, key=lambda x: [chroms_order[x[0]], x[1]]):
            chrom, pos, ref, alt, gt, score = record
            prob = np.round(1 - (10**(-float(score) / 10)), 4)
            filter_ = "REJECT"
            if prob >= pass_threshold:
                filter_ = "PASS"
            elif prob >= lowqual_threshold:
                filter_ = "LowQual"
            o_f.write("\t".join([chrom, pos, ".", ref, alt,
                                 "{:.4f}".format(
                                     float(score)), filter_, "SCORE={:.4f}".format(prob),
                                 "GT", "0/1"]) + "\n")

if __name__ == '__main__':

    FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    logging.basicConfig(level=logging.INFO, format=FORMAT)
    logger = logging.getLogger(__name__)

    parser = argparse.ArgumentParser(description='merge pred vcfs')
    parser.add_argument(
        '--ref', help='reference fasta filename', required=True)
    parser.add_argument('--resolved_vcf', help='resolved_vcf', required=True)
    parser.add_argument('--no_resolve_vcf',
                        help='no resolve vcf', required=True)
    parser.add_argument('--out_vcf', help='output vcf', required=True)
    args = parser.parse_args()
    logger.info(args)
    try:
        merge_post_vcfs(args.ref, args.resolved_vcf,
                        args.no_resolve_vcf, args.out_vcf)
    except Exception as e:
        logger.error(traceback.format_exc())
        logger.error("Aborting!")
        logger.error(
            "merge_post_vcfs.py failure on arguments: {}".format(args))
        raise e
