#-------------------------------------------------------------------------
# merge_post_vcfs.py
# Merge resolved variants and other predicted variants and output the final NeuSomatic .vcf
#-------------------------------------------------------------------------
import argparse
import traceback
import fileinput
import logging

import numpy as np

from utils import get_chromosomes_order

from _version import __version__

FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logFormatter = logging.Formatter(FORMAT)
logger = logging.getLogger(__name__)
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
logger.addHandler(consoleHandler)
logging.getLogger().setLevel(logging.INFO)


def merge_post_vcfs(ref, resolved_vcf, no_resolve_vcf, target_vcf, out_vcf,
                    pass_threshold, lowqual_threshold):

    logger.info("-----------------------------------------------------------")
    logger.info("Merge vcfs")
    logger.info("-----------------------------------------------------------")

    chroms_order = get_chromosomes_order(reference=ref)

    good_records = []
    # Use fileinput to stream as if the two files were concatenated
    for line in fileinput.input([no_resolve_vcf, resolved_vcf]):
        if line[0] == "#":
            continue
        chrom, pos, _, ref, alt, score, _, info, format_, gt = line.strip().split()
        good_records.append([chrom, pos, ref, alt, gt, score])

    with open(out_vcf, "w") as o_f:
        o_f.write("##fileformat=VCFv4.2\n")
        o_f.write("##NeuSomatic Version={}\n".format(__version__))
        o_f.write(
            "##FORMAT=<ID=SCORE,Number=1,Type=Float,Description=\"Prediction probability score\">\n")
        o_f.write("##FILTER=<ID=PASS,Description=\"Accept as a higher confidence somatic mutation calls with probability score value at least {}\">\n".format(
            pass_threshold))
        o_f.write("##FILTER=<ID=LowQual,Description=\"Less confident somatic mutation calls with probability score value at least {}\">\n".format(
            lowqual_threshold))
        o_f.write("##FILTER=<ID=REJECT,Description=\"Rejected as a confident somatic mutation with probability score value below {}\">\n".format(
            lowqual_threshold))
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
    parser = argparse.ArgumentParser(description='merge pred vcfs')
    parser.add_argument(
        '--ref', help='reference fasta filename', required=True)
    parser.add_argument('--resolved_vcf', help='resolved_vcf', required=True)
    parser.add_argument('--no_resolve_vcf',
                        help='no resolve vcf', required=True)
    parser.add_argument(
        '--target_vcf', help='resolve target vcf', required=True)
    parser.add_argument('--out_vcf', help='output vcf', required=True)
    args = parser.parse_args()
    logger.info(args)
    try:
        merge_post_vcfs(args.ref, args.resolved_vcf,
                        args.no_resolve_vcf, args.target_vcf, args.out_vcf)
    except:
        traceback.print_exc()
        logger.error("Aborting!")
        raise Exception(
            "merge_post_vcfs.py failure on arguments: {}".format(args))
