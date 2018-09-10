#-------------------------------------------------------------------------
# resolve_score.py
# resolve prediction scores for realigned variants
#-------------------------------------------------------------------------

import argparse
import logging
import traceback

import pybedtools
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


def resolve_scores(input_bam, ra_vcf, target_vcf, output_vcf):

    logger.info("-----------------------------------------------------------")
    logger.info("Resolve Prediction Scores for Realigned Variants")
    logger.info("-----------------------------------------------------------")

    ra_out = pybedtools.BedTool(ra_vcf)
    ra_target = pybedtools.BedTool(target_vcf)

    final_intervals = []
    for interval in ra_out.window(ra_target, w=5, v=True):
        interval[5] = "0.5"
        final_intervals.append(interval)

    intervals_dict = {}
    for interval in ra_out.window(ra_target, w=5):
        id_ = "{}-{}-{}-{}".format(interval[0],
                                   interval[1], interval[3], interval[4])
        if id_ not in intervals_dict:
            intervals_dict[id_] = []
        intervals_dict[id_].append(interval)

    for id_, intervals in intervals_dict.iteritems():
        if len(intervals) == 1:
            score = intervals[0][15]
            interval = intervals[0][:10]
            interval[5] = score
            interval[7] = "SCORE={:.4f}".format(
                np.round(1 - (10**(-float(score) / 10)), 4))
        else:
            len_ = (len(intervals[0][4]) - len(intervals[0][3]))
            pos_ = int(intervals[0][1])
            len_diff = map(lambda x: abs(
                (len(x[14]) - len(x[13])) - len_), intervals)
            min_len_diff = min(len_diff)
            intervals = filter(lambda x: abs(
                (len(x[14]) - len(x[13])) - len_) == min_len_diff, intervals)
            pos_diff = map(lambda x: abs(int(x[11]) - pos_), intervals)
            min_pos_diff = min(pos_diff)
            intervals = filter(lambda x: abs(
                int(x[11]) - pos_) == min_pos_diff, intervals)
            score = "{:.4f}".format(
                np.round(max(map(lambda x: float(x[15]), intervals)), 4))
            interval = intervals[0][:10]
            interval[5] = score
            interval[7] = "SCORE={:.4f}".format(
                np.round(1 - (10**(-float(score) / 10)), 4))
        final_intervals.append(interval)

    chroms_order = get_chromosomes_order(bam=input_bam)

    out_variants = sorted(final_intervals, key=lambda x: [
                          chroms_order[x[0]], int(x[1])])
    with open(output_vcf, "w") as o_f:
        o_f.write("##fileformat=VCFv4.2\n")
        o_f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
        for var in out_variants:
            o_f.write("\t".join(var) + "\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Resolve scores')
    parser.add_argument('--input_bam', type=str,
                        help='input bam', required=True)
    parser.add_argument('--ra_vcf', type=str,
                        help='realigned vcf', required=True)
    parser.add_argument('--target_vcf', type=str,
                        help='target vcf', required=True)
    parser.add_argument('--output_vcf', type=str,
                        help='output_vcf', required=True)
    args = parser.parse_args()

    try:
        resolve_scores(args.input_bam, args.ra_vcf,
                       args.target_vcf, args.output_vcf)
    except:
        traceback.print_exc()
        logger.error("Aborting!")
        raise Exception(
            "resolve_scores.py failure on arguments: {}".format(args))
