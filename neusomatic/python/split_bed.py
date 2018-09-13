#-------------------------------------------------------------------------
# split_bed.py
# split bed file to multiple sub-regions
#-------------------------------------------------------------------------
import os
import argparse
import traceback
from random import shuffle
import logging

import pybedtools


def split_region(work, region_bed_file, num_splits, max_region=1000000, min_region=20, shuffle_intervals=False):

    logger = logging.getLogger(split_region.__name__)

    logger.info("-----------------------------------------------------------")
    logger.info("Split region")
    logger.info("-----------------------------------------------------------")

    regions_bed = pybedtools.BedTool(region_bed_file).sort().merge(d=0)
    intervals = []
    for region in regions_bed:
        chrom, start, end = region.chrom, region.start, region.end
        if region.length + 1 > max_region:
            for i in range(start, end + 1, max_region):
                intervals.append(pybedtools.Interval(
                    chrom, i, min(end, i + max_region - 1)))
        else:
            intervals.append(region)
    if shuffle_intervals:
        shuffle(intervals)
    regions_bed = pybedtools.BedTool(intervals)
    total_len = sum(map(lambda x: int(x[2]) - int(x[1]) + 1, regions_bed))
    logger.info("Total length: {}".format(total_len))
    split_len = total_len / num_splits
    split_regions = []
    current_regions = []
    sofar_len = 0
    current_len = 0
    split_lens = []
    for region in regions_bed:
        chrom, start, end = region[0:3]
        start, end = int(start), int(end)
        s = start
        e = -1
        while(current_len < split_len):
            s = max(s, e + 1)
            e = min(s + split_len - current_len - 1, end)
            if (e - s + 1) < 2 * min_region:
                e = min(s + 2 * min_region - 1, end)
            if (end - e) < 2 * min_region:
                e = end
            current_regions.append(pybedtools.Interval(chrom, s, e))
            current_len += e - s + 1
            if (current_len >= split_len):
                sofar_len += current_len
                split_lens.append(current_len)
                current_len = 0
                split_regions.append(current_regions)
                current_regions = []
                if split_len < (total_len - sofar_len) < 1.5 * split_len:
                    split_len = total_len - sofar_len
            if e >= end:
                break
    if current_regions:
        split_lens.append(current_len)
        split_regions.append(current_regions)
    split_region_files = []
    for i, split_region_ in enumerate(split_regions):
        split_region_file = os.path.join(work, "region_{}.bed".format(i))
        pybedtools.BedTool(split_region_).saveas(split_region_file)
        logger.info("Split {}: {}".format(i, split_lens[i]))
        split_region_files.append(split_region_file)
    sum_len = sum(split_lens)
    logger.info("Total splitted length: {}".format(sum_len))
    return split_region_files

if __name__ == '__main__':

    FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    logging.basicConfig(level=logging.INFO, format=FORMAT)
    logger = logging.getLogger(__name__)

    parser = argparse.ArgumentParser(
        description='Split bedfile to multiple beds')
    parser.add_argument('--region_bed', type=str,
                        help='region bed', required=True)
    parser.add_argument('--num_splits', type=int,
                        help='number of splits', required=True)
    parser.add_argument('--output', type=str,
                        help='output directory', required=True)
    parser.add_argument('--max_region', type=int,
                        help='max region size in the bed (for shuffling purpose)', default=1000000)
    parser.add_argument('--min_region', type=int,
                        help='min region size for spliting', default=20)
    args = parser.parse_args()

    if not os.path.exists(args.output):
        os.mkdir(args.output)
    work = os.path.join(args.output, "splits")
    if not os.path.exists(work):
        os.mkdir(work)

    try:
        split_region_files = split_region(
            work, args.region_bed, args.num_splits, args.max_region, args.min_region)
    except Exception as e:
        logger.error(traceback.format_exc())
        logger.error("Aborting!")
        logger.error(
            "split_bed.py failure on arguments: {}".format(args))
        raise e
