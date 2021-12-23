import multiprocessing
import argparse
import os
import shutil
import traceback
import logging

import pybedtools

import sys, os, re, pysam
import scipy.stats as stats
import numpy as np

def exclude(input_vcf, exclude_vcf, output_vcf):
    ex_ids=[]
    for x in pybedtools.BedTool(exclude_vcf):
        chrom,pos,_,ref,alt=x[0:5]
        id_="-".join(map(str,[chrom,pos,ref,alt]))
        ex_ids.append(id_)

    with open(output_vcf,"w") as o_f:
        o_f.write("##fileformat=VCFv4.2\n")
        o_f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
        for x in pybedtools.BedTool(input_vcf):
            chrom,pos,_,ref,alt=x[0:5]
            id_="-".join(map(str,[chrom,pos,ref,alt]))
            if id_ not in ex_ids:
                o_f.write("\t".join(map(str,x))+"\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='exclude from vcf')
    parser.add_argument('--input', type=str, help='input vcf (or vcf.gz)',
                        required=True)
    parser.add_argument('--ex', type=str,
                        help='exclude vcf (or vcf.gz)', required=True)
    parser.add_argument('--output', type=str,
                        help='output vcf', required=True)
    args = parser.parse_args()
    print args

    try:
        exclude(args.input, args.ex, args.output)
    except Exception as e:
        print traceback.format_exc()
        print "Aborting!"
        print "preprocess.py failure on arguments: {}".format(args)
        raise e
