#!/usr/bin/env python
#-------------------------------------------------------------------------
# filter_candidates.py
# filter raw candidates extracted by 'scan_alignments.py' using min_af and other cut-offs
#-------------------------------------------------------------------------
import os
import argparse
import traceback
import logging
import multiprocessing

import pysam
import numpy as np

from utils import safe_read_info_dict, run_bedtools_cmd, vcf_2_bed, write_tsv_file, bedtools_sort, get_tmp_file, skip_empty
from defaults import VCF_HEADER


def filter_candidates(candidate_record):
    candidates_vcf, filtered_candidates_vcf, reference, dbsnp, min_dp, max_dp, good_ao, \
        min_ao, snp_min_af, snp_min_bq, snp_min_ao, ins_min_af, del_min_af,  \
        del_merge_min_af, ins_merge_min_af, merge_r = candidate_record
    thread_logger = logging.getLogger(
        "{} ({})".format(filter_candidates.__name__, multiprocessing.current_process().name))
    try:
        thread_logger.info(
            "---------------------Filter Candidates---------------------")

        if dbsnp:
            if not dbsnp.endswith("vcf.gz"):
                thread_logger.error("Aborting!")
                raise Exception(
                    "The dbSNP file should be a tabix indexed file with .vcf.gz format")
            if not os.path.exists(dbsnp + ".tbi"):
                thread_logger.error("Aborting!")
                raise Exception(
                    "The dbSNP file should be a tabix indexed file with .vcf.gz format. No {}.tbi file exists.".format(dbsnp))

        records = {}
        with open(candidates_vcf) as v_f:
            for line in skip_empty(v_f):
                if len(line.strip().split()) != 10:
                    raise RuntimeError(
                        "Bad VCF line (<10 fields): {}".format(line))
                chrom, pos, _, ref, alt, _, _, info_, _, info = line.strip().split()
                pos = int(pos)
                loc = "{}.{}".format(chrom, pos)
                dp, ro, ao = list(map(int, info.split(":")[1:4]))
                info_dict = dict(map(lambda x: x.split(
                    "="), filter(None, info_.split(";"))))
                mq_ = safe_read_info_dict(info_dict, "MQ", int, -100)
                bq_ = safe_read_info_dict(info_dict, "BQ", int, -100)
                nm_ = safe_read_info_dict(info_dict, "NM", int, -100)
                as_ = safe_read_info_dict(info_dict, "AS", int, -100)
                xs_ = safe_read_info_dict(info_dict, "XS", int, -100)
                pr_ = safe_read_info_dict(info_dict, "PR", int, -100)
                cl_ = safe_read_info_dict(info_dict, "CL", int, -100)
                st_ = safe_read_info_dict(info_dict, "ST", str, "-100,-100")
                ls_ = safe_read_info_dict(info_dict, "LS", int, -100)
                rs_ = safe_read_info_dict(info_dict, "RS", int, -100)

                if ao < min(ro, min_ao):
                    continue
                if dp >= max_dp:
                    continue

                if loc not in records:
                    records[loc] = []
                if ref == "N" or "\t".join(line.split()[0:5]) \
                        not in map(lambda x: "\t".join(x[-1].split()[0:5]), records[loc]):
                    records[loc].append([chrom, pos, ref, alt, dp, ro, ao, mq_, bq_, st_, ls_, rs_,
                                         nm_, as_, xs_, pr_, cl_, line])
                elif "\t".join(line.split()[0:5]) \
                        in map(lambda x: "\t".join(x[-1].split()[0:5]), records[loc]):
                    for i, x in enumerate(records[loc]):
                        if "\t".join(line.split()[0:5]) == "\t".join(x[-1].split()[0:5]) \
                                and ao / float(ro + 0.0001) > x[6] / float(x[5] + 0.0001):
                            records[loc][i] = [chrom, pos, ref, alt, dp, ro, ao, mq_, bq_, st_, ls_,
                                               rs_, nm_, as_, xs_, pr_, cl_, line]
                            break
        fasta_file = pysam.Fastafile(reference)
        good_records = []
        dels = []
        for loc, rs in sorted(records.items(), key=lambda x: x[1][0:2]) + \
                [["", [["", 0, "", "", 0, 0, 0, ""]]]]:
            ins = list(filter(lambda x: x[2] == "N", rs))
            if len(ins) > 1:
                # emit ins
                afs = list(map(lambda x: x[6] / float(x[5] + x[6]), ins))
                max_af = max(afs)
                ins = list(filter(lambda x: x[6] / float(x[5] +
                                                         x[6]) >= (max_af * merge_r), ins))
                chrom, pos, ref = ins[0][0:3]
                dp = max(map(lambda x: x[4], ins))
                ro = max(map(lambda x: x[5], ins))
                ao = max(map(lambda x: x[6], ins))
                mq_ = max(map(lambda x: x[7], ins))
                bq_ = max(map(lambda x: x[8], ins))
                st_ = "{},{}".format(max(map(lambda x: int(x[9].split(",")[0]), ins)),
                                     max(map(lambda x: int(x[9].split(",")[1]), ins)))
                ls_ = max(map(lambda x: x[10], ins))
                rs_ = max(map(lambda x: x[11], ins))
                nm_ = max(map(lambda x: x[12], ins))
                as_ = max(map(lambda x: x[13], ins))
                xs_ = max(map(lambda x: x[14], ins))
                pr_ = max(map(lambda x: x[15], ins))
                cl_ = max(map(lambda x: x[16], ins))
                alt = "".join(map(lambda x: x[3], ins))
                if (max_af >= ins_merge_min_af) or (ao >= good_ao):
                    ins = [[chrom, pos, ref, alt, dp, ro, ao, mq_, bq_, st_, ls_,
                            rs_, nm_, as_, xs_, pr_, cl_]]
                else:
                    ins = []
            elif len(ins) == 1:
                # emit 1-base ins
                dp, ro, ao = ins[0][4:7]
                if (ao / float(ro + ao) < (ins_min_af) and ao < good_ao) or dp <= 5:
                    ins = []
                else:
                    ins = [ins[0][:-1]]
            good_records.extend(ins)
            if dels and (ins or list(filter(lambda x: x[3] != "N" and x[2] != "N", rs))):
                # emit del
                if len(dels) == 1:
                    ro = dels[0][5]
                    ao = dels[0][6]
                    chrom, pos, ref = dels[0][0:3]
                    if ao / float(ro + ao) >= ((del_min_af)) or ao >= good_ao:
                        good_records.extend(dels)

                else:
                    afs = list(map(lambda x: x[6] / float(x[5] + x[6]), dels))
                    max_af = max(afs)
                    merge_r_thr = merge_r * max_af
                    dels = list(filter(
                        lambda x: x[6] / float(x[5] + x[6]) >= merge_r_thr, dels))
                    chrom, pos = dels[0][0:2]
                    dp = max(map(lambda x: x[4], dels))
                    ro = max(map(lambda x: x[5], dels))
                    ao = max(map(lambda x: x[6], dels))
                    mq_ = max(map(lambda x: x[7], dels))
                    bq_ = max(map(lambda x: x[8], dels))
                    st_ = "{},{}".format(max(map(lambda x: int(x[9].split(",")[0]), dels)),
                                         max(map(lambda x: int(x[9].split(",")[1]), dels)))
                    ls_ = max(map(lambda x: x[10], dels))
                    rs_ = max(map(lambda x: x[11], dels))
                    nm_ = max(map(lambda x: x[12], dels))
                    as_ = max(map(lambda x: x[13], dels))
                    xs_ = max(map(lambda x: x[14], dels))
                    pr_ = max(map(lambda x: x[15], dels))
                    cl_ = max(map(lambda x: x[16], dels))
                    ref = "".join(map(lambda x: x[2], dels))
                    alt = "N"
                    good_records.append([chrom, pos, ref, alt, dp, ro, ao, mq_, bq_, st_, ls_,
                                         rs_, nm_, as_, xs_, pr_, cl_])
                dels = []
            if not loc:
                continue

            for record in rs:
                dp = record[4]
                if dp <= min_dp:
                    continue
                if dp >= max_dp:
                    continue
                ro, ao = record[5:7]
                if record[2] != "N" and record[3] != "N" and record[2] != record[3]:
                    bq = record[8]
                    if (ao / float(ro + ao) >= (snp_min_af) or ao >= snp_min_ao) and bq >= snp_min_bq:
                        # emit SNP
                        good_records.append(record[:-1])
                elif record[2] != "N" and record[3] == "N":
                    if ao / float(ro + ao) >= (del_merge_min_af) or ao >= good_ao:
                        chrom, pos = record[0:2]
                        if dels and pos - dels[-1][1] != 1:
                            # emit del
                            if len(dels) == 1:
                                ro = dels[0][5]
                                ao = dels[0][6]
                                chrom, pos, ref = dels[0][0:3]
                                pos = int(pos)
                                if ao / float(ro + ao) >= ((del_min_af)):
                                    good_records.extend(dels)
                            else:
                                afs = list(map(lambda x: x[6] /
                                               float(x[5] + x[6]), dels))
                                max_af = max(afs)
                                merge_r_thr = merge_r * max_af
                                dels = list(filter(
                                    lambda x: x[6] / float(x[5] + x[6]) >= merge_r_thr, dels))
                                chrom, pos = dels[0][0:2]
                                dp = max(map(lambda x: x[4], dels))
                                ro = max(map(lambda x: x[5], dels))
                                ao = max(map(lambda x: x[6], dels))
                                mq_ = max(map(lambda x: x[7], dels))
                                bq_ = max(map(lambda x: x[8], dels))
                                st_ = "{},{}".format(max(map(lambda x: int(x[9].split(",")[0]),
                                                             dels)),
                                                     max(map(lambda x: int(x[9].split(",")[1]),
                                                             dels)))
                                ls_ = max(map(lambda x: x[10], dels))
                                rs_ = max(map(lambda x: x[11], dels))
                                nm_ = max(map(lambda x: x[12], dels))
                                as_ = max(map(lambda x: x[13], dels))
                                xs_ = max(map(lambda x: x[14], dels))
                                pr_ = max(map(lambda x: x[15], dels))
                                cl_ = max(map(lambda x: x[16], dels))
                                ref = "".join(map(lambda x: x[2], dels))
                                alt = "N"
                                good_records.append([chrom, pos, ref, alt, dp, ro, ao, mq_, bq_,
                                                     st_, ls_, rs_, nm_, as_, xs_, pr_, cl_])
                            dels = []
                        # accumulate dels
                        dels.append(record[:-1])

        final_records = []
        dels = []
        for i, record in enumerate(good_records):
            chrom, pos, ref, alt, dp, ro, ao, mq_, bq_, st_, ls_, rs_, nm_, as_, xs_, pr_, cl_ = record
            ref = ref.upper()
            alt = alt.upper()
            info_str = ""
            if st_ != "-100,-100":
                info_str += ";ST={}".format(st_)
            if ls_ != -100:
                info_str += ";LS={}".format(ls_)
            if rs_ != -100:
                info_str += ";RS={}".format(rs_)
            if nm_ != -100:
                info_str += ";NM={}".format(nm_)
            if as_ != -100:
                info_str += ";AS={}".format(as_)
            if xs_ != -100:
                info_str += ";XS={}".format(xs_)
            if pr_ != -100:
                info_str += ";PR={}".format(pr_)
            if cl_ != -100:
                info_str += ";CL={}".format(cl_)
            if mq_ != -100:
                info_str += ";MQ={}".format(mq_)
            if bq_ != -100:
                info_str += ";BQ={}".format(bq_)

            af = np.round(ao / float(ao + ro), 4)
            info_str += ";AF={}".format(af)
            if ref != "N" and alt != "N":
                line = "\t".join([chrom, str(pos), ".", ref, alt, "100", ".",
                                  "DP={};RO={};AO={}".format(
                                      dp, ro, ao) + info_str,
                                  "GT:DP:RO:AO:AF", "0/1:{}:{}:{}:{}".format(dp, ro, ao, af)])
                final_records.append([chrom, pos, ref, alt, line])
            elif alt == "N":
                ref = fasta_file.fetch(
                    chrom, pos - 2, pos + len(ref) - 1).upper()
                alt = fasta_file.fetch(chrom, pos - 2, pos - 1).upper()
                line = "\t".join([chrom, str(pos - 1), ".", ref, alt, "100", ".",
                                  "DP={};RO={};AO={}".format(
                                      dp, ro, ao) + info_str,
                                  "GT:DP:RO:AO:AF", "0/1:{}:{}:{}:{}".format(dp, ro, ao, af)])
                final_records.append([chrom, pos - 1, ref, alt, line])
            elif ref == "N":
                ref = fasta_file.fetch(chrom, pos - 2, pos - 1).upper()
                alt = ref + alt
                line = "\t".join([chrom, str(pos - 1), ".", ref, alt, "100", ".",
                                  "DP={};RO={};AO={}".format(
                                      dp, ro, ao) + info_str,
                                  "GT:DP:RO:AO:AF", "0/1:{}:{}:{}:{}".format(dp, ro, ao, af)])
                final_records.append([chrom, pos - 1, ref, alt, line])
        final_records = sorted(final_records, key=lambda x: x[0:2])
        if dbsnp:
            dbsnp_tb = pysam.TabixFile(dbsnp)
        with open(filtered_candidates_vcf, "w") as o_f:
            o_f.write("{}\n".format(VCF_HEADER))
            o_f.write(
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
            for record in final_records:
                if dbsnp:
                    chrom, pos, ref, alt = record[0:4]
                    var_id = "-".join(map(str,[chrom, pos, ref, alt]))
                    region = "{}:{}-{}".format(chrom, pos, pos + 1)
                    dbsnp_vars = []
                    for x in dbsnp_tb.fetch(region=region):
                        chrom_, pos_, _, ref_, alts_ = x.strip().split("\t")[
                            0:5]
                        for alt_ in alts_.split(","):
                            dbsnp_var_id = "-".join(map(str,[chrom_, pos_, ref_, alt_]))
                            dbsnp_vars.append(dbsnp_var_id)
                    if var_id in dbsnp_vars:
                        continue
                o_f.write(record[-1] + "\n")
        return filtered_candidates_vcf

    except Exception as ex:
        thread_logger.error(traceback.format_exc())
        thread_logger.error(ex)
        return None


if __name__ == '__main__':
    FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    logging.basicConfig(level=logging.INFO, format=FORMAT)
    logger = logging.getLogger(__name__)

    parser = argparse.ArgumentParser(description='filter candidates vcf')
    parser.add_argument('--candidates_vcf', type=str, help='raw candidates vcf',
                        required=True)
    parser.add_argument('--filtered_candidates_vcf', type=str, help='filtered candidates vcf',
                        required=True)
    parser.add_argument('--reference', type=str, help='reference fasta filename',
                        required=True)
    parser.add_argument('--dbsnp_to_filter', type=str,
                        help='dbsnp vcf.gz (will be used to filter candidate variants)', default=None)
    parser.add_argument('--good_ao', type=float, help='good alternate count (ignores maf)',
                        default=10)
    parser.add_argument('--min_ao', type=float,
                        help='min alternate count', default=1)
    parser.add_argument('--snp_min_af', type=float,
                        help='SNP min allele freq', default=0.05)
    parser.add_argument('--snp_min_bq', type=float,
                        help='SNP min base quality', default=10)
    parser.add_argument('--snp_min_ao', type=float,
                        help='SNP min alternate count for low AF candidates', default=3)
    parser.add_argument('--ins_min_af', type=float,
                        help='INS min allele freq', default=0.05)
    parser.add_argument('--del_min_af', type=float,
                        help='DEL min allele freq', default=0.05)
    parser.add_argument('--del_merge_min_af', type=float,
                        help='min allele freq for merging DELs', default=0)
    parser.add_argument('--ins_merge_min_af', type=float,
                        help='min allele freq for merging INSs', default=0)
    parser.add_argument('--merge_r', type=float,
                        help='merge af ratio to the max af for merging adjacent variants',
                        default=0.5)
    parser.add_argument('--min_dp', type=float, help='min depth', default=5)
    parser.add_argument('--max_dp', type=float,
                        help='max depth', default=100000)

    args = parser.parse_args()
    logger.info(args)

    try:
        output = filter_candidates((args.candidates_vcf, args.filtered_candidates_vcf,
                                    args.reference, args.dbsnp_to_filter, args.min_dp, args.max_dp,
                                    args.good_ao, args.min_ao,
                                    args.snp_min_af, args.snp_min_bq, args.snp_min_ao,
                                    args.ins_min_af, args.del_min_af,
                                    args.del_merge_min_af, args.ins_merge_min_af, args.merge_r))
        if output is None:
            raise Exception("filter_candidates failed!")
    except Exception as e:
        logger.error(traceback.format_exc())
        logger.error("Aborting!")
        logger.error(
            "filter_candidates.py failure on arguments: {}".format(args))
        raise e
