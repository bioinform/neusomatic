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

max_indel = 100


def extract_del(record):
    logger = logging.getLogger(extract_del.__name__)
    dels = []
    pos = record.pos
    cigartuples = record.cigartuples
    first_sc = 1 if cigartuples[0][0] in [
        CIGAR_SOFTCLIP, CIGAR_HARDCLIP] else 0
    last_sc = 1 if cigartuples[-1][0] in [CIGAR_SOFTCLIP,
                                          CIGAR_HARDCLIP] else 0
    for i, (C, L) in enumerate(cigartuples):
        if C in [CIGAR_SOFTCLIP, CIGAR_HARDCLIP, CIGAR_INS]:
            continue
        if C == CIGAR_DEL:
            if i > first_sc and i < len(cigartuples) - 1 - last_sc:
                L_ = min(L, max_indel)
                dels.append([record.reference_name, pos, pos + L_])
        pos += L

    return dels


def extract_ins(record):
    logger = logging.getLogger(extract_ins.__name__)
    inss = []
    pos = record.pos
    seq_pos = 0
    cigartuples = record.cigartuples
    first_sc = 1 if cigartuples[0][0] in [
        CIGAR_SOFTCLIP, CIGAR_HARDCLIP] else 0
    last_sc = 1 if cigartuples[-1][0] in [CIGAR_SOFTCLIP,
                                          CIGAR_HARDCLIP] else 0
    for i, (C, L) in enumerate(cigartuples):
        if C == CIGAR_SOFTCLIP:
            seq_pos += L
            continue
        elif C == CIGAR_HARDCLIP:
            continue
        if C == CIGAR_INS:
            if not record.seq[seq_pos:seq_pos + L]:
                logger.info([str(record).split("\t"), seq_pos,
                             L, len(record.seq), len(record.seq)])
            if i > first_sc and i < len(cigartuples) - 1 - last_sc:
                L_ = min(L, max_indel)
                inss.append([record.reference_name, pos, pos + 1,
                             record.seq[seq_pos:seq_pos + L_]])
            seq_pos += L
        else:
            if C != CIGAR_DEL:
                seq_pos += L
            pos += L
    return inss


def find_vtype(ref, alt):
    if len(alt) < len(ref):
        vtype = "DEL"
    elif len(alt) > len(ref):
        vtype = "INS"
    else:
        vtype = "SNP"
    return vtype


def push_left_var(ref_fasta, chrom, pos, ref, alt):
    logger = logging.getLogger(push_left_var.__name__)
    pos = int(pos)
    while ref[-1] == alt[-1] and pos > 1:
        prev_base = ref_fasta.fetch(chrom, pos - 2, pos - 1)
        pos -= 1
        ref = prev_base + ref[:-1]
        alt = prev_base + alt[:-1]
    while ref[0] == alt[0] and len(ref) == len(alt) and len(ref) > 1:
        pos += 1
        ref = ref[1:]
        alt = alt[1:]
    return [chrom, pos, ref, alt]


class Variant:

    def __init__(self, chrom, pos, ref, alt, gt, score, cnt, vtype):
        self.chrom = chrom
        self.pos = int(pos)
        self.ref = ref
        self.alt = alt
        self.gt = gt
        self.score = float(score)
        self.cnt = float(cnt) if cnt is not None else None
        self.vtype = vtype
        self.len = abs(len(alt) - len(ref))
        self.processed = False

    def push_left(self, ref_fasta):
        _, self.pos, self.ref, self.alt = push_left_var(
            ref_fasta, self.chrom, self.pos, self.ref, self.alt)

    def var_str(self):
        return "-".join(map(str, [self.chrom, self.pos, self.ref, self.alt, self.vtype]))

    def var_pos_vt_str(self):
        return "-".join(map(str, [self.chrom, self.pos, self.vtype]))

    def var_gt_str(self):
        return "-".join(map(str, [self.chrom, self.pos, self.ref, self.alt, self.gt, self.vtype]))

    def __str__(self):
        return "-".join(map(str, [self.chrom, self.pos, self.ref, self.alt, self.gt,
                                  self.score, self.cnt, self.vtype]))


def resolve_group(ref_fasta, variants, vars_count):
    logger = logging.getLogger(resolve_group.__name__)
    chrom = variants[0].chrom
    vars_count_ = {}
    for var_str in vars_count:
        pos, ref, alt, vtype = var_str.split("-")[-4:]
        pos = int(pos)
        v = Variant(chrom, pos, ref, alt, "0/0",
                    0, vars_count[var_str], vtype)
        v.push_left(ref_fasta)
        s = v.var_str()
        if s not in vars_count_:
            vars_count_[s] = 0
        vars_count_[s] += vars_count[var_str]
    vars_count = vars_count_

    group_vars = {}
    processed = []
    for v in variants:
        if v.pos not in group_vars:
            group_vars[v.pos] = []
        var_str = v.var_str()
        if var_str not in vars_count:
            vars_count[var_str] = 0
        cnt = vars_count[var_str]
        v.cnt = cnt
        processed.append(var_str)
        group_vars[v.pos].append(v)

    for var_str in vars_count:
        if var_str not in processed:
            pos, ref, alt, vtype = var_str.split("-")[-4:]
            pos = int(pos)
            if pos not in group_vars:
                group_vars[pos] = []
            v = Variant(chrom, pos, ref, alt, "0/0",
                        0, vars_count[var_str], vtype)
            group_vars[pos].append(v)
    for pos in group_vars:
        var_ = {}
        for v in group_vars[pos]:
            var_id = v.var_gt_str()
            if var_id not in var_:
                var_[var_id] = []
            var_[var_id].append(v)
        group_vars[pos] = []
        for var_id in var_:
            group_vars[pos].append(sorted(var_[var_id], key=lambda x: x.score, reverse=True
                                          )[0])

    out_variants_ = []
    max_target = [v.cnt for pos in group_vars for v in group_vars[
        pos] if v.score > 0 or v.len >= 3]
    if len(max_target) == 0:
        # logger.info(
        #     "No non-zero COUNT with non-zero SCORE: {}".format(list(str(x) for x in group_vars[pos])))
        return []

    # logger.info(list([pos, [str(y) for y in x]]
    #                  for pos, x in group_vars.items()))

    max_count = max(max_target)
    for pos in group_vars.keys():
        if max(map(lambda x: x.cnt, group_vars[pos])
               ) < 0.2 * max_count:
            continue
        mx = max(map(lambda x: x.cnt, group_vars[pos]))
        gts = [x.gt for x in group_vars[pos]]
        gts = set([x.gt for x in group_vars[pos]
                   if x.gt != "0/0" or x.len >= 3])
        if len(gts) == 0:
            continue
        if len(gts) > 1:
            gts_count = {"0/1": 0, "0/0": 0}
            gts_score = {"0/1": 0, "0/0": 0}
            for x in group_vars[pos]:
                if (x.gt != "0/0" or x.len >= 3) and x.cnt >= 0.4 * mx:
                    gts_count[x.gt] += x.cnt
                    gts_score[x.gt] += x.score
            priority = {"0/1": 2, "0/0": 1}
            sorted_gts = sorted(gts_count.keys(), key=lambda x: [
                                gts_count[x], gts_score[x],
                                priority[x]], reverse=True)
            gt = sorted_gts[0]
        else:
            gt = list(gts)[0]

        all_vars = sorted(group_vars[pos], key=lambda x: [
                          x.cnt, x.score, x.gt != "0/0"], reverse=True)
        vtypes = set([x.vtype for x in group_vars[pos]
                      if (x.gt != "0/0" or x.len >= 3) and x.cnt >= 0.4 * mx])
        if not vtypes:
            vtypes = set([x.vtype for x in group_vars[pos]
                          if (x.gt != "0/0" or x.len >= 3)])
        all_vars = list(
            filter(lambda x: x.vtype in vtypes, all_vars))
        if not all_vars:
            logger.info(
                "No vars: {}".format(list(str(x) for x in group_vars[pos])))
            logger.info(
                "No vars: {}".format([[list(str(x) for x in group_vars[pos_])]for pos_ in group_vars]))
            raise Exception
        score = max([v.score for v in all_vars])
        if gt == "0/0":
            nz_vars = [x for x in all_vars if x.gt !=
                       "0/0" and x.vtype == all_vars[0].vtype]
            if nz_vars:
                nz_vars = sorted(nz_vars, key=lambda x: [
                                 x.score], reverse=True)[0]
                gt = nz_vars.gt
        v = all_vars[0]
        out_variants_.append(
            [v.chrom, v.pos, v.ref, v.alt, gt, score, v.cnt])

    if len(out_variants_) == 1 and out_variants_[0][4] == "0/0" and abs(len(out_variants_[0][2]) - len(out_variants_[0][3])) >= 3:
        chrom_, pos_, ref_, alt_, gt_, score_, cnt_ = out_variants_[0]
        vtype = find_vtype(ref_, alt_)
        resolve_candids = []
        for pos in group_vars.keys():
            for y in group_vars[pos]:
                if y.vtype == vtype and y.gt != "0/0":
                    resolve_candids.append(y)
        if resolve_candids:
            resolve_candids = sorted(resolve_candids, key=lambda x: [
                x.score], reverse=True)[0]
            out_variants_ = [[chrom_, pos_, ref_, alt_,
                              resolve_candids.gt, resolve_candids.score, cnt_]]

    if len(out_variants_) > 1 and "0/0" in [x[4] for x in out_variants_]:
        nz_vars = [x for x in out_variants_ if x[4] != "0/0"]
        if nz_vars:
            nz_vtypes = [find_vtype(x[2], x[3]) for x in nz_vars]
            out_variants__ = []
            for x in out_variants_:
                if x[4] != "0/0":
                    out_variants__.append(x)
                else:
                    vtype = find_vtype(x[2], x[3])
                    if vtype not in nz_vtypes:
                        out_variants__.append(x)
            out_variants_ = out_variants__

    vars_gt = {}
    for chrom_, pos_, ref_, alt_, gt_, score_, cnt_ in out_variants_:
        if gt_ not in vars_gt:
            vars_gt[gt_] = []
        vars_gt[gt_].append(
            Variant(chrom_, pos_, ref_, alt_, gt_, score_, cnt_, ""))
    vars_gt = {gt_: sorted(vars_gt[gt_], key=lambda x: [
                           x.cnt, x.score], reverse=True) for gt_ in vars_gt}
    out_variants_ = []
    for gt_ in vars_gt:
        v0 = vars_gt[gt_][0]
        good_vs = [v0]
        for v in vars_gt[gt_][1:]:
            keep = True
            for g_v in good_vs:
                if min(v.pos + len(v.ref), g_v.pos + len(g_v.ref)) > max(v.pos, g_v.pos):
                    keep = False
                    break
            if keep:
                good_vs.append(v)
        for v in good_vs:
            out_variants_.append(
                [v.chrom, v.pos, v.ref, v.alt, v.gt, v.score])
    return out_variants_


def find_resolved_variants(input_record):
    chrom, start, end, variants, input_bam, filter_duplicate, reference = input_record
    thread_logger = logging.getLogger(
        "{} ({})".format(find_resolved_variants.__name__, multiprocessing.current_process().name))
    try:
        ref_fasta = pysam.FastaFile(reference)
        variants_ = []
        for x in variants:
            pos = int(x[1])
            ref = x[3]
            alt = x[4]
            gt = x[9].split(":")[0]
            score = x[5]
            vtype = x[-1]
            v = Variant(chrom, pos, ref, alt, gt, score, None, vtype)
            v.push_left(ref_fasta)
            variants_.append(v)
        variants = variants_
        start, end = list(map(int, [start, end]))
        region = [chrom, start, end]
        vartypes = list(map(lambda x: x.vtype, variants))
        scores = list(map(lambda x: x.score, variants))
        dels = []
        inss = []
        snps = []
        vars_count = {}
        with pysam.AlignmentFile(input_bam) as samfile:
            cov = 0
            dels_ = []
            inss_ = []
            snps_ = []
            for record in samfile.fetch(chrom, start, end):
                if record.is_unmapped:
                    continue
                if record.seq is None:
                    continue
                if not record.is_duplicate or not filter_duplicate:
                    cov += 1
                    if record.cigarstring and "D" in record.cigarstring:
                        dels_.extend(extract_del(record))
                    if record.cigarstring and "I" in record.cigarstring:
                        inss_.extend(extract_ins(record))
                    aligned_pairs = np.array(
                        record.get_aligned_pairs(matches_only=True))
                    near_pos = np.where((start <= aligned_pairs[:, 1]) & (
                        aligned_pairs[:, 1] <= end))[0]
                    if len(near_pos) != 0:
                        for pos_i in near_pos:
                            seq_pos, ref_pos = aligned_pairs[pos_i, :]
                            if seq_pos is not None:
                                ref_snp = ref_fasta.fetch(
                                    chrom, ref_pos, ref_pos + 1).upper()
                                alt_snp = record.seq[seq_pos]
                                if alt_snp != ref_snp:
                                    snps_.append(
                                        [chrom, ref_pos + 1, ref_snp, alt_snp])

            dels.extend([x + [1.0 / (cov)] for x in dels_])
            inss.extend([x + [1.0 / (cov)] for x in inss_])
            snps.extend([x + [1.0 / (cov)] for x in snps_])

        dels = list(filter(lambda x: (
            start <= x[1] <= end) or start <= x[2] <= end, dels))
        if dels:
            del_strs = []
            cnt_ = {}
            for x in dels:
                chrom, st, en, cnt = x
                del_str = "---".join(map(str, [chrom, st, en]))
                if del_str not in cnt_:
                    cnt_[del_str] = 0
                cnt_[del_str] += cnt
                del_strs.append(del_str)

            uniq_dels = list(set(del_strs))
            for del_ in uniq_dels:
                st, en = map(int, del_.split("---")[1:3])
                del_str = "-".join(list(map(str, [chrom, int(st), ref_fasta.fetch(chrom, st - 1, en).upper(),
                                                  ref_fasta.fetch(chrom, st - 1, st).upper(), "DEL"])))
                vars_count[del_str] = np.round(cnt_[del_], 4)
        inss = list(filter(lambda x: (
            start <= x[1] <= end) or start <= x[2] <= end, inss))
        if inss:
            cnt_ = {}
            ins_strs = []
            for x in inss:
                chrom, st, en, bases, cnt = x
                ins_str = "---".join(map(str, [chrom, st, en, bases]))
                if ins_str not in cnt_:
                    cnt_[ins_str] = 0
                cnt_[ins_str] += cnt
                ins_strs.append(ins_str)
            uniq_inss = list(set(ins_strs))
            for ins_ in uniq_inss:
                st, en, bases = ins_.split("---")[1:4]
                st, en = map(int, [st, en])
                ins_str = "-".join(list(map(str, [chrom, int(st), ref_fasta.fetch(chrom, st - 1, st).upper(),
                                                  ref_fasta.fetch(chrom, st - 1, st).upper() + bases, "INS"])))
                vars_count[ins_str] = np.round(cnt_[ins_], 4)

        if snps:
            cnt_ = {}
            snp_strs = []
            for x in snps:
                chrom, st, ref_, alt_, cnt = x
                snp_str = "---".join(map(str, [chrom, st, ref_, alt_]))
                if snp_str not in cnt_:
                    cnt_[snp_str] = 0
                cnt_[snp_str] += cnt
                snp_strs.append(snp_str)
            uniq_snps = list(set(snp_strs))
            for snp_ in uniq_snps:
                st, ref_, alt_ = snp_.split("---")[1:4]
                snp_str = "-".join(list(map(str, [chrom, st, ref_,
                                                  alt_, "SNP"])))
                vars_count[snp_str] = np.round(cnt_[snp_], 4)

        out_variants_ = resolve_group(ref_fasta, variants, vars_count)
        return out_variants_

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
            vartype = find_vtype(fields[3], fields[4])
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

    if num_threads > 1:
        try:
            n_per_bacth = min(10 * num_threads, len(map_args))
            out_variants_list = []
            i = 0
            while i < len(map_args):
                pool = multiprocessing.Pool(num_threads)
                batch_i_s = i
                batch_i_e = min(i + n_per_bacth, len(map_args))
                out_variants_list.extend(pool.map_async(
                    find_resolved_variants, map_args[batch_i_s:batch_i_e]).get())
                i = batch_i_e
                pool.close()
        except Exception as inst:
            logger.error(inst)
            pool.close()
            traceback.print_exc()
            raise Exception
    else:
        out_variants_list = [find_resolved_variants(w) for w in map_args]

    out_variants = [x for xs in out_variants_list for x in xs]
    chroms_order = get_chromosomes_order(bam=input_bam)

    out_variants = sorted(out_variants, key=lambda x: [
                          chroms_order[x[0]], int(x[1])])
    with open(resolved_vcf, "w") as o_f:
        o_f.write("{}\n".format(VCF_HEADER))
        o_f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
        done_id = set([])
        for chrom, pos, ref, alt, gt, phred_score in out_variants:
            if ref != alt:
                id_ = "-".join(list(map(str, [chrom, pos, ref, alt])))
                if id_ in done_id:
                    continue
                done_id.add(id_)
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
