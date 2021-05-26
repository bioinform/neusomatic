#!/usr/bin/env python
#-------------------------------------------------------------------------
# long_read_indelrealign.py
# Resolve variants (e.g. exact INDEL sequences) for long high-error rate
# sequences. Target variants are identified by
# 'extract_postprocess_targets.py'.
#-------------------------------------------------------------------------
import os
import traceback
import argparse
import shutil
import multiprocessing
import re
import fileinput
import logging
from random import shuffle
import csv
import functools
import tempfile

import numpy as np
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import DNAAlphabet

from utils import run_shell_command, get_chromosomes_order, run_bedtools_cmd, write_tsv_file, read_tsv_file, bedtools_sort, bedtools_merge, get_tmp_file, skip_empty

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


def nuc_to_num_convert(nuc):
    if nuc.upper() not in NUC_to_NUM.keys():
        nuc = "-"
    return NUC_to_NUM[nuc.upper()]


class Region:

    def __init__(self, fields, pad, len_chr):
        self.chrom = str(fields[0])
        self.start = max(2, int(fields[1]) - pad)
        self.end = min(int(fields[2]) + pad, len_chr - 2)
        self.pad = pad

    def span(self):
        return abs(self.end - self.start)

    def __str__(self):
        return "{}-{}-{}".format(self.chrom, self.start, self.end)


class Read_Info:

    def __init__(self, fields):
        self.query_name = fields[0]
        self.pos = int(fields[1])
        self.cigarstring = fields[2]
        self.start_idx = int(fields[3])
        self.end_idx = int(fields[4])
        self.del_start = int(fields[5])
        self.del_end = int(fields[6])
        self.pos_start = int(fields[7])
        self.pos_end = int(fields[8])
        self.end_ra_pos = 0


class Realign_Read:

    def __init__(self, query_name, chrom, pos, cigarstring):
        self.query_name = query_name
        self.chrom = chrom
        self.pos = pos
        self.cigarstring = cigarstring
        self.realignments = []
        self.end_ra_pos = 0

    def add_realignment(self, region_start, region_end, start_idx, end_idx, del_start, del_end,
                        pos_start, pos_end, new_cigar, excess_start, excess_end):
        logger = logging.getLogger(Realign_Read.add_realignment.__name__)
        pos_start = int(pos_start)
        if pos_start > self.end_ra_pos:
            self.realignments.append([int(region_start), int(region_end), int(start_idx),
                                      int(end_idx),
                                      int(del_start), int(del_end), int(
                                          pos_start), int(pos_end),
                                      new_cigar, int(excess_start), int(excess_end)])
            self.end_ra_pos = int(pos_end)

    def fix_record(self, record, ref_seq):
        logger = logging.getLogger(Realign_Read.fix_record.__name__)
        self.realignments = sorted(self.realignments, key=lambda x: x[0])
        cigartuples = record.cigartuples
        start_hc = cigartuples[0][1] if cigartuples[
            0][0] == CIGAR_HARDCLIP else 0
        end_hc = cigartuples[-1][1] if cigartuples[-1][0] == CIGAR_HARDCLIP else 0
        if start_hc:
            cigartuples = cigartuples[1:]
        if end_hc:
            cigartuples = cigartuples[:-1]

        new_cigartuples_list = []
        if start_hc:
            new_cigartuples_list.append([[CIGAR_HARDCLIP, start_hc]])
        try:
            assert self.realignments
        except:
            logger.error("Realignments are empty for {} at {}:{}".format(
                record, self.chrom, self.pos))
            raise Exception
        bias = 0
        for region_start, region_end, start_idx, end_idx, del_start, del_end, \
                pos_start, pos_end, new_cigar, excess_start, excess_end in self.realignments:
            c_array = np.array(list(map(lambda x: [x[0], x[1][1] if x[1][0]
                                                   != CIGAR_DEL else 0], enumerate(cigartuples))))
            c_map = np.repeat(c_array[:, 0], c_array[:, 1])

            c_i = c_map[start_idx - bias]
            c_e = c_map[end_idx - bias]
            begin_match = np.nonzero(c_map == c_i)[0][0]
            end_match = np.nonzero(c_map == c_e)[0][-1]

            if excess_start > 0:
                if (c_i > 0 and cigartuples[c_i - 1][0] == CIGAR_DEL and
                        cigartuples[c_i - 1][1] >= (excess_start)):
                    if excess_start > del_start:
                        del_start = excess_start
                else:
                    return record

            if excess_end > 0:
                if (c_e < (len(cigartuples) - 1) and cigartuples[c_e + 1][0] == CIGAR_DEL
                        and cigartuples[c_e + 1][1] >= (excess_end)):
                    if excess_end > del_end:
                        del_end = excess_end
                elif not(c_e == len(c_array) - 1 or cigartuples[c_e + 1][0] == CIGAR_SOFTCLIP):
                    return record

            left_cigartuple = cigartuples[:c_i]
            if del_start == 0:
                if begin_match < (start_idx - bias):
                    left_cigartuple.append(
                        [cigartuples[c_i][0], start_idx - bias - begin_match])
            else:
                try:
                    assert cigartuples[c_i - 1][0] == CIGAR_DEL
                except:
                    logger.error("Expect DEL (c_i) in positon {} at cigartuples {}".format(
                        c_i - 1, cigartuples))
                    raise Exception
                if cigartuples[c_i - 1][1] > del_start:
                    left_cigartuple = left_cigartuple[
                        :-1] + [[cigartuples[c_i - 1][0], cigartuples[c_i - 1][1] - del_start]]
                else:
                    left_cigartuple = left_cigartuple[:-1]

            new_cigartuples_list.append(left_cigartuple)
            new_cigartuples_list.append(list(cigarstring_to_tuple(new_cigar)))
            right_cigartuple = cigartuples[c_e + 1:]
            if del_end == 0:
                if end_match > (end_idx - bias):
                    right_cigartuple = [
                        [cigartuples[c_e][0], end_match - (end_idx - bias)]] + right_cigartuple
            else:
                try:
                    assert cigartuples[c_e + 1][0] == CIGAR_DEL
                except:
                    logger.info("Expect DEL (c_e) in positon {} at cigartuples {}, {}".format(
                        c_e + 1, cigartuples, self.realignments))
                    logger.info(cigartuple_to_string(cigartuples))
                    logger.info([region_start, region_end, start_idx, end_idx, del_start, del_end,
                                 pos_start, pos_end, new_cigar, excess_start, excess_end])
                    raise Exception
                if cigartuples[c_e + 1][1] > del_end:
                    right_cigartuple = [[cigartuples[
                        c_e + 1][0], cigartuples[c_e + 1][1] - del_end]] + right_cigartuple[1:]
                else:
                    right_cigartuple = right_cigartuple[1:]
            cigartuples = right_cigartuple
            bias = end_idx + 1
        new_cigartuples_list.append(right_cigartuple)

        if end_hc:
            new_cigartuples_list.append([[CIGAR_HARDCLIP, end_hc]])

        new_cigartuples = functools.reduce(
            lambda x, y: merge_cigartuples(x, y), new_cigartuples_list)
        if len(new_cigartuples) > 2 and new_cigartuples[-1][0] == CIGAR_SOFTCLIP \
                and new_cigartuples[-2][0] == CIGAR_DEL:
            new_cigartuples = new_cigartuples[:-2] + [new_cigartuples[-1]]
        elif new_cigartuples[-1][0] == CIGAR_DEL:
            new_cigartuples = new_cigartuples[:-1]
        if len(new_cigartuples) > 2 and new_cigartuples[0][0] == CIGAR_SOFTCLIP \
                and new_cigartuples[1][0] == CIGAR_DEL:
            new_cigartuples = [new_cigartuples[0]] + new_cigartuples[2:]
        elif new_cigartuples[0][0] == CIGAR_DEL:
            new_cigartuples = new_cigartuples[1:]

        try:
            assert(sum(map(lambda x: x[1] if x[0] != CIGAR_DEL else 0, new_cigartuples)) ==
                   sum(map(lambda x: x[1] if x[0] != CIGAR_DEL else 0, record.cigartuples)))
        except:
            logger.error("Old and new cigarstrings have different lengthes: {} vs {}".format(
                sum(map(lambda x: x[1] if x[0] !=
                        CIGAR_DEL else 0, new_cigartuples)),
                sum(map(lambda x: x[1] if x[0] != CIGAR_DEL else 0, record.cigartuples))))
            raise Exception

        record.cigarstring = cigartuple_to_string(new_cigartuples)
        NM = find_NM(record, ref_seq)
        record.tags = list(filter(
            lambda x: x[0] != "NM", record.tags)) + [("NM", int(NM))]
        return record


def get_cigar_stat(cigartuple, keys=[]):
    logger = logging.getLogger(get_cigar_stat.__name__)
    if not keys:
        keys = set(map(lambda x: x[0], cigartuple))
    n_key = {}
    for i in keys:
        n_key[i] = sum(map(lambda x: x[1] if x[0] == i else 0, cigartuple))
    return n_key


def find_NM(record, ref_seq):
    logger = logging.getLogger(find_NM.__name__)
    positions = np.array(list(map(lambda x: x if x else -1,
                                  (record.get_reference_positions(full_length=True)))))
    sc_start = (record.cigartuples[0][0] ==
                CIGAR_SOFTCLIP) * record.cigartuples[0][1]
    sc_end = (record.cigartuples[-1][0] ==
              CIGAR_SOFTCLIP) * record.cigartuples[-1][1]
    q_seq = record.seq
    q_seq = q_seq[sc_start:]
    positions = positions[sc_start:]
    if sc_end > 0:
        positions = positions[:-sc_end]
        q_seq = q_seq[:-sc_end]
    non_ins = np.nonzero(positions >= 0)
    non_ins_positions = positions[non_ins]
    mn, mx = min(non_ins_positions), max(non_ins_positions)
    refseq = ref_seq.get_seq(mn, mx + 1)
    ref_array = np.array(list(map(lambda x: NUC_to_NUM[x.upper()], list(refseq))))[
        non_ins_positions - mn]
    q_array = np.array(list(map(lambda x: NUC_to_NUM[x.upper()], list(q_seq))))[
        non_ins]
    cigar_stat = get_cigar_stat(record.cigartuples, [CIGAR_DEL, CIGAR_INS])
    assert ref_array.shape[0] == q_array.shape[0]
    NM = sum(abs(ref_array - q_array) > 0) + \
        cigar_stat[CIGAR_DEL] + cigar_stat[CIGAR_INS]
    return NM


def cigarstring_to_tuple(cigarstring):
    logger = logging.getLogger(cigarstring_to_tuple.__name__)
    return tuple((_CIGAR_OP_DICT[op],
                  int(length)) for length,
                 op in _CIGAR_PATTERN.findall(cigarstring))


def cigartuple_to_string(cigartuples):
    logger = logging.getLogger(cigartuple_to_string.__name__)
    return "".join(map(lambda x: "%d%s" % (x[1], _CIGAR_OPS[x[0]]), cigartuples))


def prepare_fasta(work, region, input_bam, ref_fasta_file, include_ref, split_i, ds, filter_duplicate):
    logger = logging.getLogger(prepare_fasta.__name__)
    in_fasta_file = os.path.join(
        work, region.__str__() + "_split_{}".format(split_i) + "_0.fasta")
    info_file = os.path.join(work, region.__str__() +
                             "_split_{}".format(split_i) + ".txt")
    with pysam.Fastafile(ref_fasta_file) as ref_fasta:
        with open(in_fasta_file, "w") as in_fasta:
            with open(info_file, "w") as info_txt:
                if include_ref:
                    ref_seq = ref_fasta.fetch(
                        region.chrom, region.start, region.end + 1).upper()
                    in_fasta.write(">0\n")
                    in_fasta.write("%s\n" % ref_seq.upper())
                cnt = 1
                with pysam.Samfile(input_bam, "rb") as samfile:
                    for record in samfile.fetch(region.chrom, region.start, region.end + 1):
                        if record.is_unmapped:
                            continue
                        if filter_duplicate and record.is_duplicate:
                            continue
                        if record.is_supplementary and "SA" in dict(record.tags):
                            sas = dict(record.tags)["SA"].split(";")
                            sas = list(filter(None, sas))
                            sas_cigs = list(
                                map(lambda x: x.split(",")[3], sas))
                            if record.cigarstring in sas_cigs:
                                continue
                        positions = np.array(list(map(lambda x: x if x else -1,
                                                      (record.get_reference_positions(
                                                          full_length=True)))))
                        if not record.cigartuples:
                            continue
                        if np.random.rand() > ds:
                            continue
                        sc_start = (record.cigartuples[0][0] ==
                                    CIGAR_SOFTCLIP) * record.cigartuples[0][1]
                        sc_end = (record.cigartuples[-1][0] ==
                                  CIGAR_SOFTCLIP) * record.cigartuples[-1][1]
                        positions = positions[sc_start:]
                        if sc_end > 0:
                            positions = positions[:-sc_end]
                        rstart = max(region.start, record.pos)
                        rend = min(record.aend - 1, region.end)
                        pos_start = positions[positions >= rstart][0]
                        pos_end = positions[
                            (positions <= rend) & (positions >= 0)][-1]
                        del_start = pos_start - rstart
                        del_end = rend - pos_end
                        start_idx = np.nonzero(positions == pos_start)[
                            0][0] + sc_start
                        end_idx = np.nonzero(positions == pos_end)[
                            0][0] + sc_start

                        if max(end_idx - start_idx, pos_end - pos_start) >= (region.span() * 0.75):
                            c_array = np.array(list(map(lambda x: [x[0], x[1][1] if x[1][0]
                                                                   != CIGAR_DEL else 0],
                                                        enumerate(record.cigartuples))))
                            c_map = np.repeat(c_array[:, 0], c_array[:, 1])
                            c_i = c_map[start_idx]
                            c_e = c_map[end_idx]
                            begin_match = np.nonzero(c_map == c_i)[0][0]
                            end_match = np.nonzero(c_map == c_e)[0][-1]
                            my_cigartuples = record.cigartuples[c_i:c_e + 1]
                            positions_ = positions[
                                (start_idx - sc_start):(end_idx - sc_start)]
                            non_ins = np.nonzero(positions_ >= 0)
                            refseq = ref_fasta.fetch(region.chrom, positions_[non_ins][0],
                                                     positions_[non_ins][-1] + 1).upper()
                            q_seq = record.seq[start_idx:end_idx + 1]
                            non_ins_positions = positions_[non_ins]
                            mn, mx = min(non_ins_positions), max(
                                non_ins_positions)
                            ref_array = np.array(list(map(lambda x:
                                                          NUC_to_NUM[
                                                              x.upper()],
                                                          list(refseq))))[non_ins_positions - mn]
                            q_array = np.array(list(map(lambda x: NUC_to_NUM[x.upper()], list(q_seq))))[
                                non_ins]
                            cigar_stat = get_cigar_stat(
                                my_cigartuples, [CIGAR_DEL, CIGAR_INS])
                            assert ref_array.shape[0] == q_array.shape[0]
                            NM_SNP = sum(abs(ref_array - q_array) > 0)
                            NM_INDEL = cigar_stat[
                                CIGAR_DEL] + cigar_stat[CIGAR_INS] + del_start + del_end
                            in_fasta.write(">%s\n" % cnt)
                            in_fasta.write(
                                "%s\n" % record.seq[start_idx:end_idx + 1].upper())
                            info_txt.write("\t".join(map(str, [cnt, record.query_name, record.pos,
                                                               record.cigarstring, start_idx,
                                                               end_idx,
                                                               del_start, del_end, pos_start,
                                                               pos_end, NM_SNP, NM_INDEL])) + "\n")
                            cnt += 1
    return in_fasta_file, info_file


def split_bam_to_chunks(work, region, input_bam, chunk_size=200,
                        chunk_scale=1.5, do_split=False, filter_duplicate=False):
    logger = logging.getLogger(split_bam_to_chunks.__name__)
    records = []
    with pysam.Samfile(input_bam, "rb") as samfile:
        for record in samfile.fetch(region.chrom, region.start, region.end + 1):
            if record.is_unmapped:
                continue
            if filter_duplicate and record.is_duplicate:
                continue
            if record.is_supplementary and "SA" in dict(record.tags):
                sas = dict(record.tags)["SA"].split(";")
                sas = list(filter(None, sas))
                sas_cigs = list(map(lambda x: x.split(",")[3], sas))
                if record.cigarstring in sas_cigs:
                    continue

            positions = np.array(
                list(map(lambda x: x if x else -1, (record.get_reference_positions(full_length=True)))))
            if not record.cigartuples:
                continue

            sc_start = (record.cigartuples[0][0] ==
                        CIGAR_SOFTCLIP) * record.cigartuples[0][1]
            sc_end = (record.cigartuples[-1][0] ==
                      CIGAR_SOFTCLIP) * record.cigartuples[-1][1]
            positions = positions[sc_start:]
            if sc_end > 0:
                positions = positions[:-sc_end]
            rstart = max(region.start, record.pos)
            rend = min(record.aend - 1, region.end)
            pos_start = positions[positions >= rstart][0]
            pos_end = positions[(positions <= rend) & (positions >= 0)][-1]
            start_idx = np.nonzero(positions == pos_start)[0][0] + sc_start
            end_idx = np.nonzero(positions == pos_end)[0][0] + sc_start
            q_seq = record.seq[start_idx:end_idx + 1]

            records.append([record, len(q_seq)])

    records = list(map(lambda x: x[0], sorted(records, key=lambda x: x[1])))
    if len(records) < chunk_size * chunk_scale:
        bams = [input_bam]
        lens = [len(records)]
        ds = [1]
    elif do_split:
        n_splits = int(max(6, len(records) // chunk_size))
        new_chunk_size = len(records) // n_splits
        bams = []
        lens = []
        ds = []
        n_split = (len(records) // new_chunk_size) + 1
        if 0 < (len(records) - ((n_split - 1) * new_chunk_size) + new_chunk_size) \
                < new_chunk_size * chunk_scale:
            n_split -= 1
        for i in range(n_split):
            i_start = i * new_chunk_size
            i_end = (i + 1) * \
                new_chunk_size if i < (n_split - 1) else len(records)
            split_input_bam = os.path.join(
                work, region.__str__() + "_split_{}.bam".format(i))
            with pysam.AlignmentFile(input_bam, "rb") as samfile:
                with pysam.AlignmentFile(split_input_bam, "wb",
                                         template=samfile) as out_samfile:
                    for record in records[i_start:i_end]:
                        out_samfile.write(record)
            pysam.sort("-o", "{}.sorted.bam".format(split_input_bam),
                       split_input_bam)
            shutil.move("{}.sorted.bam".format(
                split_input_bam), split_input_bam)
            pysam.index(split_input_bam)

            bams.append(split_input_bam)
            lens.append(i_end - i_start + 1)
            ds.append(1)
    else:
        bams = [input_bam]
        lens = [chunk_size * chunk_scale]
        ds = [chunk_size * chunk_scale / float(len(records))]
    return bams, lens, ds


def read_info(info_file):
    logger = logging.getLogger(read_info.__name__)
    info = {}
    with open(info_file, 'r') as csvfile:
        csvreader = csv.reader(csvfile, delimiter='\t', quotechar='|')
        for row in csvreader:
            info[int(row[0])] = Read_Info(row[1:])
    return info


def find_cigar(alignment):
    logger = logging.getLogger(find_cigar.__name__)
    SOME_BIG_NUMBER = 100
    augmented_alignment = np.append(alignment, [SOME_BIG_NUMBER])
    event_pos = np.append([-1], np.nonzero(np.diff(augmented_alignment)))
    event_len = np.diff(event_pos)
    if sum(event_len) != alignment.shape[0]:
        logger.error("event_len is different from length of alignment: {} vs {}".format(
            sum(event_len), alignment.shape[0]))
        raise Exception

    event_type = augmented_alignment[:-1][event_pos[1:]]
    cigartuple = list(zip(event_type, event_len))
    return cigartuple_to_string(cigartuple)


def extract_new_cigars(region, info_file, out_fasta_file):
    logger = logging.getLogger(extract_new_cigars.__name__)
    info = read_info(info_file)
    aligned_reads = {}
    records = SeqIO.to_dict(SeqIO.parse(out_fasta_file, "fasta"))
    if len(records) <= 1:
        return {}, {}, {}

    if set(map(int, records.keys())) ^ set(range(len(records))):
        logger.error("sequences are missing in the alignment {}".format(
            set(map(int, records.keys())) ^ set(range(len(records)))))
        raise Exception

    alignment = list(map(lambda x: x[1], sorted(map(lambda x: [int(x[0]), list(map(lambda x: 0 if x == "-"
                                                                                   else 1, x[1].seq))],
                                                    records.items()),
                                                key=lambda x: x[0])))
    ref_seq = np.array(alignment[0])
    pos_ref = np.cumsum(alignment[0]) - 1
    alignment = np.array(alignment[1:]) - ref_seq
    alignment = (alignment == 0) * (1 - ref_seq) * (-1) + (alignment == 0) * ref_seq * \
        CIGAR_MATCH + (alignment == 1) * CIGAR_INS + \
        (alignment == -1) * CIGAR_DEL
    N = alignment.shape[0]
    new_cigars = {}
    excess_start = {}
    excess_end = {}
    for i in range(N):
        core_alignment = alignment[i, :][alignment[i, :] >= 0]
        new_cigars[i + 1] = find_cigar(core_alignment)
        if CIGAR_MATCH in alignment[i, :]:
            excess_start[i + 1] = (info[i + 1].pos_start - region.start) - \
                (pos_ref[np.where((alignment[i, :] == CIGAR_MATCH))[0][0]])
            excess_end[i + 1] = (region.end - info[i + 1].pos_end) - (
                max(pos_ref) - (pos_ref[np.where((alignment[i, :] == CIGAR_MATCH))[0][-1]]))
        else:
            excess_start[i + 1] = (info[i + 1].pos_start - region.start)
            excess_end[i + 1] = (region.end - info[i + 1].pos_end)
    return new_cigars, excess_start, excess_end


def extract_consensus(region, out_fasta_file):
    logger = logging.getLogger(extract_consensus.__name__)
    aligned_reads = {}
    records = SeqIO.to_dict(SeqIO.parse(out_fasta_file, "fasta"))
    if len(records) <= 1:
        return {}, {}
    try:
        assert(not set(map(int, records.keys()))
               ^ set(range(1, len(records) + 1)))
    except:
        logger.error("sequences are missing in the alignment {}".format(
            set(map(int, records.keys())) ^ set(range(1, len(records) + 1))))
        raise Exception
    n = len(records)
    align_len = len(records["1"].seq)
    msa = [[] for i in range(n)]
    NUC_to_NUM = {"A": 1, "C": 2, "G": 3, "T": 4, "-": 0}
    NUM_to_NUC = {1: "A", 2: "C", 3: "G", 4: "T", 0: "-"}

    def nuc_to_num_convert(nuc):
        if nuc.upper() not in NUC_to_NUM.keys():
            nuc = "-"
        return NUC_to_NUM[nuc.upper()]
    for i, record in records.items():
        ii = int(i) - 1
        msa[ii] = list(map(lambda x: nuc_to_num_convert(x), record.seq))
    msa = np.array(msa, dtype=int)
    consensus = []
    for i in range(align_len):
        counts = np.histogram(msa[:, i], range(6))
        sorted_chars = np.argsort(counts[0])
        max_char = NUM_to_NUC[sorted_chars[-1]]
        if max_char == "-":
            count_gap = counts[0][sorted_chars[-1]]
            count_max_base = counts[0][sorted_chars[-2]]
            if count_max_base > 0.5 * count_gap:
                max_char = NUM_to_NUC[sorted_chars[-2]]
        consensus.append(max_char)
    consensus = "".join(consensus)
    return msa, consensus


def get_final_msa(region, msa_0, consensus, out_fasta_file_1, out_fasta_file_final):
    logger = logging.getLogger(get_final_msa.__name__)
    aligned_reads = {}
    records = SeqIO.to_dict(SeqIO.parse(out_fasta_file_1, "fasta"))
    if len(records) <= 1:
        return False
    if set(map(int, records.keys())) ^ set(range(2)):
        logger.error("sequences are missing in the alignment {}".format(
            set(map(int, records.keys())) ^ set(range(1, len(records) + 1))))
        raise Exception
    align_len = len(records["0"].seq)
    msa_1 = [[] for i in range(2)]
    for i, record in records.items():
        ii = int(i)
        msa_1[ii] = list(map(lambda x: nuc_to_num_convert(x), record.seq))
    msa_1 = np.array(msa_1, dtype=int)
    consensus_array = np.array(
        list(map(lambda x: nuc_to_num_convert(x), consensus)))
    consensus_cumsum = np.cumsum(consensus_array > 0)
    new_cols = np.where(msa_1[1, :] == 0)[0]
    new_cols -= np.arange(len(new_cols))
    inser_new_cols = []
    for col in new_cols:
        if col > 0:
            inser_new_cols.append(np.where(consensus_cumsum == col)[0][0] + 1)
        else:
            inser_new_cols.append(0)
    msa = np.insert(msa_0, inser_new_cols, 0, axis=1)
    new_consensus_array = np.insert(
        consensus_array, inser_new_cols, 100, axis=0)
    msa = np.insert(msa, 0, 0, axis=0)
    msa[0, np.where(new_consensus_array > 0)[0]] = msa_1[0]
    with open(out_fasta_file_final, "w") as out_fasta:
        for i, seq in enumerate(msa.tolist()):
            out_fasta.write(">%s\n" % i)
            out_fasta.write("%s\n" % "".join(
                map(lambda y: NUM_to_NUC[y], seq)))
    return True


def get_entries(region, info_file, new_cigars, excess_start, excess_end):
    logger = logging.getLogger(get_entries.__name__)
    info = read_info(info_file)
    N = len(new_cigars)

    if len(info) != N:
        logger.error(
            "number of items in info is different from length of new cigars: {} vs {}".format(
                len(info), N))
        raise Exception

    entries = []
    for i in range(1, N + 1):
        entries.append([region.chrom, region.start, region.end, info[i].query_name, info[i].pos,
                        info[i].start_idx, info[i].end_idx,
                        info[i].cigarstring, info[
                            i].del_start, info[i].del_end,
                        info[i].pos_start, info[
                            i].pos_end, new_cigars[i], excess_start[i],
                        excess_end[i]])
    return entries


def merge_cigartuples(tuple1, tuple2):
    logger = logging.getLogger(merge_cigartuples.__name__)
    if not tuple1:
        return tuple2
    if not tuple2:
        return tuple1
    if tuple1[-1][0] == tuple2[0][0]:
        return tuple1[:-1] + [[tuple1[-1][0], tuple1[-1][1] + tuple2[0][1]]] + tuple2[1:]
    return tuple1 + tuple2


def find_realign_dict(realign_bed_file, chrom):
    logger = logging.getLogger(find_realign_dict.__name__)

    realign_bed = get_tmp_file()
    with open(realign_bed_file) as i_f, open(realign_bed, "w") as o_f:
        for line in skip_empty(i_f):
            x = line.strip().split()
            if x[0] == chrom:
                o_f.write(line)

    realign_dict = {}
    chrom_regions = set([])
    with open(realign_bed) as r_f:
        for line in skip_empty(r_f):
            interval = line.strip().split("\t")
            chrom, start, end, query_name = interval[0:4]
            pos, start_idx, end_idx, cigarstring, del_start, del_end, pos_start, pos_end, new_cigar, \
                excess_start, excess_end = interval[6:]
            q_key = "{}_{}_{}".format(query_name, pos, cigarstring)
            if q_key not in realign_dict:
                realign_dict[q_key] = Realign_Read(
                    query_name, chrom, pos, cigarstring)
            realign_dict[q_key].add_realignment(start, end, start_idx, end_idx, del_start,
                                                del_end, pos_start, pos_end, new_cigar, excess_start,
                                                excess_end)
            chrom_regions.add("{}-{}".format(start, end))
    chrom_regions = sorted(
        map(lambda x: list(map(int, x.split("-"))), chrom_regions))
    return realign_dict, chrom_regions


def correct_bam_chrom(input_record):
    work, input_bam, realign_bed_file, ref_fasta_file, chrom = input_record
    thread_logger = logging.getLogger(
        "{} ({})".format(correct_bam_chrom.__name__, multiprocessing.current_process().name))
    try:
        fasta_file = pysam.Fastafile(ref_fasta_file)
        ref_seq = fasta_seq(fasta_file)
        ref_seq.set_chrom(chrom)
        out_sam = os.path.join(work, "output_{}.sam".format(chrom))
        with pysam.AlignmentFile(input_bam, "rb") as samfile:
            with pysam.AlignmentFile(out_sam, "w", template=samfile) as out_samfile:
                realign_dict, chrom_regions = find_realign_dict(
                    realign_bed_file, chrom)
                if chrom_regions:
                    region_cnt = 0
                    next_region = chrom_regions[0]
                    done_regions = False
                else:
                    done_regions = True
                in_active_region = False
                for record in samfile.fetch(chrom):
                    if record.is_unmapped:
                        continue
                    if not done_regions and in_active_region and record.pos > next_region[1]:
                        if region_cnt == len(chrom_regions) - 1:
                            done_regions = True
                        else:
                            region_cnt += 1
                            next_region = chrom_regions[region_cnt]
                    if done_regions or (record.aend < next_region[0]):
                        out_samfile.write(record)
                        continue
                    q_key = "{}_{}_{}".format(
                        record.query_name, record.pos, record.cigarstring)
                    if q_key not in realign_dict:
                        out_samfile.write(record)
                        continue
                    fixed_record = realign_dict[
                        q_key].fix_record(record, ref_seq)
                    out_samfile.write(fixed_record)
                    in_active_region = True
        return out_sam
    except Exception as ex:
        thread_logger.error(traceback.format_exc())
        thread_logger.error(ex)
        return None


def correct_bam_all(work, input_bam, output_bam, ref_fasta_file, realign_bed_file):
    logger = logging.getLogger(correct_bam_all.__name__)
    with pysam.AlignmentFile(input_bam, "rb") as samfile:
        with pysam.AlignmentFile(output_bam, "wb", template=samfile) as out_samfile:
            fasta_file = pysam.Fastafile(ref_fasta_file)
            ref_seq = fasta_seq(fasta_file)
            for chrom in samfile.references:
                ref_seq.set_chrom(chrom)
                realign_dict, chrom_regions = find_realign_dict(
                    realign_bed_file, chrom)
                if chrom_regions:
                    region_cnt = 0
                    next_region = chrom_regions[0]
                    done_regions = False
                else:
                    done_regions = True
                in_active_region = False
                for record in samfile.fetch(chrom):
                    if record.is_unmapped:
                        continue
                    if not done_regions and in_active_region and record.pos > next_region[1]:
                        if region_cnt == len(chrom_regions) - 1:
                            done_regions = True
                        else:
                            region_cnt += 1
                            next_region = chrom_regions[region_cnt]
                    if done_regions or (record.aend < next_region[0]):
                        out_samfile.write(record)
                        continue
                    q_key = "{}_{}_{}".format(
                        record.query_name, record.pos, record.cigarstring)
                    if q_key not in realign_dict:
                        out_samfile.write(record)
                        continue
                    fixed_record = realign_dict[
                        q_key].fix_record(record, ref_seq)
                    out_samfile.write(fixed_record)
                    in_active_region = True
    if os.path.exists(output_bam):
        pysam.index(output_bam)


def concatenate_sam_files(files, output, bam_header):
    logger = logging.getLogger(concatenate_sam_files.__name__)
    good_files = list(filter(lambda x: x and os.path.isfile(x), files))
    fin = fileinput.input(good_files)
    with open(output, "w") as merged_fd:
        with open(bam_header) as bh_fd:
            for line in bh_fd:
                merged_fd.write(line)
        for line in fin:
            if line[0] == "@":
                continue
            merged_fd.write(line)

    fin.close()

    return output


def parallel_correct_bam(work, input_bam, output_bam, ref_fasta_file, realign_bed_file,
                         num_threads):
    logger = logging.getLogger(parallel_correct_bam.__name__)
    if num_threads > 1:
        bam_header = output_bam[:-4] + ".header"
        with open(bam_header, "w") as h_f:
            h_f.write(pysam.view("-H", input_bam,))

        map_args = []
        with pysam.AlignmentFile(input_bam, "rb") as samfile:
            for chrom in samfile.references:
                map_args.append(
                    (work, input_bam, realign_bed_file, ref_fasta_file, chrom))

        try:
            with multiprocessing.Pool(num_threads) as pool:
                sams = pool.map_async(correct_bam_chrom, map_args).get()
        except Exception as inst:
            logger.error(inst)
            traceback.print_exc()
            raise Exception

        for o in sams:
            if o is None:
                raise Exception("correct_bam_chrom failed!")

        output_sam = output_bam[:-4] + ".sam"
        concatenate_sam_files(sams, output_sam, bam_header)
        if os.path.exists(output_sam):
            with pysam.AlignmentFile(output_sam, "r") as samfile:
                with pysam.AlignmentFile(output_bam, "wb",
                                         template=samfile) as out_samfile:
                    for record in samfile.fetch():
                        out_samfile.write(record)
            pysam.index(output_bam)

            for sam in [bam_header] + sams:
                os.remove(sam)
    else:
        correct_bam_all(work, input_bam, output_bam,
                        ref_fasta_file, realign_bed_file)


def run_msa(in_fasta_file, match_score, mismatch_penalty, gap_open_penalty, gap_ext_penalty,
            msa_binary):
    logger = logging.getLogger(run_msa.__name__)
    if not os.path.exists(msa_binary):
        raise IOError("File not found: {}".format(msa_binary))
    out_fasta_file = ".".join(in_fasta_file.split(".")[:-1]) + "_aligned.fasta"
    cmd = "{} -A {} -B {} -O {} -E {} -i {} -o {}".format(msa_binary, match_score, mismatch_penalty,
                                                          gap_open_penalty, gap_ext_penalty,
                                                          in_fasta_file, out_fasta_file)
    if not os.path.exists(out_fasta_file):
        run_shell_command(cmd, run_logger=logger)
    return out_fasta_file


def do_realign(region, info_file, max_realign_dp, thr_realign=0.0135):
    logger = logging.getLogger(do_realign.__name__)
    sum_nm_snp = 0
    sum_nm_indel = 0
    c = 0
    with open(info_file) as i_f:
        for line in skip_empty(i_f):
            x = line.strip().split()
            sum_nm_snp += int(x[-2])
            sum_nm_indel += int(x[-1])
            c += 1
    eps = 0.0001
    if (c < max_realign_dp) and (
            (sum_nm_snp + sum_nm_indel
             ) / float(c + eps) / float(region.span() + eps)
            > thr_realign):
        return True
    return False


def find_var(out_fasta_file, snp_min_af, del_min_af, ins_min_af, scale_maf, simplify):
    # Find variants from MSA:
    # In each column the AF is calculated
    # The low AF vars in each column are discarded and the variant is extracted
    logger = logging.getLogger(find_var.__name__)
    records = SeqIO.to_dict(SeqIO.parse(out_fasta_file, "fasta"))
    if set(map(int, records.keys())) ^ set(range(len(records))):
        logger.error("sequences are missing in the alignment {}".format(
            set(map(int, records.keys())) ^ set(range(len(records)))))
        raise Exception
    alignment = np.array(list(map(lambda x: x[1], sorted(map(lambda x: [int(x[0]), list(map(
        lambda x: NUC_to_NUM[x.upper()], x[1].seq))], records.items()),
        key=lambda x: x[0]))))
    ref_seq = alignment[0, :]
    counts = np.zeros((5, alignment.shape[1]))
    for i in range(5):
        counts[i, :] = np.sum(alignment[1:, :] == i, 0)
    counts = counts.astype(int)
    alt_seq = []
    afs = []
    i_afs = []
    for i in range(alignment.shape[1]):
        ref_base = ref_seq[i]
        if ref_base == 5:
            ref_count = counts[0, i]
        else:
            ref_count = counts[ref_base, i]
        sorted_idx = np.argsort(counts[:, i])
        alt_base = sorted_idx[-1]
        if alt_base == ref_base:
            alt_base = sorted_idx[-2]
        alt_count = counts[alt_base, i]
        af = alt_count / (alt_count + ref_count + 0.0001)
        if ((alt_base != '-' and ref_base == "-" and af > ins_min_af) or
                (alt_base != '-' and ref_base != "-" and af > snp_min_af) or
                (alt_base == '-' and ref_base != "-" and af > del_min_af)):
            alt_seq.append(alt_base)
            afs.append(af)
            i_afs.append(i)
        else:
            alt_seq.append(ref_base)
    if afs:
        afs = np.array(afs)
        # thr=np.percentile(afs,0.75)
        # if np.std(afs)>0.02:
        for ii in np.where(afs <= (np.max(afs) * 0.6))[0]:
            # for ii in np.where(afs<=thr)[0]:
            alt_seq[i_afs[ii]] = ref_seq[i_afs[ii]]
    afs = np.array(afs)
    ref_seq_ = "".join(
        map(lambda x: NUM_to_NUC[x], filter(lambda x: x > 0, ref_seq)))
    alt_seq_ = "".join(
        map(lambda x: NUM_to_NUC[x], filter(lambda x: x > 0, alt_seq)))
    if not simplify:
        variants = [[0, ref_seq_, alt_seq_, afs]]
    else:
        variants = []
        bias = 0
        current_ref = []
        current_alt = []
        current_af = []
        current_bias = 0
        is_ins = False
        is_del = False
        done = False
        for i, (r, a) in enumerate(zip(list(ref_seq) + [0], list(alt_seq) + [0])):
            if i in i_afs:
                af = afs[i_afs.index(i)]
            else:
                af = 0
            if r == a:
                done = True
            else:
                if r == 0 and a != 0:
                    if not is_ins:
                        done = True
                elif r != 0 and a == 0:
                    if not is_del:
                        done = True
                else:
                    done = True
            if done:
                if current_alt:
                    rr = "".join(map(lambda x: NUM_to_NUC[
                        x], filter(lambda x: x > 0, current_ref)))
                    aa = "".join(map(lambda x: NUM_to_NUC[
                        x], filter(lambda x: x > 0, current_alt)))
                    variants.append(
                        [current_bias, rr, aa, np.array(current_af)])
                    done = False
                    current_ref = []
                    current_alt = []
                    current_af = []
                    current_bias = bias
                    is_ins = False
                    is_del = False
                done = False
            if not done:
                current_ref.append(r)
                current_alt.append(a)
                current_af.append(af)
                is_ins = r == 0 and a != 0
                is_del = r != 0 and a == 0
            if r != 0:
                bias += 1

    return variants


def TrimREFALT(ref, alt, pos):
    logger = logging.getLogger(TrimREFALT.__name__)
    alte = len(alt)
    refe = len(ref)
    while (alte > 1 and refe > 1 and alt[alte - 1] == ref[refe - 1]):
        alte -= 1
        refe -= 1

    alt = alt[0:alte]
    ref = ref[0:refe]
    s = 0
    while (s < (len(alt) - 1) and s < (len(ref) - 1) and alt[s] == ref[s]):
        s += 1

    alt = alt[s:]
    ref = ref[s:]
    pos += s
    return ref, alt, pos


def run_realignment(input_record):
    work, ref_fasta_file, target_region, pad, chunk_size, chunk_scale, \
        snp_min_af, del_min_af, ins_min_af, len_chr, input_bam, \
        match_score, mismatch_penalty, gap_open_penalty, gap_ext_penalty, \
        max_realign_dp, \
        filter_duplicate, \
        msa_binary, get_var, do_split = input_record

    ref_fasta = pysam.Fastafile(ref_fasta_file)
    thread_logger = logging.getLogger(
        "{} ({})".format(run_realignment.__name__, multiprocessing.current_process().name))

    try:
        region = Region(target_region, pad, len_chr)
        not_realigned_region = None
        original_tempdir = tempfile.tempdir
        bed_tempdir = os.path.join(
            work, "bed_tmpdir_{}".format(region.__str__()))
        if not os.path.exists(bed_tempdir):
            os.mkdir(bed_tempdir)
        tempfile.tempdir = bed_tempdir
        variants = []
        all_entries = []
        input_bam_splits, lens_splits, ds_splits = split_bam_to_chunks(
            work, region, input_bam, chunk_size, chunk_scale, do_split or not get_var, filter_duplicate)
        new_seqs = []
        new_ref_seq = ""
        skipped = 0
        if len(input_bam_splits) <= 2:
            scale_maf = 1
        else:
            scale_maf = 2
        afss = []
        for i, i_bam in enumerate(input_bam_splits):
            in_fasta_file, info_file = prepare_fasta(
                work, region, i_bam, ref_fasta_file, True, i, ds_splits[i], filter_duplicate)
            if do_realign(region, info_file, max_realign_dp):
                out_fasta_file_0 = run_msa(
                    in_fasta_file, match_score, mismatch_penalty, gap_open_penalty,
                    gap_ext_penalty, msa_binary)
                if get_var:
                    var = find_var(
                        out_fasta_file_0, snp_min_af, del_min_af, ins_min_af, scale_maf, False)
                    assert(len(var) == 1)
                    _, ref_seq_, alt_seq_, afs = var[0]
                    afss.append(afs)
                    new_ref_seq = ref_seq_
                    new_seqs.append(alt_seq_)
                new_cigars, excess_start, excess_end = extract_new_cigars(
                    region, info_file, out_fasta_file_0)
                if new_cigars:
                    entries = get_entries(
                        region, info_file, new_cigars, excess_start, excess_end)
                    all_entries.extend(entries)
            else:
                skipped += 1
        if get_var:
            if new_seqs:
                for i in range(skipped):
                    new_seqs = [new_ref_seq] + new_seqs
                new_seqs = [new_ref_seq] + new_seqs
                consensus_fasta = os.path.join(
                    work, region.__str__() + "_consensus.fasta")
                with open(consensus_fasta, "w") as output_handle:
                    for i, seq in enumerate(new_seqs):
                        record = SeqRecord(
                            Seq(seq, DNAAlphabet.letters), id=str(i), description="")
                        SeqIO.write(record, output_handle, "fasta")
                consensus_fasta_aligned = run_msa(
                    consensus_fasta, match_score, mismatch_penalty, gap_open_penalty,
                    gap_ext_penalty, msa_binary)
                vars_ = find_var(
                    consensus_fasta_aligned, snp_min_af, del_min_af, ins_min_af, 1, True)
                for var in vars_:
                    pos_, ref_seq, alt_seq, afs = var
                    if ref_seq != alt_seq:
                        ref, alt, pos = ref_seq, alt_seq, int(
                            region.start) + 1 + pos_
                        if pos > 1:
                            num_add_before = min(40, pos - 1)
                            before = ref_fasta.fetch(
                                region.chrom, pos - num_add_before, pos - 1).upper()
                            pos -= num_add_before - 1
                            ref = before + ref
                            alt = before + alt
                        ref, alt, pos = TrimREFALT(
                            ref, alt, pos)
                        a = int(np.ceil(np.max(afs) * len(afss)))
                        af = sum(sorted(map(lambda x:
                                            np.max(x) if x.shape[0] > 0 else 0,
                                            afss))[-a:]) / float(len(afss))
                        dp = int(sum(lens_splits))
                        ao = int(af * dp)
                        ro = dp - ao
                        if ref == "" and pos > 1:
                            pos -= 1
                            r_ = ref_fasta.fetch(
                                region.chrom, pos - 1, pos).upper()
                            ref = r_ + ref
                            alt = r_ + alt
                        if alt == "" and pos > 1:
                            pos -= 1
                            r_ = ref_fasta.fetch(
                                region.chrom, pos - 1, pos).upper()
                            ref = r_ + ref
                            alt = r_ + alt
                        variants.append(
                            [region.chrom, pos, ref, alt, dp, ro, ao])
            else:
                if skipped > 0:
                    not_realigned_region = target_region

        shutil.rmtree(bed_tempdir)
        tempfile.tempdir = original_tempdir
        return all_entries, variants, not_realigned_region
    except Exception as ex:
        thread_logger.error(traceback.format_exc())
        thread_logger.error(ex)
        return None


class fasta_seq:

    def __init__(self, fasta_pysam):
        self.fasta_pysam = fasta_pysam
        self.chrom = ""

    def set_chrom(self, chrom):
        self.chrom = chrom

    def get_seq(self, start, end=[]):
        if not end:
            end = start + 1
        return self.fasta_pysam.fetch(self.chrom, start, end).upper()


def extend_regions_hp(region_bed_file, extended_region_bed_file, ref_fasta_file,
                      chrom_lengths, pad):
    # If boundaries of regions are in the middle of a homopolymer, this function extends the region
    # to fully include the homopolymer
    logger = logging.getLogger(extend_regions_hp.__name__)
    with pysam.Fastafile(ref_fasta_file) as ref_fasta:
        intervals = []
        with open(region_bed_file) as r_f:
            for line in skip_empty(r_f):
                interval = line.strip().split("\t")
                chrom, start, end = interval[0:3]
                start, end = int(start), int(end)
                s_base = ref_fasta.fetch(
                    chrom, start - pad, start - pad + 1).upper()
                e_base = ref_fasta.fetch(
                    chrom, end + pad, end + pad + 1).upper()
                new_start = start
                i = start - pad - 1
                while True:
                    base = ref_fasta.fetch(chrom, i, i + 1).upper()
                    if base == s_base:
                        new_start -= 1
                    else:
                        break
                    i -= 1
                    if i <= 3:
                        break
                new_end = end
                i = end + pad + 1
                while True:
                    base = ref_fasta.fetch(chrom, i, i + 1).upper()
                    if base == e_base:
                        new_end += 1
                    else:
                        break
                    i += 1
                    if i >= chrom_lengths[chrom] - 3:
                        break
                if ref_fasta.fetch(chrom, new_end + pad, new_end + pad + 1
                                   ).upper() == ref_fasta.fetch(chrom, new_end - 1 + pad,
                                                                new_end + pad).upper():
                    new_end += 1
                if ref_fasta.fetch(chrom, new_start - pad, new_start - pad + 1
                                   ).upper() == ref_fasta.fetch(chrom, new_start - pad + 1,
                                                                new_start - pad + 2).upper():
                    new_start -= 1
                seq, new_seq = ref_fasta.fetch(chrom, start - pad, end + pad + 1).upper(
                ), ref_fasta.fetch(chrom, new_start - pad, new_end + pad + 1).upper()
                intervals.append([chrom, new_start, new_end])

        tmp_ = get_tmp_file()
        write_tsv_file(tmp_, intervals)
        bedtools_sort(tmp_, output_fn=extended_region_bed_file,
                      run_logger=logger)


def check_rep(ref_seq, left_right, w):
    logger = logging.getLogger(check_rep.__name__)
    if len(ref_seq) < 2 * w:
        return False
    if left_right == "left":
        return ref_seq[0:w] == ref_seq[w:2 * w] and len(set(ref_seq[0:2 * w])) > 1
    elif left_right == "right":
        return ref_seq[-w:] == ref_seq[-2 * w:-w] and len(set(ref_seq[-2 * w:])) > 1
    else:
        logger.error("Wrong left/right value: {}".format(left_right))
        raise Exception


def extend_regions_repeat(region_bed_file, extended_region_bed_file, ref_fasta_file,
                          chrom_lengths, pad):
    logger = logging.getLogger(extend_regions_repeat.__name__)
    with pysam.Fastafile(ref_fasta_file) as ref_fasta:
        intervals = []
        with open(region_bed_file) as r_f:
            for line in skip_empty(r_f):
                interval = line.strip().split("\t")
                chrom, start, end = interval[0:3]
                start, end = int(start), int(end)
                w = 3
                new_start = max(start - pad - w, 1)
                new_end = min(end + pad + w, chrom_lengths[chrom] - 2)
                ref_seq = ref_fasta.fetch(
                    chrom, new_start, new_end + 1).upper()
                cnt_s = 0
                while check_rep(ref_seq, "left", 2):
                    new_start -= 2
                    ref_seq = ref_fasta.fetch(
                        chrom, new_start, new_end + 1).upper()
                    cnt_s += 2
                if cnt_s == 0:
                    while check_rep(ref_seq, "left", 3):
                        new_start -= 3
                        ref_seq = ref_fasta.fetch(
                            chrom, new_start, new_end + 1).upper()
                        cnt_s += 3
                if cnt_s == 0:
                    while check_rep(ref_seq, "left", 4):
                        new_start -= 4
                        ref_seq = ref_fasta.fetch(
                            chrom, new_start, new_end + 1).upper()
                        cnt_s += 4
                if cnt_s == 0:
                    new_start += w
                    ref_seq = ref_fasta.fetch(
                        chrom, new_start, new_end + 1).upper()
                if cnt_s == 0:
                    while check_rep(ref_seq, "left", 2):
                        new_start -= 2
                        ref_seq = ref_fasta.fetch(
                            chrom, new_start, new_end + 1).upper()
                        cnt_s += 2
                if cnt_s == 0:
                    while check_rep(ref_seq, "left", 3):
                        new_start -= 3
                        ref_seq = ref_fasta.fetch(
                            chrom, new_start, new_end + 1).upper()
                        cnt_s += 3
                if cnt_s == 0:
                    while check_rep(ref_seq, "left", 4):
                        new_start -= 4
                        ref_seq = ref_fasta.fetch(
                            chrom, new_start, new_end + 1).upper()
                        cnt_s += 4
                cnt_e = 0
                while check_rep(ref_seq, "right", 2):
                    new_end += 2
                    ref_seq = ref_fasta.fetch(
                        chrom, new_start, new_end + 1).upper()
                    cnt_e += 2
                if cnt_e == 0:
                    while check_rep(ref_seq, "right", 3):
                        new_end += 3
                        ref_seq = ref_fasta.fetch(
                            chrom, new_start, new_end + 1).upper()
                        cnt_e += 3
                if cnt_e == 0:
                    while check_rep(ref_seq, "right", 4):
                        new_end += 4
                        ref_seq = ref_fasta.fetch(
                            chrom, new_start, new_end + 1).upper()
                        cnt_e += 4

                if cnt_e == 0:
                    new_end -= w
                    ref_seq = ref_fasta.fetch(
                        chrom, new_start, new_end + 1).upper()
                if cnt_e == 0:
                    while check_rep(ref_seq, "right", 2):
                        new_end += 2
                        ref_seq = ref_fasta.fetch(
                            chrom, new_start, new_end + 1).upper()
                        cnt_e += 2
                if cnt_e == 0:
                    while check_rep(ref_seq, "right", 3):
                        new_end += 3
                        ref_seq = ref_fasta.fetch(
                            chrom, new_start, new_end + 1).upper()
                        cnt_e += 3
                if cnt_e == 0:
                    while check_rep(ref_seq, "right", 4):
                        new_end += 4
                        ref_seq = ref_fasta.fetch(
                            chrom, new_start, new_end + 1).upper()
                        cnt_e += 4
                intervals.append([chrom, new_start + pad, new_end - pad])

        tmp_ = get_tmp_file()
        write_tsv_file(tmp_, intervals, add_fields=[".", ".", "."])
        bedtools_sort(tmp_, output_fn=extended_region_bed_file,
                      run_logger=logger)


def long_read_indelrealign(work, input_bam, output_bam, output_vcf, output_not_realigned_bed,
                           region_bed_file,
                           ref_fasta_file, num_threads, pad,
                           chunk_size, chunk_scale, snp_min_af, del_min_af, ins_min_af,
                           match_score, mismatch_penalty, gap_open_penalty, gap_ext_penalty,
                           max_realign_dp,
                           do_split,
                           filter_duplicate,
                           msa_binary):
    logger = logging.getLogger(long_read_indelrealign.__name__)

    logger.info("-----------Resolve variants for INDELS (long-read)---------")
    if not os.path.exists(work):
        os.mkdir(work)

    if not output_bam and not output_vcf:
        logger.error(
            "At least one of --output_bam or --output_vcf should be provided.")
        raise Exception

    chrom_lengths = {}
    chroms_order = get_chromosomes_order(bam=input_bam)
    with pysam.AlignmentFile(input_bam, "rb") as samfile:
        lens = []
        for length in samfile.lengths:
            lens.append(length)
        chroms = []
        for chrom in samfile.references:
            chroms.append(chrom)
        chrom_lengths = dict(zip(chroms, lens))

    extended_region_bed_file = os.path.join(work, "regions_extended.bed")
    extend_regions_repeat(
        region_bed_file, extended_region_bed_file, ref_fasta_file, chrom_lengths, pad)
    region_bed_file = extended_region_bed_file

    region_bed_merged = region_bed_file
    len_merged = 0
    with open(region_bed_merged) as r_b:
        for line in skip_empty(r_b):
            len_merged += 1
    while True:
        region_bed_merged_tmp = bedtools_merge(
            region_bed_merged, args=" -d {}".format(pad * 2), run_logger=logger)
        len_tmp = 0
        with open(region_bed_merged_tmp) as r_b:
            for line in skip_empty(r_b):
                len_tmp += 1
        if len_tmp == len_merged:
            break
        region_bed_merged = region_bed_merged_tmp
        len_merged = len_tmp
    shutil.copyfile(region_bed_merged, os.path.join(
        work, "regions_merged.bed"))

    target_regions = read_tsv_file(region_bed_merged, fields=range(3))
    target_regions = list(
        map(lambda x: [x[0], int(x[1]), int(x[2])], target_regions))

    get_var = True if output_vcf else False
    map_args = []
    for target_region in target_regions:
        map_args.append((work, ref_fasta_file, target_region, pad, chunk_size,
                         chunk_scale, snp_min_af, del_min_af, ins_min_af,
                         chrom_lengths[target_region[0]], input_bam,
                         match_score, mismatch_penalty, gap_open_penalty, gap_ext_penalty,
                         max_realign_dp, filter_duplicate,
                         msa_binary, get_var, do_split))

    shuffle(map_args)
    try:
        if num_threads == 1:
            realign_output = [run_realignment(w) for w in map_args]
        else:
            with multiprocessing.Pool(num_threads) as pool:
                realign_output = pool.map_async(run_realignment, map_args).get()
    except Exception as inst:
        logger.error(inst)
        traceback.print_exc()
        raise Exception

    for o in realign_output:
        if o is None:
            raise Exception("run_realignment failed!")

    realign_entries = list(map(lambda x: x[0], realign_output))

    realign_variants = list(map(lambda x: x[1], realign_output))
    realign_variants = [v for var in realign_variants for v in var]
    realign_variants = list(filter(None, realign_variants))
    not_realigned_regions = list(map(lambda x: x[2], realign_output))
    not_realigned_regions = list(filter(None, not_realigned_regions))

    if get_var:
        with open(output_vcf, "w") as o_f:
            o_f.write("#" + "\t".join(["CHROM", "POS", "ID", "REF",
                                       "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]) + "\n")
            realign_variants = sorted(realign_variants, key=lambda x:
                                      [chroms_order[x[0]], x[1]])
            for variant in realign_variants:
                if variant:
                    chrom, pos, ref, alt, dp, ro, ao = variant
                    line = "\t".join([chrom, str(pos), ".", ref, alt, "100", ".",
                                      "DP={};RO={};AO={}".format(dp, ro, ao),
                                      "GT:DP:RO:AO", "0/1:{}:{}:{}".format(
                                          dp, ro, ao), ])
                    o_f.write(line + "\n")
    with open(output_not_realigned_bed, "w") as o_f:
        for x in not_realigned_regions:
            o_f.write("\t".join(map(str, x)) + "\n")

    original_tempdir = tempfile.tempdir
    bed_tempdir = os.path.join(work, "bed_tmpdir")
    if not os.path.exists(bed_tempdir):
        os.mkdir(bed_tempdir)
    tempfile.tempdir = bed_tempdir

    if realign_entries:
        realign_entries = functools.reduce(lambda x, y: x + y, realign_entries)
    realign_bed_file = os.path.join(work, "realign.bed")
    realign_entries.sort()
    with open(realign_bed_file, "w") as o_f:
        for x in realign_entries:
            o_f.write("\t".join(map(str, x[0:4] + [".", "."] + x[4:])) + "\n")

    if output_bam:
        parallel_correct_bam(work, input_bam, output_bam, ref_fasta_file,
                             realign_bed_file, num_threads)

    shutil.rmtree(bed_tempdir)
    tempfile.tempdir = original_tempdir

    logger.info("Done")


if __name__ == '__main__':

    FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    logging.basicConfig(level=logging.INFO, format=FORMAT)
    logger = logging.getLogger(__name__)

    parser = argparse.ArgumentParser(description='realign indels using MSA')
    parser.add_argument('--input_bam', type=str, help='input bam')
    parser.add_argument('--output_vcf', type=str,
                        help='output_vcf (needed for variant prediction)', default=None)
    parser.add_argument('--output_not_realigned_bed', type=str,
                        help='output_not_realigned_bed', required=True)
    parser.add_argument('--output_bam', type=str,
                        help='output_bam (needed for getting the realigned bam)', default=None)
    parser.add_argument('--region_bed', type=str,
                        help='region_bed', required=True)
    parser.add_argument('--reference', type=str,
                        help='reference fasta filename', required=True)
    parser.add_argument('--work', type=str,
                        help='work directory', required=True)
    parser.add_argument('--num_threads', type=int,
                        help='number of threads', default=1)
    parser.add_argument('--pad', type=int,
                        help='#base padding to the regions', default=1)
    parser.add_argument('--chunk_size', type=int,
                        help='chuck split size for high depth', default=600)
    parser.add_argument('--chunk_scale', type=float,
                        help='chuck scale size for high depth', default=1.5)
    parser.add_argument('--snp_min_af', type=float,
                        help='SNP min allele freq', default=0.05)
    parser.add_argument('--ins_min_af', type=float,
                        help='INS min allele freq', default=0.05)
    parser.add_argument('--del_min_af', type=float,
                        help='DEL min allele freq', default=0.05)
    parser.add_argument('--match_score', type=int,
                        help='match score', default=10)
    parser.add_argument('--mismatch_penalty', type=int,
                        help='penalty for having a mismatch', default=8)
    parser.add_argument('--gap_open_penalty', type=int,
                        help='penalty for opening a gap', default=8)
    parser.add_argument('--gap_ext_penalty', type=int,
                        help='penalty for extending a gap', default=6)
    parser.add_argument('--max_realign_dp', type=int,
                        help='max coverage for realign region', default=1000)
    parser.add_argument('--do_split',
                        help='Split bam for high coverage regions (in variant-calling mode).',
                        action="store_true")
    parser.add_argument('--filter_duplicate',
                        help='filter duplicate reads in analysis',
                        action="store_true")
    parser.add_argument('--msa_binary', type=str,
                        help='MSA binary', default="../bin/msa")
    args = parser.parse_args()
    logger.info(args)

    try:
        processor = long_read_indelrealign(args.work, args.input_bam, args.output_bam,
                                           args.output_vcf, args.output_not_realigned_bed,
                                           args.region_bed, args.reference,
                                           args.num_threads, args.pad, args.chunk_size,
                                           args.chunk_scale, args.snp_min_af, args.del_min_af,
                                           args.ins_min_af, args.match_score,
                                           args.mismatch_penalty, args.gap_open_penalty,
                                           args.gap_ext_penalty,
                                           args.gap_ext_penalty,
                                           args.max_realign_dp,
                                           args.do_split,
                                           args.filter_duplicate,
                                           args.msa_binary)
    except Exception as e:
        logger.error(traceback.format_exc())
        logger.error("Aborting!")
        logger.error(
            "long_read_indelrealign.py failure on arguments: {}".format(args))
        raise e
