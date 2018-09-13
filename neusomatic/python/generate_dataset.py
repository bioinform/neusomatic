#-------------------------------------------------------------------------
# generate_dataset.py
# Use the input filtered candidates to prepare and extracted features to generate datasets to
# be used by the NeuSomatic network.
#-------------------------------------------------------------------------
import argparse
import base64
import multiprocessing
import os
import traceback
import pickle
import zlib
import logging

import numpy as np
import pybedtools
import pysam
from scipy.misc import imresize

from split_bed import split_region
from utils import concatenate_vcfs, get_chromosomes_order


NUC_to_NUM_hp = {"A": 1, "C": 2, "G": 3, "T": 4, "N": 5}
NUC_to_NUM = {"a": 1, "c": 2, "g": 3, "t": 4, "-": 5, " ": 0, "A": 1, "C": 2, "G": 3, "T": 4, "N": 5,
              "b": 1, "d": 2, "h": 3, "u": 4}
snp = {"a": "b", "c": "d", "g": "h", "t": "u"}
NUM_to_NUC = {v: k for k, v in NUC_to_NUM.iteritems()}
NUM_to_NUC_hp = {v: k for k, v in NUC_to_NUM_hp.iteritems()}


def get_type(ref, alt):
    logger = logging.getLogger(get_type.__name__)
    len_diff = len(ref) - len(alt.split(",")[0])
    if len_diff > 0:
        return "DEL"
    elif len_diff < 0:
        return "INS"
    else:
        return "SNP"


def get_variant_matrix_tabix(ref_file, count_bed, record, matrix_base_pad):
    logger = logging.getLogger(get_variant_matrix_tabix.__name__)
    chrom, pos, ref, alt = record[0:4]
    fasta_file = pysam.Fastafile(ref_file)
    tb = pysam.TabixFile(count_bed, parser=pysam.asTuple())
    tabix_records = tb.fetch(
        chrom, pos - matrix_base_pad, pos + matrix_base_pad)

    NUC_to_NUM_tabix = {"A": 1, "C": 2, "G": 3, "T": 4, "-": 0}
    matrix_ = []
    bq_matrix_ = []
    mq_matrix_ = []
    st_matrix_ = []
    lsc_matrix_ = []
    rsc_matrix_ = []
    tag_matrices_ = [[], [], [], [], []]
    ref_array = []
    col_pos_map = {}
    cnt = 0
    curr_pos = max(1, pos - matrix_base_pad)
    for rec in tabix_records:
        pos_ = int(rec[1])
        if pos_ > pos + matrix_base_pad:
            continue
        ref_base = rec[3]
        if ref_base.upper() not in "ACGT-":
            ref_base = "-"
        if pos_ in col_pos_map and ref_base != "-":
            continue
        if pos_ > (curr_pos):
            refs = fasta_file.fetch(
                chrom, curr_pos - 1, pos_ - 1).upper().replace("N", "-")
            for i in range(curr_pos, pos_):
                ref_base_ = refs[i - curr_pos]
                if ref_base_.upper() not in "ACGT-":
                    ref_base_ = "-"
                matrix_.append([0, 0, 0, 0, 0])
                bq_matrix_.append([0, 0, 0, 0, 0])
                mq_matrix_.append([0, 0, 0, 0, 0])
                st_matrix_.append([0, 0, 0, 0, 0])
                lsc_matrix_.append([0, 0, 0, 0, 0])
                rsc_matrix_.append([0, 0, 0, 0, 0])
                for iii in range(len(tag_matrices_)):
                    tag_matrices_[iii].append([0, 0, 0, 0, 0])
                ref_array.append(NUC_to_NUM_tabix[ref_base_])
                # if ref_base_ != "-" and i not in col_pos_map:
                if i not in col_pos_map:
                    col_pos_map[i] = cnt
                cnt += 1
            curr_pos = pos_
        if pos_ == (curr_pos) and ref_base == "-" and pos_ not in col_pos_map:
            ref_base_ = fasta_file.fetch(
                chrom, pos_ - 1, pos_).upper().replace("N", "-")
            if ref_base_.upper() not in "ACGT-":
                ref_base_ = "-"
            matrix_.append([0, 0, 0, 0, 0])
            bq_matrix_.append([0, 0, 0, 0, 0])
            mq_matrix_.append([0, 0, 0, 0, 0])
            st_matrix_.append([0, 0, 0, 0, 0])
            lsc_matrix_.append([0, 0, 0, 0, 0])
            rsc_matrix_.append([0, 0, 0, 0, 0])
            for iii in range(len(tag_matrices_)):
                tag_matrices_[iii].append([0, 0, 0, 0, 0])
            ref_array.append(NUC_to_NUM_tabix[ref_base_])
            # if ref_base_ != "-" and pos_ not in col_pos_map:
            if pos_ not in col_pos_map:
                col_pos_map[pos_] = cnt
                cnt += 1
            curr_pos = pos_ + 1
        matrix_.append(map(int, rec[4].split(":")))
        bq_matrix_.append(map(int, rec[5].split(":")))
        mq_matrix_.append(map(int, rec[6].split(":")))
        st_matrix_.append(map(int, rec[7].split(":")))
        lsc_matrix_.append(map(int, rec[8].split(":")))
        rsc_matrix_.append(map(int, rec[9].split(":")))
        for iii in range(len(tag_matrices_)):
            tag_matrices_[iii].append(map(int, rec[10 + iii].split(":")))
        ref_array.append(NUC_to_NUM_tabix[ref_base])
        if ref_base != "-" and pos_ not in col_pos_map:
            col_pos_map[pos_] = cnt
        cnt += 1
        curr_pos = pos_ + 1

    if curr_pos < pos + matrix_base_pad + 1:
        refs = fasta_file.fetch(
            chrom, curr_pos - 1, pos + matrix_base_pad).upper().replace("N", "-")
        for i in range(curr_pos, pos + matrix_base_pad + 1):
            ref_base_ = refs[i - curr_pos]
            if ref_base_.upper() not in "ACGT-":
                ref_base_ = "-"
            matrix_.append([0, 0, 0, 0, 0])
            bq_matrix_.append([0, 0, 0, 0, 0])
            mq_matrix_.append([0, 0, 0, 0, 0])
            st_matrix_.append([0, 0, 0, 0, 0])
            lsc_matrix_.append([0, 0, 0, 0, 0])
            rsc_matrix_.append([0, 0, 0, 0, 0])
            for iii in range(len(tag_matrices_)):
                tag_matrices_[iii].append([0, 0, 0, 0, 0])
            ref_array.append(NUC_to_NUM_tabix[ref_base_])
            # if ref_base_ != "-" and i not in col_pos_map:
            if i not in col_pos_map:
                col_pos_map[i] = cnt
            cnt += 1
        curr_pos = pos + matrix_base_pad + 1

    matrix_ = np.array(matrix_).transpose()
    bq_matrix_ = np.array(bq_matrix_).transpose()
    mq_matrix_ = np.array(mq_matrix_).transpose()
    st_matrix_ = ((np.array(st_matrix_).transpose() / 100.0)) * matrix_
    lsc_matrix_ = np.array(lsc_matrix_).transpose()
    rsc_matrix_ = np.array(rsc_matrix_).transpose()
    for iii in range(len(tag_matrices_)):
        tag_matrices_[iii] = np.array(tag_matrices_[iii]).transpose()

    ref_array = np.array(ref_array)
    return matrix_, bq_matrix_, mq_matrix_, st_matrix_, lsc_matrix_, rsc_matrix_, tag_matrices_, ref_array, col_pos_map


def align_tumor_normal_matrices(record, tumor_matrix_, bq_tumor_matrix_, mq_tumor_matrix_, st_tumor_matrix_,
                                lsc_tumor_matrix_, rsc_tumor_matrix_,
                                tag_tumor_matrices_, tumor_ref_array, tumor_col_pos_map, normal_matrix_,
                                bq_normal_matrix_, mq_normal_matrix_, st_normal_matrix_,
                                lsc_normal_matrix_, rsc_normal_matrix_, tag_normal_matrices_,
                                normal_ref_array, normal_col_pos_map):
    logger = logging.getLogger(align_tumor_normal_matrices.__name__)
    if not tumor_col_pos_map:
        logger.error("record: {}".format(record))
        raise(RuntimeError("tumor_col_pos_map is empty."))

    tumor_col_pos_map[max(tumor_col_pos_map.keys()) +
                      1] = tumor_matrix_.shape[1]
    normal_col_pos_map[max(normal_col_pos_map.keys()) +
                       1] = normal_matrix_.shape[1]

    if set(tumor_col_pos_map.keys()) ^ set(normal_col_pos_map.keys()):
        logger.error("record: {}".format(record))
        logger.error("normal_col_pos_map: {}".format(normal_col_pos_map))
        logger.error("tumor_col_pos_map: {}".format(tumor_col_pos_map))
        raise(RuntimeError(
            "tumor_col_pos_map  and normal_col_pos_map have different keys."))

    pT = map(lambda x: tumor_col_pos_map[x], sorted(tumor_col_pos_map.keys()))
    pN = map(lambda x: normal_col_pos_map[
             x], sorted(normal_col_pos_map.keys()))

    if pT[0] != pN[0]:
        logger.error("record: {}".format(record))
        logger.error("pT, pN: {}, {}".format(pT, pN))
        raise(RuntimeError(
            "pT[0] != pN[0]"))

    min_i = pT[0]
    cols_T = np.ones(tumor_matrix_.shape[1] + 1, int)
    cols_N = np.ones(normal_matrix_.shape[1] + 1, int)
    current_col_T = min_i + 1
    current_col_N = min_i + 1
    for i in range(1, len(pT)):
        current_col_T += pT[i] - pT[i - 1]
        current_col_N += pN[i] - pN[i - 1]
        cols_T[pT[i]] += max(0, current_col_N - current_col_T)
        cols_N[pN[i]] += max(0, current_col_T - current_col_N)
        current_col_T += max(0, current_col_N - current_col_T)
        current_col_N += max(0, current_col_T - current_col_N)
    if current_col_T != current_col_N:
        logger.error("record: {}".format(record))
        raise(RuntimeError(
            "current_col_T != current_col_N"))

    del tumor_col_pos_map[max(tumor_col_pos_map.keys())]
    del normal_col_pos_map[max(normal_col_pos_map.keys())]
    new_tumor_matrix_ = np.zeros((5, current_col_T - 1))
    new_bq_tumor_matrix_ = np.zeros((5, current_col_T - 1))
    new_mq_tumor_matrix_ = np.zeros((5, current_col_T - 1))
    new_st_tumor_matrix_ = np.zeros((5, current_col_T - 1))
    new_lsc_tumor_matrix_ = np.zeros((5, current_col_T - 1))
    new_rsc_tumor_matrix_ = np.zeros((5, current_col_T - 1))
    new_tag_tumor_matrices_ = [
        np.zeros((5, current_col_T - 1)) for i in range(len(tag_tumor_matrices_))]
    new_tumor_matrix_[0, :] = max(tumor_matrix_.sum(0))
    new_bq_tumor_matrix_[0, :] = max(bq_tumor_matrix_[0, :])
    new_mq_tumor_matrix_[0, :] = max(mq_tumor_matrix_[0, :])
    new_st_tumor_matrix_[0, :] = max(st_tumor_matrix_[0, :])
    for iii in range(len(tag_tumor_matrices_)):
        new_tag_tumor_matrices_[iii][0, :] = max(
            tag_tumor_matrices_[iii][0, :])

    new_normal_matrix_ = np.zeros((5, current_col_N - 1))
    new_bq_normal_matrix_ = np.zeros((5, current_col_T - 1))
    new_mq_normal_matrix_ = np.zeros((5, current_col_T - 1))
    new_st_normal_matrix_ = np.zeros((5, current_col_T - 1))
    new_lsc_normal_matrix_ = np.zeros((5, current_col_T - 1))
    new_rsc_normal_matrix_ = np.zeros((5, current_col_T - 1))
    new_tag_normal_matrices_ = [
        np.zeros((5, current_col_T - 1)) for i in range(len(tag_normal_matrices_))]
    new_normal_matrix_[0, :] = max(normal_matrix_.sum(0))
    new_bq_normal_matrix_[0, :] = max(bq_normal_matrix_[0, :])
    new_mq_normal_matrix_[0, :] = max(mq_normal_matrix_[0, :])
    new_st_normal_matrix_[0, :] = max(st_normal_matrix_[0, :])
    for iii in range(len(tag_normal_matrices_)):
        new_tag_normal_matrices_[iii][0, :] = max(
            tag_normal_matrices_[iii][0, :])

    map_T = (np.cumsum(cols_T) - 1)[:-1]
    map_N = (np.cumsum(cols_N) - 1)[:-1]
    new_tumor_matrix_[:, map_T] = tumor_matrix_
    new_bq_tumor_matrix_[:, map_T] = bq_tumor_matrix_
    new_mq_tumor_matrix_[:, map_T] = mq_tumor_matrix_
    new_st_tumor_matrix_[:, map_T] = st_tumor_matrix_
    new_lsc_tumor_matrix_[:, map_T] = lsc_tumor_matrix_
    new_rsc_tumor_matrix_[:, map_T] = rsc_tumor_matrix_
    for iii in range(len(tag_tumor_matrices_)):
        new_tag_tumor_matrices_[iii][:, map_T] = tag_tumor_matrices_[iii]
    new_normal_matrix_[:, map_N] = normal_matrix_
    new_bq_normal_matrix_[:, map_N] = bq_normal_matrix_
    new_mq_normal_matrix_[:, map_N] = mq_normal_matrix_
    new_st_normal_matrix_[:, map_N] = st_normal_matrix_
    new_lsc_normal_matrix_[:, map_N] = lsc_normal_matrix_
    new_rsc_normal_matrix_[:, map_N] = rsc_normal_matrix_
    for iii in range(len(tag_normal_matrices_)):
        new_tag_normal_matrices_[iii][:, map_N] = tag_normal_matrices_[iii]

    new_tumor_ref_array = np.zeros((current_col_T - 1), int)
    new_normal_ref_array = np.zeros((current_col_N - 1), int)
    new_tumor_ref_array[map_T] = tumor_ref_array
    new_normal_ref_array[map_N] = normal_ref_array

    new_tumor_col_pos_map = {k: map_T[v]
                             for k, v in tumor_col_pos_map.iteritems()}
    new_normal_col_pos_map = {k: map_N[v]
                              for k, v in normal_col_pos_map.iteritems()}

    if sum(new_normal_ref_array - new_tumor_ref_array) != 0:
        logger.error("record: {}".format(record))
        logger.error("new_normal_ref_array, new_tumor_ref_array: {}, {}".format(
            new_normal_ref_array, new_tumor_ref_array))
        logger.error("new_normal_ref_array - new_tumor_ref_array: {}".format(
            new_normal_ref_array - new_tumor_ref_array))
        raise Exception

    for k in new_tumor_col_pos_map:
        assert(new_tumor_col_pos_map[k] == new_normal_col_pos_map[k])
    return [new_tumor_matrix_, new_bq_tumor_matrix_, new_mq_tumor_matrix_, new_st_tumor_matrix_,
            new_lsc_tumor_matrix_, new_rsc_tumor_matrix_, new_tag_tumor_matrices_, new_normal_matrix_,
            new_bq_normal_matrix_, new_mq_normal_matrix_, new_st_normal_matrix_,
            new_lsc_normal_matrix_, new_rsc_normal_matrix_,
            new_tag_normal_matrices_, new_tumor_ref_array,
            new_tumor_col_pos_map]


def prepare_info_matrices_tabix(ref_file, tumor_count_bed, normal_count_bed, record, rlen, rcenter, matrix_base_pad, matrix_width, min_ev_frac_per_col, min_cov):
    logger = logging.getLogger(prepare_info_matrices_tabix.__name__)

    chrom, pos, ref, alt = record[0:4]

    tumor_matrix_, bq_tumor_matrix_, mq_tumor_matrix_, st_tumor_matrix_, lsc_tumor_matrix_, rsc_tumor_matrix_, tag_tumor_matrices_, tumor_ref_array, tumor_col_pos_map = get_variant_matrix_tabix(
        ref_file, tumor_count_bed, record, matrix_base_pad)
    normal_matrix_, bq_normal_matrix_, mq_normal_matrix_, st_normal_matrix_, lsc_normal_matrix_, rsc_normal_matrix_, tag_normal_matrices_, normal_ref_array, normal_col_pos_map = get_variant_matrix_tabix(
        ref_file, normal_count_bed, record, matrix_base_pad)

    if not tumor_col_pos_map:
        logger.warning("Skip {} for all N reference".format(record))
        return None

    bq_tumor_matrix_[0, np.where(tumor_matrix_.sum(0) == 0)[
        0]] = np.max(bq_tumor_matrix_)
    bq_normal_matrix_[0, np.where(normal_matrix_.sum(0) == 0)[
        0]] = np.max(bq_normal_matrix_)
    mq_tumor_matrix_[0, np.where(tumor_matrix_.sum(0) == 0)[
        0]] = np.max(mq_tumor_matrix_)
    mq_normal_matrix_[0, np.where(normal_matrix_.sum(0) == 0)[
        0]] = np.max(mq_normal_matrix_)
    st_tumor_matrix_[0, np.where(tumor_matrix_.sum(0) == 0)[
        0]] = np.max(st_tumor_matrix_)
    st_normal_matrix_[0, np.where(normal_matrix_.sum(0) == 0)[
        0]] = np.max(st_normal_matrix_)
    lsc_tumor_matrix_[0, np.where(tumor_matrix_.sum(0) == 0)[
        0]] = np.max(lsc_tumor_matrix_)
    lsc_normal_matrix_[0, np.where(normal_matrix_.sum(0) == 0)[
        0]] = np.max(lsc_normal_matrix_)
    rsc_tumor_matrix_[0, np.where(tumor_matrix_.sum(0) == 0)[
        0]] = np.max(rsc_tumor_matrix_)
    rsc_normal_matrix_[0, np.where(normal_matrix_.sum(0) == 0)[
        0]] = np.max(rsc_normal_matrix_)
    for iii in range(len(tag_tumor_matrices_)):
        tag_tumor_matrices_[iii][0, np.where(tumor_matrix_.sum(0) == 0)[
            0]] = np.max(tag_tumor_matrices_[iii])
    for iii in range(len(tag_normal_matrices_)):
        tag_normal_matrices_[iii][0, np.where(normal_matrix_.sum(0) == 0)[
            0]] = np.max(tag_normal_matrices_[iii])

    tumor_matrix_[0, np.where(tumor_matrix_.sum(0) == 0)[
        0]] = max(np.sum(tumor_matrix_, 0))
    normal_matrix_[0, np.where(normal_matrix_.sum(0) == 0)[
        0]] = max(np.sum(normal_matrix_, 0))
    if max(np.sum(normal_matrix_, 0)) == 0:
        normal_matrix_[0, :] = np.max(np.sum(tumor_matrix_, 0))
        bq_normal_matrix_[0, :] = np.max(bq_tumor_matrix_)
        mq_normal_matrix_[0, :] = np.max(mq_tumor_matrix_)
        st_normal_matrix_[0, :] = np.max(st_tumor_matrix_)
        lsc_normal_matrix_[0, :] = np.max(lsc_tumor_matrix_)
        rsc_normal_matrix_[0, :] = np.max(rsc_tumor_matrix_)
        for iii in range(len(tag_normal_matrices_)):
            tag_normal_matrices_[iii][0, :] = np.max(tag_normal_matrices_[iii])

    tumor_matrix_, bq_tumor_matrix_, mq_tumor_matrix_, st_tumor_matrix_, lsc_tumor_matrix_, rsc_tumor_matrix_, \
        tag_tumor_matrices_, normal_matrix_, bq_normal_matrix_, mq_normal_matrix_, st_normal_matrix_, \
        lsc_normal_matrix_, rsc_normal_matrix_, tag_normal_matrices_, \
        ref_array, col_pos_map = align_tumor_normal_matrices(
            record, tumor_matrix_, bq_tumor_matrix_, mq_tumor_matrix_, st_tumor_matrix_, lsc_tumor_matrix_, rsc_tumor_matrix_,
            tag_tumor_matrices_, tumor_ref_array, tumor_col_pos_map, normal_matrix_,
            bq_normal_matrix_, mq_normal_matrix_, st_normal_matrix_, lsc_normal_matrix_, rsc_normal_matrix_, tag_normal_matrices_,
            normal_ref_array, normal_col_pos_map)

    tw = int(matrix_width)
    count_column = sum(tumor_matrix_[1:], 0)
    n_col = tumor_matrix_.shape[1]
    n_row = np.max(np.sum(tumor_matrix_, 0))

    if n_row < min_cov:
        logger.warning("Skip {} for low cov {}<{}".format(
            record, int(n_row), min_cov))
        return None
    cols_not_to_del = []
    largest_block = []
    if len(alt) > len(ref):
        z_ref_array = np.where((ref_array == 0))[0]
        if z_ref_array.shape[0] > 0:
            max_count = np.max(count_column[z_ref_array])
            cols_not_to_del = np.where(np.logical_and(count_column >= max_count * 0.7,
                                                      (ref_array == 0)))[0]
            a = [-1000] + sorted(np.where((ref_array == 0))[0]) + [1000]
            b = np.diff((np.diff(a) == 1).astype(int))
            a = a[1:-1]
            blocks = map(lambda x: [a[x[0]], a[x[1]]], zip(
                np.where(b == 1)[0], np.where(b == -1)[0]))
            if blocks:
                largest_block = sorted(
                    blocks, key=lambda x: x[1] - x[0] + 1)[-1]
                if np.max(count_column[range(largest_block[0], largest_block[1] + 1)]) > max_count * 0.05:
                    if (largest_block[1] - largest_block[0] + 1) > 2:
                        cols_not_to_del = sorted(list(set(cols_not_to_del) | set(
                            range(largest_block[0], largest_block[1] + 1))))
                else:
                    largest_block = []

    cols_to_del = sorted(list(set(np.where(np.logical_and(count_column <= (min_ev_frac_per_col * n_row), (ref_array == 0)))[0]
                                  ) - set(cols_not_to_del)))
    if n_col - len(cols_to_del) > tw:
        mn = min(count_column[np.argsort(count_column)][
                 n_col - tw - 1], max(n_row / 5, min_ev_frac_per_col * n_row))
        new_cols_to_del = set(np.where(np.logical_and(count_column <= mn, (ref_array == 0)))[
                              0]) - set(cols_to_del) - set(cols_not_to_del)
        if n_col - (len(cols_to_del) + len(new_cols_to_del)) < tw and len(new_cols_to_del) > 3 and len(alt) > len(ref):
            new_cols_to_del = map(
                lambda x: [count_column[x], x], new_cols_to_del)
            new_cols_to_del = sorted(new_cols_to_del, key=lambda x: [x[0], x[1]])[
                0:len(new_cols_to_del) - 4]
            new_cols_to_del = map(lambda x: x[1], new_cols_to_del)
        cols_to_del = sorted(set(list(new_cols_to_del) + list(cols_to_del)))
    cols_to_del = list(set(cols_to_del))
    if n_col - len(cols_to_del) > tw:
        for i in range(n_col / 10):
            if i not in cols_to_del and sum(rcenter) > -3:
                ref_b = ref_array[i]
                if tumor_matrix_[ref_b, i] > .8 * n_row:
                    cols_to_del.append(i)
            if n_col - len(cols_to_del) <= tw:
                break
            ii = n_col - i - 1
            if ii not in cols_to_del and sum(rcenter) < 3:
                ref_b = ref_array[ii]
                if tumor_matrix_[ref_b, ii] > .8 * n_row:
                    cols_to_del.append(ii)
            if n_col - len(cols_to_del) <= tw:
                break
    if n_col - len(cols_to_del) > tw and len(largest_block) > 0:
        block_len = (largest_block[1] - largest_block[0] + 1)
        if block_len > 2:
            cols_to_del = sorted(list(set(cols_to_del) |
                                      set(range(largest_block[0], largest_block[0] + min(n_col - len(cols_to_del) - tw, block_len - 3)))))
        else:
            cols_to_del = sorted(list(set(cols_to_del) | (set(np.where(count_column <= (min_ev_frac_per_col * n_row))[0]) &
                                                          set(range(largest_block[0], largest_block[1] + 1)))))

    cols_to_del.sort()
    for i, v in col_pos_map.iteritems():
        col_pos_map[i] -= sum(np.array(cols_to_del) < v)

    tumor_matrix = np.delete(tumor_matrix_, cols_to_del, 1)
    bq_tumor_matrix = np.delete(bq_tumor_matrix_, cols_to_del, 1)
    mq_tumor_matrix = np.delete(mq_tumor_matrix_, cols_to_del, 1)
    st_tumor_matrix = np.delete(st_tumor_matrix_, cols_to_del, 1)
    lsc_tumor_matrix = np.delete(lsc_tumor_matrix_, cols_to_del, 1)
    rsc_tumor_matrix = np.delete(rsc_tumor_matrix_, cols_to_del, 1)
    tag_tumor_matrices = []
    for iii in range(len(tag_tumor_matrices_)):
        tag_tumor_matrices.append(
            np.delete(tag_tumor_matrices_[iii], cols_to_del, 1))
    normal_matrix = np.delete(normal_matrix_, cols_to_del, 1)
    bq_normal_matrix = np.delete(bq_normal_matrix_, cols_to_del, 1)
    mq_normal_matrix = np.delete(mq_normal_matrix_, cols_to_del, 1)
    st_normal_matrix = np.delete(st_normal_matrix_, cols_to_del, 1)
    lsc_normal_matrix = np.delete(lsc_normal_matrix_, cols_to_del, 1)
    rsc_normal_matrix = np.delete(rsc_normal_matrix_, cols_to_del, 1)
    tag_normal_matrices = []
    for iii in range(len(tag_normal_matrices_)):
        tag_normal_matrices.append(
            np.delete(tag_normal_matrices_[iii], cols_to_del, 1))
    ref_array = np.delete(ref_array, cols_to_del, 0)

    ref_matrix = np.zeros(tumor_matrix.shape)
    ref_matrix[ref_array, range(ref_matrix.shape[1])] = np.sum(tumor_matrix, 0)

    ncols = tumor_matrix.shape[1]
    if matrix_width >= ncols:
        col_pos_map = {i: v + (matrix_width - ncols) /
                       2 for i, v in col_pos_map.iteritems()}
        tumor_count_matrix = np.zeros((5, matrix_width))
        tumor_count_matrix[:, (matrix_width - ncols) /
                           2:(matrix_width - ncols) / 2 + ncols] = tumor_matrix
        bq_tumor_count_matrix = np.zeros((5, matrix_width))
        bq_tumor_count_matrix[
            :, (matrix_width - ncols) / 2:(matrix_width - ncols) / 2 + ncols] = bq_tumor_matrix
        mq_tumor_count_matrix = np.zeros((5, matrix_width))
        mq_tumor_count_matrix[
            :, (matrix_width - ncols) / 2:(matrix_width - ncols) / 2 + ncols] = mq_tumor_matrix
        st_tumor_count_matrix = np.zeros((5, matrix_width))
        st_tumor_count_matrix[
            :, (matrix_width - ncols) / 2:(matrix_width - ncols) / 2 + ncols] = st_tumor_matrix
        lsc_tumor_count_matrix = np.zeros((5, matrix_width))
        lsc_tumor_count_matrix[
            :, (matrix_width - ncols) / 2:(matrix_width - ncols) / 2 + ncols] = lsc_tumor_matrix
        rsc_tumor_count_matrix = np.zeros((5, matrix_width))
        rsc_tumor_count_matrix[
            :, (matrix_width - ncols) / 2:(matrix_width - ncols) / 2 + ncols] = rsc_tumor_matrix
        tag_tumor_count_matrices = []
        for iii in range(len(tag_tumor_matrices)):
            tag_tumor_count_matrices.append(np.zeros((5, matrix_width)))
            tag_tumor_count_matrices[iii][
                :, (matrix_width - ncols) / 2:(matrix_width - ncols) / 2 + ncols] = tag_tumor_matrices[iii]
        normal_count_matrix = np.zeros((5, matrix_width))
        normal_count_matrix[
            :, (matrix_width - ncols) / 2:(matrix_width - ncols) / 2 + ncols] = normal_matrix
        bq_normal_count_matrix = np.zeros((5, matrix_width))
        bq_normal_count_matrix[
            :, (matrix_width - ncols) / 2:(matrix_width - ncols) / 2 + ncols] = bq_normal_matrix
        mq_normal_count_matrix = np.zeros((5, matrix_width))
        mq_normal_count_matrix[
            :, (matrix_width - ncols) / 2:(matrix_width - ncols) / 2 + ncols] = mq_normal_matrix
        st_normal_count_matrix = np.zeros((5, matrix_width))
        st_normal_count_matrix[
            :, (matrix_width - ncols) / 2:(matrix_width - ncols) / 2 + ncols] = st_normal_matrix
        lsc_normal_count_matrix = np.zeros((5, matrix_width))
        lsc_normal_count_matrix[
            :, (matrix_width - ncols) / 2:(matrix_width - ncols) / 2 + ncols] = lsc_normal_matrix
        rsc_normal_count_matrix = np.zeros((5, matrix_width))
        rsc_normal_count_matrix[
            :, (matrix_width - ncols) / 2:(matrix_width - ncols) / 2 + ncols] = rsc_normal_matrix
        tag_normal_count_matrices = []
        for iii in range(len(tag_normal_matrices)):
            tag_normal_count_matrices.append(np.zeros((5, matrix_width)))
            tag_normal_count_matrices[iii][
                :, (matrix_width - ncols) / 2:(matrix_width - ncols) / 2 + ncols] = tag_normal_matrices[iii]
        ref_count_matrix = np.zeros((5, matrix_width))
        ref_count_matrix[:, (matrix_width - ncols) /
                         2:(matrix_width - ncols) / 2 + ncols] = ref_matrix
    else:
        col_pos_map = {i: int(round(v / float(ncols) * matrix_width))
                       for i, v in col_pos_map.iteritems()}
        tumor_count_matrix = imresize(
            tumor_matrix, (5, matrix_width)).astype(int)
        bq_tumor_count_matrix = imresize(
            bq_tumor_matrix, (5, matrix_width)).astype(int)
        mq_tumor_count_matrix = imresize(
            mq_tumor_matrix, (5, matrix_width)).astype(int)
        st_tumor_count_matrix = imresize(
            st_tumor_matrix, (5, matrix_width)).astype(int)
        lsc_tumor_count_matrix = imresize(
            lsc_tumor_matrix, (5, matrix_width)).astype(int)
        rsc_tumor_count_matrix = imresize(
            rsc_tumor_matrix, (5, matrix_width)).astype(int)
        tag_tumor_count_matrices = []
        for iii in range(len(tag_tumor_matrices)):
            tag_tumor_count_matrices.append(
                imresize(tag_tumor_matrices[iii], (5, matrix_width)).astype(int))
        normal_count_matrix = imresize(
            normal_matrix, (5, matrix_width)).astype(int)
        bq_normal_count_matrix = imresize(
            bq_normal_matrix, (5, matrix_width)).astype(int)
        mq_normal_count_matrix = imresize(
            mq_normal_matrix, (5, matrix_width)).astype(int)
        st_normal_count_matrix = imresize(
            st_normal_matrix, (5, matrix_width)).astype(int)
        lsc_normal_count_matrix = imresize(
            lsc_normal_matrix, (5, matrix_width)).astype(int)
        rsc_normal_count_matrix = imresize(
            rsc_normal_matrix, (5, matrix_width)).astype(int)
        tag_normal_count_matrices = []
        for iii in range(len(tag_normal_matrices)):
            tag_normal_count_matrices.append(
                imresize(tag_normal_matrices[iii], (5, matrix_width)).astype(int))
        ref_count_matrix = imresize(ref_matrix, (5, matrix_width)).astype(int)

    if int(pos) + rcenter[0] not in col_pos_map:
        center = min(col_pos_map.values()) + rcenter[0] - 1 + rcenter[1]
    else:
        center = col_pos_map[int(pos) + rcenter[0]] + rcenter[1]

    if center > ref_count_matrix.shape[1] - 1:
        center = min(max(0, min(col_pos_map.values()) +
                         rcenter[0] - 1 + rcenter[1]), ref_count_matrix.shape[1] - 1)

    return [tumor_matrix_, tumor_matrix, normal_matrix_, normal_matrix, ref_count_matrix, tumor_count_matrix,
            bq_tumor_count_matrix, mq_tumor_count_matrix, st_tumor_count_matrix, lsc_tumor_count_matrix, rsc_tumor_count_matrix,
            tag_tumor_count_matrices, normal_count_matrix, bq_normal_count_matrix, mq_normal_count_matrix,
            st_normal_count_matrix, lsc_normal_count_matrix, rsc_normal_count_matrix,
            tag_normal_count_matrices, center, rlen, col_pos_map]


def prep_data_single_tabix((ref_file, tumor_count_bed, normal_count_bed, record, vartype, rlen, rcenter, ch_order,
                            matrix_base_pad, matrix_width, min_ev_frac_per_col, min_cov, ann)):
    thread_logger = logging.getLogger(
        "{} ({})".format(prep_data_single_tabix.__name__, multiprocessing.current_process().name))
    try:
        chrom, pos, ref, alt = record[:4]
        pos = int(pos)
        matrices_info = prepare_info_matrices_tabix(ref_file=ref_file,
                                                    tumor_count_bed=tumor_count_bed,
                                                    normal_count_bed=normal_count_bed, record=record, rlen=rlen, rcenter=rcenter,
                                                    matrix_base_pad=matrix_base_pad, matrix_width=matrix_width,
                                                    min_ev_frac_per_col=min_ev_frac_per_col,
                                                    min_cov=min_cov)
        if matrices_info:
            tumor_matrix_, tumor_matrix, normal_matrix_, normal_matrix, ref_count_matrix, tumor_count_matrix, \
                bq_tumor_count_matrix, mq_tumor_count_matrix, st_tumor_count_matrix, lsc_tumor_count_matrix, rsc_tumor_count_matrix, \
                tag_tumor_count_matrices, normal_count_matrix, bq_normal_count_matrix, mq_normal_count_matrix, st_normal_count_matrix, \
                lsc_normal_count_matrix, rsc_normal_count_matrix, tag_normal_count_matrices, center, rlen, col_pos_map = matrices_info
        else:
            return []

        candidate_mat = np.zeros((tumor_count_matrix.shape[0], tumor_count_matrix.shape[
                                 1], 13 + (len(tag_tumor_count_matrices) * 2)))
        candidate_mat[:, :, 0] = ref_count_matrix
        candidate_mat[:, :, 1] = tumor_count_matrix
        candidate_mat[:, :, 2] = normal_count_matrix
        candidate_mat[:, :, 3] = bq_tumor_count_matrix
        candidate_mat[:, :, 4] = bq_normal_count_matrix
        candidate_mat[:, :, 5] = mq_tumor_count_matrix
        candidate_mat[:, :, 6] = mq_normal_count_matrix
        candidate_mat[:, :, 7] = st_tumor_count_matrix
        candidate_mat[:, :, 8] = st_normal_count_matrix
        candidate_mat[:, :, 9] = lsc_tumor_count_matrix
        candidate_mat[:, :, 10] = lsc_normal_count_matrix
        candidate_mat[:, :, 11] = rsc_tumor_count_matrix
        candidate_mat[:, :, 12] = rsc_normal_count_matrix
        for iii in range(len(tag_tumor_count_matrices)):
            candidate_mat[:, :, 13 + (iii * 2)] = tag_tumor_count_matrices[iii]
            candidate_mat[:, :, 13 + (iii * 2) +
                          1] = tag_normal_count_matrices[iii]
        tumor_cov = int(round(max(np.sum(tumor_count_matrix, 0))))
        normal_cov = int(round(max(np.sum(normal_count_matrix, 0))))

        candidate_mat[:, :, 0] = candidate_mat[
            :, :, 0] / (max(np.max(ref_count_matrix), np.max(tumor_count_matrix)) + 0.00001) * 255
        candidate_mat[:, :, 1] = candidate_mat[:, :, 1] / \
            (np.max(tumor_count_matrix) + 0.00001) * 255
        candidate_mat[:, :, 2] = candidate_mat[:, :, 2] / \
            (np.max(normal_count_matrix) + 0.00001) * 255
        candidate_mat[:, :, 3] = candidate_mat[:, :, 3] / \
            (max(np.max(bq_tumor_count_matrix), 41.0)) * 255
        candidate_mat[:, :, 4] = candidate_mat[:, :, 4] / \
            (max(np.max(bq_normal_count_matrix), 41.0)) * 255
        candidate_mat[:, :, 5] = candidate_mat[:, :, 5] / \
            (max(np.max(mq_tumor_count_matrix), 70.0)) * 255
        candidate_mat[:, :, 6] = candidate_mat[:, :, 6] / \
            (max(np.max(mq_normal_count_matrix), 70.0)) * 255
        candidate_mat[:, :, 7] = candidate_mat[:, :, 7] / \
            (np.max(tumor_count_matrix) + 0.00001) * 255
        candidate_mat[:, :, 8] = candidate_mat[:, :, 8] / \
            (np.max(normal_count_matrix) + 0.00001) * 255
        candidate_mat[:, :, 9] = candidate_mat[:, :, 9] / \
            (np.max(tumor_count_matrix) + 0.00001) * 255
        candidate_mat[:, :, 10] = candidate_mat[:, :, 10] / \
            (np.max(normal_count_matrix) + 0.00001) * 255
        candidate_mat[:, :, 11] = candidate_mat[:, :, 11] / \
            (np.max(tumor_count_matrix) + 0.00001) * 255
        candidate_mat[:, :, 12] = candidate_mat[:, :, 12] / \
            (np.max(normal_count_matrix) + 0.00001) * 255
        for iii in range(len(tag_tumor_count_matrices)):
            candidate_mat[:, :, 13 + (iii * 2)] = candidate_mat[:, :, 13 + (iii * 2)] / (
                max(np.max(tag_tumor_count_matrices[iii]), 100.0)) * 255
            candidate_mat[:, :, 13 + (iii * 2) + 1] = candidate_mat[:, :, 13 + (
                iii * 2) + 1] / (max(np.max(tag_normal_count_matrices[iii]), 100.0)) * 255

        candidate_mat = np.maximum(0, np.minimum(
            candidate_mat, 255)).astype(np.uint8)
        tag = "{}.{}.{}.{}.{}.{}.{}.{}.{}".format(ch_order, pos, ref[0:55], alt[
                                                  0:55], vartype, center, rlen, tumor_cov, normal_cov)
        candidate_mat = base64.b64encode(zlib.compress(candidate_mat))
        return tag, candidate_mat, record[0:4], ann
    except Exception as ex:
        thread_logger.error(traceback.format_exc())
        thread_logger.error(ex)
        return None


def push_lr(fasta_file, record, left_right_both):
    logger = logging.getLogger(push_lr.__name__)
    record[0] = str(record[0])
    eqs = [record]
    if "," not in record[3]:
        if record[2] != record[3] and left_right_both in [0, 2]:
            chrom, pos, ref, alt = record[0:4]
            new_pos = pos
            new_ref = ref
            new_alt = alt
            while(new_pos > 1):
                l_base = fasta_file.fetch(
                    (chrom), new_pos - 2, new_pos - 1).upper()
                new_ref = l_base + new_ref
                new_alt = l_base + new_alt
                new_pos -= 1
                while(len(new_alt) > 1 and len(new_ref) > 1):
                    if new_alt[-1] == new_ref[-1]:
                        new_alt = new_alt[:-1]
                        new_ref = new_ref[:-1]
                    else:
                        break
                if len(new_alt) > len(alt):
                    new_ref = new_ref[1:]
                    new_alt = new_alt[1:]
                    new_pos += 1
                    break
                else:
                    eqs.append([chrom, new_pos, new_ref, new_alt] + record[4:])
            record = [chrom, new_pos, new_ref, new_alt] + record[4:]
        if record[2] != record[3] and left_right_both:
            chrom, pos, ref, alt = record[0:4]
            new_pos = pos + len(ref)
            new_ref = ref
            new_alt = alt
            max_pos = fasta_file.lengths[fasta_file.references.index(chrom)]
            while(new_pos < max_pos):
                r_base = fasta_file.fetch(
                    (chrom), new_pos - 1, new_pos).upper()
                new_ref = new_ref + r_base
                new_alt = new_alt + r_base
                new_pos += 1
                while(len(new_alt) > 1 and len(new_ref) > 1):
                    if new_alt[0] == new_ref[0] and new_alt[1] == new_ref[1]:
                        new_alt = new_alt[1:]
                        new_ref = new_ref[1:]
                    else:
                        break
                if len(new_alt) > len(alt):
                    new_ref = new_ref[:-1]
                    new_alt = new_alt[:-1]
                    new_pos -= 1
                    break
                else:
                    eqs.append([chrom, new_pos - len(new_ref),
                                new_ref, new_alt] + record[4:])
            record = [chrom, new_pos -
                      len(new_ref), new_ref, new_alt] + record[4:]
        eqs = map(lambda x: eqs[x], dict(
            map(lambda x: ["_".join(map(str, x[1])), x[0]], enumerate(eqs))).values())

        for eq in eqs:
            c, p, r, a = eq[0:4]
            assert(fasta_file.fetch((c), p - 1, p - 1 + len(r)).upper() == r)

    return record, eqs


def merge_records(fasta_file, records):
    logger = logging.getLogger(merge_records.__name__)
    if len(set(map(lambda x: x[0], records))) != 1:
        return None
    if len(records) == 1:
        return records[0][0:4]
    chrom = records[0][0]
    pos_s = map(lambda x: x[1], records)
    pos_m = min(pos_s) - 1
    b = pos_m
    ref2_ = ""
    alt2_ = ""
    for record in sorted(records, key=lambda x: x[1]):
        chrom_, pos_, ref_, alt_ = record[0:4]
        if len(ref_) == len(alt_):
            if b > pos_ - 1:
                return None
            sequence = fasta_file.fetch(chrom, b, pos_ - 1).upper()
            ref2_ += sequence + ref_
            alt2_ += sequence + alt_
        else:
            if b <= (pos_ - 1):
                sequence = fasta_file.fetch(chrom, b, pos_ - 1).upper()
                ref2_ += sequence + ref_
                alt2_ += sequence + alt_
            elif b == (pos_):
                sequence = fasta_file.fetch(chrom, b, pos_).upper()
                ref2_ += sequence + ref_[1:]
                alt2_ += sequence + alt_[1:]
            else:
                return None

        b = pos_ - 1 + len(ref_)
    while(len(alt2_) > 1 and len(ref2_) > 1):
        if alt2_[0] == ref2_[0] and (alt2_[1] == ref2_[1] or len(alt2_) == len(ref2_)):
            alt2_ = alt2_[1:]
            ref2_ = ref2_[1:]
            pos_m += 1
        else:
            break
    while(len(alt2_) > 1 and len(ref2_) > 1):
        if alt2_[-1] == ref2_[-1]:
            alt2_ = alt2_[:-1]
            ref2_ = ref2_[:-1]
        else:
            break

    return [str(chrom), pos_m + 1, ref2_, alt2_]


def is_part_of(record1, record2):
    logger = logging.getLogger(is_part_of.__name__)
    chrom1, pos1, ref1, alt1 = record1[0:4]
    chrom2, pos2, ref2, alt2 = record2[0:4]
    if chrom1 != chrom2:
        return False
    vartype1 = get_type(ref1, alt1)
    vartype2 = get_type(ref2, alt2)
    if vartype1 == "SNP" and vartype2 == "DEL":
        if pos2 < pos1 < pos2 + len(ref2):
            return True
    elif vartype2 == "SNP" and vartype1 == "DEL":
        if pos1 < pos2 < pos1 + len(ref1):
            return True
    elif vartype1 == vartype2:
        if pos1 == pos2:
            return True
        elif vartype1 == "DEL" and set(range(pos1 + 1, pos1 + len(ref1))) & set(range(pos2 + 1, pos2 + len(ref2))):
            return True
    return False


def find_i_center(ref, alt):
    logger = logging.getLogger(find_i_center.__name__)
    i_center = 0
    if len(alt) != len(ref):
        while(min(len(alt), len(ref)) > i_center and alt[i_center] == ref[i_center]):
            i_center += 1
    return [0, i_center] if (len(ref) < len(alt)) else [i_center, 0]


def find_len(ref, alt):
    logger = logging.getLogger(find_len.__name__)
    i_ = 0
    while(min(len(alt), len(ref)) > i_ and alt[i_] == ref[i_]):
        i_ += 1
    ref_ = ref[i_:]
    alt_ = alt[i_:]
    i_ = 0
    while (min(len(alt_), len(ref_)) > i_ and alt_[len(alt_) - i_ - 1] == ref_[len(ref_) - i_ - 1]):
        i_ += 1
    if i_ > 0:
        ref_ = ref_[:-i_]
        alt_ = alt_[:-i_]
    return max(len(ref_), len(alt_))


def find_records((work, split_region_file, truth_vcf_file, pred_vcf_file, ref_file, ensemble_bed, work_index)):
    thread_logger = logging.getLogger(
        "{} ({})".format(find_records.__name__, multiprocessing.current_process().name))
    try:
        thread_logger.info(
            "Start find_records for worker {}".format(work_index))

        split_bed = pybedtools.BedTool(
            split_region_file).slop(g=ref_file + ".fai", b=5)
        split_truth_vcf_file = os.path.join(
            work, "truth_{}.vcf".format(work_index))
        split_pred_vcf_file = os.path.join(
            work, "pred_{}.vcf".format(work_index))
        split_ensemble_bed_file = os.path.join(
            work, "ensemble_{}.bed".format(work_index))
        split_missed_ensemble_bed_file = os.path.join(
            work, "missed_ensemble_{}.bed".format(work_index))
        split_pred_with_missed_file = os.path.join(
            work, "pred_with_missed_{}.bed".format(work_index))
        split_in_ensemble_bed = os.path.join(
            work, "in_ensemble_{}.bed".format(work_index))
        pybedtools.BedTool(truth_vcf_file).intersect(
            split_bed, u=True).saveas(split_truth_vcf_file)
        pybedtools.BedTool(pred_vcf_file).intersect(
            split_bed, u=True).saveas(split_pred_vcf_file)
        if ensemble_bed:
            pybedtools.BedTool(ensemble_bed).intersect(
                split_bed, u=True).saveas(split_ensemble_bed_file)
            pybedtools.BedTool(split_ensemble_bed_file).window(split_pred_vcf_file, w=5, v=True).each(
                lambda x: pybedtools.Interval(x[0], int(x[1]), int(x[1]) + 1, x[3], x[4], ".", otherfields=[".", ".", ".", "."])).saveas(split_missed_ensemble_bed_file)
            concatenate_vcfs(
                [split_pred_vcf_file, split_missed_ensemble_bed_file], split_pred_with_missed_file)
            pred_with_missed = pybedtools.BedTool(split_pred_with_missed_file).each(
                lambda x: pybedtools.create_interval_from_list(
                    [x[0], str(x[1]), ".", x[3], x[4], ".", ".", ".", ".", "."])).sort().saveas(
                split_pred_with_missed_file)
            not_in_ensemble_bed = pybedtools.BedTool(
                pred_with_missed).window(split_ensemble_bed_file, w=1, v=True)
            in_ensemble_bed = pybedtools.BedTool(pred_with_missed).window(
                split_ensemble_bed_file, w=1).saveas(split_in_ensemble_bed)

        records = []
        i = 0
        anns = {}
        if ensemble_bed:
            for record in not_in_ensemble_bed:
                chrom, pos, ref, alt = [str(record[0]), int(
                    record[1]), record[3], record[4]]
                r_ = []
                if len(ref) == len(alt) and len(ref) > 1:
                    for ii in range(len(ref)):
                        ref_ = ref[ii]
                        alt_ = alt[ii]
                        if ref_ != alt_:
                            r_.append([chrom, pos + ii, ref_, alt_])
                else:
                    r_ = [[chrom, pos, ref, alt]]
                for rr in r_:
                    records.append(rr + [str(i)])
                    anns[i] = [0] * 93
                    i += 1

            curren_pos_records = []
            emit_flag = False
            len_in_ensemble_bed = len(in_ensemble_bed)
            for j in range(len_in_ensemble_bed + 1):

                if j == (len_in_ensemble_bed):
                    emit_flag = True
                else:
                    record = in_ensemble_bed[j]
                    if curren_pos_records:
                        if record[0] == curren_pos_records[0][0] and record[1] == curren_pos_records[0][1]  \
                                and record[3] == curren_pos_records[0][3] and record[4] == curren_pos_records[0][4]:
                            curren_pos_records.append(record)
                        else:
                            emit_flag = True
                    else:
                        curren_pos_records.append(record)

                if emit_flag:
                    if curren_pos_records:
                        rrs = []
                        for record_ in curren_pos_records:
                            chrom, pos, ref, alt = [str(record_[0]), int(
                                record_[1]), record_[3], record_[4]]
                            r_ = []
                            if len(ref) == len(alt) and len(ref) > 1:
                                for ii in range(len(ref)):
                                    ref_ = ref[ii]
                                    alt_ = alt[ii]
                                    if ref_ != alt_:
                                        r_.append(
                                            [chrom, pos + ii, ref_, alt_])
                            else:
                                r_ = [[chrom, pos, ref, alt]]

                            ann = [0] * 93
                            if record_[1] == record_[11]:
                                if record_[3] == record_[13] and record_[4] == record_[14]:
                                    ann = record_[15:]
                                elif len(record_[3]) > len(record_[4]) and len(record_[13]) > len(record_[14]) and \
                                        (record_[4]) == (record_[14]):
                                    if (len(record_[3]) > len(record_[13]) and record_[3][0:len(record_[13])] == record_[13]) or (
                                            len(record_[13]) > len(record_[3]) and record_[13][0:len(record_[3])] == record_[3]):
                                        ann = record_[15:]
                                elif len(record_[3]) < len(record_[4]) and len(record_[13]) < len(record_[14]) and \
                                        (record_[3]) == (record_[13]):
                                    if (len(record_[4]) > len(record_[14]) and record_[4][0:len(record_[14])] == record_[14]) or (
                                            len(record_[14]) > len(record_[4]) and record_[14][0:len(record_[4])] == record_[4]):
                                        ann = record_[15:]
                            if ann:
                                ann = map(float, ann)
                            rrs.append([r_, ann])
                        max_ann = max(map(lambda x: sum(x[1]), rrs))
                        if max_ann > 0:
                            rrs = filter(lambda x: sum(x[1]) > 0, rrs)
                        elif max_ann == 0:
                            rrs = rrs[0:1]
                        for r_, ann in rrs:
                            for rr in r_:
                                records.append(rr + [str(i)])
                                anns[i] = ann
                                i += 1
                    emit_flag = False
                    if j < len_in_ensemble_bed:
                        curren_pos_records = [record]
        else:
            with open(split_pred_vcf_file, 'r') as vcf_reader:
                for line in vcf_reader:
                    if line[0] == "#":
                        continue
                    record = line.strip().split()
                    chrom, pos, ref, alt = [record[0], int(
                        record[1]), record[3], record[4]]
                    r_ = []
                    if len(ref) == len(alt) and len(ref) > 1:
                        for ii in range(len(ref)):
                            ref_ = ref[ii]
                            alt_ = alt[ii]
                            if ref_ != alt_:
                                r_.append([chrom, pos + ii, ref_, alt_])
                    else:
                        r_ = [[chrom, pos, ref, alt]]

                    for rr in r_:
                        records.append(rr + [str(i)])
                        i += 1
        records_bed = pybedtools.BedTool(map(lambda x: pybedtools.Interval(
            x[0], x[1], x[1] + len(x[2]), x[2], x[3], x[4]), records))

        truth_records = []
        i = 0
        with open(split_truth_vcf_file, 'r') as vcf_reader:
            for line in vcf_reader:
                if line[0] == "#":
                    continue
                record = line.strip().split()
                truth_records.append(
                    [record[0], int(record[1]), record[3], record[4], str(i)])
                i += 1

        truth_bed = pybedtools.BedTool(map(lambda x: pybedtools.Interval(
            x[0], x[1], x[1] + len(x[2]), x[2], x[3], x[4]), truth_records))
        none_records_0 = records_bed.window(truth_bed, w=5, v=True)
        none_records_ids = map(lambda x: int(x[5]), none_records_0)
        other_records = records_bed.window(truth_bed, w=5)
        map_pred_2_truth = {}
        map_truth_2_pred = {}
        for record in other_records:
            id_pred = int(record[5])
            id_truth = int(record[11])
            if id_pred not in map_pred_2_truth:
                map_pred_2_truth[id_pred] = []
            map_pred_2_truth[id_pred].append(id_truth)
            if id_truth not in map_truth_2_pred:
                map_truth_2_pred[id_truth] = []
            map_truth_2_pred[id_truth].append(id_pred)

        record_center = {}

        chroms_order = get_chromosomes_order(reference=ref_file)
        fasta_file = pysam.Fastafile(ref_file)

        good_records = {"INS": [], "DEL": [], "SNP": []}
        vtype = {}
        record_len = {}
        for i, js in map_truth_2_pred.iteritems():
            truth_record = truth_records[i]
            for j in js:
                record = records[j]
                if record[0] == truth_record[0] and record[1] == truth_record[1] and record[2] == truth_record[2] and record[3] == truth_record[3]:
                    assert int(record[4]) == j
                    vartype = get_type(record[2], record[3])
                    if j not in good_records[vartype]:
                        ref, alt = truth_record[2:4]
                        record_center[j] = find_i_center(ref, alt)
                        record_len[j] = find_len(ref, alt)
                        good_records[vartype].append(j)
                        vtype[j] = vartype

        good_records_idx = [i for w in good_records.values() for i in w]
        remained_idx = sorted(set(range(len(records))) -
                              (set(good_records_idx) | set(none_records_ids)))
        done_js = list(good_records_idx)
        for i, js in map_truth_2_pred.iteritems():
            truth_record = truth_records[i]
            if set(js) & set(good_records_idx) or set(js) & set(done_js):
                continue
            i_s = sorted(set([ii for j in js for ii in map_pred_2_truth[j]]))
            done_is_ = []
            done_js_ = []
            done = False
            for idx_ii, ii in enumerate(i_s):
                if done:
                    break
                for n_merge_i in range(0, len(i_s) - idx_ii):
                    if done:
                        break
                    t_i = i_s[idx_ii:idx_ii + n_merge_i + 1]
                    t_ = [truth_records[iii] for iii in t_i]
                    mt = merge_records(fasta_file, t_)
                    if mt:
                        mt2, eqs2 = push_lr(fasta_file, mt, 2)
                        eqs2 = map(lambda x: "_".join(map(str, x[0:4])), eqs2)
                        for idx_jj, jj in enumerate(js):
                            for n_merge_j in range(0, len(js) - idx_jj):
                                r_j = js[idx_jj:idx_jj + n_merge_j + 1]
                                if set(r_j) & set(done_js_):
                                    continue
                                r_ = [records[jjj] for jjj in r_j]
                                mr = merge_records(fasta_file, r_)
                                if mr:
                                    mr2, eqs_r2 = push_lr(fasta_file, mr, 2)
                                    record_str = "_".join(map(str, mr2[0:4]))
                                    if record_str in eqs2:
                                        for j in r_j:
                                            record = records[j]
                                            vartype = get_type(
                                                record[2], record[3])
                                            pos, ref, alt = record[1:4]
                                            record_center[
                                                j] = find_i_center(ref, alt)
                                            record_len[j] = find_len(ref, alt)
                                            good_records[vartype].append(j)
                                            vtype[j] = vartype
                                            done_js.append(j)
                                            done_js_.append(j)
                                        done_is_.extend(t_i)
                                        if not set(js) - set(done_js_):
                                            done = True

        perfect_idx = [i for w in good_records.values() for i in w]
        good_records_idx = [i for w in good_records.values() for i in w]
        remained_idx = sorted(set(range(len(records))) -
                              (set(good_records_idx) | set(none_records_ids)))
        for j in remained_idx:
            record = records[j]
            pos, ref, alt = record[1:4]
            i_s = map_pred_2_truth[j]
            done = False
            for i in i_s:
                truth_record = truth_records[i]
                tr, eqs = push_lr(fasta_file, truth_record, 2)
                for eq in eqs:
                    if is_part_of(eq, record):
                        ref_t, alt_t = truth_record[2:4]
                        vartype_t = get_type(ref_t, alt_t)
                        record_center[j] = find_i_center(ref, alt)
                        record_len[j] = find_len(ref_t, alt_t)
                        good_records[vartype_t].append(j)
                        vtype[j] = vartype_t
                        done = True
                        break
                if done:
                    break
            if not done:
                p_s = set([jj for i in i_s for jj in map_truth_2_pred[i]]) & set(
                    perfect_idx)
                for p in p_s:
                    ref_p, alt_p = records[p][2:4]
                    tr, eqs = push_lr(fasta_file, records[p], 2)
                    for eq in eqs:
                        if is_part_of(eq, record):
                            vartype = vtype[p]
                            record_center[j] = find_i_center(ref, alt)
                            record_len[j] = find_len(ref_p, alt_p)
                            good_records[vartype].append(j)
                            vtype[j] = vartype
                            done = True
                            break
                    if done:
                        break

        good_records_idx = [i for w in good_records.values() for i in w]
        remained_idx = sorted(set(range(len(records))) -
                              (set(good_records_idx) | set(none_records_ids)))
        for i, js in map_truth_2_pred.iteritems():
            truth_record = truth_records[i]
            if set(js) & set(good_records_idx):
                continue
            pos_t, ref_t, alt_t = truth_record[1:4]
            vartype_t = get_type(ref_t, alt_t)
            rct = find_i_center(ref_t, alt_t)
            if len(js) == 2 and vartype_t == "SNP":
                vartype0 = get_type(records[js[0]][2], records[js[0]][3])
                vartype1 = get_type(records[js[1]][2], records[js[1]][3])
                if (vartype0 == "DEL" and vartype1 == "INS") or (vartype1 == "DEL" and vartype0 == "INS"):
                    for j in js:
                        record = records[j]
                        pos, ref, alt = record[1:4]
                        vartype = get_type(ref, alt)
                        record_center[j] = find_i_center(ref, alt)
                        record_len[j] = find_len(ref, alt)
                        good_records[vartype].append(j)
                        vtype[j] = vartype

        good_records_idx = [i for w in good_records.values() for i in w]
        remained_idx = sorted(set(range(len(records))) -
                              (set(good_records_idx) | set(none_records_ids)))
        for i, js in map_truth_2_pred.iteritems():
            truth_record = truth_records[i]
            if set(js) & set(good_records_idx):
                continue
            pos_t, ref_t, alt_t = truth_record[1:4]
            vartype_t = get_type(ref_t, alt_t)
            rct = find_i_center(ref_t, alt_t)
            for j in js:
                record = records[j]
                vartype = get_type(record[2], record[3])
                pos, ref, alt = record[1:4]
                rc = find_i_center(ref, alt)
                if vartype_t == vartype and pos_t == pos:
                    good_records[vartype_t].append(j)
                    vtype[j] = vartype_t
                    record_len[j] = find_len(ref_t, alt_t)
                    record_center[j] = rc

        good_records_idx = [i for w in good_records.values() for i in w]
        remained_idx = sorted(set(range(len(records))) -
                              (set(good_records_idx) | set(none_records_ids)))
        for i, js in map_truth_2_pred.iteritems():
            truth_record = truth_records[i]

            if set(js) & set(good_records_idx):
                continue
            pos_t, ref_t, alt_t = truth_record[1:4]
            vartype_t = get_type(ref_t, alt_t)
            rct = find_i_center(ref_t, alt_t)
            for j in js:
                if j not in remained_idx:
                    continue
                record = records[j]
                vartype = get_type(record[2], record[3])
                pos, ref, alt = record[1:4]
                rc = find_i_center(ref, alt)
                if pos_t + rct[0] + rct[1] == pos + rc[0] + rc[1]:
                    if (vartype_t == "INS" and vartype == "SNP") or (vartype == "INS" and vartype_t == "SNP"):
                        good_records[vartype_t].append(j)
                        vtype[j] = vartype_t
                        record_len[j] = find_len(ref_t, alt_t)
                        record_center[j] = rc

        good_records_idx = [i for w in good_records.values() for i in w]
        remained_idx = sorted(set(range(len(records))) -
                              (set(good_records_idx) | set(none_records_ids)))
        for i, js in map_truth_2_pred.iteritems():
            truth_record = truth_records[i]
            if set(js) & set(good_records_idx):
                continue
            pos_t, ref_t, alt_t = truth_record[1:4]
            vartype_t = get_type(ref_t, alt_t)
            for j in js:
                record = records[j]
                pos, ref, alt = record[1:4]
                vartype = get_type(record[2], record[3])
                if (vartype == vartype_t) and vartype_t != "SNP" and abs(pos - pos_t) < 2:
                    good_records[vartype_t].append(j)
                    vtype[j] = vartype_t
                    record_center[j] = find_i_center(ref, alt)
                    record_len[j] = find_len(ref_t, alt_t)

        good_records_idx = [i for w in good_records.values() for i in w]
        remained_idx = sorted(set(range(len(records))) -
                              (set(good_records_idx) | set(none_records_ids)))
        for i, js in map_truth_2_pred.iteritems():
            truth_record = truth_records[i]

            if set(js) & set(good_records_idx):
                continue

        good_records_idx = [i for w in good_records.values() for i in w]
        remained_idx = sorted(set(range(len(records))) -
                              (set(good_records_idx) | set(none_records_ids)))
        for j in remained_idx:
            none_records_ids.append(j)

        for j in none_records_ids:
            record = records[j]
            pos, ref, alt = record[1:4]
            record_center[j] = find_i_center(ref, alt)
            vtype[j] = "NONE"
            record_len[j] = 0

        good_records_idx = [i for w in good_records.values() for i in w]

        records_r = [records[x] for k, w in good_records.items() for x in w]
        records_r_bed = pybedtools.BedTool(map(lambda x: pybedtools.Interval(
            x[0], x[1], x[1] + len(x[2]), x[2], x[3], x[4]), records_r))

        N_none = len(none_records_ids)
        thread_logger.info("N_none: {} ".format(N_none))
        none_records = map(lambda x: records[x], none_records_ids)
        none_records = sorted(none_records, key=lambda x: [x[0], int(x[1])])

        return records_r, none_records, vtype, record_len, record_center, chroms_order, anns
    except Exception as ex:
        thread_logger.error(traceback.format_exc())
        thread_logger.error(ex)
        return None


def extract_ensemble(work, ensemble_tsv):
    logger = logging.getLogger(extract_ensemble.__name__)
    ensemble_data = []
    ensemble_pos = []
    header = []
    header_pos = []
    COV = 50
    with open(ensemble_tsv) as s_f:
        for line in s_f:
            if line[0:5] == "CHROM":
                header_pos = line.strip().split()[0:5]
                header = line.strip().split()[5:105]
                continue
            fields = line.strip().split()
            fields[2] = str(int(fields[1]) + len(fields[3]))
            ensemble_pos.append(fields[0:5])
            ensemble_data.append(map(lambda x: float(
                x.replace("False", "0").replace("True", "1")), fields[5:105]))
    ensemble_data = np.array(ensemble_data)

    cov_features = map(lambda x: x[0], filter(lambda x: x[1] in [
        "Consistent_Mates", "Inconsistent_Mates", "N_DP",
        "nBAM_REF_NM", "nBAM_ALT_NM", "nBAM_REF_Concordant", "nBAM_REF_Discordant", "nBAM_ALT_Concordant", "nBAM_ALT_Discordant",
        "N_REF_FOR", "N_REF_REV", "N_ALT_FOR", "N_ALT_REV", "nBAM_REF_Clipped_Reads", "nBAM_ALT_Clipped_Reads",  "nBAM_MQ0", "nBAM_Other_Reads", "nBAM_Poor_Reads",
        "nBAM_REF_InDel_3bp", "nBAM_REF_InDel_2bp", "nBAM_REF_InDel_1bp", "nBAM_ALT_InDel_3bp", "nBAM_ALT_InDel_2bp",
        "nBAM_ALT_InDel_1bp",
        "T_DP", "tBAM_REF_NM", "tBAM_ALT_NM", "tBAM_REF_Concordant", "tBAM_REF_Discordant", "tBAM_ALT_Concordant", "tBAM_ALT_Discordant",
        "T_REF_FOR", "T_REF_REV", "T_ALT_FOR", "T_ALT_REV",
        "tBAM_REF_Clipped_Reads", "tBAM_ALT_Clipped_Reads",
        "tBAM_MQ0", "tBAM_Other_Reads", "tBAM_Poor_Reads", "tBAM_REF_InDel_3bp", "tBAM_REF_InDel_2bp",
        "tBAM_REF_InDel_1bp", "tBAM_ALT_InDel_3bp", "tBAM_ALT_InDel_2bp", "tBAM_ALT_InDel_1bp",
    ], enumerate(header)))
    mq_features = map(lambda x: x[0], filter(lambda x: x[1] in [
                      "nBAM_REF_MQ", "nBAM_ALT_MQ", "tBAM_REF_MQ", "tBAM_ALT_MQ"], enumerate(header)))
    bq_features = map(lambda x: x[0], filter(lambda x: x[1] in [
                      "nBAM_REF_BQ", "nBAM_ALT_BQ", "tBAM_REF_BQ", "tBAM_ALT_BQ"], enumerate(header)))
    nm_diff_features = map(lambda x: x[0], filter(
        lambda x: x[1] in ["nBAM_NM_Diff", "tBAM_NM_Diff"], enumerate(header)))
    ranksum_features = map(lambda x: x[0], filter(lambda x: x[1] in ["nBAM_Z_Ranksums_MQ", "nBAM_Z_Ranksums_BQ",
                                                                     "nBAM_Z_Ranksums_EndPos", "tBAM_Z_Ranksums_BQ",  "tBAM_Z_Ranksums_MQ", "tBAM_Z_Ranksums_EndPos", ], enumerate(header)))
    zero_to_one_features = map(lambda x: x[0], filter(lambda x: x[1] in ["if_MuTect", "if_VarScan2", "if_SomaticSniper", "if_VarDict",
                                                                         "MuSE_Tier", "if_Strelka"] + ["nBAM_Concordance_FET", "nBAM_StrandBias_FET", "nBAM_Clipping_FET",
                                                                                                       "tBAM_Concordance_FET", "tBAM_StrandBias_FET", "tBAM_Clipping_FET"] + ["if_dbsnp", "COMMON"] + ["M2_STR"], enumerate(header)))
    stralka_scor = map(lambda x: x[0], filter(
        lambda x: x[1] in ["Strelka_Score"], enumerate(header)))
    stralka_qss = map(lambda x: x[0], filter(
        lambda x: x[1] in ["Strelka_QSS"], enumerate(header)))
    stralka_tqss = map(lambda x: x[0], filter(
        lambda x: x[1] in ["Strelka_TQSS"], enumerate(header)))
    varscan2_score = map(lambda x: x[0], filter(
        lambda x: x[1] in ["VarScan2_Score"], enumerate(header)))
    vardict_score = map(lambda x: x[0], filter(
        lambda x: x[1] in ["VarDict_Score"], enumerate(header)))
    m2_lod = map(lambda x: x[0], filter(lambda x: x[1] in [
                 "M2_NLOD", "M2_TLOD"], enumerate(header)))
    sniper_score = map(lambda x: x[0], filter(
        lambda x: x[1] in ["Sniper_Score"], enumerate(header)))
    m2_ecent = map(lambda x: x[0], filter(
        lambda x: x[1] in ["M2_ECNT"], enumerate(header)))
    sor = map(lambda x: x[0], filter(
        lambda x: x[1] in ["SOR"], enumerate(header)))
    msi = map(lambda x: x[0], filter(
        lambda x: x[1] in ["MSI"], enumerate(header)))
    msilen = map(lambda x: x[0], filter(
        lambda x: x[1] in ["MSILEN"], enumerate(header)))
    shift3 = map(lambda x: x[0], filter(
        lambda x: x[1] in ["SHIFT3"], enumerate(header)))
    MaxHomopolymer_Length = map(lambda x: x[0], filter(
        lambda x: x[1] in ["MaxHomopolymer_Length"], enumerate(header)))
    SiteHomopolymer_Length = map(lambda x: x[0], filter(
        lambda x: x[1] in ["SiteHomopolymer_Length"], enumerate(header)))
    InDel_Length = map(lambda x: x[0], filter(
        lambda x: x[1] in ["InDel_Length"], enumerate(header)))

    min_max_features = [[cov_features, 0, 2 * COV],
                        [mq_features, 0, 70],
                        [bq_features, 0, 41],
                        [nm_diff_features, -2 * COV, 2 * COV],
                        [zero_to_one_features, 0, 1],
                        [ranksum_features, -30, 30],
                        [stralka_scor, 0, 40],
                        [stralka_qss, 0, 200],
                        [stralka_tqss, 0, 4],
                        [varscan2_score, 0, 60],
                        [vardict_score, 0, 120],
                        [m2_lod, 0, 100],
                        [sniper_score, 0, 120],
                        [m2_ecent, 0, 40],
                        [sor, 0, 100],
                        [msi, 0, 100],
                        [msilen, 0, 10],
                        [shift3, 0, 100],
                        [MaxHomopolymer_Length, 0, 50],
                        [SiteHomopolymer_Length, 0, 50],
                        [InDel_Length, -30, 30],
                        ]
    selected_features = sorted([i for f in min_max_features for i in f[0]])
    selected_features_tags = map(lambda x: header[x], selected_features)
    for i_s, mn, mx in min_max_features:
        s = ensemble_data[:, np.array(i_s)]
        s = np.maximum(np.minimum(s, mx), mn)
        s = (s - mn) / (mx - mn)
        ensemble_data[:, np.array(i_s)] = s
    ensemble_data = ensemble_data[:, selected_features]
    ensemble_data = ensemble_data.tolist()
    ensemble_bed = os.path.join(work, "ensemble.bed")
    with open(ensemble_bed, "w")as f_:
        f_.write(
            "#" + "\t".join(map(str, header_pos + selected_features_tags)) + "\n")
        for i, s in enumerate(ensemble_data):
            f_.write("\t".join(map(str, ensemble_pos[i] + s)) + "\n")
    return ensemble_bed


def generate_dataset(work, truth_vcf_file, mode,  tumor_pred_vcf_file, region_bed_file, tumor_count_bed, normal_count_bed, ref_file,
                     matrix_width, matrix_base_pad, min_ev_frac_per_col, min_cov, num_threads, ensemble_tsv,
                     ensemble_bed, tsv_batch_size):
    logger = logging.getLogger(generate_dataset.__name__)

    logger.info("---------------------Generate Dataset----------------------")

    if not os.path.exists(work):
        os.mkdir(work)

    if mode == "train" and not truth_vcf_file:
        raise(RuntimeError("--truth_vcf is needed for 'train' mode"))

    if mode == "call":
        truth_vcf_file = os.path.join(work, "empty.vcf")
        with open(truth_vcf_file, "w") as o_f:
            o_f.write("##fileformat=VCFv4.2\n")
            o_f.write(
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")

    split_batch_size = 10000
    if ensemble_tsv and not ensemble_bed:
        ensemble_bed = extract_ensemble(work, ensemble_tsv)

    len_candids = len(pybedtools.BedTool(
        tumor_pred_vcf_file).intersect(region_bed_file, u=True))
    if ensemble_bed:
        len_candids += len(pybedtools.BedTool(ensemble_bed)
                           .intersect(region_bed_file, u=True))
    logger.info("len_candids: {}".format(len_candids))
    num_splits = max(len_candids / split_batch_size, num_threads)
    split_region_files = split_region(
        work, region_bed_file, num_splits, shuffle_intervals=True)

    pool = multiprocessing.Pool(num_threads)
    map_args = []
    for i, split_region_file in enumerate(split_region_files):
        map_args.append((work, split_region_file, truth_vcf_file,
                         tumor_pred_vcf_file, ref_file, ensemble_bed, i))
    try:
        records_data = pool.map_async(find_records, map_args).get()
        pool.close()
    except Exception as inst:
        logger.error(inst)
        pool.close()
        traceback.print_exc()
        raise Exception

    for o in records_data:
        if o is None:
            raise Exception("find_records failed!")

    none_vcf = "{}/none.vcf".format(work)
    var_vcf = "{}/var.vcf".format(work)
    if not os.path.exists(work):
        os.mkdir("{}".format(work))

    total_ims = 0
    for records_r, none_records, vtype, record_len, record_center, chroms_order, anns in records_data:
        total_ims += len(records_r) + len(none_records)

    candidates_split = int(total_ims / tsv_batch_size) + 1
    is_split = total_ims / candidates_split
    with open(var_vcf, "w") as vv, open(none_vcf, "w") as nv:
        is_current = 0
        is_ = -1
        while is_ < candidates_split:
            is_ += 1
            cnt = -1
            if is_ < candidates_split - 1:
                is_end = is_current + is_split
            else:
                is_end = total_ims
            candidates_tsv_file = "{}/candidates_{}.tsv".format(work, is_)
            logger.info("Write {}/{} split to {} for cnts ({}..{})/{}".format(
                is_ + 1, candidates_split, candidates_tsv_file, is_current,
                is_end, total_ims))
            pool = multiprocessing.Pool(num_threads)
            map_args_records = []
            map_args_nones = []
            for records_r, none_records, vtype, record_len, record_center, chroms_order, anns in records_data:
                if len(records_r) + cnt < is_current:
                    cnt += len(records_r)
                else:
                    for record in records_r:
                        cnt += 1
                        if is_current <= cnt < is_end:
                            vartype = vtype[int(record[-1])]
                            rlen = record_len[int(record[-1])]
                            rcenter = record_center[int(record[-1])]
                            ch_order = chroms_order[record[0]]
                            ann = anns[
                                int(record[-1])] if ensemble_bed else []
                            map_args_records.append((ref_file, tumor_count_bed, normal_count_bed, record, vartype, rlen, rcenter, ch_order,
                                                     matrix_base_pad, matrix_width, min_ev_frac_per_col, min_cov, ann))
                        if cnt >= is_end:
                            break
                    if cnt >= is_end:
                        break
                if cnt >= is_end:
                    break

                if len(none_records) + cnt < is_current:
                    cnt += len(none_records)
                else:
                    for record in none_records:
                        cnt += 1
                        if is_current <= cnt < is_end:
                            rcenter = record_center[int(record[-1])]
                            ch_order = chroms_order[record[0]]
                            ann = anns[
                                int(record[-1])] if ensemble_bed else []
                            map_args_nones.append((ref_file, tumor_count_bed, normal_count_bed, record, "NONE",
                                                   0, rcenter, ch_order,
                                                   matrix_base_pad, matrix_width, min_ev_frac_per_col, min_cov, ann))
                        if cnt >= is_end:
                            break
                    if cnt >= is_end:
                        break
                if cnt >= is_end:
                    break
            try:
                records_done = pool.map_async(
                    prep_data_single_tabix, map_args_records).get()
                pool.close()
            except Exception as inst:
                logger.error(inst)
                pool.close()
                traceback.print_exc()
                raise Exception

            for o in records_done:
                if o is None:
                    raise Exception("prep_data_single_tabix failed!")

            pool = multiprocessing.Pool(num_threads)
            try:
                none_records_done = pool.map_async(
                    prep_data_single_tabix, map_args_nones).get()
                pool.close()
            except Exception as inst:
                logger.error(inst)
                pool.close()
                traceback.print_exc()
                raise Exception

            for o in none_records_done:
                if o is None:
                    raise Exception("prep_data_single_tabix failed!")

            cnt_ims = 0
            tsv_idx = []
            with open(candidates_tsv_file, "w") as b_o:
                for x in records_done:
                    if x:
                        tag, compressed_candidate_mat, record, ann = x
                        vv.write("\t".join([record[0], str(record[1]), ".", record[2], record[
                                 3], ".", ".", "TAG={};".format(tag), ".", "."]) + "\n")
                        tsv_idx.append(b_o.tell())
                        b_o.write("\t".join([str(cnt_ims), "1", tag, compressed_candidate_mat] + map(
                            lambda x: str(np.round(x, 4)), ann)) + "\n")
                        cnt_ims += 1
                for x in none_records_done:
                    if x:
                        tag, compressed_candidate_mat, record, ann = x
                        nv.write("\t".join([record[0], str(record[1]), ".", record[2], record[
                                 3], ".", ".", "TAG={};".format(tag), ".", "."]) + "\n")
                        tsv_idx.append(b_o.tell())
                        b_o.write("\t".join([str(cnt_ims), "1", tag, compressed_candidate_mat] + map(
                            lambda x: str(np.round(x, 4)), ann)) + "\n")
                        cnt_ims += 1
                tsv_idx.append(b_o.tell())
            pickle.dump(tsv_idx, open(candidates_tsv_file + ".idx", "w"))
            is_current = is_end
            if is_current >= total_ims:
                break
    done_flag = "{}/done.txt".format(work)
    with open(done_flag, "w") as d_f:
        d_f.write("Done")


if __name__ == '__main__':
    FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    logging.basicConfig(level=logging.INFO, format=FORMAT)
    logger = logging.getLogger(__name__)

    parser = argparse.ArgumentParser(
        description='Generate dataset for train/call candidate variants on CNN')
    parser.add_argument('--mode', type=str, help='train/call mode',
                        choices=["train", "call"], required=True)
    parser.add_argument('--truth_vcf', type=str,
                        help='truth vcf (required for train mode)', default=None)
    parser.add_argument('--tumor_pred_vcf', type=str,
                        help='tumor candidate variants vcf file', required=True)
    parser.add_argument('--region_bed', type=str,
                        help='region bed', required=True)
    parser.add_argument('--tumor_count_bed', type=str,
                        help='tumor count bed.gz tabix file', required=True)
    parser.add_argument('--normal_count_bed', type=str,
                        help='normal count bed.gz tabix file', required=True)
    parser.add_argument('--reference', type=str,
                        help='reference fasta filename', required=True)
    parser.add_argument('--work', type=str,
                        help='work directory', required=True)
    parser.add_argument('--tsv_batch_size', type=int,
                        help='output files batch size', default=50000)
    parser.add_argument('--matrix_window_size', type=int,
                        help='target window width', default=32)
    parser.add_argument('--matrix_base_pad', type=int,
                        help='number of bases to pad around the candidate variant', default=7)
    parser.add_argument('--min_ev_frac_per_col', type=float,
                        help='minimum frac cov per column to keep columm', default=0.06)
    parser.add_argument('--min_cov', type=int, help='minimum cov', default=5)
    parser.add_argument('--num_threads', type=int,
                        help='number of threads', default=1)
    parser.add_argument('--ensemble_tsv', type=str,
                        help='Ensemble annotation tsv file (only for short read)', default=None)
    parser.add_argument('--ensemble_bed', type=str,
                        help='Ensemble annotation bed file (only for short read)', default=None)
    args = parser.parse_args()
    logger.info(args)

    mode = args.mode
    truth_vcf_file = args.truth_vcf
    tumor_pred_vcf_file = args.tumor_pred_vcf
    region_bed_file = args.region_bed
    tumor_count_bed = args.tumor_count_bed
    normal_count_bed = args.normal_count_bed
    ref_file = args.reference
    work = args.work
    matrix_width = args.matrix_window_size
    matrix_base_pad = args.matrix_base_pad
    min_ev_frac_per_col = args.min_ev_frac_per_col
    min_cov = args.min_cov
    num_threads = args.num_threads
    ensemble_tsv = args.ensemble_tsv
    ensemble_bed = args.ensemble_bed
    tsv_batch_size = args.tsv_batch_size

    try:
        generate_dataset(work, truth_vcf_file, mode, tumor_pred_vcf_file, region_bed_file, tumor_count_bed, normal_count_bed, ref_file,
                         matrix_width, matrix_base_pad, min_ev_frac_per_col, min_cov, num_threads, ensemble_tsv,
                         ensemble_bed, tsv_batch_size)
    except Exception as e:
        logger.error(traceback.format_exc())
        logger.error("Aborting!")
        logger.error(
            "generate_dataset.py failure on arguments: {}".format(args))
        raise e
