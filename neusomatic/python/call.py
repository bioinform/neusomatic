#!/usr/bin/env python
#-------------------------------------------------------------------------
# call.py
# Call variants using model trained by NeuSomatic network
#-------------------------------------------------------------------------

import os
import traceback
import argparse
import logging
import shutil
import pickle
import multiprocessing
import copy

import pysam
import numpy as np
from imageio import imwrite, imread
import torch
from torch.autograd import Variable
import torch.nn as nn
import torch.nn.functional as F
from torchvision import transforms
import torchvision

from network import NeuSomaticNet
from dataloader import NeuSomaticDataset, matrix_transform
from utils import get_chromosomes_order, prob2phred
from merge_tsvs import merge_tsvs

import torch._utils
try:
    torch._utils._rebuild_tensor_v2
except AttributeError:
    def _rebuild_tensor_v2(storage, storage_offset, size, stride, requires_grad, backward_hooks):
        tensor = torch._utils._rebuild_tensor(
            storage, storage_offset, size, stride)
        tensor.requires_grad = requires_grad
        tensor._backward_hooks = backward_hooks
        return tensor
    torch._utils._rebuild_tensor_v2 = _rebuild_tensor_v2


def get_type(ref, alt):
    logger = logging.getLogger(get_type.__name__)
    len_diff = len(ref) - len(alt.split(",")[0])
    if len_diff > 0:
        return "DEL"
    elif len_diff < 0:
        return "INS"
    else:
        return "SNP"


def call_variants(net, vartype_classes, call_loader, out_dir, model_tag, use_cuda):
    logger = logging.getLogger(call_variants.__name__)
    net.eval()
    nclasses = len(vartype_classes)
    final_preds = {}
    none_preds = {}
    true_path = {}

    final_preds = {}
    none_preds = {}
    loader_ = call_loader

    iii = 0
    j = 0
    for data in loader_:
        (matrices, labels, var_pos_s, var_len_s,
         non_transformed_matrices), (paths) = data

        paths_ = copy.deepcopy(paths)
        del paths
        paths = paths_

        matrices = Variable(matrices)
        iii += 1
        j += len(paths[0])
        if use_cuda:
            matrices = matrices.cuda()

        outputs, _ = net(matrices)
        [outputs1, outputs2, outputs3] = outputs

        _, predicted = torch.max(outputs1.data.cpu(), 1)
        pos_pred = outputs2.data.cpu().numpy()
        _, len_pred = torch.max(outputs3.data.cpu(), 1)
        len_pred = len_pred.numpy()

        preds = {}
        for i, path_ in enumerate(paths[0]):
            path = path_.split("/")[-1]
            preds[i] = [vartype_classes[predicted[i]], pos_pred[i], len_pred[i]]
            if vartype_classes[predicted[i]] != "NONE":
                file_name = "{}/matrices_{}/{}.png".format(
                    out_dir, model_tag, path)
                if not os.path.exists(file_name):
                    imwrite(file_name, np.array(
                        non_transformed_matrices[i, :, :, 0:3]))
                true_path[path] = file_name
                final_preds[path] = [vartype_classes[predicted[i]], pos_pred[i], len_pred[i],
                                     list(map(lambda x: round(x, 4), F.softmax(
                                         outputs1[i, :], 0).data.cpu().numpy())),
                                     list(map(lambda x: round(x, 4), F.softmax(
                                         outputs3[i, :], 0).data.cpu().numpy())),
                                     list(map(lambda x: round(x, 4),
                                              outputs1.data.cpu()[i].numpy())),
                                     list(map(lambda x: round(x, 4),
                                              outputs3.data.cpu()[i].numpy()))]
            else:
                none_preds[path] = [vartype_classes[predicted[i]], pos_pred[i], len_pred[i],
                                    list(map(lambda x: round(x, 4), F.softmax(
                                        outputs1[i, :], 0).data.cpu().numpy())),
                                    list(map(lambda x: round(x, 4), F.softmax(
                                        outputs3[i, :], 0).data.cpu().numpy())),
                                    list(map(lambda x: round(x, 4),
                                             outputs1.data.cpu()[i].numpy())),
                                    list(map(lambda x: round(x, 4),
                                             outputs3.data.cpu()[i].numpy()))]
        if (iii % 10 == 0):
            logger.info("Called {} candidates in this batch.".format(j))
    logger.info("Called {} candidates in this batch.".format(j))
    return final_preds, none_preds, true_path


def pred_vcf_records_path(record):
    path, true_path_, pred_all, chroms, vartype_classes, ref_file = record
    thread_logger = logging.getLogger(
        "{} ({})".format(pred_vcf_records_path.__name__, multiprocessing.current_process().name))
    try:
        fasta_file = pysam.FastaFile(ref_file)
        ACGT = "ACGT"
        I = imread(true_path_) / 255.0
        vcf_record = []
        Ih, Iw, _ = I.shape
        zref_pos = np.where((np.argmax(I[:, :, 0], 0) == 0) & (
            sum(I[:, :, 0], 0) > 0))[0]
        nzref_pos = np.where(
            (np.argmax(I[:, :, 0], 0) > 0) & (sum(I[:, :, 0], 0) > 0))[0]
        # zref_pos_0 = np.where((I[0, :, 0] > 0) & (sum(I[:, :, 0], 0) > 0))[0]
        # nzref_pos_0 = np.where((I[0, :, 0] == 0) & (sum(I[:, :, 0], 0) > 0))[0]
        # assert(len(set(zref_pos_0)^set(zref_pos))==0)
        # assert(len(set(nzref_pos_0)^set(nzref_pos))==0)

        chrom, pos, ref, alt, _, center, _, _, _ = path.split(
            ".")
        center = int(center)
        pos = int(pos)

        ins_no_zref_pos = False
        infered_type_pred = "NONE"
        pred = list(pred_all)
        min_acceptable_probmax = 0.05
        center_dist_roundback = 0.8
        too_far_center = 3
        neigh_bases = 7
        while True:
            pred_probs = pred[3]
            if sum(pred_probs) < min_acceptable_probmax:
                break
            amx_prob = np.argmax(pred_probs)
            type_pred = vartype_classes[amx_prob]
            if type_pred == "NONE":
                break
            center_pred = min(max(0, pred[1][0]), Iw - 1)
            if abs(center_pred - center) < center_dist_roundback:
                center_ = center
            else:
                center_ = int(round(center_pred))
            if type_pred != "INS" and center_ in zref_pos:
                if center in nzref_pos:
                    center_ = center
                elif len(nzref_pos) > 0:
                    center_ = nzref_pos[
                        np.argmin(abs(nzref_pos - center_pred))]
                else:
                    break
            elif type_pred == "INS" and center_ in nzref_pos:
                if len(zref_pos) == 0:
                    if len(alt) < len(ref) + neigh_bases:
                        break
                    center_ = center
                    ins_no_zref_pos = True
                else:
                    if center in zref_pos:
                        center_ = center
                    else:
                        center_ = zref_pos[
                            np.argmin(abs(zref_pos - center_pred))]
            if abs(center_ - center) > too_far_center:
                pred[3][amx_prob] = 0
                # thread_logger.warning("Too far center: path:{}, pred:{}".format(path, pred))
            else:
                pred[0] = type_pred
                infered_type_pred = type_pred
                break
        if infered_type_pred == "NONE":
            return vcf_record
        type_pred = infered_type_pred
        len_pred = pred[2]
        vartype_candidate = get_type(ref, alt)
        col_2_pos = {}
        if vartype_candidate == "DEL":
            anchor = [pos + 1, center]
        elif vartype_candidate == "INS":
            anchor = [pos, center - 1]
        elif vartype_candidate == "SNP":
            anchor = [pos, center]
        cnt = 0
        for i in nzref_pos:
            col_2_pos[i] = cnt
            cnt += 1
        if anchor[1] not in col_2_pos:
            # print "NNN",path,pred
            return vcf_record

        b = (anchor[0] - col_2_pos[anchor[1]])
        for i in nzref_pos:
            col_2_pos[i] += b
        pos_2_col = {v: k for k, v in col_2_pos.items()}

        if abs(center_pred - center) < too_far_center:
            if type_pred == "SNP":
                pos_ = col_2_pos[center_]
                ref_ = ""
                alt_ = ""
                for i in range(len_pred):
                    nzp = nzref_pos[nzref_pos >= (center_ + i)]
                    if len(nzp) > 0:
                        center__ = nzp[np.argmin(abs(nzp - (center_ + i)))]
                        rb = np.argmax(I[1:, center__, 0])
                        ref_ += ACGT[rb]
                        II = I.copy()
                        II[rb + 1, center__, 1] = 0
                        alt_ += ACGT[np.argmax(II[1:, center__, 1])]
                        if sum(I[1:, center__, 1]) == 0:
                            break
                if not ref_:
                    # print "SSS",path,pred
                    return vcf_record
            elif type_pred == "INS":
                if ins_no_zref_pos:
                    pos_, ref_, alt_ = pos, ref.upper(), alt.upper()
                else:
                    pos_ = -1
                    i_ = center_ - 1
                    for i_ in range(center_ - 1, 0, -1):
                        if i_ in nzref_pos:
                            pos_ = col_2_pos[i_]
                            break
                    if pos_ == -1:
                        # print "PPP-1",path,pred
                        return vcf_record
                    if (sum(I[1:, i_, 1]) == 0):
                        # path,pred,i_,nzref_pos,col_2_pos,I[1:,i_,1],true_path[path]
                        return vcf_record
                    ref_ = ACGT[np.argmax(I[1:, i_, 0])]
                    alt_ = ref_
                    if len_pred == 3:
                        len_pred = max(len(alt) - len(ref), len_pred)
                    for i in range(i_ + 1, Iw):
                        if i in zref_pos:
                            alt_ += ACGT[np.argmax(I[1:, i, 1])]
                        else:
                            break
                        if (len(alt_) - len(ref_)) >= len_pred:
                            break
            elif type_pred == "DEL":
                pos_ = col_2_pos[center_] - 1
                if pos_ not in pos_2_col:
                    # print "DDDDD2",path,pred,pos_2_col,pos_
                    return vcf_record
                ref_ = ACGT[np.argmax(I[1:, pos_2_col[pos_], 0])]
                alt_ = ref_
                if len_pred == 3:
                    len_pred = max(len(ref) - len(alt), len_pred)
                for i in range(center_, Iw):
                    if i in nzref_pos:
                        ref_ += ACGT[np.argmax(I[1:, i, 0])]
                    if (len(ref_) - len(alt_)) >= len_pred:
                        break
                if (len(ref_) - len(alt_)) < len_pred:
                    pos_, ref_, alt_ = pos, ref.upper(), alt.upper()
            chrom_ = chroms[int(chrom)]
            if fasta_file.fetch(chrom_, pos_ - 1, pos_ +
                                len(ref_) - 1).upper() != ref_.upper():
                # print "AAAA"
                return vcf_record
            if ref_ == alt_:
                return vcf_record

            pred = list(pred_all)
            if type_pred == "SNP":
                prob = pred[3][3] * (pred[4][1])
            elif type_pred == "INS":
                if not ins_no_zref_pos:
                    prob = pred[3][1] * (1 - pred[4][0])
                else:
                    prob = pred[3][1]
            else:
                prob = pred[3][0] * (1 - pred[4][0])
            vcf_record = [path, [chrom_, pos_,
                                 ref_, alt_, prob, [path, pred]]]
        else:
            return vcf_record
        return vcf_record
    except Exception as ex:
        thread_logger.error(traceback.format_exc())
        thread_logger.error(ex)
        return None


def pred_vcf_records(ref_file, final_preds, true_path, chroms, vartype_classes, num_threads):
    logger = logging.getLogger(pred_vcf_records.__name__)
    logger.info(
        "Prepare VCF records for predicted somatic variants in this batch.")
    map_args = []
    for path in final_preds.keys():
        map_args.append([path, true_path[path], final_preds[path],
                         chroms, vartype_classes, ref_file])

    pool = multiprocessing.Pool(num_threads)
    try:
        all_vcf_records = pool.map_async(pred_vcf_records_path, map_args).get()
        pool.close()
    except Exception as inst:
        logger.error(inst)
        pool.close()
        traceback.print_exc()
        raise Exception

    for o in all_vcf_records:
        if o is None:
            raise Exception("pred_vcf_records_path failed!")

    all_vcf_records = list(filter(None, all_vcf_records))

    return all_vcf_records


def pred_vcf_records_none(none_preds, chroms):
    logger = logging.getLogger(pred_vcf_records_none.__name__)
    logger.info(
        "Prepare VCF records for predicted non-somatic variants in this batch.")
    all_vcf_records = {}
    for path in none_preds.keys():
        pred = none_preds[path]
        chrom, pos, ref, alt, _, _, _, _, _ = path.split(
            ".")
        if len(ref) == len(alt):
            prob = pred[3][3] * (pred[4][1])
        elif len(ref) < len(alt):
            prob = pred[3][1] * (1 - pred[4][0])
        else:
            prob = pred[3][0] * (1 - pred[4][0])
        if prob > 0.005:
            all_vcf_records[path] = [
                chroms[int(chrom)], pos, ref, alt, prob, [path, pred]]
    return all_vcf_records.items()


def get_vcf_records(all_vcf_records):
    logger = logging.getLogger(get_vcf_records.__name__)
    vcf_records = []
    for path in all_vcf_records:
        if all_vcf_records[path]:
            vcf_records.append(all_vcf_records[path][:-1])
    return vcf_records


def write_vcf(vcf_records, output_vcf, chroms_order, pass_threshold, lowqual_threshold):
    logger = logging.getLogger(write_vcf.__name__)
    vcf_records = list(filter(lambda x: len(x) > 0, vcf_records))
    vcf_records = sorted(vcf_records, key=lambda x: [chroms_order[x[0]], x[1]])
    lines = []
    with open(output_vcf, "w") as ov:
        for chrom_, pos_, ref_, alt_, prob in vcf_records:
            if ref_ == alt_:
                continue
            filter_ = "REJECT"
            if prob >= pass_threshold:
                filter_ = "PASS"
            elif prob >= lowqual_threshold:
                filter_ = "LowQual"
            line = "\t".join([chrom_, str(pos_), ".", ref_, alt_,
                              "{:.4f}".format(np.round(prob2phred(prob), 4)),
                              filter_, "SCORE={:.4f}".format(
                                  np.round(prob, 4)),
                              "GT", "0/1"]) + "\n"
            if line not in lines:
                ov.write(line)
                lines.append(line)


def call_neusomatic(candidates_tsv, ref_file, out_dir, checkpoint, num_threads,
                    batch_size, max_load_candidates, pass_threshold, lowqual_threshold,
                    use_cuda):
    logger = logging.getLogger(call_neusomatic.__name__)

    logger.info("-----------------Call Somatic Mutations--------------------")

    logger.info("PyTorch Version: {}".format(torch.__version__))
    logger.info("Torchvision Version: {}".format(torchvision.__version__))
    if not use_cuda:
        torch.set_num_threads(num_threads)

    chroms_order = get_chromosomes_order(reference=ref_file)
    with pysam.FastaFile(ref_file) as rf:
        chroms = rf.references

    vartype_classes = ['DEL', 'INS', 'NONE', 'SNP']
    data_transform = matrix_transform((0.5, 0.5, 0.5), (0.5, 0.5, 0.5))

    ensemble = False
    with open(candidates_tsv[0]) as i_f:
        for line in i_f:
            x = line.strip().split()
            if len(x) == 97:
                ensemble = True
            break

    num_channels = 119 if ensemble else 26
    logger.info("Number of channels: {}".format(num_channels))
    net = NeuSomaticNet(num_channels)
    if use_cuda:
        logger.info("GPU calling!")
        net.cuda()
    else:
        logger.info("CPU calling!")

    if torch.cuda.device_count() > 1:
        logger.info("We use {} GPUs!".format(torch.cuda.device_count()))
        net = nn.DataParallel(net)

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    logger.info("Load pretrained model from checkpoint {}".format(checkpoint))
    pretrained_dict = torch.load(
        checkpoint, map_location=lambda storage, loc: storage)
    pretrained_state_dict = pretrained_dict["state_dict"]
    model_tag = pretrained_dict["tag"]
    logger.info("tag: {}".format(model_tag))

    matrices_dir = "{}/matrices_{}".format(out_dir, model_tag)
    if os.path.exists(matrices_dir):
        logger.warning("Remove matrices directory: {}".format(matrices_dir))
        shutil.rmtree(matrices_dir)
    os.mkdir(matrices_dir)
    coverage_thr = pretrained_dict["coverage_thr"]
    if "normalize_channels" in pretrained_dict:
        normalize_channels = pretrained_dict["normalize_channels"]
    else:
        normalize_channels = False

    model_dict = net.state_dict()

    # 1. filter out unnecessary keys
    # pretrained_state_dict = {
    #     k: v for k, v in pretrained_state_dict.items() if k in model_dict}
    if "module." in list(pretrained_state_dict.keys())[0] and "module." not in list(model_dict.keys())[0]:
        pretrained_state_dict = {k.split("module.")[1]: v for k, v in pretrained_state_dict.items(
        ) if k.split("module.")[1] in model_dict}
    elif "module." not in list(pretrained_state_dict.keys())[0] and "module." in list(model_dict.keys())[0]:
        pretrained_state_dict = {
            ("module." + k): v for k, v in pretrained_state_dict.items()
            if ("module." + k) in model_dict}
    else:
        pretrained_state_dict = {k: v for k,
                                 v in pretrained_state_dict.items() if k in model_dict}

    # 2. overwrite entries in the existing state dict
    model_dict.update(pretrained_state_dict)
    # 3. load the new state dict
    net.load_state_dict(pretrained_state_dict)

    new_split_tsvs_dir = os.path.join(out_dir, "split_tsvs")
    if os.path.exists(new_split_tsvs_dir):
        logger.warning(
            "Remove split candidates directory: {}".format(new_split_tsvs_dir))
        shutil.rmtree(new_split_tsvs_dir)
    os.mkdir(new_split_tsvs_dir)
    Ls = []
    candidates_tsv_ = []
    split_i = 0
    for candidate_file in candidates_tsv:
        idx = pickle.load(open(candidate_file + ".idx", "rb"))
        if len(idx) > max_load_candidates / 2:
            logger.info("Splitting {} of lenght {}".format(
                candidate_file, len(idx)))
            new_split_tsvs_dir_i = os.path.join(
                new_split_tsvs_dir, "split_{}".format(split_i))
            if os.path.exists(new_split_tsvs_dir_i):
                logger.warning("Remove split candidates directory: {}".format(
                    new_split_tsvs_dir_i))
                shutil.rmtree(new_split_tsvs_dir_i)
            os.mkdir(new_split_tsvs_dir_i)
            candidate_file_splits = merge_tsvs(input_tsvs=[candidate_file],
                                               out=new_split_tsvs_dir_i,
                                               candidates_per_tsv=max(
                                                   1, max_load_candidates / 2),
                                               max_num_tsvs=100000,
                                               overwrite_merged_tsvs=True,
                                               keep_none_types=True)
            for candidate_file_split in candidate_file_splits:
                idx_split = pickle.load(
                    open(candidate_file_split + ".idx", "rb"))
                candidates_tsv_.append(candidate_file_split)
                Ls.append(len(idx_split) - 1)
            split_i += 1
        else:
            candidates_tsv_.append(candidate_file)
            Ls.append(len(idx) - 1)

    current_L = 0
    candidate_files = []
    all_vcf_records = []
    all_vcf_records_none = []
    for i, (candidate_file, L) in enumerate(sorted(zip(candidates_tsv_, Ls), key=lambda x: x[1])):
        current_L += L
        candidate_files.append(candidate_file)
        if current_L > max_load_candidates / 10 or i == len(candidates_tsv_) - 1:
            logger.info("Run for candidate files: {}".format(candidate_files))
            call_set = NeuSomaticDataset(roots=candidate_files,
                                         max_load_candidates=max_load_candidates,
                                         transform=data_transform, is_test=True,
                                         num_threads=num_threads,
                                         coverage_thr=coverage_thr,
                                         normalize_channels=normalize_channels)
            call_loader = torch.utils.data.DataLoader(call_set,
                                                      batch_size=batch_size,
                                                      shuffle=True, pin_memory=True,
                                                      num_workers=num_threads)

            current_L = 0
            candidate_files = []

            logger.info("N_dataset: {}".format(len(call_set)))
            if len(call_set) == 0:
                logger.warning(
                    "Skip {} with 0 candidates".format(candidate_file))
                continue

            final_preds_, none_preds_, true_path_ = call_variants(
                net, vartype_classes, call_loader, out_dir, model_tag, use_cuda)
            all_vcf_records.extend(pred_vcf_records(
                ref_file, final_preds_, true_path_, chroms, vartype_classes, num_threads))
            all_vcf_records_none.extend(
                pred_vcf_records_none(none_preds_, chroms))

    all_vcf_records = dict(all_vcf_records)
    all_vcf_records_none = dict(all_vcf_records_none)

    if os.path.exists(new_split_tsvs_dir):
        logger.warning(
            "Remove split candidates directory: {}".format(new_split_tsvs_dir))
        shutil.rmtree(new_split_tsvs_dir)

    logger.info("Prepare Output VCF")
    output_vcf = "{}/pred.vcf".format(out_dir)
    var_vcf_records = get_vcf_records(all_vcf_records)
    write_vcf(var_vcf_records, output_vcf, chroms_order,
              pass_threshold, lowqual_threshold)

    logger.info("Prepare Non-Somatics VCF")
    output_vcf_none = "{}/none.vcf".format(out_dir)
    vcf_records_none = get_vcf_records(all_vcf_records_none)
    write_vcf(vcf_records_none, output_vcf_none,
              chroms_order, pass_threshold, lowqual_threshold)

    if os.path.exists(matrices_dir):
        logger.warning("Remove matrices directory: {}".format(matrices_dir))
        shutil.rmtree(matrices_dir)

    logger.info("Calling is Done.")

    return output_vcf

if __name__ == '__main__':

    FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    logging.basicConfig(level=logging.INFO, format=FORMAT)
    logger = logging.getLogger(__name__)

    parser = argparse.ArgumentParser(
        description='simple call variants from bam')
    parser.add_argument('--candidates_tsv', nargs="*",
                        help=' test candidate tsv files', required=True)
    parser.add_argument('--reference', type=str,
                        help='reference fasta filename', required=True)
    parser.add_argument('--out', type=str,
                        help='output directory', required=True)
    parser.add_argument('--checkpoint', type=str,
                        help='network model checkpoint path', required=True)
    parser.add_argument('--num_threads', type=int,
                        help='number of threads', default=1)
    parser.add_argument('--batch_size', type=int,
                        help='batch size', default=1000)
    parser.add_argument('--max_load_candidates', type=int,
                        help='maximum candidates to load in memory', default=100000)
    parser.add_argument('--pass_threshold', type=float,
                        help='SCORE for PASS (PASS for score => pass_threshold)', default=0.7)
    parser.add_argument('--lowqual_threshold', type=float,
                        help='SCORE for LowQual (PASS for lowqual_threshold <= score < pass_threshold)',
                        default=0.4)
    args = parser.parse_args()

    use_cuda = torch.cuda.is_available()
    logger.info("use_cuda: {}".format(use_cuda))

    try:
        output_vcf = call_neusomatic(args.candidates_tsv, args.reference, args.out,
                                     args.checkpoint,
                                     args.num_threads, args.batch_size, args.max_load_candidates,
                                     args.pass_threshold, args.lowqual_threshold,
                                     use_cuda)
    except Exception as e:
        logger.error(traceback.format_exc())
        logger.error("Aborting!")
        logger.error("call.py failure on arguments: {}".format(args))
        raise e
