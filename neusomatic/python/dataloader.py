#-------------------------------------------------------------------------
# dataloader.py
# Data loader used by NeuSomatic network for datasets created by 'generate_dataset.py'
#-------------------------------------------------------------------------
import multiprocessing
import pickle
import zlib
import glob
import logging
import base64
import traceback

from numpy import random
import numpy as np
import torch
import resource

FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger(__name__)

type_class_dict = {"DEL": 0, "INS": 1, "NONE": 2, "SNP": 3}


def extract_zlib(zlib_compressed_im):
    return np.fromstring(zlib.decompress(zlib_compressed_im), dtype="uint8").reshape((5, 32, 23))


def candidate_loader_tsv(tsv, open_tsv, idx, i):
    if open_tsv:
        i_f = open_tsv
    else:
        i_f = open(tsv, "r")
    i_f.seek(idx[i])
    fields = i_f.read(idx[i + 1] - idx[i]).strip().split()
    tag = fields[2]
    im = extract_zlib(base64.b64decode(fields[3]))
    if len(fields) > 4:
        anns = map(float, fields[4:])
    else:
        anns = []
    label = type_class_dict[tag.split(".")[4]]
    if not open_tsv:
        i_f.close()
    return tag, im, anns, label


def extract_info_tsv((i_b, tsv, idx, L, max_load_candidates, nclasses_t, nclasses_l)):
    thread_logger = logging.getLogger(
        "{} ({})".format(extract_info_tsv.__name__, multiprocessing.current_process().name))
    try:
        n_none = 0
        with open(tsv, "r") as i_f:
            for line in i_f:
                tag = line.strip().split()[2]
                n_none += (1 if "NONE" in tag else 0)
        n_var = L - n_none

        max_load_candidates_var = min(n_var, max_load_candidates)
        max_load_candidates_none = max_load_candidates - max_load_candidates_var

        j = 0
        matrices = []
        none_ids = []
        var_ids = []
        count_class_t = [0] * nclasses_t
        count_class_l = [0] * nclasses_l
        data = []
        cnt_none = 0
        cnt_var = 0
        with open(tsv, "r") as i_f:
            i = 0
            for line in i_f:
                fields = line.strip().split()
                ii = int(fields[0])
                assert ii == i
                tag = fields[2]
                matrices.append([i_b, i])
                if "NONE" in tag:
                    none_ids.append(j)
                else:
                    var_ids.append(j)
                j += 1
                _, _, _, _, vartype, _, length, _, _ = tag.split(".")
                count_class_t[type_class_dict[vartype]] += 1
                count_class_l[min(int(length), 3)] += 1
                if ((cnt_var < max_load_candidates_var) and ("NONE" not in tag)) or (
                        (cnt_none < max_load_candidates_none) and ("NONE" in tag)):
                    im = extract_zlib(base64.b64decode(fields[3]))
                    label = type_class_dict[tag.split(".")[4]]
                    if len(fields) > 4:
                        anns = map(float, fields[4:])
                    else:
                        anns = []
                    data.append([tag, im, anns, label])
                    if "NONE" in tag:
                        cnt_none += 1
                    else:
                        cnt_var += 1
                else:
                    data.append([])
                i += 1
            assert i == L
        thread_logger.info("Loaded {} candidates for {}".format(
            len(matrices), tsv))
        return matrices, data, none_ids, var_ids, count_class_t, count_class_l
    except Exception as ex:
        thread_logger.error(traceback.format_exc())
        thread_logger.error(ex)
        return None


class NeuSomaticDataset(torch.utils.data.Dataset):

    def __init__(self, roots, max_load_candidates, transform=None,
                 loader=candidate_loader_tsv, is_test=False,
                 num_threads=1, disable_ensemble=False, data_augmentation=False,
                 nclasses_t=4, nclasses_l=4, coverage_thr=100,
                 max_opended_tsv=-1):

        soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
        logger.info(resource.getrlimit(resource.RLIMIT_NOFILE))
        new_soft = max(soft, hard / 2)
        resource.setrlimit(resource.RLIMIT_NOFILE, (hard, hard))
        logger.info(resource.getrlimit(resource.RLIMIT_NOFILE))
        if max_opended_tsv == -1:
            max_opended_tsv = new_soft
        else:
            max_opended_tsv = min(max_opended_tsv, new_soft)
        self.da_shift_p = 0.3
        self.da_base_p = 0.05
        self.da_rev_p = 0.1
        self.da_cov_p = 0.2
        self.da_cov_e = 0.2

        self.matrices = []
        self.count_class_t = [1] * nclasses_t
        self.count_class_l = [1] * nclasses_l
        self.none_ids = []
        self.var_ids = []

        self.tsvs = []
        self.open_tsvs = [[] for t in range(num_threads)]
        opened_tsvs = 0
        self.num_threads = num_threads
        self.Ls = []
        self.idxs = []
        for root in roots:
            for tsv in glob.glob(root):
                self.tsvs.append(tsv)
                for t in range(num_threads):
                    if opened_tsvs < max_opended_tsv - 1:
                        self.open_tsvs[t].append(open(tsv))
                        opened_tsvs += 1
                    else:
                        self.open_tsvs[t].append("")
                idx = pickle.load(open(tsv + ".idx"))
                self.idxs.append(idx)
                self.Ls.append(len(idx) - 1)
        logger.info("Opened {}/{} tsv's".format(opened_tsvs,
                                                num_threads * len(self.tsvs)))
        self.data = []
        total_L = sum(self.Ls)
        batches = []
        new_batch = []
        for i_b, L in enumerate(self.Ls):
            new_batch.append([i_b, L])
            if sum(map(lambda x: x[1], new_batch)) > 200000 or i_b == len(self.Ls) - 1 \
                    or len(new_batch) > num_threads:
                batches.append(new_batch)
                new_batch = []

        records_done = [[] for i in range(len(batches))]
        for i, batch in enumerate(batches):
            map_args = []
            Ls_ = []
            for i_b, _ in batch:
                tsv = self.tsvs[i_b]
                max_load_ = self.Ls[i_b] * max_load_candidates / \
                    total_L if total_L > 0 else 0
                map_args.append([i_b, tsv, self.idxs[i_b], self.Ls[i_b],
                                 max_load_, nclasses_t, nclasses_l])
                Ls_.append(self.Ls[i_b])
            logger.info("Len's of tsv files in this batch: {}".format(Ls_))
            if len(map_args) == 1:
                records_ = [extract_info_tsv(map_args[0])]
            else:
                pool = multiprocessing.Pool(num_threads)
                try:
                    records_ = pool.map_async(
                        extract_info_tsv, map_args).get()
                    pool.close()
                except Exception as inst:
                    pool.close()
                    logger.error(inst)
                    traceback.print_exc()
                    raise Exception

            for o in records_:
                if o is None:
                    raise Exception("extract_info_tsv failed!")

            records_done[i] = records_

        j = 0
        for records_ in records_done:
            for matrices, data, none_ids, var_ids, count_class_t, count_class_l in records_:
                self.matrices += matrices
                self.data += data
                self.none_ids += map(lambda x: x + j, none_ids)
                self.var_ids += map(lambda x: x + j, var_ids)
                for k in range(nclasses_t):
                    self.count_class_t[k] += count_class_t[k]
                for k in range(nclasses_l):
                    self.count_class_l[k] += count_class_l[k]
                j += len(matrices)

        self.roots = roots
        self.transform = transform
        self.loader = loader
        self.is_test = is_test
        self.disable_ensemble = disable_ensemble
        self.data_augmentation = data_augmentation
        self.coverage_thr = coverage_thr

    def __getitem__(self, index):
        if len(self.data[index]) == 0:
            i_b, i = self.matrices[index]
            if multiprocessing.current_process()._identity:
                path, matrix, anns, label = candidate_loader_tsv(self.tsvs[i_b],
                                                             self.open_tsvs[
                int(multiprocessing.current_process()._identity[0]
                    ) % self.num_threads][i_b],
                self.idxs[i_b], i)
            else:
                path, matrix, anns, label = candidate_loader_tsv(self.tsvs[i_b],
                                                             self.open_tsvs[0][i_b],
                self.idxs[i_b], i)

        else:
            path, matrix, anns, label = self.data[index]

        if self.disable_ensemble:
            anns = []

        tag = path.split("/")[-1]
        _, _, _, _, vartype, center, length, tumor_cov, normal_cov = tag.split(
            ".")
        tumor_cov = int(tumor_cov)
        normal_cov = int(normal_cov)
        center = int(center)
        length = int(length)

        h, w, _ = matrix.shape
        orig_matrix = matrix.copy()
        orig_center = int(center)

        h, w, _ = matrix.shape
        far_center = False
        if (((center - 2) * 2 / 3) >= (center - 2)) or (((w - center - 2) * 2 / 3)
                                                        >= (w - center - 2)):
            far_center = True

        # Data augmentaion by shifting left or right
        if self.data_augmentation and (not self.is_test) and (random.rand() < self.da_shift_p
                                                              and (not far_center)):
            h, w, c = matrix.shape
            r = random.rand()
            if r < 0.6:
                x_left = random.randint((center - 2) * 2 / 3, center - 2)
            else:
                x_left = 0
            if r > 0.4:
                x_right = random.randint(
                    (w - center - 2) * 2 / 3, w - center - 2)
            else:
                x_right = 0
            if x_left > 0:
                matrix[:, 0:w - x_left, :] = matrix[:, x_left:, :]
                matrix[:, -x_left:, :] = -1
                center -= x_left
            if x_right > 0:
                matrix[:, x_right:, :] = matrix[:, 0:w - x_right, :]
                matrix[:, 0:x_right, :] = -1
                center += x_right

        # Data augmentaion by switch bases
        if self.data_augmentation and (not self.is_test) and random.rand() < self.da_base_p \
                and (vartype != "NONE"):
            [i, j] = random.permutation(range(1, 5))[0:2]
            a = matrix[i, :, :]
            matrix[i, :, :] = matrix[j, :, :]
            matrix[j, :, :] = a

        # Data augmentaion by random flip
        try:
            nt_matrix = matrix.copy()
            nt_center = int(center)
            if self.data_augmentation and (not self.is_test) and random.rand() < self.da_rev_p \
                    and (vartype not in ["DEL"]):
                h, w, c = matrix.shape
                refbase = np.nonzero(matrix[:, center, 0])[0]
                if len(refbase) > 1:
                    logger.warning("{} {}".format(path, refbase))
                if len(refbase) == 1:
                    refbase = refbase[0]
                    b = center
                    e = center + 1
                    if refbase != 0:
                        for i in range(center - 1, 0, -1):
                            if matrix[refbase, i, 0] == 1:
                                b -= 1
                            else:
                                break
                        for i in range(center + 1, w):
                            if matrix[refbase, i, 0] == 1:
                                e += 1
                            else:
                                break
                    elif refbase == 0:
                        hp_base = 0
                        for i in range(center + 1, w):
                            if sum(matrix[:, i, 0]) == 0:
                                break
                            base = np.nonzero(matrix[:, i, 0])[0][0]
                            if base != 0:
                                hp_base = base
                                break
                        if matrix[hp_base, center, 1] >= np.max(matrix[1:, center, 1]):
                            for i in range(center + 1, w):
                                if sum(matrix[:, i, 0]) == 0:
                                    break
                                base = np.nonzero(matrix[:, i, 0])[0][0]
                                if base != 0:
                                    if not matrix[hp_base, i, 0] == 1:
                                        break
                                    h += 1
                                else:
                                    mx_1 = np.max(matrix[:, i, 1])
                                    if matrix[0, i, 1] < mx_1 and matrix[hp_base, i, 1] < mx_1:
                                        e = center + 1
                                        break
                                e += 1
                        if h == 1:
                            e = center + 1
                    if (e - b) > 1:
                        matrix[:, b:e, :] = matrix[:, e - 1:b - 1:-1, :].copy()
                        center = e - 1 - (center - b)
                matrix = matrix[:, ::-1, :].copy()
                center = w - center - 1
        except:
            logger.warning(
                "Failed random flip center={} tag={}".format(center, tag))
            matrix = nt_matrix
            center = nt_center

        # Data augmentaion by changing coverage
        if self.data_augmentation and (not self.is_test) and random.rand() < self.da_cov_p:
            r_cov = (1 - self.da_cov_e) + (random.rand() * 2 * self.da_cov_e)
            tumor_cov *= r_cov
            normal_cov *= r_cov

        # add COV channel
        matrix_ = np.zeros((matrix.shape[0], matrix.shape[1], 26 + len(anns)))
        matrix_[:, :, 0:23] = matrix
        matrix = matrix_
        matrix[:, center, 23] = np.max(matrix[:, :, 0])
        matrix[:, :, 24] = (min(tumor_cov, self.coverage_thr) /
                            float(self.coverage_thr)) * 255.0
        matrix[:, :, 25] = (
            min(normal_cov, self.coverage_thr) / float(self.coverage_thr)) * 255.0
        for i, a in enumerate(anns):
            matrix[:, :, 26 + i] = a * 255.0

        if self.is_test:
            orig_matrix_ = np.zeros(
                (orig_matrix.shape[0], orig_matrix.shape[1], 3))
            orig_matrix_[:, :, 0:2] = orig_matrix[:, :, 0:2]
            orig_matrix_[:, orig_center, 2] = np.max(orig_matrix[:, :, 0])
            orig_matrix = orig_matrix_
            non_transformed_matrix = np.array(orig_matrix).astype(np.uint8)
        else:
            non_transformed_matrix = []

        if self.transform is not None:
            matrix = self.transform(matrix)

        var_pos = [length, center]
        var_pos = torch.Tensor(var_pos)
        varlen_label = min(length, 3)

        return (matrix, label, var_pos, varlen_label, non_transformed_matrix), [path, label]

    def __len__(self):
        return len(self.matrices)

    def get_none_indices(self):
        return self.none_ids

    def get_var_indices(self):
        return self.var_ids
