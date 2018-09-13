#-------------------------------------------------------------------------
# train.py
# Train NeuSomatic network
#-------------------------------------------------------------------------

import os
import traceback
import argparse
import datetime
import logging

import numpy as np
import torch
from torch.autograd import Variable
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torchvision import transforms

from network import NeuSomaticNet
from dataloader import NeuSomaticDataset

type_class_dict = {"DEL": 0, "INS": 1, "NONE": 2, "SNP": 3}
vartype_classes = ['DEL', 'INS', 'NONE', 'SNP']


def make_weights_for_balanced_classes(count_class_t, count_class_l, nclasses_t, nclasses_l,
                                      none_count=None):
    logger = logging.getLogger(make_weights_for_balanced_classes.__name__)

    w_t = [0] * nclasses_t
    w_l = [0] * nclasses_l

    count_class_t = list(count_class_t)
    count_class_l = list(count_class_l)
    if none_count:
        count_class_t[type_class_dict["NONE"]] = none_count
        count_class_l[0] = none_count

    logger.info("count type classes: {}".format(
        zip(vartype_classes, count_class_t)))
    N = float(sum(count_class_t))
    for i in range(nclasses_t):
        w_t[i] = (1 - (float(count_class_t[i]) / float(N))) / float(nclasses_t)
    w_t = np.array(w_t)
    logger.info("weight type classes: {}".format(zip(vartype_classes, w_t)))

    logger.info("weight length classes: {}".format(
        zip(range(nclasses_l), count_class_l)))
    N = float(sum(count_class_l))
    for i in range(nclasses_l):
        w_l[i] = (1 - (float(count_class_l[i]) / float(N))) / float(nclasses_l)
    w_l = np.array(w_l)
    logger.info("weight length classes: {}".format(
        zip(range(nclasses_l), w_l)))
    return w_t, w_l


def test(net, epoch, validation_loader, use_cuda):
    logger = logging.getLogger(test.__name__)
    net.eval()
    nclasses = len(vartype_classes)
    class_correct = list(0. for i in range(nclasses))
    class_total = list(0. for i in range(nclasses))
    class_p_total = list(0. for i in range(nclasses))

    len_class_correct = list(0. for i in range(4))
    len_class_total = list(0. for i in range(4))
    len_class_p_total = list(0. for i in range(4))

    falses = []
    for data in validation_loader:
        (matrices, labels, _, var_len_s, _), (paths) = data

        matrices = Variable(matrices)
        if use_cuda:
            matrices = matrices.cuda()

        outputs, _ = net(matrices)
        [outputs1, outputs2, outputs3] = outputs

        _, predicted = torch.max(outputs1.data.cpu(), 1)
        pos_pred = outputs2.data.cpu().numpy()
        _, len_pred = torch.max(outputs3.data.cpu(), 1)
        preds = {}
        for i, _ in enumerate(paths[0]):
            preds[i] = [vartype_classes[predicted[i]], pos_pred[i], len_pred[i]]

        compare_labels = (predicted == labels).squeeze()
        false_preds = np.where(compare_labels.numpy() == 0)[0]
        if len(false_preds) > 0:
            for i in false_preds:
                falses.append([paths[0][i], vartype_classes[predicted[i]], pos_pred[i], len_pred[i],
                               list(
                                   np.round(F.softmax(outputs1[i, :], 0).data.cpu().numpy(), 4)),
                               list(
                                   np.round(F.softmax(outputs3[i, :], 0).data.cpu().numpy(), 4))])

        for i in range(len(labels)):
            label = labels[i]
            class_correct[label] += compare_labels[i]
            class_total[label] += 1
        for i in range(len(predicted)):
            label = predicted[i]
            class_p_total[label] += 1

        compare_len = (len_pred == var_len_s).squeeze()
        for i in range(len(var_len_s)):
            len_ = var_len_s[i]
            len_class_correct[len_] += compare_len[i]
            len_class_total[len_] += 1
        for i in range(len(len_pred)):
            len_ = len_pred[i]
            len_class_p_total[len_] += 1

    for i in range(nclasses):
        SN = 100 * class_correct[i] / (class_total[i] + 0.0001)
        PR = 100 * class_correct[i] / (class_p_total[i] + 0.0001)
        F1 = 2 * PR * SN / (PR + SN + 0.0001)
        logger.info('Epoch {}: Type Accuracy of {:>5} ({}) : {:.2f}  {:.2f} {:.2f}'.format(
            epoch,
            vartype_classes[i], class_total[i],
            SN, PR, F1))
    logger.info('Epoch {}: Type Accuracy of the network on the {} test candidates: {:.4f} %'.format(
        epoch, sum(class_total), (
            100 * sum(class_correct) / float(sum(class_total)))))

    for i in range(4):
        SN = 100 * len_class_correct[i] / (len_class_total[i] + 0.0001)
        PR = 100 * len_class_correct[i] / (len_class_p_total[i] + 0.0001)
        F1 = 2 * PR * SN / (PR + SN + 0.0001)
        logger.info('Epoch {}: Type Accuracy of {:>5} ({}) : {:.2f}  {:.2f} {:.2f}'.format(
            epoch, i, len_class_total[i],
            SN, PR, F1))
    logger.info('Epoch {}: Type Accuracy of the network on the {} test candidates: {:.4f} %'.format(
        epoch, sum(len_class_total), (
            100 * sum(len_class_correct) / float(sum(len_class_total)))))

    net.train()


class SubsetNoneSampler(torch.utils.data.sampler.Sampler):

    def __init__(self, none_indices, var_indices, none_count):
        self.none_indices = none_indices
        self.var_indices = var_indices
        self.none_count = none_count
        self.current_none_id = 0

    def __iter__(self):
        if self.current_none_id > (len(self.none_indices) - self.none_count):
            this_round_nones = self.none_indices[self.current_none_id:]
            self.none_indices = map(lambda i: self.none_indices[i],
                                    torch.randperm(len(self.none_indices)).tolist())
            self.current_none_id = self.none_count - len(this_round_nones)
            this_round_nones += self.none_indices[0:self.current_none_id]
        else:
            this_round_nones = self.none_indices[
                self.current_none_id:self.current_none_id + self.none_count]
            self.current_none_id += self.none_count

        current_indices = this_round_nones + self.var_indices
        return iter(map(lambda i: current_indices[i], torch.randperm(len(current_indices))))

    def __len__(self):
        return len(self.var_indices) + self.none_count


def train_neusomatic(candidates_tsv, validation_candidates_tsv, out_dir, checkpoint,
                     num_threads, batch_size, max_epochs, learning_rate, lr_drop_epochs,
                     lr_drop_ratio, momentum, boost_none, none_count_scale,
                     max_load_candidates, coverage_thr, save_freq, use_cuda):
    logger = logging.getLogger(train_neusomatic.__name__)

    logger.info("-----------------------------------------------------------")
    logger.info("Train NeuSomatic Network")
    logger.info("-----------------------------------------------------------")

    data_transform = transforms.Compose(
        [transforms.ToTensor(),
         transforms.Normalize((0.5, 0.5, 0.5), (0.5, 0.5, 0.5))]
    )
    num_channels = 119 if args.ensemble else 26
    net = NeuSomaticNet(num_channels)
    if use_cuda:
        net.cuda()

    if torch.cuda.device_count() > 1:
        logger.info("We use {} GPUs!".format(torch.cuda.device_count()))
        net = nn.DataParallel(net)

    if not os.path.exists("{}/models/".format(out_dir)):
        os.mkdir("{}/models/".format(out_dir))

    if checkpoint:
        logger.info(
            "Load pretrained model from checkpoint {}".format(checkpoint))
        pretrained_dict = torch.load(
            checkpoint, map_location=lambda storage, loc: storage)
        pretrained_state_dict = pretrained_dict["state_dict"]
        tag = pretrained_dict["tag"]
        sofar_epochs = pretrained_dict["epoch"]
        logger.info(
            "sofar_epochs from pretrained checkpoint: {}".format(sofar_epochs))
        coverage_thr = pretrained_dict["coverage_thr"]
        logger.info(
            "Override coverage_thr from pretrained checkpoint: {}".format(coverage_thr))
        prev_epochs = sofar_epochs + 1
        model_dict = net.state_dict()
        # 1. filter out unnecessary keys
        # pretrained_state_dict = {
        # k: v for k, v in pretrained_state_dict.items() if k in model_dict}
        if "module." in pretrained_state_dict.keys()[0] and "module." not in model_dict.keys()[0]:
            pretrained_state_dict = {k.split("module.")[1]: v for k, v in pretrained_state_dict.items(
            ) if k.split("module.")[1] in model_dict}
        elif "module." not in pretrained_state_dict.keys()[0] and "module." in model_dict.keys()[0]:
            pretrained_state_dict = {
                ("module." + k): v for k, v in pretrained_state_dict.items() if ("module." + k) in model_dict}
        else:
            pretrained_state_dict = {k: v for k,
                                     v in pretrained_state_dict.items() if k in model_dict}
        # 2. overwrite entries in the existing state dict
        model_dict.update(pretrained_state_dict)
        # 3. load the new state dict
        net.load_state_dict(pretrained_state_dict)
    else:
        prev_epochs = 0
        time_now = datetime.datetime.now().strftime("%y-%m-%d-%H-%M-%S")
        tag = "neusomatic_{}".format(time_now)
    logger.info("tag: {}".format(tag))

    train_set = NeuSomaticDataset(roots=candidates_tsv,
                                  max_load_candidates=max_load_candidates,
                                  transform=data_transform, is_test=False,
                                  num_threads=num_threads, coverage_thr=coverage_thr)
    none_indices = train_set.get_none_indices()
    var_indices = train_set.get_var_indices()
    if none_indices:
        none_indices = map(lambda i: none_indices[i],
                           torch.randperm(len(none_indices)).tolist())
    logger.info("Non-somatic candidates: {}".format(len(none_indices)))
    if var_indices:
        var_indices = map(lambda i: var_indices[i],
                          torch.randperm(len(var_indices)).tolist())
    logger.info("Somatic candidates: {}".format(len(var_indices)))
    none_count = min(len(none_indices), len(var_indices) * none_count_scale)
    logger.info("Non-somatic considered in each epoch: {}".format(none_count))
    sampler = SubsetNoneSampler(none_indices, var_indices, none_count)

    train_loader = torch.utils.data.DataLoader(train_set,
                                               batch_size=batch_size,
                                               num_workers=num_threads, pin_memory=True,
                                               sampler=sampler)
    logger.info("#Train cadidates: {}".format(len(train_set)))

    if validation_candidates_tsv:
        validation_set = NeuSomaticDataset(roots=validation_candidates_tsv,
                                           max_load_candidates=max_load_candidates,
                                           transform=data_transform, is_test=True,
                                           num_threads=num_threads, coverage_thr=coverage_thr)
        validation_loader = torch.utils.data.DataLoader(validation_set,
                                                        batch_size=batch_size, shuffle=True,
                                                        num_workers=num_threads, pin_memory=True)
        logger.info("#Validation candidates: {}".format(len(validation_set)))

    weights_type, weights_length = make_weights_for_balanced_classes(
        train_set.count_class_t, train_set.count_class_l, 4, 4, none_count)

    weights_type[2] *= boost_none
    weights_length[0] *= boost_none

    logger.info("weights_type:{}, weights_length:{}".format(
        weights_type, weights_length))

    loss_s = []
    gradients = torch.FloatTensor(weights_type)
    gradients2 = torch.FloatTensor(weights_length)
    if use_cuda:
        gradients = gradients.cuda()
        gradients2 = gradients2.cuda()
    criterion_crossentropy = nn.CrossEntropyLoss(gradients)
    criterion_crossentropy2 = nn.CrossEntropyLoss(gradients2)
    criterion_smoothl1 = nn.SmoothL1Loss()
    optimizer = optim.SGD(
        net.parameters(), lr=learning_rate, momentum=momentum)

    net.train()
    len_train_set = none_count + len(var_indices)
    logger.info("Number of candidater per epoch: {}".format(len_train_set))
    print_freq = max(1, int(len_train_set / float(batch_size) / 4.0))
    curr_epoch = int(round(len(loss_s) / float(len_train_set)
                           * batch_size)) + prev_epochs
    torch.save({"state_dict": net.state_dict(),
                "tag": tag,
                "epoch": curr_epoch,
                "coverage_thr": coverage_thr},
               '{}/models/checkpoint_{}_epoch{}.pth'.format(out_dir, tag, curr_epoch))
    # loop over the dataset multiple times
    for epoch in range(max_epochs - prev_epochs):
        running_loss = 0.0
        for i, data in enumerate(train_loader, 0):
            # get the inputs
            (inputs, labels, var_pos_s, var_len_s, _), _ = data
            # wrap them in Variable
            inputs, labels, var_pos_s, var_len_s = Variable(inputs), Variable(
                labels), Variable(var_pos_s), Variable(var_len_s)
            if use_cuda:
                inputs, labels, var_pos_s, var_len_s = inputs.cuda(
                ), labels.cuda(), var_pos_s.cuda(), var_len_s.cuda()

            # zero the parameter gradients
            optimizer.zero_grad()

            outputs, _ = net(inputs)
            [outputs_classification, outputs_pos, outputs_len] = outputs
            var_len_labels = Variable(
                torch.LongTensor(var_len_s.cpu().data.numpy()))
            if use_cuda:
                var_len_labels = var_len_labels.cuda()
            loss = criterion_crossentropy(outputs_classification, labels) + 1 * criterion_smoothl1(
                outputs_pos, var_pos_s[:, 1]
            ) + 1 * criterion_crossentropy2(outputs_len, var_len_labels)

            loss.backward()
            optimizer.step()
            loss_s.append(loss.data[0])

            running_loss += loss.data[0]
            if i % print_freq == print_freq - 1:
                logger.info('epoch: {}, i: {:>5}, lr: {}, loss: {:.5f}'.format(
                            epoch + 1 + prev_epochs, i + 1,
                            learning_rate, running_loss / print_freq))
                running_loss = 0.0
        curr_epoch = int(
            round(len(loss_s) / float(len_train_set) * batch_size)) + prev_epochs
        if curr_epoch % save_freq == 0:
            torch.save({"state_dict": net.state_dict(),
                        "tag": tag,
                        "epoch": curr_epoch,
                        "coverage_thr": coverage_thr,
                        }, '{}/models/checkpoint_{}_epoch{}.pth'.format(out_dir, tag, curr_epoch))
            if validation_candidates_tsv:
                test(net, curr_epoch, validation_loader, use_cuda)
        if curr_epoch % lr_drop_epochs == 0:
            learning_rate *= lr_drop_ratio
            optimizer = optim.SGD(
                net.parameters(), lr=learning_rate, momentum=momentum)
    logger.info('Finished Training')

    curr_epoch = int(round(len(loss_s) / float(len_train_set)
                           * batch_size)) + prev_epochs
    torch.save({"state_dict": net.state_dict(),
                "tag": tag,
                "epoch": curr_epoch,
                "coverage_thr": coverage_thr}, '{}/models/checkpoint_{}_epoch{}.pth'.format(
        out_dir, tag, curr_epoch))
    if validation_candidates_tsv:
        test(net, curr_epoch, validation_loader, use_cuda)
    logger.info("Total Epochs: {}".format(curr_epoch))
    return '{}/models/checkpoint_{}_epoch{}.pth'.format(out_dir, tag, curr_epoch)

if __name__ == '__main__':

    FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    logging.basicConfig(level=logging.INFO, format=FORMAT)
    logger = logging.getLogger(__name__)

    parser = argparse.ArgumentParser(
        description='simple call variants from bam')
    parser.add_argument('--candidates_tsv', nargs="*",
                        help=' train candidate tsv files', required=True)
    parser.add_argument('--out', type=str,
                        help='output directory', required=True)
    parser.add_argument('--checkpoint', type=str,
                        help='pretrianed network model checkpoint path', default=None)
    parser.add_argument('--validation_candidates_tsv', nargs="*",
                        help=' validation candidate tsv files', default=[])
    parser.add_argument('--ensemble',
                        help='Enable training for ensemble mode',
                        action="store_true")
    parser.add_argument('--num_threads', type=int,
                        help='number of threads', default=1)
    parser.add_argument('--batch_size', type=int,
                        help='batch size', default=1000)
    parser.add_argument('--max_epochs', type=int,
                        help='maximum number of training epochs', default=1000)
    parser.add_argument('--lr', type=float, help='learning_rate', default=0.01)
    parser.add_argument('--lr_drop_epochs', type=int,
                        help='number of epochs to drop learning rate', default=400)
    parser.add_argument('--lr_drop_ratio', type=float,
                        help='learning rate drop scale', default=0.1)
    parser.add_argument('--momentum', type=float,
                        help='SGD momentum', default=0.9)
    parser.add_argument('--boost_none', type=float,
                        help='the amount to boost none-somatic classification weight', default=10)
    parser.add_argument('--none_count_scale', type=float,
                        help='ratio of the none/somatic canidates to use in each training epoch \
                        (the none candidate will be subsampled in each epoch based on this ratio',
                        default=2)
    parser.add_argument('--max_load_candidates', type=int,
                        help='maximum candidates to load in memory', default=1000000)
    parser.add_argument('--save_freq', type=int,
                        help='the frequency of saving checkpoints in terms of # epochs', default=50)
    parser.add_argument('--coverage_thr', type=int,
                        help='maximum coverage threshold to be used for network input \
                              normalization. \
                              Will be overridden if pretrained model is provided\
                              For ~50x WGS, coverage_thr=100 should work. \
                              For higher coverage WES, coverage_thr=300 should work.', default=100)
    args = parser.parse_args()

    logger.info(args)

    use_cuda = torch.cuda.is_available()
    logger.info("use_cuda: {}".format(use_cuda))

    try:
        checkpoint = train_neusomatic(args.candidates_tsv, args.validation_candidates_tsv,
                                      args.out, args.checkpoint, args.num_threads, args.batch_size,
                                      args.max_epochs,
                                      args.lr, args.lr_drop_epochs, args.lr_drop_ratio, args.momentum,
                                      args.boost_none, args.none_count_scale,
                                      args.max_load_candidates, args.coverage_thr, args.save_freq, use_cuda)
    except Exception as e:
        logger.error(traceback.format_exc())
        logger.error("Aborting!")
        logger.error(
            "traina.py failure on arguments: {}".format(args))
        raise e
