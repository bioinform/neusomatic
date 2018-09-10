#-------------------------------------------------------------------------
# network.py
# Defines the architecture of NeuSomatic network
#-------------------------------------------------------------------------
import logging

import torch.nn as nn
import torch.nn.functional as F
import numpy as np


FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logFormatter = logging.Formatter(FORMAT)
logger = logging.getLogger(__name__)
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
logger.addHandler(consoleHandler)
logging.getLogger().setLevel(logging.INFO)


class DSBlock(nn.Module):

    def __init__(self, dim, ks_1=3, ks_2=3, dl_1=1, dl_2=1, mp_ks=3, mp_st=1):
        super(DSBlock, self).__init__()
        self.dim = dim
        self.conv_r1 = nn.Conv2d(
            dim, dim, kernel_size=ks_1, dilation=dl_1, padding=(dl_1 * (ks_1 - 1)) / 2)
        self.bn_r1 = nn.BatchNorm2d(dim)
        self.conv_r2 = nn.Conv2d(
            dim, dim, kernel_size=ks_2, dilation=dl_2, padding=(dl_2 * (ks_2 - 1)) / 2)
        self.bn_r2 = nn.BatchNorm2d(dim)
        self.pool_r2 = nn.MaxPool2d((1, mp_ks), padding=(
            0, (mp_ks - 1) / 2), stride=(1, mp_st))

    def forward(self, x):
        y1 = (F.relu(self.bn_r1(self.conv_r1(x))))
        y2 = (self.bn_r2(self.conv_r2(y1)))
        y3 = x + y2
        z = self.pool_r2(y3)
        return z


class NeuSomaticNet(nn.Module):

    def __init__(self, num_channels):
        super(NeuSomaticNet, self).__init__()
        dim = 64
        self.conv1 = nn.Conv2d(num_channels, dim, kernel_size=(
            1, 3), padding=(0, 1), stride=1)
        self.bn1 = nn.BatchNorm2d(dim)
        self.pool1 = nn.MaxPool2d((1, 3), padding=(0, 1), stride=(1, 1))
        self.resblocks = [
            [3, 5, 1, 1, 3, 1],
            [3, 5, 1, 1, 3, 2],
            [3, 5, 2, 1, 3, 2],
            [3, 5, 4, 2, 3, 2],
        ]
        res_layers = []
        for ks_1, ks_2, dl_1, dl_2, mp_ks, mp_st in self.resblocks:
            rb = DSBlock(dim, ks_1, ks_2, dl_1, dl_2, mp_ks, mp_st)
            res_layers.append(rb)
        self.res_layers = nn.Sequential(*res_layers)
        ds = np.prod(map(lambda x: x[5], self.resblocks))
        self.fc_dim = dim * 32 * 5 / ds
        self.fc1 = nn.Linear(self.fc_dim, 240)
        self.fc2 = nn.Linear(240, 4)
        self.fc3 = nn.Linear(240, 1)
        self.fc4 = nn.Linear(240, 4)

    def forward(self, x):
        x = self.pool1(F.relu(self.bn1(self.conv1(x))))
        internal_outs = [x]

        x = self.res_layers(x)
        internal_outs.append(x)
        x2 = x.view(-1, self.fc_dim)
        x3 = F.relu(self.fc1(x2))
        internal_outs.extend([x2, x3])
        o1 = self.fc2(x3)
        o2 = self.fc3(x3)
        o3 = self.fc4(x3)
        return [o1, o2, o3], internal_outs
