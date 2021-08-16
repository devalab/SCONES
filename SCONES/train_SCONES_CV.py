import os
import sys
import argparse
import random
import logging

from importlib import reload
from datetime import datetime

import torch
import numpy as np

from datasets.SCONES import SCONESDataset
from model.scones import SCONESNetwork

device = 'cpu' #'cuda:0' if torch.cuda.is_available() else 'cpu'

model_opt = {
    # Dataset
    "radius_cutoff" : 8,

    # Optimizer and Scheduler
    "lr" : 0.005,
    "weight_decay" : 0.004,
    "batch_size" : 32,
    "epochs" : 40,
    "scheduler" : {
        "type" : "ReduceLROnPlateau",
        "factor" : 0.1,
        "patience" : 5,
        "cooldown" : 3,
        "min_lr" : 1e-5,
        "threshold" : 1e-3
    },
    "best_model_threshold" : 0.001,

    # Loss
    "symmetry_loss_scale" : 1.0,

    # More options
    "embed_freeze_epochs" : 10
}

def train_model(opt, train_dataset, val_dataset):
    logging.shutdown()
    reload(logging)

    logging.basicConfig(
        format="%(asctime)s %(message)s",
        datefmt="[%H:%M:%S]",
        level=logging.INFO,
    )

    logger = logging.getLogger('scones')

    model = SCONESNetwork(opt, logger).to(device=device)
    model.init_train()
    model.set_criterion()
    model.set_optimizer_and_scheduler()
    model.set_data_loader(train_dataset, val_dataset)
    model.learn(True)

def cross_validation(opt, dataset, K):
    indices = np.arange(len(dataset))
    np.random.shuffle(indices)
    split_stride = len(indices)//K
    for k in range(K):
        val_start = k * split_stride
        val_end = val_start + split_stride
        if k == K - 1:
            val_end = len(indices)

        val_indices = np.full(len(indices), False)
        val_indices[val_start : val_end] = True

        train_indices = np.logical_not(val_indices)

        train_dataset = torch.utils.data.dataset.Subset(dataset, indices[train_indices])
        val_dataset = torch.utils.data.dataset.Subset(dataset, indices[val_indices])
        train_model(opt, train_dataset, val_dataset)

def main():
    parser = argparse.ArgumentParser(description='Train SCONES')
    parser.add_argument('--train_dataset_path', dest='train_dataset_path', type=str, help='path to dataset in CSV', required=True, default=None)
    parser.add_argument('--dir_prefix', dest='dir_prefix', type=str, help='prefix for training session directory', required=True, default="")

    args = parser.parse_args()
    if args.train_dataset_path is None:
        parser.print_help()
        sys.exit()

    train_dataset_path = args.train_dataset_path
    dir_prefix = args.dir_prefix

    train_dataset = SCONESDataset()
    train_dataset.deserialize(train_dataset_path)

    # random number to minimize risk of multiple training sessions using the same directory
    model_training_dir = os.path.join(dir_prefix, "cv5-" + datetime.now().strftime("%Y-%m-%d %H_%M_%S") + "-%d" % random.randint(0, 10000))
    opt = {
        "device" : device,
        "model_category" : "SCONES",
        "model_subdir" : model_training_dir,
        "config" : model_opt
    }

    cross_validation(opt, train_dataset, 5)

if __name__ == "__main__":
    main()
