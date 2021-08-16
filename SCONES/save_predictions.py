import os
import sys
import argparse
import logging

import pandas as pd
import numpy as np

from datasets.SCONES import SCONESDataset
from model.dataset import PreprocessedDataset
from model.scones import SCONESNetwork

def get_ensemble_predictions(opt, model_paths, dataset_path, results_path):
    dataset = SCONESDataset()
    dataset.deserialize(dataset_path)

    pp_dataset = PreprocessedDataset(radius_cutoff=opt["config"]["radius_cutoff"], device=opt["device"])
    pp_dataset.add_dataset(dataset, training_filters=True)
    pp_dataset.summary()

    predictions = {}
    for sample in pp_dataset:
        predictions[sample["id"]] = []

    for path in model_paths:
        print(path)
        logging.basicConfig(
            format="%(asctime)s %(message)s",
            datefmt="[%H:%M:%S]",
            level=logging.INFO,
        )

        logger = logging.getLogger('scones')

        model = SCONESNetwork(opt, logger)
        model.init_eval()
        model.set_criterion()
        model.load_model(os.path.join(path, "best_val_model.pth"))
        model.eval()

        for sample in pp_dataset:
            sample_id = sample["id"]
            if sample_id not in predictions:
                continue
            _, _, prediction = model.predict(sample)
            predictions[sample_id].append(prediction.item())

    predictions = [{  "id" : key, "ddG_SCONES" : round(np.mean(preds), 2), "ddG_SCONES_stddev" : round(np.std(preds), 2) } for key, preds in predictions.items()]
    df = pd.DataFrame(predictions)
    df.to_csv(results_path, index=False)

def main():
    parser = argparse.ArgumentParser(description='Save predictions from SCONES training sessions')
    parser.add_argument('--dataset_path', dest='dataset_path', type=str, help='path to dataset in CSV', required=True, default=None)
    parser.add_argument('--results_path', dest='results_path', type=str, help='path to results CSV', required=True, default=None)
    parser.add_argument('--multicv_model_dir', dest='multicv_model_dir', type=str, help='path to cross-validation training directories', required=False, default=None)
    parser.add_argument('--cv_model_dir', dest='cv_model_dir', type=str, help='path to cross-validation training directory', required=False, default=None)
    parser.add_argument('--model_dir', dest='model_dir', type=str, help='path to a single training directory', required=False, default=None)

    args = parser.parse_args()
    if args.dataset_path is None:
        parser.print_help()
        sys.exit()

    if args.multicv_model_dir is None and args.cv_model_dir is None and args.model_dir is None:
        parser.print_help()
        print("Please provide (multiple) cross-validation training directory or a single model training directory")
        sys.exit()

    dataset_path = args.dataset_path
    results_path = args.results_path
    multicv_model_dir = args.multicv_model_dir
    cv_model_dir = args.cv_model_dir
    model_dir = args.model_dir

    model_paths = []
    if multicv_model_dir:
        directory_list = [f.path for f in os.scandir(multicv_model_dir) if f.is_dir()]
        for cv_dir in directory_list:
            model_paths += [f.path for f in os.scandir(cv_dir) if f.is_dir()]
    elif cv_model_dir:
        model_paths = [f.path for f in os.scandir(cv_model_dir) if f.is_dir()]
    else:
        model_paths = [model_dir]

    opt = {
        "model_category" : "SCONES",
        "device" : 'cpu',
        "config" : {
            "radius_cutoff" : 8,
        }
    }

    get_ensemble_predictions(opt, model_paths, dataset_path, results_path)

if __name__ == "__main__":
    main()
