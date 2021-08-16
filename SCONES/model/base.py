import os
import logging

import torch
import torch.nn as nn
import numpy as np

class BaseNetwork(nn.Module):
    def __init__(self, opt, logger):
        super(BaseNetwork, self).__init__()
        self.opt = opt
        self.device = torch.device(self.opt["device"])
        self.logger = logger

    def init_train(self):
        self.train()
        
        model_category = self.opt["model_category"]
        model_subdir = self.opt.get("model_subdir", "")

        from datetime import datetime
        timestamp = datetime.now().strftime("%Y-%m-%d %H_%M_%S")
        self.output_dir = os.path.join("output", model_category, model_subdir, timestamp)
        os.makedirs(self.output_dir)
        
        file_handler = logging.FileHandler(filename=os.path.join(self.output_dir, "log.txt"))
        self.logger.addHandler(file_handler)

        import json
        opt_prettyprint_str = json.dumps(self.opt, indent=4)
        self.logger.info("Model: %s", model_category)
        self.logger.info("Output directory: %s", self.output_dir)
        self.logger.info(("Options:\n" + opt_prettyprint_str).replace("\n", "\n           "))

        from distutils.dir_util import copy_tree
        copy_tree("model/", os.path.join(self.output_dir, "code"))

        self.logs = {
            "train_loss" : open(os.path.join(self.output_dir, "train_loss.txt"), "w", buffering=1),
            "val_loss" : open(os.path.join(self.output_dir, "val_loss.txt"), "w", buffering=1),
            "learning_rate" : open(os.path.join(self.output_dir, "learning_rate.txt"), "w", buffering=1)
        }

    def init_eval(self):
        self.eval()
        self.to(device=self.device)

    def generate_path(self, *rpath):
        return os.path.join(self.output_dir, *rpath)

    def save_checkpoint(self, path):
        torch.save({
            'state_dict': self.state_dict(),
            'optimizer' : self.optimizer.state_dict(),
        }, path)

    def load_checkpoint(self, path):
        checkpoint = torch.load(path)
        self.load_state_dict(checkpoint['state_dict'])
        self.optimizer.load_state_dict(checkpoint['optimizer'])

    def save_model(self, path):
        torch.save(self.state_dict(), path)

    def load_model(self, path):
        self.load_state_dict(torch.load(path))

    def set_criterion(self):
        class RobustLoss(nn.Module):
            def __init__(self, device):
                nn.Module.__init__(self)
                from robust_loss_pytorch.adaptive import AdaptiveLossFunction
                self.lossfn = AdaptiveLossFunction(num_dims = 1, float_dtype=np.float32, device=device)

            def forward(self, predicted, target):
                return torch.mean(self.lossfn.lossfun(predicted - target))

        self.criterion = RobustLoss(device=self.device)

    def set_optimizer_and_scheduler(self):
        config = self.opt["config"]
        self.optimizer = torch.optim.Adam(self.parameters(), lr=config['lr'], weight_decay=config["weight_decay"])
        self.logger.info("[Adam Optimizer] lr = %f, weight_decay = %f" % (config["lr"], config["weight_decay"]))

        self.scheduler = None
        if config["scheduler"]:
            if config["scheduler"]["type"] == "ReduceLROnPlateau":
                factor = config["scheduler"]["factor"]
                patience = config["scheduler"]["patience"]
                cooldown = config["scheduler"]["cooldown"]
                min_lr = config["scheduler"]["min_lr"]
                threshold = config["scheduler"]["threshold"]
                self.scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(self.optimizer, 'min', factor=factor, patience=patience, cooldown=cooldown, min_lr=min_lr, threshold=threshold, verbose=True)
                self.logger.info("[ReduceLROnPlateau] factor = %f, patience = %d, cooldown = %d, min_lr = %.e, threshold = %.e" % (factor, patience, cooldown, min_lr, threshold))

        total_params = sum(p.numel() for p in self.parameters())
        train_params = sum(p.numel() for p in self.parameters() if p.requires_grad)
        self.logger.info("Total Parameters: %d", total_params)
        self.logger.info("Trainable Parameters: %d", train_params)
        self.logger.info(self)

    def set_data_loader(self, trainset, valset):
        config = self.opt["config"]

        self.logger.info("Started preprocessing dataset")
        from .dataset import PreprocessedDataset
        self.train_dataset = PreprocessedDataset(radius_cutoff=config["radius_cutoff"], device=self.device, logger=self.logger)
        self.train_dataset.add_dataset(trainset, training_filters=True)
        self.train_dataset.summary()

        self.validation_dataset = PreprocessedDataset(radius_cutoff=config["radius_cutoff"], device=self.device, logger=self.logger)
        self.validation_dataset.add_dataset(valset, training_filters=False)
        self.validation_dataset.summary()

        from torch.utils.data import DataLoader
        self.train_loader = DataLoader(dataset=self.train_dataset, batch_size=config["batch_size"], collate_fn=self.__collate_fn, shuffle=True)
        self.validation_loader = DataLoader(dataset=self.validation_dataset, batch_size=1, collate_fn=self.__collate_fn, shuffle=False)
        
        self.logger.info("Number of training samples: %d" % (len(self.train_dataset)))
        self.logger.info("Number of validation samples: %d" % (len(self.validation_dataset)))
        self.logger.info("Finished preprocessing dataset")

    @staticmethod
    def __collate_fn(batch):
        collated_batch = {
            "primary_rt" : [sample["primary_rt"] for sample in batch],
            "primary_props_rt" : [sample["primary_props_rt"] for sample in batch],
            "primary_mt" : [sample["primary_mt"] for sample in batch],
            "primary_props_mt" : [sample["primary_props_mt"] for sample in batch],
            "position" : [sample["position"] for sample in batch],
            "neighbors" : [sample["neighbors"] for sample in batch],
            "edge_features" : [sample["edge_features"] for sample in batch],
            "edge_features_rev" : [sample["edge_features_rev"] for sample in batch],
            "target" : torch.stack([sample["target"] for sample in batch])
        }
        return collated_batch

    def generate_input_and_target(self, sample):
        input_data = {
            "primary_rt" : sample["primary_rt"],
            "primary_props_rt" : sample["primary_props_rt"],
            "primary_mt" : sample["primary_mt"],
            "primary_props_mt" : sample["primary_props_mt"],
            "position" : sample["position"],
            "neighbors" : sample["neighbors"],
            "edge_features" : sample["edge_features"],
            "edge_features_rev" : sample["edge_features_rev"],
        }
        return input_data, sample["target"]

    def learn(self, verbose=False):
        self.to(device=self.device)

        config = self.opt["config"]
        best_val_stats, best_selection_value = None, np.Inf
        for epoch_no in range(config["epochs"]):
            accumlated_train_loss = 0.0
            accumulated_additive_loss = 0.0
            num_samples = 0
            
            if epoch_no < config["embed_freeze_epochs"]:
                self.embeddingNet.disable_training()
            else:
                self.embeddingNet.enable_training()

            self.train() # set to training mode (validation might have set to eval mode)
            for data in self.train_loader:
                indata, target = self.generate_input_and_target(data)
                _, additive_loss, predicted = self.forward(**indata)

                batch_size = target.shape[0]
                loss = self.criterion(predicted, target)
                accumlated_train_loss += loss.item() * batch_size

                if additive_loss is not None:
                    additive_loss = torch.mean(additive_loss)
                    accumulated_additive_loss += additive_loss.item()
                    loss += config["symmetry_loss_scale"] * additive_loss
              
                num_samples += batch_size

                self.optimizer.zero_grad()
                loss.backward()
                self.optimizer.step()

            avg_epoch_train_loss = accumlated_train_loss / num_samples
            avg_epoch_additive_loss = accumulated_additive_loss / num_samples

            # validate
            validation_stats = self.evaluate(self.validation_loader)
            avg_epoch_validation_loss = validation_stats["avg_loss"]
            selection_value = validation_stats["ddG_mae"]

            if self.scheduler is not None:
                self.scheduler.step(selection_value)

            LRs = [param_group['lr'] for param_group in self.optimizer.param_groups]
            self.logs['learning_rate'].write(", ".join([str(lr) for lr in LRs]) + '\n')
            self.logs['train_loss'].write(str(avg_epoch_train_loss) + '\n')
            self.logs['val_loss'].write(str(avg_epoch_validation_loss) + '\n')

            is_best_model = best_selection_value - selection_value > config["best_model_threshold"]
            if is_best_model:
                best_selection_value = selection_value
                best_val_stats = validation_stats
                self.save_model(self.generate_path("best_val_model.pth"))

            if verbose:
                self.logger.info("[%d] Epoch completed", epoch_no)
                self.logger.info("\tTraining Stats:")
                self.logger.info("\t\taverage loss: %.4f, avg. additive loss: %.2f", avg_epoch_train_loss, avg_epoch_additive_loss)
                self.logger.info("\tValidation Stats:")
                if is_best_model:
                    self.logger.info("\t\tCURRENT BEST MODEL")
                self.print_stats(validation_stats)
            self.save_checkpoint(self.generate_path("checkpoint.ckpt"))

        self.logger.info("[best model] selection value: %f", best_selection_value)
        self.print_stats(best_val_stats)
        self.save_model(self.generate_path("final_model.pth"))

    def predict(self, sample):
        indata, _ = self.generate_input_and_target(self.__collate_fn([sample]))
        return self.forward(**indata)

    def print_stats(self, stats):
        self.logger.info("\t\tavg loss: %.4f", stats["avg_loss"])
        self.logger.info("\t\tddG PCC: %.2f, ddG r2: %.2f, ddG MAE: %.2f, ddG RMSE: %.2f"
            % (stats["ddG_pcc"], stats["ddG_r2"], stats["ddG_mae"], stats["ddG_rmse"]))
        self.logger.info("\t\t[two class] accuracy: %.2f, mcc: %.2f" % (stats["overall_2acc"], stats["overall_2mcc"]))
        self.logger.info("\t\tstab accuracy: %.2f, stab pcc: %.2f, stab r2: %.2f" % (stats["stab_2acc"], stats["stab_2pcc"], stats["stab_2r2"]))
        self.logger.info("\t\tdestab accuracy: %.2f, destab pcc: %.2f, destab r2: %.2f" % (stats["destab_2acc"], stats["destab_2pcc"], stats["destab_2r2"]))
        self.logger.info("\t\thighly stab auc: %.2f, highly destab auc: %.2f" % (stats["highly_stab_auc"], stats["highly_destab_auc"]))
        self.logger.info("\t\t[three class] accuracy: %.2f, mcc: %.2f" % (stats["overall_3acc"], stats["overall_3mcc"]))
        self.logger.info("\t\tstab accuracy: %.2f, stab pcc: %.2f, stab r2: %.2f" % (stats["stab_3acc"], stats["stab_3pcc"], stats["stab_3r2"]))
        self.logger.info("\t\tneutral accuracy: %.2f, neutral pcc: %.2f, neutral r2: %.2f" % (stats["neutral_3acc"], stats["neutral_3pcc"], stats["neutral_3r2"]))
        self.logger.info("\t\tdestab accuracy: %.2f, destab pcc: %.2f, destab r2: %.2f" % (stats["destab_3acc"], stats["destab_3pcc"], stats["destab_3r2"]))

    @staticmethod
    def evaluate_base(predicted_values, target_values):
        predicted_values = np.asarray(predicted_values)
        target_values = np.asarray(target_values)
        num_samples = len(target_values)
        assert(len(predicted_values) == len(target_values))

        from sklearn.metrics import matthews_corrcoef, accuracy_score, r2_score

        # overall metrics
        overall_ddG_mae = np.mean(np.abs(target_values - predicted_values))
        overall_ddG_rmse = np.sqrt(np.mean((target_values - predicted_values)**2))
        overall_pcc = np.corrcoef(target_values, predicted_values)[0][1]
        overall_r2 = r2_score(target_values, predicted_values)

        # two class metrics
        target_labels = target_values > 0
        predicted_labels = predicted_values > 0

        overall_2acc = accuracy_score(target_labels, predicted_labels)
        overall_2mcc = matthews_corrcoef(target_labels, predicted_labels)

        stab_indices = target_labels
        destab_indices = np.invert(stab_indices)

        stab_2acc = accuracy_score(target_labels[stab_indices], predicted_labels[stab_indices])
        stab_2pcc = np.corrcoef(target_values[stab_indices], predicted_values[stab_indices])[0][1]
        stab_2r2 = r2_score(target_values[stab_indices], predicted_values[stab_indices])
        destab_2acc = accuracy_score(target_labels[destab_indices], predicted_labels[destab_indices])
        destab_2pcc = np.corrcoef(target_values[destab_indices], predicted_values[destab_indices])[0][1]
        destab_2r2 = r2_score(target_values[destab_indices], predicted_values[destab_indices])

        def roc_auc(predictions, target_labels):
            from sklearn.metrics import roc_curve, auc
            fpr, tpr, _ = roc_curve(target_labels, predictions)
            auc = auc(fpr, tpr)
            return auc

        highly_destab_auc = roc_auc(predicted_values, target_values > 1)
        highly_stab_auc = roc_auc(-predicted_values, target_values < -1)

        # three class metrics
        STAB_LABEL = -1
        NEUTRAL_LABEL = 0
        DESTAB_LABEL = 1

        class3_truth = np.zeros(num_samples, dtype=np.int8)
        class3_truth[target_values < -0.5] = STAB_LABEL
        class3_truth[target_values > 0.5] = DESTAB_LABEL

        class3_predicted = np.zeros(num_samples, dtype=np.int8)
        class3_predicted[predicted_values < -0.5] = STAB_LABEL
        class3_predicted[predicted_values > 0.5] = DESTAB_LABEL

        overall_3acc = accuracy_score(class3_truth, class3_predicted)
        overall_3mcc = matthews_corrcoef(class3_truth, class3_predicted)

        stab_indices = class3_truth == STAB_LABEL
        stab_3acc = accuracy_score(class3_truth[stab_indices], class3_predicted[stab_indices])
        stab_3pcc = np.corrcoef(predicted_values[stab_indices], target_values[stab_indices])[0][1]
        stab_3r2 = r2_score(target_values[stab_indices], predicted_values[stab_indices])

        neutral_indices = class3_truth == NEUTRAL_LABEL
        if np.sum(neutral_indices) > 0:
            neutral_3acc = accuracy_score(class3_truth[neutral_indices], class3_predicted[neutral_indices])
            neutral_3pcc = np.corrcoef(predicted_values[neutral_indices], target_values[neutral_indices])[0][1]
            neutral_3r2 = r2_score(target_values[neutral_indices], predicted_values[neutral_indices])
        else:
            neutral_3acc = neutral_3pcc = neutral_3r2 = np.nan

        destab_indices = class3_truth == DESTAB_LABEL
        destab_3acc = accuracy_score(class3_truth[destab_indices], class3_predicted[destab_indices])
        destab_3pcc =  np.corrcoef(predicted_values[destab_indices], target_values[destab_indices])[0][1]
        destab_3r2 = r2_score(target_values[destab_indices], predicted_values[destab_indices])

        stats = {}
        stats["ddG_pcc"] = overall_pcc
        stats["ddG_r2"] = overall_r2
        stats["ddG_mae"] = overall_ddG_mae
        stats["ddG_rmse"] = overall_ddG_rmse
        stats["overall_2acc"] = overall_2acc * 100
        stats["overall_2mcc"] = overall_2mcc
        stats["stab_2acc"] = stab_2acc * 100
        stats["stab_2pcc"] = stab_2pcc
        stats["stab_2r2"] = stab_2r2
        stats["highly_stab_auc"] = highly_stab_auc
        stats["destab_2acc"] = destab_2acc * 100
        stats["destab_2pcc"] = destab_2pcc
        stats["destab_2r2"] = destab_2r2
        stats["highly_destab_auc"] = highly_destab_auc
        stats["overall_3acc"] = overall_3acc * 100
        stats["overall_3mcc"] = overall_3mcc
        stats["stab_3acc"] = stab_3acc * 100
        stats["stab_3pcc"] = stab_3pcc
        stats["stab_3r2"] = stab_3r2
        stats["neutral_3acc"] = neutral_3acc * 100
        stats["neutral_3pcc"] = neutral_3pcc
        stats["neutral_3r2"] = neutral_3r2
        stats["destab_3acc"] = destab_3acc * 100
        stats["destab_3pcc"] = destab_3pcc
        stats["destab_3r2"] = destab_3r2
        return stats

    def evaluate(self, dataset):
        self.eval()

        from torch.utils.data import DataLoader
        if not isinstance(dataset, DataLoader):
            dataset = DataLoader(dataset=dataset, batch_size=1, collate_fn=self.__collate_fn, shuffle=False)
        
        accumulated_loss = 0
        num_samples = len(dataset)
        predicted_values = np.zeros(num_samples)
        target_values = np.zeros(num_samples)
        for idx, data in enumerate(dataset):
            indata, target = self.generate_input_and_target(data)
            _, _, predicted = self.forward(**indata)
            accumulated_loss += self.criterion(predicted, target)

            predicted_values[idx] = predicted.item()
            target_values[idx] = target.item()

        stats = self.evaluate_base(predicted_values, target_values)
        stats["avg_loss"] = accumulated_loss.item() / len(dataset)
        return stats        

    