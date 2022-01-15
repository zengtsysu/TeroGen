#!/usr/bin/env python
#  coding=utf-8

"""
Script to train a model
"""

import argparse
import os.path
import glob
import itertools as it

import torch
import torch.utils.tensorboard as tbx

import collect_stats_from_model as csfm

import models.model as mm
import models.actions as ma

import utils.chem as uc
import utils.log as ul


class TrainModelPostEpochHook(ma.TrainModelPostEpochHook):

    WRITER_CACHE_EPOCHS = 25

    def __init__(self, output_prefix_path, epochs, validation_sets, lr_scheduler, collect_stats_params,
                 lr_params, collect_stats_frequency, save_frequency, logger=None):
        ma.TrainModelPostEpochHook.__init__(self, logger)

        self.validation_sets = validation_sets
        self.lr_scheduler = lr_scheduler

        self.output_prefix_path = output_prefix_path
        self.save_frequency = save_frequency
        self.epochs = epochs
        self.log_path = collect_stats_params["log_path"]

        self.collect_stats_params = collect_stats_params
        self.collect_stats_frequency = collect_stats_frequency

        self.lr_params = lr_params

        self._writer = None
        if self.collect_stats_frequency > 0:
            self._reset_writer()

    def __del__(self):
        self._close_writer()

    def run(self, model, training_set, epoch):
        if self.collect_stats_frequency > 0 and epoch % self.collect_stats_frequency == 0:
            validation_set = next(self.validation_sets)
            other_values = {"lr": self.get_lr()}

            ma.CollectStatsFromModel(
                model=model, epoch=epoch, training_set=training_set,
                validation_set=validation_set, writer=self._writer, other_values=other_values, logger=self.logger,
                sample_size=self.collect_stats_params["sample_size"]
            ).run()

        self.lr_scheduler.step(epoch=epoch)

        lr_reached_min = (self.get_lr() < self.lr_params["min"])
        if lr_reached_min or self.epochs == epoch \
                or (self.save_frequency > 0 and (epoch % self.save_frequency == 0)):
            model.save(self._model_path(epoch))

        if self._writer and (epoch % self.WRITER_CACHE_EPOCHS == 0):
            self._reset_writer()

        return not lr_reached_min

    def get_lr(self):
        return self.lr_scheduler.optimizer.param_groups[0]["lr"]

    def _model_path(self, epoch):
        return "{}.{}".format(self.output_prefix_path, epoch)

    def _reset_writer(self):
        self._close_writer()
        self._writer = tbx.SummaryWriter(log_dir=self.log_path)

    def _close_writer(self):
        if self._writer:
            self._writer.close()


def main():
    """Main function."""
    params = parse_args()
    lr_params = params["learning_rate"]
    cs_params = params["collect_stats"]
    params = params["other"]

    # ut.set_default_device("cuda")

    model = mm.DecoratorModel.load_from_file(params["input_model_path"])
    optimizer = torch.optim.Adam(model.network.parameters(), lr=lr_params["start"])
    lr_scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=lr_params["step"], gamma=lr_params["gamma"])

    training_sets = load_sets(params["training_set_path"])
    validation_sets = []
    if params["collect_stats_frequency"] > 0:
        validation_sets = load_sets(cs_params["validation_set_path"])

    post_epoch_hook = TrainModelPostEpochHook(
        params["output_model_prefix_path"], params["epochs"], validation_sets, lr_scheduler,
        cs_params, lr_params, collect_stats_frequency=params["collect_stats_frequency"],
        save_frequency=params["save_every_n_epochs"], logger=LOG
    )

    epochs_it = ma.TrainModel(model, optimizer, training_sets, params["batch_size"], params["clip_gradients"],
                              params["epochs"], post_epoch_hook, logger=LOG).run()

    for num, (total, epoch_it) in enumerate(epochs_it):
        for _ in ul.progress_bar(epoch_it, total=total, desc="#{}".format(num)):
            pass  # we could do sth in here, but not needed :)


def load_sets(set_path):
    file_paths = [set_path]
    if os.path.isdir(set_path):
        file_paths = sorted(glob.glob("{}/*.smi".format(set_path)))

    for path in it.cycle(file_paths):  # stores the path instead of the set
        yield list(uc.read_csv_file(path, num_fields=2))


SUBCATEGORIES = ["collect_stats", "learning_rate"]


def parse_args():
    """Parses input arguments."""
    parser = argparse.ArgumentParser(
        description="Train a model on a SMILES file.")

    _add_base_args(parser)
    _add_lr_args(parser)
    csfm.add_stats_args(parser, with_prefix=True, with_required=False)

    args = {k: {} for k in ["other", *SUBCATEGORIES]}
    for arg, val in vars(parser.parse_args()).items():
        done = False
        for prefix in SUBCATEGORIES:
            if arg.startswith(prefix):
                arg_name = arg[len(prefix) + 1:]
                args[prefix][arg_name] = val
                done = True
        if not done:
            args["other"][arg] = val

    # special case
    args["other"]["collect_stats_frequency"] = args["collect_stats"]["frequency"]
    del args["collect_stats"]["frequency"]

    return args


def _add_lr_args(parser):
    parser.add_argument("--learning-rate-start", "--lrs",
                        help="Starting learning rate for training. [DEFAULT: 1E-4]", type=float, default=1E-4)
    parser.add_argument("--learning-rate-min", "--lrmin",
                        help="Minimum learning rate, when reached the training stops. [DEFAULT: 1E-6]",
                        type=float, default=1E-6)
    parser.add_argument("--learning-rate-gamma", "--lrg",
                        help="Ratio which the learning change is changed. [DEFAULT: 0.95]", type=float, default=0.95)
    parser.add_argument("--learning-rate-step", "--lrt",
                        help="Number of epochs until the learning rate changes. [DEFAULT: 1]",
                        type=int, default=1)


def _add_base_args(parser):
    parser.add_argument("--input-model-path", "-i", help="Input model file", type=str, required=True)
    parser.add_argument("--output-model-prefix-path", "-o",
                        help="Prefix to the output model (may have the epoch appended).", type=str, required=True)
    parser.add_argument("--training-set-path", "-s", help="Path to a file with (scaffold, decoration) tuples \
        or a directory with many of these files to be used as training set.", type=str, required=True)
    parser.add_argument("--save-every-n-epochs", "--sen",
                        help="Save the model after n epochs. [DEFAULT: 1]", type=int, default=1)
    parser.add_argument("--epochs", "-e", help="Number of epochs to train. [DEFAULT: 100]", type=int, default=100)
    parser.add_argument("--batch-size", "-b",
                        help="Number of molecules processed per batch. [DEFAULT: 128]", type=int, default=128)
    parser.add_argument("--clip-gradients",
                        help="Clip gradients to a given norm. [DEFAULT: 1.0]", type=float, default=1.0)
    parser.add_argument("--collect-stats-frequency", "--csf",
                        help="Collect statistics every n epochs. [DEFAULT: 0]", type=int, default=0)


if __name__ == "__main__":
    LOG = ul.get_logger(name="train_model")
    main()
