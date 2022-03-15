#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
The main starting point of the Male pedigree toolbox pipeline
Author: Bram van Wersch

NOTE:
IMPORTANT: Python3.6 or higher needs to be used because some of the modules rely on sorted dictionaries
"""

# library imports
import argparse
import logging
import random
from typing import Union
import time
import datetime
import sys
import os
from sys import argv
from pathlib import Path
import re

# own imports
from MalePedigreeToolbox import utility
from MalePedigreeToolbox import thread_termination

if sys.version_info <= (3, 6):
    raise SystemExit(f"Minimal python version required is 3.6, detected version"
                     f" {'.'.join(map(str, list(sys.version_info)[:3]))}")

LOG: Union[logging.Logger, None] = None
FORCE: bool = False
RUNNING_GUI: bool = True


class CustomArgParser(argparse.ArgumentParser):
    """User as custom argument parser to make sure errors are called with the gui"""

    def error(self, message):
        if RUNNING_GUI:
            raise argparse.ArgumentError(None, f'{self.prog}: error: {message}\n')
        else:
            super().error(message)


def command_line_main():
    """Called from the if __name__ == __main__"""
    # make sure to set certain values that are only applicable when calling script directly
    thread_termination.disable_thread_termination()
    global RUNNING_GUI
    RUNNING_GUI = False
    main()


@thread_termination.ThreadTerminable
def main(*arguments, is_gui=False):
    """Main entry point"""
    global RUNNING_GUI

    if not is_gui:
        thread_termination.disable_thread_termination()
        RUNNING_GUI = False

    parser = get_argument_parser()

    # for testing purposes without popen and gui
    arguments = argv[1:] if len(arguments) == 0 else arguments
    if '-f' in arguments or '--force' in arguments:
        global FORCE
        FORCE = True
    if "--version" in arguments:
        with open(Path(__file__).parent / "__init__.py") as f:
            version_match = re.search(
                r"^__version__ = ['\"]([^'\"]*)['\"]", f.read(), re.M
            )
            if version_match:
                print(f"MalePedigreeToolbox {version_match.group(1)}")
                return
            raise RuntimeError("Failed to find version string")
    try:
        name_space = parser.parse_args(arguments)
    except utility.MalePedigreeToolboxError as e:
        parser.error(str(e))
        return

    if name_space.subcommand is None:
        parser.error("No subcommand provided")  # noqa
    log_file_name = get_log_file_name(name_space.log_name, name_space.outdir)

    start_time = time.time()
    # needs to be after parsing!!!!
    setup_logger(start_time, name_space.log_level, log_file_name, name_space.outdir)
    full_command = ' '.join(argv)
    LOG.info(f"Running MalePedigreeToolbox with command '{full_command}' from"
             f" {'command line ' if not RUNNING_GUI else 'GUI'}.")

    # all imports are locally to make sure libraries get only loaded when a certain command is chosen
    if not name_space.subcommand:
        parser.print_help()
        raise SystemExit("Please provide subcommand")
    elif name_space.subcommand == "distances":
        LOG.info("Loading libraries...")
        from MalePedigreeToolbox import distances
        distances.main(name_space)
    elif name_space.subcommand == "mut_diff":
        LOG.info("Loading libraries...")
        from MalePedigreeToolbox import mutation_diff
        mutation_diff.main(name_space)
    elif name_space.subcommand == "draw_pedigrees":
        LOG.info("Loading libraries...")
        from MalePedigreeToolbox import draw_pedigrees
        draw_pedigrees.main(name_space)
    elif name_space.subcommand == "ped_mut_graph":
        LOG.info("Loading libraries...")
        from MalePedigreeToolbox import infer_pedigree_mutations
        infer_pedigree_mutations.main(name_space)
    elif name_space.subcommand == "all":
        LOG.info("Loading libraries...")
        from MalePedigreeToolbox import run_it_all
        run_it_all.main(name_space)
    elif name_space.subcommand == "simulate":
        LOG.info("Loading libraries...")
        from MalePedigreeToolbox.prediction_code import simulate_mutations
        simulate_mutations.main(name_space)
    elif name_space.subcommand == "make_models":
        LOG.info("Loading libraries...")
        from MalePedigreeToolbox.prediction_code import make_models
        make_models.main(name_space)
    elif name_space.subcommand == "predict":
        LOG.info("Loading libraries...")
        from MalePedigreeToolbox.prediction_code import classifier_predict
        classifier_predict.main(name_space)
    else:
        # should never reach here unless mistake was made by programmer
        raise SystemExit("Invalid subcommand, fix it!")
    LOG.info(f"The log file can be found at '{name_space.outdir / log_file_name}'")


def get_log_file_name(user_provided_name: str, directory: Path) -> str:
    if user_provided_name != "no_value":
        return user_provided_name
    base_name = "run"
    name = f"{base_name}.log"
    all_present_files = set(p.name for p in Path(directory).glob("*"))
    count = 1
    while name in all_present_files:
        name = f"{base_name}{count}.log"
        count += 1
    return name


def setup_logger(start_time: float, level: str, log_file_name: str, output_directory: Path):
    """
    instantiate the logger
    """
    log_file_loc = output_directory / log_file_name
    # make sure that a potentially existing logfile is emptied
    if os.path.exists(log_file_loc):
        open(log_file_loc, "w").close()

    if level == "basic":
        level = logging.INFO
    elif level == "debug":
        level = logging.DEBUG
    else:
        level = logging.WARNING

    logger = logging.getLogger("mpt")
    handlers = logger.handlers[:]
    for handler in handlers:
        handler.close()
        logger.removeHandler(handler)

    logger.propagate = False
    logger.setLevel(level)

    handler = logging.StreamHandler(sys.stdout)
    formatter = MyFormatter('%(levelname)s %(currentTime)s (%(passedTime)s sec) - %(message)s',
                            starttime=start_time)
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    file_handler = logging.FileHandler(filename=log_file_loc)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    global LOG
    LOG = logger

    LOG.debug('Logger created')


class MyFormatter(logging.Formatter):
    """
    Copied from MultiGeneBlast (my code)
    """
    def __init__(self, fmt, starttime=time.time()):
        logging.Formatter.__init__(self, fmt)
        self._start_time = starttime

    def format(self, record):
        """
        Overwrite of the format function that prints the passed time and adds
        current time to the existing format
        :See: logging.Formatter.format()
        """
        # difference = datetime.datetime.now() - self._start_time
        record.passedTime = "{:.3f}".format(time.time() - self._start_time)
        record.currentTime = datetime.datetime.now().time()
        return super(MyFormatter, self).format(record)


def get_argument_parser() -> argparse.ArgumentParser:
    parser = CustomArgParser(description="Welcome to the male pedigree toolbox kit a toolbox for investigating "
                                         "pedigrees and estimating age based on y-chromosmal markers.",
                             epilog="Created by: Arwin Ralf, Diego Montiel Gonzalez and Bram van Wersch, 2021",)
    parser.add_argument("-ll", "--log_level",
                        help="Configure the logger for printing a certain level of information. 'basic': just simple "
                             "'output', debug: all output including issues, 'silent': no output except warnings",
                        choices=["debug", "basic", "silent"], default='basic', metavar="STRING")
    parser.add_argument("-f", "--force", help="Force file/ directory creations regardless if they exist or not.",
                        action="store_true")  # the presence -f or --force will be checked before parcing
    parser.add_argument("-ln", "--log_name",
                        help="The name of the log file. If no name is provided a unique name will be chosen "
                             "automatically", default="no_value", metavar="STRING")

    # add subparsers for the different program functionalities into subcommand
    subparsers = parser.add_subparsers(dest="subcommand")
    add_pairwise_distance_subparser(subparsers)
    add_mutation_differentation_subparser(subparsers)
    add_infer_pedigree_mutations_parser(subparsers)
    add_dendogram_parser(subparsers)
    add_all_parser(subparsers)
    add_simulation_parser(subparsers)
    add_make_models_parser(subparsers)
    add_predict_parser(subparsers)
    return parser


def add_pairwise_distance_subparser(subparsers):
    pairwise_dist_parser = subparsers.add_parser("distances",
                                                 help="Calultate the pairwise distance for all .tgf files in a folder")
    pairwise_dist_parser.add_argument("-t", "--tgf_folder",  required=True,
                                      help="Folder name containing at least 1 .tgf file.", metavar="PATH",
                                      type=utility.check_tgf_folder)
    pairwise_dist_parser.add_argument("-o", "--outdir", help="Output folder name",
                                      required=True, metavar="PATH",
                                      type=lambda path: utility.check_create_out_folder(path, force=FORCE))


def add_mutation_differentation_subparser(subparsers):
    mutation_parser = subparsers.add_parser("mut_diff",
                                            help="Calculate the number of mutations between all alleles of markers in "
                                                 "a pedigree. Additionaly calculate mutation rates if distances between"
                                                 " individuals in a pedigree are provided.")
    mutation_parser.add_argument("-af", "--allele_file", help="File containing allele frequencies in CSV format",
                                 required=True, metavar="FILE", type=utility.check_in_file)
    mutation_parser.add_argument("-df", "--dist_file", help="File of distances between samples.This can be "
                                                            "generated with the help of the 'distances' command",
                                 metavar="FILE(OPTIONAL)", type=utility.check_in_file)
    mutation_parser.add_argument("-pf", "--prediction_file", help="Add this option if you want a file to be generated"
                                                                  " that can be used by the predict command in order to"
                                                                  " predict generational distance between individuals "
                                                                  "in a pedigree.",
                                 action="store_true")
    mutation_parser.add_argument("-o", "--outdir", help="Output dir for all the resulting files",
                                 type=lambda path: utility.check_create_out_folder(path, force=FORCE),
                                 metavar="PATH", required=True)


def add_infer_pedigree_mutations_parser(subparsers):
    pedigree_parser = subparsers.add_parser("ped_mut_graph",
                                            help="Infer mutation events from pedigrees with incomplete allele data, "
                                                 "in addition draw these pedigrees.")
    pedigree_parser.add_argument("-af", "--allele_file", help="File containing allele frequencies in CSV format",
                                 required=True, metavar="FILE", type=utility.check_in_file)
    pedigree_parser.add_argument("-t", "--tgf_folder",  required=True,
                                 help="Folder name containing at least 1 .tgf file.", metavar="PATH",
                                 type=utility.check_tgf_folder)
    pedigree_parser.add_argument("-mm", "--minimum_mutations",
                                 help="Minimum mutations that need to be present in order for the pedigree to be drawn"
                                      " (default=1).",
                                 default=1, type=int, metavar="INT(OPTIONAL)")
    pedigree_parser.add_argument("-o", "--outdir", required=True,
                                 help="Folder path to store all outputs", metavar="PATH",
                                 type=lambda path: utility.check_create_out_folder(path, force=FORCE))


def add_dendogram_parser(subparsers):
    dendogram_parser = subparsers.add_parser("draw_pedigrees",
                                             help="Draw plots based on mutation distances in a pedigree, this can be "
                                                  "dendrograms or multi-dimensional scaling plots.")
    dendogram_parser.add_argument("-fm", "--full_marker_csv", help="The file containing full mutation "
                                                                   "differentiations (-fo) generated with the mut_diff"
                                                                   " command",
                                  required=True, metavar="FILE", type=utility.check_in_file)
    dendogram_parser.add_argument("-mr", "--marker_rates", help="File with mutation rates of all markers present in "
                                                                "full marker file. The expected format is a csv file "
                                                                "with 2 columns 1. marker 2. rate. This will give more "
                                                                "accurate dendrograms. Leave this field empty to assume"
                                                                " the same mutation rate for all markers.",
                                  metavar="FILE", type=utility.check_in_file)
    dendogram_parser.add_argument("-t", "--type", help="The plot type you want. You can choose 'dendrogram' or 'MDS' or"
                                                       "'both'. The default is 'dendrogram'",
                                  metavar="NAME(OPTIONAL)", choices=['dendrogram', 'MDS', 'both'], default='dendrogram')
    dendogram_parser.add_argument("-o", "--outdir", help="Folder path to store all outputs",
                                  required=True, metavar="PATH",
                                  type=lambda path: utility.check_create_out_folder(path, force=FORCE))
    dendogram_parser.add_argument("-c", "--clusters", help="The expected number of clusters for all pedigrees. This can"
                                                           " be a single value to get the same number of clusters for"
                                                           " all pedigrees or a list that will assign that number of "
                                                           "clusters to a pedigree until the list is empty. The input "
                                                           "can be a file with space seperated integers or typed out "
                                                           "on the command line. If no value is provided the optimal "
                                                           "clustering is calculated based on silhouette score.",
                                  metavar="INT / FILE(OPTIONAL)", nargs="+", type=utility.check_file_int)
    dendogram_parser.add_argument("-rs", "--random_state",
                                  help="Innitialize the MDS and clustering algorithms with a random state for "
                                       "consistent figures between runs.", default=random.randint(0, 1_000_000),
                                  type=int, metavar="INT(OPTIONAL)")
    dendogram_parser.add_argument("-md", "--min_dist",
                                  help="The minimum distance an arm of the dendrogram is drawn. Since distances can be "
                                       "often 0 dendograms can look a bit strange. (default = 0.1)", default=0.1,
                                  metavar="FLOAT(OPTIONAL)", type=float)


def add_all_parser(subparsers):
    all_parser = subparsers.add_parser("all", help="Performs all of the commands in order, first calculating distance,"
                                                   "than mutation differentiation, followed by dendograms for "
                                                   "mutation differentiation rates and finaly caclulates mutation "
                                                   "rates from pedigrees.")
    all_parser.add_argument("-t", "--tgf_folder",  required=True,
                            help="Folder name containing at least 1 .tgf file.", metavar="PATH",
                            type=utility.check_tgf_folder)
    all_parser.add_argument("-af", "--allele_file", help="File containing allele frequencies in CSV format",
                            required=True, metavar="FILE", type=utility.check_in_file)
    all_parser.add_argument("-o", "--outdir", help="Output directory where all information will end up",
                            required=True, metavar="PATH",
                            type=lambda path: utility.check_create_out_folder(path, force=FORCE))
    all_parser.add_argument("-pf", "--prediction_file", help="Add this option if you want a file to be generated"
                                                             " that can be used by the predict command in order to"
                                                             " predict generational distance between individuals "
                                                             "in a pedigree.", action="store_true")
    all_parser.add_argument("-md", "--min_dist",
                            help="The minimum distance an arm of the dendrogram is drawn. Since distances can be "
                                 "often 0 dendograms can look a bit strange. (default = 0.1)", default=0.1,
                            metavar="FLOAT(OPTIONAL)", type=float)
    all_parser.add_argument("-tp", "--type", help="The plot type you want. You can choose 'dendrogram' or 'MDS' or"
                                                  "'both'. The default is 'dendrogram'",
                            metavar="NAME(OPTIONAL)", choices=['dendrogram', 'MDS', 'both'], default='dendrogram')
    all_parser.add_argument("-rs", "--random_state",
                            help="Innitialize the MDS and clustering algorithms with a random state for "
                                   "consistent figures between runs.", default=random.randint(0, 1_000_000),
                            type=int, metavar="INT(OPTIONAL)")
    all_parser.add_argument("-mr", "--marker_rates", help="File with mutation rates of all markers present in "
                                                          "full marker file. The expected format is a csv file "
                                                          "with 2 columns 1. marker 2. rate. This will give more "
                                                          "accurate dendrograms. Leave this field empty to assume"
                                                          " the same mutation rate for all markers.",
                            metavar="FILE(OPTIONAL)", type=utility.check_in_file)

    all_parser.add_argument("-c", "--clusters", help="The expected number of clusters for all pedigrees. This can"
                                                     " be a single value to get the same number of clusters for"
                                                     " all pedigrees or a list that will assign that number of "
                                                     "clusters to a pedigree until the list is empty. The input "
                                                     "can be a file with space seperated integers or typed out "
                                                     "on the command line. If no value is provided the optimal "
                                                     "clustering is calculated based on silhouette score.",
                                  metavar="INT / FILE(OPTIONAL)", nargs="+", type=utility.check_file_int)

    all_parser.add_argument("-mm", "--minimum_mutations",
                            help="Minimum mutations that need to be present in order for the pedigree to be drawn"
                                 " (default=1).",
                            default=1, type=int, metavar="INT(OPTIONAL)")


def add_simulation_parser(subparsers):
    simulation_parser = subparsers.add_parser("simulate",
                                              help="Simulate data for creating classification models based on mutation"
                                                   " rates of markers that are provided in a file. Data is simulated"
                                                   " for a given number of generations and samples. All generations"
                                                   " are independently simulated.")
    simulation_parser.add_argument("-i", "--input", dest="input",
                                   help="Excel or csv file containing mutation rates", metavar="FILE", required=True,
                                   type=utility.check_in_file)
    simulation_parser.add_argument("-o", "--outdir", dest="outdir",
                                   help="Folder for table of simulated mutations in csv format", metavar="PATH",
                                   required=True,
                                   type=lambda path: utility.check_create_out_folder(path, force=FORCE))

    simulation_parser.add_argument("-n", "--num_sim",
                                   help="Number of simulations/samples simulated per generation",
                                   metavar="INT", required=True, type=int)

    simulation_parser.add_argument("-g", "--generations", help="Number of generations to simulate",
                                   metavar="INT", required=True, type=int)


def add_make_models_parser(subparsers):
    make_models_parser = subparsers.add_parser("make_models",
                                               help="Create one or more machine learning models for mutation rates."
                                                    " Models can be tested with a number of different hyper parameter "
                                                    "Keep in mind that making these models can take quite alot of time"
                                                    " depending on the complexity of the model and the size of the "
                                                    "input data set.")
    make_models_parser.add_argument("-i", "--input", dest="input",
                                    help="Csv file formatted in the way that the 'simulate' command outputs its files. "
                                         "Keep in mind that files that are too small (< 10.000 observations) will leed"
                                         " to very poor models.",
                                    metavar="FILE", required=True, type=utility.check_in_file)
    make_models_parser.add_argument("-o", "--outdir", dest="outdir",
                                    help="Folder where the final models, evaluation parameters and evaluation figures"
                                         "are places", metavar="FOLDER", required=True,
                                    type=lambda path: utility.check_create_out_folder(path, force=FORCE))
    make_models_parser.add_argument("-mt", "--model_types",
                                    help="The different types of models to create a predictor for. The predictor is "
                                         "always a classifier.",
                                    choices=["KNN", "LDA", "logistic", "QDA", "RF", "Gaussian", "MLP", "SVM"],
                                    required=True, nargs="+", metavar="STRING(OPTIONAL)")
    make_models_parser.add_argument("-cv", "--cv_splits",
                                    help="Cross validation splits. Leaf one out cross validation is not supported. In "
                                         "total as many final models as cv splits will be created. If 1 is provided a "
                                         "single 80/20 split will be used. Default = 5.",
                                    default=5, type=int, metavar="INT(OPTIONAL)")
    make_models_parser.add_argument("-hpc", "--hyper_parameter_choices",
                                    help="The total number of times random parameters will be rolled from a predefined"
                                         " parameter space. Default = 1000.", default=1000, type=int,
                                    metavar="INT(OPTIONAL)")
    make_models_parser.add_argument("-pcv", "--parameter_cv",
                                    help="Cross validation used for hyper parameter tuning of each model. A randomized"
                                         " search is used for the tuning. Default = 2.", default=2, type=int,
                                    metavar="INT(OPTIONAL)")
    make_models_parser.add_argument("-c", "--cpus",
                                    help="The number of cpu's to use when calculating models. Use -1 to use all "
                                         "available cpu's. Default = 1", default=1, type=int,
                                    metavar="INT(OPTIONAL)")


def add_predict_parser(subparsers):
    predict_parser = subparsers.add_parser("predict",
                                           help="Predict the generational distance between 2 individuals based on the "
                                                "number of mutations between them.")
    predict_parser.add_argument("-i", "--input",
                                help="Mutation rates table to predict. Columns are marker names, Rows as"
                                     " samples in CSV or TSV format. This file can be generated from pedigrees trough "
                                     "the mut_diff command "
                                     "instead.", metavar="FILE", required=True,
                                type=utility.check_in_file)
    predict_parser.add_argument("-m", "--model",
                                help="A path leading to a joblib dumped model. Instead of this parameter you can "
                                     "also choose one of the precalculated models see --predefined_models",
                                metavar="PATH/STR", type=utility.check_in_file, default=None)
    predict_parser.add_argument("-tf", "--training_file",
                                help="The file that was used to train the data, to make sure that the order of the"
                                     " input file is the same as the order expected by the model. This only needs to"
                                     " be specified when using the --model option.",
                                metavar="FILE", type=utility.check_in_file, default=None)
    predict_parser.add_argument("-pm", "--predefined_model",
                                help="Choose an already calculated model based on a 5 million datapoints dataset. "
                                     "The names of the models refer to sets of markers used when creating the models."
                                     " If you want to know if your set applies to a model try it out and an error will"
                                     " be raised if not all required markers are present.",
                                choices=["RMPLEX", "PPY23", "YFP", "PPY23_RMPLEX", "YFP_RMPLEX"], default=None)
    predict_parser.add_argument("-p", "--plots", help="Specify this argument if all comparissons in the input file "
                                                      "should be plotted, showing the age range at 85, 95 and 99 %% "
                                                      "confidence. If a lot of values are predicted a lot of plots are "
                                                      "generated, so be warned, it can take a while.",
                                action="store_true")
    predict_parser.add_argument("-o", "--outdir", help="Output directory", metavar="PATH", required=True,
                                type=lambda path: utility.check_create_out_folder(path, force=FORCE))


if __name__ == '__main__':
    # entry point for the tool
    command_line_main()
