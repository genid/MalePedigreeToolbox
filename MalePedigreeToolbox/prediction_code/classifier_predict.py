# -*- coding: utf-8 -*-
"""
Predict one or more new samples based on a provided or pre-computed model

@author: Bram
"""

from typing import List, TYPE_CHECKING, Union
import logging
import os
import sys
from pathlib import Path

import joblib
import pandas
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# own imports
from MalePedigreeToolbox import utility
from MalePedigreeToolbox import thread_termination

if TYPE_CHECKING:
    import argparse
    from sklearn.model_selection import RandomizedSearchCV


matplotlib.use('Agg')


if getattr(sys, 'frozen', False):
    data_path = Path(sys._MEIPASS)
elif __file__:
    data_path = Path(os.path.dirname(__file__))
else:
    raise SystemExit("Can not find application path.")


LOG = logging.getLogger("mpt")


# marker orders for the pre-made models
RMPLEX_MARKER_LIST = ['DYS713', 'DYF1001', 'DYS626', 'DYS711', 'DYS1010', 'DYS612', 'DYF1002', 'DYS724', 'DYS1005',
                      'DYF403S1b', 'DYS449', 'DYS1012', 'DYS518', 'DYF404S1', 'DYS442', 'DYF403S1a', 'DYS576', 'DYS547',
                      'DYS712', 'DYS570', 'DYS1007', 'DYF393S1', 'DYF1000', 'DYS627', 'DYS1003', 'DYF399S1', 'DYS1013',
                      'DYF387S1', 'DYR88', 'DYS526b']
PPY23_MARKER_LIST = ['DYS392', 'DYS643', 'DYS549', 'DYS438', 'DYS458', 'DYS385', 'DYS19', 'DYS481', 'DYS635', 'DYS576',
                     'DYS391', 'DYS448', 'DYS389I', 'DYS390', 'DYS533', 'DYS456', 'DYS570', 'DYS439', 'YGATAH4',
                     'DYS393', 'DYS437', 'DYS389II']
PPY23_RMPLEX_MARKER_LIST = ['DYS576', 'DYS449', 'DYS391', 'DYS713', 'DYS712', 'DYS1010', 'DYF403S1b', 'DYS1005',
                            'DYR88', 'DYS533', 'DYS393', 'DYS1013', 'DYS392', 'DYS570', 'DYS711', 'DYS1007', 'DYS481',
                            'DYS643', 'DYS442', 'DYS547', 'DYF1001', 'DYS635', 'DYF1000', 'DYS438', 'DYS458',
                            'YGATAH4', 'DYF404S1', 'DYS626', 'DYS448', 'DYS390', 'DYS437', 'DYS385', 'DYS1003',
                            'DYS549', 'DYF1002', 'DYS389I', 'DYF399S1', 'DYS518', 'DYS1012', 'DYS612', 'DYF387S1',
                            'DYS439', 'DYS19', 'DYS526b', 'DYS389II', 'DYS456', 'DYS724', 'DYF403S1a', 'DYF393S1',
                            'DYS627']
YFP_MARKER_LIST = ['DYS392', 'DYS570', 'DYS438', 'DYS456', 'DYF387S1', 'DYS19', 'DYS449', 'DYS576', 'DYS439', 'DYS460',
                   'DYS448', 'DYS389II', 'DYS391', 'DYS437', 'DYS389I', 'DYS393', 'DYS635', 'DYS627', 'DYS458',
                   'DYS518', 'YGATAH4', 'DYS481', 'DYS533', 'DYS390', 'DYS385']
YFP_RMPLEX_MARKER_LIST = ['DYS392', 'DYS460', 'DYF1002', 'DYS456', 'DYS627', 'DYR88', 'DYS391', 'DYS626', 'DYS1013',
                          'DYS1007', 'DYS1003', 'DYS389I', 'DYS442', 'DYS481', 'YGATAH4', 'DYF393S1', 'DYF403S1b',
                          'DYS448', 'DYS547', 'DYS713', 'DYF387S1', 'DYS635', 'DYS711', 'DYS19', 'DYS533', 'DYS438',
                          'DYF404S1', 'DYS390', 'DYS449', 'DYF403S1a', 'DYS389II', 'DYF399S1', 'DYS1010', 'DYS724',
                          'DYS437', 'DYS570', 'DYS518', 'DYS612', 'DYF1000', 'DYF1001', 'DYS526b', 'DYS1005', 'DYS576',
                          'DYS1012', 'DYS393', 'DYS458', 'DYS439', 'DYS712', 'DYS385']

MARKER_MAPPING = {"RMPLEX": RMPLEX_MARKER_LIST, "PPY23": PPY23_MARKER_LIST, "YFP": YFP_MARKER_LIST,
                  "PPY23_RMPLEX": PPY23_RMPLEX_MARKER_LIST, "YFP_RMPLEX": YFP_RMPLEX_MARKER_LIST}


try:
    MODELS = {"RMPLEX": joblib.load(data_path / "models" / "MLPClassifier_RMPlex.joblib"),
              "PPY23": joblib.load(data_path / "models" / "MLPClassifier_PPY23.joblib"),
              "PPY23_RMPLEX": joblib.load(data_path / "models" / "MLPClassifier_PPY23_RMplex.joblib"),
              "YFP": joblib.load(data_path / "models" / "MLPClassifier_YFP.joblib"),
              "YFP_RMPLEX": joblib.load(data_path / "models" / "MLPClassifier_YFP_RMplex.joblib")}
except FileNotFoundError:
    MODELS = {}

# speeds up calculation of prediction ranges a lot, since a lot of the same results are present
PREDICTION_RANGE_CACHE = {}


@thread_termination.ThreadTerminable
def main(name_space: "argparse.Namespace"):
    LOG.info("Started with predicting models")
    df_test = read_input_file(name_space.input)
    model_path = name_space.model
    predefined_model_name = name_space.predefined_model
    training_input_file = name_space.training_file
    output_dir = Path(name_space.outdir)
    user_wants_plots = name_space.plots

    if model_path is not None and predefined_model_name is not None:
        LOG.error("You can not both define a custom model and predefined model")
        raise utility.MalePedigreeToolboxError("You can not both define a custom model and predefined model")
    if model_path is None and predefined_model_name is None:
        LOG.error("Either a custom model has to be specified or a predefined model should be chosen.")
        raise utility.MalePedigreeToolboxError("Either a custom model has to be specified or a predefined model should"
                                               " be chosen.")
    if model_path is not None and training_input_file is None:
        LOG.error("When providing a custom model the training data has to be provided to guarantee correct input"
                  " value order")
        raise utility.MalePedigreeToolboxError("When providing a custom model the trianing data has to be provided to "
                                               "guarantee correct input value order")

    column_order = get_column_order(model_path, predefined_model_name, training_input_file)
    LOG.info("Reading input model")
    model = get_model(model_path, predefined_model_name)
    LOG.info("Making predictions")
    predict_model(df_test, model, column_order, output_dir, user_wants_plots)
    LOG.info("Finished with predicting input")


@thread_termination.ThreadTerminable
def read_input_file(input_file: Path) -> pd.DataFrame:
    ext = os.path.splitext(input_file)[1]
    if "csv" in ext:
        return pd.read_csv(input_file, index_col=0)
    elif "tsv" in ext:
        return pd.read_csv(input_file, index_col=0, sep="\t")
    else:
        raise utility.MalePedigreeToolboxError("Either supply a .csv or .tsv file.")


@thread_termination.ThreadTerminable
def get_column_order(
    model_path: Union[Path, None],
    predefined_model_name: Union[str, None],
    training_input_file: Union[Path, None]
) -> List[str]:
    if model_path is not None:
        with open(training_input_file) as f:
            header_line = f.readline()
        return header_line.split(",")
    else:
        return MARKER_MAPPING[predefined_model_name]


@thread_termination.ThreadTerminable
def get_model(
    model_path: Union[Path, None],
    predefined_model_name: Union[str, None]
) -> "RandomizedSearchCV":
    if model_path is not None:
        return joblib.load(model_path)
    else:
        # models failed to load because executing from executable
        if len(MODELS) == 0:
            LOG.error("Failed to load pre-compiled models. Make sure that these models are present in the executable "
                      "directory in a folder called 'models'.")
            raise utility.MalePedigreeToolboxError("Failed to load pre-compiled models. Make sure that these models "
                                                   "are present in the executable directory in a folder called"
                                                   " 'models'.")
        return MODELS[predefined_model_name]


@thread_termination.ThreadTerminable
def predict_model(
    df_test: pandas.DataFrame,
    model: "RandomizedSearchCV",
    column_order: List[str],
    output_dir: Path,
    user_wants_plots: bool
):
    # any markers not used by the model are dropped
    df_test.drop(df_test.columns.difference(column_order), 1, inplace=True)
    try:
        df_test = df_test[column_order]
    except KeyError as e:
        LOG.error(f"Not all required markers are present in the prediction file. {str(e)}")
        raise utility.MalePedigreeToolboxError(f"Not all required markers are present in the prediction file. {str(e)}")

    x_test = df_test.values

    try:
        model_name = type(model.best_estimator_).__name__
    except AttributeError:
        LOG.error("Expected a 'RandomizedSearchCV' object. If you are unsure please use the make_models parser to"
                  " get the correct joblib object.")
        raise utility.MalePedigreeToolboxError("Expected a 'RandomizedSearchCV' object. If you are unsure please use"
                                               " the make_models parser to get the correct joblib object.")
    y_pred_proba = model.predict_proba(x_test)
    pred_df = pd.DataFrame(y_pred_proba)
    pred_df = pred_df.rename(columns={index: index + 1 for index in range(len(pred_df.columns))})

    if user_wants_plots:
        LOG.info("Started with generating plots")
        create_plots(output_dir, pred_df, df_test.index)

    LOG.info("Estimating probability ranges")
    _99_age_ranges, _95_age_ranges, _85_age_ranges = caclulate_probability_ranges(pred_df)

    pred_df["85 prob. range"] = _85_age_ranges
    pred_df["95 prob. range"] = _95_age_ranges
    pred_df["99 prob. range"] = _99_age_ranges
    pred_df.insert(loc=0, column="sample", value=list(df_test.index))

    LOG.info("Started writing final prediction table.")
    pred_df.to_csv((output_dir / f"{model_name}_predictions.csv"))


@thread_termination.ThreadTerminable
def caclulate_probability_ranges(pred_df):
    _99_age_ranges = []
    _95_age_ranges = []
    _85_age_ranges = []
    prev_total = 0

    for index, (_, prediction) in enumerate(pred_df.iterrows()):
        _85_indexes, _85_total, _95_indexes, _95_total, _99_indexes, _99_total = \
            get_ranges(prediction, len(pred_df.columns))
        _99_age_ranges.append(f"{_99_indexes[0] + 1}-{_99_indexes[1]}({_99_total})")
        _95_age_ranges.append(f"{_95_indexes[0] + 1}-{_95_indexes[1]}({_95_total})")
        _85_age_ranges.append(f"{_85_indexes[0] + 1}-{_85_indexes[1]}({_85_total})")

        # update user periodicall
        total, remainder = divmod(index / len(pred_df.index), 0.05)
        if total != prev_total:
            LOG.info(f"Calculation progress: {round(5 * total)}%...")
            prev_total = total
    return _99_age_ranges, _95_age_ranges, _85_age_ranges


@thread_termination.ThreadTerminable
def create_plots(outdir: Path, pred_df: pd.DataFrame, samples: List[str]):
    x_values = list(range(1, len(pred_df.columns) + 1))

    prev_total = 0
    with PdfPages(outdir / "plots.pdf") as pdf:
        for index, (_, prediction) in enumerate(pred_df.iterrows()):

            _85_indexes, _85_total, _95_indexes, _95_total, _99_indexes, _99_total = \
                get_ranges(prediction, x_values[-1] - 1)
            fig = plt.figure(num=1, clear=True)
            plt.plot(x_values, prediction)

            try:
                pedigree, name1, name2 = samples[index].split("_")
            except ValueError:
                # in case the name is in a different format
                pedigree = samples[index]
                name1 = "name1"
                name2 = "name2"

            plt.fill(
                [_99_indexes[0] + 1, *list(range(_99_indexes[0] + 1, _99_indexes[1] + 2)), _99_indexes[1] + 1],
                [0, *prediction[_99_indexes[0]:_99_indexes[1] + 1], 0],
                label=f"{_99_total}% ci. ({_99_indexes[0] + 1}-{_99_indexes[1] + 1})")
            plt.fill(
                [_95_indexes[0] + 1, *list(range(_95_indexes[0] + 1, _95_indexes[1] + 2)), _95_indexes[1] + 1],
                [0, *prediction[_95_indexes[0]:_95_indexes[1] + 1], 0],
                label=f"{_95_total}% ci. ({_95_indexes[0] + 1}-{_95_indexes[1] + 1})")
            plt.fill(
                [_85_indexes[0] + 1, *list(range(_85_indexes[0] + 1, _85_indexes[1] + 2)), _85_indexes[1] + 1],
                [0, *prediction[_85_indexes[0]:_85_indexes[1] + 1], 0],
                label=f"{_85_total}% ci. ({_85_indexes[0] + 1}-{_85_indexes[1] + 1})")

            plt.title(f"Likelyhood of generational distance between {name1} and {name2} for pedigree {pedigree}")
            plt.xlabel("Number of generations apart")
            plt.ylabel("Probability")
            plt.legend()
            plt.tight_layout()
            pdf.savefig()

            # update user periodically
            total, remainder = divmod(index / len(pred_df.index), 0.01)
            if total != prev_total:
                LOG.info(f"Plotting progress: {round(1 * total)}%...")
                prev_total = total


@thread_termination.ThreadTerminable
def get_ranges(predictions, max_x):
    (_85_indexes, _85_total), (_95_indexes, _95_total), (_99_indexes, _99_total) = \
        get_accumulated_prob_ranges(predictions, [0.99, 0.95, 0.85])
    _85_indexes = (_85_indexes[0], min(_85_indexes[1], max_x))
    _95_indexes = (_95_indexes[0], min(_95_indexes[1], max_x))
    _99_indexes = (_99_indexes[0], min(_99_indexes[1], max_x))
    return _85_indexes, _85_total, _95_indexes, _95_total, _99_indexes, _99_total


@thread_termination.ThreadTerminable
def get_accumulated_prob_ranges(predictions, tresholds):
    key = tuple(predictions)
    if key in PREDICTION_RANGE_CACHE:
        return PREDICTION_RANGE_CACHE[key]
    final_data = []
    for n in range(1, len(predictions) + 1, 1):
        dict_ = {}
        for i in range(0, len(predictions), 1):
            dict_[(i, n+i)] = sum(predictions[i:i + n])
        max_key = max(dict_, key=dict_.get)  # get key with max value
        for index in range(len(tresholds) - 1, -1, -1):
            treshold = tresholds[index]
            if dict_[max_key] >= treshold:
                final_data.append((max_key, round(dict_[max_key] * 100, 2)))
                del tresholds[index]
                if len(tresholds) == 0:
                    PREDICTION_RANGE_CACHE[key] = final_data
                    return final_data
    PREDICTION_RANGE_CACHE[key] = final_data
    return final_data
