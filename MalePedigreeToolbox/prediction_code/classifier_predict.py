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
    data_path = Path(sys._MEIPASS)  # noqa
elif __file__:
    data_path = Path(os.path.dirname(__file__))
else:
    raise SystemExit("Can not find application path.")


LOG = logging.getLogger("mpt")


# marker orders for the pre-made models
RMPLEX_MARKER_LIST = ['DYF403S1a', 'DYS724', 'DYS711', 'DYS1010', 'DYF1000', 'DYS1012', 'DYF399S1', 'DYS626', 'DYS712',
                      'DYS547', 'DYS518', 'DYS442', 'DYF1002', 'DYF387S1', 'DYS576', 'DYS1013', 'DYF404S1', 'DYF403S1b',
                      'DYF393S1', 'DYS526b', 'DYR88', 'DYF1001', 'DYS1007', 'DYS1005', 'DYS1003', 'DYS627', 'DYS612',
                      'DYS713', 'DYS449', 'DYS570']
PPY23_MARKER_LIST = ['DYS643', 'DYS549', 'DYS438', 'DYS456', 'DYS570', 'DYS481', 'DYS533', 'DYS448', 'DYS437', 'DYS458',
                     'DYS385', 'DYS576', 'DYS390', 'DYS392', 'YGATAH4', 'DYS389II', 'DYS439', 'DYS389I', 'DYS391',
                     'DYS393', 'DYS19', 'DYS635']
PPY23_RMPLEX_MARKER_LIST = ['DYS570', 'DYS439', 'DYS626', 'DYF403S1a', 'DYF404S1', 'DYF1000', 'DYS635', 'DYS389I',
                            'DYS724', 'DYF387S1', 'DYS643', 'DYS19', 'DYS1003', 'DYF393S1', 'DYS391', 'DYS393',
                            'DYS526b', 'YGATAH4', 'DYS385', 'DYF1001', 'DYS1010', 'DYS627', 'DYS711', 'DYS481',
                            'DYF1002', 'DYS547', 'DYS458', 'DYS1013', 'DYS1012', 'DYS549', 'DYS712', 'DYS392',
                            'DYS713', 'DYS448', 'DYS1005', 'DYR88', 'DYS442', 'DYS438', 'DYS518', 'DYS390', 'DYF399S1',
                            'DYS437', 'DYS449', 'DYS389II', 'DYS456', 'DYS576', 'DYS612', 'DYF403S1b', 'DYS1007',
                            'DYS533']
YFP_MARKER_LIST = ['DYS392', 'DYS448', 'DYS393', 'DYS627', 'DYS437', 'DYS439', 'DYS19', 'DYS456', 'DYS458', 'DYS391',
                   'DYS576', 'DYS481', 'DYS389II', 'DYS518', 'YGATAH4', 'DYS570', 'DYS635', 'DYS390', 'DYS385',
                   'DYS389I', 'DYS449', 'DYF387S1', 'DYS438', 'DYS533', 'DYS460']
YFP_RMPLEX_MARKER_LIST = ['DYS456', 'YGATAH4', 'DYS724', 'DYS19', 'DYS437', 'DYS393', 'DYS712', 'DYF403S1a',
                          'DYF404S1', 'DYS460', 'DYS481', 'DYF1002', 'DYF387S1', 'DYS1003', 'DYS612', 'DYS1005',
                          'DYF399S1', 'DYF403S1b', 'DYS448', 'DYS389II', 'DYS390', 'DYS442', 'DYF1000', 'DYS1012',
                          'DYS526b', 'DYS385', 'DYS1007', 'DYF1001', 'DYS439', 'DYS626', 'DYS438', 'DYS711', 'DYS533',
                          'DYS570', 'DYS627', 'DYF393S1', 'DYS391', 'DYS392', 'DYS713', 'DYS1013', 'DYS547', 'DYS389I',
                          'DYS518', 'DYS449', 'DYS635', 'DYS458', 'DYS576', 'DYR88', 'DYS1010']
YFORGEN_MARKER_LIST = ['DYS438', 'DYS389II', 'DYS385', 'DYS570', 'DYS455', 'DYS459', 'DYS533', 'DYS391', 'DYS607',
                       'DYS724', 'DYS518', 'DYS392', 'DYS643', 'DYS447', 'DYS460', 'DYS426', 'DYS627', 'DYS388',
                       'DYS464', 'DYF387S1', 'DYS439', 'DYS456', 'DYS390', 'DYS635', 'DYS19', 'DYS454', 'DYS442',
                       'YCAII', 'DYS576', 'DYS549', 'DYS389I', 'YGATAH4', 'DYS437', 'DYS448', 'DYS393', 'DYS449',
                       'DYS458', 'DYS481']
YFORGEN_RMPLEX_LIST = ['DYS439', 'DYS724', 'DYS455', 'DYS449', 'DYS518', 'DYF403S1b', 'DYS549', 'DYS390', 'DYS389II',
                       'DYS713', 'DYS1013', 'DYS437', 'DYS1007', 'YGATAH4', 'DYS447', 'DYS456', 'DYS612', 'DYS712',
                       'DYS1012', 'DYS533', 'DYS464', 'DYS643', 'DYS627', 'DYS547', 'DYS481', 'DYS458', 'DYS454',
                       'DYS576', 'DYF1002', 'DYS392', 'DYS607', 'DYS19', 'DYS459', 'DYS389I', 'DYS385', 'DYS1005',
                       'DYF399S1', 'DYS626', 'DYR88', 'DYF404S1', 'DYS711', 'DYF403S1a', 'DYS635', 'DYS388', 'DYS526b',
                       'DYS391', 'DYS393', 'DYF393S1', 'DYF387S1', 'DYS570', 'DYS426', 'DYF1000', 'YCAII', 'DYS438',
                       'DYS1003', 'DYS442', 'DYS448', 'DYS460', 'DYS1010', 'DYF1001']

MARKER_MAPPING = {"RMPLEX": RMPLEX_MARKER_LIST, "PPY23": PPY23_MARKER_LIST, "YFP": YFP_MARKER_LIST,
                  "PPY23_RMPLEX": PPY23_RMPLEX_MARKER_LIST, "YFP_RMPLEX": YFP_RMPLEX_MARKER_LIST}


# speeds up calculation of prediction ranges a lot, since a lot of the same results are present
PREDICTION_RANGE_CACHE = {}


def load_model_paths():
    models = {"RMPLEX": data_path / "models" / "RMplex_model.joblib",
              "PPY23": data_path / "models" / "PPY23_model.joblib",
              "PPY23_RMPLEX": data_path / "models" / "PPY23_RMplex_model.joblib",
              "YFP": data_path / "models" / "YFP_model.joblib",
              "YFP_RMPLEX": data_path / "models" / "YFP_RMplex_model.joblib",
              "YFORGEN": data_path / "models" / "YForgen_model.joblib",
              "YFORGEN_RMPLEX": data_path / "models" / "YForgen_RMplex_model.joblib"}
    # make sure to not change size during itteration
    for name, path in list(models.items()):
        if not path.exists():
            LOG.warning(f"Failed to find pre-defined model for {name}")
            del models[name]
    return models


MODELS = load_model_paths()


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
            header_line = f.readline().strip()
        return header_line.split(",")[1:-1]
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
        try:
            return joblib.load(MODELS[predefined_model_name])
        except KeyError:
            LOG.error("Failed to load pre-compiled model. The program can not find the pre-computed models.")
            raise utility.MalePedigreeToolboxError("Failed to load pre-compiled model. The program can not find"
                                                   " the pre-computed models.")


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
