#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Make a number of different models in a simple automated way. Keep in mind that to large training sets or a large number
of more complex models will take an extended amount of time. It is advised to run these programs with a large number
of processes

@author: Bram and Diego
"""

from typing import TYPE_CHECKING, List, Tuple, Dict, Any
import time
import logging
import numpy as np
import random
import warnings

import pandas as pd
import matplotlib.pyplot as plt
from joblib import dump

from sklearn.svm import SVC
from sklearn.model_selection import RandomizedSearchCV
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis, LinearDiscriminantAnalysis
from sklearn.metrics import f1_score, accuracy_score, precision_score, recall_score, \
    matthews_corrcoef, confusion_matrix
from sklearn.metrics import ConfusionMatrixDisplay
from sklearn.model_selection import StratifiedKFold
from sklearn.exceptions import DataConversionWarning, UndefinedMetricWarning, FitFailedWarning

# own imports
from MalePedigreeToolbox import utility


if TYPE_CHECKING:
    from pathlib import Path
    import argparse


# since sklearn can go a little wild with warnings lets filter some
warnings.filterwarnings(action='ignore', category=DataConversionWarning)
warnings.filterwarnings(action='ignore', category=UndefinedMetricWarning)
warnings.filterwarnings(action='ignore', category=FitFailedWarning)
warnings.filterwarnings(action='ignore', category=UserWarning)
warnings.filterwarnings(action='ignore', category=RuntimeWarning)

LOG = logging.getLogger("mpt")


def main(name_space: "argparse.Namespace"):
    LOG.info("Starting with making models")
    df = pd.read_csv(name_space.input)
    chosen_models = name_space.model_types
    cv_splits = name_space.cv_splits
    hyper_paramater_draws = name_space.hyper_parameter_choices
    param_choice_cv = name_space.parameter_cv
    nr_jobs = name_space.cpus
    out_folder = name_space.outdir

    learn_models(df, chosen_models, cv_splits, hyper_paramater_draws, param_choice_cv, nr_jobs, out_folder)
    LOG.info("Finished making all models")


def learn_models(
    df: pd.DataFrame,
    chosen_models: List[str],
    outer_splits: int,
    total_parameter_choices: int,
    parameter_choice_cv: int,
    n_jobs: int,
    outfolder: "Path"
):
    try:
        response_variable = df["Dist"]
        response_variable = np.array(list(response_variable))
        df.drop(columns=["Dist", "sample"], inplace=True)
    except AttributeError:
        LOG.error("Missing column names. Expected a 'Dist' and 'sample' column.")
        raise utility.MalePedigreeToolboxError("Missing column names. Expected a 'Dist' and 'sample' column.")
    learning_data = df.values

    if outer_splits == 1:
        indexes = list(range(len(df.index)))
        random.shuffle(indexes)
        splits = [[indexes[:int(0.8 * len(df.index))], indexes[int(0.8 * len(df.index)):]]]
    else:
        kfold = StratifiedKFold(n_splits=outer_splits, shuffle=True)
        splits = kfold.split(learning_data, response_variable)

    LOG.info("Created output_scores file in which all output scores will be written")
    with open(outfolder / "output_scores.tsv", "w") as f:
        f.write("name\taccuracy\tprecission\trecall\tf1\tmcc\ttime\tbest params\n")

    for model_name in chosen_models:
        model, param_grid = ESTIMATORS[model_name]
        model_name = model.__name__
        LOG.info(f"Starting with predictions for model {model_name}")
        for fold, (index_train, index_test) in enumerate(splits):
            random.shuffle(index_train)
            start_time = time.time()
            LOG.info(f"Outer cross validation run #{fold + 1}")
            x_train = learning_data[index_train]
            y_train = response_variable[index_train]
            x_test = learning_data[index_test]
            y_test = response_variable[index_test]
            if parameter_choice_cv == 1:
                indexes = list(range(len(y_train)))
                random.shuffle(indexes)
                inner_cv = [[indexes[:int(0.8 * len(y_train))], indexes[int(0.8 * len(y_train)):]]]
            else:
                inner_cv = StratifiedKFold(n_splits=parameter_choice_cv, shuffle=True)
            try:
                best_algo = RandomizedSearchCV(estimator=model(probability=True),
                                               param_distributions=param_grid, cv=inner_cv,
                                               verbose=1, n_jobs=n_jobs, n_iter=total_parameter_choices)
            except TypeError:
                best_algo = RandomizedSearchCV(estimator=model(),
                                               param_distributions=param_grid, cv=inner_cv,
                                               verbose=1, n_jobs=n_jobs, n_iter=total_parameter_choices)
            best_algo.fit(x_train, y_train)
            y_pred = best_algo.best_estimator_.predict(x_test)
            f1, acc, prec, rec, mcc = get_testing_results(y_pred, y_test)
            try:
                draw_confusion_matrix(y_pred, y_test, outfolder / f"{model_name}_confusion_matrix_{fold + 1}.png")
            except Exception as e:
                LOG.warning(f"Failed to draw confusion matrix with message {str(e)}")

            end_time = time.time() - start_time

            dump(best_algo, outfolder / f"{model_name}_model_{fold + 1}.joblib")
            with open(outfolder / "output_scores.tsv", "a") as f:
                f.write(f"{model_name}\t{acc}\t{prec}\t{rec}\t{f1}\t{mcc}\t{end_time}"
                        f"\t{str(best_algo.best_params_)}\n")


def get_testing_results(
    y_pred: List[float],
    y_test: List[float]
) -> Tuple[float, float, float, float, float]:
    f1 = f1_score(y_test, y_pred, average="weighted", labels=list(set(y_pred)))
    acc = accuracy_score(y_test, y_pred)
    prec = precision_score(y_test, y_pred, average="weighted", labels=list(set(y_pred)))
    rec = recall_score(y_test, y_pred, average="weighted", labels=list(set(y_pred)))
    mcc = matthews_corrcoef(y_test, y_pred)
    return f1, acc, prec, rec, mcc


def draw_confusion_matrix(
    y_pred: List[float],
    y_test: List[float],
    outfile: "Path"
):
    matrix = confusion_matrix(y_test, y_pred)
    matrix = matrix.astype('float') / matrix.sum(axis=1)[:, np.newaxis]
    matrix = matrix.round(decimals=1)

    plt.figure()
    plt.clf()
    cm = ConfusionMatrixDisplay(matrix)

    fig, ax = plt.subplots(figsize=(15, 15))
    cm.plot(ax=ax)
    plt.tight_layout()
    plt.savefig(outfile)


def knn_params() -> Dict[str, Any]:
    return {
        'n_neighbors': list(range(5, 101)),
        'weights': ["uniform", "distance"],
        'algorithm': ["auto", "kd_tree", "ball_tree", "brute"],
        'leaf_size': list(range(30, 61)),
        'p': [1, 2]
    }


def logistic_params() -> Dict[str, Any]:
    return {
        'penalty': ['l1'],
        'C': np.logspace(-4, 4, 20),
        'max_iter': [1000],
        'solver': ['liblinear']
    }


def rf_params() -> Dict[str, Any]:
    return {
        'n_estimators': range(100, 1001, 50),
        'max_features': ["auto", "sqrt"],
        'min_samples_split': range(2, 11),
        'min_samples_leaf': range(1, 5),
        'bootstrap': [True, False],
        'max_depth': [10]
    }


def mlp_params2() -> Dict[str, Any]:
    return {
        "hidden_layer_sizes": [(nr1,) for nr1 in range(50, 100)],
        "learning_rate": ['constant', 'invscaling', 'adaptive'],
        "solver": ["sgd", "adam"],
        "activation": ["identity", "logistic", "tanh", "relu"],
        "beta_1": np.linspace(0.5, 0.99, 50),
        "beta_2": np.linspace(0.99, 0.9999, 1000),
        "alpha": np.linspace(0.00001, 0.1, 1000),
        "max_iter": [1000]
    }


def lm_params() -> Dict[str, Any]:
    return {
        "solver": ["svd", "lsqr", "eigen"],
        "tol": np.linspace(0.0001, 0.1, 1000)
    }


def guassian_params():
    return {
        "var_smoothing": np.linspace(0.000000001, 0.0001, 1000)
    }


def svm_params():
    return {
        'C': np.linspace(0.0001, 10, 100),
        'gamma': np.linspace(1, 10, 100),
        'kernel': ["rbf"]
    }


def quadratic_params():
    return {
        "tol": np.linspace(0.0001, 0.1, 1000)
    }


ESTIMATORS = {"KNN": (KNeighborsClassifier, knn_params()),
              "LDA": (LinearDiscriminantAnalysis, lm_params()),
              "logistic": (LogisticRegression, logistic_params()),
              "QDA": (QuadraticDiscriminantAnalysis, quadratic_params()),
              "RF": (RandomForestClassifier, rf_params()),
              "Gaussian": (GaussianNB, guassian_params()),
              "MLP": (MLPClassifier, mlp_params2()),
              "SVM": (SVC, svm_params())}
