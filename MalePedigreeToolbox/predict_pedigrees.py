#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Module for drawing dendograms based on mutation differentiation.

Author: Bram van Wersch

NOTE:
IMPORTANT: Python3.7 or higher needs to be used because this script relies on sorted dictionaries
"""

# library imports
from pathlib import Path
import numpy as np
import logging
from math import log10
import os
import warnings
from typing import List, Tuple, TYPE_CHECKING, Union, Dict
from collections import defaultdict

from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score  # calinski_harabasz_score, davies_bouldin_score other options

import matplotlib
import matplotlib.pyplot as plt

# own imports
from MalePedigreeToolbox import utility
from MalePedigreeToolbox import thread_termination


if TYPE_CHECKING:
    import argparse


warnings.filterwarnings(action='ignore', category=UserWarning)

matplotlib.use('Agg')

LOG: logging.Logger = logging.getLogger("mpt")

LINKAGE_METHOD = "complete"


@thread_termination.ThreadTerminable
def main(
    name_space: "argparse.Namespace"
):
    """Main entry point for this module"""
    LOG.info("Starting with creating dendograms based on mutation differentiation")
    input_file = name_space.full_marker_csv
    marker_rates_file = name_space.mutation_rates
    outdir = Path(name_space.outdir)
    min_distance = name_space.min_dist
    expected_nr_of_clusters = name_space.clusters

    with open(input_file) as f:
        lines = f.readlines()

    lines = adjust_fo_file_values(lines, marker_rates_file)

    covered_lines = 0
    prev_total = 0
    for index, pedigree_lines in enumerate(get_pedigree_lines(lines)):
        pedigree_name = pedigree_lines[0][0]

        # decide number of clusters as defined by user
        dist_mat, colnames = get_distance_matrix(pedigree_lines, min_distance)
        pedigree_dir = create_pedigree_dir(outdir, pedigree_name)

        # write the dendrogram
        LOG.debug(f"Writing dendogram for {pedigree_name}")
        draw_dendrogram(dist_mat, colnames, pedigree_dir, pedigree_name, expected_nr_of_clusters)
        LOG.debug(f"Finished making dendrogram plot for pedigree {pedigree_name}")

        covered_lines += len(pedigree_lines)

        # update user periodically
        total, remainder = divmod(covered_lines / len(lines), 0.05)
        if total != prev_total:
            LOG.info(f"Calculation progress: {round(5 * total)}%...")
            prev_total = total

    LOG.info("Finished drawing dendograms for all pedigrees that were present")


@thread_termination.ThreadTerminable
def get_requested_cluster_list(
    requested_clusters: Union[str, None]
) -> Union[None, List[int]]:
    """Get the number of clusters requested for all the different pedigrees that need to be drawn"""
    if requested_clusters is None:
        return None
    requested_cluser_numbers = []
    for value in requested_clusters:
        if isinstance(value, int):
            requested_cluser_numbers.append(value)
        else:  # in case of file
            with open(value) as f:
                try:
                    requested_cluser_numbers.extend(map(int, f.read().split()))
                except ValueError:
                    LOG.error(f"Invalid integer provided for the number of clusters in file {value}. Make sure that"
                              f"the cluster numbers are seperated by spaces.")
                    raise utility.MalePedigreeToolboxError(f"Invalid integer provided for the number of clusters in"
                                                           f" file {value}. Make sure that the cluster numbers are"
                                                           f" seperated by spaces.")
    return requested_cluser_numbers


@thread_termination.ThreadTerminable
def adjust_fo_file_values(
    lines: List[str],
    marker_rates_file: Union[Path, None]
) -> List[List[Union[str, float]]]:
    """Wheigh the line values based on mutation rates to better represent the importance of all the values.

    If no marker rates file is provided the default weight is 1 preserving the original values
    """
    if marker_rates_file is not None:
        marker_rate_dict = read_marker_rate_file(marker_rates_file)
    else:
        marker_rate_dict = defaultdict(lambda: 1)

    all_summary_lines = []
    last_values = lines[1].strip().split(",")
    summary_dict = defaultdict(int)
    for line in lines[1:]:
        values = line.strip().split(",")
        _, pedigree, first_id, second_id, marker = values[:5]
        mutations = int(values[-1])
        # a new dendrogram begins
        if pedigree != last_values[1]:
            for (id1, id2), mutation_wheight in summary_dict.items():
                all_summary_lines.append([last_values[1], id1, id2, mutation_wheight])
            summary_dict = defaultdict(int)
        try:
            summary_dict[(first_id, second_id)] += mutations * get_weight(marker_rate_dict[values[4]])
        except KeyError:
            LOG.error(f"Missing rate for marker {values[4]} in marker rates file.")
            raise utility.MalePedigreeToolboxError(f"Missing rate for marker {values[4]} in marker rates file.")
        last_values = values
    # make sure to do the last one
    for (id1, id2), mutation_wheight in summary_dict.items():
        all_summary_lines.append([last_values[1], id1, id2, mutation_wheight])
    return all_summary_lines


def get_weight(
    value: float
) -> float:
    """Weight calculation function"""
    if value == 1 or value == 0:
        return 1
    return - log10(value)


@thread_termination.ThreadTerminable
def read_marker_rate_file(
    file: Path
) -> Dict[str, float]:
    """Read a marker rate file in order to weight the number of mutations of all alleles"""
    rate_dict = {}
    with open(file) as f:
        f.readline()  # ignore the header
        for line in f:
            try:
                marker, rate = line.strip().split(",")
                rate = float(rate)
            except Exception:
                LOG.error("Failed to read marker rates file. The format of the file is not as expected. The file needs"
                          " a header and 2 columns. The first one being the marker names and the second one being the "
                          "marker rates.")
                raise utility.MalePedigreeToolboxError("Failed to reda marker rates file. The format of the file is "
                                                       "not as expected. The file needs a header and 2 columns. The "
                                                       "first one being the marker names and the second one being the "
                                                       "marker rates.")
            rate_dict[marker] = rate
    return rate_dict


@thread_termination.ThreadTerminable
def create_pedigree_dir(
    outdir: Path,
    pedigree_name: str
) -> Path:
    """Create the directory containing the files generated by this module"""
    pedigree_name = utility.sanitize_string(pedigree_name)
    pedigree_dir = outdir / pedigree_name
    # make directories for each
    try:
        os.mkdir(pedigree_dir)
    except FileExistsError:
        pass
    except IOError as e:
        LOG.error(f"Failed to create subfolder in output folder; {outdir}. Error code: {str(e)}")
        raise utility.MalePedigreeToolboxError(f"Failed to create subfolder in output folder; {outdir}."
                                               f" Error code: {str(e)}")
    return pedigree_dir


def get_pedigree_lines(
    lists: List[List[Union[str, float]]]
) -> List[List[List[str]]]:
    """Splits the input file based on pedigrees and returns connections of nodes and respective mutation
    differentiation"""
    dendogram_lines = [lists[0]]
    for list_ in lists[1:]:
        # a new dendogram begins
        if list_[0] != dendogram_lines[-1][0]:
            yield dendogram_lines
            dendogram_lines = [list_]
        else:
            dendogram_lines.append(list_)
    yield dendogram_lines


@thread_termination.ThreadTerminable
def get_distance_matrix(
    pedigree_lines: List[List[str]],
    min_dist: float
) -> Tuple[np.ndarray, List[str]]:
    """Get a distance matrix based connections of nodes in a pedigree"""
    pedigree_dict = {}
    for values in pedigree_lines:
        node1 = values[1]
        node2 = values[2]
        mutation_distance = float(values[3])
        if node1 not in pedigree_dict:
            pedigree_dict[node1] = {node2: mutation_distance}
        else:
            pedigree_dict[node1][node2] = mutation_distance

        if node2 not in pedigree_dict:
            pedigree_dict[node2] = {node1: mutation_distance}
        else:
            pedigree_dict[node2][node1] = mutation_distance
    dist_mat = np.zeros((len(pedigree_dict), len(pedigree_dict)), float)
    for r_index, key1 in enumerate(pedigree_dict):
        for c_index, key2 in enumerate(pedigree_dict):
            if key1 == key2:
                continue
            else:
                mutation_distance = pedigree_dict[key1][key2]
                if mutation_distance == 0:
                    mutation_distance = min_dist
                dist_mat[r_index, c_index] = mutation_distance
    return dist_mat, list(pedigree_dict.keys())


@thread_termination.ThreadTerminable
def draw_dendrogram(
    dist_mat: np.ndarray,
    names: List[str],
    outdir: Path,
    pedigree_name: str,
    expected_nr_of_clusters: str
):
    """Draws the dendograms"""
    fig = plt.figure(num=1, clear=True)
    dists = squareform(dist_mat)
    linkage_matrix = linkage(dists, LINKAGE_METHOD)
    if expected_nr_of_clusters == "opt":
        best_nr_clusters = max(find_optimal_clustering(dist_mat)) + 1
        ct = linkage_matrix[-(best_nr_clusters-1), 2]
    else:
        try:
            expected_nr = int(expected_nr_of_clusters)
        except ValueError:
            LOG.error(f"Invalid integer provided for the number of clusters: '{expected_nr_of_clusters}'.")
            raise utility.MalePedigreeToolboxError(f"Invalid integer provided for the number of clusters:"
                                                   f" '{expected_nr_of_clusters}'.")
        if expected_nr > len(dist_mat):
            raise utility.MalePedigreeToolboxError(f"The number of clusters can not be bigger than the number of"
                                                   f" datapoints. Got {expected_nr} clusters for pedigree"
                                                   f" {pedigree_name} with {len(dist_mat)} datapoints.")
        if expected_nr == 0:
            ct = 0
        else:
            ct = linkage_matrix[-(expected_nr-1), 2]
    if ct > 0:
        dendogram_dct = dendrogram(linkage_matrix, labels=names, leaf_rotation=90, color_threshold=ct)  # noqa
    else:
        dendogram_dct = dendrogram(linkage_matrix, labels=names, leaf_rotation=90, color_threshold=0,  # noqa
                                   above_threshold_color='k')

    plt.title(f"Pedigree {pedigree_name}")
    out_path = outdir / f"{pedigree_name}_predicted_dendogram.png"
    plt.tight_layout()
    fig.savefig(out_path)

    # write the dendogram clusters
    leave_order = dendogram_dct["ivl"]
    # old version of scipy does not have this feature
    if 'leaves_color_list' in dendogram_dct:
        leave_colors = dendogram_dct['leaves_color_list']
    else:
        # this does not always work exactly
        leave_colors = get_cluster_colors(dendogram_dct)
    current_color = leave_colors[0]
    cluster_nr = 1
    all_cluster_info = []
    cluster_info = f"cluster {cluster_nr}:\n"

    for index, color in enumerate(leave_colors):
        if current_color == color:
            cluster_info += f"{leave_order[index]}\n"
        else:
            current_color = color
            all_cluster_info.append(cluster_info + "\n")
            cluster_nr += 1
            cluster_info = f"cluster {cluster_nr}:\n{leave_order[index]}\n"
    all_cluster_info.append(cluster_info)
    with open(outdir / f"{pedigree_name}_dendogram_clusters.txt", "w") as f:
        f.write(''.join(all_cluster_info))


def get_cluster_colors(dendrogram_dct):
    # solution that works for python 3.6, can be slightly wrong unfortunately
    # https://stackoverflow.com/questions/61959602/retrieve-leave-colors-from-scipy-dendrogram/61964297#61964297
    points = dendrogram_dct['leaves']
    colors = ['none'] * len(points)
    for xs, c in zip(dendrogram_dct['icoord'], dendrogram_dct['color_list']):
        for xi in xs:
            if xi % 10 == 5:
                colors[(int(xi) - 5) // 10] = c
    return colors


@thread_termination.ThreadTerminable
def find_optimal_clustering(
    points: List[float]
) -> List[int]:
    """Calculate the optimal clustering for a set of datapoints for the cluster with the highest silhoute score.

    Return a cluster number for all the provided labels."""
    if len(points) == 2:
        return [0, 1]

    hierarchical = AgglomerativeClustering(n_clusters=2, linkage=LINKAGE_METHOD).fit(points)
    best_score = silhouette_score(points, hierarchical.labels_) # noqa
    best_labels = hierarchical.labels_  # noqa
    for k in range(3, len(points)):
        hierarchical = AgglomerativeClustering(n_clusters=k, linkage=LINKAGE_METHOD).fit(points)
        score = silhouette_score(points, hierarchical.labels_) # noqa
        if score > best_score:
            best_score = score
            best_labels = hierarchical.labels_  # noqa
    return best_labels
