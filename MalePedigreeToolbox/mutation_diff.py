#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Module for drawing dendograms based on mutation differentiation.

(If time allows it, it makes sense to rewrite the this module)

Author: Diego
Changed by: Bram
"""

import pandas as pd
import numpy as np
import logging
import math
from statsmodels.stats.proportion import proportion_confint
from typing import Union

# own imports
from MalePedigreeToolbox import utility
from MalePedigreeToolbox import thread_termination

LOG = logging.getLogger("mpt")

SUMMARY_OUT = "summary_out.csv"
FULL_OUT = "full_out.csv"
DIFFERENTIATION_OUT = "differentiation_out.csv"
PREDICT_OUT = "predict_out.csv"

SCORE_CACHE = {}  # used to store already computed scores to potentially speed up calculations


def get_score(value1, value2):
    return math.ceil(abs(value1 - value2))


def get_decimal(number: Union[int, float]) -> int:
    # make sure that the full number is returned and no strange float rounding occurs
    nr_frac = str(number).split(".")
    if len(nr_frac) == 1:
        return 0
    return int(nr_frac[1])


def get_matrix_scores(a, b):
    matrix_scores = np.ndarray((len(a), len(b)))
    match_decimal_matrix = np.ndarray((len(a), len(b)))

    for row in range(matrix_scores.shape[0]):
        for col in range(matrix_scores.shape[1]):
            matrix_scores[row][col] = get_score(a[row][1], b[col][1])
            match_decimal_matrix[row][col] = True if get_decimal(a[row][1]) == get_decimal(b[col][1]) else False
    return matrix_scores, match_decimal_matrix


def sort_allele(allele):
    allele = [(index, value) for index, value in enumerate(allele)]
    allele.sort(key=lambda x: (x[1], x[0]))
    return allele


@thread_termination.ThreadTerminable
def get_mutation_diff(alleles1, alleles2, expected_size):
    """Get the difference between 2 rows of numbers. Matching decimals is first priority then smallest difference.

    Decimals have priority because we assume that a mutation from 1.1 to 3.1 is more likely than 1.1 to 2.0 even though
    the score is lower"""
    # use frozen set to cache independant of order
    cache_key = frozenset([frozenset(alleles1), frozenset(alleles2), expected_size])
    if cache_key in SCORE_CACHE:
        return SCORE_CACHE[cache_key]

    # remove all 0 signals
    alleles1 = [value for value in alleles1 if value != 0]
    alleles2 = [value for value in alleles2 if value != 0]
    if len(alleles1) == 0:
        alleles1 = [0]
    if len(alleles2) == 0:
        alleles2 = [0]

    # sort allele values from smallest to largest to make sure that all assumptions are met, each allele value now has
    # an additional original index attached
    alleles1 = sort_allele(alleles1)
    alleles2 = sort_allele(alleles2)

    # ensure the biggest allele is always allele 1
    if len(alleles2) > len(alleles1):
        temp = alleles1
        alleles1 = alleles2
        alleles2 = temp
        scores = get_assymetrical_difference(alleles1, alleles2)
    elif len(alleles1) > len(alleles2):
        scores = get_assymetrical_difference(alleles1, alleles2)
    else:
        scores = get_symmetrical_difference(alleles1, alleles2)

    # revert back to original order
    scores.sort(key=lambda x: x[0])
    # remove all the index values
    scores = [score[1] for score in scores]
    if len(alleles1) < expected_size:
        # make sure that when both alleles are shorter than the real allele that the scores are still accurate
        missing_difference = expected_size - len(alleles1)
        min_score = min(scores)
        for _ in range(missing_difference):
            scores.append(min_score)
    SCORE_CACHE[cache_key] = scores
    return scores


@thread_termination.ThreadTerminable
def get_assymetrical_difference(alleles1, alleles2):
    # in case alleles1 > alleles2
    matrix_scores, match_decimal_matrix = get_matrix_scores(alleles1, alleles2)
    final_scores = []
    for row_index in range(len(alleles1)):
        any_decimal_match = any(match_decimal_matrix[row_index])
        lowest_score = 1_000_000
        lowest_orig_index = 0  # will be changed if used
        added_score = False
        for column_index in range(len(alleles2)):
            if any_decimal_match:
                if match_decimal_matrix[row_index][column_index] is False:
                    continue
                score = matrix_scores[row_index][column_index]
                if score == 0:
                    # add original index
                    final_scores.append((alleles1[row_index][0], score))
                    added_score = True
                    break
                if score < lowest_score:
                    lowest_score = score
                    lowest_orig_index = alleles1[row_index][0]
            else:
                matrix_row = [(index, value) for index, value in enumerate(matrix_scores[row_index])]
                final_scores.append(min(matrix_row, key=lambda x: x[1]))
                added_score = True
                break
        if added_score is False:
            final_scores.append([lowest_orig_index, lowest_score])
    return final_scores


@thread_termination.ThreadTerminable
def get_symmetrical_difference(alleles1, alleles2):

    # if the lenght of both alleles is equal
    matrix_scores, match_decimal_matrix = get_matrix_scores(alleles1, alleles2)
    final_scores = []
    alleles2_sorted_indexes = [index for index in range(len(alleles2))]
    for row_index in range(len(alleles1)):
        any_decimal_match = any(match_decimal_matrix[row_index])
        lowest_score = 1_000_000
        lowest_index = 0  # will be changed if used
        lowest_orig_index = 0  # will be changed if used
        added_score = False
        for index, column_index in enumerate(alleles2_sorted_indexes):
            if any_decimal_match:
                if match_decimal_matrix[row_index][column_index] is False:
                    continue
                score = matrix_scores[row_index][column_index]
                if score == 0:
                    final_scores.append((alleles1[row_index][0], score))
                    del alleles2_sorted_indexes[index]
                    added_score = True
                    break
                if score < lowest_score:
                    lowest_score = score
                    lowest_index = index
                    lowest_orig_index = alleles1[row_index][0]
            else:
                matrix_row = [(index, value) for index, value in enumerate(matrix_scores[row_index])]
                final_scores.append(min(matrix_row, key=lambda x: x[1]))
                added_score = True
                del alleles2_sorted_indexes[index]
                break
        if added_score is False:
            final_scores.append([lowest_orig_index, lowest_score])
            del alleles2_sorted_indexes[lowest_index]
    return final_scores


@thread_termination.ThreadTerminable
def write_differentiation_rates(mutation_dict_list, distance_dict, outfile):
    # the rate of differentiation given a certain distance of a 2 subjects in a pair
    meiosis_dict = {}
    covered_pairs = set()
    mutated_pairs = set()
    warned_pedigrees = set()  # make sure the log is less spammy
    for dictionary in mutation_dict_list:
        differentiated = dictionary["Total"] != 0
        pedigree = dictionary["Pedigree"]
        pair = dictionary["From"] + dictionary["To"]
        reverse_pair = dictionary["To"] + dictionary["From"]

        if pedigree not in distance_dict:
            if pedigree not in warned_pedigrees:
                LOG.warning(f"Can not include pedigree {pedigree} in differentiation rate calculation since they are"
                            f" not present in the distance file.")
            warned_pedigrees.add(pedigree)
            continue

        if pair in distance_dict[pedigree]:
            distance = distance_dict[pedigree][pair]
        elif reverse_pair in distance_dict[pedigree]:
            distance = distance_dict[pedigree][reverse_pair]
        else:
            LOG.warning(f"Can not include pair {dictionary['To']}-{dictionary['From']} in the differentiation rate "
                        f"calculation since they are not present in the distance file.")
            continue
        if distance in meiosis_dict:
            if pair not in covered_pairs:
                meiosis_dict[distance][0] += 1
        else:
            meiosis_dict[distance] = [1, 0]
        if pair not in mutated_pairs and differentiated:
            meiosis_dict[distance][1] += 1
            mutated_pairs.add(pair)
        covered_pairs.add(pair)
    meiosis_list = []
    for key, values in meiosis_dict.items():
        ci = [str(round(x * 100, 2)) for x in proportion_confint(values[1], values[0], method='beta')]
        meiosis_list.append((key, *values, round(values[1] / values[0] * 100, 2), *ci))
    meiosis_list.sort(key=lambda x: x[0])

    final_text = "Meioses,Pairs,Differentiated,Differentiation_rate(%),Clopper-Pearson CI lower bound, " \
                 "Clopper-Pearson CI upper bound\n"
    for values in meiosis_list:
        final_text += ",".join(map(str, values)) + "\n"

    with open(outfile, "w") as f:
        f.write(final_text)


@thread_termination.ThreadTerminable
def sample_combinations(samples):
    combinations = []
    for index, sample in enumerate(samples):
        for inner_index in range(index + 1, len(samples)):
            combinations.append([sample, samples[inner_index]])
    return combinations


@thread_termination.ThreadTerminable
def read_distance_file(distance_file):
    # read the distance file into a quickly accesible dictionary
    distance_dict = {}
    with open(distance_file) as f:
        f.readline()  # skip header
        for line in f:
            values = line.strip().split(",")
            pedigree_name = values[0]
            sample1 = values[1]
            sample2 = values[2]
            distance = int(values[3])
            pair = f"{sample1}{sample2}"
            if pedigree_name in distance_dict:
                distance_dict[pedigree_name][pair] = distance
            else:
                distance_dict[pedigree_name] = {pair: distance}
    return distance_dict


@thread_termination.ThreadTerminable
def main(name_space):
    LOG.info("Starting with calculating differentiation rates")

    alleles_file = name_space.allele_file
    distance_file = name_space.dist_file
    outdir = name_space.outdir
    include_predict_file = name_space.prediction_file
    if alleles_file.suffix == ".xlsx":
        alleles_df = pd.read_excel(alleles_file, dtype={'Pedigree': str, 'Sample': str, 'Marker': str,
                                                        'Allele_1': np.float64, 'Allele_2': np.float64,
                                                        'Allele_3': np.float64, 'Allele_4': np.float64,
                                                        'Allele_5': np.float64, 'Allele_6': np.float64})
    elif alleles_file.suffix == ".csv":
        alleles_df = pd.read_csv(alleles_file, dtype={'Pedigree': str, 'Sample': str, 'Marker': str,
                                                      'Allele_1': np.float64, 'Allele_2': np.float64,
                                                      'Allele_3': np.float64, 'Allele_4': np.float64,
                                                      'Allele_5': np.float64, 'Allele_6': np.float64})
    else:
        LOG.error(f"Unsupported file type .{alleles_file.suffix} for the alleles file.")
        raise utility.MalePedigreeToolboxError(f"Unsupported file type .{alleles_file.suffix}"
                                               f" for the alleles file.")
    run(alleles_df, distance_file, outdir, include_predict_file)


@thread_termination.ThreadTerminable
def run(alleles_df, distance_file, outdir, include_predict_file):

    alleles_list_dict = alleles_df.to_dict('records')

    if len(alleles_list_dict) == 0:
        LOG.error("Empty alleles file provided")
        raise utility.MalePedigreeToolboxError("Empty alleles file provided")
    total_alleles_specified = len(alleles_list_dict[0]) - 3

    # pre-sort pedigree information
    LOG.debug("Pre-sorting alleles information")
    grouped_alleles_dict = {}
    longest_allele_per_marker = {}
    for dictionary in alleles_list_dict:
        try:
            pedigree_name = dictionary.pop("Pedigree")
            sample_name = dictionary.pop("Sample")
            marker = dictionary.pop("Marker")
        except KeyError:
            LOG.error("Incorrect alleles file. The following three column names are required: 'Pedigree', 'Sample', "
                      "'Marker'.")
            raise utility.MalePedigreeToolboxError("Incorrect alleles file. The following three column names are "
                                                   "required: 'Pedigree', 'Sample', 'Marker'.")
        allele = [x for x in dictionary.values() if not math.isnan(x)]
        if marker not in longest_allele_per_marker:
            longest_allele_per_marker[marker] = len(allele)
        else:
            longest_allele_per_marker[marker] = max(longest_allele_per_marker[marker], len(allele))
        if pedigree_name in grouped_alleles_dict:
            if sample_name in grouped_alleles_dict[pedigree_name]:
                grouped_alleles_dict[pedigree_name][sample_name][marker] = allele
            else:
                grouped_alleles_dict[pedigree_name][sample_name] = {marker: allele}
        else:
            grouped_alleles_dict[pedigree_name] = {sample_name: {marker: allele}}
    LOG.info("Finished reading both input files")
    markers = set(alleles_df.Marker)

    LOG.info(f"In total there are {len(markers)} markers being analysed.")
    mutation_dict = []
    total_mutation_dict = []
    predict_pedigrees_list = []

    prev_total = 0
    for index, (pedigree, pedigree_data) in enumerate(grouped_alleles_dict.items()):
        sample_names = list(pedigree_data.keys())
        sample_combs = sample_combinations(sample_names)
        predict_samples_list = []
        LOG.info(f"Comparing {len(sample_combs)} allele combinations for pedigree {pedigree}")
        for sample1, sample2 in sample_combs:
            sample1_data = pedigree_data[sample1]
            sample2_data = pedigree_data[sample2]
            total_mutations = 0
            marker_values = {name: 0 for name in markers}
            for marker in markers:
                count_mutation = 0

                if marker not in sample1_data or marker not in sample2_data:
                    LOG.warning(f"Marker ({marker}) is not present in {sample1} and {sample2}. The comparisson will be"
                                f" skipped.")
                    continue
                marker_data1 = sample1_data[marker]
                marker_data2 = sample2_data[marker]
                mutations = get_mutation_diff(marker_data1, marker_data2, longest_allele_per_marker[marker])
                [mutations.append(0.0) for _ in range(total_alleles_specified - len(mutations))]

                count_mutation += np.sum(mutations)
                total_mutations += count_mutation
                if marker in marker_values and include_predict_file:
                    marker_values[marker] = count_mutation
                mutation_dict.append({"Marker": marker, "Pedigree": pedigree, "From": sample1, "To": sample2,
                                      **{f"Allele_{index + 1}": mutations[index] for index in range(total_alleles_specified)},
                                      "Total": count_mutation})

            # for predicting the generational distances
            if include_predict_file:
                predict_samples_list.append([f"{pedigree}_{sample1}_{sample2}",
                                             *[marker_values[name] for name in markers]])

            total_mutation_dict.append({"Pedigree": pedigree, "From": sample1, "To": sample2, "Total": total_mutations})

            LOG.debug(f"Finished calculating differentiation for {index} out of {len(grouped_alleles_dict)}")
            total, remainder = divmod(index / len(grouped_alleles_dict), 0.05)

            if total != prev_total:
                LOG.info(f"Calculation progress: {round(5 * total)}%...")
                prev_total = total
        if include_predict_file:
            predict_pedigrees_list.append(predict_samples_list)

    mutation_df = pd.DataFrame(mutation_dict)
    total_mutation_df = pd.DataFrame(total_mutation_dict)

    mutation_df_cols = ['Pedigree', 'From', 'To', 'Marker',
                        *[f"Allele_{index + 1}" for index in range(total_alleles_specified)], 'Total']
    total_mutations_cols = ['Pedigree', 'From', 'To', 'Total']

    LOG.info("Starting with writing mutation differentiation information to files")
    mutation_df = mutation_df[mutation_df_cols]
    total_mutation_df = total_mutation_df[total_mutations_cols]
    mutation_df.to_csv(outdir / FULL_OUT)

    total_mutation_df.to_csv(outdir / SUMMARY_OUT)

    if distance_file is not None:
        # read the distance file
        LOG.info("Started with summarising and writing meiosis differentiation rates to file")
        distance_dict = read_distance_file(distance_file)
        write_differentiation_rates(mutation_dict, distance_dict, outdir / DIFFERENTIATION_OUT)

    if include_predict_file:
        predict_text_list = [f"sample,{','.join(markers)}"]
        for pedigree_list in predict_pedigrees_list:
            for sample_list in pedigree_list:
                predict_text_list.append(','.join(list(map(str, sample_list))))
        with open(outdir / PREDICT_OUT, "w") as f:
            f.write('\n'.join(predict_text_list))

    LOG.info("Finished calculating differentiation rates.")


if __name__ == '__main__':
    # all these should be true

    func = get_mutation_diff
    l1 = [1.0, 0.0, 0.0, 0.0, 0.0]
    l2 = [4.0, 0.0, 0.0, 0.0, 0.0]
    assert func(l1, l2, 1) == [3.0]
    l1 = [12, 13, 18]
    l2 = [12, 18]
    assert func(l1, l2, 3) == [0.0, 1.0, 0.0]
    l1 = [12, 13, 14]
    l2 = [12]
    assert func(l1, l2, 3) == [0.0, 1.0, 2.0]
    l1 = [12, 13, 14]
    l2 = [12, 14, 18]
    assert func(l1, l2, 3) == [0.0, 1.0, 4.0]
    l1 = [12, 13, 18, 19]
    l2 = [12, 18]
    assert func(l1, l2, 3) == [0.0, 1.0, 0.0, 1.0]
    l1 = [12, 13, 16]
    l2 = [13, 12, 16]
    assert func(l1, l2, 3) == [0.0, 0.0, 0.0]
    l1 = [13, 12.1, 16]
    l2 = [12.1, 13, 16]
    assert func(l1, l2, 3) == [0.0, 0.0, 0.0]
    l1 = [12.1, 13.1, 16]
    l2 = [13, 12.1, 16]
    assert func(l1, l2, 3) == [0.0, 1.0, 0.0]
    l1 = [12.1, 14, 16]
    l2 = [13, 12.1, 16]
    assert func(l1, l2, 3) == [0.0, 1.0, 0.0]
    l0 = [12, 13]
    l1 = [12]
    l2 = [13]
    assert func(l1, l2, 2) == [1.0, 1.0]
    l0 = [12, 13, 15] # noqa
    l1 = [12]
    l2 = [13]
    assert func(l1, l2, 3) == [1.0, 1.0, 1.0]
    l0 = [12, 13, 15]  # noqa
    l1 = [12, 15]
    l2 = [13, 15]
    assert func(l1, l2, 3) == [1.0, 0.0, 0.0]
    l1 = [12.1, 13.0]
    l2 = [12.0, 13.1]
    assert func(l1, l2, 2) == [1.0, 1.0]
    l1 = [13]
    l2 = [12.1, 13]
    assert func(l1, l2, 2) == [1.0, 0.0]
    l1 = [12.1, 11]
    l2 = [11.1, 12.1, 11, 12]
    assert func(l1, l2, 4) == [1.0, 0.0, 0.0, 1.0]
    l1 = [16.2, 19.2, 0.0, 0.0]
    l2 = [16.2, 18.2, 19.2, 0.0]
    assert func(l1, l2, 3) == [0.0, 1.0, 0.0]
    l1 = [0.0]
    l2 = [0.0]
    assert func(l1, l2, 1) == [0.0]
    l1 = [10.0]
    l2 = [0.0]
    assert func(l1, l2, 1) == [10.0]
