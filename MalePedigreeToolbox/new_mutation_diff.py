from typing import Union, Tuple, List, Dict, FrozenSet, Any, TYPE_CHECKING, Set
from math import isclose, ceil
from collections import defaultdict
import pandas as pd
import numpy as np
import logging
import math
from statsmodels.stats.proportion import proportion_confint

# own imports
from MalePedigreeToolbox import utility
from MalePedigreeToolbox import thread_termination

if TYPE_CHECKING:
    from pathlib import Path


LOG: logging.Logger = logging.getLogger("mpt")
SUMMARY_OUT: str = "summary_out.csv"
FULL_OUT: str = "full_out.csv"
DIFFERENTIATION_OUT: str = "differentiation_out.csv"
PREDICT_OUT: str = "predict_out.csv"

SCORE_CACHE: Dict[FrozenSet, List[int]] = {}  # used to store already computed scores, speed up


class Allele:

    def __init__(self, components: List[float]):
        self._components = [comp for comp in components if comp != 0]
        if len(self._components) == 0:
            self._components = [0]
        self._decimals = [self._get_decimal(comp) for comp in self._components]

    def _get_decimal(
        self,
        number: Union[int, float]
    ) -> int:
        # make sure that the full number is returned and no strange float rounding occurs
        nr_frac = str(number).split(".")
        if len(nr_frac) == 1:
            return 0
        return int(nr_frac[1])

    def duplicate_component(self, index):
        new_value = self._components[index]
        self._components.append(new_value)
        new_decimal = self._decimals[index]
        self._decimals.append(new_decimal)

    def get_equalizable_decimals(self, other_allele: "Allele"):
        equalizables = set()
        for decimal in other_allele._decimals:
            if decimal in self._decimals:
                equalizables.add(decimal)
        if len(equalizables) == 1:
            return set()
        return equalizables

    def get_decimal_difference(self, other_allele: "Allele", equalizables: Set[int]):
        # compare if there are the same amount of decimals for components in this and the other allele
        decimal_count_dict = defaultdict(int)
        for decimal in self._decimals:
            if decimal not in equalizables:
                continue
            decimal_count_dict[decimal] += 1
        for decimal in other_allele._decimals:
            if decimal not in equalizables:
                continue
            decimal_count_dict[decimal] -= 1
            if decimal_count_dict[decimal] == 0:
                del decimal_count_dict[decimal]
        return decimal_count_dict

    def get_indexes_with_decimal(self, wanted_decimal):
        indexes = []
        for index, decimal in enumerate(self._decimals):
            if decimal == wanted_decimal:
                indexes.append(index)
        return indexes

    def __iter__(self):
        return iter(self._components)

    def __getitem__(self, item):
        return self._components[item]

    def __str__(self):
        return str(self._components)

    def __len__(self):
        return len(self._components)


class DifferenceMatrix:
    NO_MATCHING_DECIMAL_PENALTY: int = 1000

    allele1: Allele
    allele2: Allele

    def __init__(self, allele1, allele2, expected_size):
        if len(allele1) == 0 or len(allele2) == 0:
            raise ValueError("Cannot compute distance against empty allele")

        if len(allele2) > len(allele1):
            temp = allele1
            allele1 = allele2
            allele2 = temp
        self.allele1 = allele1
        self.allele2 = allele2

        self._rows = []
        self._create_matrix(expected_size)

    def _create_matrix(self, expected_size):
        # first fill the matrix based on the current alleles
        for allele_component1 in self.allele1:
            matrix_row = []
            for allele_component2 in self.allele2:
                difference = abs(allele_component1 - allele_component2)

                # make sure that float rounding does not interfere with proper scoring
                if not isclose(difference, round(difference)):
                    difference += self.NO_MATCHING_DECIMAL_PENALTY
                difference = ceil(difference)
                matrix_row.append(difference)
            self._rows.append(matrix_row)

        self._equalize_decimals()

        # in case both rows and columns need duplications we can already add some
        if self.nr_rows < expected_size:
            self._duplicate_simultanious(expected_size)

    def _equalize_decimals(self):
        equalizables = self.allele1.get_equalizable_decimals(self.allele2)
        decimal_difference_dict = self.allele1.get_decimal_difference(self.allele2, equalizables)
        # in case decimals are not matching assume that there are duplications
        if len(decimal_difference_dict) == 0:
            return
        for decimal, difference in decimal_difference_dict.items():
            if difference < 0:
                indexes = self.allele1.get_indexes_with_decimal(decimal)
                min_row_index = self._min_row(indexes)
                self.allele1.duplicate_component(min_row_index)
                self._rows.append(self._rows[min_row_index].copy())
            else:
                indexes = self.allele2.get_indexes_with_decimal(decimal)
                min_col_index = self._min_column(indexes)
                self.allele2.duplicate_component(min_col_index)

                for row in self._rows:
                    row.append(row[min_col_index])

    def _min_row(self, indexes):
        # row index with the lowest sum of values for a given number of indexes
        min_sum = sum(self._rows[indexes[0]])
        min_index = indexes[0]
        for index in indexes[1:]:
            row_sum = sum(self._rows[index])
            if row_sum > min_sum:
                continue
            min_sum = row_sum
            min_index = index
        return min_index

    def _min_column(self, indexes):
        # column index with the lowest sum of values
        min_sum = sum([row[indexes[0]] for index, row in enumerate(self._rows)])
        min_index = indexes[0]
        for index in indexes[1:]:
            col_values = [row[index] for index, row in enumerate(self._rows)]
            col_sum = sum(col_values)
            if col_sum > min_sum:
                continue
            min_sum = col_sum
            min_index = index
        return min_index

    def _duplicate_simultanious(self, expected_size):
        for _ in range(expected_size - self.nr_rows):
            rindex, cindex = self._get_minimum_coordinate()
            self._rows.append(self._rows[rindex].copy())
            self.allele1.duplicate_component(rindex)
            self.allele2.duplicate_component(cindex)
            for row in self._rows:
                row.append(row[cindex])

    def _get_minimum_coordinate(self):
        # get the coordinate with the minum value in the matrix
        lowest_value = self._rows[0][0]
        lowest_coordinate = (0, 0)
        for rindex, row in enumerate(self._rows):
            for cindex, value in enumerate(row):
                if value < lowest_value:
                    lowest_value = value
                    lowest_coordinate = (rindex, cindex)
        return lowest_coordinate

    def _get_sorted_score_coordinates(self):
        coordinates = [(i, j) for j in range(self.nr_columns) for i in range(self.nr_rows)]
        return sorted(coordinates, key=lambda coord: self._rows[coord[0]][coord[1]])

    def get_optimal_score(self):
        scores = []
        score_rows = []
        sorted_score_coordinates = self._get_sorted_score_coordinates()

        # select the optimal scores
        covered_rows = set()
        covered_columns = set()
        score_index = 0
        while len(scores) < len(self._rows[0]):
            row, column = sorted_score_coordinates[score_index]
            score_index += 1
            if row in covered_rows or column in covered_columns:
                continue
            covered_rows.add(row)
            covered_columns.add(column)
            scores.append(self._rows[row][column])
            score_rows.append(row)

        # if there are more rows add allele components by duplicating existing components
        if len(covered_rows) < self.nr_rows:
            more_scores, more_rows = self._duplicate_columns(covered_rows)
            scores.extend(more_scores)
            score_rows.extend(more_rows)

        # sort scores based on the order of the longest allele (allele1)
        zipped_values = zip(scores, score_rows)
        scores = [value[0] for value in sorted(zipped_values, key=lambda x: x[1])]

        # remove the no matching penalty when returning scores
        for index in range(len(scores)):
            if scores[index] > self.NO_MATCHING_DECIMAL_PENALTY:
                scores[index] -= self.NO_MATCHING_DECIMAL_PENALTY
        return scores

    def _duplicate_columns(self, covered_rows):
        unused_rows = [i for i in range(len(self._rows)) if i not in covered_rows]
        scores = []
        for row_index in unused_rows:
            min_value = self._rows[row_index][0]
            min_index = 0
            for col_index, value in enumerate(self._rows[row_index]):
                if value < min_value:
                    min_value = value
                    min_index = col_index
            # making matrix complete, should go at some point
            for row in self._rows:
                row.append(row[min_index])
            scores.append(min_value)
            # make sure to add the duplication to the allele
            self.allele2.duplicate_component(min_index)
        return scores, unused_rows

    @property
    def nr_rows(self):
        return len(self._rows)

    @property
    def nr_columns(self):
        return len(self._rows[0])

    def __str__(self):
        # just for visual
        longest_column_values = []
        for index in range(self.nr_columns):
            longest_column_values.append(max(max([len(str(row[index])) for row in self._rows]),
                                             len(str(self.allele2[index]))))
        longest_row_column_value = max([len(str(value)) for value in self.allele1])
        col_names = [f"{'':<{longest_row_column_value}}"] + \
                    [f"{value: <{longest_column_values[index]}}" for index, value in enumerate(self.allele2)]
        final_str_list = [" ".join(col_names)]
        for rindex, row in enumerate(self._rows):
            formatted_row_values = [f"{self.allele1[rindex]: <{longest_row_column_value}}"]
            for cindex, value in enumerate(row):
                formatted_row_values.append(f"{value: <{longest_column_values[cindex]}}")
            final_str_list.append(" ".join(formatted_row_values))
        return '\n'.join(final_str_list)


def get_mutation_diff(
    parent_allele: Allele,
    child_allele: Allele,
    expected_size: int
) -> List[float]:
    cache_key = frozenset([frozenset(str(parent_allele)), frozenset(str(child_allele)), expected_size])
    if cache_key in SCORE_CACHE:
        return SCORE_CACHE[cache_key]

    # will permanently modify alleles in order to include duplicates
    matrix = DifferenceMatrix(parent_allele, child_allele, expected_size)
    score = matrix.get_optimal_score()

    # alleles can change so make sure to remake key
    cache_key = frozenset([frozenset(str(parent_allele)), frozenset(str(child_allele)), expected_size])
    SCORE_CACHE[cache_key] = score
    return score


@thread_termination.ThreadTerminable
def write_differentiation_rates(
    mutation_dict_list: List[Dict[str, Any]],
    distance_dict: Dict[str, Dict[str, int]],
    outfile: "Path"
):
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
                            f" not present in the distance file. This is likely caused by different names in the TGF"
                            f" files and alleles file.")
            warned_pedigrees.add(pedigree)
            continue

        if pair in distance_dict[pedigree]:
            distance = distance_dict[pedigree][pair]
        elif reverse_pair in distance_dict[pedigree]:
            distance = distance_dict[pedigree][reverse_pair]
        else:
            LOG.warning(f"Can not include pair {dictionary['To']}-{dictionary['From']} in the differentiation rate "
                        f"calculation since they are not present in the distance file.  This is likely caused by "
                        f"different names in the TGF files and alleles file.")
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
def sample_combinations(
    samples: List[str]
) -> List[Tuple[str, str]]:
    # get unique pairs of all combinations
    combinations = []
    for index, sample in enumerate(samples):
        for inner_index in range(index + 1, len(samples)):
            combinations.append((sample, samples[inner_index]))
    return combinations


@thread_termination.ThreadTerminable
def read_distance_file(
    distance_file: "Path"
) -> Dict[str, Dict[str, int]]:
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


def sort_pedigree_information(
    alleles_list_dict: List[Dict[str, Any]]
) -> Tuple[Dict[str, Dict[str, Dict[str, List[float]]]], Dict[str, Dict[str, int]]]:
    grouped_alleles_dict = {}
    longest_allele_per_pedigree_marker = {}
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
        if pedigree_name not in longest_allele_per_pedigree_marker:
            longest_allele_per_pedigree_marker[pedigree_name] = {}
        if marker not in longest_allele_per_pedigree_marker[pedigree_name]:
            longest_allele_per_pedigree_marker[pedigree_name][marker] = len(allele)
        else:
            longest_allele_per_pedigree_marker[pedigree_name][marker] = \
                max(longest_allele_per_pedigree_marker[pedigree_name][marker], len(allele))
        if pedigree_name in grouped_alleles_dict:
            if sample_name in grouped_alleles_dict[pedigree_name]:
                grouped_alleles_dict[pedigree_name][sample_name][marker] = allele
            else:
                grouped_alleles_dict[pedigree_name][sample_name] = {marker: allele}
        else:
            grouped_alleles_dict[pedigree_name] = {sample_name: {marker: allele}}
    return grouped_alleles_dict, longest_allele_per_pedigree_marker


@thread_termination.ThreadTerminable
def run(
    alleles_df: pd.DataFrame,
    distance_file: "Path",
    outdir: "Path",
    include_predict_file: bool
):

    alleles_list_dict = alleles_df.to_dict('records')

    if len(alleles_list_dict) == 0:
        LOG.error("Empty alleles file provided")
        raise utility.MalePedigreeToolboxError("Empty alleles file provided")
    total_alleles_specified = len(alleles_list_dict[0]) - 3

    # pre-sort pedigree information for quick retrieval of information
    LOG.debug("Pre-sorting alleles information")
    grouped_alleles_dict, longest_allele_per_pedigree_marker = sort_pedigree_information(alleles_list_dict)

    LOG.info("Finished reading both input files")
    markers = set(alleles_df.Marker)

    LOG.info(f"In total there are {len(markers)} markers being analysed.")
    mutation_dict = []
    total_mutation_dict = []
    predict_pedigrees_list = []

    prev_total = 0
    # make comparssons within each pedigree
    for index, (pedigree, pedigree_data) in enumerate(grouped_alleles_dict.items()):
        sample_names = list(pedigree_data.keys())
        sample_combs = sample_combinations(sample_names)
        predict_samples_list = []
        LOG.info(f"Comparing {len(sample_combs)} allele combinations for pedigree {pedigree}")
        # for each combination of samples
        for sample1, sample2 in sample_combs:
            sample1_data = pedigree_data[sample1]
            sample2_data = pedigree_data[sample2]
            total_mutations = 0
            marker_values = {name: 0 for name in markers}
            # for each individual marker
            for marker in markers:
                count_mutation = 0

                if marker not in sample1_data or marker not in sample2_data:
                    LOG.warning(f"Marker ({marker}) is not present in {sample1} and {sample2}. The comparisson will be"
                                f" skipped.")
                    continue
                marker_data1 = sample1_data[marker]
                marker_data2 = sample2_data[marker]
                mutations = get_mutation_diff(marker_data1, marker_data2,
                                              longest_allele_per_pedigree_marker[pedigree][marker])
                [mutations.append(0.0) for _ in range(total_alleles_specified - len(mutations))]

                count_mutation += np.sum(mutations)
                total_mutations += count_mutation
                if marker in marker_values and include_predict_file:
                    marker_values[marker] = count_mutation
                mutation_dict.append({"Marker": marker, "Pedigree": pedigree, "From": sample1, "To": sample2,
                                      **{f"Allele_{index + 1}": mutations[index] for index in
                                         range(total_alleles_specified)}, "Total": count_mutation})

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
    l1 = Allele([12.1, 13.1, 16])
    l2 = Allele([13, 12.1, 16])

    matrix_ = DifferenceMatrix(l2, l1, 3)
    score_ = matrix_.get_optimal_score()
    print(matrix_)
    print(score_)
