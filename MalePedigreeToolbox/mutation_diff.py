from typing import Union, Tuple, List, Dict, Any, TYPE_CHECKING, Set
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


class Allele:

    def __init__(self, components: Union[List[float], Tuple[float]]):
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
        new_decimal = self._get_decimal(new_value)
        self._decimals.append(new_decimal)

    def get_equalizable_decimals(self, other_allele: "Allele"):
        # decimals present in both alleles and more then one decimal
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
        decimal_count_dict = {key: value for key, value in decimal_count_dict.items() if value != 0}
        return decimal_count_dict

    def get_indexes_with_decimal(self, wanted_decimal):
        indexes = []
        for index, decimal in enumerate(self._decimals):
            if decimal == wanted_decimal:
                indexes.append(index)
        return indexes

    @property
    def components(self):
        return sorted(self._components)

    def __iter__(self):
        return iter(self._components)

    def __getitem__(self, item: int):
        return self._components[item]

    def __str__(self):
        return str(sorted(self._components))

    def __len__(self):
        return len(self._components)

    def __hash__(self):
        return hash(tuple(self.components))

    def __eq__(self, other: "Allele"):
        if not isinstance(other, Allele):
            # this is an acceptable alternative
            if other is None:
                return False
            raise ValueError(f"Only compare alleles. Cannot compare against {other}")
        return sorted(self.components) == sorted(other.components)


class DifferenceMatrix:
    NO_MATCHING_DECIMAL_PENALTY: int = 1000

    allele1: Allele
    allele2: Allele

    def __init__(self, allele1: Allele, allele2: Allele, expected_size):
        if len(allele1) == 0 or len(allele2) == 0:
            raise ValueError("Cannot compute distance against empty allele")

        if len(allele2) > len(allele1):
            temp = allele1
            allele1 = allele2
            allele2 = temp
        self.allele1 = allele1
        self.allele2 = allele2
        self.expected_size = expected_size

        # in case of equal alleles we already know the score
        if self.allele1 != self.allele2:

            self._rows = []
            self._create_matrix()
            self.score = None

    def _create_matrix(self):
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

        # equalize the lengths of both alleles
        if len(self.allele1) != len(self.allele2):
            self._equalize_allele_lengths()
        # make alleles to the expected size
        if self.nr_rows < self.expected_size:
            self._duplicate_simultanious()

    def _equalize_allele_lengths(self):
        # find the current optimal path and calculate optimal duplication depending on remaining open components

        # prioritize duplication of imbalanced decimals
        self._equalize_decimals((0, len(self.allele1) - len(self.allele2)))
        if len(self.allele1) == len(self.allele2):
            return

        # if still a difference duplicate such that the unaligned allele gets the optimal match
        sorted_score_coordinates = self._get_sorted_score_coordinates()

        # select the optimal mutations
        covered_rows = set()
        covered_columns = set()
        score_index = 0
        while len(covered_rows) < len(self._rows[0]):
            row, column = sorted_score_coordinates[score_index]
            score_index += 1
            if row in covered_rows or column in covered_columns:
                continue
            covered_rows.add(row)
            covered_columns.add(column)

        self._duplicate_in_allele2(covered_rows)

    def _duplicate_in_allele2(self, covered_rows):
        # add components to allele2 (shortest allele) in order to equalize the lengths
        unused_rows = [i for i in range(len(self._rows)) if i not in covered_rows]
        for row_index in unused_rows:
            min_value = self._rows[row_index][0]
            min_index = 0
            for col_index, value in enumerate(self._rows[row_index]):
                if value < min_value:
                    min_value = value
                    min_index = col_index
            # making matrix complete
            for row in self._rows:
                row.append(row[min_index])
            # make sure to add the duplication to the allele
            self.allele2.duplicate_component(min_index)

    def _duplicate_simultanious(self):
        # duplicate both alleles when both of them are equally to short of the expected size
        self._equalize_decimals((self.expected_size - self.nr_rows, self.expected_size - self.nr_rows))

        for _ in range(self.expected_size - self.nr_rows):
            rindex, cindex = self._min_row(list(range(self.nr_rows))), self._min_column(list(range(self.nr_columns)))

            # if we add the same allele resulting in an addition of 0 score, dont add it since there are other
            # possibilities
            if self.allele1[rindex] == self.allele2[cindex]:
                continue
            self._rows.append(self._rows[rindex].copy())
            self.allele1.duplicate_component(rindex)
            self.allele2.duplicate_component(cindex)
            for row in self._rows:
                row.append(row[cindex])

    def _equalize_decimals(self, max_changes: Tuple[int, int] = None):
        # equalize the differences in decimals for both alleles
        # max changes is the maximum changes for each allele
        equalizables = self.allele1.get_equalizable_decimals(self.allele2)
        decimal_difference_dict = self.allele1.get_decimal_difference(self.allele2, equalizables)
        if len(decimal_difference_dict) == 0:
            return

        total_changes = [0, 0]
        for decimal, difference in decimal_difference_dict.items():
            # for each difference
            for _ in range(abs(difference)):
                if difference < 0 and (max_changes is None or max_changes[0] > total_changes[0]):
                    indexes = self.allele1.get_indexes_with_decimal(decimal)
                    min_row_index = self._min_row(indexes)
                    self.allele1.duplicate_component(min_row_index)
                    self._rows.append(self._rows[min_row_index].copy())
                    total_changes[0] += 1
                elif difference > 0 and (max_changes is None or max_changes[1] > total_changes[1]):
                    indexes = self.allele2.get_indexes_with_decimal(decimal)
                    min_col_index = self._min_column(indexes)
                    self.allele2.duplicate_component(min_col_index)

                    for row in self._rows:
                        row.append(row[min_col_index])
                    total_changes[1] += 1
                elif total_changes[0] >= max_changes[0] and total_changes[1] >= max_changes[1]:
                    # in case the max changes are reached
                    return

    def _min_row(self, indexes):
        # row index with the lowest sum of values for a given number of indexes
        min_sum = sum(self._rows[indexes[0]])
        min_value = min(self._rows[indexes[0]])
        min_index = indexes[0]
        for index in indexes[1:]:
            row_sum = sum(self._rows[index])
            lowest_value = min(self._rows[index])
            if row_sum > min_sum:
                continue
            elif row_sum == min_sum and lowest_value >= min_value:
                continue
            min_sum = row_sum
            min_index = index
            min_value = lowest_value
        return min_index

    def _min_column(self, indexes):
        # column index with the lowest sum of values
        min_sum = sum([row[indexes[0]] for row in self._rows])
        min_value = min([row[indexes[0]] for row in self._rows])
        min_index = indexes[0]
        for index in indexes[1:]:
            col_values = [row[index] for row in self._rows]
            col_sum = sum(col_values)
            lowest_value = min(col_values)
            if col_sum > min_sum:
                continue
            elif col_sum == min_sum and lowest_value >= min_value:
                continue
            min_sum = col_sum
            min_index = index
            min_value = lowest_value
        return min_index

    def _get_sorted_score_coordinates(self):
        coordinates = [(i, j) for j in range(self.nr_columns) for i in range(self.nr_rows)]
        return sorted(coordinates, key=lambda coord: self._rows[coord[0]][coord[1]])

    def calculate_mutations(self) -> List[int]:
        if self.allele1 == self.allele2:
            return [0 for _ in range(self.expected_size)]
        mutations = []
        score_rows = []
        sorted_score_coordinates = self._get_sorted_score_coordinates()

        # select the optimal mutations
        covered_rows = set()
        covered_columns = set()
        score_index = 0
        while len(mutations) < len(self._rows[0]):
            row, column = sorted_score_coordinates[score_index]
            score_index += 1
            if row in covered_rows or column in covered_columns:
                continue
            covered_rows.add(row)
            covered_columns.add(column)
            mutations.append(self._rows[row][column])
            score_rows.append(row)

        # sort mutations based on the order of the longest allele (allele1)
        zipped_values = zip(mutations, score_rows)
        mutations = [value[0] for value in sorted(zipped_values, key=lambda x: x[1])]

        # remove the no matching penalty when returning mutations
        for index in range(len(mutations)):
            if mutations[index] > self.NO_MATCHING_DECIMAL_PENALTY:
                mutations[index] -= self.NO_MATCHING_DECIMAL_PENALTY
        # add 0's at the end in case of ambigious duplications that are left out
        for _ in range(self.expected_size - len(mutations)):
            mutations.append(0)
        return mutations

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


@thread_termination.ThreadTerminable
def get_mutation_diff(
    parent_allele: Allele,
    child_allele: Allele,
    expected_size: int
) -> List[float]:
    # will permanently modify alleles in order to include duplicates
    matrix = DifferenceMatrix(parent_allele, child_allele, expected_size)
    score = matrix.calculate_mutations()
    return score


@thread_termination.ThreadTerminable
def get_optimal_nr_mutations(
    all_allele_pairs: List[Tuple[str, str]],
    allele_name_mapping: Dict[str, List[float]],
    expected_size: int,
) -> Tuple[List[List[float]], Dict[str, Allele]]:
    best_mutations, best_total_mutations, best_allele_mapping = \
        _calculate_mutations(all_allele_pairs, allele_name_mapping, expected_size)

    allowed_over_duplications = 1
    while True:
        mutations, total_mutations, allele_mapping = \
            _calculate_mutations(all_allele_pairs, allele_name_mapping, expected_size + allowed_over_duplications)
        if total_mutations < best_total_mutations:
            best_mutations = mutations
            best_total_mutations = total_mutations
            best_allele_mapping = allele_mapping
        else:
            break
    return best_mutations, best_allele_mapping
        

@thread_termination.ThreadTerminable
def _calculate_mutations(
    allele_pairs: List[Tuple[str, str]],
    allele_name_mapping: Dict[str, List[float]],
    expected_size: int
):
    mutations = []
    total_mutations = 0
    allele_mapping = {}

    # create Allele objects for every allele defined as a List of floats, alleles can be modified when calculating which
    # is important information for other comparissons
    for name1, name2 in allele_pairs:
        allele_mapping[name1] = Allele(allele_name_mapping[name1])
        allele_mapping[name2] = Allele(allele_name_mapping[name2])

    for name1, name2 in allele_pairs:
        allele1 = allele_mapping[name1]
        allele2 = allele_mapping[name2]
        score = get_mutation_diff(allele1, allele2, expected_size)
        mutations.append(score)
        total_mutations += sum(score)
    return mutations, total_mutations, allele_mapping


@thread_termination.ThreadTerminable
def write_differentiation_rates(
    all_information_list: List[List[Any]],
    distance_dict: Dict[str, Dict[str, int]],
    outfile: "Path"
):
    # the rate of differentiation given a certain distance of a 2 subjects in a pair
    meiosis_dict = {}
    covered_pairs = set()
    mutated_pairs = set()
    warned_pedigrees = set()  # make sure the log is less spammy
    for lst in all_information_list:
        differentiated = lst[-1] != 0
        pedigree = lst[1]
        pair = lst[2] + lst[3]
        reverse_pair = lst[3] + lst[2]

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
            LOG.warning(f"Can not include pair {lst[3]}-{lst[2]} in the differentiation rate "
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


@thread_termination.ThreadTerminable
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

    # pre-sort pedigree information for quick retrieval of information
    LOG.debug("Pre-sorting alleles information")
    grouped_alleles_dict, longest_allele_per_pedigree_marker = sort_pedigree_information(alleles_list_dict)

    LOG.info("Finished reading both input files")
    markers = list(set(alleles_df.Marker))

    LOG.info(f"In total there are {len(markers)} markers being analysed.")
    all_information_list = []
    total_mutation_dict = {}
    predict_samples_dict = {}
    longest_allele = 0

    prev_total = 0
    # make comparssons within each pedigree
    for index, (pedigree, pedigree_data) in enumerate(grouped_alleles_dict.items()):
        sample_names = list(pedigree_data.keys())
        sample_combs = sample_combinations(sample_names)
        LOG.info(f"Comparing {len(sample_combs)} allele combinations for pedigree {pedigree}")
        for marker in markers:
            all_allele_pairs = []
            pair_name_mapping = {}
            for sample1, sample2 in sample_combs:
                sample1_data = pedigree_data[sample1]
                sample2_data = pedigree_data[sample2]

                if marker not in sample1_data or marker not in sample2_data:
                    LOG.warning(f"Marker ({marker}) is not present in {sample1} and {sample2}. The comparisson will be"
                                f" skipped.")
                    continue
                marker_data1 = sample1_data[marker]
                marker_data2 = sample2_data[marker]
                all_allele_pairs.append((sample1, sample2))
                pair_name_mapping[sample1] = marker_data1
                pair_name_mapping[sample2] = marker_data2
            if marker not in longest_allele_per_pedigree_marker[pedigree]:
                continue
            optimal_mutations, _ = get_optimal_nr_mutations(all_allele_pairs, pair_name_mapping,
                                                            longest_allele_per_pedigree_marker[pedigree][marker])
            for mutations, (sample1, sample2) in zip(optimal_mutations, sample_combs):
                if len(mutations) > longest_allele:
                    longest_allele = len(mutations)
                sum_mutations = sum(mutations)
                all_information_list.append([len(all_information_list), pedigree, sample1, sample2, marker,
                                            mutations, sum_mutations])
                key = (pedigree, sample1, sample2)
                if key not in total_mutation_dict:
                    total_mutation_dict[key] = 0
                    predict_samples_dict[key] = {}
                total_mutation_dict[key] += sum_mutations
                if include_predict_file:
                    predict_samples_dict[key][marker] = sum_mutations

        LOG.debug(f"Finished calculating differentiation for {index} out of {len(grouped_alleles_dict)}")
        total, remainder = divmod(index / len(grouped_alleles_dict), 0.05)

        if total != prev_total:
            LOG.info(f"Calculation progress: {round(5 * total)}%...")
            prev_total = total

    with open(outdir / FULL_OUT, "w") as f:
        allele_header_section = ','.join([f"Allele_{index + 1}" for index in range(longest_allele)])
        header = f"Pedigree,From,To,Marker,{allele_header_section},Total\n"
        f.write(header)
        for lst in all_information_list:
            # mutations
            lst[5] += [0 for _ in range(longest_allele - len(lst[5]))]
            lst[5] = ','.join(map(str, lst[5]))
            f.write(','.join(map(str, lst)) + "\n")

    with open(outdir / SUMMARY_OUT, "w") as f:
        f.write(f",Pedigree,From,To,Total\n")
        for index, (key, value) in enumerate(total_mutation_dict.items()):
            pedigree, from_, to = key
            f.write(f"{index},{pedigree},{from_},{to},{value}\n")

    if distance_file is not None:
        # read the distance file
        LOG.info("Started with summarising and writing meiosis differentiation rates to file")
        distance_dict = read_distance_file(distance_file)
        write_differentiation_rates(all_information_list, distance_dict, outdir / DIFFERENTIATION_OUT)

    if include_predict_file:
        with open(outdir / PREDICT_OUT, "w") as f:
            f.write(f"sample,{','.join(markers)}\n")
            for key, marker_dict in predict_samples_dict.items():
                f.write('_'.join(key) + ",")
                f.write(','.join(str(marker_dict[marker]) if marker in marker_dict else "0" for marker in markers))
                f.write("\n")
    LOG.info("Finished calculating differentiation rates.")
