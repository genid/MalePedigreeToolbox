from typing import List, Union, Set
from math import isclose, ceil
from collections import defaultdict


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
            for _ in range(expected_size - self.nr_rows):
                rindex, cindex = self._get_minimum_coordinate()
                self._rows.append(self._rows[rindex].copy())
                self.allele1.duplicate_component(rindex)
                self.allele2.duplicate_component(cindex)
                for row in self._rows:
                    row.append(row[cindex])

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

            unused_rows = [i for i in range(len(self._rows)) if i not in covered_rows]
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
                score_rows.append(row_index)
                # make sure to add the duplication to the allele
                self.allele2.duplicate_component(min_index)

        # sort scores based on the order of the longest allele (allele1)
        zipped_values = zip(scores, score_rows)
        scores = [value[0] for value in sorted(zipped_values, key=lambda x: x[1])]

        # remove the no matching penalty when returning scores
        for index in range(len(scores)):
            if scores[index] > self.NO_MATCHING_DECIMAL_PENALTY:
                scores[index] -= self.NO_MATCHING_DECIMAL_PENALTY
        return scores

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
    # will permanently modify alleles in order to include duplicates
    matrix = DifferenceMatrix(parent_allele, child_allele, expected_size)
    score = matrix.get_optimal_score()
    return score


if __name__ == '__main__':
    l1 = Allele([12.1, 13.1, 16])
    l2 = Allele([13, 12.1, 16])

    matrix_ = DifferenceMatrix(l2, l1, 3)
    score_ = matrix_.get_optimal_score()
    print(matrix_)
    print(score_)
