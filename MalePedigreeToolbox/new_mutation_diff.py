from typing import List
from math import isclose


class Marker:
    NO_MATCHING_DECIMAL_PENALTY: int = 1000

    def __init__(self, alleles: List[float]):
        self._alleles = alleles

    def score(self, other_marker: "Marker") -> "DifferenceMatrix":

        matrix = DifferenceMatrix()
        for allele1 in self._alleles:
            matrix_row = []
            for allele2 in other_marker._alleles:
                difference = abs(allele1 - allele2)
                if not isclose(difference, round(difference)):
                    difference += self.NO_MATCHING_DECIMAL_PENALTY
                matrix_row.append(difference)
            if len(other_marker) > len(self):
                matrix.add_row(matrix_row)
            else:
                matrix.add_column(matrix_row)

        return matrix

    def __len__(self):
        return len(self._alleles)


class DifferenceMatrix:

    def __init__(self):
        self._rows = []

    def add_row(self, row: List[float]):
        self._rows.append(row)

    def add_column(self, column: List[float]):
        if len(self._rows) == 0:
            for value in column:
                self._rows.append([value])
        else:
            for index, value in enumerate(column):
                self._rows[index].append(value)

    def __str__(self):
        return '\n'.join(map(str, self._rows))


def get_mutation_diff(
    marker_parent_alleles: List[float],
    marker_children_alleles: List[List[float]],
    expected_size: int
) -> List[List[float]]:
    parent_marker = Marker(marker_parent_alleles)
    children_markers = [Marker(alleles) for alleles in marker_children_alleles]
    for marker in children_markers:
        print(parent_marker.score(marker))


if __name__ == '__main__':
    marker1 = [1.0, 1.1, 2.0]
    marker2 = [1, 2.1]
    get_mutation_diff(marker1, [marker2], 2)
