
from unittest import TestCase

from MalePedigreeToolbox import new_mutation_diff as mutation_diff


class TestMutationDiff(TestCase):
    # test a lot of different possible combinations of alleles and expected outcomes

    def test_get_optimal_nr_mutations1(self):
        l2 = [48, 66.1]
        l1 = [48, 66.1, 67.1]
        self.assertTrue(mutation_diff.get_optimal_nr_mutations([(l1, l2)], 4) == [[0, 0, 1, 0]])

    def test_get_optimal_nr_mutations2(self):
        l2 = [48, 66.1]
        l1 = [48, 66.1, 67.1]
        self.assertTrue(mutation_diff.get_optimal_nr_mutations([(l1, l2)], 5) == [[0, 0, 1, 0, 0]])

    def test_get_optimal_nr_mutations3(self):
        l1 = [12]
        l2 = [13]
        self.assertTrue(mutation_diff.get_optimal_nr_mutations([(l1, l2)], 2) == [[1, 1]])

    def test_get_optimal_nr_mutations4(self):
        l1 = [55, 63.1, 67.1]
        l2 = [54, 55, 63.1]
        self.assertTrue(mutation_diff.get_optimal_nr_mutations([(l1, l2)], 4) == [[0, 0, 4, 1]])

    def test_get_optimal_nr_mutations5(self):
        l1 = [55, 63.1, 67.1]
        l2 = [54, 55, 63.1]
        # in this case it makes more sense if there are 2 duplicated alleles
        self.assertTrue(mutation_diff.get_optimal_nr_mutations([(l1, l2)], 3) == [[0, 0, 4, 1]])

    def test_get_optimal_nr_mutations6(self):
        l1 = [1, 0, 0, 0, 0]
        l2 = [4, 0, 0, 0, 0]
        self.assertTrue(mutation_diff.get_optimal_nr_mutations([(l1, l2)], 1) == [[3]])

    def test_get_optimal_nr_mutations7(self):
        l1 = [12, 13, 18]
        l2 = [12, 18]
        self.assertTrue(mutation_diff.get_optimal_nr_mutations([(l1, l2)], 3) == [[0, 1, 0]])

    def test_get_optimal_nr_mutations8(self):
        l1 = [12, 13, 14]
        l2 = [12]
        self.assertTrue(mutation_diff.get_optimal_nr_mutations([(l1, l2)], 3) == [[0, 1, 2]])

    def test_get_optimal_nr_mutations9(self):
        l1 = [12, 13, 14]
        l2 = [12, 14, 18]
        # both of these are considered correct the first one is returned
        diff = mutation_diff.get_optimal_nr_mutations([(l1, l2)], 3)
        self.assertTrue(diff == [[0, 5, 0]] or diff == [[0, 1, 4]])

    def test_get_optimal_nr_mutations10(self):
        l1 = [12, 13, 18, 19]
        l2 = [12, 18]
        self.assertTrue(mutation_diff.get_optimal_nr_mutations([(l1, l2)], 4) == [[0, 1, 0, 1]])

    def test_get_optimal_nr_mutations11(self):
        l1 = [12, 13, 16]
        l2 = [13, 12, 16]
        self.assertTrue(mutation_diff.get_optimal_nr_mutations([(l1, l2)], 3) == [[0, 0, 0]])

    def test_get_optimal_nr_mutations12(self):
        l1 = [13, 12.1, 16]
        l2 = [12.1, 13, 16]
        self.assertTrue(mutation_diff.get_optimal_nr_mutations([(l1, l2)], 3) == [[0, 0, 0]])

    def test_get_optimal_nr_mutations13(self):
        l1 = [12.1, 13.1, 16]
        l2 = [13, 12.1, 16]
        self.assertTrue(mutation_diff.get_optimal_nr_mutations([(l1, l2)], 3) == [[0, 1, 0]])

    def test_get_optimal_nr_mutations14(self):
        l1 = [12.1, 14, 16]
        l2 = [13, 12.1, 16]
        self.assertTrue(mutation_diff.get_optimal_nr_mutations([(l1, l2)], 3) == [[0, 1, 0]])

    def test_get_optimal_nr_mutations15(self):
        l1 = [12]
        l2 = [13]
        self.assertTrue(mutation_diff.get_optimal_nr_mutations([(l1, l2)], 3) == [[1, 1, 1]])

    def test_get_optimal_nr_mutations16(self):
        l1 = [12, 15]
        l2 = [13, 15]
        self.assertTrue(mutation_diff.get_optimal_nr_mutations([(l1, l2)], 3) == [[1, 0, 0]])

    def test_get_optimal_nr_mutations17(self):
        l1 = [12.1, 13]
        l2 = [12, 13.1]
        self.assertTrue(mutation_diff.get_optimal_nr_mutations([(l1, l2)], 2) == [[1, 1]])

    def test_get_optimal_nr_mutations18(self):
        l1 = [13]
        l2 = [12.1, 13]
        self.assertTrue(mutation_diff.get_optimal_nr_mutations([(l1, l2)], 2) == [[1, 0]])

    def test_get_optimal_nr_mutations19(self):
        l1 = [12.1, 11]
        l2 = [11.1, 12.1, 11, 12]
        self.assertTrue(mutation_diff.get_optimal_nr_mutations([(l1, l2)], 4) == [[1, 0, 0, 1]])

    def test_get_optimal_nr_mutations20(self):
        l1 = [16.2, 19.2, 0, 0]
        l2 = [16.2, 18.2, 19.2, 0]
        self.assertTrue(mutation_diff.get_optimal_nr_mutations([(l1, l2)], 3) == [[0, 1, 0]])

    def test_get_optimal_nr_mutations21(self):
        l1 = [0]
        l2 = [0]
        self.assertTrue(mutation_diff.get_optimal_nr_mutations([(l1, l2)], 1) == [[0]])

    def test_get_optimal_nr_mutations22(self):
        l1 = [10]
        l2 = [0]
        self.assertTrue(mutation_diff.get_optimal_nr_mutations([(l1, l2)], 1) == [[10]])

    def test_get_optimal_nr_mutations23(self):
        l1 = [21]
        l2 = [22]
        self.assertTrue(mutation_diff.get_optimal_nr_mutations([(l1, l2)], 1) == [[1]])

    def test_get_optimal_nr_mutations24(self):
        l1 = [47, 48, 66.1, 67.1]
        l2 = [48, 66.1]
        self.assertTrue(mutation_diff.get_optimal_nr_mutations([(l1, l2)], 4) == [[1, 0, 0, 1]])

    def test_get_optimal_nr_mutations25(self):
        l1 = [48, 66.1]
        l2 = [47, 48, 66.1, 67.1]
        self.assertTrue(mutation_diff.get_optimal_nr_mutations([(l1, l2)], 4) == [[1, 0, 0, 1]])

    def test_get_optimal_nr_mutations26(self):
        l2 = [48, 66.1]
        l1 = [48, 66.1]
        self.assertTrue(mutation_diff.get_optimal_nr_mutations([(l1, l2)], 4) == [[0, 0, 0, 0]])

    def test_get_optimal_nr_mutations27(self):
        l1 = [16.2, 19.2]
        l2 = [16.2, 18.2, 19.2]
        self.assertTrue(mutation_diff.get_optimal_nr_mutations([(l1, l2)], 4) == [[0, 1, 0, 0]])

    def test_get_optimal_nr_mutations28(self):
        l1 = [55, 63, 67]
        l2 = [54, 55, 63]
        # in this case it makes more sense if there are 2 duplicated alleles
        self.assertTrue(mutation_diff.get_optimal_nr_mutations([(l1, l2)], 3) == [[0, 0, 4, 1]])

    def test_get_optimal_nr_mutations29(self):
        l1 = [12, 16]
        l2 = [13, 15]
        self.assertTrue(mutation_diff.get_optimal_nr_mutations([(l1, l2)], 3) == [[1, 1, 1]])

    def test_get_optimal_nr_mutations30(self):
        l1 = [12, 13, 18]
        l2 = [12, 13]
        print(mutation_diff.get_optimal_nr_mutations([(l1, l2)], 4))
        self.assertTrue(mutation_diff.get_optimal_nr_mutations([(l1, l2)], 4) == [[0, 0, 5, 0]])

    def test_get_optimal_nr_mutations31(self):
        l1 = [12, 14]
        l2 = [12, 15]
        self.assertTrue(mutation_diff.get_optimal_nr_mutations([(l1, l2)], 3) == [[0, 1, 0]])

    def test_get_optimal_nr_mutations32(self):
        l1 = [10, 11]
        l2 = [10, 11, 12]
        l3 = [10]
        self.assertTrue(mutation_diff.get_optimal_nr_mutations([(l1, l2), (l1, l3)], 3) == [[0, 0, 1],
                                                                                            [0, 1, 1]])

    def test_get_optimal_nr_mutations33(self):
        l1 = [10.1, 11.1]
        l2 = [10, 11, 12]
        self.assertTrue(mutation_diff.get_optimal_nr_mutations([(l1, l2)], 3) == [[1, 1, 1]])

    def test_get_optimal_nr_mutations34(self):
        l1 = [10, 11.1, 11, 12]
        l2 = [10, 11.1, 12.1, 13.1]
        print(mutation_diff.get_optimal_nr_mutations([(l1, l2)], 4))
        self.assertTrue(mutation_diff.get_optimal_nr_mutations([(l1, l2)], 4) == [[0, 0, 1, 2, 1, 2]])

    def test_get_optimal_nr_mutations35(self):
        l1 = [14, 12]
        l2 = [12, 15]
        self.assertTrue(mutation_diff.get_optimal_nr_mutations([(l1, l2)], 3) == [[0, 1, 0]])
