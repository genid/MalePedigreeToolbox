
from unittest import TestCase

from MalePedigreeToolbox import mutation_diff as mutation_diff


class TestMutationDiff(TestCase):
    # test a lot of different possible combinations of alleles and expected outcomes

    def test_get_optimal_nr_mutations1(self):
        l2 = [48, 66.1]
        l1 = [48, 66.1, 67.1]
        mapping = {"l1": l1, "l2": l2}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 4)[0],  [[0, 0, 1, 0]])

    def test_get_optimal_nr_mutations2(self):
        l2 = [48, 66.1]
        l1 = [48, 66.1, 67.1]
        mapping = {"l1": l1, "l2": l2}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 5)[0],  [[0, 0, 1, 0, 0]])

    def test_get_optimal_nr_mutations3(self):
        l1 = [12]
        l2 = [13]
        mapping = {"l1": l1, "l2": l2}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 2)[0],  [[1, 1]])

    def test_get_optimal_nr_mutations4(self):
        l1 = [55, 63.1, 67.1]
        l2 = [54, 55, 63.1]
        mapping = {"l1": l1, "l2": l2}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 4)[0],  [[0, 0, 4, 1]])

    def test_get_optimal_nr_mutations5(self):
        l1 = [55, 63.1, 67.1]
        l2 = [54, 55, 63.1]
        mapping = {"l1": l1, "l2": l2}
        # in this case it makes more sense if there are 2 duplicated alleles
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 3)[0],  [[0, 0, 4, 1]])

    def test_get_optimal_nr_mutations6(self):
        l1 = [1, 0, 0, 0, 0]
        l2 = [4, 0, 0, 0, 0]
        mapping = {"l1": l1, "l2": l2}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 1)[0],  [[3]])

    def test_get_optimal_nr_mutations7(self):
        l1 = [12, 13, 18]
        l2 = [12, 18]
        mapping = {"l1": l1, "l2": l2}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 3)[0],  [[0, 1, 0]])

    def test_get_optimal_nr_mutations8(self):
        l1 = [12, 13, 14]
        l2 = [12]
        mapping = {"l1": l1, "l2": l2}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 3)[0],  [[0, 1, 2]])

    def test_get_optimal_nr_mutations9(self):
        l1 = [12, 13, 14]
        l2 = [12, 14, 18]
        mapping = {"l1": l1, "l2": l2}
        # both of these are considered correct the first one is returned
        diff = mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 3)[0]
        self.assertListEqual(diff,  [[0, 1, 4]])

    def test_get_optimal_nr_mutations10(self):
        l1 = [12, 13, 18, 19]
        l2 = [12, 18]
        mapping = {"l1": l1, "l2": l2}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 4)[0],  [[0, 1, 0, 1]])

    def test_get_optimal_nr_mutations11(self):
        l1 = [12, 13, 16]
        l2 = [13, 12, 16]
        mapping = {"l1": l1, "l2": l2}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 3)[0],  [[0, 0, 0]])

    def test_get_optimal_nr_mutations12(self):
        l1 = [13, 12.1, 16]
        l2 = [12.1, 13, 16]
        mapping = {"l1": l1, "l2": l2}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 3)[0],  [[0, 0, 0]])

    def test_get_optimal_nr_mutations13(self):
        l1 = [12.1, 13.1, 16]
        l2 = [13, 12.1, 16]
        mapping = {"l1": l1, "l2": l2}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 3)[0],  [[0, 1, 0]])

    def test_get_optimal_nr_mutations14(self):
        l1 = [12.1, 14, 16]
        l2 = [13, 12.1, 16]
        mapping = {"l1": l1, "l2": l2}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 3)[0],  [[0, 1, 0]])

    def test_get_optimal_nr_mutations15(self):
        l1 = [12]
        l2 = [13]
        mapping = {"l1": l1, "l2": l2}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 3)[0],  [[1, 1, 1]])

    def test_get_optimal_nr_mutations16(self):
        l1 = [12, 15]
        l2 = [13, 15]
        mapping = {"l1": l1, "l2": l2}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 3)[0],  [[1, 0, 0]])

    def test_get_optimal_nr_mutations17(self):
        l1 = [12.1, 13]
        l2 = [12, 13.1]
        mapping = {"l1": l1, "l2": l2}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 2)[0],  [[1, 1]])

    def test_get_optimal_nr_mutations18(self):
        l1 = [13]
        l2 = [12.1, 13]
        mapping = {"l1": l1, "l2": l2}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 2)[0],  [[1, 0]])

    def test_get_optimal_nr_mutations19(self):
        l1 = [12.1, 11]
        l2 = [11.1, 12.1, 11, 12]
        mapping = {"l1": l1, "l2": l2}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 4)[0],  [[1, 0, 0, 1]])

    def test_get_optimal_nr_mutations20(self):
        l1 = [16.2, 19.2, 0, 0]
        l2 = [16.2, 18.2, 19.2, 0]
        mapping = {"l1": l1, "l2": l2}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 3)[0],  [[0, 1, 0]])

    def test_get_optimal_nr_mutations21(self):
        l1 = [0]
        l2 = [0]
        mapping = {"l1": l1, "l2": l2}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 1)[0],  [[0]])

    def test_get_optimal_nr_mutations22(self):
        l1 = [10]
        l2 = [0]
        mapping = {"l1": l1, "l2": l2}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 1)[0],  [[10]])

    def test_get_optimal_nr_mutations23(self):
        l1 = [21]
        l2 = [22]
        mapping = {"l1": l1, "l2": l2}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 1)[0],  [[1]])

    def test_get_optimal_nr_mutations24(self):
        l1 = [47, 48, 66.1, 67.1]
        l2 = [48, 66.1]
        mapping = {"l1": l1, "l2": l2}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 4)[0],  [[1, 0, 0, 1]])

    def test_get_optimal_nr_mutations25(self):
        l1 = [48, 66.1]
        l2 = [47, 48, 66.1, 67.1]
        mapping = {"l1": l1, "l2": l2}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 4)[0],  [[1, 0, 0, 1]])

    def test_get_optimal_nr_mutations26(self):
        l2 = [48, 66.1]
        l1 = [48, 66.1]
        mapping = {"l1": l1, "l2": l2}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 4)[0],  [[0, 0, 0, 0]])

    def test_get_optimal_nr_mutations27(self):
        l1 = [16.2, 19.2]
        l2 = [16.2, 18.2, 19.2]
        mapping = {"l1": l1, "l2": l2}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 4)[0],  [[0, 1, 0, 0]])

    def test_get_optimal_nr_mutations28(self):
        l1 = [55, 63, 67]
        l2 = [54, 55, 63]
        mapping = {"l1": l1, "l2": l2}
        # in this case it makes more sense if there are 2 duplicated alleles
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 3)[0],  [[0, 0, 4, 1]])

    def test_get_optimal_nr_mutations29(self):
        l1 = [12, 16]
        l2 = [13, 15]
        mapping = {"l1": l1, "l2": l2}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 3)[0],  [[1, 1, 1]])

    def test_get_optimal_nr_mutations30(self):
        l1 = [12, 13, 18]
        l2 = [12, 13]
        mapping = {"l1": l1, "l2": l2}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 4)[0],  [[0, 0, 5, 0]])

    def test_get_optimal_nr_mutations31(self):
        l1 = [12, 14]
        l2 = [12, 15]
        mapping = {"l1": l1, "l2": l2}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 3)[0],  [[0, 1, 0]])

    def test_get_optimal_nr_mutations32(self):
        l1 = [10, 11]
        l2 = [10, 11, 12]
        l3 = [10]
        mapping = {"l1": l1, "l2": l2, "l3": l3}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2"), ("l1", "l3")], mapping, 3)[0],
                             [[0, 0, 1], [0, 1, 1]])

    def test_get_optimal_nr_mutations33(self):
        l1 = [10.1, 11.1]
        l2 = [10, 11, 12]
        mapping = {"l1": l1, "l2": l2}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 3)[0],  [[1, 1, 1]])

    def test_get_optimal_nr_mutations34(self):
        # how far do we go in matching decimals
        l1 = [10, 11.1, 11, 12]
        l2 = [10, 11.1, 12.1, 13.1]
        mapping = {"l1": l1, "l2": l2}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 4)[0],  [[0, 0, 2, 2]])

    def test_get_optimal_nr_mutations35(self):
        l1 = [14, 12]
        l2 = [12, 15]
        mapping = {"l1": l1, "l2": l2}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 3)[0],  [[1, 0, 0]])

    def test_get_optimal_nr_mutations36(self):
        l1 = [10, 11]
        l2 = [10, 11, 12]
        l3 = [10]
        l4 = [11]
        mapping = {"l1": l1, "l2": l2, "l3": l3, "l4": l4}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2"), ("l3", "l4")], mapping, 3)[0],
                             [[0, 0, 1], [1, 1, 1]])

    def test_get_optimal_nr_mutations37(self):
        l1 = [10, 11]
        l2 = [10, 11, 12]
        l3 = [10, 11]
        l4 = [10]
        mapping = {"l1": l1, "l2": l2, "l3": l3, "l4": l4}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2"), ("l3", "l4")], mapping, 3)[0],
                             [[0, 0, 1], [0, 1, 0]])

    def test_get_optimal_nr_mutations38(self):
        l1 = [10, 11, 20]
        l2 = [10, 20, 21]
        l3 = [10]
        l4 = [18]
        mapping = {"l1": l1, "l2": l2, "l3": l3, "l4": l4}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2"), ("l3", "l4")], mapping, 3)[0],
                             [[0, 9, 1], [8, 8, 8]])

    def test_get_optimal_nr_mutations39(self):
        l1 = [10, 11, 20]
        l2 = [10, 20, 21]
        l3 = [10]
        l4 = [17]
        mapping = {"l1": l1, "l2": l2, "l3": l3, "l4": l4}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2"), ("l3", "l4")], mapping, 3)[0],
                             [[0, 1, 0, 1], [7, 7, 7, 7]])

    def test_get_optimal_nr_mutations40(self):
        l1 = [55, 56, 65.1]
        l2 = [55, 56, 64.1, 65.1]
        mapping = {"l1": l1, "l2": l2}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 4)[0],  [[0, 0, 1, 0]])

    def test_get_optimal_nr_mutations41(self):
        # order of inputs matters
        l1 = [13, 17.2, 18]
        l2 = [13, 16, 17.2]
        l3 = [13, 16, 17.2, 18.2]
        mapping = {"l1": l1, "l2": l2, "l3": l3}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2"), ("l1", "l3"), ("l2", "l3")],
                                                                    mapping, 4)[0],
                             [[0, 0, 2, 0], [0, 2, 0, 1], [0, 0, 0, 1]])

    def test_get_optimal_nr_mutations42(self):
        l1 = [13, 16, 17.2]
        l2 = [13, 17.2, 18]
        mapping = {"l1": l1, "l2": l2}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 4)[0],  [[0, 2, 0, 0]])

    def test_get_optimal_nr_mutations43(self):
        l1 = [13, 17]
        l2 = [16, 18]
        mapping = {"l1": l1, "l2": l2}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 2)[0],  [[3, 1]])

    def test_get_optimal_nr_mutations44(self):
        l1 = [58.2, 61, 64]
        l2 = [58.2, 61, 63, 64]
        mapping = {"l1": l1, "l2": l2}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 4)[0],  [[0, 0, 1, 0]])

    def test_get_optimal_nr_mutations45(self):
        l1 = [58.2, 61, 64]
        l2 = [58.2, 61, 63, 64]
        l3 = [58.2, 60.2, 64]
        mapping = {"l1": l1, "l2": l2, "l3": l3}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l3"), ("l1", "l2")],
                                                                    mapping, 4)[0],  [[0, 3, 0, 2, 0], [0, 0, 0, 0, 1]])

    def test_get_optimal_nr_mutations46(self):
        l1 = [58.2, 61, 64]
        l2 = [58.2, 60.2, 64]
        mapping = {"l1": l1, "l2": l2}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 4)[0],  [[0, 3, 0, 2]])

    def test_get_optimal_nr_mutations47(self):
        l1 = [58.2, 60.2, 64]
        l2 = [58.2, 61, 64]
        mapping = {"l1": l1, "l2": l2}
        self.assertListEqual(mutation_diff.get_optimal_nr_mutations([("l1", "l2")], mapping, 4)[0],  [[0, 2, 0, 3]])
