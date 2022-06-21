
from unittest import TestCase

from MalePedigreeToolbox import new_mutation_diff as mutation_diff


class TestMutationDiff(TestCase):
    # test a lot of different possible combinations of alleles and expected outcomes

    def test_get_mutation_diff1(self):
        l2 = mutation_diff.Allele([48, 66.1])
        l1 = mutation_diff.Allele([48, 66.1, 67.1])
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 4) == [0.0, 0.0, 1.0, 0.0])

    def test_get_mutation_diff2(self):
        l2 = mutation_diff.Allele([48, 66.1])
        l1 = mutation_diff.Allele([48, 66.1, 67.1])
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 5) == [0.0, 0.0, 1.0, 0.0, 0.0])

    def test_get_mutation_diff3(self):
        l1 = mutation_diff.Allele([12])
        l2 = mutation_diff.Allele([13])
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 2) == [1.0, 1.0])

    def test_get_mutation_diff4(self):
        l1 = mutation_diff.Allele([55, 63.1, 67.1])
        l2 = mutation_diff.Allele([54, 55, 63.1])
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 4) == [0.0, 0.0, 4.0, 1.0])

    def test_get_mutation_diff5(self):
        l1 = mutation_diff.Allele([55, 63.1, 67.1])
        l2 = mutation_diff.Allele([54, 55, 63.1])
        # in this case it makes more sense if there are 2 duplicated alleles
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 3) == [0.0, 0.0, 4.0, 1.0])

    def test_get_mutation_diff6(self):
        l1 = mutation_diff.Allele([1.0, 0.0, 0.0, 0.0, 0.0])
        l2 = mutation_diff.Allele([4.0, 0.0, 0.0, 0.0, 0.0])
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 1) == [3.0])

    def test_get_mutation_diff7(self):
        l1 = mutation_diff.Allele([12, 13, 18])
        l2 = mutation_diff.Allele([12, 18])
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 3) == [0.0, 1.0, 0.0])

    def test_get_mutation_diff8(self):

        l1 = mutation_diff.Allele([12, 13, 14])
        l2 = mutation_diff.Allele([12])
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 3) == [0.0, 1.0, 2.0])

    def test_get_mutation_diff9(self):
        l1 = mutation_diff.Allele([12, 13, 14])
        l2 = mutation_diff.Allele([12, 14, 18])
        # both of these are considered correct the first one is returned
        diff = mutation_diff.get_mutation_diff(l1, l2, 3)
        self.assertTrue(diff == [0.0, 5.0, 0.0] or diff == [0.0, 1.0, 4.0])

    def test_get_mutation_diff10(self):
        l1 = mutation_diff.Allele([12, 13, 18, 19])
        l2 = mutation_diff.Allele([12, 18])
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 4) == [0.0, 1.0, 0.0, 1.0])

    def test_get_mutation_diff11(self):
        l1 = mutation_diff.Allele([12, 13, 16])
        l2 = mutation_diff.Allele([13, 12, 16])
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 3) == [0.0, 0.0, 0.0])

    def test_get_mutation_diff12(self):
        l1 = mutation_diff.Allele([13, 12.1, 16])
        l2 = mutation_diff.Allele([12.1, 13, 16])
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 3) == [0.0, 0.0, 0.0])

    # def test_get_mutation_diff13(self):
    #     # TODO ask what to do here
    #     l1 = mutation_diff.Allele([12.1, 13.1, 16])
    #     l2 = mutation_diff.Allele([13, 12.1, 16])
    #     self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 3) == [0.0, 1.0, 0.0])

    def test_get_mutation_diff14(self):
        l1 = mutation_diff.Allele([12.1, 14, 16])
        l2 = mutation_diff.Allele([13, 12.1, 16])
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 3) == [0.0, 1.0, 0.0])

    def test_get_mutation_diff15(self):
        l1 = mutation_diff.Allele([12])
        l2 = mutation_diff.Allele([13])
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 3) == [1.0, 1.0, 1.0])

    def test_get_mutation_diff16(self):
        l1 = mutation_diff.Allele([12, 15])
        l2 = mutation_diff.Allele([13, 15])
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 3) == [1.0, 0.0, 0.0])

    def test_get_mutation_diff17(self):
        l1 = mutation_diff.Allele([12.1, 13.0])
        l2 = mutation_diff.Allele([12.0, 13.1])
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 2) == [1.0, 1.0])

    def test_get_mutation_diff18(self):
        l1 = mutation_diff.Allele([13])
        l2 = mutation_diff.Allele([12.1, 13])
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 2) == [1.0, 0.0])

    def test_get_mutation_diff19(self):
        l1 = mutation_diff.Allele([12.1, 11])
        l2 = mutation_diff.Allele([11.1, 12.1, 11, 12])
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 4) == [1.0, 0.0, 0.0, 1.0])

    def test_get_mutation_diff20(self):
        l1 = mutation_diff.Allele([16.2, 19.2, 0.0, 0.0])
        l2 = mutation_diff.Allele([16.2, 18.2, 19.2, 0.0])
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 3) == [0.0, 1.0, 0.0])

    def test_get_mutation_diff21(self):
        l1 = mutation_diff.Allele([0.0])
        l2 = mutation_diff.Allele([0.0])
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 1) == [0.0])

    def test_get_mutation_diff22(self):
        l1 = mutation_diff.Allele([10.0])
        l2 = mutation_diff.Allele([0.0])
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 1) == [10.0])

    def test_get_mutation_diff23(self):
        l1 = mutation_diff.Allele([21.0])
        l2 = mutation_diff.Allele([22])
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 1) == [1.0])

    def test_get_mutation_diff24(self):
        l1 = mutation_diff.Allele([47, 48, 66.1, 67.1])
        l2 = mutation_diff.Allele([48, 66.1])
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 4) == [1.0, 0.0, 0.0, 1.0])

    def test_get_mutation_diff25(self):
        l1 = mutation_diff.Allele([48, 66.1])
        l2 = mutation_diff.Allele([47, 48, 66.1, 67.1])
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 4) == [1.0, 0.0, 0.0, 1.0])

    def test_get_mutation_diff26(self):
        l2 = mutation_diff.Allele([48, 66.1])
        l1 = mutation_diff.Allele([48, 66.1])
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 4) == [0.0, 0.0, 0.0, 0.0])

    def test_get_mutation_diff27(self):
        l1 = mutation_diff.Allele([16.2, 19.2])
        l2 = mutation_diff.Allele([16.2, 18.2, 19.2])
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 4) == [0.0, 1.0, 0.0, 0.0])

    def test_get_mutation_diff28(self):
        l1 = mutation_diff.Allele([55, 63, 67])
        l2 = mutation_diff.Allele([54, 55, 63])
        # in this case it makes more sense if there are 2 duplicated alleles
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 3) == [0.0, 0.0, 4.0, 1.0])

    def test_get_mutation_diff29(self):
        l1 = mutation_diff.Allele([12, 16])
        l2 = mutation_diff.Allele([13, 15])
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 3) == [1.0, 1.0, 1.0])

    def test_get_mutation_diff30(self):
        l1 = mutation_diff.Allele([12, 13, 18])
        l2 = mutation_diff.Allele([12, 13])
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 4) == [0.0, 0.0, 5.0, 0.0])

    def test_get_mutation_diff31(self):
        l1 = mutation_diff.Allele([12, 14])
        l2 = mutation_diff.Allele([12, 15])
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 3) == [0.0, 1.0, 0.0])

    def test_get_mutation_diff32(self):
        l1 = mutation_diff.Allele([10, 11])
        l2 = mutation_diff.Allele([10, 11, 12])
        l3 = mutation_diff.Allele([10])
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 3) == [0.0, 1.0, 2.0])
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l3, 3) == [0.0, 1.0, 2.0])

    def test_get_mutation_diff33(self):
        l1 = mutation_diff.Allele([10.1, 11.1])
        l2 = mutation_diff.Allele([10, 11, 12])
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 3) == [0.0, 1.0, 2.0])
