# TODO add missing tests

from unittest import TestCase

from MalePedigreeToolbox import mutation_diff as mutation_diff


class Test(TestCase):

    def test_get_score1(self):
        self.assertEqual(4, mutation_diff.get_score(5, 1))

    def test_get_score2(self):
        self.assertEqual(7, mutation_diff.get_score(5, -2), )

    def test_get_score3(self):
        self.assertEqual(6, mutation_diff.get_score(5.9, 0.89))

    def test_get_score4(self):
        self.assertEqual(5, mutation_diff.get_score(5.9, 0.91))

    def test_get_score5(self):
        self.assertEqual(1, mutation_diff.get_score(9, 10))

    def test_get_decimal1(self):
        self.assertEqual(1, mutation_diff.get_decimal(1.1))

    def test_get_decimal2(self):
        self.assertEqual(0, mutation_diff.get_decimal(1))

    def test_get_decimal3(self):
        self.assertEqual(0, mutation_diff.get_decimal(1.0))

    def test_get_decimal4(self):
        self.assertEqual(12345678, mutation_diff.get_decimal(1.12345678))

    def test_get_matrix_score1(self):
        ms, mdm = mutation_diff.get_matrix_scores([(1, 2.0), (3, 1.0)], [(1, 2.0), (3, 1.0)])
        self.assertEqual(ms.tolist(), [[0, 1], [1, 0]])
        self.assertEqual(mdm.tolist(), [[1, 1], [1, 1]])

    def test_get_matrix_score2(self):
        ms, mdm = mutation_diff.get_matrix_scores([(1, 2.0), (3, 1.0)], [(1, 2.0)])
        self.assertEqual(ms.tolist(), [[0], [1]])
        self.assertEqual(mdm.tolist(), [[1], [1]])

    def test_get_matrix_score3(self):
        ms, mdm = mutation_diff.get_matrix_scores([(1, 1.1), (3, 1.0)], [(1, 2.0)])
        self.assertEqual(ms.tolist(), [[1], [1]])
        self.assertEqual(mdm.tolist(), [[0], [1]])

    def test_get_matrix_score4(self):
        ms, mdm = mutation_diff.get_matrix_scores([(1, 1.12314), (3, 2.0)], [(100, 2.0)])
        self.assertEqual(ms.tolist(), [[1], [0]])
        self.assertEqual(mdm.tolist(), [[0], [1]])

    def test_get_matrix_score5(self):
        ms, mdm = mutation_diff.get_matrix_scores([(1, 2.1), (3, 1.0)], [(1, 2.0), (3, 1.1), (3, 1.0)])
        self.assertEqual(ms.tolist(), [[1, 1, 2], [1, 1, 0]])
        self.assertEqual(mdm.tolist(), [[0, 1, 0], [1, 0, 1]])

    def test_sort_allele(self):
        self.assertEqual([(0, 2), (2, 3), (1, 4)], mutation_diff.sort_allele([2, 4, 3]))


class TestMutationDiff(TestCase):
    # test a lot of different possible combinations of alleles and expected outcomes

    def test_get_mutation_diff1(self):
        l2 = [48, 66.1]
        l1 = [48, 66.1, 67.1]
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 4) == [0.0, 0.0, 1.0, 0.0])

    def test_get_mutation_diff2(self):
        l2 = [48, 66.1]
        l1 = [48, 66.1, 67.1]
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 5) == [0.0, 0.0, 1.0, 0.0, 0.0])

    def test_get_mutation_diff3(self):
        l1 = [12]
        l2 = [13]
        print(mutation_diff.get_mutation_diff(l1, l2, 2))
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 2) == [1.0, 1.0])

    def test_get_mutation_diff4(self):
        l1 = [55, 63.1, 67.1]
        l2 = [54, 55, 63.1]
        print(mutation_diff.get_mutation_diff(l1, l2, 4))
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 4) == [0.0, 0.0, 4.0, 1.0])

    def test_get_mutation_diff5(self):
        l1 = [55, 63.1, 67.1]
        l2 = [54, 55, 63.1]
        # in this case it makes more sense if there are 2 duplicated alleles
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 3) == [0.0, 0.0, 4.0, 1.0])

    def test_get_mutation_diff6(self):
        l1 = [1.0, 0.0, 0.0, 0.0, 0.0]
        l2 = [4.0, 0.0, 0.0, 0.0, 0.0]
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 1) == [3.0])

    def test_get_mutation_diff7(self):
        l1 = [12, 13, 18]
        l2 = [12, 18]
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 3) == [0.0, 1.0, 0.0])

    def test_get_mutation_diff8(self):

        l1 = [12, 13, 14]
        l2 = [12]
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 3) == [0.0, 1.0, 2.0])

    def test_get_mutation_diff9(self):
        l1 = [12, 13, 14]
        l2 = [12, 14, 18]
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 3) == [0.0, 1.0, 4.0])

    def test_get_mutation_diff10(self):
        l1 = [12, 13, 18, 19]
        l2 = [12, 18]
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 4) == [0.0, 1.0, 0.0, 1.0])

    def test_get_mutation_diff11(self):
        l1 = [12, 13, 16]
        l2 = [13, 12, 16]
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 3) == [0.0, 0.0, 0.0])

    def test_get_mutation_diff12(self):
        l1 = [13, 12.1, 16]
        l2 = [12.1, 13, 16]
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 3) == [0.0, 0.0, 0.0])

    def test_get_mutation_diff13(self):
        l1 = [12.1, 13.1, 16]
        l2 = [13, 12.1, 16]
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 3) == [0.0, 1.0, 0.0])

    def test_get_mutation_diff14(self):
        l1 = [12.1, 14, 16]
        l2 = [13, 12.1, 16]
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 3) == [0.0, 1.0, 0.0])

    def test_get_mutation_diff15(self):
        l1 = [12]
        l2 = [13]
        print(mutation_diff.get_mutation_diff(l1, l2, 3))
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 3) == [1.0, 1.0, 1.0])

    def test_get_mutation_diff16(self):
        l1 = [12, 15]
        l2 = [13, 15]
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 3) == [1.0, 0.0, 0.0])

    def test_get_mutation_diff17(self):
        l1 = [12.1, 13.0]
        l2 = [12.0, 13.1]
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 2) == [1.0, 1.0])

    def test_get_mutation_diff18(self):
        l1 = [13]
        l2 = [12.1, 13]
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 2) == [1.0, 0.0])

    def test_get_mutation_diff19(self):
        l1 = [12.1, 11]
        l2 = [11.1, 12.1, 11, 12]
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 4) == [1.0, 0.0, 0.0, 1.0])

    def test_get_mutation_diff20(self):
        l1 = [16.2, 19.2, 0.0, 0.0]
        l2 = [16.2, 18.2, 19.2, 0.0]
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 3) == [0.0, 1.0, 0.0])

    def test_get_mutation_diff21(self):
        l1 = [0.0]
        l2 = [0.0]
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 1) == [0.0])

    def test_get_mutation_diff22(self):
        l1 = [10.0]
        l2 = [0.0]
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 1) == [10.0])

    def test_get_mutation_diff23(self):
        l1 = [21.0]
        l2 = [22]
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 1) == [1.0])

    def test_get_mutation_diff24(self):
        l1 = [47, 48, 66.1, 67.1]
        l2 = [48, 66.1]
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 4) == [1.0, 0.0, 0.0, 1.0])

    def test_get_mutation_diff25(self):
        l1 = [48, 66.1]
        l2 = [47, 48, 66.1, 67.1]
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 4) == [1.0, 0.0, 0.0, 1.0])

    def test_get_mutation_diff26(self):
        l2 = [48, 66.1]
        l1 = [48, 66.1]
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 4) == [0.0, 0.0, 0.0, 0.0])

    def test_get_mutation_diff27(self):
        l1 = [16.2, 19.2]
        l2 = [16.2, 18.2, 19.2]
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 4) == [0.0, 1.0, 0.0, 0.0])

    def test_get_mutation_diff28(self):
        l1 = [55, 63, 67]
        l2 = [54, 55, 63]
        # in this case it makes more sense if there are 2 duplicated alleles
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 3) == [0.0, 0.0, 4.0, 1.0])

    def test_get_mutation_diff29(self):
        l1 = [12, 16]
        l2 = [13, 15]
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 3) == [1.0, 1.0, 1.0])

    def test_get_mutation_diff30(self):
        l1 = [12, 13, 18]
        l2 = [12, 13]
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 4) == [0.0, 0.0, 5.0, 0.0])

    def test_get_mutation_diff31(self):
        l1 = [12, 14]
        l2 = [12, 15]
        print(mutation_diff.get_mutation_diff(l1, l2, 3))
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 3) == [0.0, 1.0, 0.0])
