# TODO add missing tests

import io
from unittest import TestCase
import logging

from MalePedigreeToolbox.tests.testing_utility import *
from MalePedigreeToolbox import mutation_diff


class Test(TestCase):

    def setUp(self) -> None:
        # ensure the stdout refers to the same stdout
        logger = logging.getLogger("mpt")
        logger.level = logging.WARNING
        self.log_capture = io.StringIO()
        self.stream_handler = logging.StreamHandler(self.log_capture)
        logger.addHandler(self.stream_handler)

    def tearDown(self) -> None:
        logger = logging.getLogger("mpt")
        logger.removeHandler(self.stream_handler)

    def _check_warning_messages(self, expected_messages):
        captured_log = self.log_capture.getvalue()
        self.log_capture.close()
        lines = captured_log.strip().split("\n")
        if len(lines) != len(expected_messages):
            self.fail("Not enough / to much warning messages generated")
        for line in lines:
            if line not in expected_messages:
                self.fail(f"Unexpected warning message recieved: {line}")


class TestMutationDiff(TestCase):

    def test_get_mutation_diff1(self):
        # a lot of cases that are expected to be true
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
        self.assertTrue(mutation_diff.get_mutation_diff(l1, l2, 2) == [1.0, 1.0])

    def test_get_mutation_diff4(self):
        l1 = [55, 63.1, 67.1]
        l2 = [54, 55, 63.1]
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