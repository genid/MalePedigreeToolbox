import io
from unittest import TestCase
import logging

from MalePedigreeToolbox.tests.testing_utility import *
from MalePedigreeToolbox import distances


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

    def test_read_graph(self):
        graph, id_name_link = distances.read_graph(TEST_FILE_DIR / "distances_tgfs" / "distance_test.tgf")
        # check warning messages
        expected_messages = {"File distance_test.tgf contains an edge with an unknown node.",
                             "File distance_test.tgf contains an edge between 3 node(s). An edge should be "
                             "between 2 nodes only."}
        self._check_warning_messages(expected_messages)

        # check output
        self.assertDictEqual(graph, {'1': {'2', '3'}, '2': {'1', '3'}, '3': {'1', '2'}, '10': set(), '100': set()})
        # the loose nodes not connected to the graph are not notified. I think that is fine
        self.assertDictEqual(id_name_link, {'1': 'a', '2': 'b', '100': '5'})

    def test_get_distance1(self):
        dists = distances.get_distances(
            {'1': {'2', '3'}, '2': {'1', '3'}, '3': {'1', '2'}, '10': set(), '100': set()},
            {'1': 'a', '2': 'b', '100': '5'})

        expected_messages = {"Encountered a circular relation. Aborting distance calculation."}
        self._check_warning_messages(expected_messages)
        self.assertEqual(dists, None)

    def test_get_distance2(self):
        dists = distances.get_distances(
            {'1': {'2'}, '2': {'1'}, '10': set(), '100': set()},
            {'1': 'a', '2': 'b', '100': '5'})
        self._check_warning_messages({""})
        self.assertEqual(dists, [('a', 'b', 1)])
