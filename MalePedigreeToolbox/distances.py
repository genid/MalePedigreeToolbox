#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Code for calculating the distances between all named nodes in a tgf file

Author: Bram van Wersch
"""


from pathlib import Path
import logging
from typing import Dict, Tuple, Set, List, TYPE_CHECKING, Union

# own imports
from MalePedigreeToolbox import thread_termination

if TYPE_CHECKING:
    import argparse


LOG = logging.getLogger("mpt")
DISTANCE_FILE_NAME = "distances.csv"


@thread_termination.ThreadTerminable
def main(name_space: "argparse.Namespace"):
    """Main entry point"""
    LOG.info("Started with calculating pairwise distances.")
    files = Path(name_space.tgf_folder).absolute().glob("*.tgf")
    outdir = Path(name_space.outdir).absolute()
    all_distances = []
    for file in files:
        graph, nodes_of_interest = read_graph(file)
        distances = get_distances(graph, nodes_of_interest)
        if distances is None:
            continue
        else:
            all_distances.append((distances, file.name))
    write_results(all_distances, outdir)


@thread_termination.ThreadTerminable
def read_graph(
    file: Path
) -> Tuple[Dict[str, Set[str]], Dict[str, str]]:
    """Read the Trivial Graph Format (tgf) file into a dictionary of nodes as keys and a list of connected nodes as
    value.

    return a dictionary keyed on nodes with connected nodes as values and a set linking node ID's to sample names
    if applicable for a certain ID
    """
    with open(file) as f:
        lines = f.readlines()
    # read the graph into a dictionary with keys as nodes and a list of connected nodes as values
    graph = {}
    id_name_link = {}

    reading_edges = False
    for line in lines:
        # all values under this indicate edges
        if line.startswith("#"):
            reading_edges = True
            continue

        values = line.strip().split(" ")
        # reading in nodes
        if not reading_edges:
            # save all the nodes with ID's these are the nodes that we want the distances of
            if len(values) > 1:
                id_name_link[values[0]] = values[1]
            graph[values[0]] = set()
        else:
            if len(values) != 2:
                # empty line ignore and continue
                if len(values) != 1:
                    LOG.warning(f"File {file.name} contains an edge between {len(values)} node(s). "
                                f"An edge should be between 2 nodes only.")
                continue

            id1, id2 = values
            # in case an edge is referencing a non existing node ignore it
            if id1 not in graph or id2 not in graph:
                LOG.warning(f"File {file.name} contains an edge with an unknown node.")
                continue
            graph[id1].add(id2)
            graph[id2].add(id1)
    return graph, id_name_link


@thread_termination.ThreadTerminable
def get_distances(
    graph: Dict[str, Set[str]],
    nodes_of_interest: Dict[str, str]
) -> Union[List[Tuple[str, str, int]], None]:
    """Get the distance between all named nodes in the graph.

    :param graph:  a dictionary keyed on nodes with connected nodes as values
    :param nodes_of_interest: a set linking node ID's to sample names if applicable for a certain ID
    :return: a list of connecting named nodes and respective distance between them
    """
    distances = []
    covered_pairs = set()
    for id_ in nodes_of_interest:
        covered_connected = {id_}
        connected = {(i, 1) for i in graph[id_]}

        while len(connected) > 0:
            connected_id, distance = connected.pop()
            if connected_id in nodes_of_interest:
                if (connected_id, id_) not in covered_pairs:
                    distances.append((nodes_of_interest[id_], nodes_of_interest[connected_id], distance))
                    covered_pairs.add((id_, connected_id))
            covered_connected.add(connected_id)

            covered_node_ids = {v[0] for v in connected}
            for node_id in graph[connected_id]:
                if node_id in covered_node_ids:
                    LOG.warning("Encountered a circular relation. Aborting distance calculation.")
                    return None
                if node_id not in covered_connected:
                    connected.add((node_id, distance + 1))
                    covered_node_ids.add(node_id)
    return distances


def write_results(all_distances, outdir):
    distance_text = "File,From,To,Distance\n"
    for distances, file in all_distances:
        for distance in distances:
            distance_text += f"{file.replace('.tgf', '')},{distance[0]},{distance[1]},{distance[2]}\n"

    with open(outdir / DISTANCE_FILE_NAME, "w") as f:
        f.write(distance_text[:-1])  # remove the last newline
    LOG.info("Finished calculating pairwise distances")
