#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Infer mutation events from pedigrees with incomplete allele data, in addition draw these pedigrees.

general approach:

1. estimate alleles based on known information for all individuals missing allele data
2. Test all starting alleles and count minum needed number of mutations
3. choose the optimal starting allele
4. draw pedigree

Author: Diego and Bram
"""

# library imports
import logging
from pathlib import Path
import graphviz
import sys
import os
from statsmodels.stats.proportion import proportion_confint
from typing import List, Union, Dict, Tuple, Any
from collections import defaultdict
import numpy as np
import pandas as pd

# own imports
from MalePedigreeToolbox import thread_termination
from MalePedigreeToolbox import utility
from MalePedigreeToolbox import mutation_diff

# determine if application is a script file or frozen exe
if getattr(sys, 'frozen', False):
    application_path = Path(sys._MEIPASS)  # noqa needed in order for pyinstaller to find packaged files
elif __file__:
    application_path = os.path.dirname(__file__)
else:
    raise SystemExit("Can not find application path.")

os.environ["PATH"] += os.pathsep + f"{application_path}{os.sep}Graphviz{os.sep}bin"

LOG = logging.getLogger("mpt")
# color to be assign to a mutation --> this lists forces a consistent order of color use starting from the back
COLORS = ['darkorchid', 'coral2', 'azure3', 'bisque4', 'brown3', 'darkseagreen4', 'darkseagreen1', 'bisque1',
          'darkolivegreen4', 'cornsilk', 'brown1', 'cornsilk1', 'darksalmon', 'aquamarine1', 'bisque', 'cornsilk3',
          'darkslategray3', 'darkolivegreen2', 'coral1', 'darkslategray4', 'darkorange1', 'brown2', 'darkorchid2',
          'brown4', 'antiquewhite1', 'darkgoldenrod', 'darkkhaki', 'chocolate3', 'chocolate1', 'darkolivegreen',
          'darkred', 'aqua', 'darkcyan', 'coral3', 'cyan', 'aquamarine4' 'darkolivegreen3', 'aliceblue', 'darkorchid3',
          'darkgoldenrod4', 'darkgrey', 'cornsilk2', 'darkseagreen3', 'darkslategray1', 'brown', 'darkturquoise',
          'cadetblue1', 'cornsilk4', 'crimson', 'blanchedalmond', 'darkslateblue', 'burlywood', 'chocolate4', 'bisque3',
          'darkseagreen2', 'chocolate2', 'darkorange', 'bisque2', 'darkorange4', 'burlywood2', 'darkslategrey',
          'coral', 'aquamarine2', 'burlywood3', 'darkorange3', 'coral4', 'aquamarine3',
          'darkseagreen', 'darkolivegreen1', 'beige', 'cornflowerblue', 'darkorchid4', 'darkmagenta', 'cadetblue4',
          'burlywood4', 'darkorange2', 'darkorchid1', 'cadetblue', 'chartreuse4', 'chocolate',
          'burlywood1', 'darkgreen', 'aquamarine', 'darkgoldenrod2', 'plum', 'orange', 'beige', 'brown',
          'yellow', 'gray', 'blue', 'pink', 'green', 'purple', 'red']

# local testing specific values
TESTING = False
MARKER = "marker1"
PEDIGREE = "correct_split"


@thread_termination.ThreadTerminable
def main(name_space):
    """Starting point of the script"""
    LOG.info("Start with caclulating mutations from pedigrees")

    alleles_file = name_space.allele_file
    tgf_dir = name_space.tgf_folder
    minimum_mutations = name_space.minimum_mutations
    output_folder = Path(name_space.outdir)

    alleles_df = pd.read_csv(alleles_file, dtype={'Pedigree': str, 'Sample': str, 'Marker': str,
                                                  'Allele_1': np.float64, 'Allele_2': np.float64,
                                                  'Allele_3': np.float64, 'Allele_4': np.float64,
                                                  'Allele_5': np.float64, 'Allele_6': np.float64})
    pedigrees = list(set(alleles_df.Pedigree.values))  # all unique pedigree names

    # for testing specific pedigree
    if TESTING:
        pedigrees = [PEDIGREE]  # temporary testing

    files = Path(tgf_dir).absolute().glob("*.tgf")
    tgf_files = {}
    for tgf in files:
        tgf_files[tgf.stem] = tgf

    all_marker_mutation_dict = {}
    prev_total = 0
    for pedigree_index, pedigree in enumerate(pedigrees):
        if pedigree not in tgf_files:
            LOG.warning(f"Pedigree {pedigree} from alleles file not found in the tgf files. Skipping...")
            continue
        output_path = get_output_path(output_folder, pedigree)
        infer_pedigree_mutations(pedigree, alleles_df, tgf_files, all_marker_mutation_dict, output_path,
                                 minimum_mutations)

        total, remainder = divmod(pedigree_index / len(pedigrees), 0.05)
        if total != prev_total:
            LOG.info(f"Calculation progress: {round(5 * total)}%...")
            prev_total = total
    write_marker_mutations(all_marker_mutation_dict, output_folder / f"total_mutations.csv")

    LOG.info("Finished calculating mutations from pedigrees")


@thread_termination.ThreadTerminable
def parse_tgf(
    tgf_text: str
) -> "PedigreeMarkerGraph":
    """Parse a tgf file into a PedigreeDict object"""
    marker_graph = PedigreeMarkerGraph()
    in_sample_mapping_section = True
    node_mapping = {}
    sample_mapping = {}
    for elem in tgf_text.split("\n"):
        elem = elem.rstrip().split(" ")
        if elem[0] == '#' or elem[0] == "":
            in_sample_mapping_section = False
        else:
            if len(elem) <= 1:
                continue  # no associated sample
            if elem[0] not in node_mapping:
                node_mapping[elem[0]] = PedigreeNode(elem[0])
            parent_node = node_mapping[elem[0]]
            if not in_sample_mapping_section:
                if elem[1] not in node_mapping:
                    node_mapping[elem[1]] = PedigreeNode(elem[1])
                child_node = node_mapping[elem[1]]
                marker_graph.add_connection(parent_node, child_node)
            else:
                sample_mapping[parent_node] = elem[1]
    for node, sample_name in sample_mapping.items():
        marker_graph.add_sample(node, sample_name)
    return marker_graph


class FullPedgreeGraph:
    """Save information of all markers with all alleles for making a complete plot"""
    def __init__(self, pedigree_marker_graph: "PedigreeMarkerGraph"):
        self.store, self.id_node_mapping = self._innitiate_store(pedigree_marker_graph)
        self.base_node = next(iter(self.store.keys()))  # assumes sorted dicts of python > 3.6

    def _innitiate_store(
        self,
        pedigree_marker_graph: "PedigreeMarkerGraph"
    ) -> Tuple[Dict["FullPedigreeNode", List["FullPedigreeNode"]], Dict[str, "FullPedigreeNode"]]:
        next_nodes = [pedigree_marker_graph.get_base_node()]
        self_next_nodes = [FullPedigreeNode(next_nodes[0].id, next_nodes[0].sample)]
        new_store = {}
        id_mapping = {}
        while len(next_nodes) > 0:
            new_next_nodes = []
            new_self_next_nodes = []
            # add nodes
            for index, parent_node in enumerate(next_nodes):
                new_parent_node = self_next_nodes[index]
                id_mapping[new_parent_node.id] = new_parent_node
                child_nodes = pedigree_marker_graph.get_children(parent_node)
                new_child_nodes = []
                for child_node in child_nodes:
                    new_child_node = FullPedigreeNode(child_node.id, child_node.sample)
                    new_child_nodes.append(new_child_node)
                    id_mapping[new_child_node.id] = new_child_node
                new_store[new_parent_node] = new_child_nodes
                new_next_nodes.extend(child_nodes)
                new_self_next_nodes.extend(new_child_nodes)
            self_next_nodes = new_self_next_nodes
            next_nodes = new_next_nodes
        return new_store, id_mapping

    def merge_pedigree_marker_dict(
        self,
        pedigree_marker_graph: "PedigreeMarkerGraph",
        marker: str
    ):
        """Add a PedigreeMarkerGraph to this object"""
        # if this fails the merge should never happen
        assert self.base_node.id == pedigree_marker_graph.get_base_node().id
        next_nodes = [pedigree_marker_graph.get_base_node()]

        while len(next_nodes) > 0:
            new_next_nodes = []
            # add nodes
            for parent_node in next_nodes:
                parent_full_node = self.id_node_mapping[parent_node.id]
                parent_full_node.add_allele(parent_node.allele, marker)
                child_nodes = pedigree_marker_graph.get_children(parent_node)
                new_next_nodes.extend(child_nodes)
            next_nodes = new_next_nodes

    def get_children(
        self,
        parent_node: "FullPedigreeNode"
    ) -> List["FullPedigreeNode"]:
        try:
            return self.store[parent_node]
        except KeyError:
            return []


class FullPedigreeNode:
    def __init__(self, id_, sample):
        self.id = id_
        self.sample = sample
        self.alleles = {}

    def add_allele(self, allele, marker):
        if isinstance(allele, list):
            allele = mutation_diff.Locus(allele)
        self.alleles[marker] = allele

    def allele_id(self):
        # create an id that is consistent for pedigrees with the same alleles
        markers = sorted(self.alleles.keys())
        allele_marker_string_list = []
        for marker in markers:
            if self.alleles[marker] is None:
                continue
            allele_marker_string_list.append(marker + str(self.alleles[marker]))
        if len(allele_marker_string_list) == 0:
            return None
        return hash(tuple(allele_marker_string_list))

    def __repr__(self):
        return f"<id: {self.id}, sample: {self.sample}, alleles: {self.alleles}>"


class PedigreeMarkerGraph:
    """
    Store a tgf file in a easily navigable object. Containing the graph stored in both directions for quick look ups.
    This object is meant to store allele information for one marker.
    """

    id_node_mapping: Dict[str, "PedigreeNode"]

    def __init__(self):
        self.store = {}
        self.id_node_mapping = {}  # link id to node
        self.sample_node_mapping = {}  # link sample to node --> only when sample is present
        self.reverse_store = {}
        self.alleles = set()
        self._longest_allele = 0

    def check_alleles(self):
        # check alleles based on the current named alleles and the best way these alleles can be distributed
        name_node_mapping = {}
        for id_, node in self.id_node_mapping.items():
            if node.allele is not None:
                name_node_mapping[id_] = node.allele.components.copy()
        all_id_combinations = mutation_diff.sample_combinations(list(name_node_mapping.keys()))
        _, optimal_alleles = mutation_diff.get_optimal_nr_mutations(all_id_combinations, name_node_mapping,
                                                                    self._longest_allele)
        # set all the optimal alleles
        for id_ in optimal_alleles:
            self.id_node_mapping[id_].set_allele(optimal_alleles[id_])

    def add_connection(self, parent_node: "PedigreeNode", child_node: "PedigreeNode"):
        if parent_node.allele is not None and len(parent_node.allele) > self._longest_allele:
            self._longest_allele = len(parent_node.allele)
        if child_node.allele is not None and len(child_node.allele) > self._longest_allele:
            self._longest_allele = len(child_node.allele)
        self.id_node_mapping[parent_node.id] = parent_node
        self.id_node_mapping[child_node.id] = child_node
        if parent_node in self.store:
            self.store[parent_node].append(child_node)
        else:
            self.store[parent_node] = [child_node]
        # because it is a one way graph this will work
        self.reverse_store[child_node] = parent_node

    def add_sample(self, node: "PedigreeNode", sample_name):
        node.set_sample(sample_name)

        self.sample_node_mapping[sample_name] = node

    def remove_sample(self, sample_name: str):
        node = self.sample_node_mapping[sample_name]
        node.set_sample(None)
        del self.sample_node_mapping[sample_name]

    def add_allele(self, node: "PedigreeNode", allele: Union[List[int], mutation_diff.Locus]):
        node.set_allele(allele)
        self.alleles.add(node.allele)

    def total_included_nodes(self):
        # get all nodes that have an allele that is not None
        return len([node for node in self.id_node_mapping.values() if node.allele is not None])

    def get_base_node(self):
        base_nodes = set(self.id_node_mapping.values()) - set(self.reverse_store.keys())
        if len(base_nodes) != 1:
            LOG.error("Invalid tgf file provided. The file does not contain a pedigree and or the graph is circular.")
            raise utility.MalePedigreeToolboxError("Invalid tgf file provided. The file does not contain a pedigree "
                                                   "and or the graph is circular.")
        return base_nodes.pop()

    def get_sample_names(self) -> List[str]:
        return list(self.sample_node_mapping.keys())

    def get_sample_nodes(self) -> List["PedigreeNode"]:
        return list(self.sample_node_mapping.values())

    def get_children(self, parent_node: "PedigreeNode") -> List["PedigreeNode"]:
        try:
            return self.store[parent_node]
        except KeyError:
            return []

    def get_node_by_id(self, id_):
        return self.id_node_mapping[id_]

    def get_node_by_sample(self, sample):
        return self.sample_node_mapping[sample]

    def get_parent(self, child_node: "PedigreeNode") -> Union["PedigreeNode", None]:
        try:
            return self.reverse_store[child_node]
        except KeyError:
            return None

    def get_final_nodes(self):
        # nodes at the bottom of the pedigree without children
        bottom_nodes = []
        for node in iter(self.id_node_mapping.values()):
            if node not in self.store:
                bottom_nodes.append(node)
        return bottom_nodes

    def __iter__(self):
        return iter(self.id_node_mapping.values())

    def __repr__(self):
        return f"{repr(self.store)}\n{repr(self.reverse_store)}"


class PedigreeNode:
    """tracks information for a node in the pedigree graph"""

    allele: Union[mutation_diff.Locus, None]

    def __init__(self, id_):
        self.id = id_
        self.sample = None
        self.allele = None

    def set_sample(self, sample):
        self.sample = sample

    def set_allele(self, allele: Union[List[int], mutation_diff.Locus]):
        if isinstance(allele, list) or isinstance(allele, tuple):
            self.allele = mutation_diff.Locus(allele)
        else:
            self.allele = allele

    def __repr__(self):
        return f"<id: {self.id}, sample:{self.sample}, allele:{self.allele}>"


@thread_termination.ThreadTerminable
def get_output_path(
    output_folder: Path,
    pedigree: str
) -> Path:
    """Get/create the ouptut folder for a given pedigree"""
    output_path = None
    try:
        pedigree = utility.sanitize_string(pedigree)
        output_path = output_folder / pedigree
        os.mkdir(output_path)
    except FileExistsError:
        pass
    except IOError as e:
        LOG.error(f"Failed to create subfolder in output folder; {output_folder}. Error code: {str(e)}")
        raise utility.MalePedigreeToolboxError("Failed to create subfolder in output folder; {output_folder}."
                                               " Error code: {str(e)}")
    LOG.debug("Output directory has been created")
    return output_path


@thread_termination.ThreadTerminable
def infer_pedigree_mutations(
    pedigree: str,
    alleles_df: pd.DataFrame,
    tgf_files: Dict[str, Path],
    all_marker_mutation_dict: Dict[str, Dict[str, int]],
    output_path: Path,
    minimum_mutations: int
):
    """Process a pedigree and all markers and estimate the location of mutations in the pedigree"""
    LOG.info(f"Processing pedigree {pedigree}")

    df_pedigree = alleles_df.loc[alleles_df["Pedigree"] == pedigree]  # not great is very slow with large files
    full_pedigree_graph = None

    with open(tgf_files[pedigree]) as f:
        tgf = f.read()
    # all nodes found in the pedigree_marker_graph
    markers = set(df_pedigree.Marker)  # all markers found in the alleles file

    # for testing specific marker --> used for testing this module only
    if TESTING:
        if MARKER != "":
            markers = {MARKER}  # temporary for testing
    marker_mutation_dict = {}
    for marker in markers:
        LOG.debug(f"Starting analysing marker {marker} for pedigree {pedigree}")
        pedigree_marker_graph = parse_tgf(tgf)

        # get all allele lines with specific marker
        marker_rows = alleles_df.loc[(alleles_df["Pedigree"] == pedigree) & (alleles_df["Marker"] == marker)]
        if marker_rows.size <= 0:
            LOG.warning(f"Marker {marker} not found for pedigree {pedigree}")
            continue
        # remove samples names not matching between alleles and tgf file
        keep_samples = set(pedigree_marker_graph.get_sample_names()) & set(marker_rows["Sample"].values)
        allele_remove_samples = set(marker_rows["Sample"].values) - keep_samples
        tgf_remove_samples = set(pedigree_marker_graph.get_sample_names()) - keep_samples
        # filter the marker rows
        marker_rows = marker_rows.loc[~marker_rows["Sample"].isin(allele_remove_samples)]
        for sample_name in tgf_remove_samples:
            pedigree_marker_graph.remove_sample(sample_name)
        marker_rows.fillna(0, inplace=True)
        if len(marker_rows) == 0:
            LOG.warning(f"No known alleles left after filtering for marker {marker}. Skipping...")
            continue
        for _, row in marker_rows.iterrows():
            row_values = list(row.values)
            sample_name = row_values[1]
            allele = tuple(row_values[3:])
            # set the known alleles in the graph
            node = pedigree_marker_graph.get_node_by_sample(sample_name)
            pedigree_marker_graph.add_allele(node, allele)

        # assign the optimal alleles based on the currently known information in the pedigree
        pedigree_marker_graph.check_alleles()

        if full_pedigree_graph is None:
            full_pedigree_graph = FullPedgreeGraph(pedigree_marker_graph)

        # estimate the number of allele changes needed for all unknown alleles
        allele_counts = estimate_allele_changes(pedigree_marker_graph)

        # use these counts to estimate alleles of unknowns
        estimate_unknown_alleles(pedigree_marker_graph, allele_counts)
        outfile = output_path / f"pedigree_{pedigree}_marker_{marker}"
        total_mutations = plot_pedigree(pedigree_marker_graph, allele_counts, outfile, pedigree, marker,
                                        minimum_mutations)

        # save values for csv file
        add_marker_mutations(marker_mutation_dict, marker, pedigree_marker_graph.total_included_nodes(),
                             total_mutations)
        add_marker_mutations(all_marker_mutation_dict, marker, pedigree_marker_graph.total_included_nodes(),
                             total_mutations)
        full_pedigree_graph.merge_pedigree_marker_dict(pedigree_marker_graph, marker)
    write_marker_mutations(marker_mutation_dict, output_path / f"{pedigree}_mutations.csv")
    if full_pedigree_graph is not None:
        plot_full_pedigree(full_pedigree_graph, output_path / f"pedigree_{pedigree}_all_markers", pedigree)
    else:
        LOG.warning(f"Not making a full graph for pedigree {pedigree} since no markers had known alleles.")


@thread_termination.ThreadTerminable
def estimate_allele_changes(
    pedigree_marker_graph: "PedigreeMarkerGraph"
) -> Dict["PedigreeNode", Dict[List[int], int]]:
    """Estimate the number of mutations individuals needs to make to accomodate the alleles of its children

    This function traverses the pedigree_marker_graph starting from the final nodes back to the base node. For each
    node of which the alleles are unknown it is tracked how allele changes are neccesairy starting from that node
    downwards
    """
    # get mutation numbers for children in all unknown nodes
    allele_counts = defaultdict(lambda: defaultdict(int))
    covered_sample_nodes = set()
    for child_node in pedigree_marker_graph.get_final_nodes():
        allele = child_node.allele
        while True:
            parent_node = pedigree_marker_graph.get_parent(child_node)
            if parent_node is None:
                break  # reached the top of the graph
            elif allele is None:  # unknown allele due to difference in tgf and alleles file --> try the parent
                allele = parent_node.allele
            elif child_node in covered_sample_nodes:
                break
            elif parent_node.sample is None:
                # parent has no known allele
                allele_counts[parent_node][allele] += 1
                total_children = len(pedigree_marker_graph.get_children(parent_node))
                allele_counts[parent_node][allele] = min(total_children, allele_counts[parent_node][allele])
                if allele_counts[parent_node][allele] > 1:
                    break
            else:
                # parent and child have both known alleles
                if child_node.allele != parent_node.allele:
                    # encountered a mutation in known alleles. Disregard collected information and reset to the parent
                    #  allele
                    allele = parent_node.allele
            if child_node.sample is not None:
                covered_sample_nodes.add(child_node)
            child_node = parent_node
    return allele_counts


@thread_termination.ThreadTerminable
def estimate_unknown_alleles(
    pedigree_marker_graph: "PedigreeMarkerGraph",
    needed_allele_changes: Dict["PedigreeNode", Dict[mutation_diff.Locus, int]]
):
    """Estimate allelels based on the number of estimated allele changes for all possible alleles for a node without
    known allele

    First all alleleles for the starting node are tested, then the alleles are set a final time for the starting allele
    with the lowest number of mutations
    """
    if len(needed_allele_changes) == 0:  # we know the full graph so skip this step
        return
    base_node = pedigree_marker_graph.get_base_node()

    best_score = 1_000_000  # is close to the total number of mutations in the graph
    best_allele = None
    expected_allele_number = longest_allele_length(pedigree_marker_graph)

    if base_node in needed_allele_changes:
        alleles = needed_allele_changes[base_node]
    else:
        alleles = [base_node.allele]

    # check all possible beginning alleles to see what gives the least score
    for max_allele in alleles:
        pedigree_marker_graph.add_allele(base_node, max_allele)
        score = set_alleles(pedigree_marker_graph, needed_allele_changes, base_node, expected_allele_number)
        if score < best_score:
            best_score = score
            best_allele = max_allele

    pedigree_marker_graph.add_allele(base_node, best_allele)
    score = set_alleles(pedigree_marker_graph, needed_allele_changes, base_node, expected_allele_number)
    LOG.debug(f"At least {score} mutations are needed for this pedigree")


@thread_termination.ThreadTerminable
def set_alleles(
    pedigree_marker_graph: "PedigreeMarkerGraph",
    needed_allele_changes: Dict["PedigreeNode", Dict[mutation_diff.Locus, int]],
    base_node: "PedigreeNode",
    expected_allele_number: int
) -> int:
    """Set alleles for all nodes. For a node with unknown allele the following aproach is used to choose the allle:

    Get the difference between the parent allelle and all possible child alleleles saved in needed_allele_changes. Then
    simply choose the allele that needs the most allele changes minus the difference. We choose this allele to avoid
    having to make all the allele changes that we would have to make staying on the parent allele.

    A score is returned that indicates the minum number of mutations that need to be made given the alleles that
    where set.
    """
    score = 0
    next_nodes = [base_node]
    while len(next_nodes) > 0:
        new_next_nodes = []
        for parent_node in next_nodes:
            parent_allele = parent_node.allele
            child_nodes = pedigree_marker_graph.get_children(parent_node)
            for child_node in child_nodes:
                if child_node not in needed_allele_changes:
                    if child_node.allele is not None:
                        difference = sum(mutation_diff.get_mutation_diff(parent_allele, child_node.allele,
                                                                         expected_allele_number))
                        score += difference
                    continue
                id_entry_child = needed_allele_changes[child_node]
                path_scores = {}
                all_zero = True  # if this remains true there is no difference from staying on the current allele
                for child_allele in id_entry_child:
                    difference = sum(mutation_diff.get_mutation_diff(parent_allele, child_allele,
                                                                     expected_allele_number))
                    # number of mutations that are avoided choosing this one
                    avoided_mutations = id_entry_child[child_allele] - difference
                    if avoided_mutations != 0:
                        all_zero = False
                    path_scores[child_allele] = avoided_mutations, difference
                if not all_zero or len(path_scores) == 1:  # if only one possible allele might as well switch right now
                    max_allele, take_score = max(path_scores.items(), key=lambda x: x[1][0])
                    score += take_score[1]
                    pedigree_marker_graph.add_allele(child_node, max_allele)
                else:
                    pedigree_marker_graph.add_allele(child_node, parent_allele)
            new_next_nodes.extend(child_nodes)
        next_nodes = new_next_nodes
    return score


@thread_termination.ThreadTerminable
def longest_allele_length(
    pedigree_marker_graph: "PedigreeMarkerGraph"
) -> int:
    """Get the lenght longest allele for all alleles in the graph"""
    longest = 0
    for allele in pedigree_marker_graph.alleles:
        allele = np.array(allele)
        size = len(allele[allele != 0])
        if size > longest:
            longest = size
    return longest


@thread_termination.ThreadTerminable
def plot_pedigree(
    pedigree_marker_graph: "PedigreeMarkerGraph",
    needed_allele_changes: Dict["PedigreeNode", Dict[List[int], int]],
    output_file: Path,
    pedigree: str,
    marker: str,
    minimum_mutations: int
):
    """Use graphiz to draw a Digraph of the pedigree with mutations colored in it."""
    dot = graphviz.Digraph(comment='Pedigree_')
    colors = COLORS.copy()

    next_nodes = [pedigree_marker_graph.get_base_node()]
    allele_color_mapping = {tuple(next_nodes[0].allele.components): "white"}
    expected_allele_number = longest_allele_length(pedigree_marker_graph)

    # track this in order to draw edges the same color that could have had a certain mutation as well
    mutated_parents = []
    total_mutations = 0  # not the greatest place but is most efficient
    while len(next_nodes) > 0:
        new_next_nodes = []
        new_mutated_parents = []
        # add nodes
        for parent_node in next_nodes:
            color = get_allele_color(colors, allele_color_mapping, tuple(parent_node.allele.components))
            if parent_node.sample is not None:
                dot.node(parent_node.id, parent_node.sample, shape='box', style='filled', fillcolor=color)
            else:
                dot.node(parent_node.id, "", shape='box')
            child_nodes = [node for node in pedigree_marker_graph.get_children(parent_node) if node.allele is not None]
            # add edges
            for child_node in child_nodes:
                if child_node.allele != parent_node.allele:
                    edge_color = get_edge_color(colors, allele_color_mapping, tuple(child_node.allele.components))

                    distance = sum(mutation_diff.get_mutation_diff(parent_node.allele, child_node.allele,
                                                                   expected_allele_number))
                    new_mutated_parents.append(child_node)
                    total_mutations += distance
                    dot.edge(parent_node.id, child_node.id, label=str(int(distance)), color=edge_color)
                elif parent_node in needed_allele_changes and len(needed_allele_changes[parent_node]) == 1:
                    edge_color = get_edge_color(colors, allele_color_mapping, tuple(parent_node.allele.components))
                    if parent_node in mutated_parents:
                        dot.edge(parent_node.id, child_node.id, color=edge_color)
                        new_mutated_parents.append(child_node)
                    else:
                        dot.edge(parent_node.id, child_node.id)
                else:
                    dot.edge(parent_node.id, child_node.id)
            new_next_nodes.extend(child_nodes)
        next_nodes = new_next_nodes
        mutated_parents = new_mutated_parents
    dot.attr(label=f"Pedigree {pedigree} for marker {marker}", labelloc="t", fontsize="30")

    # only draw when needded
    if total_mutations >= minimum_mutations:
        with dot.subgraph(name="cluster_legend",) as c:
            prev_allele = None
            c.attr(rankdir="LR", rank="same",  label="legend")
            for allele, color in allele_color_mapping.items():
                str_allele = format_allele(allele)
                # add this to make sure there is no way that an id and allele interfere
                c.node(str_allele + " ", label=f"allele: {str_allele}", shape='box', style='filled', fillcolor=color)
                if prev_allele is not None:
                    c.edge(prev_allele + " ", str_allele + " ", style="invis", weight="10", minlen="1")
                prev_allele = str_allele

        dot.render(output_file, view=False)
        # remove the graphiz dot file because it creates extra clutter
        try:
            os.remove(output_file)
        except FileNotFoundError:
            pass
    return int(total_mutations)


def get_allele_color(
    colors: List[str],
    allele_color_mapping: Dict[Any, str],
    allele: Union[Tuple[int], int]
) -> str:
    """Get the color for an allele. Return a color from the mapping if already present otherwise map a new color"""
    if allele is None:  # special cases of values present in tgf but not allele file
        return "white"
    if allele not in allele_color_mapping:
        allele_color_mapping[allele] = colors.pop()
    return allele_color_mapping[allele]


def get_edge_color(
    colors: List[str],
    allele_color_mapping: Dict[Any, str],
    allele: Union[Tuple[int], int]
) -> str:
    edge_color = get_allele_color(colors, allele_color_mapping, allele)
    if edge_color == "white":
        edge_color = "black"
    return edge_color


@thread_termination.ThreadTerminable
def plot_full_pedigree(
    full_pedigree_marker_graph: "FullPedgreeGraph",
    output_file: Path,
    pedigree: str
):
    """Plot the pedigree of all markers combined"""
    dot = graphviz.Digraph(comment='Pedigree_')
    colors = COLORS.copy()

    next_nodes = [full_pedigree_marker_graph.base_node]
    if next_nodes[0].sample is not None:
        dot.node(next_nodes[0].id, next_nodes[0].sample, shape='box', style='filled', fillcolor="white")
    else:
        dot.node(next_nodes[0].id, "", shape='box', style='filled', fillcolor="white")

    allele_color_mapping = {next_nodes[0].allele_id(): "white"}
    edge_change_dict = {}

    while len(next_nodes) > 0:
        new_next_nodes = []
        # add nodes
        for parent_node in next_nodes:
            child_nodes = [node for node in full_pedigree_marker_graph.get_children(parent_node) if
                           node.allele_id() is not None]
            parent_alleles_dict = parent_node.alleles
            # add edges
            for child_node in child_nodes:
                edge_changes = []
                for marker, child_allele in child_node.alleles.items():
                    parent_allele = parent_alleles_dict[marker]
                    if child_allele is not None and child_allele != parent_allele:
                        edge_changes.append(f"{marker},{format_allele(parent_allele)},{format_allele(child_allele)}")

                # dont show label if no mutations
                if len(edge_changes) == 0:
                    dot.edge(parent_node.id, child_node.id)
                else:
                    edge_id = get_edge_id(len(edge_change_dict))
                    edge_change_dict[edge_id] = edge_changes
                    dot.edge(parent_node.id, child_node.id, label=f"{edge_id}: {len(edge_changes)}")

                if len(edge_changes) == 0:
                    # if nothing changed give same color. The markers can be different because of differences between
                    # allele file and tgf file.
                    color = get_allele_color(colors, allele_color_mapping, parent_node.allele_id())
                else:
                    color = get_allele_color(colors, allele_color_mapping, child_node.allele_id())
                if child_node.sample is not None:
                    dot.node(child_node.id, child_node.sample, shape='box', style='filled', fillcolor=color)
                else:
                    dot.node(child_node.id, "", shape='box', style='filled', fillcolor=color)

            new_next_nodes.extend(child_nodes)
        next_nodes = new_next_nodes
    dot.attr(label=f"Pedigree {pedigree} for all markers", labelloc="t", fontsize="30")
    dot.render(output_file, view=False)
    # remove the graphiz dot file because it creates extra clutter
    try:
        os.remove(output_file)
    except FileNotFoundError:
        pass

    write_full_pedigree_edge_information(edge_change_dict,
                                         output_file.parent / f"pedigree_{pedigree}_all_marker_edge_info.csv")


def format_allele(
    allele: Union[Tuple[float], mutation_diff.Locus]
) -> str:
    if isinstance(allele, mutation_diff.Locus):
        allele = allele.components
    allele_string_list = []
    for number in allele:
        number_decimal = str(number).split(".")
        if len(number_decimal) == 1:
            allele_string_list.append(str(int(number)))
        elif int(number_decimal[1]) == 0:
            allele_string_list.append(str(int(number)))
        else:
            allele_string_list.append(f"{number_decimal[0]}.{number_decimal[1]}")
    return '/'.join(allele_string_list)


def get_edge_id(
    current_id_nr: int
) -> str:
    """Convert a number into a two letter sequence. Starting from aa. This function can maximum handle 26 * 26 edges
     this should be more than enough"""
    first_id_letter = chr(int(current_id_nr / 26) + 97)
    second_id_letter = chr((current_id_nr % 26) + 97)
    return first_id_letter + second_id_letter


def add_marker_mutations(
    mutation_dict: Dict[str, Dict[str, int]],
    marker: str,
    total_individuals: int,
    total_mutations: int
):
    """Record the number of mutations in a mutation_dict"""
    if marker in mutation_dict:
        mutation_dict[marker]["m"] += total_individuals - 1
        mutation_dict[marker]["tm"] += total_mutations

    else:
        mutation_dict[marker] = {"m": total_individuals - 1, "tm": total_mutations}


@thread_termination.ThreadTerminable
def write_marker_mutations(
    mutation_dict: Dict[str, Dict[str, int]],
    outfile: Path
):
    """Write a mutation dict to a csv file"""
    text_list = ["Marker,Total_meioses,Total_mutations,Mutation_rate,CI (95% lower bound), CI (95% upper bound)"]
    for marker, nr_dict in mutation_dict.items():
        total_individuals = nr_dict["m"]
        total_mutations = nr_dict["tm"]
        # in case of no mutations in the whole pedigree
        if total_individuals == total_mutations == 0:
            text_list.append(f"{marker},{total_individuals},{total_mutations},nan,nan,nan")
            continue
        if total_individuals == 0:
            ratio = 0
        else:
            ratio = round(total_mutations / total_individuals, 5)
        ci = proportion_confint(total_mutations, total_individuals, method='beta')
        ci = list(map(lambda x: round(x, 5), ci))
        text_list.append(f"{marker},{total_individuals},{total_mutations},{ratio},{ci[0]},{ci[1]}")

    with open(outfile, "w") as f:
        f.write('\n'.join(text_list))


@thread_termination.ThreadTerminable
def write_full_pedigree_edge_information(
    edge_change_dict: Dict[str, List[str]],
    outfile: Path
):
    """Write information about the edges that are numbered in the combined pedigree plot"""
    full_text_list = ["edge_id,marker,parent_allele,child_allele"]
    for edge_id, marker_changes in edge_change_dict.items():
        for marker_change in marker_changes:
            full_text_list.append(f"{edge_id},{marker_change}")
    with open(outfile, "w") as f:
        f.write('\n'.join(full_text_list))
