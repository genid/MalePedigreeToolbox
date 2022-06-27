#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Module for running all other modules in tandem
Author: Bram van Wersch

NOTE:
IMPORTANT: Python3.6 or higher needs to be used because some of the modules rely on sorted dictionaries
"""

# library imports
import logging
from typing import TYPE_CHECKING
from pathlib import Path

# own imports
from MalePedigreeToolbox import distances
from MalePedigreeToolbox import mutation_diff
from MalePedigreeToolbox import predict_pedigrees
from MalePedigreeToolbox import infer_pedigree_mutations
from MalePedigreeToolbox import thread_termination

if TYPE_CHECKING:
    import argparse


LOG: logging.Logger = logging.getLogger("mpt")


@thread_termination.ThreadTerminable
def main(
    name_space: "argparse.Namespace"
):
    LOG.info("Running all modules in tandem...")
    main_out_folder = Path(name_space.outdir).resolve()

    # set the variables for the distance module in the namespace
    distance_out = main_out_folder / distances.DISTANCE_FILE_NAME
    LOG.info("")
    LOG.info("Step 1/4")
    distances.main(name_space)

    # set the variables for the mutation module in the namespace
    name_space.dist_file = distance_out
    LOG.info("")
    LOG.info("Step 2/4")
    mutation_diff.main(name_space)

    # set the variables for the dendogram module in the namespace
    name_space.full_marker_csv = main_out_folder / mutation_diff.FULL_OUT
    name_space.outdir = main_out_folder
    LOG.info("")
    LOG.info("Step 3/4")
    predict_pedigrees.main(name_space)

    # set the variables for the pedigree graphs module in the namespace --> all values are already set
    LOG.info("")
    LOG.info("Step 4/4")
    infer_pedigree_mutations.main(name_space)

    LOG.info("Finished running all modules")
