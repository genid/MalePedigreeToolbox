from typing import List, Dict, TYPE_CHECKING

import numpy as np
import pandas as pd
from tqdm import tqdm
import logging

from MalePedigreeToolbox import utility

if TYPE_CHECKING:
    from pathlib import Path


LOG = logging.getLogger("mpt")
SIMULATE_FILE = "simulation.csv"


# threading individual generations is possible --> not realy needed

def main(name_space):
    """Main starting function"""
    LOG.info("started with simulating mutations.")
    mut_rates = name_space.input
    outdir = name_space.outdir
    nr_simulated = int(name_space.num_sim)
    nr_generations = int(name_space.generations)
    if mut_rates.suffix in (".xlsx", ".xls"):
        mut_rates_df = pd.read_excel(mut_rates)
    elif mut_rates.suffix == ".csv":
        mut_rates_df = pd.read_csv(mut_rates)
    else:
        LOG.error("Invalid file type provided for mutation rate file. Either .xlsx/.xls or"
                  " .csv files are required.")
        raise ValueError("Invalid file type provided for mutation rate file. Either .xlsx/.xls or"
                         " .csv files are required.")
    generation_list = simulate_marker_mutation(nr_generations, nr_simulated, mut_rates_df)
    LOG.info("Writing simulated data to file")
    write_simulation_results(generation_list, nr_simulated, outdir)
    LOG.info("Finished simulating mutations")


def simulate_marker_mutation(
    nr_generations: int,
    individuals_per_generation: int,
    mut_rates_df: pd.DataFrame
) -> List[Dict[str, List[int]]]:
    """Generate marker mutations for each individual per generation based on predefined rates. The mutations are ran
    independant for each generation"""
    generation_list = []
    for current_generation_number in tqdm(range(1, nr_generations + 1)):
        generation_dict = {}
        for row in mut_rates_df.iterrows():
            try:
                marker = row[1].Marker
                mutation_rate_0 = row[1].Prob_0
                mutation_rate_1 = row[1].Prob_1
                mutation_rate_2 = row[1].Prob_2
            except AttributeError:
                LOG.error("Missing column names. Expected a column 'Marker' containing marker names and 3 columns "
                          "'Prob_0', 'Prob_1' and 'Prob_2' giving the probability for no mutation, a one step mutation,"
                          " and a 2 step mutation.")
                raise utility.MalePedigreeToolboxError("Missing column names. Expected a column 'Marker' containing "
                                                       "marker names and 3 columns 'Prob_0', 'Prob_1' and 'Prob_2' "
                                                       "giving the probability for no mutation, a one step mutation,"
                                                       " and a 2 step mutation.")
            mutation_options = [0, 1, 2]
            mutation_chances = [mutation_rate_0, mutation_rate_1, mutation_rate_2]
            mutation_list = [0 for _ in range(individuals_per_generation)]
            # re-simulate for each new generation to make events independant
            for _ in range(current_generation_number):
                mutation_amnt = list(np.random.choice(mutation_options, individuals_per_generation, p=mutation_chances))
                for index, amnt in enumerate(mutation_amnt):
                    # only apply mutation if actually mutated
                    if amnt == 0:
                        continue
                    mutate_amnt = mutate(amnt)
                    mutation_list[index] += mutate_amnt
                    mutation_list[index] = abs(mutation_list[index])
            generation_dict[marker] = mutation_list
        generation_list.append(generation_dict)
    return generation_list


def mutate(amnt):
    mutation_options = [-amnt, amnt]
    mutation_chances = [0.5, 0.5]
    return np.random.choice(mutation_options, 1, p=mutation_chances)[0]


def write_simulation_results(
    generation_list: List[Dict[str, List[int]]],
    individuals_per_generation: int,
    outdir: "Path"
):
    """Writes all information to a file. Marker order will be random. For large simulations this can take a while"""
    all_markers = set(generation_list[0].keys())

    copy_markers = {}
    for marker in all_markers:
        if "_" in marker:
            non_copy_name = marker.split("_")[0]
            if non_copy_name in copy_markers:
                copy_markers[non_copy_name].append(marker)
            else:
                copy_markers[non_copy_name] = [marker]

    for marker, coppied_names in copy_markers.items():
        for c_marker in coppied_names:
            all_markers.remove(c_marker)
        all_markers.add(marker)

    output_text = [f"sample,{','.join(all_markers)},Dist\n"]

    for generation, dictionary in enumerate(generation_list):
        # first de-duplicate copys
        for marker, coppied_names in copy_markers.items():
            total_values = [0 for _ in range(individuals_per_generation)]
            for copy_marker in coppied_names:
                values = dictionary.pop(copy_marker)
                for index, value in enumerate(values):
                    total_values[index] += value
            dictionary[marker] = total_values

        # then record the markers in the file
        line = ""
        for index in range(individuals_per_generation):
            line += f"{generation + 1}_{index + 1},"
            for marker in all_markers:
                line += f"{dictionary[marker][index]},"
            line += f"{generation + 1}\n"
        output_text.append(line)
    with open(outdir / SIMULATE_FILE, "w") as f:
        f.write(''.join(output_text))
