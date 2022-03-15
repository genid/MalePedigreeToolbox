import PySimpleGUI as sg

from MalePedigreeToolbox.gui.gui_parts import TextLabel, Frame
from MalePedigreeToolbox.gui.gui_constants import LINE_LENGTH, HALFWAY_START_NR


sg.theme("Lightgrey1")


all_frame = Frame(
    "Run it all",
    layout=[
        [sg.Text(
            "Performs all of the commands in order, first calculating distance,than mutation differentiation, "
            "followed by dendograms for mutation differentiation rates and finaly caclulates mutation rates from "
            "pedigrees.",
            size=(LINE_LENGTH, 3)
        )],
        [TextLabel("Tgf folder"),
         sg.InputText(size=(HALFWAY_START_NR, 1), key="tgf_folder_all"),
         sg.FolderBrowse(key="tgf_folder_all")],
        [sg.Text(
            "Folder containing at least 1 .tgf file",
            size=(LINE_LENGTH, 1)
        )],
        [TextLabel("Allele file"),
         sg.InputText(key="allele_all", size=(HALFWAY_START_NR, 1)),
         sg.FileBrowse(key="allele_all")],
        [sg.Text(
            "File containing allele frequencies in CSV format.",
            size=(LINE_LENGTH, 1)
        )],
        [TextLabel("Marker rate file (optional)"),
         sg.InputText(size=(HALFWAY_START_NR, 1), key="marker_rate_all"),
         sg.FileBrowse(key="marker_rate_all")],
        [sg.Text(
            "File with mutation rates of all markers present in full marker file. The expected format is a csv file "
            "with 2 columns 1. marker 2. rate. This will give more accurate dendrograms. Leave this field empty to "
            "assume the same mutation rate for all markers.",
            size=(LINE_LENGTH, 3)
        )],
        [TextLabel("Plot choice"),
         sg.Combo(values=["dendrogram", "MDS", 'both'], key="plot_choice_all", readonly=True,
                  default_value='dendrogram')],
        [sg.Text(
            "The plot type you want.",
            size=(LINE_LENGTH, 1)
        )],
        [TextLabel("Nr. of clusters (optional)"),
         sg.InputText(size=(HALFWAY_START_NR, 1), key="clusters_all"),
         sg.FileBrowse(key="clusters_all")],
        [sg.Text(
            "The expected number of clusters for all pedigrees. This can be a single value to get the same number of"
            " clusters for all pedigrees or a text file containing space separated positive integers. If"
            " no value is provided the optimal clustering is calculated based on silhouette score.",
            size=(LINE_LENGTH, 4)
        )],
        [TextLabel("Random state (optional)"),
         sg.InputText(size=(HALFWAY_START_NR, 1), key="random_state_all")],
        [sg.Text(
            "An integer representing a random start state for the MDS plot. This will ensure that consecutive runs on"
            " the same data provide the same plot.",
            size=(LINE_LENGTH, 2)
        )],
        [TextLabel("Minimum mut. (optional)"),
         sg.InputText(size=(HALFWAY_START_NR, 1), key="minimum_mutations_all")],
        [sg.Text(
            "The minimum number of mutations for a pedigree to be drawn.",
            size=(LINE_LENGTH, 1)
        )],
        [TextLabel("Include predict file"),
         sg.Checkbox(
             "",
             key=f'predict_file_all',
             enable_events=True)],
        [sg.Text(
            "Add this option if you want a file to be generated that can be used by the predict command in order to"
            " predict generational distance between individuals in a pedigree.",
            size=(LINE_LENGTH, 2)
        )],
        [TextLabel("Output folder"),
         sg.InputText(key="output_all", size=(HALFWAY_START_NR, 1)),
         sg.FolderBrowse(key="output_all")],
        [sg.Text(
            "Folder path to store all outputs",
            size=(LINE_LENGTH, 1)
        )]
    ],
)

layout = [[all_frame]]
