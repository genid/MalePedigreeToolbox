import PySimpleGUI as sg

from MalePedigreeToolbox.gui.gui_parts import TextLabel, Frame
from MalePedigreeToolbox.gui.gui_constants import LINE_LENGTH, HALFWAY_START_NR


sg.theme("Lightgrey1")


all_frame = Frame(
    "Full analysis",
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
            "File containing the pedigree, individual and genotipic information.",
            size=(LINE_LENGTH, 1)
        )],
        [TextLabel("Mutation rate file (optional)"),
         sg.InputText(size=(HALFWAY_START_NR, 1), key="marker_rate_all"),
         sg.FileBrowse(key="marker_rate_all")],
        [sg.Text(
            "A .csv file with two columns containing the marker name and mutation rate of each Y-STR.",
            size=(LINE_LENGTH, 2)
        )],
        [TextLabel("Nr. of clusters (optional)"),
         sg.InputText(size=(HALFWAY_START_NR, 1), key="clusters_all")],
        [sg.Text(
            "The number of clusters you want each dendogram to be divided in or type 'opt' for the optimal clustering "
            "based on silhoute score.",
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
        [TextLabel("Output folder"),
         sg.InputText(key="output_all", size=(HALFWAY_START_NR, 1)),
         sg.FolderBrowse(key="output_all")],
        [sg.Text(
            "Output directory for all files.",
            size=(LINE_LENGTH, 1)
        )]
    ],
)

layout = [[all_frame]]
