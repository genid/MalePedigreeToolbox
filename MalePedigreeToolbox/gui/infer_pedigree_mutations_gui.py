import PySimpleGUI as sg

from MalePedigreeToolbox.gui.gui_parts import TextLabel, Frame
from MalePedigreeToolbox.gui.gui_constants import LINE_LENGTH, HALFWAY_START_NR


sg.theme("Lightgrey1")


ped_mut_graph_frame = Frame(
    "Pedigree mutation graphs",
    layout=[
        [sg.Text("This module will use both the genotypic information as the pedigree information to estimate the total"
                 " number of mutations that have occured within the pedigree. The module will visualize the infered"
                 " mutations for each pedigree. Additional it will create a file with estimations of the"
                 " locus-specific mutation rates based on all provided pedigrees together.", size=(LINE_LENGTH, 5))],
        [TextLabel("Allele file"),
         sg.InputText(key="allele_pmg", size=(HALFWAY_START_NR, 1)),
         sg.FileBrowse(key="allele_pmg")],
        [sg.Text("File containing the pedigree, individual and genotypic information.", size=(LINE_LENGTH, 1))],
        [TextLabel("Tgf folder"),
         sg.InputText(size=(HALFWAY_START_NR, 1), key="tgf_folder_pmg"),
         sg.FolderBrowse(key="tgf_folder_pmg")],
        [sg.Text(
            "Folder containing at least 1 .tgf file",
            size=(LINE_LENGTH, 1)
        )],
        [TextLabel("Minimum mut. (optional)"),
         sg.InputText(size=(HALFWAY_START_NR, 1), key="minimum_mutations_pmg")],
        [sg.Text(
            "The minimum number of mutations for a pedigree to be drawn.",
            size=(LINE_LENGTH, 1)
        )],
        [TextLabel("Output folder"),
         sg.InputText(key="output_pmg", size=(HALFWAY_START_NR, 1)),
         sg.FolderBrowse(key="output_pmg")],
        [sg.Text(
            "Output directory for all files.",
            size=(LINE_LENGTH, 1)
        )]
    ],
)

layout = [[ped_mut_graph_frame]]
