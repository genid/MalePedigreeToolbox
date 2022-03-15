import PySimpleGUI as sg

from MalePedigreeToolbox.gui.gui_parts import TextLabel, Frame
from MalePedigreeToolbox.gui.gui_constants import LINE_LENGTH, HALFWAY_START_NR


sg.theme("Lightgrey1")


ped_mut_graph_frame = Frame(
    "Pedigree mutation graphs",
    layout=[
        [sg.Text(
            "Infer mutation events from pedigrees with incomplete allele data, in addition draw these pedigrees.",
            size=(LINE_LENGTH, 3)
        )],
        [TextLabel("Allele file"),
         sg.InputText(key="allele_pmg", size=(HALFWAY_START_NR, 1)),
         sg.FileBrowse(key="allele_pmg")],
        [sg.Text(
            "File containing allele frequencies in CSV format.",
            size=(LINE_LENGTH, 1)
        )],
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
            "Folder path to store all outputs",
            size=(LINE_LENGTH, 1)
        )]
    ],
)

layout = [[ped_mut_graph_frame]]
