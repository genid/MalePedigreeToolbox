import PySimpleGUI as sg

from MalePedigreeToolbox.gui.gui_parts import TextLabel, Frame
from MalePedigreeToolbox.gui.gui_constants import LINE_LENGTH, HALFWAY_START_NR


sg.theme("Lightgrey1")


distance_frame = Frame(
    "Pairwise distance calculation",
    layout=[
        [sg.Text("This module calculates the pairwise distances (in meioses) of all individuals in a predefined "
                 "pedigree. The pedigree structure needs to be provided in TFG format as can be done easily using, "
                 "e.g., yEd. Multiple TGF files in a single folder can be analyzed together. Individuals with genotypic"
                 " data available need to receive a label, while individuals without genotypic data should remain "
                 "unlabeled.", size=(LINE_LENGTH, 5))],
        [TextLabel("Tgf folder"),
         sg.InputText(size=(HALFWAY_START_NR, 1), key="tgf_folder_d"),
         sg.FolderBrowse(key="tgf_folder_d")],
        [sg.Text(
            "Folder containing at least 1 .tgf file",
            size=(LINE_LENGTH, 1)
        )],
        [TextLabel("Output folder"),
         sg.InputText(key="output_d", size=(HALFWAY_START_NR, 1)),
         sg.FolderBrowse(key="output_d")],
        [sg.Text(
            "Output directory for all files.",
            size=(LINE_LENGTH, 1)
        )],
    ],
)

layout = [[distance_frame]]
