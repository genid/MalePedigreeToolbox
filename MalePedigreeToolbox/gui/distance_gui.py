import PySimpleGUI as sg

from MalePedigreeToolbox.gui.gui_parts import TextLabel, Frame
from MalePedigreeToolbox.gui.gui_constants import LINE_LENGTH, HALFWAY_START_NR


sg.theme("Lightgrey1")


distance_frame = Frame(
    "Tgf distance calculation",
    layout=[
        [sg.Text(
            "The distance command allows to calculate the distance for all .tgf files in a folder",
            size=(LINE_LENGTH, 2)
        )],
        [TextLabel("Tgf folder"),
         sg.InputText(size=(HALFWAY_START_NR, 1), key="tgf_folder_d"),
         sg.FolderBrowse(key="tgf_folder_d")],
        [sg.Text(
            "Folder containing at least 1 .tgf file",
            size=(LINE_LENGTH, 1)
        )],
        [TextLabel("Output directory"),
         sg.InputText(key="output_d", size=(HALFWAY_START_NR, 1)),
         sg.FolderBrowse(key="output_d")],
        [sg.Text(
            "Output directory.",
            size=(LINE_LENGTH, 1)
        )],
    ],
)

layout = [[distance_frame]]
