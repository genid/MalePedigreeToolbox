import PySimpleGUI as sg

from MalePedigreeToolbox.gui.gui_parts import TextLabel, Frame
from MalePedigreeToolbox.gui.gui_constants import LINE_LENGTH, HALFWAY_START_NR


sg.theme("Lightgrey1")


mut_diff_frame = Frame(
    "Mutation differentiation calculation",
    layout=[
        [sg.Text(
            "Calculate the number of mutations between all alleles of markers in a pedigree. Additionaly calculate "
            "mutation rates if distances between individuals in a pedigree are provided.",
            size=(LINE_LENGTH, 2)
        )],
        [TextLabel("Allele file"),
         sg.InputText(key="allele_md", size=(HALFWAY_START_NR, 1)),
         sg.FileBrowse(key="allele_md")],
        [sg.Text(
            "File containing allele frequencies in CSV format.",
            size=(LINE_LENGTH, 1)
        )],
        [TextLabel("Distance file (optional)"),
         sg.InputText(key="distance_md", size=(HALFWAY_START_NR, 1)),
         sg.FileBrowse(key="distance_md")],
        [sg.Text(
            "Distance file location. Can be generated in the distance tab. This file is optional, it will provide a way"
            " of seeing average mutations per generation plus confidence intervals.",
            size=(LINE_LENGTH, 2)
        )],
        [TextLabel("Include predict file"),
         sg.Checkbox(
             "",
             key=f'predict_file_md',
             enable_events=True)],
        [sg.Text(
            "Add this option if you want a file to be generated that can be used by the predict command in order to"
            " predict generational distance between individuals in a pedigree.",
            size=(LINE_LENGTH, 2)
        )],
        [TextLabel("Outdir"),
         sg.InputText(key="output_md", size=(HALFWAY_START_NR, 1)),
         sg.FolderBrowse(key="output_md")],
        [sg.Text(
            "Output directory for all files.",
            size=(LINE_LENGTH, 1)
        )],
    ],
)

layout = [[mut_diff_frame]]
