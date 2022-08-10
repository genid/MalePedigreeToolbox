import PySimpleGUI as sg

from MalePedigreeToolbox.gui.gui_parts import TextLabel, Frame
from MalePedigreeToolbox.gui.gui_constants import LINE_LENGTH, HALFWAY_START_NR


sg.theme("Lightgrey1")


mut_diff_frame = Frame(
    "Mutation differentiation calculation",
    layout=[
        [sg.Text("This module estimates the number of mutations that may have occured between each pair that belong "
                 "to the same pedigree. When combined with the output file from the pairwise distance calculation this"
                 " module will also calculate the differentiation rates. Lastly this module can produce an input file "
                 "for the prediction module.", size=(LINE_LENGTH, 5))],
        [TextLabel("Allele file"),
         sg.InputText(key="allele_md", size=(HALFWAY_START_NR, 1)),
         sg.FileBrowse(key="allele_md")],
        [sg.Text(
            "File containing the pedigree, individual and genotipic information.",
            size=(LINE_LENGTH, 1)
        )],
        [TextLabel("Pairwise distance file (opt)"),
         sg.InputText(key="distance_md", size=(HALFWAY_START_NR, 1)),
         sg.FileBrowse(key="distance_md")],
        [TextLabel("Generate prediction file"),
         sg.Checkbox(
             "",
             key=f'predict_file_md',
             enable_events=True)],
        [TextLabel("Output folder"),
         sg.InputText(key="output_md", size=(HALFWAY_START_NR, 1)),
         sg.FolderBrowse(key="output_md")],
        [sg.Text(
            "Output directory for all files.",
            size=(LINE_LENGTH, 1)
        )],
    ],
)

layout = [[mut_diff_frame]]
