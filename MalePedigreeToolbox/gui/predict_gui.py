import PySimpleGUI as sg

from MalePedigreeToolbox.gui.gui_parts import TextLabel, Frame
from MalePedigreeToolbox.gui.gui_constants import LINE_LENGTH, HALFWAY_START_NR


sg.theme("Lightgrey1")


predict_frame = Frame(
    "Predict",
    layout=[
        [sg.Text(
            "Predict the generational distance between 2 individuals based on the "
            "number of mutations between them.",
            size=(LINE_LENGTH, 2)
        )],
        [TextLabel("Input file"),
         sg.InputText(key="input_pr", size=(HALFWAY_START_NR, 1)),
         sg.FileBrowse(key="input_pr")],
        [sg.Text(
            "Mutation rates table to predict. Columns are marker names, Rows as"
            " samples in CSV or TSV format. This file can be generated from pedigrees trough "
            "the mut_diff command instead.",
            size=(LINE_LENGTH, 2)
        )],
        [TextLabel("Model choice"),
         sg.Combo(values=["RMPLEX", "PPY23", "YFP", "PPY23_RMPLEX", "YFP_RMPLEX"], key="model_choice_pr", readonly=True,
                  default_value='RMPLEX')],
        [sg.Text(
            "Choose an already calculated model based on a 5 million datapoints dataset. "
            "The names of the models refer to sets of markers used when creating the models."
            " If you want to know if your set applies to a model try it out and an error will"
            " be raised if not all required markers are present.",
            size=(LINE_LENGTH, 3)
        )],
        [TextLabel("Include probability plots"),
         sg.Checkbox(
             "",
             key='plots_pr',
             enable_events=True)],
        [sg.Text(
            "Specify this argument if all comparissons in the input file "
            "should be plotted, showing the age range at 85, 95 and 99 % "
            "confidence. If a lot of values are predicted a lot of plots are "
            "generated, so be warned, it can take a while.",
            size=(LINE_LENGTH, 3)
        )],
        [TextLabel("Outdir"),
         sg.InputText(key="output_pr", size=(HALFWAY_START_NR, 1)),
         sg.FolderBrowse(key="output_pr")],
        [sg.Text(
            "Output directory for all files.",
            size=(LINE_LENGTH, 1)
        )],
    ],
)

layout = [[predict_frame]]
