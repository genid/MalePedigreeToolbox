import PySimpleGUI as sg

from MalePedigreeToolbox.gui.gui_parts import TextLabel, Frame
from MalePedigreeToolbox.gui.gui_constants import LINE_LENGTH, HALFWAY_START_NR


sg.theme("Lightgrey1")


predict_frame = Frame(
    "Predict generations",
    layout=[
        [sg.Text(
            "Predict the generational distance between 2 individuals based on the "
            "number of mutations between them.",
            size=(LINE_LENGTH, 3)
        )],
        [TextLabel("Input file"),
         sg.InputText(key="input_pr", size=(HALFWAY_START_NR, 1)),
         sg.FileBrowse(key="input_pr")],
        [sg.Text(
            "Mutation rates table to predict. Columns are marker names, Rows as"
            " samples in CSV or TSV format. This file can be generated from pedigrees trough "
            "the pairwise_mutation command instead.",
            size=(LINE_LENGTH, 2)
        )],
        [TextLabel("Prediction model (opt)"),
         sg.Combo(values=["RMPLEX", "PPY23", "YFP", "PPY23_RMPLEX", "YFP_RMPLEX", "YFORGEN", "YFORGEN_RMPLEX"],
                  key="model_choice_pr", readonly=True, default_value='RMPLEX', size=(18, 1))],
        [sg.Text(
            "Select a pre-made prediction model, or select a custom model below.", size=(LINE_LENGTH, 2))],
        [TextLabel("Custom model (optional)"),
            sg.InputText(key="custom_model_pr", size=(HALFWAY_START_NR, 1)),
            sg.FileBrowse(key="custom_model_pr")],
        [sg.Text(
            "A path leading to a joblib dumped model. If a custom model is provided the model choice is ignored",
            size=(LINE_LENGTH, 1)
        )],
        [TextLabel("Training file (optional)"),
         sg.InputText(key="training_file_pr", size=(HALFWAY_START_NR, 1)),
         sg.FileBrowse(key="training_file_pr")],
        [sg.Text(
            "The file that was used to train the data, to make sure that the order of the"
            " input file is the same as the order expected by the model. This only needs to"
            " be specified when using a custom model",
            size=(LINE_LENGTH, 2)
        )],
        [TextLabel("Include plots"),
         sg.Checkbox(
             "",
             key='plots_pr',
             enable_events=True)],
        [sg.Text("Warning! If many comparisons are made, selecting this option will use a lot of resouces and it will "
                 "significantly extend the running time.", size=(LINE_LENGTH, 3))],
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
