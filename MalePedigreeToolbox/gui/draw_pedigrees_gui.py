import PySimpleGUI as sg

from MalePedigreeToolbox.gui.gui_parts import TextLabel, Frame
from MalePedigreeToolbox.gui.gui_constants import LINE_LENGTH, HALFWAY_START_NR


sg.theme("Lightgrey1")


draw_pedigree_frame = Frame(
    "Draw Pedigrees",
    layout=[
        [sg.Text(
            "Draw plots based mutation distances in a pedigree, this can be dendograms or multi-dimensional"
            " scaling plots.",
            size=(LINE_LENGTH, 3)
        )],
        [TextLabel("Full marker file"),
         sg.InputText(size=(HALFWAY_START_NR, 1), key="full_marker_dp"),
         sg.FileBrowse(key="full_marker_dp")],
        [sg.Text(
            "The file containing full mutation differentiations this file can be generated in the mutation "
            "differentiation tab.",
            size=(LINE_LENGTH, 2)
        )],
        [TextLabel("Marker rate file (optional)"),
         sg.InputText(size=(HALFWAY_START_NR, 1), key="marker_rate_dp"),
         sg.FileBrowse(key="marker_rate_dp")],
        [sg.Text(
            "File with mutation rates of all markers present in full marker file. The expected format is a csv file "
            "with 2 columns 1. marker 2. rate. This will give more accurate dendrograms. Leave this field empty to "
            "assume the same mutation rate for all markers.",
            size=(LINE_LENGTH, 3)
        )],
        [TextLabel("Plot choice"),
         sg.Combo(values=["dendrogram", "MDS", 'both'], key="plot_choice_dp", readonly=True,
                  default_value='dendrogram')],
        [sg.Text(
            "The plot type you want.",
            size=(LINE_LENGTH, 1)
        )],
        [TextLabel("Nr. of clusters (optional)"),
         sg.InputText(size=(HALFWAY_START_NR, 1), key="clusters_dp"),
         sg.FileBrowse(key="clusters_dp")],
        [sg.Text(
            "The expected number of clusters for all pedigrees. This can be a single value to get the same number of"
            " clusters for all pedigrees or a text file containing space separated positive integers. If"
            " no value is provided the optimal clustering is calculated based on silhouette score.",
            size=(LINE_LENGTH, 4)
        )],
        [TextLabel("Random state (optional)"),
         sg.InputText(size=(HALFWAY_START_NR, 1), key="random_state_dp")],
        [sg.Text(
            "An integer representing a random start state for the MDS plot. This will ensure that consecutive runs on"
            " the same data provide the same plot.",
            size=(LINE_LENGTH, 2)
        )],
        [TextLabel("Output folder"),
         sg.InputText(key="output_dp", size=(HALFWAY_START_NR, 1)),
         sg.FolderBrowse(key="output_dp")],
        [sg.Text(
            "Folder path to store all outputs",
            size=(LINE_LENGTH, 1)
        )]
    ],
)

layout = [[draw_pedigree_frame]]
