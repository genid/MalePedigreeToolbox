# library imports
import tempfile

import PySimpleGUI as sg
from pathlib import Path
import sys
import os
from threading import Thread
import traceback

# own imports
from MalePedigreeToolbox.gui import distance_gui
from MalePedigreeToolbox.gui import mutation_differentiation_gui
from MalePedigreeToolbox.gui import infer_pedigree_mutations_gui
from MalePedigreeToolbox.gui import predict_pedigrees_gui
from MalePedigreeToolbox.gui import all_gui
from MalePedigreeToolbox.gui import predict_gui

from MalePedigreeToolbox.gui import gui_parts
from MalePedigreeToolbox.gui.gui_constants import LINE_LENGTH

from MalePedigreeToolbox import thread_termination
from MalePedigreeToolbox import main
from MalePedigreeToolbox import __version__

sg.theme("Lightgrey1")

# in case runnign from the exe
if getattr(sys, 'frozen', False):
    user_execute_file_path = Path(os.path.dirname(sys.executable))
elif __file__:
    user_execute_file_path = Path(os.getcwd())
else:
    raise SystemExit("Can not find application path.")


TESTING = False


def mpt_gui():
    layout = [
        [sg.Text("Male Pedigree Toolbox", font="Arial 18 bold", pad=(0, 0))],
        [sg.Text(f"Version: {__version__}", font="Arial 10", pad=(0, 0))],
        [sg.Text("Created by: Arwin Ralf, Diego Montiel Gonzalez and Bram van Wersch, 2022", font="Arial 10",
                 pad=(0, 0))],
        [sg.Text("", font="Arial 10",
                 pad=(0, 0))],
        [sg.TabGroup([
            [sg.Tab("Pairwise distance", [[gui_parts.Column(distance_gui.layout, scrollable=True)]])],
            [sg.Tab("Pairwise mutation",
                    [[gui_parts.Column(mutation_differentiation_gui.layout, scrollable=True)]])],
            [sg.Tab("Pedigree mutation",
                    [[gui_parts.Column(infer_pedigree_mutations_gui.layout, scrollable=True)]])],
            [sg.Tab("Dendrograms",
                    [[gui_parts.Column(predict_pedigrees_gui.layout, scrollable=True)]])],
            [sg.Tab("Full analysis", [[gui_parts.Column(all_gui.layout, scrollable=True)]])],
            [sg.Tab("Prediction", [[gui_parts.Column(predict_gui.layout, scrollable=True)]])],
        ], enable_events=True, key="mpt_tabs"
        )],
        [sg.Button("Start", key="start_button", button_color=("white", "green")),
         sg.Button("Exit", key="exit_button", button_color=("white", "red"))]
    ]

    window = sg.Window(
        "Male Pedigree Toolbox",
        layout,
        size=(700, 800),
        element_padding=(5, 5),
        element_justification="center",
        finalize=True
    )

    while True:
        event, values = window.read()

        if event in (None, "exit_button"):
            break
        if event:
            if event == "start_button":
                run_main_command(values)
    window.close()


class CommandThread(Thread):
    def __init__(self, arguments):
        super().__init__(daemon=True)
        self._arguments = arguments
        self.final_error = ""

    def run(self) -> None:
        try:
            main.main(*self._arguments, is_gui=True)
        except BaseException as e:  # use base exception to catch system-exit calls
            if TESTING:
                self.final_error = traceback.format_exc()
            else:
                self.final_error = f"Execution failed with message: {str(e)}."


def run_main_command(values):
    thread_termination.clear_semaphore()

    arguments, log_file_loc = get_command(values)

    temp = tempfile.NamedTemporaryFile(delete=False)
    if log_file_loc is None:
        # since we will fail because the output is not specified, make a temp file that can satisfy some other processes
        # before the application inevitably crashes
        log_file_loc = temp.name

    command_thread = CommandThread(arguments)
    command_thread.start()

    layout = [[sg.Multiline("Welcome to the Male Pedigree Toolbox!", key="multiline", size=(LINE_LENGTH, 25),
                            disabled=True, autoscroll=True)],
              [sg.Button("OK", key="ok_button", disabled=True),
               sg.Button("ABORT", key="abort_button")]]
    window = sg.Window("Command run window", layout, modal=True, element_justification="center")
    thread_is_finished = False  # make sure to not continiously trigger updates while user decides to click OK
    while True:
        event, values = window.read(timeout=1)

        if event == "Exit" or event == sg.WIN_CLOSED or event == "ok_button" or event == "abort_button":
            thread_termination.request_exit()
            # try notify the user of their doing
            try:
                with open(log_file_loc) as f:
                    text = window.Element("multiline")
                    text.update(f.read() + os.linesep + command_thread.final_error)
            except Exception as e:
                print(f"failed to update error message with: {str(e)}")
            break
        elif not command_thread.is_alive() and not thread_is_finished:
            # try notify the user of their doing
            try:
                with open(log_file_loc) as f:
                    text = window.Element("multiline")
                    text.update(f.read() + os.linesep + command_thread.final_error)
            except Exception as e:
                print(f"failed to update error message with: {str(e)}")
            window.Element("ok_button").update(disabled=False)
            window.Element("abort_button").update(disabled=True)
            thread_is_finished = True

        # otherwise the file is closed
        elif command_thread.is_alive():
            with open(log_file_loc) as f:
                text = window.Element("multiline")
                text.update(f.read())
    try:
        temp.close()
        os.unlink(temp.name)
    except Exception as e:
        print(f"failed to delete temp file with: {str(e)}")
        pass
    window.close()


def get_command(values):
    arguments = ["-f"]

    if values["mpt_tabs"] == "Pairwise distance":
        arguments += ["distances", "-t", values["tgf_folder_d"], "-o", values["output_d"]]
        output_dir = values["output_d"]

    elif values["mpt_tabs"] == "Pairwise mutation":
        arguments += ["pairwise_mutation", "-af", values["allele_md"], "-o", values["output_md"]]
        if values["predict_file_md"] is True:
            arguments.append("-pf")
        if values["distance_md"] != '':
            arguments.extend(["-df", values["distance_md"]])
        output_dir = values["output_md"]

    elif values["mpt_tabs"] == "Pedigree mutation":

        arguments += ["pedigree_mutation", "-af", values['allele_pmg'], "-t", values['tgf_folder_pmg'], "-o",
                      values['output_pmg']]

        mm_arg = values["minimum_mutations_pmg"]
        if mm_arg != '':
            arguments += ["-mm", mm_arg]

        output_dir = values['output_pmg']

    elif values["mpt_tabs"] == "Dendrograms":
        arguments += ["dendrograms", "-fm", values['full_marker_dp'], '-o', values['output_dp']]
        marker_rate_file = values['marker_rate_dp']
        if marker_rate_file != '':
            arguments += ['-mr', values['marker_rate_dp']]

        if values['clusters_dp'] != "":
            arguments += ["-c", values['clusters_dp']]
        output_dir = values['output_dp']

    elif values["mpt_tabs"] == "Full analysis":
        arguments += ["all", "-t", values['tgf_folder_all'], '-af', values['allele_all'], '-o', values['output_all']]
        marker_rate_file = values['marker_rate_all']
        if marker_rate_file != '':
            arguments += ['-mr', values['marker_rate_all']]

        mm_arg = values["minimum_mutations_all"]
        if mm_arg != '':
            arguments += ["-mm", mm_arg]

        if values['clusters_all'] != "":
            arguments += ["-c", values['clusters_dp']]

        if values["predict_file_all"] is True:
            arguments.append("-pf")
        output_dir = values['output_all']
    elif values["mpt_tabs"] == "Prediction":
        arguments += ["prediction", "-i", values["input_pr"], "-o", values["output_pr"]]
        if values["custom_model_pr"] != '':
            if values["training_file_pr"] != '':
                arguments.extend(["-tf", values["training_file_pr"]])
            arguments.extend(["-m", values["custom_model_pr"]])
        else:
            arguments.extend(["-pm", values["model_choice_pr"]])
        if values["plots_pr"] is True:
            arguments.append("-p")
        output_dir = values["output_pr"]

    else:
        raise ValueError(f"Got tab: {values['mpt_tabs']}. There is no command for this tab, bad programmer")

    # request unique log name that is know at execution time
    log_name = main.get_log_file_name("no_value", output_dir)

    arguments.insert(1, log_name)
    arguments.insert(1, "-ln")

    if output_dir != '':
        output_dir = Path(output_dir)
        log_file = output_dir / log_name

        try:
            os.mkdir(output_dir)
        except FileExistsError:
            pass
        open(log_file, "w").close()
    else:
        log_file = None

    return arguments, log_file


if __name__ == '__main__':
    mpt_gui()
