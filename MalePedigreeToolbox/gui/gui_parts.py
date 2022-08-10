import PySimpleGUI as sg


def SectionLabel(text):
    return sg.Text(text, justification="l", font="Arial 12 bold")


def TextLabel(text, size=(22, 1)):
    return sg.Text(
        f"{text}:",
        justification="l",
        size=size,
        font="Arial 10 bold"
    )


def Column(layout, scrollable=False):
    return sg.Column(
        layout,
        scrollable=scrollable,
        size=(649, 580),
        vertical_scroll_only=True
    )


def Frame(text, layout):
    return sg.Frame(
        text,
        layout=layout,
        font="Arial 10 bold",
        title_color="blue",
        relief="flat",
    )