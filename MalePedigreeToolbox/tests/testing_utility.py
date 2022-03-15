import os
import shutil
from pathlib import Path
from collections import defaultdict

from MalePedigreeToolbox import main

CURRENT_DIR = Path(__file__).resolve().parent
TEST_FILE_DIR = CURRENT_DIR / 'test_files'
TEMP_OUT_DIR = CURRENT_DIR / 'temp_out'

try:
    os.mkdir(TEMP_OUT_DIR)
except IOError:
    pass


def clean_temp_out():
    for file in Path(TEMP_OUT_DIR / 'test_outdir').rglob("*"):
        if file.suffix == ".log" or file.is_dir():
            continue
        else:
            os.remove(file)


def clean_log_files():
    # make sure all files are closed and we can safely delete
    handlers = main.LOG.handlers[:]
    for handler in handlers:
        handler.close()
        main.LOG.removeHandler(handler)
    shutil.rmtree(TEMP_OUT_DIR / 'test_outdir')


def create_temp_out():
    if not TEMP_OUT_DIR.exists():
        os.mkdir(TEMP_OUT_DIR)


def run_command(command):
    in_wrapped = False
    arguments = []
    current_value = ""
    for char in command:
        if char == '"':
            if in_wrapped:
                if current_value != "":
                    arguments.append(current_value)
                    current_value = ""
                in_wrapped = False
            else:
                in_wrapped = True
        elif char == ' ' and not in_wrapped:
            if current_value != "":
                arguments.append(current_value)
                current_value = ""
        else:
            current_value += char
    arguments = arguments[1:]
    main.main(*arguments)


def confirm_files_exist(*expected_files):
    for file in expected_files:
        if not Path(file).exists():
            print(f"{file} not in expected_files")
            return False
    return True


def confirm_lines_equal(created_file, expected_file):
    # compare 2 files, line order does not have to be guaranteed
    expected_dict = defaultdict(int)
    with open(expected_file) as f:
        for line in f:
            expected_dict[line.strip()] += 1
    created_dict = defaultdict(int)
    with open(created_file) as f:
        for line in f:
            created_dict[line.strip()] += 1

    for entry, value in expected_dict.items():
        if entry not in created_dict:
            print(f"{entry} not present in file")
            return False
        if created_dict[entry] != value:
            print(f"{entry} is repeated {created_dict[entry]}, while expected {value}")
            return True
        created_dict.pop(entry)
    return len(created_dict) == 0


def confirm_columns_equal(created_file, expected_file, columns, sep=","):
    # confirm that specified columns in 2 files are the same as expected or all specified columns are excluded from
    # comparisson the order of lines in the file is not important
    expected_dict = defaultdict(int)
    with open(expected_file) as f:
        for line in f:
            compare_line = _get_requested_columns(line.strip(), sep, columns)
            expected_dict[compare_line] += 1
    created_dict = defaultdict(int)
    with open(created_file) as f:
        for line in f:
            compare_line = _get_requested_columns(line.strip(), sep, columns)
            created_dict[compare_line] += 1

    for entry, value in expected_dict.items():
        if entry not in created_dict:
            print(f"{entry} not present in file")
            return False
        if created_dict[entry] != value:
            print(f"{entry} is repeated {created_dict[entry]}, while expected {value}")
            return True
        created_dict.pop(entry)
    return len(created_dict) == 0


def _get_requested_columns(line, separator, ignore_columns):
    values = line.split(separator)
    ignore_columns = set(ignore_columns)
    wanted_values = []
    for index, value in enumerate(values):
        if index in ignore_columns:
            continue
        wanted_values.append(value)
    return f'{separator}'.join(wanted_values)
