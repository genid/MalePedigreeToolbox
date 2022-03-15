#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Some utility functions

author: Bram
"""
import os
from sys import exit
from pathlib import Path
from typing import Union
import uuid


class MalePedigreeToolboxError(Exception):
    pass


def check_tgf_folder(
    folder_path: str
) -> Path:
    """Check if a given folder exists"""
    folder_path = _verify_str_path(folder_path)
    path = Path(folder_path)
    if not path.exists():
        raise MalePedigreeToolboxError(f"Provided tgf folder ({path}) does not exist.")
    if not path.is_dir():
        raise MalePedigreeToolboxError("Provided output directory is a file, please provide a directory.")
    all_tgf_files = list(path.glob("*.tgf"))
    if len(all_tgf_files) < 1:
        raise MalePedigreeToolboxError(f"At least 1 tgf files should be present in {path}. "
                                       f"But no tgf files where found.")
    return path


def check_file_create(
    file_path: str
) -> Path:
    """Check if a file can be created"""
    file_path = _verify_str_path(file_path)
    path = Path(file_path)
    if path.is_dir():
        raise MalePedigreeToolboxError("Provided file is a directory, please provide a file.")
    if not path.exists():
        # create empty file to check if there are writing permissons and the path is correct
        try:
            open(path, "w").close()
        except Exception as e:
            raise MalePedigreeToolboxError(f"Can not create file at {file_path}. Error code {str(e)}.")
    return path


def check_create_out_folder(
    path: str,
    force: bool = False
) -> Path:
    """Check that the output folder can be created and ask user to overwite if force is not used"""
    path = _verify_str_path(path)
    folder = Path(path).resolve()
    if folder.exists() and not force:
        print("The specified folder already exists. Are you sure you want to potentially overwrite its contents?"
              " Y / N")
        while True:
            val = input()
            val = val.upper()
            if val == "Y":
                if not folder.is_dir():
                    raise MalePedigreeToolboxError("Provided output directory is a file, please provide a directory.")
                return folder
            elif val == "N":
                exit("No valid output folder selected")
            else:
                print("Pleas enter Y or N")
    try:
        os.mkdir(folder)
    except Exception as e:
        # when forcing just ignore this error and write to the requested folder
        if force:
            return folder
        raise MalePedigreeToolboxError(f"Failed to create specified directory folder {folder}. Error code: {str(e)}")
    return folder


def check_in_file(
    file_path: str
) -> Path:
    """Check that the input file is valid and not a directory"""
    file_path = _verify_str_path(file_path)
    path = Path(file_path)
    if path.is_dir():
        raise MalePedigreeToolboxError("Provided file is a directory, please provide a file.")
    if not path.exists():
        raise MalePedigreeToolboxError(f"Specified input file does not exist. Can not find path: {path}")
    return path


def check_file_int(
    value: str
) -> Union[int, Path]:
    """See if the string is an int"""
    try:
        return int(value)
    except ValueError:
        pass
    return check_in_file(value)


def _verify_str_path(
    path: str
) -> str:
    """Verify that a path is not empty to not cause unexpected behaviour when wrapping a string with pathlib"""
    if path == '':
        raise MalePedigreeToolboxError(f"Empty path was provided. One required argument is missing a path.")
    return path


def sanitize_string(
    string: str
) -> str:
    """Sanitize strings to make sure that they are allowed as filenames"""
    new_name = ""
    for char in string:
        if char in ["#", "?", "/", "\\", "'", '"', ";", ":", "%"]:
            continue
        new_name += char
    if len(new_name) == 0:
        return str(uuid.uuid4())
    return new_name
