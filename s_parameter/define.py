#!/usr/bin/env python3

import pathlib
import re

from format_multiple_errors import format_multiple_errors
from uncertainties import UFloat


substitutions = {
    "Nf12": "NfTwelve",
    "10": "Ten",
    "0": "Zero",
    "1": "One",
    "2": "Two",
    "3": "Three",
    "4": "Four",
    "5": "Five",
    "6": "Six",
    "7": "Seven",
    "8": "Eight",
    "9": "Nine",
    ".": "Point",
    "_": "",
    " ": "",
}


def get_filename(file_obj):
    filename = pathlib.Path(file_obj.name).name
    if filename.endswith(".tex"):
        filename = filename[:-4]
    return filename


def sanitize(key):
    """
    LaTeX names can only contain alphanumeric characters.
    As such,
    we must strip out spaces,
    and replace symbols and numerals with words.
    """
    original_key = key
    for from_token, to_token in substitutions.items():
        key = key.replace(from_token, to_token)
    if bad_char := re.search("[^A-Za-z0-9]", key):
        message = (
            f"Unable to sanitize key {original_key}; character {bad_char.group()} "
            f"remains at index {bad_char.start()} of {key}"
        )
        raise ValueError(message)
    return key


def reformat(value, **flags):
    if isinstance(value, UFloat):
        return f"{value:.02uSL}"
    if (
        isinstance(value, tuple)
        and isinstance(value[0], float)
        and all(
            isinstance(subvalue, float) or len(subvalue) == 2 for subvalue in value[1:]
        )
    ):
        return format_multiple_errors(*value, abbreviate=True, latex=True, **flags)
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return f"{value:.03e}"
    if isinstance(value, str):
        return value

    raise ValueError(f"I don't know how to format a {type(value)}.")


def define(name, value, **flags):
    return f"\\newcommand \\{sanitize(name)} {{{reformat(value, **flags)}}}"


def define_many(*values, **flags):
    return "\n".join(
        define(name, value, **flags) for name, value in dict(values).items()
    )


def name_ensemble(datum):
    return f"Nf{datum['Nf']}_mf{datum['mass']}_L{datum['Nx']}T{datum['Nt']}"
