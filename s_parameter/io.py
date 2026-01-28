#!/usr/bin/env python3

"""
Tools for getting data in and out of files.
"""

import json

import numpy as np


def serialise_value(value):
    match value:
        case np.ndarray():
            return value.tolist()
        case tuple():
            return list(map(serialise_value, value))
        case _:
            return value


def dump_numpy(data, file_object):
    """
    JSON doesn't support Numpy arrays.
    Serialise them as nested lists.
    """
    json.dump(
        {key: serialise_value(value) for key, value in data.items()},
        file_object,
    )


def read_Z_A(filename):
    with open(filename, "r") as file_object:
        data = json.load(file_object)

    data["Z_A_error"] = data["Z_A"][1]
    data["Z_A"] = data["Z_A"][0]
    data["Z_A_eff_error"] = np.array(data["Z_A_eff"][1])
    data["Z_A_eff"] = np.array(data["Z_A_eff"][0])
    data["Z_A_samples"] = np.array(data["Z_A_samples"])
    return data
