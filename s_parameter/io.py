#!/usr/bin/env python3

"""
Tools for getting data in and out of files.
"""

import json

import numpy as np


def _serialise_value(value):
    """
    Helper function for dump_numpy:
    convert Numpy arrays to lists,
    and recurse into other data structures to find them.
    """
    match value:
        case np.ndarray():
            return value.tolist()
        case tuple():
            return list(map(_serialise_value, value))
        case dict():
            return {
                key: _serialise_value(inner_value) for key, inner_value in value.items()
            }
        case list():
            return list(map(_serialise_value, value))
        case _:
            return value


def dump_numpy(data, file_object):
    """
    JSON doesn't support Numpy arrays.
    Serialise them as nested lists.
    """
    json.dump(
        {key: _serialise_value(value) for key, value in data.items()},
        file_object,
    )


def to_numpy_recursive(data):
    """
    Given a nested dict,
    find any lists and turn them into Numpy arrays.
    """
    result = {}
    for key, value in data.items():
        match value:
            case list():
                result[key] = np.array(value)
            case dict():
                result[key] = to_numpy_recursive(value)
            case _:
                result[key] = value

    return result


def convert_types(items):
    return {
        int(k) if k.isdigit() else k: np.array(v) if isinstance(v, list) else v
        for k, v in items
    }


def read_numpy(input_filename):
    with open(input_filename, "r") as file_object:
        data = json.load(file_object)

    return to_numpy_recursive(data)


def read_numpy_optional(filename):
    if filename is None:
        return None
    return read_numpy(filename)


def read_Z_A(filename):
    with open(filename, "r") as file_object:
        data = json.load(file_object)

    data["Z_A_error"] = data["Z_A"][1]
    data["Z_A"] = data["Z_A"][0]
    data["Z_A_eff_error"] = np.array(data["Z_A_eff"][1])
    data["Z_A_eff"] = np.array(data["Z_A_eff"][0])
    data["Z_A_samples"] = np.array(data["Z_A_samples"])
    return data


def get_samples(data, ensemble_mass, key):
    result = [
        datum["fit_result_samples"][key]
        for datum in data
        if key in datum["fit_result_samples"] and datum["mass"] == ensemble_mass
    ]
    if len(result) != 1:
        raise ValueError("Missing or duplicate data")

    return result[0]


def update_without_overwrite(dest, source):
    for key in source:
        if key in dest:
            if isinstance(dest[key], dict) and isinstance(source[key], dict):
                update_without_overwrite(dest[key], source[key])
            else:
                assert dest[key] == source[key]
        else:
            dest[key] = source[key]


def collate_ensembles(data):
    collated_data = {}
    for datum in data:
        descriptor_keys = ["Nt", "Nx", "Ny", "Nz", "mass"]
        descriptor = tuple(datum[key] for key in descriptor_keys)
        if descriptor in collated_data:
            update_without_overwrite(collated_data[descriptor], datum)
        else:
            collated_data[descriptor] = datum

    return list(collated_data.values())
