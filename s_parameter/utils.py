#!/usr/bin/env python3


def nested_get(datum, keys):
    """
    Where a key has multiple parts,
    treat each level as a key to nested dictionary.
    """
    # For convenience, allow a single level to be given without packaging it in a list
    if isinstance(keys, str):
        keys = [keys]

    for key in keys:
        datum = datum[key]

    return datum
