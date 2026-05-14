#!/usr/bin/env python3

from argparse import ArgumentParser, FileType

import numpy as np

from ..io import read_numpy, dump_numpy


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_files", nargs="+")
    parser.add_argument("--output_file", default="-", type=FileType("w"))
    return parser.parse_args()


def get(datum, channel):
    return datum["pade_fit_result"][channel]


def add_systematics(result, max2, max3):
    for channel in ["OneLink", "Conserved"]:
        for key in "S", "b0", "b1", "c1", "c2", "params":
            systematic = np.abs(get(max2, channel)[key][0] - get(max3, channel)[key][0])
            target = get(result, channel)
            target[key] = np.array([*target[key], systematic])


def check_consistency(data):
    keys = ["Nt", "Nx", "Ny", "Nz", "Nf", "mass"]
    assert all(all(data[0][key] == datum[key] for datum in data) for key in keys)


def upper_bound(datum):
    return datum["pade_fit_result"]["max_momentum_squared"]


def combine(raw_data):
    data = {datum["upper_bound"]: datum for datum in raw_data}

    assert all(key in data for key in ["max2", "max3"])

    if "1" in data and upper_bound(data["1"]) < upper_bound(data["max2"]):
        result = data["1"]
    else:
        result = data["max2"]

    add_systematics(result, data["max2"], data["max3"])
    return result


def main():
    args = get_args()
    data = [read_numpy(filename) for filename in args.input_files]
    check_consistency(data)
    result = combine(data)

    dump_numpy(result, args.output_file)


if __name__ == "__main__":
    main()
