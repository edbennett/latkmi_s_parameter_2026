#!/usr/bin/env python3

"""
Extract momentum correlators and PBP from a full log.
"""

from argparse import ArgumentParser
from compression import zstd
from collections import defaultdict
import re
import json
import sys


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_filename")
    parser.add_argument("--output_file", default=None)
    parser.add_argument("--use_complex", action="store_true")
    parser.add_argument("--Nt", type=int, default=None)
    parser.add_argument("--Nx", type=int, default=None)
    parser.add_argument("--Ny", type=int, default=None)
    parser.add_argument("--Nz", type=int, default=None)
    parser.add_argument("--Nf", type=int, default=None)
    return parser.parse_args()


def get_direction(operator_name):
    if operator_name[-1] in "0123":
        return int(operator_name[-1])
    else:
        return None


def recursive_defaultdict():
    return defaultdict(recursive_defaultdict)


def read_data(file_object, use_complex=False):
    # Read momentum correlators into a nested dict,
    # data["correlator"],
    # with index levels
    # - channel (V or A)
    # - current (OneLink or Conserved)
    # - source position (as an index into a separate list)
    # - source direction (denoted in the paper as nu)
    # - sink direction (denoted in the paper as mu)
    # Leaf nodes are lists of lists, with index levels
    # - Trajectory/configuration
    # - Momentum
    # Additional data in result dict gives:
    # - source positions
    # - trajectory indices
    # - momentum units
    # - probe mass

    if use_complex:
        raise NotImplementedError("complex number support not yet implemented")

    data = recursive_defaultdict()
    data["source_positions"] = []
    data["trajectory_indices"] = []
    data["1LPBP"] = {}
    mass = None

    reading_header_block = False
    reading_data_block = False

    target_sources = [
        f"OneLinkCurrent{channel}{nu}" for channel in ["V", "A"] for nu in range(4)
    ]
    target_sinks = [
        f"{current}Current{channel}"
        for current in ["Conserved", "OneLink"]
        for channel in ["V", "A"]
    ]

    for line in file_object:
        if line.startswith("start trajectory "):
            trajectory_index = int(line.split()[2])
            data["trajectory_indices"].append(trajectory_index)

        if line.startswith("SOURCE:"):
            source_name = line.split()[1]
            if source_name in target_sources:
                reading_header_block = True
                source_direction = get_direction(source_name)
                source_channel = source_name[-2]
            if source_name == "OneLinkCurrent":
                reading_header_block = True
                source_direction = None
                source_channel = None

        if line.startswith("MASSES:"):
            current_mass = float(line.split()[1])
            if mass is not None:
                if current_mass != mass:
                    raise ValueError("Inconsistent masses!")
            mass = current_mass

        if reading_header_block:
            if line.startswith("END"):
                reading_header_block = False
                reading_data_block = False
            if line.startswith("SOURCE_POS:"):
                source_position = tuple(map(int, line.split()[1:]))
                if source_position not in data["source_positions"]:
                    data["source_positions"].append(source_position)

                source_position_index = data["source_positions"].index(source_position)
                if source_position_index not in data["1LPBP"]:
                    data["1LPBP"][source_position_index] = []

            if line.startswith("MOMENTUM_UNITS"):
                momentum_units = list(map(float, line.split()[1:]))
                if "momentum_units" in data:
                    assert momentum_units == data["momentum_units"]
                else:
                    data["momentum_units"] = momentum_units

            if line.startswith("SINK:"):
                sink = line.split()[1]
                if sink in target_sinks:
                    assert sink[-1] == source_channel
                    reading_data_block = True
                    reading_header_block = False

                    (current,) = re.match("(.*)Current", sink).groups()
                    momentum_direction = []

                    datum = data["correlator"][source_channel][current][
                        source_position_index
                    ][source_direction]
                    for mu in range(4):
                        if mu not in datum:
                            datum[mu] = []
                        datum[mu].append([])
                    continue

            if line.startswith("SINKS:") and line.split()[1] == "1LPBP":
                sink = "1LPBP"
                reading_data_block = True
                reading_header_block = False
                continue

        if reading_data_block:
            if line.startswith("END") or line.startswith("SINK:"):
                if "momentum_direction" in data:
                    assert data["momentum_direction"] == momentum_direction
                else:
                    data["momentum_direction"] = momentum_direction

            if line.startswith("END"):
                reading_header_block = False
                reading_data_block = False
                continue
            if line.startswith("SINK:"):
                sink = line.split()[1]
                if sink not in target_sinks:
                    reading_data_block = False
                    reading_header_block = True
                else:
                    assert sink[-1] == source_channel
                    (current,) = re.match("(.*)Current", sink).groups
                    momentum_direction = []
                    datum = data["correlator"][source_channel][current][
                        source_position_index
                    ][source_direction]
                    for mu in range(4):
                        if mu not in datum:
                            datum[mu] = []
                        datum[mu].append([])
                continue

            split_line = line.split()

            if sink == "1LPBP":
                data["1LPBP"][source_position_index].append(
                    list(map(float, split_line[1::2]))
                )
                continue

            momentum_direction.append(list(map(int, split_line[:4])))
            for sink_direction, correlator_value in enumerate(split_line[4::2]):
                datum[sink_direction][-1].append(float(correlator_value))

    return {
        "data": data,
        "mass": mass,
    }


def main():
    args = get_args()
    with zstd.open(args.input_filename, "rt") as file_object:
        data = read_data(file_object, args.use_complex)

    # TODO: Check consistency between correlator length and given Nt
    data.update(
        {"Nt": args.Nt, "Nx": args.Nx, "Ny": args.Ny, "Nz": args.Nz, "Nf": args.Nf}
    )

    if args.output_file is None:
        json.dump(data, sys.stdout)
    else:
        with zstd.open(args.output_file, "wt") as output_file:
            json.dump(data, output_file)


if __name__ == "__main__":
    main()
