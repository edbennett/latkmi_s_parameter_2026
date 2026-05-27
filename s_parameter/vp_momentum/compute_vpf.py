#!/usr/bin/env python3

r"""
Compute the vacuum polarisation function $\Pi^{V-A}(q^2)$
following the prescription in Section 2 of the paper.
"""

from argparse import ArgumentParser, FileType
from compression import zstd
from functools import partial
import itertools
import json

import numpy as np

from ..io import dump_numpy, convert_types
from ..stats import jackknife_mean_variance, sample_jackknife, bin_data


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_file")
    parser.add_argument("--output_file", type=FileType("w"), default="-")
    parser.add_argument("--bin_size", type=int, default=1)
    return parser.parse_args()


def empty_tensor_dict():
    return {
        source_direction: {sink_direction: [] for sink_direction in range(4)}
        for source_direction in range(4)
    }


def reorder_momentum(momentum):
    abs_momentum = np.abs(momentum)
    return np.concatenate(
        [np.sort(abs_momentum[:, :3], axis=1), abs_momentum[:, 3:]],
        axis=1,
    )


def project(data):
    momentum = data["momentum_direction"] * data["momentum_units"]
    # We will index on momentum_squared,
    # so coerce all extremely close numbers to be the same
    # In principle this introduces a very small systematic effect,
    # but at least ten orders of magnitude smaller than the signal
    # If this affects our results,
    # we are working far too close to machine precision to trust them either way.
    momentum_squared = (momentum**2).sum(axis=1).round(15)

    sort_index = np.argsort(momentum_squared)
    unique_momentum_squared, momentum_groups, momentum_group_count = np.unique(
        momentum_squared[sort_index],
        return_index=True,
        return_counts=True,
    )
    unique_momentum = reorder_momentum(momentum[sort_index][momentum_groups])

    projection = {}

    def project_single(channel, current, source_position_index):
        # Use a closure to reduce indentation level
        # because we need so much surrounding state
        source_position = tuple(data["source_positions"][source_position_index])
        projection[channel][current][source_position] = (dest_data := {})
        source_data = data["correlator"][channel][current][source_position_index]

        for key in ["PiA", "PiB"]:
            for subtraction in [""] + (
                ["_subtracted"] if current == "Conserved" else []
            ):
                dest_data[f"{key}{subtraction}"] = 0

        # mu: sink  direction; nu: source direction
        for mu, nu in itertools.product(range(4), range(4)):
            momentum_factor = momentum[:, mu] * momentum[:, nu] / momentum_squared
            pi_b = (
                np.add.reduceat(
                    (source_data[nu][mu] * momentum_factor)[:, sort_index],
                    momentum_groups,
                    axis=1,
                )
                / momentum_group_count
            )
            dest_data["PiB"] += pi_b

            if current == "Conserved":
                dest_data["PiB_subtracted"] += pi_b
                if mu == nu:
                    source_pbp = data["1LPBP"][source_position_index][:, mu]
                    dest_data["PiB_subtracted"] += (
                        np.add.reduceat(
                            source_pbp[:, np.newaxis] * momentum_factor,
                            momentum_groups,
                            axis=1,
                        )
                        / momentum_group_count
                    )

        # Trace
        for nu in range(4):
            pi_a = (
                np.add.reduceat(
                    source_data[nu][nu][:, sort_index], momentum_groups, axis=1
                )
                / momentum_group_count
            )
            dest_data["PiA"] += pi_a
            if current == "Conserved":
                dest_data["PiA_subtracted"] += (
                    pi_a + data["1LPBP"][source_position_index][:, nu, np.newaxis]
                )

    def normalise_single(channel, current, source_position):
        dest_data = projection[channel][current][source_position]
        for observable in ["PiA", "PiB", "PiA_subtracted", "PiB_subtracted"]:
            if observable in dest_data:
                dest_data[observable] /= 4

    for channel in ["V", "A"]:
        projection[channel] = {}
        for current in ["OneLink", "Conserved"]:
            projection[channel][current] = {}
            for source_position_index, source_position in enumerate(
                data["source_positions"]
            ):
                project_single(channel, current, source_position_index)
                normalise_single(channel, current, tuple(source_position.tolist()))

    return projection, unique_momentum, unique_momentum_squared


def subtract_V_A(data):
    result = {}
    for current in ["OneLink", "Conserved"]:
        result[current] = {}
        for source_position in data["V"][current].keys():
            result[current][source_position] = {}
            for current_type in ["PiA", "PiB"]:
                result[current][source_position][current_type] = (
                    data["V"][current][source_position][current_type]
                    - data["A"][current][source_position][current_type]
                )
    return result


def map_currents_sources(func, data):
    return {
        current: {
            source_position: func(source_data)
            for source_position, source_data in current_data.items()
        }
        for current, current_data in data.items()
    }


def compute_vpf(datum, momentum_squared):
    momentum_filter = momentum_squared != 0
    return (datum["PiA"][:, momentum_filter] - datum["PiB"][:, momentum_filter]) / 3


def combine_sources(data):
    return np.array(list(data.values())).swapaxes(0, 1)


def map_dict(func, data):
    return {key: func(value) for key, value in data.items()}


def get_vpf(data, bin_size=1):
    projected_data, momentum, momentum_squared = project(data)
    v_minus_a = subtract_V_A(projected_data)

    vpf_raw = map_currents_sources(
        partial(compute_vpf, momentum_squared=momentum_squared), v_minus_a
    )
    vpf_combined = map_dict(combine_sources, vpf_raw)
    vpf_bins = map_dict(partial(bin_data, bin_size=bin_size), vpf_combined)
    vpf_samples = map_dict(sample_jackknife, vpf_bins)
    vpf = map_dict(jackknife_mean_variance, vpf_samples)

    return {
        "vpf": vpf,
        "vpf_samples": vpf_samples,
        "momentum_squared": momentum_squared[momentum_squared != 0],
        "reordered_momentum": momentum[momentum_squared != 0],
        "momentum_units": data["momentum_units"],
    }


def main():
    args = get_args()
    with zstd.open(args.input_file, "rt") as input_file:
        full_data = json.load(input_file, object_pairs_hook=convert_types)

    result = get_vpf(full_data["data"], args.bin_size)
    dump_numpy(
        {
            **result,
            **{key: full_data[key] for key in ["mass", "Nf", "Nt", "Nx", "Ny", "Nz"]},
            "bin_size": args.bin_size,
        },
        args.output_file,
    )


if __name__ == "__main__":
    main()
