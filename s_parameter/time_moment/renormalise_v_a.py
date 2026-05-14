#!/usr/bin/env python3

r"""
Multiply the vector and axial correlation functions
by appropriate factors of the renormalisation factor Z_A.
"""

from argparse import ArgumentParser, FileType
import json

import numpy as np

from ..io import dump_numpy, convert_types
from ..stats import jackknife_mean_variance, bin_data, sample_jackknife


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_file_v_a")
    parser.add_argument("input_file_Z_A")
    parser.add_argument("--output_file", type=FileType("w"), default="-")
    return parser.parse_args()


def structure_data(full_data, current):
    """
    Rearrange the nested dictionaries as parsed out of the input file
    into a format allowing easier downstream processing.
    """
    return {
        channel: np.array(
            [
                [
                    [
                        source_data[f"{current}Current{channel}{mu}"]
                        for source_data in cfg_data
                        if source_data["_source"] == f"OneLinkCurrent{channel}{mu}"
                        if f"{current}Current{channel}{mu}" in source_data
                    ]
                    for _, cfg_data in sorted(
                        full_data.items(), key=lambda x: int(x[0])
                    )
                ]
                for mu in range(4)
            ]
        )
        for channel in ["V", "A"]
    }


def fold(correlator):
    """
    Symmetrise a correlator about the midpoint in the temporal direction.
    """
    correlator_length = correlator.shape[-1]
    assert correlator_length % 2 == 0
    folded_length = correlator_length // 2
    folded_correlator = correlator[:, : folded_length + 1]
    folded_correlator[:, 1:folded_length] += correlator[:, -1:-folded_length:-1]
    folded_correlator[:, 1:folded_length] /= 2
    return folded_correlator


def multiply(vpf, Z_A):
    # Correlator has dimensions [N_samples, N_timeslices],
    # Z_A has dimensions [N_samples],
    # so transpose to get the correct broadcasting
    return (vpf.T * Z_A).T


def renormalise(v_a, Z_A, current):
    """
    Multiply each one-link current by a factor Z_A.
    (I.e. OneLink-OneLink gets Z_A^2,
    while OneLink-Conserved gets Z_A.)
    """
    powers = {"OneLink": 2, "Conserved": 1}
    return multiply(v_a, Z_A["Z_A_samples"] ** powers[current])


def process(raw_v_a, Z_A):
    """
    Extract the zero-momentum V and A correlators and renormalise them.
    """
    results = {}
    for current in ["OneLink", "Conserved"]:
        result = {}
        current_data = structure_data(raw_v_a["data"], current)

        for channel in ["V", "A"]:
            zero_momentum_correlator = current_data[channel][:3, :, :, :].mean(
                axis=(0, 2)
            )
            binned_correlator = bin_data(zero_momentum_correlator, Z_A["bin_size"])
            correlator_samples = sample_jackknife(binned_correlator)
            renormalised_samples = renormalise(correlator_samples, Z_A, current)
            folded_renormalised_samples = fold(renormalised_samples)

            for samples, label in [
                # (renormalised_samples, "_renormalised"),
                (folded_renormalised_samples, "_renormalised"),
                (correlator_samples, ""),
            ]:
                result[f"{channel}{label}_samples"] = samples
                result[f"{channel}{label}"] = jackknife_mean_variance(samples)

        v_minus_a_samples = (
            result["V_renormalised_samples"] - result["A_renormalised_samples"]
        )
        result["V-A_renormalised_samples"] = v_minus_a_samples
        result["V-A_renormalised"] = jackknife_mean_variance(v_minus_a_samples)

        results[current] = result

    return results


def check_consistent_metadata(v_a, Z_A, keys):
    for key in keys:
        assert v_a[key] == Z_A[key]


def main():
    args = get_args()
    with open(args.input_file_v_a, "r") as input_file_v_a:
        v_a = json.load(input_file_v_a, object_pairs_hook=convert_types)
    with open(args.input_file_Z_A, "r") as input_file_Z_A:
        Z_A = json.load(input_file_Z_A, object_pairs_hook=convert_types)

    common_metadata_keys = ["mass", "Nf", "Nt", "Nx", "Ny", "Nz"]
    check_consistent_metadata(v_a, Z_A, common_metadata_keys)

    result = process(v_a, Z_A)

    dump_numpy(
        {
            **result,
            **{key: Z_A[key] for key in common_metadata_keys + ["bin_size"]},
        },
        args.output_file,
    )


if __name__ == "__main__":
    main()
