#!/usr/bin/env python3

"""
Compute the renormalization factor $Z_A$
following the prescription in Section 2 of the paper.
"""

from argparse import ArgumentParser, FileType
import json

import numpy as np

from ..io import dump_numpy
from ..stats import jackknife_mean_variance, sample_jackknife_ratio, bin_data


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_file")
    parser.add_argument("--output_file", type=FileType("w"), default="-")
    parser.add_argument("--tmin", type=int, default=0)
    parser.add_argument("--tmax", type=int, default=-1)
    parser.add_argument("--bin_size", type=int, default=1)
    return parser.parse_args()


def structure_data(data, channels):
    """
    Pull out the specified channels from the input data,
    and reshape them into Numpy arrays suitable for downstream processing.
    """
    result = {key: [] for key in channels}
    result["configurations"] = []
    for configuration_index, configuration_data in data.items():
        result["configurations"].append(int(configuration_index))
        for channel in channels:
            configuration_result = []
            for source in configuration_data:
                configuration_result.append(source[channel])

            result[channel].append(configuration_result)

    return {key: np.array(value) for key, value in result.items()}


def get_Z_A(data, tmin, tmax, bin_size):
    """
    Given data structured as read from the raw correlator output,
    restructure, bin, and jackknife,
    using the relevant channels to compute the effective Z_A,
    and from this fit the plateau region to obtain Z_A with uncertainties.
    """

    structured_data = structure_data(data, ["A4_1LINK", "ConservedA4"])
    binned_conserved_current = bin_data(structured_data["ConservedA4"], bin_size)
    binned_one_link_current = bin_data(structured_data["A4_1LINK"], bin_size)
    Z_A_eff_samples = sample_jackknife_ratio(
        binned_conserved_current, binned_one_link_current
    )

    # Include lower and upper bound, and allow negative indices
    Nt = Z_A_eff_samples.shape[-1]
    t_upper_bound = tmax % Nt + 1
    Z_A_eff = jackknife_mean_variance(Z_A_eff_samples)
    Z_A_samples = np.average(
        Z_A_eff_samples[:, tmin:t_upper_bound],
        axis=-1,
        weights=1 / Z_A_eff[1][tmin:t_upper_bound] ** 2,
    )
    return {
        "Z_A_samples": Z_A_samples,
        "Z_A_eff": Z_A_eff,
        "Z_A": jackknife_mean_variance(Z_A_samples),
    }


def main():
    args = get_args()
    with open(args.input_file, "r") as input_file:
        full_data = json.load(input_file)

    result = get_Z_A(full_data["data"], args.tmin, args.tmax, args.bin_size)
    dump_numpy(
        {
            **result,
            **{key: full_data[key] for key in ["mass", "Nf", "Nt", "Nx", "Ny", "Nz"]},
            "plateau_start": args.tmin,
            "plateau_end": args.tmax,
            "bin_size": args.bin_size,
        },
        args.output_file,
    )


if __name__ == "__main__":
    main()
