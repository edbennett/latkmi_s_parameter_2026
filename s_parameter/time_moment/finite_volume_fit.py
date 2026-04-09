#!/usr/bin/env python3

from argparse import ArgumentParser, FileType

import numpy as np
from scipy.optimize import curve_fit

from ..io import read_numpy, dump_numpy


METADATA_KEYS = ["mass", "Nt", "Nx", "Ny", "Nz", "bin_size"]


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_files", metavar="INPUT_FILE", nargs="+")
    parser.add_argument("--output_file", default="-", type=FileType("w"))
    return parser.parse_args()


def get_fit_form(data):
    masses = []
    mass_indices = []

    for datum in data:
        if datum["mass"] not in masses:
            masses.append(datum["mass"])
        mass_indices.append(masses.index(datum["mass"]))

    mass_indices = np.array(mass_indices)

    def fit_form(volume_factor, C, *S_infinity):
        return C * volume_factor + np.array(S_infinity)[mass_indices]

    return masses, fit_form


def get_samples(data):
    breakpoint()
    (num_samples,) = set(len(datum["delta_fv_S_samples"]) for datum in data)
    errors = [datum["S_infinite_t_samples"].std() for datum in data]
    for sample_idx in range(num_samples):
        yield (
            [datum["delta_fv_S_samples"][sample_idx] for datum in data],
            [datum["S_infinite_t_samples"][sample_idx] for datum in data],
            errors,
        )


def fit(data):
    mass_ordering, fit_form = get_fit_form(data)
    starting_guess = [1.0] + [1.0 for _ in mass_ordering]

    fit_mean, fit_covariance = curve_fit(
        fit_form,
        [datum["delta_fv_S"][0] for datum in data],
        [datum["S_infinite_t"][0] for datum in data],
        sigma=[datum["S_infinite_t"][0] for datum in data],
        p0=starting_guess,
    )
    const_coefficient, *S_infinite_volume = zip(fit_mean, fit_covariance.diagonal())
    return {
        "masses": mass_ordering,
        "C": const_coefficient,
        "S_infinite_volume": S_infinite_volume,
    }


def main():
    args = get_args()
    data = [read_numpy(input_file) for input_file in args.input_files]
    result = fit(data)
    dump_numpy(
        {
            "fit_result": result,
            "source_metadata": [
                {key: datum[key] for key in METADATA_KEYS} for datum in data
            ],
        },
        args.output_file,
    )


if __name__ == "__main__":
    main()
