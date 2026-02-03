#!/usr/bin/env python3

"""
Fit the vacuum polarisation function with a Padé fit
via the Ansatz f(q^2) = (b0 + b1 q^2)/(1 + c1 q^2 + c2 q^4);
see Eq. (9) of the paper.
"""

from argparse import ArgumentParser, FileType
import json

import numpy as np
from scipy.optimize import curve_fit

from ..io import dump_numpy, convert_types
from .pade import pade
from ..stats import jackknife_mean_variance


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_file")
    parser.add_argument("--output_file", type=FileType("w"), default="-")
    parser.add_argument("--renormalised", action="store_true")
    return parser.parse_args()


def max_q(length):
    # See below Eq. (9) of the paper
    return 2 * (2 * np.pi) / length


def get_s_parameter(param_samples):
    b0, b1, c1, _ = param_samples.T
    s_parameter_samples = b1 - b0 * c1
    return jackknife_mean_variance(s_parameter_samples)


def fit_single_pade(momentum_squared, vpf):
    popt, *_ = curve_fit(pade, momentum_squared, vpf)
    return popt


def fit_samples_pade(momentum_squared, vpf_samples):
    vpf_values, vpf_errors = jackknife_mean_variance(vpf_samples)

    result_samples = np.array(
        [fit_single_pade(momentum_squared, sample) for sample in vpf_samples]
    )
    fit_values, fit_errors = jackknife_mean_variance(result_samples)
    chisquare = (
        ((pade(momentum_squared, *fit_values) - vpf_values) / vpf_errors) ** 2
    ).sum()
    b0, b1, c1, c2 = zip(fit_values, fit_errors)

    return {
        "params": np.array([fit_values, fit_errors]),
        "b0": b0,
        "b1": b1,
        "c1": c1,
        "c2": c2,
        "chisquare": chisquare,
        "dof": len(vpf_samples) - 4,
        "S": get_s_parameter(result_samples),
    }


def fit_pade(data, key="renormalised_vpf_samples"):
    momentum = data["reordered_momentum"] * data["momentum_units"]
    momentum_filter = (
        (momentum[:, 0] <= max_q(data["Nx"]))
        & (momentum[:, 1] <= max_q(data["Ny"]))
        & (momentum[:, 2] <= max_q(data["Nz"]))
        & (momentum[:, 3] <= max_q(data["Nt"]))
        & (data["momentum_squared"] < 1)
    )
    filtered_momentum_squared = data["momentum_squared"][momentum_filter]
    result = {
        current: fit_samples_pade(
            filtered_momentum_squared, data[key][current][:, momentum_filter]
        )
        for current in ["OneLink", "Conserved"]
    }
    result["max_momentum_squared"] = filtered_momentum_squared.max()
    return result


def main():
    args = get_args()
    with open(args.input_file, "r") as input_file:
        data = json.load(input_file, object_pairs_hook=convert_types)

    keys = {True: "renormalised_vpf_samples", False: "vpf_samples"}
    result = fit_pade(data, keys[args.renormalised])

    dump_numpy(
        {
            "pade_fit_result": result,
            **{key: data[key] for key in ["mass", "Nt", "Nx", "Ny", "Nz"]},
        },
        args.output_file,
    )


if __name__ == "__main__":
    main()
