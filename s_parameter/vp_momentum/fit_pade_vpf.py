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
from .pade import pade, get_momentum_filter
from ..stats import jackknife_mean_variance


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_file")
    parser.add_argument("--output_file", type=FileType("w"), default="-")
    parser.add_argument("--renormalised", action="store_true")
    parser.add_argument("--upper_bound", choices=["1", "max2", "max3", "final"])
    return parser.parse_args()


def max_q(length):
    # See below Eq. (9) of the paper
    return 2 * (2 * np.pi) / length


def get_s_parameter(param_samples):
    b0, b1, c1, _ = param_samples.T
    s_parameter_samples = (b1 - b0 * c1) * 2 * np.pi
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
    chisquare_samples = (
        (
            (pade(momentum_squared[:, np.newaxis], *result_samples.T).T - vpf_values)
            / vpf_errors
        )
        ** 2
    ).sum(axis=1)
    chisquare = jackknife_mean_variance(chisquare_samples)
    b0, b1, c1, c2 = zip(fit_values, fit_errors)

    return {
        "params": np.array([fit_values, fit_errors]),
        "b0": b0,
        "b1": b1,
        "c1": c1,
        "c2": c2,
        "chisquare_samples": chisquare_samples,
        "chisquare": chisquare,
        "dof": len(vpf_samples) - 4,
        "S": get_s_parameter(result_samples),
    }


def fit_pade(data, key="renormalised_vpf_samples", momentum_squared_upper_bound=1):
    momentum_filter = get_momentum_filter(
        data["reordered_momentum"],
        [data[key] for key in ["Nx", "Ny", "Nz", "Nt"]],
        momentum_squared_upper_bound,
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


def get_upper_bound(data, key):
    """
    Implements Eq. (B3) of the paper.
    """
    Nt, Nx, Ny, Nz = data["Nt"], data["Nx"], data["Ny"], data["Nz"]
    assert Nx == Ny and Nx == Nz
    common_factor = 3 * (2 * np.pi / Nx) ** 2 + (2 * np.pi / Nt) ** 2
    q_max2 = 2**2 * common_factor
    q_max3 = 3**2 * common_factor

    return {
        "1": 1,
        "max2": q_max2,
        "max3": q_max3,
        "final": min(1, q_max2),
    }[key]


def main():
    args = get_args()
    with open(args.input_file, "r") as input_file:
        data = json.load(input_file, object_pairs_hook=convert_types)

    keys = {True: "renormalised_vpf_samples", False: "vpf_samples"}
    fit_upper_bound = get_upper_bound(data, args.upper_bound)
    result = fit_pade(data, keys[args.renormalised], fit_upper_bound)

    dump_numpy(
        {
            "pade_fit_result": result,
            "upper_bound": args.upper_bound,
            **{key: data[key] for key in ["mass", "Nt", "Nx", "Ny", "Nz"]},
        },
        args.output_file,
    )


if __name__ == "__main__":
    main()
