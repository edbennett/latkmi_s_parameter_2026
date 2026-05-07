#!/usr/bin/env python3

"""
Fit Z_A as a function of fermion mass using linear and quadratic Ansätze
"""

from argparse import ArgumentParser, FileType
import json

import numpy as np
from scipy.optimize import curve_fit

from ..fits import fit_form_linear, fit_form_quadratic


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_files", nargs="+", metavar="input_file")
    parser.add_argument("--output_file", type=FileType("w"), default="-")
    return parser.parse_args()


def read_single(filename):
    with open(filename, "r") as file_object:
        data = json.load(file_object)

    data["Z_A_error"] = data["Z_A"][1]
    data["Z_A"] = data["Z_A"][0]
    data["Z_A_eff_error"] = np.array(data["Z_A_eff"][1])
    data["Z_A_eff"] = np.array(data["Z_A_eff"][0])
    data["Z_A_samples"] = np.array(data["Z_A_samples"])
    return data


def read(filenames):
    return [read_single(filename) for filename in filenames]


def fit_single(data, fit_form):
    mf_values = np.array([ensemble["mass"] for ensemble in data])
    Z_A_values = np.array([ensemble["Z_A"] for ensemble in data])
    Z_A_errors = np.array([ensemble["Z_A_error"] for ensemble in data])
    popt, pcov, info, _, _ = curve_fit(
        fit_form, mf_values, Z_A_values, sigma=Z_A_errors, full_output=True
    )
    return {
        "Z_A_0": (popt[0], pcov[0][0] ** 0.5),
        "c": (popt[1], pcov[1][1] ** 0.5),
        "chisquare": (info["fvec"] ** 2).sum(),
        "dof": len(mf_values) - len(popt),
        "mf_max": mf_values.max(),
    }


def fit(data):
    results = {
        "linear": fit_single(data, fit_form_linear),
        "quadratic": fit_single(
            sorted(data, key=lambda d: d["mass"])[:3],
            fit_form_quadratic,
        ),
    }
    Z_A_0_linear = results["linear"]["Z_A_0"][0]
    Z_A_0_quadratic = results["quadratic"]["Z_A_0"][0]
    central_result = (Z_A_0_linear + Z_A_0_quadratic) / 2
    central_statistical = 0.5 * (
        abs(results["linear"]["Z_A_0"][1]) + abs(results["quadratic"]["Z_A_0"][1])
    )
    central_error = abs(Z_A_0_linear - Z_A_0_quadratic) / 2
    results["central"] = {"Z_A_0": (central_result, central_statistical, central_error)}
    return results


def main():
    args = get_args()
    data = read(args.input_files)
    result = fit(data)
    json.dump(result, args.output_file)


if __name__ == "__main__":
    main()
