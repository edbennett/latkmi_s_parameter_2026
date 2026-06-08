#!/usr/bin/env python3

from argparse import ArgumentParser, FileType

import numpy as np
from scipy.special import gamma

from ..io import read_numpy, dump_numpy
from ..stats import jackknife_mean_error


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_file")
    parser.add_argument("--output_file", type=FileType("w"), default="-")
    return parser.parse_args()


def compute_S_parameter_contribution(correlator, start_time=0):
    """
    Given the folded, zero-momentum projected V-A correlator,
    or jackknife samples thereof,
    compute the contribution of each correlator point to the S parameter.

    C.f. Eq. (20) of the paper,
    but we do not yet perform the summation.
    """

    # 1.0/4.0 = 2.0/(2.0*4.0)
    # numerator 2.0 = redundant counting of the following (1) and (2)
    # denominaator 8.0 = to cancel the following (1) and (2)
    # (1) 4-flavor --> 2-flavor (one doublet)
    # (2) staggered 4-source meas.
    # c.f. -1/2! is taken account in the following and not included in fac
    normalisation_factor = 4.0 * np.pi / 4.0

    moment_index = 1
    times = start_time + np.arange(correlator.shape[-1])
    return (
        normalisation_factor
        * ((-1) ** moment_index)
        / gamma(2 * moment_index + 1)
        * times ** (2 * moment_index)
        * correlator
    )


def get_S_parameter(data):
    results = {}

    for current in ["OneLink", "Conserved"]:
        v_minus_a = data[current]["V-A_renormalised_samples"]

        # It's convenient to have this in the same file as S_eff for later plots
        result = {
            "V-A_renormalised_samples": v_minus_a,
            "V-A_renormalised": data[current]["V-A_renormalised"],
        }

        # Contribution to the total S parameter from each time slice
        S_parameter_contrib = compute_S_parameter_contribution(v_minus_a)
        result["S_parameter_contrib_samples"] = S_parameter_contrib
        result["S_parameter_contrib"] = jackknife_mean_error(S_parameter_contrib)

        # Effective S parameter; plateaus to S at large t
        S_parameter_eff = np.cumsum(S_parameter_contrib, axis=1)
        result["S_parameter_eff_samples"] = S_parameter_eff
        result["S_parameter_eff"] = jackknife_mean_error(S_parameter_eff)

        results[current] = result

    return results


def main():
    args = get_args()
    data = read_numpy(args.input_file)

    result = get_S_parameter(data)
    dump_numpy(
        {
            **result,
            **{
                key: data[key]
                for key in ["mass", "Nt", "Nx", "Ny", "Nz", "Nf", "bin_size"]
            },
        },
        args.output_file,
    )


if __name__ == "__main__":
    main()
