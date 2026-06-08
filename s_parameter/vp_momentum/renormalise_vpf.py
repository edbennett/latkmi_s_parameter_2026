#!/usr/bin/env python3

r"""
Multiply the vacuum polarisation function $\Pi^{V-A}(q^2)$
by the renormalisation factor Z_A,
per e.g. Eq. (5) of the paper.
"""

from argparse import ArgumentParser, FileType
import json

from ..io import dump_numpy, convert_types
from ..stats import jackknife_mean_error


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_file_vpf")
    parser.add_argument("input_file_Z_A")
    parser.add_argument("--output_file", type=FileType("w"), default="-")
    return parser.parse_args()


def multiply(vpf, Z_A):
    # VPF has dimensions [N_samples, N_momenta],
    # Z_A has dimensions [N_samples],
    # so transpose to get the correct broadcasting
    return (vpf.T * Z_A).T


def renormalise(vpf, Z_A):
    """
    Multiply each one-link current by a factor Z_A.
    (I.e. OneLink-OneLink gets Z_A^2,
    while OneLink-Conserved gets Z_A.)
    """
    powers = {
        "OneLink": 2,
        "Conserved": 1,
    }
    renormalised_vpf_samples = {
        current: multiply(samples, Z_A["Z_A_samples"] ** powers[current])
        for current, samples in vpf["vpf_samples"].items()
    }
    renormalised_vpf = {
        current: jackknife_mean_error(samples)
        for current, samples in renormalised_vpf_samples.items()
    }
    return {
        "renormalised_vpf_samples": renormalised_vpf_samples,
        "renormalised_vpf": renormalised_vpf,
    }


def main():
    args = get_args()
    with open(args.input_file_vpf, "r") as input_file_vpf:
        vpf = json.load(input_file_vpf, object_pairs_hook=convert_types)
    with open(args.input_file_Z_A, "r") as input_file_Z_A:
        Z_A = json.load(input_file_Z_A, object_pairs_hook=convert_types)

    result = renormalise(vpf, Z_A)
    dump_numpy(
        {
            **result,
            **{
                key: vpf[key]
                for key in [
                    "mass",
                    "Nt",
                    "Nx",
                    "Ny",
                    "Nz",
                    "Nf",
                    "bin_size",
                    "momentum_squared",
                    "reordered_momentum",
                    "momentum_units",
                ]
            },
        },
        args.output_file,
    )


if __name__ == "__main__":
    main()
