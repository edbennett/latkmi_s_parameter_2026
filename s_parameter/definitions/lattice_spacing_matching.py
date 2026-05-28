#!/usr/bin/env python3

from argparse import ArgumentParser, FileType

import pandas as pd

from ..define import define
from ..stats import add_quadrature
from ..comparison.plot_infinite_volume_S_lsd import get_rho_mass


def get_args():
    parser = ArgumentParser()
    parser.add_argument("--chiral_fit_result_latkmi", required=True)
    parser.add_argument("--chiral_spectrum_lsd", required=True)
    parser.add_argument("--output_definitions", type=FileType("w"), default="-")
    return parser.parse_args()


def get_definitions(chiral_fit_result_latkmi, chiral_spectrum_lsd):
    value_latkmi_mass, error_latkmi_mass, systematic_error_latkmi_mass = get_rho_mass(
        chiral_fit_result_latkmi
    )
    lsd_subset = chiral_spectrum_lsd.query("Nf == 8")
    assert len(lsd_subset) == 1
    value_lsd_mass = lsd_subset.iloc[0].value_rho_mass
    error_lsd_mass = lsd_subset.iloc[0].error_rho_mass
    systematic_error_lsd_mass = lsd_subset.iloc[0].systematic_error_rho_mass

    value_ratio = value_lsd_mass / value_latkmi_mass
    error_ratio = value_ratio * add_quadrature(
        (value_lsd_mass, error_lsd_mass),
        (value_lsd_mass, systematic_error_lsd_mass),
        (value_latkmi_mass, error_latkmi_mass),
        (value_latkmi_mass, systematic_error_latkmi_mass),
    )

    return define("Lattice_Spacing_Ratio_LSD_LatKMI", (value_ratio, error_ratio))


def main():
    args = get_args()
    chiral_fit_result_latkmi = pd.read_csv(args.chiral_fit_result_latkmi, comment="#")
    chiral_spectrum_lsd = pd.read_csv(args.chiral_spectrum_lsd, comment="#")
    print(
        get_definitions(chiral_fit_result_latkmi, chiral_spectrum_lsd),
        file=args.output_definitions,
    )


if __name__ == "__main__":
    main()
