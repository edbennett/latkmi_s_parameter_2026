#!/usr/bin/env python3

from argparse import ArgumentParser, FileType

import pandas as pd

from ..define import define


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_data")
    parser.add_argument("--output_definitions", type=FileType("w"), default="-")
    return parser.parse_args()


def get_lightest_mrho(spectrum):
    Nf8_spectrum = spectrum.query("Nf == 8")
    lightest_mass = min(Nf8_spectrum.mf)
    lightest_results = Nf8_spectrum.query(f"mf == {lightest_mass}")
    assert len(lightest_results) == 1

    lightest_spectrum = lightest_results.iloc[0]
    ensemble_name = (
        f"Nf{lightest_spectrum['Nf']}_mf{lightest_spectrum['mf']}_"
        f"L{lightest_spectrum['L']}T{lightest_spectrum['T']}"
    )
    return define(
        f"Rho_Mass_{ensemble_name}",
        (lightest_spectrum["value_mrho"], lightest_spectrum["error_mrho"]),
    )


def get_mpi_over_mrho_ratio_bounds(spectrum):
    Nf8_spectrum = spectrum.query("Nf == 8")
    mpi_over_mrho = Nf8_spectrum["value_mpi"] / Nf8_spectrum["value_mrho"]
    return define(
        "Mpi_Over_Mrho_Range", f"{min(mpi_over_mrho):.02g} - {max(mpi_over_mrho):.02g}"
    )


def main():
    args = get_args()
    spectrum = pd.read_csv(args.input_data, comment="#")
    print(get_lightest_mrho(spectrum), file=args.output_definitions)
    print(get_mpi_over_mrho_ratio_bounds(spectrum), file=args.output_definitions)


if __name__ == "__main__":
    main()
