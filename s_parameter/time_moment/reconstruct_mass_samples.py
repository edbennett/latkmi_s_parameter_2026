#!/usr/bin/env python3

from argparse import ArgumentParser, FileType

import pandas as pd

from ..io import read_numpy, dump_numpy
from ..stats import generate_jackknife


def get_args():
    parser = ArgumentParser()
    parser.add_argument(
        "match_data",
        help="Data from ensemble to produce sample for",
    )
    parser.add_argument(
        "--spectrum_data",
        help="CSV of spectrum from previous work to base samples on",
        required=True,
    )
    parser.add_argument("--output_file", type=FileType("w"), default="-")
    return parser.parse_args()


def generate_mass_samples(match_datum, spectrum_data):
    result = {}
    keys = {"ma1": "a_1_mass", "mrho": "rho_mass"}
    matching_spectra = (
        spectrum_data.query(
            "Nf == 8 & beta == 3.8 & mf == {mass}".format(**match_datum)
        )
        .dropna(subset=[f"value_{key}" for key in keys])
        .sort_values(by="L", ascending=False)
    )
    if len(matching_spectra) == 0:
        # No matching data; write null data into the file
        return {key: None for key in keys.values()}

    ensemble_spectrum = matching_spectra.iloc[0]

    for old_key, new_key in keys.items():
        value = ensemble_spectrum[f"value_{old_key}"]
        error = ensemble_spectrum[f"error_{old_key}"]
        result[new_key] = generate_jackknife(
            value,
            error,
            match_datum,
            match_datum["fit_result_samples"]["osc_mass"].shape,
        )

    return result


def main():
    args = get_args()
    match_datum = read_numpy(args.match_data)
    spectrum_data = pd.read_csv(args.spectrum_data)
    result = generate_mass_samples(match_datum, spectrum_data)

    dump_numpy(
        {
            **result,
            **{
                key: match_datum[key]
                for key in ["mass", "Nt", "Nx", "Ny", "Nz", "bin_size"]
            },
            "source": args.spectrum_data,
        },
        args.output_file,
    )


if __name__ == "__main__":
    main()
