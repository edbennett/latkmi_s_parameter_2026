#!/usr/bin/env python3

from argparse import ArgumentParser, FileType
from collections import defaultdict

from format_multiple_errors import format_multiple_errors
import pandas as pd

from ..io import read_numpy


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_files", metavar="input_file", nargs="+")
    parser.add_argument("--output_file", type=FileType("w"), default="-")
    return parser.parse_args()


def sort_values_by_key(a_dict):
    return [a_dict[key] for key in sorted(a_dict)]


def tabulate(data):
    dataframe_source = defaultdict(dict)

    for datum in data:
        datum_key = (datum["mass"], datum["Nx"], datum["Nt"])
        structured_datum = dataframe_source[datum_key]
        assert datum["Nx"] == datum["Ny"] and datum["Nx"] == datum["Nz"]
        structured_datum["$L$"] = datum["Nx"]
        structured_datum["$T$"] = datum["Nt"]
        structured_datum["$m_f$"] = str(datum["mass"])
        for channel, channel_label in [("rho", r"\rho"), ("a_1", "a_1")]:
            for observable, observable_label in [("mass", "M"), ("decay_const", "F")]:
                target_key = f"{channel}_{observable}"
                if target_key in datum["fit_result"]:
                    structured_datum[f"${observable_label}_{{{channel_label}}}$"] = (
                        format_multiple_errors(
                            *datum["fit_result"][f"{channel}_{observable}"],
                            abbreviate=True,
                            latex=True,
                        )
                    )

    column_order = [
        "$L$",
        "$T$",
        "$m_f$",
        r"$F_{\rho}$",
        r"$F_{a_1}$",
        r"$M_{\rho}$",
        r"$M_{a_1}$",
    ]
    return pd.DataFrame(sort_values_by_key(dataframe_source))[column_order].to_latex(
        index=False
    )


def main():
    args = get_args()
    data = [read_numpy(input_file) for input_file in args.input_files]
    print(tabulate(data), file=args.output_file)


if __name__ == "__main__":
    main()
