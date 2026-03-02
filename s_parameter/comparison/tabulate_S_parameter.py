#!/usr/bin/env python3

from argparse import ArgumentParser, FileType
from functools import partial

from format_multiple_errors import format_multiple_errors

from ..io import read_numpy


def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_files", metavar="input_file", nargs="+")
    parser.add_argument("--output_file", type=FileType("w"), default="-")
    return parser.parse_args()


def get_with_attribute(data, spatial_size, temporal_size, mass, attribute):
    subset = [
        datum
        for datum in data
        if datum["Nx"] == spatial_size
        and datum["Ny"] == spatial_size
        and datum["Nz"] == spatial_size
        and datum["Nt"] == temporal_size
        and datum["mass"] == mass
        and attribute in datum
    ]
    if len(subset) == 0:
        raise ValueError("Datum not found")
    if len(subset) > 1:
        breakpoint()
        raise ValueError("Multiple results found")
    return subset[0]


def get_row(data, spatial_size, temporal_size, mass):
    vacuum_polarisation_datum = get_with_attribute(
        data, spatial_size, temporal_size, mass, "pade_fit_result"
    )
    time_moment_datum = get_with_attribute(
        data, spatial_size, temporal_size, mass, "S_infinite_t"
    )

    formatter = partial(
        format_multiple_errors,
        abbreviate=True,
        length_control="central",
        significant_figures=3,
    )
    vacuum_polarisation_S = formatter(
        *vacuum_polarisation_datum["pade_fit_result"]["Conserved"]["S"]
    )
    time_moment_S = formatter(*time_moment_datum["S_infinite_t"])

    row_data = [spatial_size, temporal_size, mass, vacuum_polarisation_S, time_moment_S]
    return " & ".join(map(str, row_data)) + r" \\"


def tabulate(data):
    ensembles = sorted(
        sorted(
            set([(datum["Nx"], datum["Nt"], datum["mass"]) for datum in data]),
            key=lambda x: x[0],  # Reverse sort by volume first
            reverse=True,
        ),
        key=lambda x: x[2],  # Then sort by mass
    )
    header = "\n".join(
        [
            r"\begin{tabular}{lllll}",
            r"\toprule",
            " & ".join(
                ["$L$", "$T$", "$m_f$", r"$S_{\textrm{VP-Mom}}$", r"$S_{\textrm{TM}}$"]
            )
            + r" \\",
            r"\midrule",
        ]
    )
    footer = "\n".join([r"\bottomrule", r"\end{tabular}"])

    content = [get_row(data, *ensemble) for ensemble in ensembles]
    return "\n".join([header, *content, footer])


def main():
    args = get_args()
    data = [read_numpy(input_file) for input_file in args.input_files]
    print(tabulate(data), file=args.output_file)


if __name__ == "__main__":
    main()
